/* Ergo, version 3.6, a program for linear scaling electronic structure
 * calculations.
 * Copyright (C) 2017 Elias Rudberg, Emanuel H. Rubensson, Pawel Salek,
 * and Anastasia Kruchinina.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Primary academic reference:
 * Kohn−Sham Density Functional Theory Electronic Structure Calculations 
 * with Linearly Scaling Computational Time and Memory Usage,
 * Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek,
 * J. Chem. Theory Comput. 7, 340 (2011),
 * <http://dx.doi.org/10.1021/ct100611z>
 * 
 * For further information about Ergo, see <http://www.ergoscf.org>.
 */

/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/

/** @file xc_matrix_sparse.cc

    @brief The sparse XC matrix evaluator.

    @author: Pawel Salek <em>responsible</em>

    (c) Pawel Salek, pawsa@theochem.kth.se.
    2002.04.05

    This module evaluates DFT contribution KS matrix.
*/

#define WITH_PTHREAD 1
#if defined(WITH_PTHREAD)
#include <pthread.h>
static pthread_mutex_t dft_prop_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t dft_hicu_grid_init_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

#include "aos.h"
#include "integrator.h"
#include "sparse_matrix.h"
#include "xc_matrix_sparse.h"
#include "dft_common.h"
#include "grid_reader.h"
#include "output.h"
#include "utilities.h"
#include "xc_evaluators.h"

/* FIXME: remove this dependency */
#include "grid_hicu.h"

/* restrict hints should not be necessary... */
#if !defined(restrict)
#define restrict
#endif

BEGIN_NAMESPACE(Dft)


class XCEvaluator {
protected:
    const BasisInfoStruct& bisOrig;
    const IntegralInfo& integralInfo;
    const Molecule& mol;
    const Dft::GridParams& gss;
    std::vector<int> const & permutationHML;
    int* aoMap;
    BasisInfoStruct* bisPermuted;
    Dft::SparsePattern* pattern;
public:
    XCEvaluator(const BasisInfoStruct& bisOrig_,
                const IntegralInfo& integralInfo_,
                const Molecule& mol_,
                const Dft::GridParams& gss_,
                std::vector<int> const & permutationHML_,
                const symmMatrix& dens);
    ~XCEvaluator();
};

XCEvaluator::XCEvaluator(const BasisInfoStruct& bisOrig_,
                         const IntegralInfo& integralInfo_,
                         const Molecule& mol_,
                         const Dft::GridParams& gss_,
                         std::vector<int> const & permutationHML_,
                         const symmMatrix& dens)
    : bisOrig(bisOrig_), integralInfo(integralInfo_),
      mol(mol_), gss(gss_), permutationHML(permutationHML_),
      aoMap(new int[bisOrig.noOfBasisFuncs])
{
    /* We need to create our own permutation of shells. The
       permutation used in the matrix library permutes basis functions
       even in the middle of the basis function shell, breaking in
       this way simple mapping between shells and basis function
       blocks. We therefore generate here own permutation of
       shells. This is used to create own BasisInfoStruct used
       throughout. The final conversion from own format to the "Common
       Matrix Format" is done at the final assignment phase: the basis
       function indices are permuted with the inverse
       transformation. */

    int *shellMap = new int[bisOrig.noOfShells];

    Dft::setupShellMap(bisOrig, shellMap, aoMap);
    bisPermuted = bisOrig.permuteShells(shellMap, integralInfo);
    delete []shellMap;

    
    /* Force grid creation so that the sparse pattern is
       available... Maybe it is getSparsePattern's problem. */
    ErgoMolInfo molInfo(*bisPermuted, mol);
    pattern = new SparsePattern(*bisPermuted);

    /* When using HiCu grid generation we need a sparse density
       matrix. In order to create a sparse density matrix we need a
       SparsePattern. However, this is only needed when grid file does
       not exist, and only for the HiCu case. */
    SparseMatrix* dMatForGridGen = NULL;
    /* Use mutex lock here to make sure only one thread does this
       (there is anyway threading inside the HiCu code). */
    pthread_mutex_lock(&dft_hicu_grid_init_mutex);
    if(!grid_is_ready() && gss.gridType == Dft::GridParams::TYPE_HICU) {
      Dft::SparsePattern patternTmp(*bisPermuted);
      /* Generate sparse pattern. */
      grid_generate_sparse_pattern(*bisPermuted, gss.hicuParams.maxError, 
				   gss.hicuParams.box_size, gss.hicuParams.start_box_size_debug, patternTmp);
      /* Create dMat already here so it can be used by HiCu grid generator! */
      dMatForGridGen = new SparseMatrix(patternTmp, dens, aoMap, permutationHML);
    }
    else {
      // We do not really need a sparse density matrix in this case; we just create a dummy.
      Dft::SparsePattern patternTmp(*bisPermuted);
      dMatForGridGen = new SparseMatrix(patternTmp);
    }
    pthread_mutex_unlock(&dft_hicu_grid_init_mutex);

    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "getXC calling grid_open_full().");
    Dft::Matrix *mat = createGridMatrix(*dMatForGridGen);
    DftGridReader* rawgrid = grid_open_full(&molInfo, gss, pattern, 
                                            mat, *bisPermuted);
    delete mat;
    delete dMatForGridGen;

    grid_close(rawgrid);    
}

XCEvaluator::~XCEvaluator()
{
    delete bisPermuted;
    delete pattern;
    delete []aoMap;
}

class XCEvaluatorRestricted : public XCEvaluator {
  SparseMatrix *densityMatrix;
public:
  XCEvaluatorRestricted(const BasisInfoStruct& bisOrig_,
                        const IntegralInfo& integralInfo_,
                        const Molecule& mol_,
                        const Dft::GridParams& gss_,
                        std::vector<int> const & permutationHML_,
                        const symmMatrix& density)
    : XCEvaluator(bisOrig_, integralInfo_, mol_, gss_, permutationHML_,
                  density),
      densityMatrix(NULL)
  {
    densityMatrix = new SparseMatrix(*pattern, density, aoMap, permutationHML);
  }
  ~XCEvaluatorRestricted()
  {
    delete densityMatrix;
  }
      
  real getXC(int nElectrons, symmMatrix& xcm, real* xcEnergy,
             int nThreads) const;
};

/** Computes Fock matrix xcm corresponding to given density matrix dmat.
   fast version - uses memory bandwidth-efficient algorithm.
*/
real
XCEvaluatorRestricted::getXC(int nElectrons, symmMatrix& xcm,
                             real* xcEnergy, int nThreads) const
{
    real electrons;
    Util::TimeMeter tm;

    *xcEnergy = 0;
    sync_threads(false, nThreads);
    SparseMatrix excmat(*pattern);

    KsData<Dft::SparseMatrix> ds(&excmat, DFT_MAX_BLLEN);

    void (*cblda)(DftIntegratorBl* grid, real * restrict tmp, 
                  int bllen, int blstart, int blend,
                  KsData<Dft::SparseMatrix>* data)
      = xcCallbackLdaR<Dft::SparseMatrix,XCDistributorLda<Dft::SparseMatrix> >;
    void (*cbgga)(DftIntegratorBl* grid, real * restrict tmp, 
                  int bllen, int blstart, int blend,
                  KsData<Dft::SparseMatrix>* data)
      = xcCallbackGgaR<Dft::SparseMatrix,XCDistributorGga<Dft::SparseMatrix> >;

    electrons = integrate(1, &densityMatrix, *bisPermuted, mol, gss, nThreads,
                          (DftBlockCallback)
                          (selected_func->is_gga() ? cbgga : cblda),
                          &ds);

    pthread_mutex_lock(&dft_prop_mutex);
    *xcEnergy +=ds.energy;
    excmat.addSymmetrizedTo(xcm, aoMap, permutationHML);
    pthread_mutex_unlock(&dft_prop_mutex);


    if(nThreads<=1) {
        do_output(LOG_CAT_INFO, LOG_AREA_DFT, 
                  "Electrons: %11.7f %7.1g: xc energy %f (serial)", 
                  (double)electrons,
                  (double)((electrons-nElectrons)/nElectrons), 
                  (double)ds.energy);
        tm.print(LOG_AREA_DFT, __func__);
    }
    return electrons;
}

struct XcData {
    const XCEvaluatorRestricted* xcEvaluator;
    int nElectrons;
    symmMatrix* xcm;
    real xcEnergy;
    real el;
    int nThreads;
};

static void*
xcWorker(void *data)
{
  static const int XCWORKER_ERROR =0;
  struct XcData *d = (XcData*)data;
  try {
      d->el = d->xcEvaluator->getXC(d->nElectrons, *d->xcm, &d->xcEnergy,
                                    d->nThreads);
  } catch(const char *s) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, 
              "xcWorker thread caught an exception '%s'", s);
    return (void*)&XCWORKER_ERROR;
  } catch(const std::bad_alloc & e) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, 
              "xcWorker thread caught an exception '%s'", e.what());
    return (void*)&XCWORKER_ERROR;
  } catch(const std::runtime_error & e) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, 
              "xcWorker thread caught an exception '%s'", e.what());
    return (void*)&XCWORKER_ERROR;
  }  catch(...) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, 
              "xcWorker thread caught unexpected exception.");
    return (void*)&XCWORKER_ERROR;
  }
  return NULL;
}

real
getXC_mt(const BasisInfoStruct& bis, const IntegralInfo& integralInfo,
         const Molecule& mol, const Dft::GridParams& gss, int nElectrons,
         const symmMatrix& dens, symmMatrix& xcm, real* xcEnergy,
         std::vector<int> const & permutationHML)
{
    int i;
    real electrons;
    Util::TimeMeter tm;
    
    int nThreads = dft_get_num_threads();
    std::vector<XcData> data(nThreads);
    std::vector<pthread_t> pids(nThreads);

    XCEvaluatorRestricted xcEvaluator(bis, integralInfo, mol, gss,
                                      permutationHML, dens);

    if(nThreads == 1) {
        /* Do not create any threads at all to avoid stack allocation. */
        *xcEnergy = 0.0;
        electrons = xcEvaluator.getXC(nElectrons, xcm, xcEnergy, 1);
    } else {
        for(i=0; i<nThreads; i++) {
            data[i].xcEvaluator = &xcEvaluator;
            data[i].nElectrons = nElectrons;
            data[i].xcm   = &xcm;
            data[i].xcEnergy = 0.0;
            data[i].nThreads = nThreads;
            if(pthread_create(&pids[i], NULL, xcWorker, &data[i])) {
              do_output(LOG_CAT_ERROR, LOG_AREA_DFT, 
                        "Creation of thread # %d failed\n", i);
              if (i==0)
                throw "No worker threads could be started";
              else 
                break;
            }
        }
        *xcEnergy = 0;
        electrons = 0;
        while (--i >= 0) {
            pthread_join(pids[i], NULL);
            *xcEnergy += data[i].xcEnergy;
            electrons += data[i].el;
        }
    }
    if(nThreads>1) {
        do_output(LOG_CAT_INFO, LOG_AREA_DFT, 
                  "Electrons: %11.7f %7.1g: KS Energy %f (mt)", 
                  (double)electrons,
                  (double)((electrons-nElectrons)/nElectrons), 
                  (double)*xcEnergy);
        tm.print(LOG_AREA_DFT, __func__);
    }
    return electrons;
}

real
getXC_seq(const BasisInfoStruct& bis, const IntegralInfo& integralInfo,
          const Molecule& mol, const Dft::GridParams& gss, int nElectrons,
          const symmMatrix& dens, symmMatrix& xcm, real* xcEnergy,
          std::vector<int> const & permutationHML)
{
  XCEvaluatorRestricted xcEvr(bis, integralInfo, mol, gss,
                              permutationHML, dens);
  printf("%p\n", (void*)(&xcEvr));
  return xcEvr.getXC(nElectrons, xcm, xcEnergy, 1);
}

/* =================================================================== */
/* Unrestricted sparse code. */
class XCEvaluatorUnrestricted : public XCEvaluator {
  SparseMatrix *dMat[2];
public:
  XCEvaluatorUnrestricted(const BasisInfoStruct& bisOrig_,
                          const IntegralInfo& integralInfo_,
                          const Molecule& mol_,
                          const Dft::GridParams& gss_,
                          std::vector<int> const & permutationHML_,
                          const symmMatrix& densA,
                          const symmMatrix& densB)
    : XCEvaluator(bisOrig_, integralInfo_, mol_, gss_, permutationHML_, densA)
  {
    dMat[0] = new SparseMatrix(*pattern, densA, aoMap, permutationHML);
    dMat[1] = new SparseMatrix(*pattern, densB, aoMap, permutationHML);
  }
  ~XCEvaluatorUnrestricted()
  {
    delete dMat[0];
    delete dMat[1];
  }

  real getXC(int nElectrons, symmMatrix& xcA, symmMatrix& xcB,
             real* xcEnergy, int nThreads) const;
};

real
XCEvaluatorUnrestricted::getXC(int nElectrons,
                               symmMatrix& xca, symmMatrix& xcb,
                               real* xcEnergy, int nThreads) const
{
    real electrons;

    Util::TimeMeter tm;
    bool isGGA = selected_func->is_gga();

    *xcEnergy = 0.0;
    sync_threads(false, nThreads);
    SparseMatrix mata(*pattern), matb(*pattern);
    struct UksData<Dft::SparseMatrix> ds(&mata, &matb, DFT_MAX_BLLEN);

    void (*cblda)(DftIntegratorBl* grid, real * restrict tmp, 
                  int bllen, int blstart, int blend,
                  UksData<SparseMatrix>* data)
      = xcCallbackLdaU<SparseMatrix,XCDistributorLda<SparseMatrix> >;
    void (*cbgga)(DftIntegratorBl* grid, real * restrict tmp, 
                  int bllen, int blstart, int blend,
                  UksData<SparseMatrix>* data)
      = xcCallbackGgaU<SparseMatrix,XCDistributorGgaU<SparseMatrix> >;

    electrons = integrate(2, dMat, *bisPermuted, mol, gss,
                          nThreads,
                          DftBlockCallback(isGGA ? cbgga : cblda),
                          &ds);
    
    pthread_mutex_lock(&dft_prop_mutex);
    *xcEnergy +=ds.energy;
    pthread_mutex_unlock(&dft_prop_mutex);
    mata.addSymmetrizedTo(xca, aoMap, permutationHML);
    matb.addSymmetrizedTo(xcb, aoMap, permutationHML);
    pthread_mutex_unlock(&dft_prop_mutex);

    if(nThreads <= 1) {
      do_output(LOG_CAT_INFO, LOG_AREA_DFT,
                "Electrons: %11.7f %7.1g:  U-xc energy %f (serial)",
                (double)electrons,
                (double)((electrons-nElectrons)/nElectrons), 
                (double)ds.energy);
      tm.print(LOG_AREA_DFT, __func__);
    }
    return electrons;
}

/* multithreaded interface TBW... */

real
getUXC_seq(const BasisInfoStruct& bis, const IntegralInfo& integralInfo,
          const Molecule& mol, const Dft::GridParams& gss, int nElectrons,
          const symmMatrix& densA, const symmMatrix& densB,
          symmMatrix& xcA, symmMatrix& xcB, real* xcEnergy,
          std::vector<int> const & permutationHML)
{
    XCEvaluatorUnrestricted xcEvaluator(bis, integralInfo, mol, gss,
                                        permutationHML, densA, densB);
    return xcEvaluator.getXC(nElectrons, xcA, xcB, xcEnergy, 1);
}


END_NAMESPACE(Dft)
