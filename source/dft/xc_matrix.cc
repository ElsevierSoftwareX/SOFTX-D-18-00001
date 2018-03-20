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
/** @file xc_matrix.cc The XC matrix evaluator.
   (c) Pawel Salek, pawsa@theochem.kth.se.
   2002.04.05

   This module evaluates DFT contribution KS matrix.
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>

#define WITH_PTHREAD 1
#if defined(WITH_PTHREAD)
#include <pthread.h>
static pthread_mutex_t dft_prop_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

#define __CVERSION__
#include "aos.h"
#include "integrator.h"
#include "functionals.h"
#include "dft_common.h"

#include "mat_gblas.h"
#include "output.h"
#include "utilities.h"
#include "matrix_utilities.h"
#include "grid_matrix.h"
#include "xc_evaluators.h"

/* restrict hints should not be necessary... */
#if !defined(restrict)
#define restrict
#endif

const static int KOHNSH_DEBUG = 0;

void lrao2mo_(const real* cmo, const int *ksymop, 
              const real*res, real* fmat, real* work, int*lw);

#if defined(VAR_MPI)
#include <mpi.h>
#define MASTER_NO 0
#endif
#if 0 && defined(VAR_MPI)
#include <mpi.h>
#define MASTER_NO 0

/* dft_kohn_sham_slave:
   this is a slave driver. It's task is to allocate memory needed by
   the main property evaluator (dft_kohn_sham in this case) and call it.
*/
void
dft_kohn_sham_slave(real* work, int* lwork, const int* iprint)
{
    real* dmat = malloc(inforb_.n2basx*sizeof(real));
    real* ksm  = calloc(inforb_.n2basx,sizeof(real));
    int iprfck = 0;
    dft_kohn_sham_(dmat, ksm, work, lwork, &iprfck);
    free(dmat);
    free(ksm);
}

static __inline__ void
dft_kohn_sham_sync_slaves(real* dmat)
{
    MPI_Bcast(dmat, inforb_.n2basx,MPI_DOUBLE,
	      MASTER_NO, MPI_COMM_WORLD);
}

static __inline__ void
dft_kohn_sham_collect_info(real*ksm, real* energy, real* work)
{
    real tmp = *energy;
    dcopy_(&inforb_.n2basx, ksm,&ONEI, work, &ONEI);
    MPI_Reduce(work, ksm, inforb_.n2basx, MPI_DOUBLE, MPI_SUM, 
	       MASTER_NO, MPI_COMM_WORLD);
    MPI_Reduce(&tmp, energy, 1, MPI_DOUBLE, MPI_SUM, 
	       MASTER_NO, MPI_COMM_WORLD);
}

#else /* VAR_MPI */
#define dft_kohn_sham_sync_slaves(dmat)
#define dft_kohn_sham_collect_info(myksm, ksm, energy)
#endif /* VAR_MPI */

/* =================================================================== */
/*                    BLOCKED PROPERTY EVALUATORS                      */
/* =================================================================== */

struct XCDistributorLdaBlas {
  static void distribute(DftIntegratorBl *grid,
                         int bllen, int blstart, int blend,
                         real * restrict tmp, const real *restrict dR,
                         Dft::FullMatrix& excmat);
};

void
XCDistributorLdaBlas::distribute(DftIntegratorBl *grid,
                                 int bllen, int blstart, int blend,
                                 real * restrict tmp, const real *restrict dR,
                                 Dft::FullMatrix& mat)
{
    int isym, jbl, j, ibl, k;
    const real * const aos = grid->atv;
    real * restrict excmat = mat.mat;

    for(isym=0; isym<grid->nsym; isym++) {
        int (*restrict blocks)[2] = BASBLOCK(grid,isym);
        int bl_cnt = grid->bas_bl_cnt[isym];

        for(jbl=0; jbl<bl_cnt; jbl++)
            for(j=blocks[jbl][0]; j<blocks[jbl][1]; j++) { 
                int joff = j*bllen;
                for(k=blstart; k<blend; k++)
                    tmp[k+joff] = aos[k+joff]*dR[k];
            }

        for(jbl=0; jbl<bl_cnt; jbl++) {
	    static const ergo_real HALF = 0.5;
	    static const ergo_real QUARTER = 0.25;
	    for(ibl=0; ibl<jbl; ibl++) {
		int cRows = blocks[ibl][1]-blocks[ibl][0];
		int cCols = blocks[jbl][1]-blocks[jbl][0];
		static const ergo_real ONER = 1.0;
		mat::gemm("T","N", &cRows, &cCols, &bllen, &HALF, 
			  aos+blocks[ibl][0]*bllen, &bllen,
			  tmp+blocks[jbl][0]*bllen, &bllen, &ONER,
			  excmat+blocks[jbl][0]*grid->nbast+blocks[ibl][0],
			  &grid->nbast);
	    }
	    /* This will double-count diagonal elements, need to
	       correct for it later. Or maybe not? */
	    int cRows = blocks[jbl][1]-blocks[jbl][0];
	    int cCols = blocks[jbl][1]-blocks[jbl][0];
	    mat::gemm("T","N", &cRows, &cCols, &bllen, &QUARTER, 
		      aos+blocks[jbl][0]*bllen, &bllen,
		      tmp+blocks[jbl][0]*bllen, &bllen, &ONER,
		      excmat+blocks[jbl][0]*grid->nbast+blocks[jbl][0],
		      &grid->nbast);
	}
    }
}

struct XCDistributorGgaBlas {
  static void distribute(DftIntegratorBl *grid,
                         int bllen, int blstart, int blend,
                         real * restrict tmp,
                         const real *dR, const real *dZ,
                         Dft::FullMatrix& mat);
};

void
XCDistributorGgaBlas::distribute(DftIntegratorBl *grid,
                                 int bllen, int blstart, int blend,
                                 real * restrict tmp,
                                 const real * dR, const real * dZ,
                                 Dft::FullMatrix& mat)
{
    int isym, jbl, j, ibl, k;
    const real * restrict aox = grid->atv+bllen*grid->nbast;
    const real * restrict aoy = grid->atv+bllen*grid->nbast*2;
    const real * restrict aoz = grid->atv+bllen*grid->nbast*3;
    const real * restrict aos = grid->atv;
    real * restrict excmat = mat.mat;

    for(isym=0; isym<grid->nsym; isym++) {
        int (*restrict blocks)[2] = BASBLOCK(grid,isym);
        int nblocks = grid->bas_bl_cnt[isym];
        for(jbl=0; jbl<nblocks; jbl++)
            for(j=blocks[jbl][0]; j<blocks[jbl][1]; j++) { 
                int joff = j*bllen;
                for(k=0; k<bllen; k++)
                    tmp[k+joff] = 
                        dR[k]* aos[k+j*bllen] +
                        dZ[k]*(aox[k+j*bllen]*grid->g.rad.a[k][0]+
                               aoy[k+j*bllen]*grid->g.rad.a[k][1]+
                               aoz[k+j*bllen]*grid->g.rad.a[k][2]);
        }
        
        for(jbl=0; jbl<nblocks; jbl++) {
	    for(ibl=0; ibl<nblocks; ibl++) {
		int cRows = blocks[ibl][1]-blocks[ibl][0];
		int cCols = blocks[jbl][1]-blocks[jbl][0];
		static const ergo_real ONER = 1.0;
		mat::gemm("T","N", &cRows, &cCols, &bllen, &ONER, 
			  aos+blocks[ibl][0]*bllen, &bllen,
			  tmp+blocks[jbl][0]*bllen, &bllen, &ONER,
			  excmat+blocks[jbl][0]*grid->nbast+blocks[ibl][0],
			  &grid->nbast);
	    }
	}
    }
}


/* =================================================================== */
/*                 blocked density and KS evaluation                   */
/* =================================================================== */


#if 0
static void
printmat(int n, const ergo_real *m, const char *name)
{
  int i, j;
  printf("Printing matrix %s\n", name);
  for(i=0; i<n; i++) {
    for(j=0; j<n; j++)
      printf("%10.5f", (double)m[i + j*n]);
    puts("");
  }
}
#endif


/** computes Fock matrix ksm corresponding to given density matrix
   dmat.  fast version - uses memory bandwidth-efficient algorithm.
*/
EXTERN_C real
dft_get_xc(int nElectrons, const real* dmat, const BasisInfoStruct& bis,
           const Molecule& mol, const Dft::GridParams& gss,
	   real* ksm, real* edfty, int nThreads)
{
    int nbast2, i, j;
    real electrons;
    int nbast = bis.noOfBasisFuncs;
    Util::TimeMeter tm;
    bool isGGA = selected_func->is_gga();
    Dft::FullMatrix res(nbast);
    Dft::FullMatrix density(dmat, nbast);
    KsData<Dft::FullMatrix> ds(&res, DFT_MAX_BLLEN);
    const Dft::FullMatrix *densPtr = &density;
    nbast2   = nbast*nbast;

#if USE_BLAS_IN_XC == 1
    void (*cblda)(DftIntegratorBl* grid, real * restrict tmp, 
                  int bllen, int blstart, int blend,
                  KsData<Dft::FullMatrix>* data) =
      &xcCallbackLdaR<Dft::FullMatrix,XCDistributorLdaBlas>;
    void (*cbgga)(DftIntegratorBl* grid, real * restrict tmp, 
                  int bllen, int blstart, int blend,
                  KsData<Dft::FullMatrix>* data) =
      &xcCallbackGgaR<Dft::FullMatrix,XCDistributorGgaBlas>;
#else
    void (*cblda)(DftIntegratorBl* grid, real * restrict tmp, 
              int bllen, int blstart, int blend,
                  KsData<Dft::FullMatrix>* data) =
      &xcCallbackLdaR<Dft::FullMatrix,XCDistributorLda<Dft::FullMatrix> >;
    void (*cbgga)(DftIntegratorBl* grid, real * restrict tmp, 
                  int bllen, int blstart, int blend,
                  KsData<Dft::FullMatrix>* data) =
      &xcCallbackGgaR<Dft::FullMatrix,XCDistributorGga<Dft::FullMatrix> >;
#endif
    electrons = Dft::integrate(1, &densPtr, bis, mol, gss, nThreads,
                               (DftBlockCallback)
                               (isGGA ? cbgga : cblda),
                               &ds);

    for(i=0; i<nbast; i++) {
	int ioff = i*nbast;
	for(j=0; j<i; j++) {
	    int joff = j*nbast;
	    real averag = (res.mat[i+joff] + res.mat[j+ioff]);
	    res.mat[i+joff] = res.mat[j+ioff] = averag;
	}
#if (USE_BLAS_IN_XC  == 1)
        res.mat[i+i*nbast] *= 2.0;
#endif
    }

    pthread_mutex_lock(&dft_prop_mutex);
    *edfty +=ds.energy;
    mat::axpy(&nbast2, &ONER, res.mat, &ONEI, ksm, &ONEI);
    pthread_mutex_unlock(&dft_prop_mutex);

    if(nThreads<=1) {
        if(KOHNSH_DEBUG) {
            output_matrix(nbast, res.mat, "kohn sham matrix");
        }
        int nElectrons = mol.getNumberOfElectrons();
        do_output(LOG_CAT_INFO, LOG_AREA_DFT, 
		  "Electrons: %11.7f %7.1g: xc energy %f (serial)", 
                  (double)electrons,
                  (double)((electrons-nElectrons)/nElectrons), 
                  (double)ds.energy);
	tm.print(LOG_AREA_DFT, __func__);
    }

    return electrons;
}

/* multithreaded interface... */
struct xc_data {
    const real *dmat;
    const BasisInfoStruct *bis;
    const Molecule *mol;
    const Dft::GridParams *gss;
    real *xc, edfty;
    real el;
    int nElectrons;
    int nThreads;
};
static void*
dft_get_xc_worker(void *data)
{
  static const int XCWORKER_ERROR = 3;
    struct xc_data *d = (struct xc_data*)data;
    try {
      d->el = dft_get_xc(d->nElectrons, d->dmat, *d->bis, *d->mol, *d->gss,
			 d->xc, &d->edfty, d->nThreads);
    } catch(const char *s) {
      do_output(LOG_CAT_ERROR, LOG_AREA_DFT, 
		"dft_get_xc_worker thread caught an exception '%s'", s);
      return (void*)&XCWORKER_ERROR;
    } catch(const std::bad_alloc & e) {
      do_output(LOG_CAT_ERROR, LOG_AREA_DFT, 
		"dft_get_xc_worker thread caught an exception '%s'", e.what());
      return (void*)&XCWORKER_ERROR;
    } catch(const std::runtime_error & e) {
      do_output(LOG_CAT_ERROR, LOG_AREA_DFT, 
		"dft_get_xc_worker thread caught an exception '%s'", e.what());
      return (void*)&XCWORKER_ERROR;
    }  catch(...) {
      do_output(LOG_CAT_ERROR, LOG_AREA_DFT, 
		"dft_get_xc_worker thread caught unexpected exception.");
      return (void*)&XCWORKER_ERROR;
    }
    
    return NULL;
}

/** Computes the XC interaction matrix for given density matrix @param dmat .
    @returns the integrated number of electrons.
    @param nElectrons number of electrons.
    @param bis a structure describing the used basis set.
    @param mol a structure describing the molecule.
    @param gss a structure describing the grid settings.
    @param xc resulting XC matrix.
    @param edfty resulting XC energy.
*/
EXTERN_C real
dft_get_xc_mt(int nElectrons, const real* dmat, const BasisInfoStruct& bis,
              const Molecule& mol, const Dft::GridParams& gss,
	      real *xc, real* edfty)
{

    int i, threads;
    real electrons;
    Util::TimeMeter tm;
    
    threads = dft_get_num_threads();
    std::vector<xc_data> data(threads);
    std::vector<pthread_t> pids(threads);
    if(threads == 1) {
	/* Do not create any threads at all to avoid stack allocation. */
        *edfty = 0.0;
	electrons = dft_get_xc(nElectrons, dmat, bis, mol, gss, xc, edfty, 1);
    } else {
	for(i=0; i<threads; i++) {
	    data[i].nElectrons = nElectrons;
	    data[i].dmat = dmat;
	    data[i].xc   = xc;
	    data[i].edfty = 0.0;
	    data[i].bis = &bis;
	    data[i].mol = &mol;
	    data[i].gss = &gss;
            data[i].nThreads = threads;
	    if (pthread_create(&pids[i], NULL, dft_get_xc_worker, &data[i])) {
	      do_output(LOG_CAT_ERROR, LOG_AREA_DFT, 
			"Creation of thread # %d failed\n", i);
	      if (i==0)
		throw "No worker threads could be started";
	      else 
		break;
	    }
	}
	*edfty = 0;
	electrons = 0;
	while ( --i >= 0) {
	    pthread_join(pids[i], NULL);
	    *edfty += data[i].edfty;
	    electrons += data[i].el;
	}

        int nElectrons = mol.getNumberOfElectrons();
        do_output(LOG_CAT_INFO, LOG_AREA_DFT, 
                  "Electrons: %11.7f %7.1g: xc energy %f (mt)", 
                  (double)electrons,
                  (double)((electrons-nElectrons)/nElectrons), 
                  (double)*edfty);
        tm.print(LOG_AREA_DFT, __func__);
    }
    return electrons;
}

/* ===================================================================
   Blocked, unrestricted code
   =================================================================== */
struct uks_data {
  Dft::FullMatrix *exca, *excb;
  real* dRa, *dRb;
  real* dZa, *dZb, *dZab;
  real energy;
};

EXTERN_C real
dft_get_uxc(int nElectrons, const real* dmata, const real *dmatb,
            const BasisInfoStruct& bis, const Molecule& mol,
             const Dft::GridParams& gss,
	    real* xca, real *xcb, real* edfty, int nThreads)
{
    int nbast = bis.noOfBasisFuncs;
    int nbast2, i, j, imat;
    real electrons;
    const Dft::FullMatrix *dmat[2];
    Util::TimeMeter tm;
    bool isGGA = selected_func->is_gga();
    Dft::FullMatrix mata(nbast), matb(nbast);
    Dft::FullMatrix densa(dmata, nbast);
    Dft::FullMatrix densb(dmatb, nbast);
    dmat[0] = &densa; dmat[1] = &densb;
    nbast2   = nbast*nbast;
    UksData<Dft::FullMatrix> ds(&mata, &matb, DFT_MAX_BLLEN);

    void (*cblda)(DftIntegratorBl* grid, real * restrict tmp, 
                  int bllen, int blstart, int blend,
                  UksData<Dft::FullMatrix>* data)
#if USE_BLAS_IN_XC == 1
      = xcCallbackLdaU<Dft::FullMatrix,XCDistributorLdaBlas >;
#else
      = xcCallbackLdaU<Dft::FullMatrix,XCDistributorLda<Dft::FullMatrix> >;
#endif
    void (*cbgga)(DftIntegratorBl* grid, real * restrict tmp, 
                  int bllen, int blstart, int blend,
                  UksData<Dft::FullMatrix>* data)
      = xcCallbackGgaU<Dft::FullMatrix,XCDistributorGgaU<Dft::FullMatrix> >;

    electrons = Dft::integrate(2, dmat, bis, mol, gss, nThreads,
                               (DftBlockCallback)
                               (isGGA ? cbgga : cblda),
                               &ds);

    for(imat=0; imat<2; imat++) {
        real * e = imat ? mata.mat : matb.mat;
        for(i=0; i<nbast; i++) {
            int ioff = i*nbast;
            for(j=0; j<i; j++) {
                int joff = j*nbast;
                real averag = (e[i+joff] + e[j+ioff]);
                e[i+joff] = e[j+ioff] = averag;
            }
#if (USE_BLAS_IN_XC  == 1)
            if (!isGGA) e[i+i*nbast] *= 2.0;
#endif
        }
    }

    pthread_mutex_lock(&dft_prop_mutex);
    *edfty=ds.energy;
    mat::axpy(&nbast2, &ONER, mata.mat, &ONEI, xca, &ONEI);
    mat::axpy(&nbast2, &ONER, matb.mat, &ONEI, xcb, &ONEI);
    pthread_mutex_unlock(&dft_prop_mutex);
    if(KOHNSH_DEBUG) {
        output_matrix(nbast, mata.mat, "Unrestricted xc_alpha matrix");
        output_matrix(nbast, matb.mat, "Unrestricted xc_alpha matrix");
    }

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

/* multithreaded interface... */
struct uxc_data {
    const real* dmata, *dmatb;
    const BasisInfoStruct *bis;
    const Molecule *mol;
    const Dft::GridParams *gss;
    real* xca,   *xcb;
    real edfty, el;
    int nElectrons;
    int nThreads;
};
static void*
dft_get_uxc_worker(void *data)
{
    struct uxc_data *d = (struct uxc_data*)data;
    d->el = dft_get_uxc(d->nElectrons, d->dmata, d->dmatb, *d->bis, *d->mol,
			*d->gss, d->xca, d->xcb, &d->edfty, d->nThreads);
    return NULL;
}

EXTERN_C real
dft_get_uxc_mt(int nElectrons, const real* dmata, const real *dmatb,
               const BasisInfoStruct& bis, const Molecule& mol,
	       const Dft::GridParams& gss,
               real* xca,   real *xcb, real* edfty)
{

    int i, threads;
    real electrons = 0;

    Util::TimeMeter tm;

    threads = dft_get_num_threads();
    std::vector<uxc_data> data(threads);
    std::vector<pthread_t> pids(threads);

    *edfty = 0.0;
    if(threads == 1) {
	/* Do not create any threads at all to avoid stack allocation. */
      electrons = dft_get_uxc(nElectrons, dmata, dmatb, bis, mol,
                              gss, xca, xcb, edfty, threads);
    } else {
      for(i=0; i<threads; i++) {
        data[i].nElectrons = nElectrons;
        data[i].dmata = dmata;
        data[i].dmatb = dmatb;
        data[i].xca = xca;
        data[i].xcb = xcb;
        data[i].edfty = 0.0;
        data[i].bis = &bis;
        data[i].mol = &mol;
	data[i].gss = &gss;
        data[i].nThreads = threads;
	if (pthread_create(&pids[i], NULL, dft_get_uxc_worker, &data[i])) {
          do_output(LOG_CAT_ERROR, LOG_AREA_DFT, 
                    "Creation of thread # %d failed\n", i);
	  if (i==0)
	    throw "No worker threads could be started";
	  else 
	    break;
	}
      }
      while ( --i >= 0) {
        pthread_join(pids[i], NULL);
        *edfty += data[i].edfty;
        electrons += data[i].el;
      }
      int nElectrons = mol.getNumberOfElectrons();
      do_output(LOG_CAT_INFO, LOG_AREA_DFT, 
                "Electrons: %11.7f %7.1g: u-xc energy %f (mt)", 
                (double)electrons,
                (double)((electrons-nElectrons)/nElectrons), 
                (double)*edfty);
      tm.print(LOG_AREA_DFT, __func__);
    }

    return electrons;
}
