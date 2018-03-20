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

/** @file grid_hicu.cc

    @brief Hierarchical Cubature (HiCu) grid generation.

    @author: Elias Rudberg <em>responsible</em>
*/

/* The multi-scale cartesian cubature grid generator. The original
 * implementation described in M. Challacombe, JCP 113(22),
 * p.10037. This one is modified wrt to the paper.
 *
 *  Elias Rudberg, 2004-04
 *
 *  Modified by Elias in Feb-Apr 2010.
*/

#define __CVERSION__

#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <errno.h>
#include <memory.h>
#include <time.h>
#include <pthread.h>
#include <stdexcept>

#include "grid_hicu.h"
#include "basisinfo.h"
#include "integrals_general.h"
#include "cubature_rules.h"
#include "utilities.h"
#include "pi.h"
#include "box_system.h"
#include "integrator.h"
#include "functionals.h"
#include "aos.h"
#include "dft_common.h"
#include "rho-mat.h"
#include "units.h"


static const int CUBATURE_RULE = 3;
static const int CUBATURE_RULE_2 = 6;


const real COORD_DIFF_FOR_SAMEPOINT_CRITERION = 1.0e-11; // FIXME! hard-coded constant here.
const real DISTR_PRODUCT_THRESHOLD = 1e-12; // FIXME! hard-coded constant here.
const real RELATIVE_DENSITY_VARIATION_LIMIT = 0.5; // FIXME! hard-coded constant here.

/* Cutoff value used to decide which gaussian products should be
   ignored. A product with a lower coefficient than this will be
   thrown away. */
const real DISTR_COEFF_CUTOFF_VALUE = 1e-12; // FIXME! hard-coded constant here.

/* FIXME: Descibe the meaning of this constant. */
const real TARGET_RHO_ERROR_FACTOR = 1e-4;


/* flag for test integration. If turned on, the grid file is
reopened after it has been created, and the density is integrated
using the newly created grid. Just to check that it gives the 
same result as reported by the grid generation and by 
the dft integrator. */
//static int    global_doTestIntegration = 0; 

/* Output level. 0 means minimum output, 1 means little output,
   2 means a lot of output. */
//static int    global_outputLevel       = 1;

/* Threshold value for distributions. A gaussian is ignored in areas
where its value is below this threshold. A low value is 
computationally expensive. */
//static real   global_targetRhoError = 1.0e-11;

/* Main threshold error value for grid generation. The difference of
analytical and numerical integrals is forced below this value.
This is the most important parameter, and probably the only one
that a typical user should worry about. */
//static real   global_maxerror       = 1.0e-7;




pthread_mutex_t global_main_hicu_mutex = PTHREAD_MUTEX_INITIALIZER;




#define USE_EXP_STD
#define USE_ERF_STD
#define DO_EXTRA_ERROR_CHECKING
#define FILE_BATCH_N 1000000
#define MAX_NO_OF_POINTS_PER_BATCH 100
#define MAX_NO_OF_SHLBLOCKS 44444
#define EXPONENT_DIFF_LIMIT 1e-22
#define DISTR_CENTER_DIST_LIMIT 1e-22
#define N_BATCH_JOBS 22
#define MAX_NO_OF_POINTS_PER_WRITE 50000


#define HICU_SPARSE_MATRIX_ACCESS_ROUTINE at
//#define HICU_SPARSE_MATRIX_ACCESS_ROUTINE at_safe

const int HICU_GRID_PLOT_RESOLUTION = 50;

/*//////////////////////////////////////////////////////////////////////// */
/*/////////////////  typedef section  //////////////////////////////////// */
/*//////////////////////////////////////////////////////////////////////// */


typedef struct {
  ShellSpecStruct s;
  real extent;
} ShellSpecStructWithExtent;


typedef struct
{
  int noOfShells;
  ShellSpecStructWithExtent* shellList;
  int nbast;
  const Dft::Matrix* dmat;
  BasisFuncStruct* basisFuncList;
  int noOfDistributions;
  DistributionSpecStruct* distrList;
} DensitySpecStruct;



struct rhoTreeNode_{
  BoxStruct box;
    struct rhoTreeNode_* child1; /* NULL for leaf node */
    struct rhoTreeNode_* child2; /* NULL for leaf node */
    int distrIndex;      /* -1 for non-leaf node */
};
typedef struct rhoTreeNode_ rhoTreeNode;

struct GridGenerationParamsStruct {
  real maxerrorPerBox;
  real targetRhoError;
  bool doDoubleChecking;
  bool compareToRefined;
  bool useEnergyCriterion;
  bool useEnergyCriterionOnly;
  bool useErrorPerVolume;
  bool doVariationChecking;
  GridGenerationParamsStruct() :
    maxerrorPerBox(0), /* Must be set later. */
    targetRhoError(0), /* Must be set later. */
    doDoubleChecking(false),
    compareToRefined(false),
    useEnergyCriterion(false),
    useEnergyCriterionOnly(false),
    useErrorPerVolume(false),
    doVariationChecking(false)
  { /* Constructor body. Do nothing here. */ }
};

struct compute_grid_for_box_params_struct 
{
  const BasisInfoStruct& bis;
  DensitySpecStruct density;
  int noOfNonzeroBasisFuncs;
  int* nonZeroBasisFuncIndexList;
  int noOfNonzeroShells;
  int* nonZeroShellsIndexList;
  std::vector<real> localFullDensityMatrix;
  GridGenerationParamsStruct gridGenerationParams;
  int nShlblocks;
  int (*listShlblocks_otherformat)[2];
  DftIntegratorBl* dftIntegrator;
  real* dmagao;
  explicit compute_grid_for_box_params_struct(const BasisInfoStruct& bis_) : bis(bis_) { }
};

typedef std::vector< std::vector< std::vector<int> > > tripleVectorOfInt;

struct ComputeGridResultValuesStruct {
  real totalIntegralResultNumerical;
  real totalIntegralResultAnalytical;
  real totalIntegralResultEnergy;
  real estimatedIntegralErrorEnergy;
  real estimatedIntegralErrorDensity;
  ComputeGridResultValuesStruct() :
    totalIntegralResultNumerical(0),
    totalIntegralResultAnalytical(0),
    totalIntegralResultEnergy(0),
    estimatedIntegralErrorEnergy(0),
    estimatedIntegralErrorDensity(0)
  { /* Constructor. Do nothing here. */ }
};


struct compute_grid_thread_func_struct
{
  const BasisInfoStruct& bis;
  DensitySpecStruct* density;
  rhoTreeNode* rhoTreeRootNode;
  rhoTreeNode* rhoTreeRootNodeShells;
  GridGenerationParamsStruct gridGenerationParams;
  FILE* gridFile;
  BoxStruct* startBox;
  int Nx;
  int Ny;
  int Nz;
  int maxNoOfRelevantDistrsPerBox;
  pthread_mutex_t* fileMutex;
  pthread_mutex_t* jobMutex;
  pthread_t thread;
  int* currJobNumber;
  int noOfPoints;         /* OUTPUT */
  int noOfWrittenBatches; /* OUTPUT */
  ComputeGridResultValuesStruct resultValues;  /* OUTPUT */
  int threadNo;
  int resultCode;
  bool generateSparsePatternOnly;
  Dft::SparsePattern* sparsePattern;
  tripleVectorOfInt* counterArrForPlot;
  explicit compute_grid_thread_func_struct(const BasisInfoStruct& bis_) : bis(bis_) { }
};


/*//////////////////////////////////////////////////////////////////////// */
/*/////////////////  end of typedef section  ///////////////////////////// */
/*//////////////////////////////////////////////////////////////////////// */


/* Solid harmonics based on the table 6.3 of Molecular
 * Electronic-Structure Theory by Helgaker, Jørgensen and Olsen. */

#define solid_harmonic_s_0(x, y, z, x2, y2, z2, r2) 1

/* Elias note: changed order here from 0 1 2 to 2 0 1 which seemed to help. */
#define solid_harmonic_p_2(x, y, z, x2, y2, z2, r2) x
#define solid_harmonic_p_0(x, y, z, x2, y2, z2, r2) y
#define solid_harmonic_p_1(x, y, z, x2, y2, z2, r2) z

#define solid_harmonic_d_0(x, y, z, x2, y2, z2, r2) (x * y)
#define solid_harmonic_d_1(x, y, z, x2, y2, z2, r2) (y * z)
#define solid_harmonic_d_2(x, y, z, x2, y2, z2, r2) ((2 * z2 - x2 - y2) / (2 * template_blas_sqrt((ergo_real)3)))
#define solid_harmonic_d_3(x, y, z, x2, y2, z2, r2) (x * z)
#define solid_harmonic_d_4(x, y, z, x2, y2, z2, r2) (0.5 * (x2 - y2))

#define solid_harmonic_f_0(x, y, z, x2, y2, z2, r2) ((0.5 * template_blas_sqrt(2.5) * (3 * x2 - y2) * y) / template_blas_sqrt((ergo_real)15))
#define solid_harmonic_f_1(x, y, z, x2, y2, z2, r2) (x * y * z)
#define solid_harmonic_f_2(x, y, z, x2, y2, z2, r2) (0.5 * template_blas_sqrt((ergo_real)1.5) * (5 * z2 - r2) * y / template_blas_sqrt((ergo_real)15))
#define solid_harmonic_f_3(x, y, z, x2, y2, z2, r2) (0.5 * (5 * z2 - 3 * r2) * z / template_blas_sqrt((ergo_real)15))
#define solid_harmonic_f_4(x, y, z, x2, y2, z2, r2) (0.5 * template_blas_sqrt((ergo_real)1.5) * (5 * z2 - r2) * x / template_blas_sqrt((ergo_real)15))
#define solid_harmonic_f_5(x, y, z, x2, y2, z2, r2) (0.5 * (x2 - y2) * z)
#define solid_harmonic_f_6(x, y, z, x2, y2, z2, r2) (0.5 * template_blas_sqrt((ergo_real)2.5) * (x2 - 3 * y2) * x / template_blas_sqrt((ergo_real)15))




static void 
print_box(BoxStruct* box) {
  do_output(LOG_CAT_INFO, LOG_AREA_DFT, "print_box:");
  for(int i = 0; i < NO_OF_DIMENSIONS; i++) {
    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "min = %.11f   max = %.11f", 
              (double)box->min[i], (double)box->max[i]);
  } /* END FOR i */
}


static void
get_distribution_box(BoxStruct* box, 
                     DistributionSpecStruct* distr, 
                     real targetRhoError)
{
  real targetError = targetRhoError;
  real arg = distr->coeff / targetError;
  if(arg < 0) arg *= -1;
  if(arg < 1e-30)
    throw std::runtime_error("error in get_distribution_box: (arg < 1e-30).");
  real r1 = template_blas_log(arg);
  if(r1 < 0) r1 *= -1;
  real extent = template_blas_sqrt(r1 / distr->exponent);
  for(int i = 0; i < NO_OF_DIMENSIONS; i++) {
    box->min[i] = distr->centerCoords[i] - extent;
    box->max[i] = distr->centerCoords[i] + extent;
  } /* END FOR i */
} /* END get_distribution_box */

static void 
get_shell_box(BoxStruct* box, ShellSpecStructWithExtent* shell) {
  for(int i = 0; i < NO_OF_DIMENSIONS; i++) {
    box->min[i] = shell->s.centerCoords[i] - shell->extent;
    box->max[i] = shell->s.centerCoords[i] + shell->extent;
  } /* END FOR i */
} /* END get_shell_box */



#if 0
static real
compute_density_at_point_simple(int n, 
                                DensitySpecStruct* density,
                                const BasisInfoStruct & basisInfo,
                                ergo_real x,
                                ergo_real y,
                                ergo_real z)
{
  const Dft::Matrix* dmat = density->dmat;
  // first evaluate each basis function at the given point
  const int MAX_NO_OF_BASIS_FUNCS = 8888;
  if(n > MAX_NO_OF_BASIS_FUNCS)
    throw std::runtime_error("error: (n > MAX_NO_OF_BASIS_FUNCS)");
  ergo_real basisFuncValues[MAX_NO_OF_BASIS_FUNCS];
  for(int i = 0; i < n; i++) {
    ergo_real dx = x - basisInfo.basisFuncList[i].centerCoords[0];
    ergo_real dy = y - basisInfo.basisFuncList[i].centerCoords[1];
    ergo_real dz = z - basisInfo.basisFuncList[i].centerCoords[2];
    ergo_real r2 = dx*dx + dy*dy + dz*dz;
    int nPrims = basisInfo.basisFuncList[i].noOfSimplePrimitives;
    int primIndex = basisInfo.basisFuncList[i].simplePrimitiveIndex;
    ergo_real sum = 0;
    for(int j = 0; j < nPrims; j++) {
      DistributionSpecStruct* currPrim = &basisInfo.simplePrimitiveList[primIndex + j];
      ergo_real factor = 1;
      for(int k = 0; k < currPrim->monomialInts[0]; k++) factor *= dx;
      for(int k = 0; k < currPrim->monomialInts[1]; k++) factor *= dy;
      for(int k = 0; k < currPrim->monomialInts[2]; k++) factor *= dz;
      sum += currPrim->coeff * factor * template_blas_exp(-currPrim->exponent*r2);
    } // END FOR j
    basisFuncValues[i] = sum;
  } // END FOR i
  ergo_real sum = 0;

  // diagonal part
  for(int i = 0; i < n; i++)
    sum += dmat->at(i, i) * basisFuncValues[i] * basisFuncValues[i];
  // non-diagonal part
  for(int i = 0; i < n; i++)
    for(int j = i+1; j < n; j++)
      sum += 2 * dmat->at(i, j) * basisFuncValues[i] * basisFuncValues[j];

  return sum;
}
#endif



static real 
compute_value_at_point(
                       DensitySpecStruct* density,
                       int noOfNonzeroShells,
                       int* nonZeroShellsIndexList,
                       int noOfNonzeroBasFuncs,
                       int* nonZeroBasFuncsIndexList,
                       const real* localFullDensityMatrix,
                       real (*coor)[3],
                       real* workList)
{
  ShellSpecStruct* currShell;
  int i, j, count;
  real expFactor, result, currivalue;
  real xdiff, ydiff, zdiff;
  real x0, y0, z0;
  real x2, y2, z2, r2;
  int nbast;

  nbast = density->nbast;
  //  const Dft::Matrix* dmat = density->dmat;
  if(noOfNonzeroBasFuncs > nbast)
    throw std::runtime_error("error in compute_value_at_point: "
                             "(noOfNonzeroBasFuncs > nbast).");

  /* compute values of contracted distributions at given point */
  count = 0;
  for(i = 0; i < noOfNonzeroShells; i++)
    {
      currShell = &density->shellList[nonZeroShellsIndexList[i]].s;
      x0 = currShell->centerCoords[0];
      y0 = currShell->centerCoords[1];
      z0 = currShell->centerCoords[2];

      xdiff = coor[0][0] - x0;
      ydiff = coor[0][1] - y0;
      zdiff = coor[0][2] - z0;
      x2 = xdiff * xdiff;
      y2 = ydiff * ydiff;
      z2 = zdiff * zdiff;
      r2 = x2 + y2 + z2;

      /* compute expFactor (this is the same procedure for all shell types) */
      expFactor = 0;
      for(j = 0; j < currShell->noOfContr; j++)
          expFactor += currShell->coeffList[j] * 
            template_blas_exp(-currShell->exponentList[j] * r2);
      /* OK, expFactor computed */

      /* now there will be a different number of entries  */
      /* depending on shell type */
      switch(currShell->shellType)
        {
        case 0:
            /* 's' type shell, 1 function */
          workList[count] = expFactor * 
            solid_harmonic_s_0(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
          break;
        case 1:
            /* 'p' type shell, 3 functions */
          workList[count] = expFactor * 
            solid_harmonic_p_0(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
          workList[count] = expFactor * 
            solid_harmonic_p_1(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
          workList[count] = expFactor * 
            solid_harmonic_p_2(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
          break;
        case 2:
            /* 'd' type shell, 5 functions */
          workList[count] = expFactor * 
            solid_harmonic_d_0(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
          workList[count] = expFactor * 
            solid_harmonic_d_1(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
          workList[count] = expFactor * 
            solid_harmonic_d_2(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
          workList[count] = expFactor * 
            solid_harmonic_d_3(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
          workList[count] = expFactor * 
            solid_harmonic_d_4(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
          break;
        case 3:
            /* 'f' type shell, 7 functions */
          workList[count] = expFactor * 
            solid_harmonic_f_0(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
          workList[count] = expFactor * 
            solid_harmonic_f_1(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
          workList[count] = expFactor * 
            solid_harmonic_f_2(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
          workList[count] = expFactor * 
            solid_harmonic_f_3(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
          workList[count] = expFactor * 
            solid_harmonic_f_4(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
          workList[count] = expFactor * 
            solid_harmonic_f_5(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
          workList[count] = expFactor * 
            solid_harmonic_f_6(xdiff, ydiff, zdiff, x2, y2, z2, r2); count++;
          break;
        default:
          throw std::runtime_error("error in compute_value_at_point: "
                                   "only spdf type shells implemented");
        } /* END SWITCH shellType       */
    } /* END FOR i (for each shell) */
  
  if(count > nbast)
    throw std::runtime_error("error in compute_value_at_point: (count > nbast).");
  
  /* now use density matrix to obtain final result */
  result = 0;
  /* Diagonal part. */
  for(i = 0; i < noOfNonzeroBasFuncs; i++)
    result += localFullDensityMatrix[i*noOfNonzeroBasFuncs+i] * workList[i] * workList[i];
  /* Off-diagonal part. */
  for(i = 0; i < noOfNonzeroBasFuncs; i++) {
    currivalue = workList[i];
    for(j = i+1; j < noOfNonzeroBasFuncs; j++)
      result += 2 * localFullDensityMatrix[i*noOfNonzeroBasFuncs+j] * currivalue * workList[j];
  } /* END FOR i */
  
  return result;
} /* END compute_value_at_point */







static real 
compute_integral_from_points(const BasisInfoStruct& bis,
                             DensitySpecStruct* density,
                             int noOfNonzeroShells,
                             int* nonZeroShellsIndexList,
                             int noOfNonzeroBasFuncs,
                             int* nonZeroBasFuncsIndexList,
                             const real* localFullDensityMatrix,
                             int nPoints,
                             real (*coor)[3],
                             real* weight,
                             real* workList,
			     real & minValue,
			     real & maxValue,
			     real & maxAbsValue)
{
  real sum = 0;
  for(int i = 0; i < nPoints; i++) {
    real value = compute_value_at_point(density,
                                        noOfNonzeroShells,
                                        nonZeroShellsIndexList,
                                        noOfNonzeroBasFuncs,
                                        nonZeroBasFuncsIndexList,
                                        localFullDensityMatrix,
                                        &coor[i],
                                        workList);
    if(i == 0) {
      minValue = maxValue = value;
      maxAbsValue = template_blas_fabs(value);
    }
    if(value < minValue)
      minValue = value;
    if(value > maxValue)
      maxValue = value;
    if(template_blas_fabs(value) > maxAbsValue)
      maxAbsValue = template_blas_fabs(value);
#if 0
    if(i == 0) {
      // Verify value by computing in a differrent way.
      real value2 = compute_density_at_point_simple(bis.noOfBasisFuncs, 
                                                    density,
                                                    bis,
                                                    coor[i][0],
                                                    coor[i][1],
                                                    coor[i][2]);
      real tol = 1e-11;
      if(template_blas_fabs(value - value2) > tol)
        throw "Error in compute_integral_from_points: (template_blas_fabs(value - value2) > tol).";
    }
#endif
    sum += value * weight[i];
  } /* END FOR i */
  return sum;
} /* END compute_integral_from_points */

#if 0
static float
hicuErf(float a)
{ return erff(a); }

static double
hicuErf(double a)
{ return erf(a); }

static long double
hicuErf(long double a)
{ return erfl(a); }
#endif

template<class Treal>
Treal hicuErf(Treal a) {
  throw "Error: hicuErf default implementation does not exist.";
}

template<>
float
hicuErf(float a)
{ return erff(a); }

template<>
double
hicuErf(double a)
{ return erf(a); }

// Need to check HAVE_ERFL here, otherwise cannot compile in Cygwin.
#ifdef HAVE_ERFL
template<>
long double
hicuErf(long double a)
{ return erfl(a); }
#endif

#ifdef PRECISION_QUAD_FLT128
template<>
__float128
hicuErf(__float128 a)
{ return erfq(a); }
#endif

static real 
to_power(real x, int n) {
  real result = 1;
  for(int i = 0; i < n; i++)
    result *= x;
  return result;
}


static real 
compute_1d_gaussian_integral_recursive(real a, real b, int n, real alpha)
{
  real result, sqrtalpha, term1, term2;
  real aToPowerNminus1, bToPowerNminus1;
  if(n == 0)
    {
      sqrtalpha = template_blas_sqrt(alpha);
      result = template_blas_sqrt(pi/(4*alpha)) * (hicuErf(sqrtalpha*b)-hicuErf(sqrtalpha*a));
      return result;
    }
  if(n == 1)
    {
      result = -(1 / (2*alpha)) * (template_blas_exp(-alpha*b*b) - template_blas_exp(-alpha*a*a));
      return result;
    }
  if(n < 0)
    throw std::runtime_error("error in 1dintegral: n < 0");
  /* now we know that n >= 2 */
  term1 = (n - 1) * compute_1d_gaussian_integral_recursive(a, b, n-2, alpha);
  aToPowerNminus1 = to_power(a, n-1);
  bToPowerNminus1 = to_power(b, n-1);
  term2  = 
    bToPowerNminus1 * template_blas_exp(-alpha*b*b) - 
    aToPowerNminus1 * template_blas_exp(-alpha*a*a);
  result = (term1 - term2) / (2 * alpha);
  /*  return 0; */
  return result;
} /* END compute_1d_gaussian_integral_recursive */



static real
compute_1d_gaussian_integral(real a, real b, int n, real alpha)
{
  real result, sqrtalpha, term1, term2;
  //  return compute_1d_gaussian_integral_recursive(a, b, n, alpha);
  result = 0;
  switch(n)
    {
    case 0:
      sqrtalpha = template_blas_sqrt(alpha);
      result = template_blas_sqrt(pi/(4*alpha)) * (hicuErf(sqrtalpha*b)-hicuErf(sqrtalpha*a));
      break;
    case 1:
      result = -(1 / (2*alpha)) * (template_blas_exp(-alpha*b*b) - template_blas_exp(-alpha*a*a));
      break;
    case 2:
      sqrtalpha = template_blas_sqrt(alpha);
      term1 = 
        template_blas_sqrt(pi/(16*alpha*alpha*alpha)) * 
        (hicuErf(sqrtalpha*b) - hicuErf(sqrtalpha*a));
      term2 = -(1 / (2 * alpha)) * (b*template_blas_exp(-alpha*b*b) - a*template_blas_exp(-alpha*a*a));
      result = term1 + term2;
      break;
    case 3:
      result = -(1 / (2*alpha*alpha)) * ((1+alpha*b*b)*template_blas_exp(-alpha*b*b) - 
                                         (1+alpha*a*a)*template_blas_exp(-alpha*a*a));
      break;
    default:
      return compute_1d_gaussian_integral_recursive(a, b, n, alpha);
      break;
    } /* END SWITCH n */
  return result;
} /* END compute_1d_gaussian_integral */


static real 
compute_integral_over_box(DistributionSpecStruct* distr, BoxStruct* box)
{
  real result, a, b, alpha;
  int i, n;
  result = distr->coeff;
  alpha = distr->exponent;
  for(i = 0; i < NO_OF_DIMENSIONS; i++)
    {
      n = distr->monomialInts[i];
      a = box->min[i] - distr->centerCoords[i];
      b = box->max[i] - distr->centerCoords[i];
      result *= compute_1d_gaussian_integral(a, b, n, alpha);
    } /* END FOR i */
  return result;
} /* END compute_integral_over_box */


static int 
get_rhotree_indexes_for_box(int* resultList, int resultListMaxCount, const rhoTreeNode* node, const BoxStruct* inputBoxPtr)
{
#define MAX_DEPTH 888
  int n, i, overlap, currDepth;
  const rhoTreeNode* nodeList[MAX_DEPTH];
  int statusList[MAX_DEPTH];
  const rhoTreeNode* currNode;
  BoxStruct box;
  const BoxStruct* currBox;

  memcpy(&box, inputBoxPtr, sizeof(BoxStruct));

  n = 0;
  currDepth = 0;
  nodeList[0] = node;
  statusList[0] = 0;
  while(currDepth >= 0)
    {
      if(statusList[currDepth] == 2)
        currDepth--;
      else
        {

          currNode = nodeList[currDepth];
          currBox = &currNode->box;

          /* check for box overlap */
          overlap = 1;
          for(i = 0; i < NO_OF_DIMENSIONS; i++)
            {
              if(currBox->min[i] > box.max[i])
                overlap = 0;
              if(currBox->max[i] < box.min[i])
                overlap = 0;
            } /* END FOR i */
          if(overlap == 0)
            currDepth--;
          else
            {

              if(statusList[currDepth] == 0)
                {
                  if(currNode->distrIndex >= 0)
                    {
		      if(resultList) {
			assert(n < resultListMaxCount);
			resultList[n] = currNode->distrIndex;
		      }
                      n++;
                      currDepth--;
                    }
                  else
                    {
                      statusList[currDepth] = 1;
                      currDepth++;
                      statusList[currDepth] = 0;
                      nodeList[currDepth] = currNode->child1;
                    }
                } /* END IF status 0 */
              else
                {
                    /* status is 1 */
                  statusList[currDepth] = 2;
                  currDepth++;
                  statusList[currDepth] = 0;
                  nodeList[currDepth] = currNode->child2;             
                } /* END ELSE status 1 */
            }
        }
    } /* END WHILE (currDepth >= 0) */

  return n;
} /* END get_rhotree_indexes_for_box */



static void
callbackGga(DftIntegratorBl* grid, int bllen, real & energy) {
  FunDensProp dp = { 0 };
  assert(grid->ntypso >0);
  for(int k = 0; k < bllen; k++) {
    real weight = grid->weight[grid->curr_point+k];
    dp.grada = 0.5*template_blas_sqrt(grid->g.grad[k][0]*grid->g.grad[k][0]+
			grid->g.grad[k][1]*grid->g.grad[k][1]+
			grid->g.grad[k][2]*grid->g.grad[k][2]);
    dp. rhoa = dp.rhob = 0.5*grid->r.rho[k];
    dp.gradb  = dp.grada;
    dp.gradab = dp.grada*dp.gradb;
    if(dp.rhoa>1e-14) {
      if(dp.grada<1e-35) dp.grada = 1e-35;
      energy += selected_func->func(&dp)*weight;
    }
  }
}

static void
callbackLda(DftIntegratorBl *grid, int bllen, real & energy) {
  FunDensProp dp = { 0 };    
  assert(grid->ntypso >0);
  for(int k = 0; k < bllen; k++) {
    real weight = grid->weight[grid->curr_point+k];
    dp.rhoa = dp. rhob = 0.5*grid->r.rho[k];
    energy += selected_func->func(&dp)*weight;
  }
}


static void
integrate_density_and_energy(const BasisInfoStruct& bis,
			     DensitySpecStruct* density,
			     DftIntegratorBl* integrator,
			     real & electrons,
			     real & energy,
			     int noOfGridPoints,
			     real (*coor)[3],
			     real *weight,
			     real* dmagao) {
  // Initialize integrator.
  for(int kk = 0; kk < noOfGridPoints; kk++) {
    for(int mm = 0; mm < 3; mm++)
      integrator->coor[kk][mm] = coor[kk][mm];
    integrator->weight[kk] = weight[kk];
  }
  int ipnt = 0;
  integrator->curr_point  = ipnt;
  int len = noOfGridPoints;
  int nder = integrator->dogga ? 1 : 0;
  dft_get_orbs(len, integrator->atv, (real(*)[3]) &integrator->coor[ipnt][0],
	       integrator->shl_bl_cnt, (int(*)[2]) &integrator->shlblocks[0][0],
	       nder, bis);
  //  const real** dmatFullPtr = NULL;
  int nbast = bis.noOfBasisFuncs;
  const Dft::Matrix *dens = density->dmat;
  if (dens->isSparse()) {
    if(integrator->dogga)
      getrho_blocked_gga(nbast,
                         *dens->asSparse(), integrator->atv,
                         integrator->bas_bl_cnt,
                         integrator->basblocks, integrator->shl_bl_cnt,
                         &dmagao[0], len, integrator->r.rho,
                         integrator->g.rad.a);
    else
      getrho_blocked_lda(nbast,
                         *dens->asSparse(),
                         integrator->atv,
                         integrator->bas_bl_cnt,
                         integrator->basblocks, integrator->shl_bl_cnt,
                         &dmagao[0], len, integrator->r.rho);
  } else {
    if(integrator->dogga)
      getrho_blocked_gga(nbast,
                         dens->asFull(), integrator->atv,
                         integrator->bas_bl_cnt,
                         integrator->basblocks, integrator->shl_bl_cnt,
                         &dmagao[0], len, integrator->r.rho,
                         integrator->g.rad.a);
    else
      getrho_blocked_lda(nbast,
                         dens->asFull(),
                         integrator->atv,
                         integrator->bas_bl_cnt,
                         integrator->basblocks, integrator->shl_bl_cnt,
                         &dmagao[0], len, integrator->r.rho);
  }

  for(int j=0; j<len; j++)
    electrons += integrator->weight[ipnt+j]*integrator->r.rho[j];
  real energyTmp = 0;
  if(selected_func->is_gga()) 
    callbackGga(integrator, len, energyTmp);
  else
    callbackLda(integrator, len, energyTmp);
  energy = energyTmp;
}



static int 
compute_grid_for_box(compute_grid_for_box_params_struct* params,
                     int maxlen,
                     real (*coor)[3],
                     real *weight,
                     BoxStruct* box,
                     real analyticalIntegralValue,
                     real* workList,
		     ComputeGridResultValuesStruct & resultValues,
		     bool resolutionIsOk)
{
#define MAX_NO_OF_TEST_POINTS 1000

  int Ngrid;
  BoxStruct box1;
  BoxStruct box2;
  int bestcoord, nPoints1, nPoints2;
  real dist, maxdist, halfway;
  real IexactAbs;
  real analyticalIntegralBox1, analyticalIntegralBox2;

  real minDensityValue = 0;
  real maxDensityValue = 0;
  real maxDensityAbsValue = 0;

  int splitBox = 1;
  int noOfGridPoints = 0;

  bool resolutionIsOkForNextLevel = resolutionIsOk;

  // Compute box volume.
  ergo_real boxVolume = 1;
  for(int i = 0; i < NO_OF_DIMENSIONS; i++)
    boxVolume *= (box->max[i] - box->min[i]);
  ergo_real absErrorLimit = params->gridGenerationParams.maxerrorPerBox;
  if(params->gridGenerationParams.useErrorPerVolume) {
    // Modify absErrorLimit accordingly.
    absErrorLimit *= boxVolume;
  }

  real Iapprox = 0; // To be computed below.
  real Iexact = analyticalIntegralValue;
  real xcEnergy = 0; // To be computed below.
  real xcEnergyApproxError = 0; // To be computed below.
  real densityApproxError = 0; // To be computed below.

  if(params->gridGenerationParams.compareToRefined) {
    /* Create several different grids for comparison: first a rough
       grid where the cubature rule is only applied once for the whole
       box, then refined grids where the cubature rule is applied in
       (level)^3 sub-boxes. */
    const int NLEVELSMAX = 4;
    real integralResultList_density[NLEVELSMAX];
    /* We choose the number of levels depending on how close we are to
       the point when the density integral is itself already below
       threshold. */
    int noOfLevels = 3;
    int resultLevelIndex = 0;
    if(Iexact < absErrorLimit)
      noOfLevels = 1;
    else if(Iexact < absErrorLimit*10 && noOfLevels > 2)
      noOfLevels = 2;
    else if(Iexact < absErrorLimit*100 && noOfLevels > 3)
      noOfLevels = 3;
    if(resolutionIsOk && noOfLevels > 2)
      noOfLevels = 2;
    real integralResultList_energy [NLEVELSMAX];
    for(int levelIdx = 0; levelIdx < noOfLevels; levelIdx++) {
      real tmpCoor[MAX_NO_OF_TEST_POINTS][3];
      real tmpWeight[MAX_NO_OF_TEST_POINTS];
      int nBoxesPerDim = levelIdx+1;
      integralResultList_density[levelIdx] = 0;
      integralResultList_energy [levelIdx] = 0;
      int count = 0;
      for(int ix = 0; ix < nBoxesPerDim; ix++)
	for(int iy = 0; iy < nBoxesPerDim; iy++)
	  for(int iz = 0; iz < nBoxesPerDim; iz++) {
	    BoxStruct boxTmp;
	    boxTmp.min[0] = box->min[0] + (ix+0) * (box->max[0] - box->min[0]) / nBoxesPerDim;
	    boxTmp.max[0] = box->min[0] + (ix+1) * (box->max[0] - box->min[0]) / nBoxesPerDim;
	    boxTmp.min[1] = box->min[1] + (iy+0) * (box->max[1] - box->min[1]) / nBoxesPerDim;
	    boxTmp.max[1] = box->min[1] + (iy+1) * (box->max[1] - box->min[1]) / nBoxesPerDim;
	    boxTmp.min[2] = box->min[2] + (iz+0) * (box->max[2] - box->min[2]) / nBoxesPerDim;
	    boxTmp.max[2] = box->min[2] + (iz+1) * (box->max[2] - box->min[2]) / nBoxesPerDim;
	    int nTmp = use_cubature_rule(MAX_NO_OF_TEST_POINTS-count, 
					 &tmpCoor[count], &tmpWeight[count], &boxTmp, CUBATURE_RULE);
	    if(nTmp <= 0)
	      throw std::runtime_error("error in use_cubature_rule.");
	    real electronsTmp = 0;
	    real energyTmp = 0;
	    integrate_density_and_energy(params->bis,
					 &params->density,
					 params->dftIntegrator, 
					 electronsTmp, 
					 energyTmp,
					 nTmp, 
					 &tmpCoor[count], 
					 &tmpWeight[count],
					 params->dmagao);
	    integralResultList_density[levelIdx] += electronsTmp;
	    integralResultList_energy [levelIdx] += energyTmp;
	    count += nTmp;
	  }
      if(levelIdx == resultLevelIndex) {
	if(count > maxlen)
	  throw std::runtime_error("error in compute_grid_for_box: (count > maxlen).");
	for(int i = 0; i < count; i++) {
	  for(int k = 0; k < 3; k++)
	    coor[i][k] = tmpCoor[i][k];
	  weight[i] = tmpWeight[i];
	}
	noOfGridPoints = count;
	Iapprox = integralResultList_density[levelIdx];
	xcEnergy = integralResultList_energy[levelIdx];
      }
    }
    // Compute errors in density integrals by comparing to analytical value.
    real integralErrors_density[NLEVELSMAX];
    for(int levelIdx = 0; levelIdx < noOfLevels; levelIdx++)
      integralErrors_density[levelIdx] = template_blas_fabs(integralResultList_density[levelIdx] - analyticalIntegralValue);
    // Compute errors in energy integrals by comparing to most accurate value.
    real integralErrors_energy[NLEVELSMAX];
    for(int levelIdx = 0; levelIdx < NLEVELSMAX; levelIdx++)
      integralErrors_energy[levelIdx] = 0;
    for(int levelIdx = 0; levelIdx < noOfLevels-1; levelIdx++) {
      integralErrors_energy[levelIdx] = 
	template_blas_fabs(integralResultList_energy[levelIdx] - integralResultList_energy[noOfLevels-1]);
//      printf("integralErrors_energy[levelIdx] = %33.22f\n", integralErrors_energy[levelIdx]);
    }
    // Compute improvement factors;
    real improvementFactorList_density[NLEVELSMAX-1];
    for(int levelIdx = 0; levelIdx < NLEVELSMAX-1; levelIdx++)
      improvementFactorList_density[levelIdx] = 0;
    for(int levelIdx = 0; levelIdx < noOfLevels-1; levelIdx++)
      improvementFactorList_density[levelIdx] = 
	integralErrors_density[levelIdx] / integralErrors_density[levelIdx+1];
    real improvementFactorList_energy [NLEVELSMAX-2];
    for(int levelIdx = 0; levelIdx < NLEVELSMAX-2; levelIdx++)
      improvementFactorList_energy[levelIdx] = 0;
    for(int levelIdx = 0; levelIdx < noOfLevels-2; levelIdx++)
      improvementFactorList_energy[levelIdx] = 
	integralErrors_energy[levelIdx] / integralErrors_energy[levelIdx+1];

    ergo_real expectedErrorFromDensityEvaluations = 
      noOfGridPoints * params->gridGenerationParams.targetRhoError * boxVolume;

#if 0
    printf("resolutionIsOk = %d\n", (int)resolutionIsOk);
    printf("Improvement factors,  density: ");
    for(int levelIdx = 0; levelIdx < noOfLevels-1; levelIdx++)
      printf("%6.2f  ", improvementFactorList_density[levelIdx]);
    printf("   energy: ");
    for(int levelIdx = 0; levelIdx < noOfLevels-2; levelIdx++)
      printf("%6.2f  ", improvementFactorList_energy[levelIdx]);
    printf("\n");
    printf("integralErrors_density[resultLevelIndex]    = %33.22f\n", 
	   integralErrors_density[resultLevelIndex]);
    printf("expectedErrorFromDensityEvaluations         = %33.22f\n", expectedErrorFromDensityEvaluations);
    printf("integralErrors_energy[resultLevelIndex]     = %33.22f\n", 
	   integralErrors_energy[resultLevelIndex]);
    printf("absErrorLimit                               = %33.22f\n", absErrorLimit);
    printf("Iexact                                      = %33.22f\n", Iexact);
#endif

    ergo_real expectedImprovementFactors[5];
    expectedImprovementFactors[0] = 64.00;
    expectedImprovementFactors[1] = 11.39;
    expectedImprovementFactors[2] =  5.62;
    expectedImprovementFactors[3] =  3.81;
    expectedImprovementFactors[4] =  2.99;
    if(noOfLevels > 5)
      throw std::runtime_error("Error: (noOfLevels > 5).");

    densityApproxError = integralErrors_density[resultLevelIndex];
    xcEnergyApproxError = integralErrors_energy[resultLevelIndex];

    // TODO: use some clever splitBox criterion here.
    splitBox = 0;
    if(params->gridGenerationParams.useEnergyCriterionOnly == false) {
      if(densityApproxError > absErrorLimit)
	splitBox = 1;
    }
    if(params->gridGenerationParams.useEnergyCriterion) {
      if(xcEnergyApproxError > absErrorLimit)
	splitBox = 1;
    }
    if(!resolutionIsOk) {
      // Also check that all improvement factors are reasonably near expected values.
      resolutionIsOkForNextLevel = true;
      if(params->gridGenerationParams.useEnergyCriterionOnly == false) {
	for(int levelIdx = 0; levelIdx < noOfLevels-1; levelIdx++) {
	  if(improvementFactorList_density[levelIdx] < expectedImprovementFactors[levelIdx]*0.5)
	    resolutionIsOkForNextLevel = false;
	}
      }
      if(params->gridGenerationParams.useEnergyCriterion) {
	for(int levelIdx = 0; levelIdx < noOfLevels-2; levelIdx++) {
	  if(improvementFactorList_energy[levelIdx] < expectedImprovementFactors[levelIdx]*0.5)
	    resolutionIsOkForNextLevel = false;
	}
      }
      if(!resolutionIsOkForNextLevel)
	splitBox = 1;
    }


    if(params->gridGenerationParams.useEnergyCriterionOnly == false) {
      /* If the integral value is itself below the error limit, we are
	 happy. This happens for example if the integral is completely
	 zero, then it makes no sense to compare different cubature
	 rules. */
      if(Iexact < absErrorLimit)
	splitBox = 0;
      /* If the expected error from the numerical evaluation of the
	 density is comparable to the integral error, there is no point
	 in dividing box further. */
      if(integralErrors_density[resultLevelIndex] < expectedErrorFromDensityEvaluations*2)
	splitBox = 0;
    }

    /* If the computed errors on all levels are well below the
       threshold, we are happy. */
    bool allErrorsWellBelowThreshold = true;
    if(params->gridGenerationParams.useEnergyCriterionOnly == false) {
      for(int levelIdx = 0; levelIdx < noOfLevels; levelIdx++) {
	if(integralErrors_density[levelIdx] > absErrorLimit/100)
	  allErrorsWellBelowThreshold = false;
      }
    }
    if(params->gridGenerationParams.useEnergyCriterion) {
      for(int levelIdx = 0; levelIdx < noOfLevels-1; levelIdx++) {
	if(integralErrors_energy[levelIdx] > absErrorLimit/100)
	  allErrorsWellBelowThreshold = false;
      }
    }
    if(allErrorsWellBelowThreshold)
      splitBox = 0;

    //    printf("splitBox = %d\n\n", splitBox);
      
  } else { // old version

    /* Define Ngrid points inside box, with corresponding weights */
    /* this is where the 'cubature rule' is used */
    Ngrid = use_cubature_rule(maxlen, coor, weight, box, CUBATURE_RULE);
    if(Ngrid <= 0)
      throw std::runtime_error("error in use_cubature_rule.");
    noOfGridPoints = Ngrid;

    Iapprox = compute_integral_from_points(params->bis,
					   &params->density,
					   params->noOfNonzeroShells,
					   params->nonZeroShellsIndexList,
					   params->noOfNonzeroBasisFuncs,
					   params->nonZeroBasisFuncIndexList,
					   &params->localFullDensityMatrix[0],
					   Ngrid,
					   &coor[0],
					   weight,
					   workList,
					   minDensityValue,
					   maxDensityValue,
					   maxDensityAbsValue);
    IexactAbs = Iexact;
    if(IexactAbs < 0) IexactAbs *= -1;

    /* compute absolute error */
    densityApproxError = template_blas_fabs(Iexact - Iapprox);

    /* check if error is too large */
    splitBox = 1;

#if 0
    /* It may happen that the box is now so small that the absErrorLimit
       is so small that it approaches the accuracy with which we can
       compute the density (targetRhoError). In that case there is no
       point in using such a small absErrorLimit; we then set it to a
       value we can handle. */
    const real DENSITY_ACCURACY_COMPARISON_FACTOR = 10.0; // FIXME! hard-coded constant here.
    if(absErrorLimit < params->targetRhoError*DENSITY_ACCURACY_COMPARISON_FACTOR)
      absErrorLimit = params->targetRhoError*DENSITY_ACCURACY_COMPARISON_FACTOR;
#endif

    /* If the integral value is itself below the error limit, we are
       happy. This happens for example if the integral is completely
       zero, then it makes no sense to compare different cubature
       rules. */
    if(Iexact < absErrorLimit)
      splitBox = 0;

    if((splitBox == 1) && (densityApproxError < absErrorLimit))
      {

	if(params->gridGenerationParams.doDoubleChecking) {
	  /* it seems that the error is small enough. */
	  /* however, this could be a coincidence. */
	  /* to check, compare with denser grid */
	  real testCoor[MAX_NO_OF_TEST_POINTS][3];
	  real testWeight[MAX_NO_OF_TEST_POINTS];
	  real testIapprox;
	  int Ngrid2;

	  Ngrid2 = use_cubature_rule(MAX_NO_OF_TEST_POINTS, 
				     testCoor, testWeight, box, CUBATURE_RULE_2);
	  if(Ngrid2 <= 0)
	    throw std::runtime_error("error in use_cubature_rule");

	  real minValueDummy = 0, maxValueDummy = 0, maxAbsValueDummy = 0;
	  testIapprox = 
	    compute_integral_from_points(params->bis,
					 &params->density,
					 params->noOfNonzeroShells,
					 params->nonZeroShellsIndexList,
					 params->noOfNonzeroBasisFuncs,
					 params->nonZeroBasisFuncIndexList,
					 &params->localFullDensityMatrix[0],
					 Ngrid2,
					 &testCoor[0],
					 testWeight,
					 workList,
					 minValueDummy, 
					 maxValueDummy, 
					 maxAbsValueDummy);
	  real testAbsError = template_blas_fabs(Iexact - testIapprox);
	  /* We demand that the denser grid should also work. */
	  if(testAbsError < absErrorLimit)
	    splitBox = 0;
	}
	else
	  splitBox = 0;
      }
    if(splitBox == 0 && Iexact > absErrorLimit*2 && params->gridGenerationParams.doVariationChecking) {
      /* Check that variation of density is not too large. */
      real diff = maxDensityValue - minDensityValue;
      real relativeVariation = diff / maxDensityAbsValue;
      if(relativeVariation > RELATIVE_DENSITY_VARIATION_LIMIT)
	splitBox = 1;
    }
  } // end else old version

  if(splitBox == 1)
    {
      /* error too large, split box into box1 and box2 */
      /* first determine in which coordinate direction to do the split */
      maxdist = 0;
      bestcoord = -1;
      for(int i = 0; i < NO_OF_DIMENSIONS; i++) {
	dist = box->max[i] - box->min[i];
	if(dist > maxdist) {
	  maxdist = dist;
	  bestcoord = i;
	}
      } /* END FOR i */
      if(bestcoord < 0)
	throw std::runtime_error("error in compute_grid_for_box: (bestcoord < 0).");
      /* now create new boxes box1 and box2 */
      for(int i = 0; i < NO_OF_DIMENSIONS; i++)
        {
          if(i == bestcoord)
            {
	      /* direction of split */
              halfway = (box->max[i] + box->min[i]) / 2;
              box1.min[i] = box->min[i];
              box1.max[i] = halfway;
              box2.min[i] = halfway;
              box2.max[i] = box->max[i];
            }
          else
            {
	      /* other direction, simply copy bounds */
              box1.min[i] = box->min[i];
              box1.max[i] = box->max[i];
              box2.min[i] = box->min[i];
              box2.max[i] = box->max[i];
            }
        } /* END FOR i */
      /* now boxes box1 and box2 are now created */
      
      analyticalIntegralBox1 = 0;
      for(int i = 0; i < params->density.noOfDistributions; i++)
        analyticalIntegralBox1 += 
          compute_integral_over_box(&params->density.distrList[i], &box1);

#if 1
      analyticalIntegralBox2 = 
        analyticalIntegralValue - analyticalIntegralBox1;
#else
      analyticalIntegralBox2 = 0;
      for(i = 0; i < params->density.noOfDistributions; i++)
        analyticalIntegralBox2 += 
          compute_integral_over_box(&params->density.distrList[i], &box2);
#endif

      /* create grid points for box1 */
      nPoints1 = compute_grid_for_box(params,
                                      maxlen,
                                      coor,
                                      weight,
                                      &box1,
                                      analyticalIntegralBox1,
                                      workList,
				      resultValues,
				      resolutionIsOkForNextLevel);
      if(nPoints1 < 0)
	throw std::runtime_error("error in compute_grid_for_box: (nPoints1 < 0).");
      /* create grid points for box2 */
      nPoints2 = compute_grid_for_box(params,
                                      maxlen-nPoints1,
                                      &coor[nPoints1],
                                      &weight[nPoints1],
                                      &box2,
                                      analyticalIntegralBox2,
                                      workList,
				      resultValues,
				      resolutionIsOkForNextLevel);
      if(nPoints2 < 0)
	throw std::runtime_error("error in compute_grid_for_box: (nPoints2 < 0).");
      noOfGridPoints = nPoints1 + nPoints2;
    } /* END IF error too large */
  else
    {
        /* error acceptable,  */
        /* the computed grid points for this box are good enough. */
        /* do nothing more, just return the number of points */
      resultValues.totalIntegralResultNumerical += Iapprox;
      resultValues.totalIntegralResultAnalytical += analyticalIntegralValue;
      resultValues.totalIntegralResultEnergy += xcEnergy;
      resultValues.estimatedIntegralErrorEnergy += xcEnergyApproxError;
      resultValues.estimatedIntegralErrorDensity += densityApproxError;
    }

  return noOfGridPoints;
} /* END compute_grid_for_box */



static rhoTreeNode* 
BuildRhoTreeBranch(int noOfDistributionsTot,
                   DistributionSpecStruct* rho_alt_1,
                   ShellSpecStructWithExtent* rho_alt_2,
                   int distrIndexListN, 
                   int* distrIndexList,
                   real targetRhoError)
{
  int n1, n2, bestCoord;
  real currCoord, currDiff, maxDiff, extent1, extent2, testCoord;
  int tempInt;

  if(distrIndexListN < 1) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "Error in BuildRhoTreeBranch: (distrIndexListN < 1), "
	      "distrIndexListN = %i", distrIndexListN);
    return NULL;
  }

  rhoTreeNode* newNode = new rhoTreeNode;

  /* compute bounding box for this node */
  if(rho_alt_1 != NULL)
    get_distribution_box(&newNode->box, &rho_alt_1[distrIndexList[0]], targetRhoError);
  else
    get_shell_box(&newNode->box, &rho_alt_2[distrIndexList[0]]);
  BoxStruct tempBox;
  for(int i = 1; i < distrIndexListN; i++) {
    if(rho_alt_1 != NULL)
      get_distribution_box(&tempBox, &rho_alt_1[distrIndexList[i]], targetRhoError);
    else
      get_shell_box(&tempBox, &rho_alt_2[distrIndexList[i]]);
    for(int j = 0; j < NO_OF_DIMENSIONS; j++) {
      if(tempBox.min[j] < newNode->box.min[j]) 
	newNode->box.min[j] = tempBox.min[j];
      if(tempBox.max[j] > newNode->box.max[j]) 
	newNode->box.max[j] = tempBox.max[j];
    } /* END FOR j */
  } /* END FOR i */

  /* check if only one distr */
  if(distrIndexListN == 1)
    {
        /* OK, this becomes a leaf node */
      newNode->child1 = NULL;
      newNode->child2 = NULL;
      newNode->distrIndex = distrIndexList[0];
      return newNode;
    }

  /* There is more than one distribution */
  /* Get box that encloses all distributions */
  for(int i = 0; i < NO_OF_DIMENSIONS; i++) {
    if(rho_alt_1 != NULL) {
      tempBox.min[i] = rho_alt_1[distrIndexList[0]].centerCoords[i];
      tempBox.max[i] = rho_alt_1[distrIndexList[0]].centerCoords[i];
    }
    else {
      tempBox.min[i] = rho_alt_2[distrIndexList[0]].s.centerCoords[i];
      tempBox.max[i] = rho_alt_2[distrIndexList[0]].s.centerCoords[i];
    }
  } /* END FOR i */
  for(int i = 1; i < distrIndexListN; i++) {
    for(int j = 0; j < NO_OF_DIMENSIONS; j++) {
      if(rho_alt_1 != NULL)
	currCoord = rho_alt_1[distrIndexList[i]].centerCoords[j];
      else
	currCoord = rho_alt_2[distrIndexList[i]].s.centerCoords[j];
      if(tempBox.min[j] > currCoord) tempBox.min[j] = currCoord;
      if(tempBox.max[j] < currCoord) tempBox.max[j] = currCoord;
    } /* END FOR j */
  } /* END FOR i */
  
  /* check if all distrs are at the same point */

  bestCoord = -1;
  maxDiff = 0;
  for(int i = 0; i < NO_OF_DIMENSIONS; i++) {
    currDiff = tempBox.max[i] - tempBox.min[i];
    if(currDiff > maxDiff)
      {
	bestCoord = i;
	maxDiff = currDiff;
      }
  } /* END FOR i */
  bool samePoint = false;
  if(bestCoord < 0)
    samePoint = true;
  else {
    if(maxDiff > COORD_DIFF_FOR_SAMEPOINT_CRITERION) {
      samePoint = false;
    }
    else
      samePoint = true;
  }

  if(samePoint) {
    /* all distrs are at the same point */
    /* sort by extent */
    /* bubble sort (this could be optimized) */
    for(int i = 0; i < (distrIndexListN-1); i++) {
      for(int j = 0; j < (distrIndexListN-1-i); j++) {
	if(rho_alt_1 != NULL) {
	  extent1 = rho_alt_1[distrIndexList[j]].extent;
	  extent2 = rho_alt_1[distrIndexList[j+1]].extent;
	}
	else {
	  extent1 = rho_alt_2[distrIndexList[j]].extent;
	  extent2 = rho_alt_2[distrIndexList[j+1]].extent;
	}
	if(extent1 > extent2) {
	  /* do switch */
	  tempInt = distrIndexList[j];
	  distrIndexList[j] = distrIndexList[j+1];
	  distrIndexList[j+1] = tempInt;
	} /* END IF SWITCH */
      } /* END FOR j bubble sort */
    } /* END FOR i bubble sort       */
      /* check sort */
    for(int i = 0; i < (distrIndexListN-1); i++) {
      if(rho_alt_1 != NULL) {
	extent1 = rho_alt_1[distrIndexList[i]].extent;
	extent2 = rho_alt_1[distrIndexList[i+1]].extent;
      }
      else {
	extent1 = rho_alt_2[distrIndexList[i]].extent;
	extent2 = rho_alt_2[distrIndexList[i+1]].extent;
      }
      if(extent1 > extent2)
	throw std::runtime_error("error in BuildRhoTreeBranch: list not sorted.");
    } /* END FOR i check sort */

      /* create 2 new boxes: small extent and large extent */
    n1 = distrIndexListN / 2;
    n2 = distrIndexListN - n1;
  }
  else {
    /* all distrs are NOT at the same point */
    /* Compute limit as midpoint between min and max for bestCoord, try to avoid rounding errors. */
    real tmpDiff = tempBox.max[bestCoord] - tempBox.min[bestCoord];
    real limit = tempBox.min[bestCoord] + tmpDiff / 2;

    std::vector<int> tempList(distrIndexListN);
    n1 = 0;
    n2 = 0;
    for(int i = 0; i < distrIndexListN; i++) {
      if(rho_alt_1 != NULL)
	testCoord = rho_alt_1[distrIndexList[i]].centerCoords[bestCoord];
      else
	testCoord = rho_alt_2[distrIndexList[i]].s.centerCoords[bestCoord];
      if(testCoord > limit) {
	tempList[n1] = distrIndexList[i];
	n1++;
      }
      else {
	tempList[distrIndexListN-1-n2] = distrIndexList[i];
	n2++;
      }
    } /* END FOR i */
    if((n1 == 0) || (n2 == 0)) {
      do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "Error in BuildRhoTreeBranch (after split): "
		"n1 = %i, n2 = %i\n", n1, n2);
      do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "maxDiff = %33.22f = %6.3g", (double)maxDiff, (double)maxDiff);
      do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "distrIndexListN = %d", distrIndexListN);
      return NULL;
    }

    memcpy(distrIndexList, &tempList[0], distrIndexListN * sizeof(int));
  }
  if((n1 == 0) || (n2 == 0)) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "Error in BuildRhoTreeBranch: n1 = %i, n2 = %i\n", n1, n2);
    return NULL;
  }
  rhoTreeNode* child1 = BuildRhoTreeBranch(noOfDistributionsTot, rho_alt_1, rho_alt_2,
					   n1, distrIndexList, targetRhoError);
  if(child1 == NULL)
    return NULL;
  rhoTreeNode* child2 = BuildRhoTreeBranch(noOfDistributionsTot, rho_alt_1, rho_alt_2,
					   n2, distrIndexList + n1, targetRhoError);
  if(child2 == NULL)
    return NULL;
  newNode->child1 = child1;
  newNode->child2 = child2;
  newNode->distrIndex = -1;

  return newNode;
} /* END  */


static rhoTreeNode* 
BuildRhoTree(int noOfDistributions,
             DistributionSpecStruct* rho_alt_1,
             ShellSpecStructWithExtent* rho_alt_2,
             real targetRhoError)
{
  rhoTreeNode* rootNode;
  int i;
  real targetError, arg, r1;
  DistributionSpecStruct* distr;

  if(rho_alt_1 != NULL)
    {
        /* compute extent for each distribution in list */
      for(i = 0; i < noOfDistributions; i++)
        {
          distr = &rho_alt_1[i];
          targetError = distr->coeff / 1e20;
          arg = distr->coeff / targetError;
          r1 = template_blas_log(arg);
          if(r1 < 0) r1 *= -1;
          distr->extent = template_blas_sqrt(r1 / distr->exponent);
        } /* END FOR i */
    }

  /* set up initial index list: all distributions included */
  std::vector<int> distrIndexList(noOfDistributions);
  for(i = 0; i < noOfDistributions; i++)
    distrIndexList[i] = i;

  rootNode = BuildRhoTreeBranch(noOfDistributions, rho_alt_1, rho_alt_2, 
                                noOfDistributions, &distrIndexList[0], 
                                targetRhoError);

  if(rootNode == NULL)
    throw std::runtime_error("error in BuildRhoTreeBranch.");

  return rootNode;
} /* END BuildRhoTree */






static void free_rho_tree_memory(rhoTreeNode* rootNode)
{
  rhoTreeNode* child1;
  rhoTreeNode* child2;
  child1 = rootNode->child1;
  child2 = rootNode->child2;
  if(child1 != NULL)
    free_rho_tree_memory(child1);
  if(child2 != NULL)
    free_rho_tree_memory(child2);
  delete rootNode;
} /* END free_rho_tree_memory */


static int round_real(real x)
{
  int x1, x2;
  real err1, err2;

  x1 = (int)x;
  x2 = x1 + 1;
  err1 = x - (real)x1;
  err2 = (real)x2 - x;
  if(err1 <= err2)
    return x1;
  else
    return x2;
}


static void getSubBox(const BoxStruct & startBox, BoxStruct & subBox,
		      int Nx, int Ny, int Nz, int ix, int iy, int iz) {
  subBox.min[0] = startBox.min[0] + (real)(ix + 0) * (startBox.max[0] - startBox.min[0]) / Nx;
  subBox.max[0] = startBox.min[0] + (real)(ix + 1) * (startBox.max[0] - startBox.min[0]) / Nx;
  subBox.min[1] = startBox.min[1] + (real)(iy + 0) * (startBox.max[1] - startBox.min[1]) / Ny;
  subBox.max[1] = startBox.min[1] + (real)(iy + 1) * (startBox.max[1] - startBox.min[1]) / Ny;
  subBox.min[2] = startBox.min[2] + (real)(iz + 0) * (startBox.max[2] - startBox.min[2]) / Nz;
  subBox.max[2] = startBox.min[2] + (real)(iz + 1) * (startBox.max[2] - startBox.min[2]) / Nz;
}


typedef real coor3DPtr[3];

static void*
compute_grid_thread_func(void* arg)
{
  try {
  int maxNoOfPoints;
  int noOfShells;
  int noOfNonzeroBasisFuncs;
  int currShellNo, prevShellNo, tempInt;
  int writeResultsToFile;
  int noOfWrittenBatches, noOfGridPoints;
  BoxStruct startBox;
  BoxStruct subBox;
  DensitySpecStruct* density;
  compute_grid_thread_func_struct* inputParams;
  rhoTreeNode* rhoTreeRootNode;
  rhoTreeNode* rhoTreeRootNodeShells;
  int m, ii, jj;
  int nFunctions, count, nPoints, nblocks, blockStarted;
  int startShellNo, NthisWrite;
  FILE* gridFile;
  int jobCount, assignedJobNumber;
  ShellSpecStruct* currShell;
  
  /* get hold of input params */
  inputParams = (compute_grid_thread_func_struct*)arg;
  inputParams->resultCode = -1; // set to zero on success.
  density = inputParams->density;
  memcpy(&startBox, inputParams->startBox, sizeof(BoxStruct));
  rhoTreeRootNode = inputParams->rhoTreeRootNode;
  rhoTreeRootNodeShells = inputParams->rhoTreeRootNodeShells;
  noOfShells = density->noOfShells;
  int Nx = inputParams->Nx;
  int Ny = inputParams->Ny;
  int Nz = inputParams->Nz;
  gridFile = inputParams->gridFile;
  bool generateSparsePatternOnly = inputParams->generateSparsePatternOnly;
  Dft::SparsePattern* sparsePattern = inputParams->sparsePattern;

  compute_grid_for_box_params_struct paramsStruct(inputParams->bis);
  std::vector<real> dmagao(density->nbast*DFT_MAX_BLLEN);
  paramsStruct.dmagao = &dmagao[0];

  writeResultsToFile = 1;

  do_output(LOG_CAT_INFO, LOG_AREA_DFT, "thread %i entering compute_grid_thread_func..", 
            inputParams->threadNo);

  /* allocate memory */
  maxNoOfPoints = FILE_BATCH_N;

  real (*coor)[3] = new real[maxNoOfPoints][3];
  std::vector<real> weight(maxNoOfPoints);
  std::vector<real> coorx(maxNoOfPoints);
  std::vector<real> coory(maxNoOfPoints);
  std::vector<real> coorz(maxNoOfPoints);
  std::vector<int> nonZeroShellIndexList(noOfShells);
  std::vector<int> nonZeroBasisFuncIndexList(density->nbast);
  std::vector<real> workList(density->nbast * MAX_NO_OF_POINTS_PER_BATCH);
  std::vector<DistributionSpecStruct> rhoForSubBox(inputParams->maxNoOfRelevantDistrsPerBox);
  std::vector<int> tempList(inputParams->maxNoOfRelevantDistrsPerBox);
  std::vector<int> listShlblocks(MAX_NO_OF_SHLBLOCKS * 2);

  int (*listShlblocks_otherformat)[2];
  listShlblocks_otherformat = new int[MAX_NO_OF_SHLBLOCKS][2];

  /* get initial assignedJobNumber */
  pthread_mutex_lock(inputParams->jobMutex);
  assignedJobNumber = *inputParams->currJobNumber;
  *inputParams->currJobNumber += N_BATCH_JOBS;
  pthread_mutex_unlock(inputParams->jobMutex);
  jobCount = 0;
  noOfWrittenBatches = 0;
  ComputeGridResultValuesStruct resultValues;
  noOfGridPoints = 0;
  for(int i = 0; i < Nx; i++) {
    for(int j = 0; j < Ny; j++) {
      for(int k = 0; k < Nz; k++) {

	jobCount++;
	if(jobCount >= (assignedJobNumber + N_BATCH_JOBS))
	  {
	    /* get new assignedJobNumber */
	    pthread_mutex_lock(inputParams->jobMutex);
	    assignedJobNumber = *inputParams->currJobNumber;
	    *inputParams->currJobNumber += N_BATCH_JOBS;
	    pthread_mutex_unlock(inputParams->jobMutex);
	  }
	if(jobCount < assignedJobNumber)
	  continue;
	if(assignedJobNumber > (Nx*Ny*Nz))
	  continue;

	/* determine current sub-box */
	getSubBox(startBox, subBox, Nx, Ny, Nz, i, j, k);

	/* get list of non-zero shells for current sub-box */
	int noOfNonzeroShells = get_rhotree_indexes_for_box(&nonZeroShellIndexList[0],
							    nonZeroShellIndexList.size(),
							    &rhoTreeRootNodeShells[0], 
							    &subBox);
	if(noOfNonzeroShells < 0)
	  throw std::runtime_error("error in get_distrs_for_box");
	if(noOfNonzeroShells == 0)
	  continue;

	/* sort list of non-zero shells (bubble sort, could be optimized) */
	for(int kk = 0; kk < (noOfNonzeroShells - 1); kk++)
	  {
	    for(jj = 0; jj < (noOfNonzeroShells - 1 - kk); jj++)
	      {
		if(nonZeroShellIndexList[jj] > 
		   nonZeroShellIndexList[jj+1])
		  {
		    tempInt = nonZeroShellIndexList[jj];
		    nonZeroShellIndexList[jj] = 
		      nonZeroShellIndexList[jj+1];
		    nonZeroShellIndexList[jj+1] = tempInt;
		  }
	      } /* END FOR jj */
	  } /* END FOR kk */

	/* translate list of nonzero shells to list of  */
	/* nonzero contracted distributions */
	noOfNonzeroBasisFuncs = 0;
	for(int kk = 0; kk < noOfNonzeroShells; kk++)
	  {
	    currShell = &density->shellList[nonZeroShellIndexList[kk]].s;
	    nFunctions = 1 + 2 * 
	      currShell->shellType;
	    for(ii = 0; ii < nFunctions; ii++)
	      {
		nonZeroBasisFuncIndexList[noOfNonzeroBasisFuncs] = 
		  currShell->startIndexInMatrix + ii;
		noOfNonzeroBasisFuncs++;
	      } /* END FOR ii      */
	  } /* END FOR kk */
	if(noOfNonzeroBasisFuncs > density->nbast)
	  throw std::runtime_error("error: (noOfNonzeroBasisFuncs > nbast)");

	/* make block-list of non-zero shells to write to file */
	nblocks = 0;
	blockStarted = 0;
	startShellNo = -1;
	prevShellNo = -1;
	for(int kk = 0; kk < noOfNonzeroShells; kk++)
	  {
	    currShellNo = nonZeroShellIndexList[kk];
	    if(blockStarted == 0)
	      {
		blockStarted = 1;
		startShellNo = currShellNo;
	      }
	    else
	      {
		if(currShellNo != (prevShellNo + 1))
		  {
		    /* register previous block */
		    listShlblocks[nblocks*2] = startShellNo; // + 1 here??
		    listShlblocks[nblocks*2+1] = prevShellNo+1; // + 1 here??
		    nblocks++;
		    startShellNo = currShellNo;
		  }
	      }
	    prevShellNo = currShellNo;
	  } /* END FOR kk */
	if(blockStarted == 1)
	  {
	    /* register previous block */
	    listShlblocks[nblocks*2] = startShellNo; // + 1 here??
	    listShlblocks[nblocks*2+1] = prevShellNo+1; // + 1 here??
	    nblocks++;
	  }
	for(int kk = 0; kk < nblocks; kk++) {
	  listShlblocks_otherformat[kk][0] = listShlblocks[kk*2+0];
	  listShlblocks_otherformat[kk][1] = listShlblocks[kk*2+1];
	}

	nPoints = 0;
	if(!generateSparsePatternOnly) {
	  /* get list of relevant distributions for sub-box */
	  count = get_rhotree_indexes_for_box(&tempList[0], tempList.size(), rhoTreeRootNode, &subBox);
	  if(count < 0)
	    throw std::runtime_error("error in get_distrs_for_box");
	  if(count == 0)
	    continue;
	  assert(count <= inputParams->maxNoOfRelevantDistrsPerBox);

	  for(m = 0; m < count; m++)
	    memcpy(&rhoForSubBox[m], 
		   &density->distrList[tempList[m]], 
		   sizeof(DistributionSpecStruct));

	  real Iexact = 0;
	  for(int kk = 0; kk < count; kk++)
	    Iexact += compute_integral_over_box(&rhoForSubBox[kk], &subBox);

	  /* create grid for sub-box */

	  memcpy(&paramsStruct.density, 
		 density, 
		 sizeof(DensitySpecStruct));
	  paramsStruct.density.noOfDistributions = count;
	  paramsStruct.density.distrList = &rhoForSubBox[0];
	  paramsStruct.gridGenerationParams = inputParams->gridGenerationParams;
	  paramsStruct.nonZeroBasisFuncIndexList = &nonZeroBasisFuncIndexList[0];
	  paramsStruct.noOfNonzeroBasisFuncs = noOfNonzeroBasisFuncs;
	  paramsStruct.nonZeroShellsIndexList = &nonZeroShellIndexList[0];
	  paramsStruct.noOfNonzeroShells = noOfNonzeroShells;
	  paramsStruct.nShlblocks = nblocks;
	  paramsStruct.listShlblocks_otherformat = listShlblocks_otherformat;

	  /* Create DFT integrator. */
	  int ndmat = 1;
	  paramsStruct.dftIntegrator = dft_integrator_bl_new(selected_func, ndmat,
							     DFT_MAX_BLLEN, false, inputParams->bis);
	  paramsStruct.dftIntegrator->shl_bl_cnt = nblocks;
	  for(int kk = 0; kk < nblocks; kk++)
	    for(int mm = 0; mm < 2; mm++)
	      paramsStruct.dftIntegrator->shlblocks[kk][mm] = listShlblocks_otherformat[kk][mm];
	  ergoShellsToOrbs(&paramsStruct.dftIntegrator->shl_bl_cnt, paramsStruct.dftIntegrator->shlblocks, 
			   paramsStruct.dftIntegrator->bas_bl_cnt, paramsStruct.dftIntegrator->basblocks,
			   inputParams->bis);

	  /* Setup local full density matrix for current box. */
	  int nnzbf = noOfNonzeroBasisFuncs;
	  paramsStruct.localFullDensityMatrix.resize(nnzbf*nnzbf);
	  for(int kk = 0; kk < nnzbf; kk++) {
	    int kkIndex = nonZeroBasisFuncIndexList[kk];
	    for(int mm = 0; mm < nnzbf; mm++) {
	      int mmIndex = nonZeroBasisFuncIndexList[mm];
	      real dmatElement;
              dmatElement = paramsStruct.density.dmat->at(kkIndex, mmIndex);
	      paramsStruct.localFullDensityMatrix[kk*nnzbf+mm] = dmatElement;
	    }
	  }

	  nPoints = compute_grid_for_box(&paramsStruct,
					 maxNoOfPoints,
					 &coor[0],
					 &weight[0],
					 &subBox,
					 Iexact,
					 &workList[0],
					 resultValues,
					 false);
	  if(nPoints < 0)
	    throw std::runtime_error("error in compute_grid_for_box");

	  dft_integrator_bl_free(paramsStruct.dftIntegrator);

	  if(nPoints == 0)
	    continue;

	  noOfGridPoints += nPoints;
	} // end if (!generateSparsePatternOnly)

	if(writeResultsToFile == 1)
	  {
	    /* set up separate x, y, z vectors for writing to file */
	    if(nPoints > maxNoOfPoints)
	      throw std::runtime_error("error in HiCu compute_grid_thread_func: (nPoints > maxNoOfPoints).");
	    for(int kk = 0; kk < nPoints; kk++)
	      {
		coorx[kk] = coor[kk][0];
		coory[kk] = coor[kk][1];
		coorz[kk] = coor[kk][2];
	      }

	    /* write grid points to file */
	    int nPointsLeft = nPoints;
	    pthread_mutex_lock(inputParams->fileMutex);
	    while(nPointsLeft > 0)
	      {
		if(nPointsLeft <= MAX_NO_OF_POINTS_PER_WRITE)
		  NthisWrite = nPointsLeft;
		else
		  NthisWrite = MAX_NO_OF_POINTS_PER_WRITE;
		fwrite(&NthisWrite, sizeof(int), 1, gridFile);
		fwrite(&nblocks, sizeof(int), 1, gridFile);
		fwrite(&listShlblocks[0], sizeof(int), 2*nblocks, gridFile);
		fwrite(&(coor[nPoints-nPointsLeft][0]),
		       sizeof(real), 3*NthisWrite, gridFile);
		fwrite(&weight[nPoints-nPointsLeft], 
		       sizeof(real), NthisWrite, gridFile);
		nPointsLeft -= NthisWrite;
		noOfWrittenBatches++;                   
	      } /* END WHILE points left */
	    /* Update counters for plot. Note that this is also protected by "fileMutex" being locked.  */
	    tripleVectorOfInt* counterArrForPlot = inputParams->counterArrForPlot;
	    for(int kk = 0; kk < nPoints; kk++) {
	      int idxList[3];
	      for(int coordIdx = 0; coordIdx < 3; coordIdx++) {
		int idx = (int)(HICU_GRID_PLOT_RESOLUTION * 
				((coor[kk][coordIdx] - startBox.min[coordIdx]) / (startBox.max[coordIdx] - startBox.min[coordIdx])));
		if(idx < 0 || idx >= HICU_GRID_PLOT_RESOLUTION)
		  throw std::runtime_error("error in HiCu compute_grid_thread_func: trouble getting indexes for plot counters.");
		idxList[coordIdx] = idx;
	      }
	      int ix = idxList[0];
	      int iy = idxList[1];
	      int iz = idxList[2];
	      (*counterArrForPlot)[ix][iy][iz] ++;
	    }
	    pthread_mutex_unlock(inputParams->fileMutex);
	  } /* end if writeResultsToFile */

	/* Add to sparsePattern if needed. */
	if(sparsePattern) {
	  pthread_mutex_lock(inputParams->fileMutex);
	  sparsePattern->add(nblocks, listShlblocks_otherformat); 
	  pthread_mutex_unlock(inputParams->fileMutex);
	}

      } /* END FOR k */
    } /* END FOR j */
  } /* END FOR i */
  
  do_output(LOG_CAT_INFO, LOG_AREA_DFT, "thread %i loops done, freeing memory..", 
            inputParams->threadNo);

  /* free memory */
  delete [] listShlblocks_otherformat;
  delete [] coor;

  do_output(LOG_CAT_INFO, LOG_AREA_DFT, "thread %i mem freed OK, setting result params..", 
            inputParams->threadNo);
  
  /* report results through input structure */
  inputParams->noOfPoints = noOfGridPoints;
  inputParams->noOfWrittenBatches = noOfWrittenBatches;
  inputParams->resultValues = resultValues;

  do_output(LOG_CAT_INFO, LOG_AREA_DFT, "thread %i exiting compute_grid_thread_func", 
            inputParams->threadNo);
  
  inputParams->resultCode = 0; // set to zero to indicate success.

  }
  catch ( std::exception & e ) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "HiCu Error: Exception caught in compute_grid_thread_func.");
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "what(): %s", e.what());
    return NULL;
  }
  catch (...) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "HiCu Error: Exception caught in compute_grid_thread_func.");
    return NULL;
  }

  return NULL;
} /* END compute_grid_thread_func */



static int compute_grid(
                 const BasisInfoStruct& bis,
                 DensitySpecStruct* density,
		 const GridGenerationParamsStruct & gridGenerationParams,
                 real boxdist,
		 real startBoxSizeDebug,
                 const char* gridFileName,
                 int noOfThreads,
                 bool generateSparsePatternOnly,
                 Dft::SparsePattern* sparsePattern
                 )
{
  BoxStruct startBox;
  BoxStruct tempBox;
  rhoTreeNode* rhoTreeRootNode;
  rhoTreeNode* rhoTreeRootNodeShells;
  real Iexact, absRelError;
  int Nxyz[3]; /* Nx Ny Nz */
  int IexactInteger;
  int noOfDistributions;
  int currJobNumber, noOfShells;

  do_output(LOG_CAT_INFO, LOG_AREA_DFT, "entering compute_grid..");

  noOfShells = density->noOfShells;
  if(noOfShells <= 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "Error in compute_grid: (noOfShells <= 0).");
    return -1;
  }

  noOfDistributions = density->noOfDistributions;
  if(noOfDistributions < 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "Error in compute_grid: (noOfDistributions < 0).");
    return -1;
  }

  do_output(LOG_CAT_INFO, LOG_AREA_DFT, "Entering compute_grid, noOfDistributions = %i, "
            "maxerrorPerBox = %9.3g, targetRhoError = %9.3g",
            noOfDistributions, (double)gridGenerationParams.maxerrorPerBox, (double)gridGenerationParams.targetRhoError);

  do_output(LOG_CAT_INFO, LOG_AREA_DFT, "nbast = %i", density->nbast);

  /* set up starting box */
  get_shell_box(&startBox, &density->shellList[0]);
  for(int i = 1; i < noOfShells; i++) {
    get_shell_box(&tempBox, &density->shellList[i]);
    for(int j = 0; j < NO_OF_DIMENSIONS; j++) {
      if(tempBox.min[j] < startBox.min[j]) 
        startBox.min[j] = tempBox.min[j];
      if(tempBox.max[j] > startBox.max[j]) 
        startBox.max[j] = tempBox.max[j];
    } /* END FOR j */
  } /* END FOR i */

  if(startBoxSizeDebug > 0) {
    for(int j = 0; j < NO_OF_DIMENSIONS; j++) {
      startBox.min[j] = 100*UNIT_one_Angstrom - startBoxSizeDebug;
      startBox.max[j] = 100*UNIT_one_Angstrom + startBoxSizeDebug;
    }
  }

  if(!generateSparsePatternOnly) {
    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "compute_grid starting box:");
    print_box(&startBox);
  }

  Iexact = 0;
  if(!generateSparsePatternOnly) {
    for(int i = 0; i < noOfDistributions; i++)
      Iexact += compute_integral_over_box(&density->distrList[i], &startBox);
    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "Analytical integral over starting box: %.22f", (double)Iexact);
    IexactInteger = round_real(Iexact);
    absRelError = template_blas_fabs((double)IexactInteger - Iexact) / (double)IexactInteger;
    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "Assuming that the correct value is %i, "
              "the relative error is %9.3g", IexactInteger, (double)absRelError);
  }

  rhoTreeRootNode = NULL;
  if(!generateSparsePatternOnly) {
    Util::TimeMeter tmRhoTree;
    rhoTreeRootNode = BuildRhoTree(noOfDistributions, 
                                   density->distrList, 
                                   NULL, 
                                   gridGenerationParams.targetRhoError);
    if(rhoTreeRootNode == NULL)
      {
        do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "error in BuildRhoTree\n");
        return -1;
      }
    tmRhoTree.print(LOG_AREA_DFT, "BuildRhoTree for distrs");
  }

  Util::TimeMeter tmRhoTreeForShells;
  rhoTreeRootNodeShells = BuildRhoTree(noOfShells, 
                                       NULL, 
                                       density->shellList, 
                                       gridGenerationParams.targetRhoError);
  if(rhoTreeRootNodeShells == NULL)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "error in BuildRhoTree.");
      return -1;
    }
  tmRhoTreeForShells.print(LOG_AREA_DFT, "BuildRhoTree for shells");

  /* compute Nx Ny Nz */
  for(int i = 0; i < 3; i++)
      Nxyz[i] = 1 + (int)((startBox.max[i] - startBox.min[i]) / boxdist);
  int Nx = Nxyz[0];
  int Ny = Nxyz[1];
  int Nz = Nxyz[2];

  do_output(LOG_CAT_INFO, LOG_AREA_DFT, "boxdist = %f, Nx = %i, Ny = %i, Nz = %i, Ntot = %i", 
            (double)boxdist, Nx, Ny, Nz, Nx*Ny*Nz);

  /* Now go through all boxes to find the largest number of relevant
     distributions for any single box, since this number is needed to
     allocate work space later. */
  int maxNoOfRelevantDistrsPerBox = 0;
  if(!generateSparsePatternOnly) {
    for(int i = 0; i < Nx; i++)
      for(int j = 0; j < Ny; j++)
	for(int k = 0; k < Nz; k++) {
	  BoxStruct subBox;
	  getSubBox(startBox, subBox, Nx, Ny, Nz, i, j, k);
	  int count = get_rhotree_indexes_for_box(NULL, 0, rhoTreeRootNode, &subBox);
	  if(count > maxNoOfRelevantDistrsPerBox)
	    maxNoOfRelevantDistrsPerBox = count;
	}
    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "maxNoOfRelevantDistrsPerBox = %9d", maxNoOfRelevantDistrsPerBox);
  }

  FILE* gridFile = NULL;
  if(!generateSparsePatternOnly) {
    /* create grid file */
    gridFile = fopen(gridFileName, "wb");
    if(gridFile == NULL)
      {
        do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "error opening grid file '%s' for writing", gridFileName);
        return -1;
      }
  }

  std::vector< std::vector< std::vector<int> > > counterArrForPlot(HICU_GRID_PLOT_RESOLUTION);
  for(int i = 0; i < HICU_GRID_PLOT_RESOLUTION; i++) {
    counterArrForPlot[i].resize(HICU_GRID_PLOT_RESOLUTION);
    for(int j = 0; j < HICU_GRID_PLOT_RESOLUTION; j++) {
      counterArrForPlot[i][j].resize(HICU_GRID_PLOT_RESOLUTION);
      for(int k = 0; k < HICU_GRID_PLOT_RESOLUTION; k++)
	counterArrForPlot[i][j][k] = 0;
    }
  }

  /* up to this point there is no parallellization */
  /* this is where we start to think about threading */

  pthread_mutex_t fileMutex = PTHREAD_MUTEX_INITIALIZER;
  pthread_mutex_t jobMutex = PTHREAD_MUTEX_INITIALIZER;

  std::vector<compute_grid_thread_func_struct*> threadParamsList(noOfThreads);
  for(int i = 0; i < noOfThreads; i++)
    threadParamsList[i] = new compute_grid_thread_func_struct(bis);

  currJobNumber = 1;

  for(int i = 0; i < noOfThreads; i++) {
    threadParamsList[i]->density = density;
    threadParamsList[i]->rhoTreeRootNode = rhoTreeRootNode;
    threadParamsList[i]->rhoTreeRootNodeShells = rhoTreeRootNodeShells;
    threadParamsList[i]->gridGenerationParams = gridGenerationParams;
    threadParamsList[i]->gridFile = gridFile;
    threadParamsList[i]->startBox = &startBox;
    threadParamsList[i]->Nx = Nx;
    threadParamsList[i]->Ny = Ny;
    threadParamsList[i]->Nz = Nz;
    threadParamsList[i]->maxNoOfRelevantDistrsPerBox = maxNoOfRelevantDistrsPerBox;
    threadParamsList[i]->fileMutex = &fileMutex;
    threadParamsList[i]->jobMutex = &jobMutex;
    threadParamsList[i]->currJobNumber = &currJobNumber;
    threadParamsList[i]->noOfPoints = -1;
    threadParamsList[i]->noOfWrittenBatches = 0;
    threadParamsList[i]->generateSparsePatternOnly = generateSparsePatternOnly;
    threadParamsList[i]->sparsePattern = sparsePattern;
    threadParamsList[i]->counterArrForPlot = &counterArrForPlot;
    threadParamsList[i]->threadNo = i;
  } /* END FOR i */

  if(noOfThreads == 1) {
    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "compute_grid: (noOfThreads == 1), no threads created.");
    compute_grid_thread_func(threadParamsList[0]);
    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "Single call to compute_grid_thread_func done.");
  }
  else {
    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "Starting %i threads.", noOfThreads);

    /* start threads */
    for(int i = 0; i < noOfThreads; i++) {
      if(pthread_create(&threadParamsList[i]->thread, 
			NULL, 
			compute_grid_thread_func, 
			threadParamsList[i]) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "Error in pthread_create for thread %i", i);
	  do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "waiting for already created threads..");
	  for(int j = 0; j < i; j++) {
	    if(pthread_join(threadParamsList[j]->thread, NULL) != 0)
	      do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "Error in pthread_join for thread %i", j);
	  } /* END FOR j */
	  do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "all threads finished, returning error code");
	  return -1;
	}
    } /* END FOR i */

    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "%i threads started OK.", noOfThreads);

    /* wait for threads to finish */
    for(int i = 0; i < noOfThreads; i++) {
      if(pthread_join(threadParamsList[i]->thread, NULL) != 0)
        do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "Error in pthread_join for thread %i", i);
    } /* END FOR i */
    
    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "all %i threads have finished:", noOfThreads);
    for(int i = 0; i < noOfThreads; i++)
      do_output(LOG_CAT_INFO, LOG_AREA_DFT, "thread %2i noOfWrittenBatches = %6i", 
                i, threadParamsList[i]->noOfWrittenBatches);
  } // end if using threads
  
  /* now all threads have finished, check for errors */
  for(int i = 0; i < noOfThreads; i++) {
    if(threadParamsList[i]->noOfPoints < 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "error in compute_grid_thread_func"
                " for thread %i\n", i);
      return -1;
    }
  } /* END FOR i */



  int noOfGridPoints = 0;
  int noOfWrittenBatches = 0;
  real totalIntegralResultNumerical = 0;
  real totalIntegralResultAnalytical = 0;
  real totalIntegralResultEnergy = 0;
  real estimatedIntegralErrorDensity = 0;
  real estimatedIntegralErrorEnergy = 0;
  for(int i = 0; i < noOfThreads; i++)
    {
      noOfGridPoints += threadParamsList[i]->noOfPoints;
      noOfWrittenBatches += threadParamsList[i]->noOfWrittenBatches;
      totalIntegralResultNumerical += threadParamsList[i]->resultValues.totalIntegralResultNumerical;
      totalIntegralResultAnalytical += threadParamsList[i]->resultValues.totalIntegralResultAnalytical;
      totalIntegralResultEnergy += threadParamsList[i]->resultValues.totalIntegralResultEnergy;
      estimatedIntegralErrorDensity += threadParamsList[i]->resultValues.estimatedIntegralErrorDensity;
      estimatedIntegralErrorEnergy += threadParamsList[i]->resultValues.estimatedIntegralErrorEnergy;
    } /* END FOR i */

  if(gridFile)
    fclose(gridFile);

  if(!generateSparsePatternOnly) {
    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "noOfWrittenBatches = %i", noOfWrittenBatches);
    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "compute_grid ending OK, noOfGridPoints = %i",
              noOfGridPoints);
    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "compute_grid maxerrorPerBox param          = %25.15f = %9.4g", 
	      (double)gridGenerationParams.maxerrorPerBox, (double)gridGenerationParams.maxerrorPerBox);
    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "compute_grid totalIntegralResultAnalytical = %25.15f", 
	      (double)totalIntegralResultAnalytical);
    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "compute_grid totalIntegralResultNumerical  = %25.15f", 
	      (double)totalIntegralResultNumerical);
    ergo_real absDiff = template_blas_fabs(totalIntegralResultAnalytical - totalIntegralResultNumerical);
    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "compute_grid numer/analy integral abs diff = %25.15f = %9.4g", 
	      (double)absDiff, (double)absDiff);
    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "compute_grid totalIntegralResultEnergy     = %25.15f", 
	      (double)totalIntegralResultEnergy);
    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "compute_grid estimatedIntegralErrorDensity = %25.15f = %9.4g", 
	      (double)estimatedIntegralErrorDensity, (double)estimatedIntegralErrorDensity);
    do_output(LOG_CAT_INFO, LOG_AREA_DFT, "compute_grid estimatedIntegralErrorEnergy  = %25.15f = %9.4g", 
	      (double)estimatedIntegralErrorEnergy, (double)estimatedIntegralErrorEnergy);

#if 0
    // Create m-file for plot.
    FILE* mfile = fopen("grid_plot_file_3d.m", "wt");
    fprintf(mfile, "M = [\n");
    for(int i = 0; i < HICU_GRID_PLOT_RESOLUTION; i++)
      for(int j = 0; j < HICU_GRID_PLOT_RESOLUTION; j++)
	for(int k = 0; k < HICU_GRID_PLOT_RESOLUTION; k++) {
	  double x = startBox.min[0] + (double)i * (startBox.max[0] - startBox.min[0]) / HICU_GRID_PLOT_RESOLUTION;
	  double y = startBox.min[1] + (double)j * (startBox.max[1] - startBox.min[1]) / HICU_GRID_PLOT_RESOLUTION;
	  double z = startBox.min[2] + (double)k * (startBox.max[2] - startBox.min[2]) / HICU_GRID_PLOT_RESOLUTION;
	  fprintf(mfile, "%15.5f %15.5f %15.5f %9d\n", x, y, z, counterArrForPlot[i][j][k]);
	}
    fprintf(mfile, "];\n");
    fclose(mfile);
    mfile = fopen("grid_plot_file_2d_z.m", "wt");
    fprintf(mfile, "M = [\n");
    for(int i = 0; i < HICU_GRID_PLOT_RESOLUTION; i++) {
      for(int j = 0; j < HICU_GRID_PLOT_RESOLUTION; j++) {
	int count = 0;
	for(int k = 0; k < HICU_GRID_PLOT_RESOLUTION; k++)
	  count += counterArrForPlot[i][j][k];
	fprintf(mfile, " %9d", count);
      }
      fprintf(mfile, "\n");
    }
    fprintf(mfile, "];\n");
    fclose(mfile);
#endif
  }
  
  if(!generateSparsePatternOnly)
    free_rho_tree_memory(rhoTreeRootNode);
  free_rho_tree_memory(rhoTreeRootNodeShells);

  for(int i = 0; i < noOfThreads; i++)
    delete threadParamsList[i];

  return noOfGridPoints;
} /* END compute_grid */


static int 
do_merge_sort_distrs(int n, 
                     DistributionSpecStruct* list, 
                     DistributionSpecStruct* workList)
{
    /* merge sort:  */
    /* first sort the first half, */
    /* then sort the second half, */
    /* then merge results to form final sorted list. */
  int n1, n2, nn, decision, i1, i2, i;
  DistributionSpecStruct* d1;
  DistributionSpecStruct* d2;

  if(n < 1) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "Error in do_merge_sort_distrs: (n < 1).");
    return -1;
  }
  if(n == 1)
    return 0;
  
  n1 = n / 2;
  n2 = n - n1;

  /* sort first half */
  if(do_merge_sort_distrs(n1, list, workList) != 0)
    return -1;

  /* sort second half */
  if(do_merge_sort_distrs(n2, &list[n1], workList) != 0)
    return -1;

  /* merge results */
  nn = 0;
  i1 = 0;
  i2 = 0;
  while(nn < n)
    {
      if((i1 < n1) && (i2 < n2))
        {
            /* compare */
          d1 = &list[i1];
          d2 = &list[n1+i2];
          decision = 0;
          for(i = 0; i < 3; i++)
            {
              if(decision == 0)
                {
                  if(d1->monomialInts[i] != d2->monomialInts[i])
                    {
                      if(d1->monomialInts[i] > d2->monomialInts[i])
                        decision = 1;
                      else
                        decision = 2;
                    }
                } /* END IF (decision == 0) */
            } /* END FOR i */
          if(decision == 0)
            {
                /* check exponents */
              if(d1->exponent > d2->exponent)
                decision = 1;
              else
                decision = 2;
            }
        }
      else
        {
          if(i1 == n1)
              decision = 2;
          else
              decision = 1;
        }
      if(decision <= 0) {
        do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "Error in do_merge_sort_distrs: (decision <= 0).");
        return -1;
      }
      if(decision == 1) {
        memcpy(&workList[nn], &list[i1], sizeof(DistributionSpecStruct));
        i1++;
      }
      else {
        memcpy(&workList[nn], &list[n1+i2], sizeof(DistributionSpecStruct));
        i2++;
      }
      nn++;
    } /* END WHILE (nn < n) */
  if(i1 != n1) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "Error in do_merge_sort_distrs: (i1 != n1).");
    return -1;
  }
  if(i2 != n2) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "Error in do_merge_sort_distrs: (i2 != n2).");
    return -1;
  }
  if(nn != n) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "Error in do_merge_sort_distrs: (nn != n).");
    return -1;
  }
  memcpy(list, workList, n * sizeof(DistributionSpecStruct));
  return 0;
} /* END do_merge_sort_distrs */



static int
compute_extent_for_shells(ShellSpecStructWithExtent* shellList, 
                          const BasisInfoStruct& bis, 
                          real targetRhoError) {
  /* We do this using existing function for getting the extent of all basis functions. */
  std::vector<real> basisFuncExtentList(bis.noOfBasisFuncs);
  ergo_real maxAbsDensityMatrixElement = 1.0; /* FIXME: use correct value here. */
  real maxAbsValue = targetRhoError / (get_max_basis_func_abs_value(bis) * maxAbsDensityMatrixElement);
  get_basis_func_extent_list(bis, &basisFuncExtentList[0], maxAbsValue);
  for(int i = 0; i < bis.noOfShells; i++) {
    ShellSpecStructWithExtent* currShell = &shellList[i];
    real largestExtent = 0;
    int startIdx = currShell->s.startIndexInMatrix;
    for(int j = 0; j < currShell->s.noOfBasisFuncs; j++) {
      ergo_real currBasisFuncExtent = basisFuncExtentList[startIdx+j];
      if(currBasisFuncExtent > largestExtent)
	largestExtent = currBasisFuncExtent;
    }
    currShell->extent = largestExtent;
  }
  return 0;
}


static int get_product_distrs(const BasisInfoStruct& bis,
                              const Dft::Matrix& dmat,
                              real targetRhoError,
                              DistributionSpecStruct* rho,  /* may be NULL. */
                              int maxCount /* only used if rho != NULL. */
                              ) {
  Util::TimeMeter tm;
  const int MAX_DISTRS_IN_TEMP_LIST = 4444;
  int nBasisFuncs = bis.noOfBasisFuncs;
  std::vector<ergo_real> basisFuncExtentList(nBasisFuncs);
  ergo_real maxAbsDensityMatrixElement = 1.0; /* FIXME: use correct value here. */
  real maxAbsValue = targetRhoError / (get_max_basis_func_abs_value(bis) * maxAbsDensityMatrixElement);
  get_basis_func_extent_list(bis, &basisFuncExtentList[0], maxAbsValue);
  ergo_real maxExtent = 0;
  for(int i = 0; i < nBasisFuncs; i++) {
    ergo_real currExtent = basisFuncExtentList[i];
    if(currExtent > maxExtent)
      maxExtent = currExtent;
  }  
  // Create box system.
  std::vector<box_item_struct> itemList(nBasisFuncs);
  for(int i = 0; i < nBasisFuncs; i++) {
    for(int j = 0; j < 3; j++)
      itemList[i].centerCoords[j] = bis.basisFuncList[i].centerCoords[j];
    itemList[i].originalIndex = i;
  }
  ergo_real toplevelBoxSize = 7.0;
  BoxSystem boxSystem;
  if(boxSystem.create_box_system(&itemList[0],
                                 nBasisFuncs,
                                 toplevelBoxSize) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "error in get_product_distrs: error creating box system.");
    return -1;
  }
  std::vector<int> orgIndexList(nBasisFuncs);
  int nn = 0;
  for(int i = 0; i < nBasisFuncs; i++) {
    // Now, instead of looping again over all nBasisFuncs basis
    // functions, we use box system to find relevant ones.
    ergo_real maxDistance = basisFuncExtentList[i] + maxExtent;
    ergo_real coords[3];
    for(int coordNo = 0; coordNo < 3; coordNo++)
      coords[coordNo] = bis.basisFuncList[i].centerCoords[coordNo];
    int nRelevant = boxSystem.get_items_near_point(&itemList[0], coords, maxDistance, &orgIndexList[0]);
    for(int jRelevant = 0; jRelevant < nRelevant; jRelevant++) {
      int j = orgIndexList[jRelevant];
      DistributionSpecStruct tempList[MAX_DISTRS_IN_TEMP_LIST];
      int nPrimitives;
      /* the matrix M is symmetric: include diagonal terms once, */
      /* and include upper off-diagonal terms multiplied by 2 */
      int symmetryFactor;
      if(i == j)
        symmetryFactor = 1;
      else
        symmetryFactor = 2;
      if(i > j)
        continue;
      nPrimitives = 
        get_product_simple_primitives(bis, i,
                                      bis, j,
                                      tempList,
                                      MAX_DISTRS_IN_TEMP_LIST,
                                      DISTR_PRODUCT_THRESHOLD);
      if(nPrimitives < 0)
        throw std::runtime_error("error in get_product_simple_primitives");
      for(int k = 0; k < nPrimitives; k++) {
        DistributionSpecStruct* currDistr = &tempList[k];
        real Mij;
        Mij = dmat.at(i, j);
        real newCoeff = currDistr->coeff * Mij * symmetryFactor;
        if(template_blas_fabs(newCoeff) > DISTR_COEFF_CUTOFF_VALUE) {
          /* add to final list */
          if(rho) {
            if(nn >= maxCount)
              throw std::runtime_error("error: (nn >= maxCount)");
            memcpy(&rho[nn], currDistr, 
                   sizeof(DistributionSpecStruct));
            rho[nn].coeff = newCoeff;
          }
          nn++;
        }
      }
    }
  }
  tm.print(LOG_AREA_DFT, "get_product_distrs");
  return nn;
}


static void
get_shell_list_with_extents(const BasisInfoStruct& bis,
                            int maxCountShellList,
                            ShellSpecStructWithExtent* shellList,
                            real targetRhoError) {
  if(maxCountShellList < bis.noOfShells)
    throw std::runtime_error("Error: (maxCountShellList < bis.noOfShells)");
  for(int i = 0; i < bis.noOfShells; i++) {
    shellList[i].s = bis.shellList[i];
    shellList[i].extent = 0; // to be computed later.
  }
  if(compute_extent_for_shells(shellList, bis, targetRhoError) != 0)
    throw std::runtime_error("Error in compute_extent_for_shells.");
}

static int
get_density(const BasisInfoStruct& bis,
            DistributionSpecStruct* rho,
            int maxCountRho,
            real targetRhoError, 
            int nbast, 
            const Dft::Matrix& dmat,
            BasisFuncStruct* basisFuncList)
{
  Util::TimeMeter tm;
  Util::TimeMeter tmFirstPart;

  do_output(LOG_CAT_INFO, LOG_AREA_DFT, "entering function get_density, targetRhoError = %22.15f", 
	    (double)targetRhoError);

  int nn = get_product_distrs(bis, dmat, targetRhoError, rho, maxCountRho);
  if(nn != maxCountRho)
    throw std::runtime_error("Error in get_density: (nn != maxCountRho).");

  memcpy(basisFuncList, bis.basisFuncList,
         bis.noOfBasisFuncs * sizeof(BasisFuncStruct));  

  do_output(LOG_CAT_INFO, LOG_AREA_DFT, "loop ended OK; list 'rho' created, nn = %i", nn);

  /* Now all distributions are stored in the list 'rho'. */
  /* The number of entries in the list is nn. */
  /* It could happen that all entries are not unique. */
  /* We want to join distributions that have the same center  */
  /* and the same exponent. */
  /* To do this, start with sorting the list by nx, ny, nz, exponent. */
  std::vector<DistributionSpecStruct> workList(nn);
  
  do_output(LOG_CAT_INFO, LOG_AREA_DFT, "calling do_merge_sort_distrs, nn = %i", nn);
  if(do_merge_sort_distrs(nn, &rho[0], &workList[0]) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "error in do_merge_sort_distrs");
    return -1;
  }
  do_output(LOG_CAT_INFO, LOG_AREA_DFT, "do_merge_sort_distrs returned OK");  

  /* check that list is sorted */
  for(int i = 0; i < (nn-1); i++) {
    if(rho[i].exponent < rho[i+1].exponent) {
      int sameYesNo = 1;
      for(int j = 0; j < 3; j++) {
        if(rho[i].monomialInts[j] != rho[i+1].monomialInts[j])
          sameYesNo = 0;
      } /* END FOR j */
      if(sameYesNo == 1) {
        do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "error: distr list NOT properly sorted.");
        return -1;
      }
    }
  } /* END FOR i */
  do_output(LOG_CAT_INFO, LOG_AREA_DFT, "sort checked OK");

  tmFirstPart.print(LOG_AREA_DFT, "get_density first part");

  std::vector<int> markList(nn);
  for(int i = 0; i < nn; i++)
    markList[i] = 0;

  /* Create box system to help finding distrs that have centers that
     are close to eachother in space. */
  std::vector<box_item_struct> itemList(nn);
  for(int i = 0; i < nn; i++) {
    for(int j = 0; j < 3; j++)
      itemList[i].centerCoords[j] = rho[i].centerCoords[j];
    itemList[i].originalIndex = i;
  }
  real toplevelBoxSize = 4.0;
  BoxSystem boxSystem;
  if(boxSystem.create_box_system(&itemList[0],
                                 nn,
                                 toplevelBoxSize) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "error in get_product_distrs: error creating box system.");
    return -1;
  }
  std::vector<int> orgIndexList(nn);

  /* now go through sorted list, joining distributions where possible */
  int icurr = 0;
  int count = 0;
  int firstIndex = 0;
  while(icurr < nn) {
    /* check if this entry has the same nx ny nz as current 'firstIndex' */
    int sameYesNo = 1;
    for(int j = 0; j < 3; j++) {
      if(rho[icurr].monomialInts[j] != rho[firstIndex].monomialInts[j])
	sameYesNo = 0;
    } /* END FOR j */
    /* check exponent */
    real absdiff = template_blas_fabs(rho[icurr].exponent - rho[firstIndex].exponent);
    if(absdiff > EXPONENT_DIFF_LIMIT)
      sameYesNo = 0;
    if(sameYesNo == 0) {
      /* Now take care of all distrs from firstIndex to icurr-1. We
	 know that all of them have identical monomialInts and
	 exponents. */
      for(int j = firstIndex; j < icurr; j++) {
        if(markList[j] == 0) {
          markList[j] = 1;
          /* join distrs that have centers within  */
          /* DISTR_CENTER_DIST_LIMIT of this one */
	  

	  ergo_real maxDistance = DISTR_CENTER_DIST_LIMIT;
	  ergo_real coords[3];
	  for(int coordNo = 0; coordNo < 3; coordNo++)
	    coords[coordNo] = rho[j].centerCoords[coordNo];
	  int nRelevant = boxSystem.get_items_near_point(&itemList[0], coords, maxDistance, &orgIndexList[0]);
          real coeffSum = rho[j].coeff;
	  for(int jRelevant = 0; jRelevant < nRelevant; jRelevant++) {
	    if(orgIndexList[jRelevant] >= j+1 && orgIndexList[jRelevant] < icurr) {
	      int k = orgIndexList[jRelevant];
	      //          for(int k = j+1; k < icurr; k++) {
	      int withinLimit = 1;
	      for(int kk = 0; kk < 3; kk++) {
		real absdiff = template_blas_fabs(rho[j].centerCoords[kk] - 
				    rho[k].centerCoords[kk]);
		if(absdiff > DISTR_CENTER_DIST_LIMIT)
		  withinLimit = 0;
	      } /* END FOR kk */
	      if(withinLimit == 1) {
		coeffSum += rho[k].coeff;
		markList[k] = 1;
	      }
	    } /* end if index within range. */
          } /* end for jRelevant */
          memcpy(&workList[count], 
                 &rho[j], 
                 sizeof(DistributionSpecStruct));
          workList[count].coeff = coeffSum;
          count++;
        } /* END IF (markList[j] == 0) */
      } /* END FOR j */
      firstIndex = icurr;
    } /* end if (sameYesNo == 0) */
    else {
      /* Do nothing here. */
    }
    icurr++;
  } /* END WHILE (icurr < nn) */
  /* take care of last part */
  for(int j = firstIndex; j < nn; j++) {
    if(markList[j] == 0) {
      markList[j] = 1;
      /* join distrs that have centers within  */
      /* DISTR_CENTER_DIST_LIMIT of this one */
      real coeffSum = rho[j].coeff;
      for(int k = j+1; k < nn; k++) {
        int withinLimit = 1;
        for(int kk = 0; kk < 3; kk++) {
          real absdiff = template_blas_fabs(rho[j].centerCoords[kk] - 
			      rho[k].centerCoords[kk]);
          if(absdiff > DISTR_CENTER_DIST_LIMIT)
            withinLimit = 0;
        } /* END FOR kk */
        if(withinLimit == 1) {
	  coeffSum += rho[k].coeff;
	  markList[k] = 1;
	}
      } /* END FOR k */
      memcpy(&workList[count], &rho[j], sizeof(DistributionSpecStruct));
      workList[count].coeff = coeffSum;
      count++;
    } /* END IF (markList[j] == 0) */
  } /* END FOR j */

  for(int j = 0; j < nn; j++) {
    if(markList[j] != 1) {
      do_output(LOG_CAT_ERROR, LOG_AREA_DFT, "Error: (markList[%i] != 1).", j);
      return -1;
    }
  } /* END FOR j */


  /* now move results back to list 'rho',  */
  /* skipping those that have too small coeff */
  int resultCount = 0;
  for(int i = 0; i < count; i++) {
      real sqrtValue = template_blas_sqrt(pi / workList[i].exponent);
      real absvalue = workList[i].coeff * sqrtValue * sqrtValue * sqrtValue;
      if(absvalue < 0) absvalue *= -1;      
      if(absvalue > DISTR_COEFF_CUTOFF_VALUE) {
	memcpy(&rho[resultCount], &workList[i], sizeof(DistributionSpecStruct));
	resultCount++;
      }
  } /* END FOR i */

  do_output(LOG_CAT_INFO, LOG_AREA_DFT, "nn          = %9i", nn);
  do_output(LOG_CAT_INFO, LOG_AREA_DFT, "count       = %9i", count);
  do_output(LOG_CAT_INFO, LOG_AREA_DFT, "resultCount = %9i", resultCount);

  tm.print(LOG_AREA_DFT, "get_density");

  return resultCount;
} /* end get_density */



int hicu_grid_generate(const char* grid_file_name, 
                       const BasisInfoStruct& bis,
		       ergo_real maxError,
		       ergo_real boxSize,
		       ergo_real startBoxSizeDebug,
		       int use_error_per_volume,
		       int do_double_checking,
		       int compare_to_refined,
		       int use_energy_criterion,
		       int use_energy_criterion_only,
		       int do_variation_checking,
                       const Dft::Matrix* dmat,
                       Dft::SparsePattern* sparsePattern, 
                       int nThreads,
                       bool generateSparsePatternOnly) {
  /* Use mutex lock to be sure only one thread at a time executes this
     call. */
  pthread_mutex_lock(&global_main_hicu_mutex);

  output_current_memory_usage(LOG_AREA_DFT, "hicu_grid_generate start");
  Util::TimeMeter tm;
  DensitySpecStruct density;
  int nGridPoints = -1;

  do_output(LOG_CAT_INFO, LOG_AREA_DFT, 
            "hicu_grid_generate checking what kind of density matrix is given.");

  int nbast = bis.noOfBasisFuncs;

  if(dmat->isSparse())
    do_output(LOG_CAT_INFO, LOG_AREA_DFT, 
              "hicu_grid_generate using sparse dmat.");
  else {
    /* check dmat */
    real maxabs = 0;
    for(int i = 0; i < nbast; i++) {
      for(int j = 0; j < nbast; j++) {

        real temp = template_blas_fabs(dmat->at(i,j));
        if(temp > maxabs)
          maxabs = temp;
      }
    }
    do_output(LOG_CAT_INFO, LOG_AREA_DFT,
              "hicu_grid_generate checking dmat: maxabs = %22.15f", 
	      (double)maxabs);
  }

  int noOfShells = bis.noOfShells;

  ergo_real targetRhoError = maxError * TARGET_RHO_ERROR_FACTOR;

  std::vector<ShellSpecStructWithExtent> shellList(noOfShells);
  get_shell_list_with_extents(bis, noOfShells, &shellList[0], targetRhoError);
  output_current_memory_usage(LOG_AREA_DFT, "hicu_grid_generate after getting shellList");

  std::vector<BasisFuncStruct> basisFuncList(nbast);
  output_current_memory_usage(LOG_AREA_DFT, "hicu_grid_generate after allocating basisFuncList");

  /* Call get_product_distrs here to get number of distrs to allocate. */
  int noOfDistributions1 = 0;
  if(!generateSparsePatternOnly) {
    noOfDistributions1 = get_product_distrs(bis, *dmat, targetRhoError,
                                            NULL, 0);
    if(noOfDistributions1 <= 0)
      throw std::runtime_error("Error in hicu_grid_generate: (noOfDistributions1 <= 0).");
  }
  std::vector<DistributionSpecStruct> rho(noOfDistributions1);
  output_current_memory_usage(LOG_AREA_DFT, "hicu_grid_generate after allocating rho");

  do_output(LOG_CAT_INFO, LOG_AREA_DFT, "Calling get_density().");

  int noOfDistributions = get_density(bis,
                                      &rho[0], 
                                      noOfDistributions1, 
                                      targetRhoError, 
                                      nbast, 
                                      *dmat,
                                      &basisFuncList[0]);
  if(noOfDistributions < 0 || noOfDistributions > noOfDistributions1)
    throw std::runtime_error("error in get_density! "
                             "(noOfDistributions < 0 || noOfDistributions > noOfDistributions1).");
  output_current_memory_usage(LOG_AREA_DFT, "hicu_grid_generate after get_density");
  
  density.noOfShells = noOfShells;
  density.shellList = &shellList[0];
  density.nbast = nbast;
  density.dmat = dmat;
  density.basisFuncList = &basisFuncList[0];
  density.noOfDistributions = noOfDistributions;
  density.distrList = &rho[0];

  GridGenerationParamsStruct gridGenerationParams;
  gridGenerationParams.maxerrorPerBox = maxError;
  gridGenerationParams.targetRhoError = targetRhoError;
  gridGenerationParams.doDoubleChecking = do_double_checking;
  gridGenerationParams.compareToRefined = compare_to_refined;
  gridGenerationParams.useEnergyCriterion = use_energy_criterion;
  gridGenerationParams.useEnergyCriterionOnly = use_energy_criterion_only;
  gridGenerationParams.useErrorPerVolume = use_error_per_volume;
  gridGenerationParams.doVariationChecking = do_variation_checking;

  /* get grid */
  nGridPoints = compute_grid(bis,
                             &density,
			     gridGenerationParams,
                             boxSize,
			     startBoxSizeDebug,
                             grid_file_name,
                             nThreads,
                             generateSparsePatternOnly,
                             sparsePattern);
  if(nGridPoints < 0)
    throw std::runtime_error("Error in compute_grid.");

  do_output(LOG_CAT_INFO, LOG_AREA_DFT, "HiCu grid generated OK, nGridPoints = %9d", nGridPoints);

  output_current_memory_usage(LOG_AREA_DFT, "hicu_grid_generate end");
  tm.print(LOG_AREA_DFT, __func__);

  pthread_mutex_unlock(&global_main_hicu_mutex);
  
  return nGridPoints;
}


void
grid_generate_sparse_pattern(const BasisInfoStruct& bis,
			     ergo_real maxError,
			     ergo_real boxSize,
			     ergo_real startBoxSizeDebug,
                             Dft::SparsePattern& sparsePattern) {
  /* Use mutex lock to be sure only one thread at a time executes this
     call. */
  pthread_mutex_lock(&global_main_hicu_mutex);

  output_current_memory_usage(LOG_AREA_DFT, "grid_generate_sparse_pattern start");
  Util::TimeMeter tm;
  DensitySpecStruct density;

  int nbast = bis.noOfBasisFuncs;

  int noOfShells = bis.noOfShells;

  ergo_real targetRhoError = maxError * TARGET_RHO_ERROR_FACTOR;

  std::vector<ShellSpecStructWithExtent> shellList(noOfShells);
  get_shell_list_with_extents(bis, noOfShells, &shellList[0], targetRhoError);
  
  density.noOfShells = noOfShells;
  density.shellList = &shellList[0];
  density.nbast = nbast;
  density.dmat = NULL;
  density.basisFuncList = NULL;
  density.noOfDistributions = 0;
  density.distrList = NULL;

  GridGenerationParamsStruct gridGenerationParams;
  gridGenerationParams.maxerrorPerBox = maxError;
  gridGenerationParams.targetRhoError = targetRhoError;

  /* get grid */
  int nGridPoints = compute_grid(bis,
                                 &density,
				 gridGenerationParams,
                                 boxSize,
				 startBoxSizeDebug,
                                 NULL,
                                 1,
                                 true,
                                 &sparsePattern);
  if(nGridPoints < 0)
    throw std::runtime_error("Error in compute_grid.");

  output_current_memory_usage(LOG_AREA_DFT, "grid_generate_sparse_pattern end");
  tm.print(LOG_AREA_DFT, __func__);

  pthread_mutex_unlock(&global_main_hicu_mutex);
}

