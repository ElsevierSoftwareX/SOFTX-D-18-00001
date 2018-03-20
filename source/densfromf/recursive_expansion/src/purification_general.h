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

/** @file purification_general.h

    @brief Recursive density matrix expansion (or density matrix
    purification).

    @author Anastasia Kruchinina <em>responsible</em>
*/


#ifndef HEADER_PURIFICATION_GENERAL
#define HEADER_PURIFICATION_GENERAL

#include <iostream>
#include <fstream>
#include <sstream>

#include "matrix_typedefs.h"  // definitions of matrix types and interval type
#include "realtype.h"         // definitions of types
#include "matrix_utilities.h"
#include "integral_matrix_wrappers.h"
#include "output.h"
#include "matrix_proxy.h"

#include "puri_info.h"
#include "constants.h"
#include "utilities.h"
#include "units.h"

#include "files_dense.h"
#include "files_sparse.h"
#include "files_sparse_bin.h"

#include "get_eigenvectors.h"

typedef ergo_real real;

/****************  DEFINE  ******************/

// number of additional iterations

/* After that stopping criterion says to stop, we save matrix X into the file,
 * perform additional iterations and read matrix X back.
 * It is done just for testing purposes in case we wish to run a few SCF iterations
 * and use results in the iterations given by the stopping criterion,
 * but at the same time we want to see what happen after stoping criterion says
 * to stop the iterations. */
#define NUM_ADDITIONAL_ITERATIONS    0

// enable more output
//#define DEBUG_PURI_OUTPUT
#define PURI_OUTPUT_NNZ


// enable printf output instead of output to ergoscf.out
// do it if you need run just recursive expansion, not the whole SCF calculations
//#define ENABLE_PRINTF_OUTPUT


//#define SAVE_MATRIX_IN_EACH_ITERATION // save F and X_i for every i



/**********************************/


extern real eucl_acc;
extern real mixed_acc;
extern real TOL_OVERLAPPING_BOUNDS;
extern real THRESHOLD_EIG_TOLERANCE;
extern int  EIG_EMPTY;
extern int  EIG_SQUARE_INT;
extern int  EIG_PROJECTION_INT;
extern int  EIG_POWER_INT;
extern int  EIG_LANCZOS_INT;


/** PurificationGeneral is an abstract class which provides an
 * interface for SP2, SP2ACC and possibly other recursive expansions.
 *
 * \tparam MatrixType Type of a matrix (ex. symmMatrix). */
template<typename MatrixType>
class PurificationGeneral
{
public:
   typedef ergo_real         real;
   typedef intervalType      IntervalType;
   typedef mat::normType     NormType;
   typedef std::vector<int>  VectorTypeInt;
   typedef std::vector<real> VectorTypeReal;
   typedef generalVector     VectorType;
   typedef MatrixType        MatrixTypeWrapper;

   PurificationGeneral()
   {
      initialized_flag         = false;
      puri_is_prepared_flag    = false;
      computed_spectrum_bounds = false;

      // eigenvectors staff
      number_of_eigenvalues = 0;
      compute_eigenvectors_in_each_iteration = false;
      eigVecHOMO  = NULL;
      eigVecLUMO  = NULL;
   }

   /** Set imporatant parameters for the recursive expansion.
    */
   virtual void initialize(const MatrixType&   F_,                          /**< [in] Effective Hamiltonian matrix. */
                           const IntervalType& lumo_bounds_,                /**< [in] Bounds for lumo eigenvalue. */
                           const IntervalType& homo_bounds_,                /**< [in] Bounds for homo eigenvalue. */
                           int           maxit_,                      /**< [in] Maximum number of recursive expansion iterations. */
                           real          error_sub_,                  /**< [in] Allowed error in invariant subspaces. */
                           real          error_eig_,                  /**< [in] Error in eigenvalues (used just in old stopping criterion). */
                           int           use_new_stopping_criterion_, /**< [in] Set if want to use new stopping criterion. */
                           NormType      norm_truncation,             /**< [in] Truncation norm. */
                           NormType      norm_stop_crit,              /**< [in] Stopping criterion norm. */
                           int           nocc_                        /**< [in] Number of occupied orbitals. */
                           );


   /** Check is function initialize() is already called.
    */
   virtual bool is_initialized() const { return initialized_flag; }

   /** Check is function prepare_to_purification() is already called.
    */
   virtual bool puri_is_prepared() const { return puri_is_prepared_flag; }

   /** Start recursive expansion.
    */
   virtual void PurificationStart();

   /** Prepare data for recursive expansion. */
   virtual void prepare_to_purification();

   /** Run recursive expansion. */
   virtual void purification_process();

   /** Estimate eigenvalues near homo-lumo gap. */
   virtual void eigenvalue_bounds_estimation();

   virtual ~PurificationGeneral() {}

   /** Clear all matrices in the class.
    *  Needed to be called if Chunks and Tasks are used,
    *  since we need to delete all ChunkIDs
    *  before exiting the program. */
   virtual void clear() { X.clear(); Xsq.clear(); }


   /** Set spectrum bounds.
    *
    *  Used if we know spectrum bounds or want to compute them outside
    *  the class. */
   void set_spectrum_bounds(real eigmin, real eigmax);

   /** Get spectrum bounds.
    *
    *  Return spectrum bounds for the matrix F. */
   void get_spectrum_bounds(real& eigmin, real& eigmax);

   int get_exact_number_of_puri_iterations();
   int get_est_number_of_puri_iterations();


   virtual real total_subspace_error(int it);

   /** Get machine epsilon. */
   static real get_epsilon()
   { return mat::getMachineEpsilon<real>(); }

   /** Get largest number. */
   static real get_max_double()
   { return std::numeric_limits<real>::max(); }

   /** Get smallest number. */
   static real get_min_double()
   { return std::numeric_limits<real>::min(); }

   /** Create MATLAB .m file which plots the idempotency error in each recursive expansion iteration.  */
   void gen_matlab_file_norm_diff(const char *filename) const;
   /** Create MATLAB .m file which plots the actual introduced error after truncation of the matrix X_i in each recursive expansion iteration.  */
   void gen_matlab_file_threshold(const char *filename) const;
   /** Create MATLAB .m file which plots the number of non-zero elements in matrices X_i and X_i^2 in each recursive expansion iteration.  */
   void gen_matlab_file_nnz(const char *filename) const;
   /** Create MATLAB .m file which plots the homo and lumo bounds in each recursive expansion iteration.  */
   void gen_matlab_file_eigs(const char *filename) const;
   /** Create MATLAB .m file which creates a bar plot presenting time spent on various parts of the iteration (such as matrix square and computation of eigenvectors) in each recursive expansion iteration.  */
   void gen_matlab_file_time(const char *filename) const;
   /** Create MATLAB .m file which plots a condition number of a problem of computing the density matrix in each recursive expansion iteration. The condition number is equal to inverse of the homo-lumo gap approximation. */
   void gen_matlab_file_cond_num(const char *filename) const;

   /** Create PYTHON .py file which plots number of non-zero elements in matrices X_i and X_i^2 in each recursive expansion iteration.  */
   void gen_python_file_nnz(const char *filename) const;



   /** Set parameters for computing eigenvectors. */
   void set_eigenvectors_params(string        eigenvectors_method_,
                                string        eigenvectors_iterative_method_,
                                real          eigensolver_accuracy_,
                                int           eigensolver_maxiter_,
                                int           scf_step_,
                                int           use_prev_vector_as_initial_guess_,
                                int           try_eigv_on_next_iteration_if_fail_,
                                VectorType    *eigVecLUMO_,
                                VectorType    *eigVecHOMO_
                                );

   void set_compute_eigenvectors_in_each_iteration() { compute_eigenvectors_in_each_iteration = true; }
   void unset_compute_eigenvectors_in_each_iteration() { compute_eigenvectors_in_each_iteration = false; }


   void compute_eigenvectors_without_diagonalization_on_F(const MatrixType& F, int eigensolver_maxiter_for_F); // for testing


   PuriInfo info;  /**< Fill in during purification with useful information. */

   MatrixType X;   /**< Matrix X. */

   MatrixType Xsq; /**< Matrix X^2. */


protected:

   MatrixType F;              /**< Matrix F. Needed for computation of eigenvectors.*/

   MatrixType X_homo, X_lumo; /**< Save matrix Xi in certain iterations for computing eigenvectors (projection method). */

   void save_matrix_now(string str);

   /** Compute spectrum bounds. */
   void compute_spectrum_bounds();

   /** Get matrix X0 by mapping spectrum of F into [0,1] in reverse
    *  order.*/
   void compute_X();

   /** Get eigenvalue bounds for X0. */
   void map_bounds_to_0_1();

   /** Check stopping criterion (obsolete).
    *
    *  Use stopping criterion based on user defined threshold
    *  values. */
   virtual void check_standard_stopping_criterion(const real XmX2_norm, int& stop);

   /** Check stopping criterion.
    *
    *  The new stopping criterion based on the order of convergence is
    *  used, see article "Parameterless Stopping Criteria for Recursive Density
    *  Matrix Expansions", J. Chem. Theory Comput., 2016, 12 (12), pp 5788–5802
    *  DOI: 10.1021/acs.jctc.6b00626 */
   virtual void check_new_stopping_criterion(const int it, const real XmX2_norm_it, const real XmX2_norm_itm2, const real XmX2_trace,
                                             int& stop, real& estim_order);

   /** Choose stopping criterion and check it.  */
   virtual void stopping_criterion(IterationInfo& iter_info, int& stop, real& estim_order);

   int get_int_eig_iter_method(string eigenvectors_iterative_method);
   int get_int_eig_method(string eigenvectors_method);

   /** Compute HOMO and LUMO eigenvalues and eigenvectors of the matrix F.
   *
   * The method uses the polynomial constructed during the recursive expansion, 
   * so no additional matrix multiplications are required. See article 
   * J. Chem. Theory Comput., Just Accepted Manuscript, 
   * DOI: 10.1021/acs.jctc.7b00968
   */
   void compute_eigenvectors_without_diagonalization(int it, IterationInfo& iter_info);
   void compute_eigenvectors_without_diagonalization_last_iter_proj();

   void compute_eigenvector(MatrixType const& M, VectorType *eigVecHOMOorLUMO, int it, bool is_homo);


   /** Get nnz of X in %. */
   double get_nnz_X(size_t& nnzX /**< Number of nz of X*/)
   { nnzX = X.nnz(); return (double)(((double)nnzX) / ((double)X.get_ncols() * X.get_nrows()) * 100); }

   /** Get nnz of X in %. */
   double get_nnz_X()
   { return (double)(((double)X.nnz()) / ((double)X.get_ncols() * X.get_nrows()) * 100); }

   /** Get nnz of X^2 in %. */
   double get_nnz_Xsq(size_t& nnzXsq /**< Number of nz of X^2*/)
   { nnzXsq = Xsq.nnz(); return (double)(((double)nnzXsq) / ((double)Xsq.get_ncols() * Xsq.get_nrows()) * 100); }

   /** Get nnz of X^2 in %. */
   double get_nnz_Xsq()
   { return (double)(((double)Xsq.nnz()) / ((double)Xsq.get_ncols() * Xsq.get_nrows()) * 100); }

   /** Get homo and lumo bounds from traces and norms of Xi-Xi^2.
    *
    * Used at the end of the recursive expansion.  See article SIAM
    * J. Sci. Comput., 36(2), B147–B170. */
   void estimate_homo_lumo(const VectorTypeReal& XmX2_norm_mixed,
                           const VectorTypeReal& XmX2_norm_frob,
                           const VectorTypeReal& XmX2_trace);


   void get_frob_norm_est(const VectorTypeReal&    XmX2_norm_frob,
                          const std::vector<real>& h_in,
                          const std::vector<real>& l_in,
                          VectorTypeReal&          YmY2_norm_frob_est);



   void get_eigenvalue_estimates(const VectorTypeReal& XmX2_norm_mixed,
                                 const VectorTypeReal& XmX2_norm_frob,
                                 const VectorTypeReal& XmX2_trace);


   void propagate_values_in_each_iter(real value_unocc, real value_occ,
                                      VectorTypeReal& out_unocc,
                                      VectorTypeReal& out_occ,
                                      int nmax);


   void determine_iteration_for_eigenvectors();

   void get_iterations_for_lumo_and_homo(int& chosen_iter_lumo,
                                         int& chosen_iter_homo);

   void check_eigenvectors_at_the_end();
   void discard_homo_eigenvector();
   void discard_lumo_eigenvector();

   void output_norms_and_traces(IterationInfo& iter_info) const;

   void get_interval_with_lambda(real& eigVal, VectorType& eigVec, bool& is_homo, bool& is_lumo, const int iter);
   void get_eigenvalue_of_F_from_eigv_of_Xi(real& eigVal, const VectorType& eigVec);

   void save_eigenvectors_to_file(bool is_homo, bool is_lumo, int it);

   void set_truncation_parameters();

   void find_truncation_thresh_every_iter();
   void find_shifts_every_iter();
   void find_eig_gaps_every_iter();


   void writeToTmpFile(MatrixType& A) const;
   void readFromTmpFile(MatrixType& A) const;

   /*
    * virtual void update_bounds(const real value);
    * real compute_chemical_potential(const int it);
    * real get_lower_bound_in_estim_homo_lumo(const int it);
    */

   virtual void set_init_params() = 0;
   virtual void truncate_matrix(real& thresh, int it) = 0;
   virtual void estimate_number_of_iterations(int& estim_num_iter) = 0;
   virtual void purify_X(const int it)      = 0;
   virtual void purify_bounds(const int it) = 0;
   virtual void save_other_iter_info(IterationInfo& iter_info, int it) = 0;
   virtual void apply_inverse_poly_vector(const int it, VectorTypeReal& bounds_from_it) = 0;
   virtual void return_constant_C(const int it, real& Cval) = 0;

   /*virtual real apply_inverse_poly(const int it, real x) = 0;*/
   virtual real apply_poly(const int it, real x) = 0;
   virtual void apply_poly_to_eigs(const int it, real& homo, real& lumo) = 0;
   virtual real compute_derivative(const int it, real x, real& DDf)      = 0;



   /* PARAMETERS */

   bool initialized_flag;
   bool puri_is_prepared_flag;

   int use_new_stopping_criterion;    /**< True for new stopping criterion */
   int additional_iterations;         /**< Do a few more iterations after convergence */

   int maxit;                         /**< Maximum number of iterations */
   int check_stopping_criterion_iter; /**< Iteration when to start to check stopping criterion. */

   int nocc;                          /**<  Number of occupied orbitals */



   NormType normPuriTrunc;  /**< Norm used for the truncation of matrices.
                             * Can be mat::frobNorm, mat::mixedNorm, or mat::euclNorm. */


   NormType normPuriStopCrit; /**< Norm used in the stopping criterion
                               * Can be mat::frobNorm, mat::mixedNorm, or mat::euclNorm. */


   real error_sub;      /**< Allowed error in invariant subspaces. */
   real error_eig;      /**< Error in eigenvalues (used just in old stopping criterion). */

   real error_per_it;   /**< Error allowed in each iteration due to truncation. */

   real constant_C;     /**< Asymptotic constant C needed for the new stopping criterion. */

   real gammaStopEstim; /**< Used on the stopping criterion for
                         * estimation of eigenvalues from
                         * purification */


   VectorTypeInt VecPoly; /**< Polynomials computed in the function
                           * estimated_number_of_iterations()
                           * VecPoly[i] = 1 if we use X=X^2
                           * VecPoly[i] = 0 if we use X=2X-X^2
                           * (or their accelerated versions) */
   VectorTypeReal VecGap; /**< Gap computed using inner homo and lumo bounds on each iteration. */


   VectorTypeReal ITER_ERROR_VEC;                             /**< (Eigenvectors) Maximum error introduced in each iteration. */
   VectorTypeReal SIGMA_HOMO_VEC, SIGMA_LUMO_VEC;             /**< (Eigenvectors) Approximation of shifts in each iteration. */
   VectorTypeReal EIG_ABS_GAP_LUMO_VEC, EIG_ABS_GAP_HOMO_VEC; /**< (Eigenvectors) Absolute and relative gap in filter for lumo eigenvalue. */
   VectorTypeReal EIG_REL_GAP_LUMO_VEC, EIG_REL_GAP_HOMO_VEC; /**< (Eigenvectors) Absolute and relative gap in filter for homo eigenvalue. */


   /*EIGENVECTORS STAFF*/

   int number_of_eigenvalues;

   IntervalType homo_bounds;    /**< (1-homo) bounds for Xi in iteration i */
   IntervalType lumo_bounds;    /**< Lumo bounds for Xi in iteration i */

   IntervalType homo_bounds_X0; /**< Initial lumo bounds for X */
   IntervalType lumo_bounds_X0; /**< Initial lumo bounds for X */

   IntervalType homo_bounds_F;  /**< Initial lumo bounds for F */
   IntervalType lumo_bounds_F;  /**< Initial homo bounds for F */

   IntervalType homo_bounds_F_new;
   IntervalType lumo_bounds_F_new;


   IntervalType spectrum_bounds; /**< Outer bounds for the whole spectrum of F/Xi. */
   bool computed_spectrum_bounds;


   /* Eigenvectors */
   int eigenvectors_method;              /**< Chosen method for computing eigenvectors. */
   int eigenvectors_iterative_method;    /**< Chosen eigensolver. */

   // accuracy of the eigenvalue problem solver
   // when residual is less then this value,stop
   real eigensolver_accuracy;
   // maximum number of iterations for eigensolver
   int eigensolver_maxiter;

   string eigenvectors_method_str;
   string eigenvectors_iterative_method_str;

   int use_prev_vector_as_initial_guess;

   bool compute_eigenvectors_in_this_SCF_cycle;
   bool try_eigv_on_next_iteration_if_fail;

   VectorType *eigVecLUMO;
   VectorType *eigVecHOMO;

   VectorType eigVecLUMORef;
   VectorType eigVecHOMORef;


   real eigValLUMO;
   real eigValHOMO;

   int iter_for_homo;
   int iter_for_lumo;

   VectorTypeInt good_iterations_homo;        /**< Iterations where homo eigenvector can be computed. */
   VectorTypeInt good_iterations_lumo;        /**< Iterations where homo eigenvector can be computed. */

   VectorTypeInt really_good_iterations_homo; /**< Iterations where homo eigenvector is actually computed. */
   VectorTypeInt really_good_iterations_lumo; /**< Iterations where lumo eigenvector is actually computed.  */

   int scf_step;

   bool compute_eigenvectors_in_each_iteration; /**< Compute homo and lumo eigenpairs in every iteration
                                                 * and save eigenvectors in txt files.
                                                 *
                                                 * NOTE: if we are computing eigenvector in every iteration, we will
                                                 * save eigenvector computed in the last iteration, and it is most
                                                 * probably the wrong eigenvector.
                                                 */
};

/**************************************/


/******************* INIT *******************/


template<typename MatrixType>
void PurificationGeneral<MatrixType>::initialize(const MatrixType&   F_,
                                                 const IntervalType& lumo_bounds_,
                                                 const IntervalType& homo_bounds_,
                                                 int           maxit_,
                                                 real          error_sub_,
                                                 real          error_eig_,
                                                 int           use_new_stopping_criterion_,
                                                 NormType      norm_truncation_,
                                                 NormType      norm_stop_crit_,
                                                 int           nocc_
                                                 )
{
   X     = F_;
   maxit = maxit_;
   assert(maxit >= 1);
   error_sub = error_sub_;
   error_eig = error_eig_;
   use_new_stopping_criterion = use_new_stopping_criterion_;
   normPuriTrunc    = norm_truncation_;
   normPuriStopCrit = norm_stop_crit_;
   nocc             = nocc_;

   initialized_flag = true;

   // save bounds for the matrix F
   lumo_bounds_F = lumo_bounds_;
   homo_bounds_F = homo_bounds_;

   /* Use this function to enable printf output of the purification
    * work if you want to run just the purification, not the whole scf calculations */
#ifdef ENABLE_PRINTF_OUTPUT
   enable_printf_output();
#endif

   size_t nnzX;
   double nnzXp = get_nnz_X(nnzX);
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Creating purification object: N = %d"
                                               " , nocc = %d , NNZ = %lu  <-> %.5lf %%",
             X.get_nrows(), nocc, nnzX, nnzXp);


   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Chosen norm for the truncation: ");
   switch (normPuriTrunc)
   {
   case mat::mixedNorm:
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "mixed");
      break;

   case mat::euclNorm:
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "eucl");
      break;

   case mat::frobNorm:
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "frob");
      break;

   default:
      throw std::runtime_error("Unknown norm in PurificationGeneral");
   }


#ifdef USE_CHUNKS_AND_TASKS
   if ((normPuriTrunc == mat::mixedNorm) || (normPuriTrunc == mat::euclNorm))
   {
      normPuriTrunc = mat::frobNorm;
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Change norm for the truncation to Frobenius.");
   }
#endif

   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Chosen norm for the stopping criterion: ");
   switch (normPuriStopCrit)
   {
   case mat::mixedNorm:
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "mixed");
      break;

   case mat::euclNorm:
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "eucl");
      break;

   case mat::frobNorm:
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "frob");
      break;

   default:
      throw std::runtime_error("Unknown norm in PurificationGeneral");
   }


#ifdef USE_CHUNKS_AND_TASKS
   if ((normPuriStopCrit == mat::mixedNorm) || (normPuriStopCrit == mat::euclNorm))
   {
      normPuriStopCrit = mat::frobNorm;
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Change norm the stopping criterion to Frobenius.");
   }
#endif

   if (this->use_new_stopping_criterion == 1)
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Chosen the NEW stopping criterion.");
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Allowed error in subspace %e", (double)error_sub);
   }
   else
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Chosen the OLD stopping criterion.");
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Allowed error in subspace %e , in eigenvalues %e", (double)error_sub, (double)error_eig);
   }

   check_stopping_criterion_iter          = -1;    // will be set in function estimate_number_of_iterations()
   compute_eigenvectors_in_this_SCF_cycle = false; // can be set to true in prepare_to_purification()

   set_init_params();

   additional_iterations   = NUM_ADDITIONAL_ITERATIONS;
   info.stopping_criterion = use_new_stopping_criterion; // 1 if new, 0 if not
   info.error_subspace     = error_sub;

   info.debug_output = 0;
#ifdef DEBUG_PURI_OUTPUT
   info.debug_output = 1;
#endif
}


/*******************************************************/


template<typename MatrixType>
void PurificationGeneral<MatrixType>::set_eigenvectors_params(string        eigenvectors_method_,
  string        eigenvectors_iterative_method_,
  real          eigensolver_accuracy_,
  int           eigensolver_maxiter_,
  int           scf_step_,
  int           use_prev_vector_as_initial_guess_,
  int           try_eigv_on_next_iteration_if_fail_,
  VectorType    *eigVecLUMO_,
  VectorType    *eigVecHOMO_
)
{
   // before this was an input parameter
   number_of_eigenvalues = 1;
   if (number_of_eigenvalues > 1)
   {
      throw std::runtime_error("Error in set_eigenvectors_params() : cannot compute more than 1 eigenpair.");
   }

   // can be NULL, empty or containing eigenvector from the previous cycle
   eigVecLUMO  = eigVecLUMO_;
   eigVecHOMO  = eigVecHOMO_;

   eigensolver_accuracy = eigensolver_accuracy_;
   eigensolver_maxiter  = eigensolver_maxiter_;
   assert(eigensolver_maxiter >= 1);
   scf_step = scf_step_;
   eigenvectors_method_str           = eigenvectors_method_;
   eigenvectors_iterative_method_str = eigenvectors_iterative_method_;
   eigenvectors_method                = get_int_eig_method(eigenvectors_method_);
   eigenvectors_iterative_method      = get_int_eig_iter_method(eigenvectors_iterative_method_);
   use_prev_vector_as_initial_guess   = use_prev_vector_as_initial_guess_;
   try_eigv_on_next_iteration_if_fail = try_eigv_on_next_iteration_if_fail_;

   info.lumo_eigenvector_is_computed = false;
   info.homo_eigenvector_is_computed = false;

   iter_for_homo = -1;
   iter_for_lumo = -1;

   // no given method for computing eigenvectors
   if (((eigVecLUMO != NULL) || (eigVecHOMO != NULL)) && (eigenvectors_method == EIG_EMPTY))
   {

      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "No given method for computing eigenvectors."
                                                  "Eigenvectors will not be computed in this SCF cycle. Set eigenvectors to NULL.");

      delete eigVecLUMO; // not NULL here
      delete eigVecHOMO; // not NULL here
      eigVecLUMO = NULL;
      eigVecHOMO = NULL;
   }
   else
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Chosen method to compute eigenvectors: %s", eigenvectors_method_str.c_str());
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Chosen iterative method to compute eigenvectors: %s", eigenvectors_iterative_method_str.c_str());
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Chosen eigensolver accuracy: %g", (double)eigensolver_accuracy);
   }

   // reuse eigenvector computed in the previous SCF cycle as an initial guess
   if (((eigVecLUMO != NULL) || (eigVecHOMO != NULL)) && use_prev_vector_as_initial_guess)
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Use eigenvectors from the previous SCF cycle as an initial guess for this SCF cycle");
   }


#ifndef USE_CHUNKS_AND_TASKS
   if (compute_eigenvectors_in_each_iteration)
   {
      /* if we compute eigenvectors in every iteration, we save initial guess in the separate vector eigVecLUMORef*/
      if ((eigVecLUMO != NULL) && use_prev_vector_as_initial_guess)
      {
         mat::SizesAndBlocks cols;
         if (X.is_empty())
         {
            throw std::runtime_error("Error in set_eigenvectors_params() : cannot save initial guess for LUMO!");
         }
         X.getCols(cols);
         eigVecLUMORef.resetSizesAndBlocks(cols);
         eigVecLUMORef = *eigVecLUMO;
      }

      if ((eigVecHOMO != NULL) && use_prev_vector_as_initial_guess)
      {
         mat::SizesAndBlocks cols;
         if (X.is_empty())
         {
            throw std::runtime_error("Error in set_eigenvectors_params() : cannot save initial guess for HOMO!");
         }
         X.getCols(cols);
         eigVecHOMORef.resetSizesAndBlocks(cols);
         eigVecHOMORef = *eigVecHOMO;
      }
   }
#endif
}


template<typename MatrixType>
int PurificationGeneral<MatrixType>::get_int_eig_method(string eigenvectors_method)
{
   if (eigenvectors_method == "square")
   {
      return EIG_SQUARE_INT;
   }
   if (eigenvectors_method == "projection")
   {
      return EIG_PROJECTION_INT;
   }
   if (eigenvectors_method == "")
   {
      return EIG_EMPTY;
   }
   throw std::runtime_error("Error in get_int_eig_method(): unknown method to compute eigenvectors");
}


template<typename MatrixType>
int PurificationGeneral<MatrixType>::get_int_eig_iter_method(string eigenvectors_iterative_method)
{
   if (eigenvectors_iterative_method == "power")
   {
      return EIG_POWER_INT;
   }
   if (eigenvectors_iterative_method == "lanczos")
   {
      return EIG_LANCZOS_INT;
   }
   if (eigenvectors_iterative_method == "")
   {
      return EIG_EMPTY;
   }
   throw std::runtime_error("Error in get_int_eig_iter_method(): unknown iterative method to compute eigenvectors");
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::output_norms_and_traces(IterationInfo& iter_info) const
{
   real XmX2_fro_norm   = iter_info.XmX2_fro_norm;
   real XmX2_eucl       = iter_info.XmX2_eucl;
   real XmX2_mixed_norm = iter_info.XmX2_mixed_norm;
   real XmX2_trace      = iter_info.XmX2_trace;
   int  it = iter_info.it;

   assert(it >= 0);

#ifndef USE_CHUNKS_AND_TASKS
   if (normPuriStopCrit == mat::euclNorm)
   {
      assert(XmX2_eucl >= 0);
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "||X-X^2||_2 = %e", (double)XmX2_eucl);
   }

   if (normPuriStopCrit == mat::mixedNorm)
   {
      assert(XmX2_fro_norm >= 0);
      assert(XmX2_mixed_norm >= 0);
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "||X-X^2||_F = %e , ||X-X^2||_mixed = %e", (double)XmX2_fro_norm, (double)XmX2_mixed_norm);
   }
#endif

   if (normPuriStopCrit == mat::frobNorm)
   {
      assert(XmX2_fro_norm >= 0);
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "||X-X^2||_F = %e", (double)XmX2_fro_norm);
   }

   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "trace(X-X^2) = %e", (double)XmX2_trace);
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "trace(X) = %e", (double)X.trace());
}


/****************** PURIFICATION_START ********************/


template<typename MatrixType>
void PurificationGeneral<MatrixType>::PurificationStart()
{
   Util::TimeMeter total_time;   // measure total time of the purification process

   prepare_to_purification();

#ifdef DEBUG_PURI_OUTPUT
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Starting recursive expansion");
#endif
   purification_process();

   eigenvalue_bounds_estimation();



   /* COMPUTE EIGENVECTORS WITH PROJECTION METHOD */
   if (info.converged == 1)
   {
      if (compute_eigenvectors_in_this_SCF_cycle && (eigenvectors_method == EIG_PROJECTION_INT))
      {
         compute_eigenvectors_without_diagonalization_last_iter_proj();
      }
   }
   else
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Cannot compute eigenvectors using projection method since the purification did not converge");
   }

   /* CHECK IF EIGENVECTORS ARE CORRECT - EITHER COMPUTED WITH SQUARE METHOD OR PROJECTION */
   check_eigenvectors_at_the_end();

   total_time.print(LOG_AREA_DENSFROMF, "Recursive expansion");
   double total_time_stop = total_time.get_elapsed_wall_seconds();
   info.total_time = total_time_stop;
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::prepare_to_purification()
{
   if (!is_initialized())
   {
      throw std::runtime_error("Error in prepare_to_purification() : function is called for an uninitialized class.");
   }

   if (!computed_spectrum_bounds)
   {
      Util::TimeMeter total_time_spectrum_bounds;
      compute_spectrum_bounds();
      total_time_spectrum_bounds.print(LOG_AREA_DENSFROMF, "compute_spectrum_bounds");
      double total_time_spectrum_bounds_stop = total_time_spectrum_bounds.get_elapsed_wall_seconds();
      info.time_spectrum_bounds = total_time_spectrum_bounds_stop;
   }

   info.set_spectrum_bounds(spectrum_bounds.low(), spectrum_bounds.upp());
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Spectrum of F: \t [ %.12lf , %.12lf ]", (double)spectrum_bounds.low(), (double)spectrum_bounds.upp());

   map_bounds_to_0_1();
   set_truncation_parameters();

   if ((eigVecLUMO != NULL) || (eigVecHOMO != NULL))
   {
      // check if we have non-overlapping homo and lumo bounds
      if (1 - homo_bounds.upp() - lumo_bounds.upp() >= TOL_OVERLAPPING_BOUNDS)
      {
         compute_eigenvectors_in_this_SCF_cycle      = true;
         info.compute_eigenvectors_in_this_SCF_cycle = compute_eigenvectors_in_this_SCF_cycle;
         F = X;    // Save original matrix F, needed for computation of the Rayleigh quotient.
         // Quotient is needed to compute eigenvalue corresponding to the eigenvector
         // and see which eigenvalue it is, homo or lumo.
         writeToTmpFile(F);
         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "calling determine_iteration_for_eigenvectors()");
         determine_iteration_for_eigenvectors();    // on which iterations we should compute eigenvectors
      }
      else
      {
         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Homo and lumo inner bounds are too close (< %g), "
                                                     "homo and lumo eigenvectors will not be computed", (double)TOL_OVERLAPPING_BOUNDS);
      }
   }

   compute_X(); // F->X, put eigenvalues to the [0,1]

   puri_is_prepared_flag = true;
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::purification_process()
{
   if (!puri_is_prepared())
   {
      throw std::runtime_error("Error in purification_process() : "
                               "first expect a call for function prepare_to_purification().");
   }

   int  it;
   int  stop = 0;
   real thresh;
   real Xsquare_time_stop = -1, total_time_stop = -1, trunc_time_stop = -1, purify_time_stop = -1;
   real frob_diff_time_stop = -1, eucl_diff_time_stop = -1, trace_diff_time_stop = -1, mixed_diff_time_stop = -1, stopping_criterion_time_stop = -1;
   int  maxit_tmp = maxit;
   real estim_order = -1;
   real XmX2_trace = -1;
   real XmX2_fro_norm = -1;
   real XmX2_mixed_norm = -1;
   real XmX2_eucl = -1;

   int already_stop_before = 0; // flag for checking stopping criterion, needed in case additional_iterations != 0

   info.Iterations.clear();
   info.Iterations.reserve(100);
   //std::cout  <<  "Cap: " << this->info.Iterations.capacity() << std::endl;

   IterationInfo iter_info; // 0-th iterations

   // 0 iteration
   it = 0;


   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "BEFORE ITERATIONS:");



#ifdef SAVE_MATRIX_IN_EACH_ITERATION
   {
      ostringstream str;
      str << it;
      save_matrix_now(str.str());
   }
#endif

   Util::TimeMeter total_time; // for this iteration

#ifdef PURI_OUTPUT_NNZ
   double nnzX = get_nnz_X();
#endif

   // truncate matrix
   Util::TimeMeter trunc_time;
   truncate_matrix(thresh, it);
   trunc_time.print(LOG_AREA_DENSFROMF, "truncate_matrix");
   trunc_time_stop = trunc_time.get_elapsed_wall_seconds();

#ifdef PURI_OUTPUT_NNZ
   double nnzXafter = get_nnz_X();
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Actual introduced error   %e , nnz before   %.2lf %% , nnz after   %.2lf %%", (double)thresh, nnzX, nnzXafter);
#else
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Actual introduced error   %e", (double)thresh);
#endif

   output_current_memory_usage(LOG_AREA_DENSFROMF, "Before X square");

   Util::TimeMeter Xsquare_time;
   Xsq = (real)1.0 * X * X;
   Xsquare_time.print(LOG_AREA_DENSFROMF, "square");
   Xsquare_time_stop = Xsquare_time.get_elapsed_wall_seconds();

#ifdef PURI_OUTPUT_NNZ
   size_t nnzXsq;
   double nnzXsqp = get_nnz_Xsq(nnzXsq);
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "NNZ Xsq = %lu <-> %.2lf %%", nnzXsq, nnzXsqp);
#endif


   // compute frob norm of X-X2
   Util::TimeMeter frob_diff_time;
   XmX2_fro_norm = MatrixType::frob_diff(X, Xsq);
   frob_diff_time.print(LOG_AREA_DENSFROMF, "Frobenius norm of X-X^2");
   frob_diff_time_stop = frob_diff_time.get_elapsed_wall_seconds();

   if (normPuriStopCrit == mat::mixedNorm)
   {
#ifndef USE_CHUNKS_AND_TASKS
      Util::TimeMeter mixed_diff_time;
      XmX2_mixed_norm = MatrixType::mixed_diff(X, Xsq, mixed_acc);
      mixed_diff_time.print(LOG_AREA_DENSFROMF, "Mixed norm of X-X^2");
      mixed_diff_time_stop = mixed_diff_time.get_elapsed_wall_seconds();
#endif
   }

   // compute trace of X-X2
   Util::TimeMeter trace_diff_time;
   XmX2_trace = X.trace() - Xsq.trace();
   trace_diff_time.print(LOG_AREA_DENSFROMF, "Trace of X-X^2");
   trace_diff_time_stop = trace_diff_time.get_elapsed_wall_seconds();

   if (normPuriStopCrit == mat::euclNorm)
   {
#ifndef USE_CHUNKS_AND_TASKS
      Util::TimeMeter eucl_diff_time;
      XmX2_eucl = MatrixType::eucl_diff(X, Xsq, eucl_acc);
      eucl_diff_time.print(LOG_AREA_DENSFROMF, "Spectral norm of X-X^2");
      eucl_diff_time_stop = eucl_diff_time.get_elapsed_wall_seconds();
#endif
   }


   iter_info.it              = it;
   iter_info.XmX2_fro_norm   = XmX2_fro_norm;
   iter_info.XmX2_eucl       = XmX2_eucl;
   iter_info.XmX2_mixed_norm = XmX2_mixed_norm;
   iter_info.XmX2_trace      = XmX2_trace;
#ifdef DEBUG_PURI_OUTPUT
   output_norms_and_traces(iter_info);
#endif


   if (compute_eigenvectors_in_this_SCF_cycle)
   {
      compute_eigenvectors_without_diagonalization(it, iter_info);
   }

   ostringstream str_out;
   str_out << "Iteration " << it;
   total_time.print(LOG_AREA_DENSFROMF, str_out.str().c_str());
   total_time_stop = total_time.get_elapsed_wall_seconds();

// SAVE INFO ABOUT ITERATION
   {
      iter_info.gap                     = 1 - homo_bounds.upp() - lumo_bounds.upp(); // = VecGap[0]
      iter_info.threshold_X             = thresh;
      iter_info.Xsquare_time            = Xsquare_time_stop;
      iter_info.trunc_time              = trunc_time_stop;
      iter_info.purify_time             = 0;
      iter_info.NNZ_X                   = get_nnz_X();
      iter_info.NNZ_X2                  = get_nnz_Xsq();
      iter_info.total_time              = total_time_stop;
      iter_info.homo_bound_low          = homo_bounds.low();
      iter_info.lumo_bound_low          = lumo_bounds.low();
      iter_info.homo_bound_upp          = homo_bounds.upp();
      iter_info.lumo_bound_upp          = lumo_bounds.upp();
      stopping_criterion_time_stop      = 0; // we are not checking stopping criterion on the 1 iteration
      iter_info.stopping_criterion_time = stopping_criterion_time_stop;
      iter_info.eucl_diff_time          = eucl_diff_time_stop;
      iter_info.frob_diff_time          = frob_diff_time_stop;
      iter_info.mixed_diff_time         = mixed_diff_time_stop;
      iter_info.trace_diff_time         = trace_diff_time_stop;

      save_other_iter_info(iter_info, it);
   }
   /**************/

   info.Iterations.push_back(iter_info); // add info about 0 iteration


   output_current_memory_usage(LOG_AREA_DENSFROMF, "Before iteration 1");
   Util::TimeMeter timeMeterWriteAndReadAll;
   std::string     sizesStr = mat::FileWritable::writeAndReadAll();
   timeMeterWriteAndReadAll.print(LOG_AREA_DENSFROMF, "FileWritable::writeAndReadAll");
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, ((std::string)"writeAndReadAll sizesStr: '" + sizesStr).c_str());


   it = 1;
   while (it <= maxit_tmp)
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "ITERATION %d :", it);

      IterationInfo iter_info;    // i-th iteration

      Util::TimeMeter total_time; // for this iteration

      Util::TimeMeter purify_time;
      purify_X(it);
      purify_time.print(LOG_AREA_DENSFROMF, "purify_X");
      purify_time_stop = purify_time.get_elapsed_wall_seconds();


#ifdef SAVE_MATRIX_IN_EACH_ITERATION
      {
         ostringstream str;
         str << it;
         save_matrix_now(str.str());
      }
#endif

#ifdef PURI_OUTPUT_NNZ
      double nnzX = get_nnz_X();
#endif

      // truncate matrix
      Util::TimeMeter trunc_time;
      truncate_matrix(thresh, it);
      trunc_time.print(LOG_AREA_DENSFROMF, "truncate_matrix");
      trunc_time_stop = trunc_time.get_elapsed_wall_seconds();

#ifdef PURI_OUTPUT_NNZ
      double nnzXafter = get_nnz_X();
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Actual introduced error   %e , nnz before   %.2lf %% , nnz after   %.2lf %%", thresh, nnzX, nnzXafter);
#else
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Actual introduced error   %e", thresh);
#endif


      output_current_memory_usage(LOG_AREA_DENSFROMF, "Before X square");

      Util::TimeMeter Xsquare_time;
      Xsq = (real)1.0 * X * X;
      Xsquare_time.print(LOG_AREA_DENSFROMF, "square");
      Xsquare_time_stop = Xsquare_time.get_elapsed_wall_seconds();
#ifdef PURI_OUTPUT_NNZ
      size_t nnzXsq;
      double nnzXsqp = get_nnz_Xsq(nnzXsq);
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "NNZ Xsq = %lu <-> %.2lf %%", nnzXsq, nnzXsqp);
#endif


      // update bounds for homo and lumo using corresponding polynomial
      purify_bounds(it);

      // compute frob norm of X-X2
      Util::TimeMeter frob_diff_time;
      XmX2_fro_norm = MatrixType::frob_diff(X, Xsq);
      frob_diff_time.print(LOG_AREA_DENSFROMF, "Frobenius norm of X-X^2");
      frob_diff_time_stop = frob_diff_time.get_elapsed_wall_seconds();

      if (normPuriStopCrit == mat::mixedNorm)
      {
#ifndef USE_CHUNKS_AND_TASKS
         Util::TimeMeter mixed_diff_time;
         XmX2_mixed_norm = MatrixType::mixed_diff(X, Xsq, mixed_acc);
         //XmX2_mixed_norm = Xsq.mixed(mixed_acc);
         mixed_diff_time.print(LOG_AREA_DENSFROMF, "Mixed norm of X-X^2");
         mixed_diff_time_stop = mixed_diff_time.get_elapsed_wall_seconds();
#endif
      }

      if (normPuriStopCrit == mat::euclNorm)
      {
#ifndef USE_CHUNKS_AND_TASKS
         Util::TimeMeter eucl_diff_time;
         XmX2_eucl = MatrixType::eucl_diff(X, Xsq, eucl_acc);
         eucl_diff_time.print(LOG_AREA_DENSFROMF, "Spectral norm of X-X^2");
         eucl_diff_time_stop = eucl_diff_time.get_elapsed_wall_seconds();
#endif
      }

      // compute trace of X-X2
      Util::TimeMeter trace_diff_time;
      XmX2_trace = X.trace() - Xsq.trace();
      trace_diff_time.print(LOG_AREA_DENSFROMF, "Trace of X-X^2");
      trace_diff_time_stop = trace_diff_time.get_elapsed_wall_seconds();


      // CHECK FOR A NEGATIVE TRACE
      if (XmX2_trace < -1e10) // here is definitively some misconvergence
      {
         throw std::runtime_error("Error in purification_process() : trace of X-X^2 is negative, seems as a"
                                  " misconvergence of the recursive expansion.");
      }

      iter_info.it              = it;
      iter_info.XmX2_fro_norm   = XmX2_fro_norm;
      iter_info.XmX2_eucl       = XmX2_eucl;
      iter_info.XmX2_mixed_norm = XmX2_mixed_norm;
      iter_info.XmX2_trace      = XmX2_trace;
#ifdef DEBUG_PURI_OUTPUT
         output_norms_and_traces(iter_info);
#endif


// SAVE INFO ABOUT ITERATION
      {
         iter_info.threshold_X  = thresh;
         iter_info.Xsquare_time = Xsquare_time_stop;
         iter_info.trunc_time   = trunc_time_stop;
         iter_info.purify_time  = purify_time_stop;

         iter_info.eucl_diff_time  = eucl_diff_time_stop;
         iter_info.frob_diff_time  = frob_diff_time_stop;
         iter_info.mixed_diff_time = mixed_diff_time_stop;
         iter_info.trace_diff_time = trace_diff_time_stop;
      }


      /* COMPUTE EIGENVECTORS WITH FOLDED SPECTRUM METHOD (SHIFT_AND_SQUARE) */
      if (compute_eigenvectors_in_this_SCF_cycle)
      {
         compute_eigenvectors_without_diagonalization(it, iter_info);
      }


      // check stopping criterion (call function on every iteration
      // larger or equal to check_stopping_criterion_iter)
      if (it >= check_stopping_criterion_iter)
      {
         Util::TimeMeter stopping_criterion_time;
         stopping_criterion(iter_info, stop, estim_order);
         stopping_criterion_time.print(LOG_AREA_DENSFROMF, "stopping_criterion");
         stopping_criterion_time_stop      = stopping_criterion_time.get_elapsed_wall_seconds();
         iter_info.stopping_criterion_time = stopping_criterion_time_stop;
      }
      else
      {
         stop        = 0;
         estim_order = -1;
      }

      // if we reach stopping iteration or maximum allowed number or iterations
      // and we are not already stop (in case we have additional_iterations != 0)
      if ((already_stop_before == 0) && ((stop == 1) || (it == maxit)))
      {
         if (stop == 1)
         {
            do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "PURIFICATION CONVERGED after %d iterations", it);
            info.converged = 1;
         }
         else  // if it == maxit
         {
            do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "NOT CONVERGED. Reached the maximum number of iterations %d", maxit);
            info.converged = 0;
         }

         assert(maxit_tmp <= (int)VecPoly.size());
         maxit_tmp           = it + additional_iterations;
         already_stop_before = 1;
      }

      ostringstream str_out;
      str_out << "Iteration " << it;
      total_time.print(LOG_AREA_DENSFROMF, str_out.str().c_str());
      total_time_stop = total_time.get_elapsed_wall_seconds();



      /******************/


      iter_info.NNZ_X  = get_nnz_X();
      iter_info.NNZ_X2 = get_nnz_Xsq();

      iter_info.homo_bound_low = homo_bounds.low();
      iter_info.homo_bound_upp = homo_bounds.upp();
      iter_info.lumo_bound_low = lumo_bounds.low();
      iter_info.lumo_bound_upp = lumo_bounds.upp();

      iter_info.total_time = total_time_stop;
      iter_info.constantC  = constant_C;
      if (use_new_stopping_criterion)
      {
         iter_info.order = estim_order;
      }

      save_other_iter_info(iter_info, it);

      /*******************/

      info.Iterations.push_back(iter_info); // add info about i-th iteration to the info

      it++;
   }  //while


   // output number of non-zeros in the obtained density matrix
   double nnzD = get_nnz_X();
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Number of non-zeros in D is  %.2lf %%", nnzD);


   output_current_memory_usage(LOG_AREA_DENSFROMF, "After the last iteration");
   Util::TimeMeter timeMeterWriteAndReadAll_end;
   sizesStr = mat::FileWritable::writeAndReadAll();
   timeMeterWriteAndReadAll_end.print(LOG_AREA_DENSFROMF, "FileWritable::writeAndReadAll");
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, ((std::string)"writeAndReadAll sizesStr: '" + sizesStr).c_str());


   info.total_it = maxit_tmp;
   info.additional_iterations = additional_iterations;

   real acc_err_sub = total_subspace_error(maxit_tmp - additional_iterations);
#ifdef DEBUG_PURI_OUTPUT
   if (acc_err_sub != -1)
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "TOTAL accumulated subspace error is %e", acc_err_sub);
   }
#endif
   info.accumulated_error_subspace = acc_err_sub;
}


/* EIGENVALUE BOUND ESTIMATION  */

template<typename MatrixType>
void PurificationGeneral<MatrixType>::eigenvalue_bounds_estimation()
{
   if (info.converged == 1)
   {
      // estimate eigenvalues of the matrix F
      VectorTypeReal norms_mixed, norms_frob, traces;
      // use mixed norm instead of the Frobenius if it is possible
      if (normPuriStopCrit == mat::mixedNorm)
      {
         info.get_vec_mixed_norms(norms_mixed);
      }
      info.get_vec_frob_norms(norms_frob);
      info.get_vec_traces(traces);
      get_eigenvalue_estimates(norms_mixed, norms_frob, traces);


      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Estimated bounds for the eigenvalues for the Fock matrix:");
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "LUMO: [ %.12lf , %.12lf ]", (double)lumo_bounds_F_new.low(), (double)lumo_bounds_F_new.upp());
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "HOMO: [ %.12lf , %.12lf ]", (double)homo_bounds_F_new.low(), (double)homo_bounds_F_new.upp());
   }
   else
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Cannot estimate eigenvalues since the purification did not converge");
   }


   info.homo_estim_low_F = homo_bounds_F_new.low();
   info.homo_estim_upp_F = homo_bounds_F_new.upp();
   info.lumo_estim_low_F = lumo_bounds_F_new.low();
   info.lumo_estim_upp_F = lumo_bounds_F_new.upp();
}


/******************************************************************************************************************************/


template<typename MatrixType>
void PurificationGeneral<MatrixType>::discard_homo_eigenvector()
{
   if ((eigVecHOMO == NULL) || eigVecHOMO->is_empty())
   {
      return;
   }
   info.eigValHOMO = 0;
   info.homo_eigenvector_is_computed = false;

   /* NOTE: clear() gives a zero vector, which may not be empty (i.e. still save some data about the structure)!
    * Thus we use clear_structure() to get an empty vector. */
   eigVecHOMO->clear_structure();
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::discard_lumo_eigenvector()
{
   if ((eigVecLUMO == NULL) || eigVecLUMO->is_empty())
   {
      return;
   }
   info.eigValLUMO = 0;
   info.lumo_eigenvector_is_computed = false;
   eigVecLUMO->clear_structure();
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::check_eigenvectors_at_the_end()
{
   if (info.converged != 1)
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Discard all computed eigenvectors since the purification did not converge");
      discard_homo_eigenvector();
      discard_lumo_eigenvector();
      return;
   }



   if (compute_eigenvectors_in_this_SCF_cycle)
   {
      int ITER_DIFF = 2;
      // check if we passed iter_for_homo iteration
      if (!eigVecHOMO->is_empty() && (info.total_it < iter_for_homo))
      {
         do_output(LOG_CAT_WARNING, LOG_AREA_DENSFROMF, "HOMO eigenvector was not computed. Iteration for homo: %d, total number of iterations: %d",
                   iter_for_homo, info.total_it);
         discard_homo_eigenvector();
      }
      else
      {
         // if yes, was it one of the last iterations?
         if (!eigVecHOMO->is_empty() && (info.total_it - iter_for_homo < ITER_DIFF) && info.homo_eigenvector_is_computed)
         {
            do_output(LOG_CAT_WARNING, LOG_AREA_DENSFROMF, "HOMO eigenvector was computed in one of the last recursive expansion iterations (%d of total %d). "
                                                           "Eigenvalues of the matrix X in this iteration probably already reached the level of numerical errors, "
                                                           "thus result may not be accurate!", iter_for_homo, info.total_it);
         }
      }

      // check if we passed iter_for_lumo iteration
      if (!eigVecLUMO->is_empty() && (info.total_it < iter_for_lumo))
      {
         do_output(LOG_CAT_WARNING, LOG_AREA_DENSFROMF, "LUMO eigenvector was not computed. Iteration for lumo: %d, total number of iterations: %d",
                   iter_for_lumo, info.total_it);
         discard_lumo_eigenvector();
      }
      else
      {
         // if yes, was it one of the last iterations?
         if (!eigVecLUMO->is_empty() && (info.total_it - iter_for_lumo < ITER_DIFF) && info.lumo_eigenvector_is_computed)
         {
            do_output(LOG_CAT_WARNING, LOG_AREA_DENSFROMF, "LUMO eigenvector was computed in one of the last recursive expansion iterations (%d of total %d). "
                                                           "Eigenvalues of the matrix X in this iteration probably already reached the level of numerical errors, "
                                                           "thus result may not be accurate!", iter_for_lumo, info.total_it);
         }
      }



      real low_lumo_F_bound  = info.lumo_estim_low_F;
      real upp_lumo_F_bound  = info.lumo_estim_upp_F;
      real low_homo_F_bound  = info.homo_estim_low_F;
      real upp_homo_F_bound  = info.homo_estim_upp_F;

      // For small cases bounds can be too good or even slightly wrong
      // due to numerical errors; thus we allow a small flexibility
      real flex_tolerance = THRESHOLD_EIG_TOLERANCE;

      if (info.homo_eigenvector_is_computed) // check if eigenvector was computed
      {
         if ((low_homo_F_bound - flex_tolerance <= eigValHOMO) && (eigValHOMO <= upp_homo_F_bound + flex_tolerance))
         {
            do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "HOMO eigenvalue is %lf , HOMO bounds are [ %lf , %lf ]",
                      (double)eigValHOMO, (double)low_homo_F_bound, (double)upp_homo_F_bound);
            info.eigValHOMO = eigValHOMO;
         }
         else
         {
            do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Computed HOMO eigenvalue is outside HOMO bounds [ %lf , %lf ],"
                                                        " discard computed HOMO eigenvector.",
                      (double)low_homo_F_bound, (double)upp_homo_F_bound);
            discard_homo_eigenvector();
         }
      }
      else
      {
         discard_homo_eigenvector();
         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "HOMO eigenvector is not computed.");
      }

      if (info.lumo_eigenvector_is_computed) // check if eigenvector was computed
      {
         if ((low_lumo_F_bound - flex_tolerance <= eigValLUMO) && (eigValLUMO <= upp_lumo_F_bound + flex_tolerance))
         {
            do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "LUMO eigenvalue is %lf , LUMO bounds are [ %lf , %lf ]",
                      (double)eigValLUMO, (double)low_lumo_F_bound, (double)upp_lumo_F_bound);
            info.eigValLUMO = eigValLUMO;
         }
         else
         {
            do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Computed LUMO eigenvalue is outside LUMO bounds [ %lf , %lf ],"
                                                        " discard computed LUMO eigenvector.",
                      (double)low_lumo_F_bound, (double)upp_lumo_F_bound);
            discard_lumo_eigenvector();
         }
      }
      else
      {
         discard_lumo_eigenvector();
         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "LUMO eigenvector is not computed.");
      }

      if (info.homo_eigenvector_is_computed && info.lumo_eigenvector_is_computed)
      {
         real gap = eigValLUMO - eigValHOMO;
         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Computed HOMO-LUMO gap is %lf = %lf eV", (double)gap, (double)(gap / UNIT_one_eV));
      }
   }
}


/***************** COMPUTE_X *****************/

template<typename MatrixType>
void PurificationGeneral<MatrixType>::compute_X()
{
   if (spectrum_bounds.empty())
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Interval spectrum_bounds is empty in compute_X().");
   }

#ifdef DEBUG_PURI_OUTPUT
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Put eigenvalues of F to the interval [0,1] in reverse order.");
#endif

   real eigmin = spectrum_bounds.low();
   real eigmax = spectrum_bounds.upp();
   X.add_identity(-eigmax);
   X *= ((real)1.0 / (eigmin - eigmax));
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::map_bounds_to_0_1()
{
   if (spectrum_bounds.empty())
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Interval spectrum_bounds is empty in map_bounds_to_0_1().");
   }

#ifdef DEBUG_PURI_OUTPUT
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Transform homo and lumo bounds...");
#endif

   real eigmin = spectrum_bounds.low();
   real eigmax = spectrum_bounds.upp();

   // Compute transformed homo and lumo eigenvalues.

   homo_bounds = homo_bounds_F;
   lumo_bounds = lumo_bounds_F;
   // homo and lumo must be in the [lmin, lmax] interval
   homo_bounds.intersect(spectrum_bounds);
   lumo_bounds.intersect(spectrum_bounds);

#ifdef DEBUG_PURI_OUTPUT
   if (homo_bounds.empty())
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Interval homo_bounds is empty.");
   }
   if (lumo_bounds.empty())
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Interval lumo_bounds is empty.");
   }
#endif

   if (!mat::Interval<real>::intersect(homo_bounds, lumo_bounds).empty())
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Bounds for homo and lumo of F are overlapping.");
   }

   // homo_bounds : (1-homo) bounds for matrix X0, later for matrix Xi
   // homo_bounds_X0 : homo bounds for matrix X0
   homo_bounds    = (homo_bounds - eigmax) / (eigmin - eigmax);
   homo_bounds_X0 = homo_bounds;
   homo_bounds    = IntervalType(1 - homo_bounds.upp(), 1 - homo_bounds.low());

   // lumo_bounds : lumo bounds for matrix X0, later for matrix Xi
   // lumo_bounds_X0 : lumo bounds for matrix X0
   lumo_bounds    = (lumo_bounds - eigmax) / (eigmin - eigmax);
   lumo_bounds_X0 = lumo_bounds;

   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "HOMO bounds of X: \t [ %.12lf , %.12lf ]", (double)(1 - homo_bounds.upp()), (double)(1 - homo_bounds.low()));
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "LUMO bounds of X: \t [ %.12lf , %.12lf ]", (double)lumo_bounds.low(), (double)lumo_bounds.upp());

}


/**************************************************************************************/


template<typename MatrixType>
void PurificationGeneral<MatrixType>::set_truncation_parameters()
{
   int estim_num_iter = 0;

   estimate_number_of_iterations(estim_num_iter);
   info.estim_total_it = estim_num_iter;

   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "ESTIMATED NUMBER OF ITERATIONS IS %d", estim_num_iter);

   if (estim_num_iter < maxit)
   {
      // Estimated number of iterations estim_num_iter is less than maxit, set maxit to estim_num_iter
      maxit = estim_num_iter;
   }

   // error due to truncation allowed in each iteration
   error_per_it = error_sub / estim_num_iter;

#ifdef DEBUG_OUTPUT
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "The total allowed subspace error is %e", error_sub);
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Then error due to truncation allowed in each iteration is %e", error_per_it);
#endif
}


/****************** SET_SPECTRUM_BOUNDS *****************/

template<typename MatrixType>
void PurificationGeneral<MatrixType>::set_spectrum_bounds(real eigmin, real eigmax)
{
   spectrum_bounds          = IntervalType(eigmin, eigmax);
   computed_spectrum_bounds = true;
}


/****************** GET_SPECTRUM_BOUNDS *****************/

template<typename MatrixType>
void PurificationGeneral<MatrixType>::get_spectrum_bounds(real& eigmin, real& eigmax)
{
   if (!computed_spectrum_bounds)
   {
      compute_spectrum_bounds();
   }
   eigmin = spectrum_bounds.low();
   eigmax = spectrum_bounds.upp();
}


/****************** COMPUTE_SPECTRUM_BOUNDS *****************/

template<typename MatrixType>
void PurificationGeneral<MatrixType>::compute_spectrum_bounds()
{
   // find approximations using Gershgorin bounds
   real eigminG, eigmaxG, eigmin, eigmax;

   Util::TimeMeter total_time_spectrum_bounds;
   X.gershgorin(eigminG, eigmaxG);
   total_time_spectrum_bounds.print(LOG_AREA_DENSFROMF, "gershgorin");


   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Gershgorin bounds: [ %.12lf , %.12lf ]", (double)eigminG, (double)eigmaxG);


   /* ELIAS NOTE 2016-02-08: Expand Gershgorin bounds by a very small
    * amount to avoid problems of misconvergence in case one of the
    * bounds is exact and the gap is zero (in such cases we want
    * convergence failure, not convergence to a solution with wrong
    * occupation). */
   real smallNumberToExpandWith = template_blas_sqrt(mat::getMachineEpsilon<real>());
   eigminG -= smallNumberToExpandWith;
   eigmaxG += smallNumberToExpandWith;
#ifdef DEBUG_PURI_OUTPUT
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "EXPANDED Gershgorin bounds: [ %.12lf , %.12lf ]", (double)eigminG, (double)eigmaxG);
#endif

   eigmin = eigminG;
   eigmax = eigmaxG;

   // Lanczos helps us to improve Gershgorin bounds
#if 1 // 0 - without Lanczos, 1 - with Lanczos
#ifndef USE_CHUNKS_AND_TASKS

   // try to impove with Lanczos algorithm
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Trying to impove bounds using Lanczos algorithm...");

   real       acc = 1e3 * template_blas_sqrt(get_epsilon()); //  this accuracy may be too high for single precision
   MatrixType Xshifted(X);
   Xshifted.add_identity((real)(-1.0) * eigminG);            // Xsh = X - eigmin*I

   int maxIter = 100;
   try
   {
      eigmax = Xshifted.eucl(acc, maxIter) + eigminG + acc;
   }
   catch (mat::AcceptableMaxIter& e)
   {

      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Lanczos failed to find extreme upper eigenvalue within maxiter... using Gershgorin bound");

      eigmax = eigmaxG;
   }

   // Now we want to create Fshifted = ( F - lambdaMaxGers*id ) but we
   // do this starting from the existing Fshifted, correcting it back
   // to F and then subtracting lambdaMaxGers*id.
   Xshifted.add_identity((real)(1.0) * eigminG); // Now Fshifted = F.
   Xshifted.add_identity((real)(-1.0) * eigmaxG);

   try
   {
      eigmin = -Xshifted.eucl(acc, maxIter) + eigmaxG - acc;
   }
   catch (mat::AcceptableMaxIter& e)
   {

      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Lanczos failed to find extreme lower eigenvalue within maxiter... using Gershgorin bound");

      eigmin = eigminG;
   }
#endif // USE_CHUNKS_AND_TASKS
#endif

   // extreme case of matrix with 1 element
   if (eigmin == eigmax)
   {
      real tol = 1e-2;
      eigmin -= tol;
      eigmax += tol;
   }

   spectrum_bounds = IntervalType(eigmin, eigmax);

   computed_spectrum_bounds = true;
}


/****************** CHECK_STOPPING_CRITERION  **********************/

template<typename MatrixType>
void PurificationGeneral<MatrixType>::stopping_criterion(IterationInfo& iter_info, int& stop, real& estim_order)
{
   int  it = iter_info.it;
   real XmX2_norm_it = -1, XmX2_norm_itm2 = -1;

   if (it < check_stopping_criterion_iter)
   {
      return;                                    // do not check the stopping criterion
   }
   if (use_new_stopping_criterion)
   {
      // if spectral norm is used for the etimation of the order
      if (normPuriStopCrit == mat::euclNorm)
      {
         XmX2_norm_it   = iter_info.XmX2_eucl;
         XmX2_norm_itm2 = info.Iterations[it - 2].XmX2_eucl;
      }
      else
      if (normPuriStopCrit == mat::frobNorm)
      {
         XmX2_norm_it   = iter_info.XmX2_fro_norm;
         XmX2_norm_itm2 = info.Iterations[it - 2].XmX2_fro_norm;
      }
      else
      if (normPuriStopCrit == mat::mixedNorm)
      {
         XmX2_norm_it   = iter_info.XmX2_mixed_norm;
         XmX2_norm_itm2 = info.Iterations[it - 2].XmX2_mixed_norm;
      }
      else
      {
         throw std::runtime_error("Error in stopping_criterion() : unknown matrix norm.");
      }

      real XmX2_trace = iter_info.XmX2_trace;
      check_new_stopping_criterion(it, XmX2_norm_it, XmX2_norm_itm2, XmX2_trace, stop, estim_order);
   }
   else // use standard stopping criterion
   {
      if (normPuriStopCrit == mat::euclNorm)
      {
         XmX2_norm_it = iter_info.XmX2_eucl;
      }
      if (normPuriStopCrit == mat::frobNorm)
      {
         XmX2_norm_it = iter_info.XmX2_fro_norm;
      }
      if (normPuriStopCrit == mat::mixedNorm)
      {
         XmX2_norm_it = iter_info.XmX2_mixed_norm;
      }

      estim_order = -1;
      check_standard_stopping_criterion(XmX2_norm_it, stop);
   }
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::check_standard_stopping_criterion(const real XmX2_norm, int& stop)
{

   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Checking standard stopping criterion...    ");

   bool homoLumo_converged = (homo_bounds.upp() < error_eig &&
                              lumo_bounds.upp() < error_eig);
   bool XmX2norm_converged = XmX2_norm < error_eig;
   if (homoLumo_converged || XmX2norm_converged)
   {
      stop = 1;
   }

#ifdef DEBUG_PURI_OUTPUT
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "stop = %d", stop);
#endif
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::check_new_stopping_criterion(const int it, const real XmX2_norm_it, const real XmX2_norm_itm2, const real XmX2_trace,
                                                                   int& stop, real& estim_order)
{
   // must do at least 2 iterations
   if (it < 2)
   {
      estim_order = -1;
      return;
   }


   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Checking stopping criterion...   ");


   real C; // constant may depend on the purification method
   return_constant_C(it, C);
   this->constant_C = C;
   if (C == -1)
   {
      estim_order = -1;
      return;
   }

   estim_order = template_blas_log(XmX2_norm_it / C) / template_blas_log(XmX2_norm_itm2); // rate of convergence - due to overflow may be Inf

   if ((VecPoly[it - 1] != VecPoly[it]) &&                                                // alternating polynomials
       (XmX2_norm_itm2 < 1) &&                                                            // assumption for frob and mixed norms, always true for the spectral norm
       (XmX2_norm_it >= C * template_blas_pow(XmX2_norm_itm2, (real)ORDER)))              // r <= 2 (or smaller value)
   {
      stop = 1;
   }


   if ((stop != 1) && (XmX2_norm_it < get_epsilon() * template_blas_sqrt(template_blas_sqrt(get_epsilon()))))
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "************************************************************************************************************");
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "The norm value went much below machine precision, therefore we stop here since n_max can be underestimated.");
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "************************************************************************************************************");
      stop = 1;
   }



#ifdef DEBUG_PURI_OUTPUT
   if ((VecPoly[it - 1] != VecPoly[it]) && (XmX2_norm_itm2 < 1))
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "e_i =%e, C*e_{i-1}^q = %e", XmX2_norm_it, C * template_blas_pow(XmX2_norm_itm2, (real)ORDER));
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Order of convergence = %lf, stop = %d", (double)estim_order, stop);
   }
   else
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Order of convergence cannot be computed");
   }
#endif
}


/********************* TOTAL_SUBSPACE_ERROR ****************/

template<typename MatrixType>
typename PurificationGeneral<MatrixType>::real
PurificationGeneral<MatrixType>::total_subspace_error(int it)
{
   assert(it <= (int)VecGap.size());

   real error = 0;
   real normE;
   for (int i = 0; i <= it; ++i)
   {
      if (VecGap[i] == -1)
      {
         return -1;                  // gap is not known
      }
      normE  = info.Iterations[i].threshold_X;
      error += normE / (VecGap[i] - normE);
   }

   return error;
}


template<typename MatrixType>
int PurificationGeneral<MatrixType>::get_est_number_of_puri_iterations()
{
   return info.estim_total_it;
}


template<typename MatrixType>
int PurificationGeneral<MatrixType>::get_exact_number_of_puri_iterations()
{
   if (info.converged == 1)
   {
      return info.total_it;
   }
   else
   {
      return -1;
   }
}


/************ GET ESTIMATE OF EIGENVALUES OF F FROM PURIFICATION **************/
template<typename MatrixType>
void PurificationGeneral<MatrixType>::propagate_values_in_each_iter(real value_unocc, real value_occ,
                                                                    VectorTypeReal& out_unocc,
                                                                    VectorTypeReal& out_occ,
                                                                    int nmax)
{
   out_occ.clear();
   out_occ.resize(nmax);
   out_unocc.clear();
   out_unocc.resize(nmax);

   out_unocc[0] = value_unocc;
   out_occ[0]   = value_occ;
   real occ, unocc;

   for (int i = 1; i < nmax; ++i) // note: i starts from 1
   {
      occ   = out_occ[i - 1];
      unocc = out_unocc[i - 1];
      apply_poly_to_eigs(i, occ, unocc);
      out_occ[i]   = occ;
      out_unocc[i] = unocc;
   }
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::get_eigenvalue_estimates(const VectorTypeReal& XmX2_norm_mixed,
                                                               const VectorTypeReal& XmX2_norm_frob,
                                                               const VectorTypeReal& XmX2_trace)
{
   estimate_homo_lumo(XmX2_norm_mixed, XmX2_norm_frob, XmX2_trace);

   VectorTypeReal h_in, l_in;

   real eigmax    = spectrum_bounds.upp();
   real eigmin    = spectrum_bounds.low();
   real homo_in_0 = 1 - (homo_bounds_F_new.low() - eigmax) / (eigmin - eigmax);
   real lumo_in_0 = (lumo_bounds_F_new.upp() - eigmax) / (eigmin - eigmax);

   int n = get_exact_number_of_puri_iterations();
   assert(n > 0);

   propagate_values_in_each_iter(lumo_in_0, homo_in_0, l_in, h_in, n);

   VectorTypeReal YmY2_norm_frob_est;
   get_frob_norm_est(XmX2_norm_frob, h_in, l_in, YmY2_norm_frob_est);
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::get_frob_norm_est(const VectorTypeReal& XmX2_norm_frob,
                                                        const VectorTypeReal& h_in,
                                                        const VectorTypeReal& l_in,
                                                        VectorTypeReal&       YmY2_norm_frob_est)
{
   int n = get_exact_number_of_puri_iterations();

   assert(n > 0);

   YmY2_norm_frob_est.clear();
   YmY2_norm_frob_est.resize(n);

   real th, tl;
   for (int i = 0; i < n; ++i)
   {
      th = h_in[i] - h_in[i] * h_in[i];
      tl = l_in[i] - l_in[i] * l_in[i];
      YmY2_norm_frob_est[i] = template_blas_sqrt(XmX2_norm_frob[i] * XmX2_norm_frob[i] - th * th - tl * tl);
   }
}



/*
 ||X-X^2||^2_F / trace(X-X^2)  <= ||X-X^2||_2 <= ||X-X^2||_m ( <= ||X-X^2||_F),
 ||X-X^2||_m - mixed norm of X-X^2
 */
template<typename MatrixType>
void PurificationGeneral<MatrixType>::estimate_homo_lumo(const VectorTypeReal& XmX2_norm_mixed,
                                                         const VectorTypeReal& XmX2_norm_frob,
                                                         const VectorTypeReal& XmX2_trace)
{
   homo_bounds_F_new = intervalType(-1e22, 1e22);
   lumo_bounds_F_new = intervalType(-1e22, 1e22);

   // do not use additional iterations in the estimation
   int total_it = info.total_it - info.additional_iterations;

   // lumo_out, lumo_in, 1-homo_out, 1-homo_in
   VectorTypeReal bounds_from_it(4);
   VectorTypeReal final_bounds(4, 1); // set all to one

   // criterion for the eligible iterations for the estimation of the bounds
   real STOP_NORM = gammaStopEstim - gammaStopEstim * gammaStopEstim;
   real vi, wi, mi;
   real temp_value;

   VectorTypeReal XmX2_norm_out;
   if (XmX2_norm_mixed.size() == XmX2_norm_frob.size())
   {
      XmX2_norm_out = XmX2_norm_mixed;
   }
   else
   {
      XmX2_norm_out = XmX2_norm_frob;
   }

   for (int it = total_it; it >= 0; it--)
   {
      vi = XmX2_norm_frob[it];
      wi = XmX2_trace[it];
      mi = XmX2_norm_out[it];

      if (vi >= STOP_NORM)
      {
         break;
      }

      if (wi <= template_blas_sqrt(get_epsilon()))
      {
         continue;
      }

      // lumo bounds
      temp_value        = vi * vi / wi;
      bounds_from_it[0] = 0.5 * (1 - template_blas_sqrt(1 - 4 * temp_value));
      bounds_from_it[1] = 0.5 * (1 - template_blas_sqrt(1 - 4 * mi));

      // bounds for 1-homo
      bounds_from_it[2] = bounds_from_it[0];
      bounds_from_it[3] = bounds_from_it[1];


      apply_inverse_poly_vector(it, bounds_from_it);


      final_bounds[0] = std::min(final_bounds[0], bounds_from_it[0]); // outer
      final_bounds[1] = std::min(final_bounds[1], bounds_from_it[1]); // inner

      final_bounds[2] = std::min(final_bounds[2], bounds_from_it[2]); // outer
      final_bounds[3] = std::min(final_bounds[3], bounds_from_it[3]); // inner
   }

   // get bounds for F
   real maxeig = spectrum_bounds.upp();
   real mineig = spectrum_bounds.low();
   lumo_bounds_F_new = IntervalType(maxeig * (1 - final_bounds[1]) + mineig * final_bounds[1],
                                    maxeig * (1 - final_bounds[0]) + mineig * final_bounds[0]);
   homo_bounds_F_new = IntervalType(mineig * (1 - final_bounds[2]) + maxeig * final_bounds[2],
                                    mineig * (1 - final_bounds[3]) + maxeig * final_bounds[3]);
}


/*************************** SAVE MATRIX **************************/

template<typename MatrixType>
void PurificationGeneral<MatrixType>::save_matrix_now(string str)
{
#ifdef USE_CHUNKS_AND_TASKS
   vector<int>  Itmp, I, Jtmp, J;
   vector<real> Vtmp, V;
   X.get_all_values(Itmp, Jtmp, Vtmp);

   size_t nnz = 0;
   // Count nonzeros
   for (size_t i = 0; i < Itmp.size(); i++)
   {
      nnz += (Vtmp[i] != 0);
   }

   I.reserve(nnz);
   J.reserve(nnz);
   V.reserve(nnz);
   // Extract nonzeros
   for (size_t i = 0; i < Itmp.size(); i++)
   {
      if (Vtmp[i] != 0)
      {
         I.push_back(Itmp[i]);
         J.push_back(Jtmp[i]);
         V.push_back(Vtmp[i]);
      }
   }

   string name = "X_" + str + ".mtx";
   if (write_matrix_to_mtx(name.c_str(), I, J, V, X.get_nrows()) == -1)
   {
      throw std::runtime_error("Error in save_matrix_now : error in write_matrix_to_mtx.\n");
   }
#endif
}


/********************************************************
 ***          COMPUTATION OF EIGENVECTORS       ***
 *********************************************************/


// FUNCTION FOR COMPARISON OF THE RESULTS
#ifdef USE_CHUNKS_AND_TASKS

template<typename MatrixType>
void PurificationGeneral<MatrixType>::compute_eigenvectors_without_diagonalization_on_F(const MatrixType& F, int eigensolver_maxiter_for_F)
{
   throw std::runtime_error("compute_eigenvectors_without_diagonalization_on_F() is not implemented for CHTMatrix.");
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::compute_eigenvectors_without_diagonalization_last_iter_proj()
{
   throw std::runtime_error("compute_eigenvectors_without_diagonalization_last_iter_proj() is not implemented for CHTMatrix.");
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::compute_eigenvectors_without_diagonalization(int it, IterationInfo& iter_info)
{
   throw std::runtime_error("compute_eigenvectors_without_diagonalization() is not implemented for CHTMatrix.");
}


#else
// Assumption: F is not in the file
template<typename MatrixType>
void PurificationGeneral<MatrixType>::compute_eigenvectors_without_diagonalization_on_F(const MatrixType& F, int eigensolver_maxiter_for_F)
{
   real start_shift = homo_bounds_F.upp();
   real end_shift   = lumo_bounds_F.low();
   real dist        = (end_shift - start_shift) / 15;
   real acc_eigv    = eigensolver_accuracy;
   real eigval;
   int  eig_num = 1;

   MatrixType Fsq;

   Fsq = (real)1.0 * F * F; // F^2

   real sigma;

   // choose different shifts
   for (sigma = start_shift; sigma < end_shift; sigma += dist)
   {
      MatrixType M(Fsq);
      M.add_identity(sigma * sigma);       // F^2 + sigma^2I
      M += ((real) - 2 * sigma) * F;       // F^2 + sigma^2I - 2sigma*F
      M  = ((real) - 1.0) * M;             // -(F-sigma*I)^2

      vector<real> eigValTmp(1);           // here will be computed eigenvalues of M
      vector<int>  num_iter_solver(1, -1); // number of iterations

      mat::SizesAndBlocks rows;
      F.getRows(rows);
      vector<VectorType> eigVec(1, VectorType(rows)); // here will be computed eigenvectors
      eigVec[0].rand();                               // initial guess
      eigvec::computeEigenvectors(M, acc_eigv, eigValTmp, eigVec, eig_num,
                                  eigenvectors_iterative_method_str, num_iter_solver,
                                  eigensolver_maxiter_for_F);

      eigval = eigvec::compute_rayleigh_quotient<real>(F, eigVec[0]);

      printf("sigma = %lf , eigval = %lf , iters = %d\n", (double)sigma, (double)eigval, num_iter_solver[0]);
   }

   sigma = end_shift;
   {
      MatrixType M(Fsq);
      M.add_identity(sigma * sigma);       // F^2 + sigma^2I
      M += ((real) - 2 * sigma) * F;       // F^2 + sigma^2I - 2sigma*F
      M  = ((real) - 1.0) * M;             // -(F-sigma*I)^2

      vector<real> eigValTmp(1);           // here will be computed eigenvalues of M
      vector<int>  num_iter_solver(1, -1); // number of iterations

      mat::SizesAndBlocks rows;
      F.getRows(rows);
      vector<VectorType> eigVec(1, VectorType(rows)); // here will be computed eigenvectors
      eigVec[0].rand();                               // initial guess
      eigvec::computeEigenvectors(M, acc_eigv, eigValTmp, eigVec, eig_num,
                                  eigenvectors_iterative_method_str, num_iter_solver,
                                  eigensolver_maxiter_for_F);

      eigval = eigvec::compute_rayleigh_quotient<real>(F, eigVec[0]);

      printf("sigma = %lf , eigval = %lf , iters = %d\n", (double)sigma, (double)eigval, num_iter_solver[0]);
   }
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::compute_eigenvectors_without_diagonalization(int it, IterationInfo& iter_info)
{
   real homo_total_time_stop, lumo_total_time_stop, homo_solver_time_stop, lumo_solver_time_stop;

   /* flags deciding on which iteration to compute eigenvectors */
   bool compute_homo_now = false;
   bool compute_lumo_now = false;


   if (compute_eigenvectors_in_each_iteration)
   {
      if (eigenvectors_method == EIG_SQUARE_INT)
      {
         // can we compute homo in this iteration?
         if (eigVecHOMO != NULL)
         {
            for (size_t i = 0; i < good_iterations_homo.size(); ++i)
            {
               if (good_iterations_homo[i] == it)
               {
                  compute_homo_now = true;
                  break;
               }
            }
         }

         // can we compute lumo in this iteration?
         if (eigVecLUMO != NULL)
         {
            for (size_t i = 0; i < good_iterations_lumo.size(); ++i)
            {
               if (good_iterations_lumo[i] == it)
               {
                  compute_lumo_now = true;
                  break;
               }
            }
         }
      }
      else // eigenvectors_method == EIG_PROJECTION_INT
      {
         // for projection method we can compute for each X_i, including X_0
         if (eigVecHOMO != NULL)
         {
            compute_homo_now = true;
         }
         if (eigVecLUMO != NULL)
         {
            compute_lumo_now = true;
         }
      }
   }
   else
   {
      // compute just in chosen iterations
      if (eigVecHOMO != NULL)
      {
         if (it == iter_for_homo)
         {
            compute_homo_now = true;
         }
      }
      if (eigVecLUMO != NULL)
      {
         if (it == iter_for_lumo)
         {
            compute_lumo_now = true;
         }
      }
   }

   if (compute_homo_now && !info.homo_eigenvector_is_computed)
   {
      if (eigenvectors_method == EIG_SQUARE_INT)
      {
         Util::TimeMeter homo_total_time;

         MatrixType M(Xsq); // M = Xsq
         writeToTmpFile(Xsq);

         real sigma = SIGMA_HOMO_VEC[it]; // take precomputed shift
         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "choose shift %lf", (double)sigma);

         /* GET MATRIX */

         M.add_identity(sigma * sigma); // X^2 + sigma^2I
         M += ((real) - 2 * sigma) * X;
         M  = ((real) - 1.0) * M;       // -(X-sigma*I)^2

         Util::TimeMeter homo_solver_time;
         compute_eigenvector(M, eigVecHOMO, it, true);
         homo_solver_time.print(LOG_AREA_DENSFROMF, "compute_eigenvector()");
         homo_solver_time_stop = homo_solver_time.get_elapsed_wall_seconds();

         homo_total_time.print(LOG_AREA_DENSFROMF, "compute homo eigenvector");
         homo_total_time_stop = homo_total_time.get_elapsed_wall_seconds();

         iter_info.homo_eig_solver_time = homo_solver_time_stop;
         iter_info.orbital_homo_time    = homo_total_time_stop;


         M.clear();
         readFromTmpFile(Xsq);
      }
      else // Projection method
      {
         if (compute_eigenvectors_in_each_iteration)
         {
            // Save matrix X_i into the mtx file. We will use it later to compute homo eigenvectors.
            ostringstream name;
            name << "homo_" << it;
            Util::TimeMeter homo_time_save_matrix;
            save_matrix_now(name.str());
            homo_time_save_matrix.print(LOG_AREA_DENSFROMF, "saving  homo matrix into mtx");
         }
         else
         {
            // Save matrix X_i. We will use it later to compute homo eigenvector.
            Util::TimeMeter homo_time_save_matrix;
            writeToTmpFile(X);
            X_homo = X;
            readFromTmpFile(X);
            homo_time_save_matrix.print(LOG_AREA_DENSFROMF, "saving homo matrix using writeToFile");
         }
      }
   }



   if (compute_lumo_now && !info.lumo_eigenvector_is_computed)
   {
      if (eigenvectors_method == EIG_SQUARE_INT)
      {
         Util::TimeMeter lumo_total_time;

         MatrixType M(Xsq); // M = Xsq
         writeToTmpFile(Xsq);

         real sigma = SIGMA_LUMO_VEC[it]; // take precomputed shift
         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "choose shift %lf", (double)sigma);

         /* GET MATRIX */

         M.add_identity(sigma * sigma); // X^2 + sigma^2I
         M += ((real) - 2 * sigma) * X;
         M  = ((real) - 1.0) * M;       // -(X-sigma*I)^2

         Util::TimeMeter lumo_solver_time;
         compute_eigenvector(M, eigVecLUMO, it, false);
         lumo_solver_time.print(LOG_AREA_DENSFROMF, "compute_eigenvector()");
         lumo_solver_time_stop = lumo_solver_time.get_elapsed_wall_seconds();

         lumo_total_time.print(LOG_AREA_DENSFROMF, "compute lumo eigenvector");
         lumo_total_time_stop = lumo_total_time.get_elapsed_wall_seconds();

         iter_info.lumo_eig_solver_time = lumo_solver_time_stop;
         iter_info.orbital_lumo_time    = lumo_total_time_stop;

         M.clear();
         readFromTmpFile(Xsq);
      }
      else // Projection method
      {
         if (compute_eigenvectors_in_each_iteration)
         {
            // Save matrix X_i into the mtx file. We will use it later to compute lumo eigenvectors.
            ostringstream name;
            name << "lumo_" << it;
            Util::TimeMeter lumo_time_save_matrix;
            save_matrix_now(name.str());
            lumo_time_save_matrix.print(LOG_AREA_DENSFROMF, "saving lumo matrix into mtx");
         }
         else
         {
            // Save matrix X_i. We will use it later to compute lumo eigenvector.
            Util::TimeMeter lumo_time_save_matrix;
            writeToTmpFile(X);
            X_lumo = X;
            readFromTmpFile(X);
            lumo_time_save_matrix.print(LOG_AREA_DENSFROMF, "saving lumo matrix using writeToFile");
         }
      }
   }


   if (compute_eigenvectors_in_each_iteration)
   {
      // compute again in the next iteration
      info.homo_eigenvector_is_computed = false;
      info.lumo_eigenvector_is_computed = false;
   }
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::compute_eigenvectors_without_diagonalization_last_iter_proj()
{
   real homo_total_time_stop, lumo_total_time_stop, homo_solver_time_stop, lumo_solver_time_stop;
   real DX_mult_time_homo_stop, DX_mult_time_lumo_stop;

   output_current_memory_usage(LOG_AREA_DENSFROMF, "Before computing eigenvectors:");
   Util::TimeMeter timeMeterWriteAndReadAll;
   std::string     sizesStr = mat::FileWritable::writeAndReadAll();
   timeMeterWriteAndReadAll.print(LOG_AREA_DENSFROMF, "FileWritable::writeAndReadAll");
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, ((std::string)"writeAndReadAll sizesStr: '" + sizesStr).c_str());

   // we can clear here X^2 !
   Xsq.clear();
   assert(eigVecHOMO != NULL);
   assert(eigVecLUMO != NULL);

   if (compute_eigenvectors_in_each_iteration)
   {
      int total_it = info.total_it;

      /*
       * Since we are computing eigenvectors in every iteration,
       * use the same matrix for homo and for lumo
       */
      for (int it = 0; it <= total_it; ++it)
      {
         // read matrix from mtx file
         MatrixType    Xi;
         ostringstream name;
         name << "X_lumo_" << it << ".mtx";
         mat::SizesAndBlocks rows;
         mat::SizesAndBlocks cols;
         X.getRows(rows);
         X.getCols(cols);
         // TODO
         // Xi.read_from_mtx(name.str(), rows, cols);

         Util::TimeMeter homo_total_time;

         MatrixType DXi(X);

         Util::TimeMeter DX_mult_time_homo; // separate timing for multiplication D*X
         MatrixType      TMP(DXi);
         // DXi = 1 * TMP * Xi + 0 * DXi
         MatrixType::ssmmUpperTriangleOnly((real)1.0, TMP, Xi, 0, DXi);
         DX_mult_time_homo.print(LOG_AREA_DENSFROMF, "computing D*X");
         DX_mult_time_homo_stop = DX_mult_time_homo.get_elapsed_wall_seconds();

         // get HOMO matrix
         MatrixType Zh(X);
         Zh -= DXi;  // D-DXi

         Util::TimeMeter homo_solver_time;
         compute_eigenvector(Zh, eigVecHOMO, it, true);
         homo_solver_time.print(LOG_AREA_DENSFROMF, "compute_eigenvector()");
         homo_solver_time_stop = homo_solver_time.get_elapsed_wall_seconds();

         homo_total_time.print(LOG_AREA_DENSFROMF, "computing homo eigenvector");
         homo_total_time_stop = homo_total_time.get_elapsed_wall_seconds();

         Util::TimeMeter lumo_total_time;

         // get LUMO matrix
         DXi -= Xi;
         DXi  = (real)(-1) * DXi; // Xi-DXi

         Util::TimeMeter lumo_solver_time;
         compute_eigenvector(DXi, eigVecLUMO, it, false);
         lumo_solver_time.print(LOG_AREA_DENSFROMF, "compute_eigenvector()");
         lumo_solver_time_stop = lumo_solver_time.get_elapsed_wall_seconds();

         lumo_total_time.print(LOG_AREA_DENSFROMF, "computing lumo eigenvector");
         lumo_total_time_stop = lumo_total_time.get_elapsed_wall_seconds();

         info.Iterations[it].DX_mult_homo_time = DX_mult_time_homo_stop;
         info.Iterations[it].DX_mult_lumo_time = 0;  // dummy in this case since we reuse matrix

         info.Iterations[it].homo_eig_solver_time = homo_solver_time_stop;
         info.Iterations[it].lumo_eig_solver_time = lumo_solver_time_stop;

         info.Iterations[it].orbital_homo_time = homo_total_time_stop;
         info.Iterations[it].orbital_lumo_time = lumo_total_time_stop;
      }

      output_current_memory_usage(LOG_AREA_DENSFROMF, "After computing eigenvectors:");
      return;
   }


   // check if we passed iter_for_homo iteration and saved the matrix
   if (info.total_it >= iter_for_homo)
   {
      // reading X_homo matrix from the bin file
      Util::TimeMeter X_homo_read;
      readFromTmpFile(X_homo);
      X_homo_read.print(LOG_AREA_DENSFROMF, "reading X matrix (for homo) using readFromFile");


      Util::TimeMeter homo_total_time; // total time for homo

      MatrixType DXi(X);

      // multiplying D*Xi
      Util::TimeMeter DX_mult_time_homo;
      MatrixType      TMP(DXi);
      // DXi = 1 * TMP * X_homo + 0 * DXi
      MatrixType::ssmmUpperTriangleOnly((real)1.0, TMP, X_homo, 0, DXi);
      DX_mult_time_homo.print(LOG_AREA_DENSFROMF, "computing D*X (for homo)");
      DX_mult_time_homo_stop = DX_mult_time_homo.get_elapsed_wall_seconds();

      // get HOMO matrix
      // note: we may need DXi for computing lumo, do not overwrite it
      MatrixType Zh(X);
      Zh -= DXi; // D-DXi
      Util::TimeMeter homo_solver_time;
      compute_eigenvector(Zh, eigVecHOMO, iter_for_homo, true);
      homo_solver_time.print(LOG_AREA_DENSFROMF, "compute_eigenvector()");
      homo_solver_time_stop = homo_solver_time.get_elapsed_wall_seconds();

      info.Iterations[iter_for_homo].homo_eig_solver_time = homo_solver_time_stop; // note: here is included just time for compute_eigenvector()
      info.Iterations[iter_for_homo].DX_mult_homo_time    = DX_mult_time_homo_stop;
      Zh.clear();

      homo_total_time.print(LOG_AREA_DENSFROMF, "computing homo eigenvector");
      homo_total_time_stop = homo_total_time.get_elapsed_wall_seconds();
      info.Iterations[iter_for_homo].orbital_homo_time = homo_total_time_stop;


      // if we are computing both homo and lumo in the same iteration
      if (iter_for_homo == iter_for_lumo)
      {
         Util::TimeMeter lumo_total_time;

         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Reuse matrix D*X_i for lumo computations");  // no matrix read

         // get LUMO matrix
         DXi -= X_homo;
         DXi  = (real)(-1) * DXi; // Xi-DXi

         Util::TimeMeter lumo_solver_time;
         compute_eigenvector(DXi, eigVecLUMO, iter_for_lumo, false);
         lumo_solver_time.print(LOG_AREA_DENSFROMF, "compute_eigenvector()");
         lumo_solver_time_stop = lumo_solver_time.get_elapsed_wall_seconds();


         lumo_total_time.print(LOG_AREA_DENSFROMF, "computing lumo eigenvector");
         lumo_total_time_stop = lumo_total_time.get_elapsed_wall_seconds();

         info.Iterations[iter_for_lumo].DX_mult_lumo_time    = 0;                     // we reuse DX matrix
         info.Iterations[iter_for_lumo].lumo_eig_solver_time = lumo_solver_time_stop; // note: here is included just time for compute_eigenvector()
         info.Iterations[iter_for_lumo].orbital_lumo_time    = lumo_total_time_stop;
      }                                                                               // LUMO
   }  // HOMO

   X_homo.clear();

   // check if we passed iter_for_lumo iteration and saved the matrix, and that it was not the same iteration
   // as iter_for_homo
   if ((info.total_it >= iter_for_lumo) && (iter_for_homo != iter_for_lumo))
   {
      Util::TimeMeter X_lumo_read;
      readFromTmpFile(X_lumo);
      X_lumo_read.print(LOG_AREA_DENSFROMF, "reading X matrix (for lumo) using readFromFile");

      Util::TimeMeter lumo_total_time;

      MatrixType DXi(X); // D

      Util::TimeMeter DX_mult_time_lumo;
      MatrixType      TMP(DXi);
      // DXi = 1 * TMP * X_lumo + 0 * DXi
      MatrixType::ssmmUpperTriangleOnly((real)1.0, TMP, X_lumo, 0, DXi);
      DX_mult_time_lumo.print(LOG_AREA_DENSFROMF, "computing D*X (for lumo)");
      DX_mult_time_lumo_stop = DX_mult_time_lumo.get_elapsed_wall_seconds();

      // get LUMO matrix
      DXi -= X_lumo;
      DXi  = (real)(-1) * DXi; // Xi-DXi

      Util::TimeMeter lumo_solver_time;
      compute_eigenvector(DXi, eigVecLUMO, iter_for_lumo, false);
      lumo_solver_time.print(LOG_AREA_DENSFROMF, "compute_eigenvector()");
      lumo_solver_time_stop = lumo_solver_time.get_elapsed_wall_seconds();
      info.Iterations[iter_for_lumo].lumo_eig_solver_time = lumo_solver_time_stop;
      info.Iterations[iter_for_lumo].DX_mult_lumo_time    = DX_mult_time_lumo_stop;

      lumo_total_time.print(LOG_AREA_DENSFROMF, "computing lumo eigenvector");
      lumo_total_time_stop = lumo_total_time.get_elapsed_wall_seconds();
      info.Iterations[iter_for_lumo].orbital_lumo_time = lumo_total_time_stop;

      X_lumo.clear();
   }  // LUMO
}


#endif // USE_CHUNKS_AND_TASKS



template<typename MatrixType>
void PurificationGeneral<MatrixType>::determine_iteration_for_eigenvectors()
{
   if (eigenvectors_method == EIG_SQUARE_INT)
   {
      get_iterations_for_lumo_and_homo(iter_for_lumo, iter_for_homo);
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Eigenvector for HOMO will be computed on the iteration %d. ", iter_for_homo);
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Eigenvector for LUMO will be computed on the iteration %d. ", iter_for_lumo);
   }
   else if (eigenvectors_method == EIG_PROJECTION_INT)
   {
      // Compute eigenvectors for X_0.
      // Projection method is used just for comparison, thus either eigenvectors
      // are computed only on the first iteration or in each iteration.
      iter_for_lumo = 0;
      iter_for_homo = 0;
   }
   else
   {
      throw std::runtime_error("Error in determine_iteration_for_eigenvectors: unknown method for computing eigenvectors.");
   }
}


/** \brief Find the best iterations for computing eigenvectors.
 *
 *  The best iteration defined by both a small error in the eigenvector
 *  and a small number of iterations of an iterative method.
 *********************************************************************/
template<typename MatrixType>
void PurificationGeneral<MatrixType>::get_iterations_for_lumo_and_homo(int& chosen_iter_lumo,
                                                                       int& chosen_iter_homo)
{
   int  maxiter = maxit;
   // Inner bounds (from the previous SCF cycle {i-1}, expanded
   // with norm of F_i - F_{i-1}).
   real homo0 = 1 - homo_bounds.upp(); // inner bounds
   real lumo0 = lumo_bounds.upp();
   real homoi = homo0, lumoi = lumo0;
   real dummy = 0;
   real Dh_homo, Dh_lumo, Dgh_homo, Dgh_lumo,
        Dgh_homo_max = get_min_double(), Dgh_lumo_max = get_min_double();

   chosen_iter_lumo = -1;
   chosen_iter_homo = -1;

   good_iterations_homo.clear();
   good_iterations_lumo.clear();

   find_shifts_every_iter();

   for (int i = 1; i <= maxiter; ++i)
   {
      homoi = apply_poly(i, homoi); // apply POLY[i] on homo
      lumoi = apply_poly(i, lumoi); // apply POLY[i] on lumo

      Dh_homo = compute_derivative(i, homo0, dummy);
      Dh_lumo = compute_derivative(i, lumo0, dummy);

      // derivative in every iteration
      Dgh_homo = 2 * (homoi - SIGMA_HOMO_VEC[i]) * Dh_homo;
      Dgh_lumo = 2 * (lumoi - SIGMA_LUMO_VEC[i]) * Dh_lumo;

      if (homoi >= SIGMA_HOMO_VEC[i])
      {
         good_iterations_homo.push_back(i);
         if (Dgh_homo >= Dgh_homo_max)
         {
            Dgh_homo_max     = Dgh_homo;
            chosen_iter_homo = i;
         }
      }   // else we cannot be sure which eigenvalue we are computing

      if (lumoi <= SIGMA_LUMO_VEC[i])
      {
         good_iterations_lumo.push_back(i);
         if (template_blas_fabs(Dgh_lumo) >= Dgh_lumo_max) // derivative for lumo is negative
         {
            Dgh_lumo_max     = template_blas_fabs(Dgh_lumo);
            chosen_iter_lumo = i;
         }
      }   // else we cannot be sure which eigenvalue we are computing
   }

   if ((chosen_iter_homo == -1) || (chosen_iter_lumo == -1))
   {
      throw "Error in get_iterations_for_lumo_and_homo() : something went wrong, cannot choose iteration to compute HOMO or LUMO eigenvector.";
   }
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::find_truncation_thresh_every_iter()
{
   int maxiter = this->maxit;

   this->ITER_ERROR_VEC.clear();
   this->ITER_ERROR_VEC.resize(maxiter + 1);

   real error = error_per_it;
   if (error_per_it == 0)
   {
      error = this->get_epsilon();
   }

   for (int i = 1; i <= maxiter; ++i)
   {
      this->ITER_ERROR_VEC[i] = (error * this->VecGap[i]) / (1 + error);
   }
}


/** /brief Find shifts sigma which will be used for construction of
 * the filtering polynomial for computing eigenvectors.
 */
template<typename MatrixType>
void PurificationGeneral<MatrixType>::find_shifts_every_iter()
{
   int maxiter = maxit;

   SIGMA_HOMO_VEC.resize(maxiter + 1);
   SIGMA_LUMO_VEC.resize(maxiter + 1);

   // Inner bounds (from the previous SCF cycle {i-1}, expanded
   // with norm of F_i - F_{i-1}).
   real homo     = 1 - homo_bounds.upp();   // inner bounds
   real lumo     = lumo_bounds.upp();
   real homo_out = 1 - homo_bounds.low();   // outer bounds
   real lumo_out = lumo_bounds.low();

   for (int i = 1; i <= maxiter; ++i)
   {
      homo = apply_poly(i, homo); // apply POLY[i] on homo
      lumo = apply_poly(i, lumo); // apply POLY[i] on lumo

      homo_out = apply_poly(i, homo_out); // apply POLY[i] on homo_out
      lumo_out = apply_poly(i, lumo_out); // apply POLY[i] on lumo_out

      SIGMA_HOMO_VEC[i] = (homo_out + lumo) / 2;
      SIGMA_LUMO_VEC[i] = (lumo_out + homo) / 2;
   }
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::find_eig_gaps_every_iter()
{
   int maxiter = maxit;

   EIG_REL_GAP_HOMO_VEC.resize(maxiter + 1);
   EIG_REL_GAP_LUMO_VEC.clear();
   EIG_REL_GAP_LUMO_VEC.resize(maxiter + 1);
   EIG_ABS_GAP_HOMO_VEC.clear();
   EIG_ABS_GAP_HOMO_VEC.resize(maxiter + 1);
   EIG_ABS_GAP_LUMO_VEC.clear();
   EIG_ABS_GAP_LUMO_VEC.resize(maxiter + 1);

   // Inner bounds (from the previous SCF cycle {i-1}, expanded
   // with norm of F_i - F_{i-1}).
   real homo0 = homo_bounds_X0.low(); // inner bounds
   real lumo0 = lumo_bounds_X0.upp();
   real one   = 1.0;
   real zero  = 0.0;

   real homo_map, lumo_map, one_map, zero_map, sigma;

   real homo = homo0, lumo = lumo0;

   for (int i = 1; i <= maxiter; ++i)
   {
      homo = apply_poly(i, homo); // apply POLY[i] on homo
      lumo = apply_poly(i, lumo); // apply POLY[i] on lumo

      sigma = SIGMA_HOMO_VEC[i];

      homo_map = (homo - sigma) * (homo - sigma);
      lumo_map = (lumo - sigma) * (lumo - sigma);
      one_map  = (one - sigma) * (one - sigma);
      zero_map = (zero - sigma) * (zero - sigma);

      EIG_ABS_GAP_HOMO_VEC[i] = min(lumo_map - homo_map, one_map - homo_map);
      EIG_REL_GAP_HOMO_VEC[i] = EIG_ABS_GAP_HOMO_VEC[i] /
                                max(zero_map - homo_map, one_map - homo_map);
   }

   homo = homo0, lumo = lumo0;

   for (int i = 1; i <= maxiter; ++i)
   {
      homo = apply_poly(i, homo); // apply POLY[i] on homo
      lumo = apply_poly(i, lumo); // apply POLY[i] on lumo

      sigma = SIGMA_LUMO_VEC[i];

      homo_map = (homo - sigma) * (homo - sigma);
      lumo_map = (lumo - sigma) * (lumo - sigma);
      zero_map = (zero - sigma) * (zero - sigma);
      one_map  = (one - sigma) * (one - sigma);

      EIG_ABS_GAP_LUMO_VEC[i] = min(homo_map - lumo_map, zero_map - lumo_map);
      EIG_REL_GAP_LUMO_VEC[i] = EIG_ABS_GAP_LUMO_VEC[i] /
                                max(zero_map - lumo_map, one_map - lumo_map);
   }

}


#ifdef USE_CHUNKS_AND_TASKS
template<typename MatrixType>
void PurificationGeneral<MatrixType>::compute_eigenvector(MatrixType const& M, VectorType *eigVecHOMOorLUMO, int it, bool is_homo)
{
   throw "compute_eigenvector() is not implemented for CHTMatrix.";
}


#else
template<typename MatrixType>
void PurificationGeneral<MatrixType>::compute_eigenvector(MatrixType const& M, VectorType *eigVecHOMOorLUMO, int it, bool is_homo)
{
  assert(eigVecHOMOorLUMO != NULL);
  real acc_eigv = eigensolver_accuracy;
  
  #ifdef DEBUG_PURI_OUTPUT
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Starting compute_eigenvector()");
  #endif
  
  
  if (compute_eigenvectors_in_each_iteration && use_prev_vector_as_initial_guess)
  {
    // copy vector from previous SCF cycle to be an initial guess
    if (is_homo)
    {
      *eigVecHOMO = eigVecHOMORef;
    }
    else
    {
      *eigVecLUMO = eigVecLUMORef;
    }
  }
  
  
  vector<real>        eigValTmp(number_of_eigenvalues); // here will be computed eigenvalues of M
  mat::SizesAndBlocks rows;
  
  
  X.getRows(rows);
  
  /*
  * Apparently the std::vector constructor calls the copy constructor which is not allowed if the data structure of VectorType is not set.
  * In VectorGeneral class were added a new constructor receiving data structure.
  */
  vector<VectorType> eigVec(number_of_eigenvalues, VectorType(rows)); // here will be computed eigenvectors
  if (use_prev_vector_as_initial_guess)                               // copy vector from previous SCF cycle to be an initial guess
  {
    use_prev_vector_as_initial_guess = 0;                            // we firstly check if we have an eigenvector from the previous SCF cycle
    if (is_homo)
    {
      if (!eigVecHOMO->is_empty())
      {
        eigVec[0] = *eigVecHOMO;
        do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Use HOMO eigenvector computed in the previous SCF cycle as initial guess");
        use_prev_vector_as_initial_guess = 1;
      }
    }
    else
    {
      if (!eigVecLUMO->is_empty())
      {
        eigVec[0] = *eigVecLUMO;
        do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Use LUMO eigenvector computed in the previous SCF cycle as initial guess");
        use_prev_vector_as_initial_guess = 1;
      }
    }
  }
  else // use random initial guess
  {
    eigVec[0].rand();
  }
  
  
  vector<int> num_iter_solver(number_of_eigenvalues, -1); // number of iterations
  
  Util::TimeMeter computeEigenvectors_time;
  // run non-deflated version
  int eig_num = 1;
  
  eigvec::computeEigenvectors(M, acc_eigv, eigValTmp, eigVec, eig_num,
    eigenvectors_iterative_method_str, num_iter_solver,
    eigensolver_maxiter);
    double eigv_elapsed_wall_sec = computeEigenvectors_time.get_elapsed_wall_seconds();
    computeEigenvectors_time.print(LOG_AREA_DENSFROMF, "eigensolver");
    
    if (num_iter_solver.empty())
    {
      throw std::runtime_error("Error in compute_eigenvector() : (num_iter_solver.empty())");
    }
    
    
    // initialize
    bool is_homo_tmp = false, is_lumo_tmp = false;
    
    if (num_iter_solver[0] == eigensolver_maxiter)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Eigensolver did not converge within maxiter = %d iterations.", eigensolver_maxiter);
    }
    else
    {
      is_homo_tmp = is_homo, is_lumo_tmp = !is_homo;
      real eigVal;
      get_interval_with_lambda(eigVal, eigVec[0], is_homo_tmp, is_lumo_tmp, it); // compute also the corresponding eigenvalue of F
      if (is_homo_tmp)
      {
        really_good_iterations_homo.push_back(it);
        *eigVecHOMO = eigVec[0];
        info.homo_eigenvector_is_computed         = true;
        info.homo_eigenvector_is_computed_in_iter = it;
        info.homo_eigensolver_iter = num_iter_solver[0];
        info.homo_eigensolver_time = eigv_elapsed_wall_sec;
        eigValHOMO = eigVal;
        do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "compute_eigenvector() for HOMO in iteration %d "
        ": %d iterations, %lf wall sec", it, num_iter_solver[0], eigv_elapsed_wall_sec);
      }
      else if (is_lumo_tmp)
      {
        really_good_iterations_lumo.push_back(it);  // iterations where we computed lumo (for any reason)
        *eigVecLUMO = eigVec[0];
        info.lumo_eigenvector_is_computed         = true;
        info.lumo_eigenvector_is_computed_in_iter = it;
        info.lumo_eigensolver_iter = num_iter_solver[0];
        info.lumo_eigensolver_time = eigv_elapsed_wall_sec;
        eigValLUMO = eigVal;
        do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "compute_eigenvector() for LUMO in iteration %d "
        ": %d iterations, %lf wall sec", it, num_iter_solver[0], eigv_elapsed_wall_sec);
      }
      else
      {
        do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "compute_eigenvector() in iteration %d "
        ": number of eigensolver iterations is %d", it, num_iter_solver[0]);
        do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Error in compute_eigenvector: wrong eigenvalue ");
      }
    }
    
    if (compute_eigenvectors_in_each_iteration)
    {
      save_eigenvectors_to_file(is_homo_tmp, is_lumo_tmp, it);
    }
  }


#endif // USE_CHUNKS_AND_TASKS


template<typename MatrixType>
void PurificationGeneral<MatrixType>::
   writeToTmpFile(MatrixType& A) const
{
#ifndef USE_CHUNKS_AND_TASKS
   A.writeToFile();
#endif
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::
   readFromTmpFile(MatrixType& A) const
{
#ifndef USE_CHUNKS_AND_TASKS
   A.readFromFile();
#endif
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::
   save_eigenvectors_to_file(bool is_homo, bool is_lumo, int it)
{
   if (is_lumo || is_homo)
   {
      std::vector<real> fullVector;
      string            eig_name;
      if (is_homo)
      {
         eigVecHOMO->fullvector(fullVector);
         eig_name = "homo";
      }
      else
      {
         eigVecLUMO->fullvector(fullVector);
         eig_name = "lumo";
      }

      // save vector to file
      ostringstream name;
      if (scf_step != -1)
      {
         name << eig_name << "_" << it << "_scf_step_" << scf_step << ".txt";
      }
      else
      {
         name << eig_name << "_" << it << ".txt";
      }
      if (write_vector(name.str().c_str(), fullVector) == -1)
      {
         throw std::runtime_error("Error in save_eigenvectors_to_file() : error in write_vector.");
      }
      name.str("");
   }
}


// FIXME: add const to vector
template<typename MatrixType>
void PurificationGeneral<MatrixType>::get_eigenvalue_of_F_from_eigv_of_Xi(real& eigVal, const VectorType& eigVec)
{
   readFromTmpFile(F);
   real approx_eig = eigvec::compute_rayleigh_quotient<real>(F, eigVec);
   writeToTmpFile(F);
   eigVal = approx_eig;
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::get_interval_with_lambda(real& eigVal, VectorType& eigVec, bool& is_homo, bool& is_lumo, const int iter)
{
   assert(is_homo || is_lumo);

   bool is_homo_init = is_homo;
   bool is_lumo_init = is_lumo;

   /* COMPUTE EIGENVALUE OF F CORRESPONDING TO THE COMPUTED EIGENVECTOR USING RAYLEIGH QUOTIENT. */

   /*
    * Note: For the square method we compute eigenvalues in the current
    * iteration during the purification. The bounds lumo_bounds and
    * homo_bounds are changing in every iteration according to the
    * polynomial (without expansion by tau). Thus we should use this
    * bounds if square method is used.  However, for the projection
    * method we should used bounds which are saved into the info object.
    */
   real low_lumo_F_bound, low_homo_F_bound;
   real upp_lumo_F_bound, upp_homo_F_bound;

   if (eigenvectors_method == EIG_SQUARE_INT)
   // for square method use bounds from the previous SCF cycle (i.e. bounds expanded with norm ||F_i-F_{i-1}||)
   {
      low_lumo_F_bound = lumo_bounds_F.low();
      upp_lumo_F_bound = lumo_bounds_F.upp();
      low_homo_F_bound = homo_bounds_F.low();
      upp_homo_F_bound = homo_bounds_F.upp();
   }
   else if (eigenvectors_method == EIG_PROJECTION_INT)
   // for projection method we can use new bounds computed in this SCF cycle
   {
      low_lumo_F_bound = info.lumo_estim_low_F;
      upp_lumo_F_bound = info.lumo_estim_upp_F;
      low_homo_F_bound = info.homo_estim_low_F;
      upp_homo_F_bound = info.homo_estim_upp_F;
   }
   else
   {
      throw std::runtime_error("Error in get_interval_with_lambda() : unexpected eigenvectors_method value.");
   }

#ifdef DEBUG_PURI_OUTPUT
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Check Rayleigh quotient...");
#endif

   readFromTmpFile(F);
   real approx_eig = eigvec::compute_rayleigh_quotient<real>(F, eigVec);
   writeToTmpFile(F);

   eigVal = approx_eig;

   real flex_tolerance = THRESHOLD_EIG_TOLERANCE;

   // it is HOMO
   if ((approx_eig <= upp_homo_F_bound + flex_tolerance) && (approx_eig >= low_homo_F_bound - flex_tolerance))
   {
      is_homo = true;
      is_lumo = false;

      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Computed HOMO eigenvalue of F is %lf, "
      "HOMO bounds are  [ %lf , %lf ]", (double)approx_eig, (double)low_homo_F_bound, (double)upp_homo_F_bound);

      iter_for_homo = iter; // We already computed homo in this iteration (in case we thought that it was lumo)
      // Do we want to recompute lumo in the next iteration?
      if (is_lumo_init && (eigenvectors_method == EIG_SQUARE_INT) && try_eigv_on_next_iteration_if_fail)
      {
         iter_for_lumo = iter_for_lumo + 1;
         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "We will try to compute LUMO in the next iteration.");
      }
   }
   // it is LUMO
   else if ((approx_eig <= upp_lumo_F_bound + flex_tolerance) && (approx_eig >= low_lumo_F_bound - flex_tolerance))
   {
      is_homo = false;
      is_lumo = true;

      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Computed LUMO eigenvalue of F is %lf, "
      "LUMO interval [ %lf , %lf ]", (double)approx_eig, (double)low_lumo_F_bound, (double)upp_lumo_F_bound);

      iter_for_lumo = iter; // We already computed lumo in this iteration (in case we thought that it was homo)
      // Do we want to recompute homo in the next iteration?
      if (is_homo_init && (eigenvectors_method == EIG_SQUARE_INT) && try_eigv_on_next_iteration_if_fail)
      {
         iter_for_homo = iter_for_homo + 1;
         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "We will try to compute HOMO in the next iteration.");
      }
   }
   else
   {
      is_homo = false;
      is_lumo = false;

      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Eigenvalue is outside of both intervals for homo and lumo.");
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Eigenvalue is %lf,  HOMO interval [ %lf , %lf ], LUMO interval [ %lf , %lf ]",
                (double)approx_eig, (double)low_homo_F_bound, (double)upp_homo_F_bound, (double)low_lumo_F_bound, (double)upp_lumo_F_bound);

      // Do we want to recompute homo (or lumo) in the next iteration?
      // We will try to compute HOMO, however, it can be LUMO.
      // If it will be LUMO, we save it as computed LUMO and continue with computing HOMO in the next iteration.
      if ((eigenvectors_method == EIG_SQUARE_INT) && try_eigv_on_next_iteration_if_fail)
      {
         iter_for_homo = iter_for_homo + 1;
         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "We will try to compute HOMO (or LUMO) in the next iteration.");
      }
   }
}



/******************************************************************/
/*********************** MATLAB FUNCTIONS *************************/
/******************************************************************/

template<typename MatrixType>
void PurificationGeneral<MatrixType>::gen_matlab_file_norm_diff(const char *filename) const
{
   std::ofstream f;
   f.open(filename, std::ios::out);
   if (!f.is_open())
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Error: cannot open file");
      return;
   }

   int it = info.total_it;
   // save POLY
   f << "POLY = [";
   for (int i = 0; i <= it; ++i)
   {
      f << VecPoly[i] << " ";
   }
   f << "];" << std::endl;

   // choose norm
   if (normPuriStopCrit == mat::euclNorm)
   {
      f << "X_norm = [";
      for (int i = 0; i <= it; ++i) // including 0 iteration = original matrix X
      {
         f << (double)info.Iterations[i].XmX2_eucl << " ";
      }
      f << "];" << std::endl;
      f << " norm_letter = '2';" << std::endl;
   }
   else if (normPuriStopCrit == mat::frobNorm)
   {
      f << "X_norm = [";
      for (int i = 0; i <= it; ++i) // including 0 iteration = original matrix X
      {
         f << (double)info.Iterations[i].XmX2_fro_norm << " ";
      }
      f << "];" << std::endl;
      f << " norm_letter = 'F';" << std::endl;
   }
   else if (normPuriStopCrit == mat::mixedNorm)
   {
      f << "X_norm = [";
      for (int i = 0; i <= it; ++i) // including 0 iteration = original matrix X
      {
         f << (double)info.Iterations[i].XmX2_mixed_norm << " ";
      }
      f << "];" << std::endl;
      f << " norm_letter = 'M';" << std::endl;
   }
   else
   {
      throw "Wrong norm in PurificationGeneral::gen_matlab_file_norm_diff()";
   }

   f << "stop_iteration = " << it - info.additional_iterations << ";" << std::endl;
   f << "it = " << it << ";" << std::endl;
   f << "plot_props = {'LineWidth', 2, 'MarkerSize', 8};" << std::endl;
   f << "fighandle = figure; clf;" << std::endl;
   f << "MARKER = ['o', '+'];" << std::endl;
   f << "semilogy(0:stop_iteration, X_norm(1:stop_iteration+1), '-b', plot_props{:});" << std::endl;
   f << "hold on" << std::endl;
   f << "for i = 1:stop_iteration+1" << std::endl;
   f << " if POLY(i) == 1" << std::endl;
   f << "  h1 = semilogy(i-1, X_norm(i), [MARKER((POLY(i) == 1) + 1) 'b'], plot_props{:});" << std::endl;
   f << " else" << std::endl;
   f << "  h2 = semilogy(i-1, X_norm(i), [MARKER((POLY(i) == 1) + 1) 'b'], plot_props{:});" << std::endl;
   f << " end" << std::endl;
   f << "end" << std::endl;
   f << "if stop_iteration ~= it" << std::endl;
   f << "h3 = semilogy(stop_iteration+1:it, X_norm(stop_iteration+2:it+1), '-.vr', plot_props{:});" << std::endl;
   f << "semilogy(stop_iteration:stop_iteration+1, X_norm(stop_iteration+1:stop_iteration+2), '-.r', plot_props{:});" << std::endl;
   f << "legend([h1 h2 h3],{'$x^2$', '$2x-x^2$', 'After stop'}, 'Interpreter','latex', 'Location','SouthWest');" << std::endl;
   f << "else" << std::endl;
   f << "legend([h1 h2],{'$x^2$', '$2x-x^2$'}, 'Interpreter','latex', 'Location','SouthWest');" << std::endl;
   f << "end" << std::endl;
   f << "xlabel('Iteration SP2', 'Interpreter','latex');" << std::endl;
   f << "ylabel({['$\\|X_i-X_i^2\\|_{' norm_letter '}$']},'interpreter','latex');" << std::endl;
   f << "grid on" << std::endl;
   f << "set(gca, 'XMinorGrid','off', 'YMinorGrid','off', 'GridAlpha',0.6, 'GridLineStyle', '--')" << std::endl;
   f << "set(gca,'FontSize',20);" << std::endl;
   f << "xlim([0 it]);" << std::endl;
   f << "ylim([-inf inf]);" << std::endl;
   f << "set(gca,'XTick',[0 5:5:it]);" << std::endl;
   f << "a = 16; S = logspace(-a, 1, a+2);" << std::endl;
   f << "set(gca,'YTick',S(1:2:end));" << std::endl;

   f << "hold off" << std::endl;

   f << "% print( fighandle, '-depsc2', 'norm_diff.eps');" << std::endl;

   f.close();
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::gen_matlab_file_threshold(const char *filename) const
{
   std::ofstream f;
   f.open(filename, std::ios::out);
   if (!f.is_open())
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Error: cannot open file");
      return;
   }

   int it = info.total_it;
   f << "Thresh = [";

   for (int i = 0; i <= it; ++i)
   {
      f << (double)info.Iterations[i].threshold_X << " ";
   }
   f << "];" << std::endl;

   f << "stop_iteration = " << it - info.additional_iterations << ";" << std::endl;
   f << "it = " << it << ";" << std::endl;
   f << "plot_props = {'LineWidth', 2, 'MarkerSize', 8};" << std::endl;
   f << "fighandle = figure; clf;" << std::endl;
   f << "semilogy(0:stop_iteration, Thresh(1:stop_iteration+1), '-vb', plot_props{:});" << std::endl;
   f << "hold on" << std::endl;
   f << "if stop_iteration ~= it" << std::endl;
   f << "semilogy(stop_iteration+1:it, Thresh(stop_iteration+2:it+1), '-^r', plot_props{:});" << std::endl;
   f << "semilogy(stop_iteration:stop_iteration+1, Thresh(stop_iteration+1:stop_iteration+2), '-r', plot_props{:});" << std::endl;
   f << "legend('before stop', 'after stop', 'Location','NorthWest');" << std::endl;
   f << "end" << std::endl;
   f << "grid on" << std::endl;
   f << "set(gca, 'XMinorGrid','off', 'YMinorGrid','off', 'GridAlpha',0.6, 'GridLineStyle', '--')" << std::endl;
   f << "set(gca,'FontSize',20);" << std::endl;
   f << "xlim([0 it]);" << std::endl;
   f << "ylim([-inf inf]);" << std::endl;
   f << "set(gca,'XTick',[0 5:5:it]);" << std::endl;
   f << "hold off" << std::endl;
   f << "xlabel('Iteration SP2', 'interpreter','latex');" << std::endl;
   f << "ylabel('Threshold value', 'interpreter','latex');" << std::endl;
   f << "% print( fighandle, '-depsc2', 'threshold.eps');" << std::endl;

   f.close();
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::gen_matlab_file_nnz(const char *filename) const
{
   std::ofstream f;
   f.open(filename, std::ios::out);
   if (!f.is_open())
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Error: cannot open file");
      return;
   }

   int it = info.total_it;
   f << "NNZ_X = [";

   for (int i = 0; i <= it; ++i)
   {
      f << (double)info.Iterations[i].NNZ_X << " ";
   }
   f << "];" << std::endl;

   f << "NNZ_X2 = [";

   for (int i = 0; i <= it; ++i)
   {
      f << (double)info.Iterations[i].NNZ_X2 << " ";
   }
   f << "];" << std::endl;



   f << "stop_iteration = " << it - info.additional_iterations << ";" << std::endl;
   f << "it = " << it << ";" << std::endl;
   f << "plot_props = {'LineWidth', 2, 'MarkerSize', 8};" << std::endl;
   f << "fighandle = figure; clf;" << std::endl;
   f << "h2 = plot(0:stop_iteration, NNZ_X(1:stop_iteration+1), '-ob', plot_props{:});" << std::endl;
   f << "hold on" << std::endl;
   f << "h1 = plot(0:stop_iteration, NNZ_X2(1:stop_iteration+1), '-vm', plot_props{:});" << std::endl;
   f << "if stop_iteration ~= it" << std::endl;
   f << "plot(stop_iteration+1:it, NNZ_X(stop_iteration+2:it+1), '-vr', plot_props{:});" << std::endl;
   f << "plot(stop_iteration+1:it, NNZ_X2(stop_iteration+2:it+1), '-*r', plot_props{:});" << std::endl;
   f << "plot(stop_iteration:stop_iteration+1, NNZ_X(stop_iteration+1:stop_iteration+2), '-r', plot_props{:});" << std::endl;
   f << "plot(stop_iteration:stop_iteration+1, NNZ_X2(stop_iteration+1:stop_iteration+2), '-r', plot_props{:});" << std::endl;
   f << "end" << std::endl;
   f << "legend([h1, h2], {'$X^2$', '$X$'},  'interpreter','latex');" << std::endl;
   f << "grid on" << std::endl;
   f << "set(gca, 'XMinorGrid','off', 'YMinorGrid','off', 'GridAlpha',0.6, 'GridLineStyle', '--')" << std::endl;
   f << "set(gca,'FontSize',20);" << std::endl;
   f << "xlim([0 it]);" << std::endl;
   f << "ylim([0 inf]);" << std::endl;
   f << "set(gca,'XTick',[0 5:5:it]);" << std::endl;
   f << "hold off" << std::endl;
   f << "xlabel('Iteration SP2', 'interpreter','latex');" << std::endl;
   f << "ylabel('NNZ [\\%]', 'interpreter','latex');" << std::endl;
   f << "% print( fighandle, '-depsc2', 'nnz.eps');" << std::endl;

   f.close();
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::gen_matlab_file_cond_num(const char *filename) const
{
   std::ofstream f;
   f.open(filename, std::ios::out);
   if (!f.is_open())
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Error: cannot open file");
      return;
   }

   int it = info.total_it;
   f << "in= [";

   for (int i = 0; i <= it; ++i)
   {
      f << (double)(1 - info.Iterations[i].homo_bound_upp - info.Iterations[i].lumo_bound_upp) << " ";
   }
   f << "];" << std::endl;

   f << "stop_iteration = " << it - info.additional_iterations << ";" << std::endl;
   f << "plot_props = {'LineWidth', 2, 'MarkerSize', 8};" << std::endl;
   f << "fighandle = figure; clf;" << std::endl;
   f << "plot(0:stop_iteration, 1./in(1:stop_iteration+1), '-*r', plot_props{:});" << std::endl;
   f << "grid on" << std::endl;
   f << "set(gca, 'XMinorGrid','off', 'YMinorGrid','off', 'GridAlpha',0.6, 'GridLineStyle', '--')" << std::endl;
   f << "set(gca,'FontSize',20);" << std::endl;
   f << "xlim([0 it]);" << std::endl;
   f << "ylim([-inf inf]);" << std::endl;
   f << "set(gca,'XTick',[0 5:5:it]);" << std::endl;
   f << "hold off" << std::endl;
   f << "xlabel('Iteration SP2', 'interpreter','latex');" << std::endl;
   f << "ylabel('$\\kappa$', 'interpreter','latex');" << std::endl;
   f << "% print( fighandle, '-depsc2', 'cond.eps');" << std::endl;

   f.close();
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::gen_matlab_file_eigs(const char *filename) const
{
   std::ofstream f;
   f.open(filename, std::ios::out);
   if (!f.is_open())
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Error: cannot open file");
      return;
   }

   int it = info.total_it;
   f << "homo_low= [";

   for (int i = 0; i <= it; ++i)
   {
      f << (double)info.Iterations[i].homo_bound_low << " ";
   }
   f << "];" << std::endl;

   f << "homo_upp= [";

   for (int i = 0; i <= it; ++i)
   {
      f << (double)info.Iterations[i].homo_bound_upp << " ";
   }
   f << "];" << std::endl;

   f << "lumo_low= [";

   for (int i = 0; i <= it; ++i)
   {
      f << (double)info.Iterations[i].lumo_bound_low << " ";
   }
   f << "];" << std::endl;

   f << "lumo_upp= [";

   for (int i = 0; i <= it; ++i)
   {
      f << (double)info.Iterations[i].lumo_bound_upp << " ";
   }
   f << "];" << std::endl;



   f << "stop_iteration = " << it - info.additional_iterations << ";" << std::endl;
   f << "plot_props = {'LineWidth', 2, 'MarkerSize', 8};" << std::endl;
   f << "fighandle = figure; clf;" << std::endl;
   f << "semilogy(0:stop_iteration, homo_upp(1:stop_iteration+1), '-^b', plot_props{:});" << std::endl;
   f << "hold on" << std::endl;
   f << "semilogy(0:stop_iteration, homo_low(1:stop_iteration+1), '-vb', plot_props{:});" << std::endl;
   f << "semilogy(0:stop_iteration, lumo_low(1:stop_iteration+1), '-vr', plot_props{:});" << std::endl;
   f << "semilogy(0:stop_iteration, lumo_upp(1:stop_iteration+1), '-^r', plot_props{:});" << std::endl;
   f << "grid on" << std::endl;
   f << "set(gca, 'XMinorGrid','off', 'YMinorGrid','off', 'GridAlpha',0.6, 'GridLineStyle', '--')" << std::endl;
   f << "set(gca,'FontSize',20);" << std::endl;
   f << "xlim([0 stop_iteration]);" << std::endl;
   f << "set(gca,'XTick',[0 5:5:stop_iteration]);" << std::endl;
   f << "ylim([-inf inf]);" << std::endl;
   f << "hold off" << std::endl;
   f << "xlabel('Iteration SP2', 'interpreter','latex');" << std::endl;
   f << "ylabel('Homo and lumo bounds', 'interpreter','latex');" << std::endl;
   f << "% print( fighandle, '-depsc2', 'eigs.eps');" << std::endl;

   f.close();
}


template<typename MatrixType>
void PurificationGeneral<MatrixType>::gen_matlab_file_time(const char *filename) const
{
   std::ofstream f;
   f.open(filename, std::ios::out);
   if (!f.is_open())
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Error: cannot open file");
      return;
   }

   int it = info.total_it;

   f << "time_total = [";
   for (int i = 0; i <= it; ++i)
   {
      f << std::setprecision(16) << (double)info.Iterations[i].total_time << " ";
   }
   f << "];" << std::endl;

   f << "time_square = [";
   for (int i = 0; i <= it; ++i)
   {
      f << std::setprecision(16) << (double)info.Iterations[i].Xsquare_time << " ";
   }
   f << "];" << std::endl;

   f << "time_trunc = [";
   for (int i = 0; i <= it; ++i)
   {
      f << std::setprecision(16) << (double)info.Iterations[i].trunc_time << " ";
   }
   f << "];" << std::endl;

   if (info.compute_eigenvectors_in_this_SCF_cycle)
   {
      f << "time_eigenvectors_homo = [";
      for (int i = 0; i <= it; ++i)
      {
         f << std::setprecision(16) << (double)info.Iterations[i].orbital_homo_time << " ";
      }
      f << "];" << std::endl;
      f << "time_eigenvectors_lumo = [";
      for (int i = 0; i <= it; ++i)
      {
         f << std::setprecision(16) << (double)info.Iterations[i].orbital_lumo_time << " ";
      }
      f << "];" << std::endl;

      f << "time_solver_homo = [";
      for (int i = 0; i <= it; ++i)
      {
         f << std::setprecision(16) << (double)info.Iterations[i].homo_eig_solver_time << " ";
      }
      f << "];" << std::endl;
      f << "time_solver_lumo = [";
      for (int i = 0; i <= it; ++i)
      {
         f << std::setprecision(16) << (double)info.Iterations[i].lumo_eig_solver_time << " ";
      }
      f << "];" << std::endl;


      if (eigenvectors_method == EIG_SQUARE_INT)
      {
         // time for X^2, truncation, eigensolver for homo and lumo, additional time for homo and lumo and other time
         f << "X = [time_square; time_trunc; time_solver_homo; time_solver_lumo; time_eigenvectors_homo-time_solver_homo;"
              " time_eigenvectors_lumo-time_solver_lumo; "
              " time_total - time_square - time_trunc - time_eigenvectors_homo - time_eigenvectors_lumo];" << std::endl;
      }
      else
      {
         f << "time_DX_homo = [";
         for (int i = 0; i <= it; ++i)
         {
            f << std::setprecision(16) << (double)info.Iterations[i].DX_mult_homo_time << " ";
         }
         f << "];" << std::endl;
         f << "time_DX_lumo = [";
         for (int i = 0; i <= it; ++i)
         {
            f << std::setprecision(16) << (double)info.Iterations[i].DX_mult_lumo_time << " ";
         }
         f << "];" << std::endl;

         // for projection total time of the iteration does not include computation of homo and lumo
         // time for X^2, eigensolver for homo and lumo, DX multiplication
         f << "X = [time_square; time_trunc; time_solver_homo; time_solver_lumo; time_DX_homo; time_DX_lumo;"
              " time_eigenvectors_homo - time_DX_homo - time_solver_homo; time_eigenvectors_lumo - time_DX_lumo - time_solver_lumo;"
              " time_total - time_square - time_trunc];" << std::endl;
      }
   }
   else
   {
      f << "X = [time_square; time_trunc; time_total - time_square - time_trunc];" << std::endl;
   }

   f << "it = " << it << ";" << std::endl;
   f << "xtick = 0:it;" << std::endl;
   f << "fighandle = figure; clf;" << std::endl;
   f << "b=bar(xtick, X', 'stacked');" << std::endl;
   f << "axis tight;" << std::endl;
   f << "grid on" << std::endl;
   f << "set(gca, 'XMinorGrid','off', 'YMinorGrid','off', 'GridAlpha',0.6, 'GridLineStyle', '--')" << std::endl;
   f << "set(gca,'FontSize',20);" << std::endl;
   f << "xlim([0 it]);" << std::endl;
   f << "set(gca,'XTick',[0 5:5:it]);" << std::endl;
   f << "ylim([-inf inf]);" << std::endl;
   f << "xlabel('Iteration SP2', 'interpreter','latex');" << std::endl;
   f << "ylabel('Time [s]', 'interpreter','latex');" << std::endl;
   if (info.compute_eigenvectors_in_this_SCF_cycle)
   {
/*
 * Legend with matlab's bar appear with default settings in reverse or order. Thus we force it to follow the order of the bar manually.
 */
      if (eigenvectors_method == EIG_SQUARE_INT)
      {
         f << "legend(flipud(b(:)), {'other', 'lumo other', 'homo other', 'lumo solver', 'homo solver', 'truncation', '$X^2$'}, 'interpreter','latex');" << std::endl;
      }
      else
      {
         f << "legend(flipud(b(:)), {'other', 'lumo other', 'homo other', 'DX (lumo)', 'DX (homo)', 'lumo solver', 'homo solver', 'truncation', '$X^2$'}, 'interpreter','latex');" << std::endl;
      }
   }
   else
   {
      f << "legend(flipud(b(:)), {'other', 'truncation', '$X^2$'}, 'interpreter','latex');" << std::endl;
   }
   f << "% print( fighandle, '-depsc2', 'time.eps');" << std::endl;

   f.close();
}



/******************************************************************/
/*********************** PYTHON FUNCTIONS *************************/
/******************************************************************/

template<typename MatrixType>
void PurificationGeneral<MatrixType>::gen_python_file_nnz(const char *filename) const
{
   std::ofstream f;
   f.open(filename, std::ios::out);
   if (!f.is_open())
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Error: cannot open file");
      return;
   }

   f << "import numpy as np" << std::endl;
   f << "import pylab as pl" << std::endl;
   f << "import matplotlib.font_manager as font_manager" << std::endl;
   f << "from matplotlib import rc" << std::endl;
   f << "rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})" << std::endl;
   f << "rc('text', usetex=True)" << std::endl;
   f << std::endl;

   int it = info.total_it;
   f << "NNZ_X = np.array([";

   for (int i = 0; i <= it; ++i)
   {
      f << (double)info.Iterations[i].NNZ_X << ", ";
   }
   f << "])" << std::endl;

   f << "NNZ_X2 = np.array([";

   for (int i = 0; i <= it; ++i)
   {
      f << (double)info.Iterations[i].NNZ_X2 << ", ";
   }
   f << "])" << std::endl;

   f << "stop_iteration = " << it - info.additional_iterations << std::endl;
   f << "it = " << it << std::endl;
   f << "font_prop = font_manager.FontProperties(size=20)" << std::endl;
   f << "prop = {'markersize':8, 'fillstyle':'none', 'linewidth':2, 'markeredgewidth':2}" << std::endl;
   f << "fig1 = pl.figure(figsize = (8, 6), num='nnz')" << std::endl;
   f << "p1, = pl.plot(range(0, stop_iteration+1), NNZ_X2[0:stop_iteration+1], '-vm', **prop);" << std::endl;
   f << "p2, = pl.plot(range(0, stop_iteration+1), NNZ_X[0:stop_iteration+1], '-ob', **prop);" << std::endl;
   f << "if stop_iteration != it:" << std::endl;
   f << "        pl.plot(range(stop_iteration+1,it+1), NNZ_X[stop_iteration+1:it+1], '-vr', **prop);" << std::endl;
   f << "        pl.plot(range(stop_iteration+1,it+1), NNZ_X2[stop_iteration+1:it+1], '-*r', **prop);" << std::endl;
   f << "        pl.plot(range(stop_iteration,stop_iteration+2), NNZ_X[stop_iteration:stop_iteration+2], '-r', **prop);" << std::endl;
   f << "        pl.plot(range(stop_iteration,stop_iteration+2), NNZ_X2[stop_iteration:stop_iteration+2], '-r', **prop);" << std::endl;
   f << std::endl;
   f << "pl.xlim(0, it);" << std::endl;
   f << "pl.ylim(ymin=0);" << std::endl;
   f << "pl.xlabel('Iteration SP2', fontproperties=font_prop);" << std::endl;
   f << "pl.ylabel('NNZ [\\%]', fontproperties=font_prop);    " << std::endl;
   f << "pl.legend((p1, p2), ('$X^2$', '$X$'), loc='best', numpoints=1)" << std::endl;
   f << "pl.xticks(range(5, it, 5))" << std::endl;
   f << "pl.tick_params(labelsize=20)" << std::endl;
   f << "pl.grid(which='major', alpha=0.5, color='k', linestyle='--', linewidth=1) " << std::endl;
   f << "pl.savefig('nnz.eps', format='eps', dpi=1000)" << std::endl;
   f << "pl.show()" << std::endl;

   f.close();
}


#endif //HEADER_PURIFICATION_GENERAL
