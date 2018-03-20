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

/** @file GetDensFromFock.h
 *
 *  @brief Routines for getting density matrix from a given Fock
 *         matrix.
 *
 *  @author Anastasia Kruchinina <em>responsible</em>
 *  @author Elias Rudberg
 */

#ifndef GETDENSFROMFOCKHEADER
#define GETDENSFROMFOCKHEADER

#include "realtype.h"
#include "matrix_typedefs.h"
#include "matrix_typedefs_chtml.h"
#include "transform.h"
#include "output.h"



/** GetDensFromFock class containing parameters and functions for computing density matrix.
 *
 * Flags are set to undefined value by default.  User should define
 * them explicitly, otherwise exception is thrown if undefined flag is
 * used.
 */
class GetDensFromFock
{
public:

   static const int UNDEF_VALUE;          // for flags
   static const int UNDEF_VALUE_UINT;
   static const ergo_real UNDEF_VALUE_REAL;
   static const std::string UNDEF_VALUE_STRING;

   static const int SET;
   static const int UNSET;

   void create_checkpoint(symmMatrix&   Finput,       /**< [in] Effective Hamiltonian matrix (written to file) */
                          symmMatrix&   F_ort_prev,   /**< [in/out]
                                                       *    Input: Previous F matrix in orthogonal basis. (written to file)
                                                       *    Output: New F matrix in orthogonal basis ( ZT*Finput*Z ). (written to file) */
                          generalVector *eigVecLUMO,  /**< [out] LUMO eigenvector */
                          generalVector *eigVecHOMO,  /**< [out] HOMO eigenvector */
                          std::string   IDstr         /**< [in] File identificator; added to the name of each file */
                          );


   static void restore_from_checkpoint(GetDensFromFock& DensFromFock,    /**< [out] Instance of GetDensFromFock class contatining all data for computing the density matrix */
                                       symmMatrix&      Finput,          /**< [out] Effective Hamiltonian matrix (written to file) */
                                       symmMatrix&      F_ort_prev,      /**< [out] F matrix in orthogonal basis ( ZT*Finput*Z ). (written to file) */
                                       generalVector    *eigVecLUMO,     /**< [out] LUMO eigenvector */
                                       generalVector    *eigVecHOMO,     /**< [out] HOMO eigenvector */
                                       std::string      checkpoint_path, /**< [out] HOMO eigenvector */
                                       std::string      IDstr,           /**< [in]  File identificator; added to the name of each file. */
                                       int              SCF_step         /**< [in]  SCF step which should be restored; added to the name of each file in given SCF cycle. */
                                       );


   //constructor
   GetDensFromFock()
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Create object from GetDensFromFock.");

      // set all variables and flags to UNDEF_VALUE

      n = UNDEF_VALUE_UINT;
      noOfOccupiedOrbs       = UNDEF_VALUE_UINT;
      factor                 = UNDEF_VALUE_UINT;
      invCholFactor_euclnorm = UNDEF_VALUE_REAL;
      maxMul                 = UNDEF_VALUE_UINT;
      plot_puri_results      = UNDEF_VALUE;
      SCF_step               = UNDEF_VALUE_UINT;

      use_diagonalization     = UNDEF_VALUE;
      use_purification        = UNDEF_VALUE;

      electronicTemperature                = UNDEF_VALUE_REAL;
      gap_expected_lower_bound             = UNDEF_VALUE_REAL;
      eigvalueErrorLimit                   = UNDEF_VALUE_REAL;
      subspaceErrorLimit                   = UNDEF_VALUE_REAL;
      puri_eig_acc_factor_for_guess        = UNDEF_VALUE;
      use_diag_on_error       = UNDEF_VALUE;
      use_diag_on_error_guess = UNDEF_VALUE;

      create_m_files     = UNDEF_VALUE;
      output_homo_and_lumo_eigenvectors    = UNDEF_VALUE;
      ignore_purification_failure          = UNDEF_VALUE;
      use_rand_perturbation_for_alleigsint = UNDEF_VALUE;
      use_acceleration                   = UNDEF_VALUE;
      use_new_stopping_criterion         = UNDEF_VALUE;
      store_all_eigenvalues_to_file      = UNDEF_VALUE;
      try_eigv_on_next_iteration_if_fail = UNDEF_VALUE;

      leavesSizeMax         = UNDEF_VALUE_UINT;
      blocksize             = UNDEF_VALUE_UINT;

      eigenvectors_method                 = UNDEF_VALUE_STRING;
      eigenvectors_iterative_method       = UNDEF_VALUE_STRING;
      use_prev_vector_as_initial_guess    = UNDEF_VALUE;
      puri_compute_eigv_in_each_iteration = UNDEF_VALUE;
      run_shift_and_square_method_on_F    = UNDEF_VALUE;
      save_permuted_F_matrix_in_bin       = UNDEF_VALUE;

      eigensolver_accuracy = UNDEF_VALUE_REAL;
      eigensolver_maxiter  = UNDEF_VALUE_UINT;

      std::string stats_prefix = ""; // default value

      clean_eigs_intervals();
      clean_puri_stats();


      filenameFinput        = "matrix_Finput";
      filenameF_ort_prev    = "matrix_F_ort_prev";
      filenameeigVecLUMO    = "vector_eigVecLUMO";
      filenameeigVecHOMO    = "vector_eigVecHOMO";
      filenameOverlap       = "matrix_Overlap";
      filenameD_ort_prev    = "matrix_D_ort_prev";
      filenameinvCholFactor = "matrix_invCholFactor";
      file_for_basic_types  = "basic_types";
   }

   /** Choose which method to use for computing the density matrix from Fock matrix.
    */
   int get_dens_from_fock(symmMatrix&   Finput,           /**< [in] Effective Hamiltonian matrix. (written to file) */
                          symmMatrix&   resultDens,       /**< [out] Density matrix. (written to file) */
                          symmMatrix&   F_ort_prev,       /**< [in/out]
                                                           *      Input: Previous F matrix in orthogonal basis. (written to file)
                                                           *      Output: New F matrix in orthogonal basis ( ZT*Finput*Z ). (written to file) */
                          generalVector *eigVecLUMO = 0,  /**< [out] LUMO eigenvector */
                          generalVector *eigVecHOMO = 0   /**< [out] HOMO eigenvector */
                          );



   /** Use recursive expansion for computing the density matrix from Fock matrix.
    */
   int get_dens_from_fock_sparse(symmMatrix&   F,                  /**< [in] Effective Hamiltonian matrix. (written to file) */
                                 symmMatrix&   resultDens,         /**< [out] Density matrix. (written to file) */
                                 symmMatrix&   F_ort_prev,         /**< [in/out]
                                                                    *      Input: Previous F matrix in orthogonal basis. (written to file)
                                                                    *      Output: New F matrix in orthogonal basis ( ZT*Finput*Z ). (written to file) */
                                 generalVector *eigVecLUMO = 0,    /**< [out] HOMO eigenvector */
                                 generalVector *eigVecHOMO = 0     /**< [out] HOMO eigenvector */
                                 );


   /** Set bounds for HOMO and LUMO eigenvalues to -/+ inf, thus remove
    *  any known bounds.
    */
   inline void clean_eigs_intervals()
   {
      homoInterval_Finput = intervalType(-1e22, 1e22);
      lumoInterval_Finput = intervalType(-1e22, 1e22);
      homoInterval_Finput = intervalType(-1e22, 1e22);
      lumoInterval_Finput = intervalType(-1e22, 1e22);

      homoInterval_F_ort_prev = intervalType(-1e22, 1e22);
      lumoInterval_F_ort_prev = intervalType(-1e22, 1e22);
      homoInterval_F_ort_prev = intervalType(-1e22, 1e22);
      lumoInterval_F_ort_prev = intervalType(-1e22, 1e22);
   }

   inline void set_SCF_step(int step /**< [in] Current SCF step */)
   { SCF_step = step; }
   inline void unset_SCF_step()
   { SCF_step = UNDEF_VALUE_UINT; }


   /** Plot figures from the  recursive expansion.
    */
   inline void set_generate_figures(std::string str = "" /**< [in] String added to each generated file. */)
   {
      assert(create_m_files != UNDEF_VALUE);
      if (create_m_files == SET)
      {
         assert(SCF_step >= 0);
         plot_puri_results     = SET;
         plot_puri_results_str = str;
      }
   }

   /** Do not plot figures from the  recursive expansion.
    */
   inline void unset_generate_figures()
   {
      assert(create_m_files != UNDEF_VALUE);
      plot_puri_results     = UNSET;
      plot_puri_results_str = "";
   }

   inline void set_general_params(const int                  n_,                   /**< [in] Number of basis functions. */
                                  mat::SizesAndBlocks const& matrixSizesAndBlocks_ /**< [in] Matrix library parameters. */
                                  )
   {
      assert(n_ >= 1);
      n = n_;
      matrixSizesAndBlocks = matrixSizesAndBlocks_;
   }

   inline void set_cht_matrix_params(const int leavesSizeMax_,  /**< [in] CHTMatrix library parameter leavesSizeMax. */
                                     const int blocksize_       /**< [in] CHTMatrix library parameter blocksize. */
                                     )
   {
      assert(leavesSizeMax_ >= 1);
      assert(blocksize_ >= 1);
      leavesSizeMax = leavesSizeMax_;
      blocksize     = blocksize_;
   }

   inline void get_SizesAndBlocks(mat::SizesAndBlocks& matrixSizesAndBlocks_   /**< [out] Matrix library parameters. */
                                  ) const
   {
      matrixSizesAndBlocks_ = matrixSizesAndBlocks;
   }

   /** Set truncation norm used in the recursive expansion.
    *  Possible norms: spectral, Frobenius or mixed.
    */
   inline void set_truncationNormPurification(mat::normType const truncationNormPurification_ /**< [in]  Norm used in truncation. */)
   { truncationNormPurification = truncationNormPurification_; }

   /** Set stopping criterion norm used in the recursive expansion.
    *  Possible norms: spectral, Frobenius or mixed.
    */
   inline void set_stopCriterionNormPurification(mat::normType const stopCriterionNormPurification_ /**< [in] Norm used in the stopping criterion. */)
   { stopCriterionNormPurification = stopCriterionNormPurification_; }


   inline void do_restricted_calculations()
   { factor = 2; }

   inline void do_unrestricted_calculations()
   { factor = 1; }

   inline void set_no_occupied_orbs(int noOfOccupiedOrbs_)
   {
      assert(noOfOccupiedOrbs_ >= 0);
      noOfOccupiedOrbs = noOfOccupiedOrbs_;
   }

   inline void clean_puri_stats()
   { puri_stats.clear(); }


   inline void set_invCholFactor(triangMatrix const& invCholFactor_,
                                 ergo_real           invCholFactor_euclnorm_)
   {
      invCholFactor = invCholFactor_;
      assert(invCholFactor_euclnorm_ >= 0);
      invCholFactor_euclnorm = invCholFactor_euclnorm_;
   }

   inline void set_gap_expected_lower_bound(ergo_real gap_expected_lower_bound_)
   {
      assert(gap_expected_lower_bound_ >= 0);
      gap_expected_lower_bound = gap_expected_lower_bound_;
   }

   /** Set maximum allowed number of iterations in recursive expansion.
    */
   inline void set_purification_maxmul(ergo_real purification_maxmul_)
   {
      assert(purification_maxmul_ > 0);
      maxMul = purification_maxmul_;
   }

   /****  SET/UNSET SECTION *****/

   inline int get_purification_create_m_files() const
   { return create_m_files == SET; }
   inline void set_purification_create_m_files()
   { create_m_files = SET; }
   inline void unset_purification_create_m_files()
   { create_m_files = UNSET; }



   inline int get_output_homo_and_lumo_eigenvectors() const
   { return output_homo_and_lumo_eigenvectors == SET; }
   inline void set_output_homo_and_lumo_eigenvectors()
   { output_homo_and_lumo_eigenvectors = SET; }
   inline void unset_output_homo_and_lumo_eigenvectors()
   { output_homo_and_lumo_eigenvectors = UNSET; }


   inline int get_purification_ignore_failure() const
   { return ignore_purification_failure == SET; }
   inline void set_purification_ignore_failure()
   { ignore_purification_failure = SET; }
   inline void unset_purification_ignore_failure()
   { ignore_purification_failure = UNSET; }


   inline int get_use_rand_perturbation_for_alleigsint() const
   { return use_rand_perturbation_for_alleigsint == SET; }
   inline void set_purification_use_rand_perturbation_for_alleigsint()
   { use_rand_perturbation_for_alleigsint = SET; }
   inline void unset_purification_use_rand_perturbation_for_alleigsint()
   { use_rand_perturbation_for_alleigsint = UNSET; }

   inline int get_use_diagonalization() const
   { return use_diagonalization == SET; }
   inline void set_use_diagonalization()
   { use_diagonalization = SET; }
   inline void unset_use_diagonalization()
   { use_diagonalization = UNSET; }


   inline int get_use_purification() const
   { return use_purification == SET; }
   inline void set_use_purification()
   { use_purification = SET; }
   inline void unset_use_purification()
   { use_purification = UNSET; }


   inline int get_use_diag_on_error_guess() const
   { return use_diag_on_error_guess == SET; }
   inline void set_use_diag_on_error_guess()
   { use_diag_on_error_guess = SET; }
   inline void unset_use_diag_on_error_guess()
   { use_diag_on_error_guess = UNSET; }


   inline int get_use_diag_on_error() const
   { return use_diag_on_error == SET; }
   inline void set_use_diag_on_error()
   { use_diag_on_error = SET; }
   inline void unset_use_diag_on_error()
   { use_diag_on_error = UNSET; }


   inline std::string get_stats_prefix() const
   { return stats_prefix; }
   inline void set_stats_prefix(std::string stats_prefix_)
   { stats_prefix = stats_prefix_; }
   inline void unset_stats_prefix()
   { stats_prefix = ""; }


   inline int get_use_acceleration() const
   { return use_acceleration == SET; }
   inline void set_use_acceleration()
   { use_acceleration = SET; }
   inline void unset_use_acceleration()
   { use_acceleration = UNSET; }

   inline int get_use_new_stopping_criterion() const
   { return use_new_stopping_criterion == SET; }
   inline void set_use_new_stopping_criterion()
   { use_new_stopping_criterion = SET; }
   inline void unset_use_new_stopping_criterion()
   { use_new_stopping_criterion = UNSET; }

   inline int get_store_all_eigenvalues_to_file() const
   { return store_all_eigenvalues_to_file == SET; }
   inline void set_store_all_eigenvalues_to_file()
   { store_all_eigenvalues_to_file = SET; }
   inline void unset_store_all_eigenvalues_to_file()
   { store_all_eigenvalues_to_file = UNSET; }

   inline int get_save_permuted_F_matrix_in_bin()
   { return save_permuted_F_matrix_in_bin; }
   inline void set_save_permuted_F_matrix_in_bin()
   { save_permuted_F_matrix_in_bin = SET; }
   inline void unset_save_permuted_F_matrix_in_bin()
   { save_permuted_F_matrix_in_bin = UNSET; }



   inline int get_puri_compute_eigv_in_each_iteration()
   { return puri_compute_eigv_in_each_iteration; }
   inline void set_puri_compute_eigv_in_each_iteration()
   { puri_compute_eigv_in_each_iteration = SET; }
   inline void unset_puri_compute_eigv_in_each_iteration()
   { puri_compute_eigv_in_each_iteration = UNSET; }


   inline int get_run_shift_and_square_method_on_F()
   { return run_shift_and_square_method_on_F; }
   inline void set_run_shift_and_square_method_on_F()
   { run_shift_and_square_method_on_F = SET; }
   inline void unset_run_shift_and_square_method_on_F()
   { run_shift_and_square_method_on_F = UNSET; }


   inline int get_try_eigv_on_next_iteration_if_fail()
   { return try_eigv_on_next_iteration_if_fail; }
   inline void set_try_eigv_on_next_iteration_if_fail()
   { try_eigv_on_next_iteration_if_fail = SET; }
   inline void unset_try_eigv_on_next_iteration_if_fail()
   { try_eigv_on_next_iteration_if_fail = UNSET; }


   inline int get_use_prev_vector_as_initial_guess()
   { return use_prev_vector_as_initial_guess; }
   inline void set_use_prev_vector_as_initial_guess()
   { use_prev_vector_as_initial_guess = SET; }
   inline void unset_use_prev_vector_as_initial_guess()
   { use_prev_vector_as_initial_guess = UNSET; }


   inline void set_diagonalization_params(ergo_real   electronicTemperature_,
                                          symmMatrix& overlapMatrix_)
   {
      set_overlapMatrix(overlapMatrix_);
      assert(electronicTemperature_ >= 0);
      electronicTemperature = electronicTemperature_;
   }

   inline void set_overlapMatrix(symmMatrix& overlapMatrix_)
   { overlapMatrix = overlapMatrix_; }


   inline void set_purification_limits(ergo_real subspaceErrorLimit_,
                                       ergo_real eigvalueErrorLimit_ = 0,
                                       ergo_real puri_eig_acc_factor_for_guess = 0)
   {
      set_eigvalueErrorLimit(eigvalueErrorLimit_);
      set_subspaceErrorLimit(subspaceErrorLimit_);
      set_puri_eig_acc_factor_for_guess(puri_eig_acc_factor_for_guess);
   }

   /** Set maximum allowed error in eigenvalues of the density matrix.
    */
   inline void set_eigvalueErrorLimit(ergo_real eigvalueErrorLimit_)
   { eigvalueErrorLimit = eigvalueErrorLimit_; }

   /** Set maximum allowed error in invariant subspaces of the density
    * matrix.
    */
   inline void set_subspaceErrorLimit(ergo_real subspaceErrorLimit_)
   { subspaceErrorLimit = subspaceErrorLimit_; }

   /** Set puri_eig_acc_factor_for_guess parameter.
    *
    * Obsolete parameter needed for the old stopping criterion for
    * creating the initial guess.
    */
   inline void set_puri_eig_acc_factor_for_guess(ergo_real puri_eig_acc_factor_for_guess_)
   { puri_eig_acc_factor_for_guess = puri_eig_acc_factor_for_guess_; }



   // get some results from the purification

   ergo_real get_result_entropy_term() const
   { return resultEntropyTerm; }

   inline void get_puri_stats(std::map<std::string, double>& puri_stats_) const
   { puri_stats_ = puri_stats; }



   // Fprev is effective Hamiltonian matrix (=Finput)
   // Intervals contain the homo and lumo eigenvalues of Fprev
   inline void set_eigs_Fprev(intervalType& homoInterval_Finput_,
                              intervalType& lumoInterval_Finput_)
   {
      homoInterval_Finput = intervalType(homoInterval_Finput_);
      lumoInterval_Finput = intervalType(lumoInterval_Finput_);
   }

   inline void get_eigs_Fprev(intervalType& homoInterval_Finput_,
                              intervalType& lumoInterval_Finput_) const
   {
      homoInterval_Finput_ = intervalType(homoInterval_Finput);
      lumoInterval_Finput_ = intervalType(lumoInterval_Finput);
   }

   // F_ort_prev is matrix in orthogonal basis
   // Intervals contain the homo and lumo eigenvalues of F_ort_prev
   inline void set_eigs_F_ort_prev(intervalType& homoInterval_F_ort_prev_,
                                   intervalType& lumoInterval_F_ort_prev_)
   {
      homoInterval_F_ort_prev = intervalType(homoInterval_F_ort_prev_);
      lumoInterval_F_ort_prev = intervalType(lumoInterval_F_ort_prev_);
   }

   inline void get_eigs_F_ort_prev(intervalType& homoInterval_F_ort_prev_,
                                   intervalType& lumoInterval_F_ort_prev_) const
   {
      homoInterval_F_ort_prev_       = intervalType(homoInterval_F_ort_prev);
      lumoInterval_F_ort_prev_       = intervalType(lumoInterval_F_ort_prev);
   }

   inline ergo_real get_eigvalueErrorLimit() const
   { return eigvalueErrorLimit; }
   inline ergo_real get_subspaceErrorLimit() const
   { return subspaceErrorLimit; }
   inline ergo_real get_puri_eig_acc_factor_for_guess() const
   { return puri_eig_acc_factor_for_guess; }



   inline void compute_eigenvectors(std::string eigenvectors_method_,     
                                    std::string eigenvectors_iterative_method_, 
                                    ergo_real eigensolver_accuracy_, 
                                    int eigensolver_maxiter_,
                                    int use_prev_vector_as_initial_guess_, 
                                    int try_eigv_on_next_iteration_if_fail_)
   {
      assert(eigenvectors_method_ == "square" || eigenvectors_method_ == "projection");
      eigenvectors_method = eigenvectors_method_;

      assert(eigenvectors_iterative_method_ == "power" || eigenvectors_iterative_method_ == "lanczos");
      eigenvectors_iterative_method = eigenvectors_iterative_method_;

      eigensolver_accuracy = eigensolver_accuracy_;
      eigensolver_maxiter  = eigensolver_maxiter_;

      if (use_prev_vector_as_initial_guess_ > 0)
      {
         set_use_prev_vector_as_initial_guess();
      }
      else
      {
         unset_use_prev_vector_as_initial_guess();
      }
      if (try_eigv_on_next_iteration_if_fail_ > 0)
      {
         set_try_eigv_on_next_iteration_if_fail();
      }
      else
      {
         unset_try_eigv_on_next_iteration_if_fail();
      }

   }

   inline void compute_eigenvectors_extra(int puri_compute_eigv_in_each_iteration_, int run_shift_and_square_method_on_F_)
   {
      if (puri_compute_eigv_in_each_iteration_ > 0)
      {
         set_puri_compute_eigv_in_each_iteration();
      }
      else
      {
         unset_puri_compute_eigv_in_each_iteration();
      }

      if (run_shift_and_square_method_on_F_ > 0)
      {
         set_run_shift_and_square_method_on_F();
      }
      else
      {
         unset_run_shift_and_square_method_on_F();
      }
   }

private:

   int SCF_step;



   int use_diagonalization;                  /**< Flag to turn on diagonalization. */

   int use_purification;                     /**< Flag to turn on purification. */

   int store_all_eigenvalues_to_file;        /**< Store eigenvalues to the file when doing diagonalization.
                                              * NOTE: works just with diagonalization */

   int try_eigv_on_next_iteration_if_fail;   /**< For square method: if eigenvector is not computed in iteration
                                              * i, try to compute it in iteration i+1 */

   ergo_real electronicTemperature;          /**< Electronic temperature */

   ergo_real gap_expected_lower_bound;       /**< Expected lower bound for the gap to be used in early iterations.  */

   ergo_real eigvalueErrorLimit;             /**< Tolerated deviation of eigenvalues from 0 and 1 in the computed density matrix.  */

   ergo_real subspaceErrorLimit;             /**< Tolerated error in the occupied subspace as measured by the sinus of the largest canonical angle. */

   ergo_real puri_eig_acc_factor_for_guess;  /**< With this number will be multiplied the tolerated deviation of
                                              * eigenvalues from 0 and 1 in the computed density matrix for the
                                              * initial guess density matrix */

   int use_diag_on_error;                    /**< Flag to fall back on diagonalization if purification fails. */
   int use_diag_on_error_guess;

   int create_m_files;                       /**< Flag to create m-files with information about the purification process.  */

   int output_homo_and_lumo_eigenvectors;    /**< Compute homo and lumo eigenvectors and write them to the file */

   int use_prev_vector_as_initial_guess;     /**< Use eigenvector from the previous SCF cycle as an initial guess in this cycle */

   int puri_compute_eigv_in_each_iteration;  /**< Compute eigenvectors in each iteration of the recursive expansion. */

   int run_shift_and_square_method_on_F;     /**< (for comparison) Run shift_and_square method to
                                              * get eigenvectors of the matrix F for various shifts. */

   int save_permuted_F_matrix_in_bin;        /**< Save sparse matrix F into bin file in the current permutation of rows and columns. */

   int ignore_purification_failure;          /**< Continue even if purification fails to converge.  */

   int use_rand_perturbation_for_alleigsint; /**< Apply a random
                                              * perturbation to (try
                                              * to) improve the
                                              * convergence speed of
                                              * Lanczos calculation of
                                              * extremal
                                              * eigenvalues.  */

   std::string stats_prefix;                 /**< Prefix to be added to statistics files. */

   int plot_puri_results;                    /**< Plot results of the purification from this function call */

   std::string plot_puri_results_str;

   int use_acceleration;                        /**< Use acceleration in the purification */

   int use_new_stopping_criterion;              /**< Use new parameterless stopping criterion */

   std::string eigenvectors_method;             /**< Method for computing eigenvectors: square or projection */

   std::string eigenvectors_iterative_method;   /**< Iterative method for computing eigenvectors: power or lanczos */

   ergo_real eigensolver_accuracy;              /**< The accuracy for the eigenvalue problem solver */

   int eigensolver_maxiter;                     /**< Maximum number of iterations for the eigenvalue problem solver */

   int n;                                       /**< System size. */

   int noOfOccupiedOrbs;                        /**< Number of occupied orbitals.  */

   ergo_real factor;                            /**< Factor to scale the resulting density matrix.
                                                 * (for restricted vs unrestricted calc) */

   symmMatrix overlapMatrix;                    /**< Overlap matrix (written to file) */

   symmMatrix D_ort_prev;                       /**< Density matrix from previous SCF cycle (written to file) */

   triangMatrix invCholFactor;                  /**< Inverse Cholesky factor (written to file) */

   ergo_real invCholFactor_euclnorm;            /**< Euclidean norm of inverse Cholesky factor. */

   mat::normType truncationNormPurification;    /**< Norm to be used for truncation.  */

   mat::normType stopCriterionNormPurification; /**< Norm to be used for stopping criterion.  */

   int maxMul;                                  /**< Maximum allowed number of matrix multiplications in the purification */

   mat::SizesAndBlocks matrixSizesAndBlocks;    /**< Information about HML matrix block sizes etc. */

   int leavesSizeMax;                           /**< Information about leavesSizeMax and blocksize for CHTMatrix */
   int blocksize;                               /**< Information about leavesSizeMax and blocksize for CHTMatrix */

   intervalType homoInterval_Finput;
   intervalType lumoInterval_Finput;

   intervalType homoInterval_F_ort_prev;
   intervalType lumoInterval_F_ort_prev;

   ergo_real resultEntropyTerm;
   std::map<std::string, double> puri_stats;


   // Names of files needed for checkpoints
   const char *filenameFinput;
   const char *filenameF_ort_prev;
   const char *filenameeigVecLUMO;
   const char *filenameeigVecHOMO;
   const char *filenameOverlap;
   const char *filenameD_ort_prev;
   const char *filenameinvCholFactor;
   const char *file_for_basic_types;
};



#endif // GETDENSFROMFOCKHEADER
