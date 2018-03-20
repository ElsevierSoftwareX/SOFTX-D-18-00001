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

/** @file scf.h

    @brief Code for classes containing various options related to
    self-consistent field (SCF) calculations.

    @author: Elias Rudberg <em>responsible</em>.
*/

#ifndef SCF_HEADER
#define SCF_HEADER

#include <string.h>

#include "molecule.h"
#include "basisinfo.h"
#include "integrals_2el.h"
#include "matrix_typedefs.h"


namespace SCF {

static const int DISTURB_ELEMENT_MAX_COUNT = 60;

struct Options {
  std::string calculation_identifier;
  std::string method_and_basis_set;
  Vector3D electric_field;
  ergo_real electronic_temperature;
  ergo_real sparse_threshold_for_S;
  ergo_real sparse_threshold_for_Z;
  ergo_real convergence_threshold;
  ergo_real step_length_giveup;
  ergo_real step_length_start;
  ergo_real puri_eig_acc_factor_for_guess;
  ergo_real purification_conv_limit;
  int create_checkpoints;
  std::string checkpoint_IDstr;
  ergo_real purification_eigvalue_err_limit; 
  ergo_real purification_subspace_err_limit;  
  int purification_with_acceleration;
  int use_new_stopping_criterion;
  ergo_real gap_expected_lower_bound;
  mat::normType purification_truncation_norm;
  mat::normType purification_stop_crit_norm;
  int cht_leavesSizeMax;
  int cht_blocksize;
  ergo_real subspace_factor_fock;
  ergo_real subspace_factor_dens;
  int use_artificial_subspace_disturbances;
  int no_of_threads_for_V;
  ergo_real box_size_for_V_and_T;
  int purification_maxmul;
  int purification_create_m_files;
  int purification_ignore_failure;
  int purification_use_rand_perturbation_for_alleigsint;
  int use_dft;
  int use_simple_starting_guess;
  int use_diag_guess_from_file;
  int write_diag_dens_to_file;
  ergo_real starting_guess_disturbance;
  int sg_disturb_specific_elements;
  int disturbedElementIndexVector[DISTURB_ELEMENT_MAX_COUNT];
  ergo_real shift_using_prev_density_matrix;
  int skip_H_core;
  int use_simple_dense_H_core;
  int break_on_energy_increase;
  int force_restricted;  /**< use a restricted determinant for open shell. */
  int force_unrestricted; /**< use an unrestricted det. for closed shell. */
  int spin_flip_atom_count;
  int starting_guess_spin_diff;
  int max_no_of_diis_matrices;
  int max_restart_count;
  int no_of_impr_req_for_diis;
  int use_diis_always;
  int do_f_thresh_verification;
  int output_statistics_mfiles;
  int no_of_careful_first_scf_steps;
  int do_report_density_diff;
  ergo_real error_maxabs_for_diis;
  int min_number_of_iterations;
  int max_number_of_iterations;
  int output_density_at_every_step;
  int output_expected_values_pos_operator;
  int output_density_images;
  int output_density_images_only;
  int write_guess_density_only;
  int compute_core_density;
  int no_of_core_electrons;
  ergo_real output_density_images_boxwidth;
  int image_view_axis;
  int save_final_potential;
  int use_diagonalization;
  int use_diag_on_error;
  int use_diag_on_error_guess;
  int store_all_eigenvalues_to_file;
  int try_eigv_on_next_iteration_if_fail;
  int puri_compute_eigv_in_each_iteration; 
  int run_shift_and_square_method_on_F;
  int save_permuted_F_matrix_in_bin;
  int write_overlap_matrix;
  int save_full_matrices_for_matlab;
  int analyze_result_after_scf;
  int do_acc_scan_J;
  int do_acc_scan_K;
  int do_acc_scan_Vxc;
  int scan_do_invcholfactor_transf;
  int scan_no_of_steps;
  ergo_real scan_start_thresh;
  ergo_real scan_step_factor;
  int create_mtx_file_S;  
  int create_mtx_file_H_core;
  int create_mtx_files_F;  
  int create_mtx_files_D;  
  int create_mtx_files_dipole;
  int create_mtx_files_S_and_quit;
  int create_2el_integral_m_file;
  int create_basis_func_coord_file;
  int use_prev_vector_as_initial_guess;
  int output_homo_and_lumo_eigenvectors;
  std::string eigenvectors_method;	
  std::string eigenvectors_iterative_method;
  ergo_real eigensolver_accuracy;
  int eigensolver_maxiter;
  int output_mulliken_pop;
  int compute_gradient_fixeddens;
  int verify_gradient_fixeddens;

  /** Initializes all the fields to sane values. */
Options() : calculation_identifier("N/A"),
    method_and_basis_set("N/A"),
    electric_field(0,0,0),
    electronic_temperature(0),
    sparse_threshold_for_S(1e-9),
    sparse_threshold_for_Z(1e-8),
    convergence_threshold(2e-7),
    step_length_giveup(0.00005),
    step_length_start(0.4),
    puri_eig_acc_factor_for_guess(1e-2),
    purification_conv_limit(0.1),
    create_checkpoints(0),
    checkpoint_IDstr(""),
    purification_eigvalue_err_limit(1e-8),
    purification_subspace_err_limit(1e-6),
    purification_with_acceleration(0),
    use_new_stopping_criterion(1),
    gap_expected_lower_bound(0.05),
    purification_truncation_norm(mat::mixedNorm),
    purification_stop_crit_norm(mat::mixedNorm),
    cht_leavesSizeMax(1024), 
    cht_blocksize(64),
    subspace_factor_fock(0.1),
    subspace_factor_dens(0.1),
    use_artificial_subspace_disturbances(0),
    no_of_threads_for_V(1),
    box_size_for_V_and_T(6.5),
    purification_maxmul(100),
    purification_create_m_files(0),
    purification_ignore_failure(0),
    purification_use_rand_perturbation_for_alleigsint(0),
    use_dft(0),
    use_simple_starting_guess(0),
    use_diag_guess_from_file(0),
    write_diag_dens_to_file(0),
    starting_guess_disturbance(0.0),
    sg_disturb_specific_elements(0),
    shift_using_prev_density_matrix(0.0),
    skip_H_core(0),
    use_simple_dense_H_core(0),
    break_on_energy_increase(0),
    force_restricted(0),
    force_unrestricted(0),
    spin_flip_atom_count(0),
    starting_guess_spin_diff(0),
    max_no_of_diis_matrices(10),
    max_restart_count(2),
    no_of_impr_req_for_diis(4),
    use_diis_always(0),
    do_f_thresh_verification(0),
    output_statistics_mfiles(0),
    no_of_careful_first_scf_steps(0),
    do_report_density_diff(1),
    error_maxabs_for_diis(0.5),
    min_number_of_iterations(),
    max_number_of_iterations(),
    output_density_at_every_step(1),
    output_expected_values_pos_operator(0),
    output_density_images(0),
    output_density_images_only(0),
    write_guess_density_only(0),
    compute_core_density(0),
    no_of_core_electrons(0),
    output_density_images_boxwidth(0.5),
    image_view_axis(),
    save_final_potential(0),
    use_diagonalization(0),
    use_diag_on_error(1),
    use_diag_on_error_guess(1),
    store_all_eigenvalues_to_file(0),
    try_eigv_on_next_iteration_if_fail(0),
    puri_compute_eigv_in_each_iteration(0), 
    run_shift_and_square_method_on_F(0),
    save_permuted_F_matrix_in_bin(0),
    write_overlap_matrix(0),
    save_full_matrices_for_matlab(0),
    analyze_result_after_scf(0),
    do_acc_scan_J(0),
    do_acc_scan_K(0),
    do_acc_scan_Vxc(0),
    scan_do_invcholfactor_transf(1),
    scan_no_of_steps(16),
    scan_start_thresh(1e-9),
    scan_step_factor(template_blas_sqrt((ergo_real)10)),
    create_mtx_file_S(0),
    create_mtx_file_H_core(0),
    create_mtx_files_F(0),
    create_mtx_files_D(0),
    create_mtx_files_dipole(0),
    create_mtx_files_S_and_quit(0),
    create_2el_integral_m_file(0),
    create_basis_func_coord_file(0),
    use_prev_vector_as_initial_guess(0),
    output_homo_and_lumo_eigenvectors(0),
    eigenvectors_method("square"),
    eigenvectors_iterative_method("lanczos"),
    eigensolver_accuracy(1e-4*template_blas_sqrt(mat::getMachineEpsilon<ergo_real>())),
    eigensolver_maxiter(200),
    output_mulliken_pop(0),
    compute_gradient_fixeddens(0),
    verify_gradient_fixeddens(0)
  { 
    memset(disturbedElementIndexVector, 0,
           sizeof(disturbedElementIndexVector));
  }
};

/** An object respresenting the configuration of the matrix
    library. All the thresholds and relevant parameters are collected
    in one object for the purposes of the input processing. */
struct MatOptions {
  mat::SizesAndBlocks size_block_info;
  std::vector<int> permutationHML;
  std::vector<int> inversePermutationHML;
  ergo_real sparse_threshold; /**< threshold value for sparse matrix
				 truncation. */
  ergo_real threshold_inch; /**< Truncation threshold in INCH function. */
  int sparse_matrix_block_size;
  int sparse_matrix_block_factor_3;
  int sparse_matrix_block_factor_2;
  int sparse_matrix_block_factor_1;
  int threads;
  int parallelLevel;
  int no_of_buffers_per_allocator;
  int use_allocator_manager;

  MatOptions() :
    sparse_threshold(1e-8),
    threshold_inch(1e-10),
    sparse_matrix_block_size(32),
    sparse_matrix_block_factor_3(8),
    sparse_matrix_block_factor_2(8),
    sparse_matrix_block_factor_1(32),
    threads(1),
    parallelLevel(1),
    /* FIXME: there should be a param to set no_of_buffers_per_allocator, for large calculations it needs to be larger, e.g. 10 x larger seems to give much better performance of matrix operations for large cases. 
       This is also connected to blocksize, maybe the best solution would be to have a param determining the number of MegaBytes per allocator or something like that.  */
    no_of_buffers_per_allocator(20000),
    use_allocator_manager(1)
  {};
  ~MatOptions() {
  }
  /** after the parameters are called, this routine is to be called
      to figure out the basis set permutation. */
  void prepare(const BasisInfoStruct& basisInfo);
};

struct OutputOptions {
  OutputOptions() 
  {}
    
};

} /* end of SCF name space */



#endif
