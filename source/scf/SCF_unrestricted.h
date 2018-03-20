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

/** @file SCF_unrestricted.h

    @brief Class for self-consistent field (SCF) procedure;
    spin-unrestricted case.

    @author: Elias Rudberg <em>responsible</em>.
*/

#ifndef SCF_UNRESTRICTED_HEADER
#define SCF_UNRESTRICTED_HEADER

#include "SCF_general.h"


class SCF_unrestricted : public SCF_general
{
 public:  

  // Constructor
  SCF_unrestricted(const Molecule& molecule_, 
		   const Molecule& extraCharges_,
		   const BasisInfoStruct & basisInfo_, 
		   const IntegralInfo & integralInfo_,
		   const char* guessDmatFileName_,
		   const JK::Params& J_K_params_,
		   const Dft::GridParams& gridParams_,
		   const SCF::Options& scfopts,
		   const SCF::MatOptions& matOpts,
		   ergo_real threshold_integrals_1el_input,
		   int alpha_beta_diff_input);

  // Destructor
  ~SCF_unrestricted();

  void get_Fock_matrices(symmMatrix & FockMatrix_a, symmMatrix & FockMatrix_b);
  void get_no_of_electrons(int & noOfElectrons_a, int & noOfElectrons_b);

 private:
  void initialize_matrices();
  void check_params();
  void get_starting_guess_density();
  void initialize_homo_lumo_limits();
  void write_matrices_to_file();
  void get_2e_part_and_energy();
  void output_sparsity_S_F_D(SCF_statistics & stats);
  void calculate_energy();
  void get_FDSminusSDF();
  void get_error_measure();
  void add_to_DIIS_list();
  void update_best_fock_so_far();
  void combine_old_fock_matrices(ergo_real stepLength);
  void use_diis_to_get_new_fock_matrix();
  void clear_diis_list();
  void clear_error_matrices();
  void save_current_fock_as_fprev();
  void get_new_density_matrix();
  void write_density_to_file();
  void save_final_potential();
  void add_random_disturbance_to_starting_guess();
  void output_expected_values_pos_operator();
  void output_density_images();
  void output_density_images_orbital(generalVector &eigVec, const std::string &filename_id);
  void write_diag_dens_to_file();
  void report_final_results();
  void save_density_as_prevdens();
  void update_subspace_diff();
  void disturb_fock_matrix(ergo_real subspaceError);
  void disturb_dens_matrix(ergo_real subspaceError);
  void do_spin_flip(int atomCount);
  void disturb_dens_matrix_exact(ergo_real subspaceError);
  void save_full_matrices_for_matlab();
  void report_density_difference();
  void create_mtx_files_F(int const scfIter);
  void create_mtx_files_D(int const scfIter);
  /* void create_homo_eigvec_file() const; */
  /* void create_lumo_eigvec_file() const; */
  void create_eigenvectors_files() const;
  void create_eigvec_file(const generalVector &eigVec_alpha, 
			  const generalVector &eigVec_beta,
			  const char *vector_name, 
			  const char *filename_id) const;
  void create_gabedit_file() const;
  void create_gabedit_file_2() const;
  void compute_dipole_moment();
  void do_mulliken_pop_stuff();
  void compute_gradient_fixeddens();

  void get_S2(ergo_real & S2_exact, ergo_real & S2);

  symmMatrix densityMatrix_alpha;
  symmMatrix densityMatrix_beta;
  symmMatrix FockMatrix_alpha;
  symmMatrix FockMatrix_beta;
  symmMatrix Fprev_alpha;
  symmMatrix Fprev_beta;
  symmMatrix Dprev_alpha;
  symmMatrix Dprev_beta;
  symmMatrix F_ort_prev_alpha; // Used by purification
  symmMatrix F_ort_prev_beta; // Used by purification
  symmMatrix D_ort_prev_alpha; // Used for computing eigenvectors
  symmMatrix D_ort_prev_beta; // Used for computing eigenvectors
  symmMatrix bestFockMatrixSoFar_alpha;
  symmMatrix bestFockMatrixSoFar_beta;
  symmMatrix bestFockMatrixSoFar2_alpha;
  symmMatrix bestFockMatrixSoFar2_beta;
  normalMatrix ErrorMatrix_alpha;
  normalMatrix ErrorMatrix_beta;
  symmMatrix G_alpha;
  symmMatrix G_beta;

  // HOMO/LUMO info
  intervalType homoInterval_F_ort_prev_alpha;
  intervalType lumoInterval_F_ort_prev_alpha;
  intervalType homoInterval_F_ort_prev_beta;
  intervalType lumoInterval_F_ort_prev_beta;
  intervalType homoInterval_Fprev_alpha;
  intervalType lumoInterval_Fprev_alpha;
  intervalType homoInterval_Fprev_beta;
  intervalType lumoInterval_Fprev_beta;

  int alpha_beta_diff;
  int noOfElectrons_alpha;
  int noOfElectrons_beta;

  generalVector eigVecLUMO_alpha;
  generalVector eigVecHOMO_alpha;
  generalVector eigVecLUMO_beta;
  generalVector eigVecHOMO_beta;
};





#endif
