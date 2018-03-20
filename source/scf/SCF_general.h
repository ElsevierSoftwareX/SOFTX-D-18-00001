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

/** @file SCF_general.h

    @brief Class for self-consistent field (SCF) procedure; base class
    that can be used for both restricted and unrestricted cases.

    @author: Elias Rudberg <em>responsible</em>.
*/

#ifndef SCF_GENERAL_HEADER
#define SCF_GENERAL_HEADER

#include "molecule.h"
#include "basisinfo.h"
#include "integrals_2el.h"
#include "grid_stream.h"
#include "scf.h"
#include "diis_general.h"
#include "SCF_statistics.h"
#include "GetDensFromFock.h"


class SCF_general
{
 public:  

  // SCF convergence routine
  void do_SCF_iterations();

  void get_overlap_matrix(symmMatrix & S);
  void get_invCholFactor_matrix(triangMatrix & invCholFactor_);
  void get_H_core_matrix(symmMatrix & H_core);
  void get_energy(ergo_real & E, ergo_real & E_nuclear);

 protected:
  // Constructor
  SCF_general(const Molecule& molecule_, 
	      const Molecule& extraCharges_,
	      const BasisInfoStruct & basisInfo_, 
	      const IntegralInfo & integralInfo_,
	      const char* guessDmatFileName_,
	      const JK::Params& J_K_params_,
	      const Dft::GridParams& gridParams_,
	      const SCF::Options& scfopts,
	      const SCF::MatOptions& matOpts,
	      ergo_real threshold_integrals_1el_input);

  // Destructor
  virtual ~SCF_general();

  const Molecule& molecule;
  const Molecule& extraCharges;
  const BasisInfoStruct & basisInfo;
  const IntegralInfo& integralInfo;
  const char* guessDmatFileName;
  const JK::Params& J_K_params;
  const Dft::GridParams& gridParams;
  const SCF::Options& scfopts;
  const SCF::MatOptions& matOpts;
  ergo_real threshold_integrals_1el;

  //integral_prep_struct* integralPrep;

  JK::ExchWeights CAM_params;  // range-separated exchange parameters

  // number of SCF cycle
  int SCF_step;

  // nuclearEnergy is nuclear repulsion energy plus contribution from external electric field.
  ergo_real nuclearEnergy;

  ergo_real energy_2el;
  ergo_real energy;

  ergo_real energy_2el_core; // only used when "core density matrix" is used
  ergo_real energy_2el_valence; // only used when "core density matrix" is used
  ergo_real energy_of_valence; // only used when "core density matrix" is used
  ergo_real energy_reference; // only used when "core density matrix" is used

  ergo_real electronicEntropyTerm;

  ergo_real errorMeasure;

  ergo_real curr_subspace_diff;

  symmMatrix S_symm;
  triangMatrix invCholFactor;
  ergo_real invCholFactor_euclnorm;
  symmMatrix H_core_Matrix;

  DIISManager* DIIS; // Must be initialized by restricted/unrestricted derived class.

  int noOfElectrons;

  SCF_statistics* curr_cycle_stats;

  GetDensFromFock DensFromFock;

  ergo_real GetEuclideanNormOfMatrix(const symmMatrix & A);

  virtual void initialize_matrices() = 0;
  virtual void check_params() = 0;
  virtual void get_starting_guess_density() = 0;
  virtual void initialize_homo_lumo_limits() = 0;
  virtual void write_matrices_to_file() = 0;
  virtual void get_2e_part_and_energy() = 0;
  virtual void output_sparsity_S_F_D(SCF_statistics & stats) = 0;
  virtual void calculate_energy() = 0;
  virtual void get_FDSminusSDF() = 0;
  virtual void get_error_measure() = 0;
  virtual void add_to_DIIS_list() = 0;
  virtual void update_best_fock_so_far() = 0;
  virtual void combine_old_fock_matrices(ergo_real stepLength) = 0;
  virtual void use_diis_to_get_new_fock_matrix() = 0;
  virtual void clear_diis_list() = 0;
  virtual void clear_error_matrices() = 0;
  virtual void save_current_fock_as_fprev() = 0;
  virtual void get_new_density_matrix() = 0;
  virtual void write_density_to_file() = 0;
  virtual void save_final_potential() = 0;
  virtual void add_random_disturbance_to_starting_guess() = 0;
  virtual void output_expected_values_pos_operator() = 0;
  virtual void output_density_images() = 0;
  virtual void write_diag_dens_to_file() = 0;
  virtual void report_final_results() = 0;
  virtual void save_density_as_prevdens() = 0;
  virtual void update_subspace_diff() = 0;
  virtual void disturb_fock_matrix(ergo_real subspaceError) = 0;
  virtual void disturb_dens_matrix(ergo_real subspaceError) = 0;
  virtual void do_spin_flip(int atomCount) = 0;
  virtual void disturb_dens_matrix_exact(ergo_real subspaceError) = 0;
  virtual void save_full_matrices_for_matlab() = 0;
  virtual void report_density_difference() = 0;
  virtual void create_mtx_files_F(int const scfIter) = 0;
  virtual void create_mtx_files_D(int const scfIter) = 0;
  virtual void create_eigenvectors_files() const = 0;
  virtual void create_gabedit_file() const = 0;
  virtual void compute_dipole_moment() = 0;
  virtual void do_mulliken_pop_stuff() = 0;
  virtual void compute_gradient_fixeddens() = 0;
};





#endif
