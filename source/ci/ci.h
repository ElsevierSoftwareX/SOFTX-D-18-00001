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

/** @file ci.h

    @brief Configuration Interaction (CI) code.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef CI_HEADER
#define CI_HEADER

#include "basisinfo.h"
#include "integrals_2el.h"


namespace CI {

struct Options {
  int use_random_orbitals;
  int use_lowdin_orbitals;
  int no_of_core_electrons;
  int use_random_starting_guess;
  ergo_real convergence_threshold;
  ergo_real initial_step_length;
  int max_no_of_iterations;
  ergo_real shift;
  int use_energy_diff_limit;
  ergo_real energy_diff_limit;
  
  /** Initializes all the fields to sane values. */
  Options() : use_random_orbitals(0),
              use_lowdin_orbitals(0),
	      no_of_core_electrons(0),
	      use_random_starting_guess(0),
	      convergence_threshold(1e-4),
	      initial_step_length(0.01),
	      max_no_of_iterations(100),
	      shift(0.0),
	      use_energy_diff_limit(0),
	      energy_diff_limit(10.0)
  {
  }
};

} /* End of CI namespace */


int do_CI(
	  const BasisInfoStruct & basisInfo, 
	  const IntegralInfo & integralInfo,
	  const CI::Options& options,
	  const Molecule & molecule,
	  const ergo_real* S,
	  const ergo_real* h_AO,
	  const ergo_real* F_a,
	  const ergo_real* F_b,
	  int n_el_a,
	  int n_el_b,
	  ergo_real nuclearEnergy,
	  ergo_real HF_energy
	  );


#endif
