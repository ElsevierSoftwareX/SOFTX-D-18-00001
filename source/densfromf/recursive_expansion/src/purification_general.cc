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

/** @file purification_general.cc

    @brief File contataining definitions of various variable used in the recursive density matrix expansion 
    (or density matrix purification).

    @author Anastasia Kruchinina <em>responsible</em>
*/


#include "purification_general.h"


real eucl_acc = 1e3*mat::getMachineEpsilon<real>();   /**< Tolerance used for computation of spectral norm.  */
real mixed_acc = 1e3*mat::getMachineEpsilon<real>();  /**< Tolerance used for computation of mixed norm. NOTE: If truncation is 0 this may not be enough, set to machine epsilon. */

real TOL_OVERLAPPING_BOUNDS = 1e-10;  /**< If the difference between inner bounds for homo and lumo 
					 is less then this tolerance, i.e. bounds are still bad, eigenvectors will not be computed.
					 (Inner bounds are used to estimate iterations for computation of homo and lumo eigenvectors.)  */ 

real THRESHOLD_EIG_TOLERANCE = 1e-5;  /**< Inner homo and lumo bounds may be too good, and it may happen that computed eigenvalue slightly outside of given intervals.
					 Thus we allow some flexibility for eigenvalue.
					 Set threshold 1e-5 since otherwise for small molecules does not work. */

int EIG_EMPTY = 0;
int EIG_SQUARE_INT = 1;
int EIG_PROJECTION_INT = 2;
int EIG_POWER_INT = 3;
int EIG_LANCZOS_INT = 4;



