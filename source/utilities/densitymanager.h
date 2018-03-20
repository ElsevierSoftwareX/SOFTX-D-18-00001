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

/** @file densitymanager.h

    @brief Functionality for working with the electron density as a
    function of space, for a given basis set and density matrix.

    @author: Elias Rudberg <em>responsible</em>
*/

#include "basisinfo.h"

/* 
get_no_of_primitives_for_density is a helper function for
get_density. Call get_no_of_primitives_for_density to find out
how long the result list needs to be.
A negative return value indicates failure.
*/
int get_no_of_primitives_for_density(ergo_real cutoff,
				     const ergo_real *dmat,
				     const BasisInfoStruct & basisInfo);

/* 
get_density creates the list resultRho using information from
basisInfo and dmat, using given threshold.
A negative return value indicates failure.
*/
int get_density(const BasisInfoStruct & basisInfo,
		const ergo_real* dmat, /* density matrix */
		ergo_real cutoff, /* threshold */
		int maxCountRho, /* maxcount for result list */
		DistributionSpecStruct* resultRho);

ergo_real integrate_density_in_box(int nPrims,
				   DistributionSpecStruct* rho,
				   ergo_real mid_x,
				   ergo_real mid_y,
				   ergo_real mid_z,
				   ergo_real box_width);

ergo_real integrate_density_in_box_2(int nPrims,
				     DistributionSpecStruct* rho,
				     ergo_real* minVect, 
				     ergo_real* maxVect, 
				     std::vector<int> monomialIntsAdd = std::vector<int>(3, 0));

