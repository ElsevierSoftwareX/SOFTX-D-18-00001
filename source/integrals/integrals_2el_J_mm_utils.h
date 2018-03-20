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

/** @file integrals_2el_J_mm_utils.h

    \brief Utility functions related to multipole method, used in
    construction of the Coulomb matrix J.

    @author: Elias Rudberg <em>responsible</em>.
*/

#ifndef INTEGRALS_2EL_J_MM_UTILS_HEADER
#define INTEGRALS_2EL_J_MM_UTILS_HEADER

#include "integral_info.h"
#include "organize_distrs.h"
#include "organize_distrs_mm.h"

int
check_if_multipoles_can_be_used(const IntegralInfo & integralInfo,
				ergo_real threshold,
				const ergo_real* boxCenterCoords_1,
				const ergo_real* boxCenterCoords_2,
				ergo_real boxWidth,
				const distr_org_struct & org_1,
				const distr_org_mm_struct & org_mm_1,
				const distr_org_struct & org_2,
				const distr_org_mm_struct & org_mm_2);

int
create_list_of_multipoles_for_box(const IntegralInfo& integralInfo,
				  const distr_org_struct & org,
				  multipole_struct_small* result_multipoleList);

#endif
