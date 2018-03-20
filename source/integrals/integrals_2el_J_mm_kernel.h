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

/** @file integrals_2el_J_mm_kernel.h

    \brief Code for multipole method computational kernel for
    computing the Coulomb matrix J.

    @author: Elias Rudberg <em>responsible</em>.
*/

#ifndef INTEGRALS_2EL_J_MM_KERNEL_HEADER
#define INTEGRALS_2EL_J_MM_KERNEL_HEADER

#include "organize_distrs_mm.h"
#include "integrals_2el_utils.h"

int
do_multipole_interaction_between_2_boxes_branches(const distr_list_description_struct & distrDescription_1,
						  const multipole_struct_large & branchMultipole,
						  const multipole_struct_small* multipoleList_1,
						  ergo_real* result_J_list, // NULL if not used
						  ResultMatContrib* resultMatContrib, // NULL if not used
						  ergo_real threshold,
						  int* largest_L_used_so_far, // optional output, NULL if not used
						  MMInteractor & interactor,
						  const MMLimitTable & mmLimitTable
						  );

#endif
