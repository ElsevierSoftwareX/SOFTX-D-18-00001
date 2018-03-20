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

/** @file grid_hicu.h

    @brief Hierarchical Cubature (HiCu) grid generation.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef _GRID_HICU_H_
#define _GRID_HICU_H_ 1

#include "grid_matrix.h"
#include "sparse_matrix.h"

int hicu_grid_generate(const char* grid_file_name,
                       const BasisInfoStruct& bis,
		       ergo_real maxError,
		       ergo_real boxSize,
		       ergo_real startBoxSizeDebug,
		       int use_error_per_volume,
		       int do_double_checking,
		       int compare_to_refined,
		       int use_energy_criterion,
		       int use_energy_criterion_only,
		       int do_variation_checking,
                       const Dft::Matrix* dmat,
                       Dft::SparsePattern* pattern,
                       int nThreads,
                       bool generateSparsePatternOnly);

void grid_generate_sparse_pattern(const BasisInfoStruct& bis,
				  ergo_real maxError,
				  ergo_real boxSize,
				  ergo_real startBoxSizeDebug,
				  Dft::SparsePattern& pattern);

#endif /* _GRID_HICU_H_ */
