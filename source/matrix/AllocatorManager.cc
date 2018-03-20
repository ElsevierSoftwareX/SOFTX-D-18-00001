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
#include "AllocatorManager.h"

/** @file AllocatorManager.cc

    \brief Code for AllocatorManager that is used to allocate memory
    for submatrices in the hierarchical matrix library.
*/

/* We need to include config.h to get macro PRECISION_QUAD_FLT128. */
#include "config.h"

namespace mat {

  template<>
  AllocatorManager<float> & AllocatorManager<float>::instance() {
    static AllocatorManager<float> theInstance;
    return theInstance;
  }

  template<>
  AllocatorManager<double> & AllocatorManager<double>::instance() {
    static AllocatorManager<double> theInstance;
    return theInstance;
  }

  template<>
  AllocatorManager<long double> & AllocatorManager<long double>::instance() {
    static AllocatorManager<long double> theInstance;
    return theInstance;
  }

#ifdef PRECISION_QUAD_FLT128
  template<>
  AllocatorManager<__float128> & AllocatorManager<__float128>::instance() {
    static AllocatorManager<__float128> theInstance;
    return theInstance;
  }
#endif

} /* end namespace mat */
