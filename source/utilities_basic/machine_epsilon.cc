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

/** @file machine_epsilon.cc

    @brief Functionality for determining an approximate value of the
    "machine epsilon" -- the smallest number that can be added to 1
    giving a result different from 1, using the ergo_real
    floating-point type.

    @author: Elias Rudberg <em>responsible</em>
*/

#include "machine_epsilon.h"

/** return machine epsilon. This routine is deprecated. Use instead:

- std::numeric_limits<double>::epsilon() in C++;

- or FLT_EPSILON for single precision, DBL_EPSILON for double
  precision or or LDBL_EPSILON in long double - in C.
*/
ergo_real 
get_machine_epsilon()
{
  ergo_real volatile x, y, z;
  x = 1;
  y = 0.1;
  while(1)
    {
      z = x + y;
      if(z == x)
        return y;
      y *= 0.98;
    }
  /* this point should never be reached */
  return 0;
}

