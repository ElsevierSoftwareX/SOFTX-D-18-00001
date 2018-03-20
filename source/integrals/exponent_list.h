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

/** @file exponent_list.h

    @brief ExponentList class for keeping track of a set of unique
    exponents in Gaussian functions, together with a max coefficient
    value for each exponent.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef EXPONENT_LIST_HEADER
#define EXPONENT_LIST_HEADER

#include "realtype.h"
#include "basisinfo.h"

const int MAX_NO_OF_UNIQUE_EXPONENTS = 222;
const ergo_real CONST_EXPONENT_DIFF_TOLERANCE = 0.0001;

typedef struct
{
  ergo_real exponent;
  ergo_real maxAbsCoeff;
} unique_exponent_struct;


class ExponentList
{
 public:
  int noOfExponents;
  unique_exponent_struct list[MAX_NO_OF_UNIQUE_EXPONENTS];
  int get_list_of_available_exponents(const BasisInfoStruct & basisInfo);
};


#endif
