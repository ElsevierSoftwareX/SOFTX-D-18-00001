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

/** @file integrals_2el_repeating.h

    @brief Functionality for keeping track of certain kinds of
    integrals that are computed repeatedly, saving the computed values
    instead of recomputing them.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef INTEGRALS_2EL_REPEATING
#define INTEGRALS_2EL_REPEATING


#include "basisinfo.h"

class IntegratorCase
{
 public:
  IntegratorCase(int Nmax_in, 
		 int noOfMonomials, 
		 ergo_real exponent_in, 
		 const ergo_real* newList);
  ~IntegratorCase();
  int Nmax;
  ergo_real exponent;
  ergo_real* list;
};

class IntegratorWithMemory
{
 public:
  explicit IntegratorWithMemory(const IntegralInfo* b);
  ~IntegratorWithMemory();
  
  ergo_real do_2e_integral(const DistributionSpecStruct* psi);

 private:
  const IntegralInfo* integralInfo;
  int noOfCases;
  IntegratorCase** caseList;
};



#endif
