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

/** @file diis_restricted.h

    @brief DIISManagerRestricted class implementing direct inversion
    in the iterative subspace (DIIS) for restricted SCF calculations.

    @author: Elias Rudberg <em>responsible</em>.
*/

#ifndef DIIS_RESTRICTED_HEADER
#define DIIS_RESTRICTED_HEADER

#include "diis_general.h"

class DIISManagerRestricted : public DIISManager
{
 public:
  DIISManagerRestricted();
  ~DIISManagerRestricted();

  int AddIterationToList(symmMatrix & F, normalMatrix & E);
  int ClearList();
  int GetCombinedFockMatrix(symmMatrix & result);
};


#endif
