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

/** @file multipole.h

    @brief Code for computing multipole moments, and multipole
    interaction and translation matrices.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef MULTIPOLE_HEADER
#define MULTIPOLE_HEADER

#include "realtype.h"
#include "integral_info.h"
#include "basisinfo.h"
#include "multipole_prep.h"

int
compute_multipole_moments(const IntegralInfo& integralInfo,
			  const DistributionSpecStruct* distr,
			  multipole_struct_small* result);

class MMTranslator {
  static const int MMDP1 = MAX_MULTIPOLE_DEGREE+1;
  ergo_real *buffer_W_cc;
  ergo_real *buffer_W_cs;
  ergo_real *buffer_W_sc;
  ergo_real *buffer_W_ss;
  const MultipolePrepManager & multipolePrep;
 public:
  MMTranslator(const MultipolePrepManager & multipolePrepManager);
  ~MMTranslator();
  int getTranslationMatrix(ergo_real dx,
                           ergo_real dy,
                           ergo_real dz,
                           int l_1,
                           int l_2,
                           ergo_real* result_W) const;
};

class MMInteractor {
  static const int MMDP1 = MAX_MULTIPOLE_DEGREE+1;
  ergo_real *buffer_T_cc;
  ergo_real *buffer_T_cs;
  ergo_real *buffer_T_sc;
  ergo_real *buffer_T_ss;
  const MultipolePrepManager & multipolePrep;
 public:
  MMInteractor(const MultipolePrepManager & multipolePrepManager);
  ~MMInteractor();
  int getInteractionMatrix(ergo_real dx,
                           ergo_real dy,
                           ergo_real dz,
                           int l_1,
                           int l_2,
                           ergo_real* result_T);
};

int
setup_multipole_maxAbsMomentList(multipole_struct_large* multipole);

#endif
