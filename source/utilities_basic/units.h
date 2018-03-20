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

/** @file units.h

    @brief Constants for conversion between units for some common
    units like Angstrom, electron-volt (eV), Kelvin etc.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef UNITS_HEADER
#define UNITS_HEADER

#include "realtype.h"

#define UNIT_one_Angstrom ((ergo_real)(1.0 / 0.5291772083))
#define UNIT_one_Debye    0.393456
#define UNIT_one_eV       0.036749322    // According to nist.gov 2014 CODATA values
#define UNIT_one_Kelvin   3.16681046e-06 // Computed from nist.gov 2014 CODATA values

#endif
