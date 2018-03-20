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

/** @file cubature_rules.h

    @brief Code for various cubature rules that can be used for
    numerical integration.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef DFT_CUBATURE_RULES_HEADER
#define DFT_CUBATURE_RULES_HEADER

#include "realtype.h"

typedef ergo_real real;

const int NO_OF_DIMENSIONS = 3;  /* Must be 3, really. */

struct BoxStruct_{
    real min[NO_OF_DIMENSIONS]; /* xmin, ymin, zmin */
    real max[NO_OF_DIMENSIONS]; /* xmax, ymax, zmax */
};
typedef struct BoxStruct_ BoxStruct;


int
use_cubature_rule(int maxlen,
		  real (*coor)[3],
		  real *weight,
		  BoxStruct* box,
		  int ruleNumber);

#endif
