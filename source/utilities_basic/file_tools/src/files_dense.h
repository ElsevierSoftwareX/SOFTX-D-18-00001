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

/** @file files_dense.h

    @brief File containing declaration of functions for reading/writing dense matrices and vectors.

    @author Anastasia Kruchinina <em>responsible</em>
*/


#ifndef FILES_DENSE_HEADER_TXT
#define FILES_DENSE_HEADER_TXT

#include "realtype.h"   // definitions of types (utilities_basic)

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string.h>


using namespace std;

typedef ergo_real real;
#define PRECISION 16



/* READ DENSE MATRIX FROM THE FILE */
int read_matrix(const char *filename, vector<real> &A, int &N, int &M, bool is_binary = false);

int read_vector(const char *filename, vector<real> &A, int &N, bool is_binary);

/* WRITE DENSE MATRIX TO THE FILE */
int write_matrix(const char *filename, const vector<real> &A, int N, int M);

/* WRITE VECTOR TO THE FILE */
int write_vector(const char *filename, const vector<real> &v);



#endif //FILES_DENSE_HEADER_TXT




