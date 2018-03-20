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

/** @file files_sparse.h

    @brief File containing declarations of functions for reading/writing sparse matrices from/to mtx (MatrixMarket format) files.

    @author Anastasia Kruchinina <em>responsible</em>
*/


#ifndef FILES_SPARSE_HEADER
#define FILES_SPARSE_HEADER

#include "realtype.h"   // definitions of types (utilities_basic)

#include <iostream>
#include <fstream>
#include <string.h>
#include <assert.h>
#include <vector>
#include <sstream>
#include <stdexcept>
#include "mmio.h"


using namespace std;

typedef ergo_real real;


/* READ SPARSE MATRIX FROM THE MATRIX MARKET FILE */
int read_matrix_from_mtx(const char* filename, vector<int> &I, vector<int> &J, vector<real> &val, int &N, int &M);


/* WRITE SPARSE MATRIX TO THE MATRIX MARKET FILE */
int write_matrix_to_mtx(const char* filename, const vector<int> &I, const vector<int> &J, const vector<real> &val, const int &N);
int write_matrix_to_mtx_nonsymm(const char* filename, const vector<int> &I, const vector<int> &J, const vector<real> &val, const int &N, const int &M);

#endif //FILES_SPARSE_HEADER




