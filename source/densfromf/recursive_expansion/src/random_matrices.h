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

/** @file random_matrices.h

    @brief Header file containing declarations of functions required for testing purposes.
    Functions include generation of the random dense matrices, random sparse symmetric matrices, initialization of the 
    hierarchical matrix structure, work with files and printing matrix to the screen.

    @author Anastasia Kruchinina 
    @sa random_matrices.cc
*/


#include "matrix_typedefs.h" // definitions of matrix types and interval type (source)
#include "realtype.h"   // definitions of types (utilities_basic)
#include "matrix_utilities.h"
#include "integral_matrix_wrappers.h"
#include "SizesAndBlocks.h"
#include "Matrix.h"
#include "Vector.h"
#include "MatrixSymmetric.h"
#include "MatrixTriangular.h"
#include "MatrixGeneral.h"
#include "VectorGeneral.h"
#include "output.h"

#include "files_dense.h"
#include "files_sparse.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>


using namespace std;



typedef intervalType IntervalType;
typedef symmMatrix MatrixTypeInner;
typedef triangMatrix TriangMatrixType;
typedef normalMatrix MatrixGeneral;

typedef std::vector<int> VectorTypeInt;

#define MAX_DOUBLE std::numeric_limits<real>::max()
#define MIN_DOUBLE std::numeric_limits<real>::min()


#define PI 3.14159265 // needed for sprandsym


void print_matrix(std::vector<ergo_real> const &A);
template<typename Matrix>
void init_matrix(Matrix &X, const int N, int blockSizesMultuple = 4);
void get_random_matrix(int N, MatrixTypeInner &X);
void get_all_eigenvalues_of_matrix(std::vector<ergo_real> & eigvalList, const MatrixTypeInner & M);
void sprandsym(int N, MatrixTypeInner &X, MatrixGeneral &Q, vector<ergo_real> &D, const double MATRIX_SPARSITY);
int get_matrix_from_sparse(char *filename, MatrixTypeInner &X);
int get_matrix_from_sparse_vec(char *filename, std::vector<int> &I, std::vector<int> &J, std::vector<real> &val);
int get_matrix_from_binary(char *filename, MatrixTypeInner &X);
int get_matrix_from_binary_vec(char *filename, std::vector<int> &I, std::vector<int> &J, std::vector<real> &val, int &N);
int get_matrix_from_full(char * filename, MatrixTypeInner &X);


/**
 *  \brief Create hierarchical matrix structure.
 *
 *  \tparam Matrix    type of the matrix (ex. symmMatrix)
 */
template<typename Matrix>
void init_matrix(Matrix &X, const int N, int blockSizesMultuple /*=4*/)
{
  /********** Initialization of SizesAndBlocks */
  int size = N;
  int nlevels = 5; //!!!
  std::vector<int> blockSizes(nlevels);
  blockSizes[nlevels - 1] = 1; // should always be one
    for (int ind = nlevels - 2; ind >= 0; ind--)
      blockSizes[ind] = blockSizes[ind + 1] * blockSizesMultuple;
  mat::SizesAndBlocks rows(blockSizes, size);
  mat::SizesAndBlocks cols(blockSizes, size);
  /********************************************/  
  X.resetSizesAndBlocks(rows,cols);
}
