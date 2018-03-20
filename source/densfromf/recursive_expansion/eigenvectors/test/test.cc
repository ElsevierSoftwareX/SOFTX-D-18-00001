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
#include <iostream>

#include "matrix_typedefs.h" // definitions of matrix types and interval type (source)
#include "realtype.h"        // definitions of types (utilities_basic)
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

#include "get_eigenvectors.h"

typedef ergo_real real;

#define SCALAR_TOL    template_blas_sqrt(mat::getMachineEpsilon<real>())


mat::SizesAndBlocks rows;
mat::SizesAndBlocks cols;

template<typename Matrix>
void init_matrix(Matrix& X, const int N, const int M)
{
   /********** Initialization of SizesAndBlocks */
   int nlevels = 5; //!!!

   std::vector<int> blockSizes(nlevels);
   blockSizes[nlevels - 1] = 1; // should always be one
   for (int ind = nlevels - 2; ind >= 0; ind--)
   {
      blockSizes[ind] = blockSizes[ind + 1] * 10;
   }
   rows = mat::SizesAndBlocks(blockSizes, N);
   cols = mat::SizesAndBlocks(blockSizes, M);
   /********************************************/
   X.resetSizesAndBlocks(rows, cols);
}


void get_matrix_from_full(std::vector<real> const& A, int N, int M, symmMatrix& X)
{
   init_matrix<symmMatrix>(X, N, M);
   assert(X.get_nrows() * X.get_ncols() == N * M);
   X.assignFromFull(A);
}


int main()
{
#ifdef _OPENMP
   int        defThreads;
   const char *env = getenv("OMP_NUM_THREADS");
   if (!(env && ((defThreads = atoi(env)) > 0)))
   {
      defThreads = 1;
   }
   mat::Params::setNProcs(defThreads);
   mat::Params::setMatrixParallelLevel(2);
#endif

   std::cout << "Start tests" << '\n';

   int n = 10;

   std::vector<real> vec(n *n);
   symmMatrix        B;
   get_matrix_from_full(vec, n, n, B);
   B.random();

   generalVector v;
   v.resetSizesAndBlocks(cols);

   std::cout << "--------------------------------" << '\n';
   std::cout << "Call computeEigenvectors " << '\n';

   std::cout << "Use Lanczos method" << '\n';

   real                       tol = 1e-12;
   std::vector<real>          eigVal(1, 0);
   std::vector<generalVector> eigVec(1, generalVector(rows));
   eigVec[0].rand();    // initial guess
   std::vector<int> num_iter(1, 0);
   real lambda1_lanczos, lambda1_power;

   try{
      eigvec::computeEigenvectors(B, tol, eigVal, eigVec, 1, "lanczos", num_iter);
      std::cout << "num_iter[0] = " << num_iter[0] << '\n';
      std::cout << "eigVal[0] = " << (double)eigVal[0] << '\n'; // (double) cast here to make it work for quad precision
      std::cout << "Call compute_rayleigh_quotient " << '\n';
      real lambda = eigvec::compute_rayleigh_quotient<real>(B, eigVec[0]);
      std::cout << "lambda = " << (double)lambda << '\n'; // (double) cast here to make it work for quad precision
      assert(template_blas_fabs(eigVal[0] - lambda) < SCALAR_TOL);
      lambda1_lanczos = lambda;
      assert(num_iter[0] <= 20); // avoid drastic growth of number of iterations
   }
   catch (std::exception& e)
   {
      std::cerr << e.what() << '\n';
   }

   std::cout << "Use power method" << '\n';
   std::vector<real>          eigValP(1, 0);
   std::vector<generalVector> eigVecP(1, generalVector(rows));
   eigVecP[0].rand();    // initial guess
   std::vector<int> num_iterP(1, 0);
   try{
      eigvec::computeEigenvectors(B, tol, eigValP, eigVecP, 1, "power", num_iterP);
      std::cout << "num_iter[0] = " << num_iterP[0] << '\n';
      std::cout << "eigVal[0] = " << (double)eigValP[0] << '\n'; // (double) cast here to make it work for quad precision
      std::cout << "Call compute_rayleigh_quotient " << '\n';
      real lambda = eigvec::compute_rayleigh_quotient<real>(B, eigVecP[0]);
      std::cout << "lambda = " << (double)lambda << '\n'; // (double) cast here to make it work for quad precision
      assert(template_blas_fabs(eigValP[0] - lambda) < SCALAR_TOL);
      lambda1_power = lambda;
      assert(num_iter[0] <= 30); // avoid drastic growth of number of iterations
   }
   catch (std::exception& e)
   {
      std::cerr << e.what() << '\n';
      std::cerr << "Ending with return -1 due to error." << std::endl;
      return -1;
   }

   // compare results obtained with power method and with Lanczos (largest eigenvalue)
   assert(template_blas_fabs(lambda1_lanczos - lambda1_power) < SCALAR_TOL);

   std::cout << "--------------------------------" << '\n';

   std::cout << "Call computeEigenvectors (deflation)" << '\n';

   try{
      eigVal.resize(2);
      std::vector<generalVector> eigVec2(2, generalVector(rows));
      eigVec2[0] = eigVec[0]; // copy computed eigenvector
      eigVec2[1].rand();      // initial guess
      num_iter.resize(2);
      eigvec::computeEigenvectors(B, tol, eigVal, eigVec2, 2, "lanczos", num_iter, 200, true);
      std::cout << "num_iter[0] = " << num_iter[0] << '\n';
      std::cout << "eigVal[0] = " << (double)eigVal[0] << '\n'; // (double) cast here to make it work for quad precision
      std::cout << "num_iter[1] = " << num_iter[1] << '\n';
      std::cout << "eigVal[1] = " << (double)eigVal[1] << '\n'; // (double) cast here to make it work for quad precision
      std::cout << "Call compute_rayleigh_quotient " << '\n';
      real lambda2 = eigvec::compute_rayleigh_quotient<real>(B, eigVec2[1]);
      std::cout << "lambda2 = " << (double)lambda2 << '\n'; // (double) cast here to make it work for quad precision
      
      assert(template_blas_fabs(lambda1_lanczos - eigVal[0]) < SCALAR_TOL);
      assert(lambda1_lanczos >= lambda2);
      assert(num_iter[0] <= 20); 
      assert(num_iter[1] <= 20); // avoid drastic growth of number of iterations
   }
   catch (std::exception& e)
   {
      std::cerr << e.what() << '\n';
      std::cerr << "Ending with return -1 due to error." << std::endl;
      return -1;
   }

   std::cout << "--------------------------------" << '\n';

   std::cout << "End tests" << '\n';
   return 0;
}
