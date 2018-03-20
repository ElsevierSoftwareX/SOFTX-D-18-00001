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

/** @file test_LanczosSeveralLargestEig.cc

    \brief Code for testing functionality for somputing several
    eigenpairs using the Lanczos method.

    @author Anastasia Kruchinina
*/

#include <fstream>  /* For ifstream */
#include <iomanip> /* For setprecision in fstream */
#include <iostream>
#include <cmath>
#include <stdio.h> /* For FILE */
#include "SizesAndBlocks.h"
#include "Matrix.h"
#include "Vector.h"
#include "MatrixSymmetric.h"
#include "MatrixTriangular.h"
#include "MatrixGeneral.h"
#include "VectorGeneral.h"
#include "mat_gblas.h"
#include "Lanczos.h"
#include "LanczosSeveralLargestEig.h"

using namespace mat;


typedef double real;

// define matrix and vector hierarchy 
typedef Matrix<real, real> Mat_1;
typedef Matrix<real, Mat_1> Mat_2;
typedef Matrix<real, Mat_2> Mat_3;
typedef Vector<real, real > Vec_1;
typedef Vector<real, Vec_1> Vec_2;
typedef Vector<real, Vec_2> Vec_3;
  
typedef Mat_3 matri;
typedef Vec_3 vect;
typedef MatrixSymmetric<real, matri> symmMatrix;
typedef MatrixTriangular<real, matri> triangMatrix;
typedef MatrixGeneral<real, matri> normalMatrix;
typedef VectorGeneral<real, vect> normalVector;


int main()
{

  srand(time(NULL));
  
#ifdef _OPENMP	
  int defThreads;
  const char *env = getenv("OMP_NUM_THREADS");
  if ( !(env && (defThreads=atoi(env)) > 0) ) {
    defThreads = 1;
  }
  
  mat::Params::setNProcs(defThreads);
  mat::Params::setMatrixParallelLevel(2);
  std::cout<<"OpenMP is used, number of threads set to "
	   <<mat::Params::getNProcs()<<". Matrix parallel level: "
	   <<mat::Params::getMatrixParallelLevel()<<"."<<std::endl;
#endif

  try
    {
      typedef arn::LanczosSeveralLargestEig<real, symmMatrix, normalVector> myLanczosType;
      real epsilon = template_blas_sqrt(mat::getRelPrecision<real>());
  

  /********** Initialization of SizesAndBlocks                            */
  int size = 51; /* Use weird size to find more bugs. */
  int nlevels = 3;
  std::vector<int> blockSizes(nlevels);
  blockSizes[nlevels - 1] = 1; // should always be one
#if 1
  blockSizes[nlevels - 2] = 1; // lowest level blocksize
  blockSizes[nlevels - 3] = 5;
#else
  for (int ind = nlevels - 2; ind >= 0; ind--)
    blockSizes[ind] = blockSizes[ind + 1] * 10;
#endif

  std::cout << "Running tests with blocksize vector: ";
  for (int ind = 0; ind < nlevels; ind++)
    std::cout << blockSizes[ind] << "  ";
  std::cout << std::endl;

  SizesAndBlocks rows(blockSizes, size);
  SizesAndBlocks cols(blockSizes, size);
    
  real ONEreal = 1.0;

  symmMatrix syA;
  syA.resetSizesAndBlocks(rows,cols);
  syA.randomZeroStructure(0.3);
  //syA.random();
  // std::vector<int> rowsA;
  // std::vector<int> colsA;
  // std::vector<real> valsA;
  // syA.get_all_values(rowsA, colsA, valsA);
  // std::cout << "Matrix:" << std::endl;
  // std::cout << rowsA.size() << std::endl;
  // for(int i = 0; i < rowsA.size(); ++i)
  //   std::cout << rowsA[i] << " " << colsA[i] << " " << valsA[i] << std::endl;

  normalVector x;
  x.resetSizesAndBlocks(rows);
  x.rand();
  int maxit = 400;

  myLanczosType lan(syA, x, 5, maxit);
  real lanEpsilon = epsilon;
  lan.setRelTol( lanEpsilon );
  lan.run();
  normalVector eigVec;
  real eigVal;
  real accuracy;

  lan.get_ith_eigenpair(1, eigVal, eigVec, accuracy);
  normalVector resVec(eigVec); // residual
  resVec *= eigVal;
  resVec += -ONEreal * syA * eigVec;

  normalVector RelResVec(eigVec); // relative residual
  RelResVec *= eigVal;
  RelResVec += -ONEreal * syA * eigVec;
  RelResVec *= 1.0 / eigVal;
  
  std::cout<<"\nLanczos several largest magnitude test : \n\n"
	   << "FIRST EIGENPAIR: \n"
	   << "Eigenvalue: " << std::setprecision(12) << eigVal <<std::setw(15)
	   <<"\n Requested relative accuracy: "
	   <<std::setprecision(10)<<std::setw(15)
	   <<lanEpsilon
    	   <<"\n Relative residual norm:       "
	   <<std::setprecision(10)<<std::setw(15)
	   <<RelResVec.eucl() 
	   <<"\n Indicated absolute error (max norm of the residual): "
	   <<std::setprecision(10)<<std::setw(15)
	   <<accuracy
	   <<"\n Residual:       "
	   <<std::setprecision(10)<<std::setw(15)
	   <<resVec.eucl() << std::endl;
  if (RelResVec.eucl() < lanEpsilon && 
      (resVec.eucl() < accuracy || 
       resVec.eucl() < mat::getRelPrecision<real>() * 100))
    std::cout<<"   OK" <<std::endl;
  else {
    std::cout<<"   ERROR" <<std::endl;
    std::exit(1);
  }    


  lan.get_ith_eigenpair(2, eigVal, eigVec, accuracy);
  resVec = eigVec; // residual
  resVec *= eigVal;
  resVec += -ONEreal * syA * eigVec;

  RelResVec = eigVec; // relative residual
  RelResVec *= eigVal;
  RelResVec += -ONEreal * syA * eigVec;
  RelResVec *= 1.0 / eigVal;

  std::cout<< "SECOND EIGENPAIR: \n"
	   << "Eigenvalue: " << std::setprecision(12) << eigVal <<std::setw(15)
	   <<"\n Requested relative accuracy: "
	   <<std::setprecision(10)<<std::setw(15)
	   <<lanEpsilon
    	   <<"\n Relative residual norm:       "
	   <<std::setprecision(10)<<std::setw(15)
	   <<RelResVec.eucl() 
	   <<"\n Indicated absolute error (max norm of the residual): "
	   <<std::setprecision(10)<<std::setw(15)
	   <<accuracy
	   <<"\n Residual:       "
	   <<std::setprecision(10)<<std::setw(15)
	   <<resVec.eucl() << std::endl;
  if (RelResVec.eucl() < lanEpsilon && 
      (resVec.eucl() < accuracy || 
       resVec.eucl() < mat::getRelPrecision<real>() * 100))
    std::cout<<"   OK" <<std::endl;
  else {
    std::cout<<"   ERROR" <<std::endl;
    std::exit(1);
  }   


  lan.get_ith_eigenpair(5, eigVal, eigVec, accuracy);
  resVec = eigVec; // residual
  resVec *= eigVal;
  resVec += -ONEreal * syA * eigVec;

  RelResVec = eigVec; // relative residual
  RelResVec *= eigVal;
  RelResVec += -ONEreal * syA * eigVec;
  RelResVec *= 1.0 / eigVal;
  
  std::cout<< "FIFTH EIGENPAIR: \n"
	   << "Eigenvalue: " << std::setprecision(12) << eigVal <<std::setw(15)
	   <<"\n Requested relative accuracy: "
	   <<std::setprecision(10)<<std::setw(15)
	   <<lanEpsilon
    	   <<"\n Relative residual norm:       "
	   <<std::setprecision(10)<<std::setw(15)
	   <<RelResVec.eucl() 
	   <<"\n Indicated absolute error (max norm of the residual): "
	   <<std::setprecision(10)<<std::setw(15)
	   <<accuracy
	   <<"\n Residual:       "
	   <<std::setprecision(10)<<std::setw(15)
	   <<resVec.eucl() << std::endl;
    if (RelResVec.eucl() < lanEpsilon && 
      (resVec.eucl() < accuracy || 
       resVec.eucl() < mat::getRelPrecision<real>() * 100))
    std::cout<<"   OK" <<std::endl;
  else {
    std::cout<<"   ERROR" <<std::endl;
    std::exit(1);
  }    

}
  catch (std::exception & e) {
  std::cout << "Exception caught: "<<e.what() << std::endl;
  std::exit(1);
}

    return 0;

}
