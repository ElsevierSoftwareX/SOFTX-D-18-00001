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

/** @file Perturb_Test.cc Test of the denaity matrix perturbation iterations
 *
 * Copyright(c) Emanuel Rubensson 2008
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date June 2008
 *
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
#include "Perturbation.h"

template<typename Treal>
class expRule {
public:
  Treal set(int const row, int const col) {
    return (rand() / (Treal)RAND_MAX)  * 
      template_blas_exp(-(template_blas_fabs(Treal(row)-Treal(col))));
  } 
};

template<typename Treal>
class setFromFullRule {
public:
  Treal* fMat;
  int const n;
  setFromFullRule(Treal* fullMat, int const n_in) :fMat(fullMat), n(n_in) {}
  Treal set(int const row, int const col) {
    return fMat[row + col*n];
  } 
};



template<typename real>
int mainFun(int argc,char* argv[]) {
  //  using namespace mat;

  try {
    typedef mat::Matrix<real, real> Mat_1;
    typedef mat::Matrix<real, Mat_1> Mat_2;
    typedef mat::Matrix<real, Mat_2> Mat_3;
    //    typedef mat::Vector<real, real > Vec_1;
    //    typedef mat::Vector<real, Vec_1> Vec_2;
    //    typedef mat::Vector<real, Vec_2> Vec_3;
  
    typedef Mat_3 matri;
    //    typedef Vec_3 vect;
    typedef mat::MatrixSymmetric<real, matri> symmMatrix;
    //    typedef mat::MatrixGeneral<real, matri> generalMatrix;
    //    typedef mat::VectorGeneral<real, vect> generalVector;
  
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

    /********** Initialization of SizesAndBlocks                            */
    int size = 7; /* Use weird size to find more bugs. */
    int nlevels = 3;
    std::vector<int> blockSizes(nlevels);
    blockSizes[nlevels - 1] = 1;
#if 1
    blockSizes[nlevels - 2] = 1;
    blockSizes[nlevels - 3] = 5;
#else
    for (int ind = nlevels - 2; ind >= 0; ind--)
      blockSizes[ind] = blockSizes[ind + 1] * 10;
#endif
    mat::SizesAndBlocks rows(blockSizes, size);
    mat::SizesAndBlocks cols(blockSizes, size);

    real fullMat[]={5.8009e-01, 2.3252e-01, 8.0704e-02, 5.4284e-03, 1.1523e-02, 1.0574e-03, 1.3265e-03, 2.3252e-01, 6.5275e-01, 1.2204e-01, 6.4472e-02, 2.7865e-02, 2.3571e-03, 1.4737e-03, 8.0704e-02, 1.2204e-01, 1.9027e-01, 1.0961e-01, 1.1393e-01, 2.9060e-02, 4.1039e-04, 5.4284e-03, 6.4472e-02, 1.0961e-01, 7.0819e-01, 4.5117e-02, 8.3039e-02, 4.0969e-02, 1.1523e-02, 2.7865e-02, 1.1393e-01, 4.5117e-02, 1.0802e-01, 1.2922e-01, 8.2088e-02, 1.0574e-03, 2.3571e-03, 2.9060e-02, 8.3039e-02, 1.2922e-01, 9.3466e-01, 3.1171e-01, 1.3265e-03, 1.4737e-03, 4.1039e-04, 4.0969e-02, 8.2088e-02, 3.1171e-01, 5.4863e-01};

    /* Eigs = 0.010467 0.209761 0.375634 0.377040 0.683334 0.905561 1.160811 */
    mat::Interval<real> gap     (0.4, 0.6);
    mat::Interval<real> allEigs (0.0  , 1.2);
    
    symmMatrix* syA = new symmMatrix;
    syA->resetSizesAndBlocks(rows,cols);
    setFromFullRule<real> sr(fullMat, size);
    syA->setElementsByRule(sr);
    std::cout<<"Norm of syA = "<<syA->eucl(1e-7)<<std::endl;
    symmMatrix* syB = new symmMatrix;
    syB->resetSizesAndBlocks(rows,cols);
    syB->random();
    (*syB) *= real(1e-3);
    std::cout<<"Norm of syB = "<<syB->eucl(1e-7)<<std::endl;
    
    std::vector<symmMatrix*> F;
    std::vector<symmMatrix*> D;
    F.push_back(syA);
    F.push_back(syB);
    //    real deltaMax = 0.1;
    //    real errorTol = 1e-8;
    //    real errorTol = 1e-25; // TRY THIS!!
#if 0
    per::Perturbation<real, symmMatrix, generalVector> 
      perturbObject(F, D, gap, allEigs, deltaMax, errorTol, mat::euclNorm, myVector);
    perturbObject.perturb();

    std::vector<real> idemErrors;
    perturbObject.checkIdempotencies(idemErrors);
    for (unsigned int ind = 0; ind < idemErrors.size();++ind)
      std::cout<<"idemErrors["<<ind<<"] = "<<idemErrors[ind]<<std::endl;

    std::vector<real> commErrors;
    generalMatrix dummyMat;
    perturbObject.checkCommutators(commErrors,dummyMat);
    for (unsigned int ind = 0; ind < commErrors.size();++ind)
      std::cout<<"commErrors["<<ind<<"] = "<<commErrors[ind]<<std::endl;
    
    real subsError;
    perturbObject.checkMaxSubspaceError(subsError);
    std::cout<<"subsError = "<<subsError<<"    errorTol = "<<errorTol<<std::endl;

#endif

  } /* end try */
  catch (mat::Failure e) {
    std::cout << "Failure caught: "<<e.what() << std::endl;
    std::exit(1);
  }
  catch (std::exception e) {
    std::cout << "Exception caught: "<<e.what() << std::endl;
    std::exit(1);
  }
  return 0;

} /*end mainFun */

int main(int argc,char* argv[]){
  
  std::cout<<"=========================================================\n"
	   <<"Testing density matrix perturbation with double precision:\n"
	   <<"=========================================================\n";
  if (!mainFun<double>(argc,argv)) {
    std::cout
      <<"Perturbation tests with double precision completed successfully!" 
      <<std::endl<<std::endl;
  }
#if 0
  std::cout<<"=========================================================\n"
	   <<"Testing  density matrix perturbation with single precision:\n"
	   <<"=========================================================\n";
  if (!mainFun<float>(argc,argv))
    std::cout
      <<"Perturbation tests with single precision completed successfully!" 
      <<std::endl<<std::endl;
#endif
  std::exit(0);
}
