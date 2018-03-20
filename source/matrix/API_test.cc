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

/** @file API_test.cc Test of the matrix library
 *
 * Copyright(c) Emanuel Rubensson 2005
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date September 2005
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

template<typename Treal>
inline bool realIsSingle() {return false;}
template<>
inline bool realIsSingle<float>() {return true;}


template<class Treal>
static Treal maxdiff(std::vector<Treal> const & f1,
		     std::vector<Treal> const & f2);
template<class Treal>
static Treal maxdiff_tri(const Treal* f1,const Treal* f2,int size);
template<class Treal>
static Treal frobdiff(const Treal* f1,const Treal* f2,int size);

using namespace mat;

template<typename Treal>
class Sum {
public:
  Treal accumulate(const Treal& a, int row, int col) {
    return a;
  }
};
template<typename Treal>
class expRule {
public:
  Treal set(int const row, int const col) {
    return (rand() / (Treal)RAND_MAX)  * 
      template_blas_exp(-(template_blas_fabs(Treal(row)-Treal(col))));
  } 
};

template<typename Treal, typename Tmatrix>
static void
testAccumulation(const Tmatrix& syFock, int size, Treal *fockfull) {
  Tmatrix f(syFock);
  Treal ref = 0;
  Sum<Treal> sumOperator;
  Treal res = mat::accumulate(f, sumOperator);
  for(int i=0; i<size; i++)
    for(int j=0; j<size; j++)
      ref += fockfull[j+i*size];
  if( template_blas_fabs(ref-res) / template_blas_fabs(ref) > template_blas_sqrt(mat::getRelPrecision<Treal>()) ) {
    std::cout << "Reference: "<< ref << " computed:" << res << 
      "\nAccumulate test failed.\n";
    std::exit(1);
  } else
    std::cout << "\nAccumulate test OK.\n";
}

template<typename Treal>
bool dotIsBroken() {
  Treal x[8];
  Treal y[8];
  x[0] = 1;
  y[0] = 1;
  int n = 1;
  int incx = 1;
  int incy = 1;
  Treal sdot_result = mat::dot(&n, x, &incx, y, &incy);
  /* Commment: return type of sdot is different in different libraries. 
   * In ATLAS, GOTO, and ACML the return type is double
   * In the 'netlib version' and in MKL it is float
   */
  return template_blas_fabs(sdot_result - (Treal)1.0) > 0.001;  
}

template<typename real>
int mainFun(int argc,char* argv[]) {
  if (dotIsBroken<real>()) {
    std::cout<<"    FAIL: BLAS dot product is broken for selected precision! "
	     <<std::endl<<"    Aborting test"<<std::endl;
    return 1;
  }

  /* different precision                      */
    
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
  typedef arn::LanczosLargestMagnitudeEig<real, symmMatrix, normalVector>
    myLanczosType;
  
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
  
  

  try {
    real epsilon = template_blas_sqrt(mat::getRelPrecision<real>());
    //    int tmpexp;
    int info;

  
    /********** Initialization of SizesAndBlocks                            */
    int size = 87; /* Use weird size to find more bugs. */
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
    
    /********** Obtain row and column index vectors. */
    std::vector<int> rowindex(size * size);
    std::vector<int> colindex(size * size);
    for (int col = 0; col < size; col++)
      for (int row = 0; row < size; row++) {
	rowindex[col * size + row] = row;
	colindex[col * size + row] = col;
      }
    
    std::vector<int> rowindexUpperTriangle((size * (size+1)) / 2);
    std::vector<int> colindexUpperTriangle((size * (size+1)) / 2);
    /* Columnwise only upper triangular storage. */ {
      int ind = 0;
      for (int col = 0; col < size; col++)
	for (int row = 0; row <= col; row++) {
	  rowindexUpperTriangle[ind] = row;
	  colindexUpperTriangle[ind] = col;
	  ++ind;
	}
    }

    
    /********** Test of matrices and vectors                      
     *  The tests following here are testing fundational behavior of the 
     *  data structure, such as constructors, assignments, etc.
     */
    
    /**** General matrices - general tests. */ {
      normalMatrix A, B;
      A.resetSizesAndBlocks(rows,cols);
      A.randomZeroStructure(0.7);
      std::cout<<"Matrix has "<<A.nnz()<<" nonzero elements out of "
	       <<size*size<<" possible.\n";
      //      A.random();
      B = A;             /* Test assignment.         */
      normalMatrix C(A); /* Test copy constructor.   */
      /* Test get_values and assignment and copy */
      std::vector<real> valuesA;
      A.get_values(rowindex, colindex, valuesA);
      std::vector<real> valuesB;
      B.get_values(rowindex, colindex, valuesB);
      std::vector<real> valuesC;
      C.get_values(rowindex, colindex, valuesC);
      assert(valuesA.size() == valuesB.size());
      assert(valuesA.size() == valuesC.size());
      for (unsigned int ind = 0; ind < valuesA.size(); ++ind) {
	assert(valuesA[ind] == valuesB[ind]);
	assert(valuesA[ind] == valuesC[ind]);
      }
      
      /* Test get_all_values */
      std::vector<int> rowindexA2;
      std::vector<int> colindexA2;
      std::vector<real> valuesA2;
      A.get_all_values(rowindexA2, colindexA2, valuesA2);
      int row = 0;
      int col = 0;
      for (unsigned int ind = 0; ind < valuesA2.size(); ++ind) {
	row = rowindexA2[ind];
	col = colindexA2[ind];
	assert(valuesA2[ind] == valuesA[col * size + row]);
      }
      
      /* Test assign_from_sparse */
      normalMatrix D;
      D.assign_from_sparse(rowindex, colindex, valuesA, rows, cols);
      std::vector<real> valuesD;
      D.get_values(rowindex, colindex, valuesD);
      assert(valuesA.size() == valuesD.size());
      for (unsigned int ind = 0; ind < valuesA.size(); ++ind) 
	assert(valuesA[ind] == valuesD[ind]);
    } /* End of general matrices test. */


    /**** General matrix - Test assignments using permutation. */ {
      /* Get a random matrix */
      normalMatrix A;
      A.resetSizesAndBlocks(rows,cols);
      A.randomZeroStructure(0.7);
      //      A.random();
      /* Get a random permutation */
      std::vector<int> permutationVec(size);
      for (int ind = 0; ind < size; ++ind )
	permutationVec[ind] = ind;
      std::random_shuffle ( permutationVec.begin(), permutationVec.end() );
      /* Get inverse permutation */
      std::vector<int> inversePermutationVec(size);
      for (int ind = 0; ind < size; ++ind) {
	inversePermutationVec[permutationVec[ind]] = ind;
      }
      
      /* Test get values with permutation */
      std::vector<real> fullMatRef(size * size);
      {
	std::vector<real> valuesARef;
	std::vector<int> permRowIndex(rowindex.size());
	std::vector<int> permColIndex(colindex.size());
	for (unsigned int ind = 0; ind < rowindex.size(); ++ind) {
	  permRowIndex[ind] = permutationVec[rowindex[ind]];
	  permColIndex[ind] = permutationVec[colindex[ind]];
	}
	A.get_values(permRowIndex, permColIndex, valuesARef);
	std::vector<real> valuesAPerm;
	A.get_values(rowindex, colindex, valuesAPerm, 
		     permutationVec, permutationVec);
	assert(valuesARef.size() == valuesAPerm.size());
	for (unsigned int ind = 0; ind < valuesARef.size(); ++ind) {
	  assert(permRowIndex[ind] == permutationVec[rowindex[ind]]);
	  assert(permColIndex[ind] == permutationVec[colindex[ind]]);
	  assert(valuesARef[ind] == valuesAPerm[ind]);
	  fullMatRef[rowindex[ind] + A.get_nrows() * colindex[ind]] =
	    valuesARef[ind];
	}
      }

      /* Test get all values with permutation */ {
	std::vector<int> rowIndex;
	std::vector<int> colIndex;
	std::vector<real>   values;
	A.get_all_values(rowIndex, colIndex, values, 
			 inversePermutationVec, inversePermutationVec);
	for (unsigned int ind = 0; ind < values.size(); ++ind) {
	  assert(values[ind] == 
		 fullMatRef[rowIndex[ind] + A.get_nrows() * colIndex[ind]] );
	}
      }

      /* Test fullMatrix */ {
	std::vector<real> fullMat;
	A.fullMatrix(fullMat,
		     inversePermutationVec, inversePermutationVec);
	for (unsigned int ind = 0; ind < fullMat.size(); ++ind)
	  assert(fullMat[ind] == fullMatRef[ind]);
      }
      
    } /* End: General matrix - Test assignments using permutation. */

    
    
#if 1
    /**** Symmetric matrices - general tests. */ {
      symmMatrix syA, syB;
      syA.resetSizesAndBlocks(rows,cols);
      syA.randomZeroStructure(0.7);
      //      syA.random();
      syA.add_identity(100);
      syB = syA;
      symmMatrix syC(syA); /* Test copy */
      
      /* Test get_values and assignment and copy */
      std::vector<real> valuesA;
      syA.get_values(rowindexUpperTriangle, colindexUpperTriangle, valuesA);
      std::vector<real> valuesB;
      syB.get_values(rowindexUpperTriangle, colindexUpperTriangle, valuesB);
      std::vector<real> valuesC;
      syC.get_values(rowindexUpperTriangle, colindexUpperTriangle, valuesC);
      assert(valuesA.size() == valuesB.size());
      assert(valuesA.size() == valuesC.size());
      for (unsigned int ind = 0; ind < valuesA.size(); ++ind) {
	assert(valuesA[ind] == valuesB[ind]);
	assert(valuesA[ind] == valuesC[ind]);
      }
      
      /* Test get_all_values */
      std::vector<int> rowindexA2;
      std::vector<int> colindexA2;
      std::vector<real> valuesA2;
      syA.get_all_values(rowindexA2, colindexA2, valuesA2);
      int row = 0;
      int col = 0;
      int offset = 0;
      for (unsigned int ind = 0; ind < valuesA2.size(); ++ind) {
	row = rowindexA2[ind];
	col = colindexA2[ind];
	offset = 0;
	for (int indOff = 1; indOff <= col; ++indOff)
	  offset = offset + indOff;
	assert(valuesA2[ind] == valuesA[offset + row]);
      }

      /* Test assign_from_sparse */
      symmMatrix syD;
      syD.assign_from_sparse(rowindexUpperTriangle, 
			     colindexUpperTriangle, 
			     valuesA, 
			     rows, cols);
      std::vector<real> valuesD;
      syD.get_values(rowindexUpperTriangle, 
		     colindexUpperTriangle, 
		     valuesD);
      assert(valuesA.size() == valuesD.size());
      for (unsigned int ind = 0; ind < valuesA.size(); ++ind) 
	assert(valuesA[ind] == valuesD[ind]);
      
      
    } /* End of symmetric matrices test. */

    /**** Regression test for assign */ {
      symmMatrix A,B;
      A.resetSizesAndBlocks(rows,cols);
      A.random();
      B = real(1.3) * A;
    }

    /**** Symmetric matrices - Test assignments using permutation. */ {
      /* Get a random matrix */
      symmMatrix A;
      A.resetSizesAndBlocks(rows,cols);
      A.randomZeroStructure(0.7);
      //      A.random();
      /* Get a random permutation */
      std::vector<int> permutationVec(size);
      std::vector<int> noPermutationVec(size);
      for (int ind = 0; ind < size; ++ind )
	permutationVec[ind] = ind;
      noPermutationVec = permutationVec;
      std::random_shuffle ( permutationVec.begin(), permutationVec.end() );
      /* Get inverse permutation */
      std::vector<int> inversePermutationVec(size);
      for (int ind = 0; ind < size; ++ind) {
	inversePermutationVec[permutationVec[ind]] = ind;
      }
      
      /* Test get values with permutation */
      std::vector<real> fullMatRef(size * size);
      {
	std::vector<real> valuesARef;
	std::vector<int> permRowIndex(rowindexUpperTriangle.size());
	std::vector<int> permColIndex(colindexUpperTriangle.size());
	for (unsigned int ind = 0; ind < rowindexUpperTriangle.size(); ++ind) {
	  permRowIndex[ind] = permutationVec[rowindexUpperTriangle[ind]];
	  permColIndex[ind] = permutationVec[colindexUpperTriangle[ind]];
	}
	A.get_values(permRowIndex, permColIndex, valuesARef,
		     noPermutationVec, noPermutationVec);
	std::vector<real> valuesAPerm;
	A.get_values(rowindexUpperTriangle, colindexUpperTriangle, valuesAPerm, 
		     permutationVec, permutationVec);
	assert(valuesARef.size() == valuesAPerm.size());
	for (unsigned int ind = 0; ind < valuesARef.size(); ++ind) {
	  assert(permRowIndex[ind] == 
		 permutationVec[rowindexUpperTriangle[ind]]);
	  assert(permColIndex[ind] == 
		 permutationVec[colindexUpperTriangle[ind]]);
	  assert(valuesARef[ind] == valuesAPerm[ind]);
	  fullMatRef[rowindexUpperTriangle[ind] + 
		     A.get_nrows() * colindexUpperTriangle[ind]] =
	    valuesARef[ind];
	  fullMatRef[colindexUpperTriangle[ind] + 
		     A.get_nrows() * rowindexUpperTriangle[ind]] =
	    valuesARef[ind];

	}
      }
      std::vector<int> rowIndex;
      std::vector<int> colIndex;
      std::vector<real>   values;
      /* Test get all values with permutation */ {
	A.get_all_values(rowIndex, colIndex, values, 
			 inversePermutationVec, inversePermutationVec);
	for (unsigned int ind = 0; ind < values.size(); ++ind) {
	  assert(values[ind] == 
		 fullMatRef[rowIndex[ind] + A.get_nrows() * colIndex[ind]]);
	}
	
	/* Test assign from sparse with permutation*/ {
	  symmMatrix syB;
	  syB.resetSizesAndBlocks(rows,cols);
	  syB.assign_from_sparse(rowIndex, colIndex, values, 
				 permutationVec, permutationVec);
	  std::vector<real> valuesCopy;
	  syB.get_values(rowIndex, colIndex, valuesCopy, 
			 permutationVec, permutationVec);
	  for (unsigned int ind = 0; ind < values.size(); ++ind) {
	    assert(values[ind] == valuesCopy[ind]);
	  }
	  
	}
	
      }

      /* Test fullMatrix */ {
	std::vector<real> fullMat;
	A.fullMatrix(fullMat,
		     inversePermutationVec, inversePermutationVec);
	for (unsigned int ind = 0; ind < fullMat.size(); ++ind)
	  assert(fullMat[ind] == fullMatRef[ind]);
      }
      
      
    } /* End: Symmetric matrices - Test assignments using permutation. */
    
    
    
    /**** Triangular matrices. */ {
      triangMatrix trA, trB;
      trA.resetSizesAndBlocks(rows,cols);
      trA.random();
      trA.add_identity(100);
      trB = trA;
      triangMatrix trC(trA); /* Test copy */
      
      /* Test get_values and assignment and copy */
      std::vector<real> valuesA;
      trA.get_values(rowindexUpperTriangle, colindexUpperTriangle, valuesA);
      std::vector<real> valuesB;
      trB.get_values(rowindexUpperTriangle, colindexUpperTriangle, valuesB);
      std::vector<real> valuesC;
      trC.get_values(rowindexUpperTriangle, colindexUpperTriangle, valuesC);
      assert(valuesA.size() == valuesB.size());
      assert(valuesA.size() == valuesC.size());
      for (unsigned int ind = 0; ind < valuesA.size(); ++ind) {
	assert(valuesA[ind] == valuesB[ind]);
	assert(valuesA[ind] == valuesC[ind]);
      }
      
      /* Test get_all_values */
      std::vector<int> rowindexA2;
      std::vector<int> colindexA2;
      std::vector<real> valuesA2;
      trA.get_all_values(rowindexA2, colindexA2, valuesA2);
      int row = 0;
      int col = 0;
      int offset = 0;
      for (unsigned int ind = 0; ind < valuesA2.size(); ++ind) {
	row = rowindexA2[ind];
	col = colindexA2[ind];
	offset = 0;
	for (int indOff = 1; indOff <= col; ++indOff)
	  offset = offset + indOff;
	assert(valuesA2[ind] == valuesA[offset + row]);
      }
      
      /* Test assign_from_sparse */
      triangMatrix trD;
      trD.assign_from_sparse(rowindexUpperTriangle, 
			     colindexUpperTriangle, 
			     valuesA, 
			     rows, cols);
      std::vector<real> valuesD;
      trD.get_values(rowindexUpperTriangle, 
		     colindexUpperTriangle, 
		     valuesD);
      assert(valuesA.size() == valuesD.size());
      for (unsigned int ind = 0; ind < valuesA.size(); ++ind) 
	assert(valuesA[ind] == valuesD[ind]);
      
    } /* End of triangular matrices test. */

    /******** Test assign with permutation (again) */ {
      normalMatrix A, B;
      real myReals[] = {3.4, 65.246, 234.23456, 655.1235};
      std::vector<real> values (myReals, myReals+4);
      int myRowInds[] = {size/4, size/3, size/2, size-1};
      std::vector<int> rowInd (myRowInds, myRowInds+4);
      int myColInds[] = {size/5, size/3, size/2, size-1};
      std::vector<int> colInd (myColInds, myColInds+4);
      std::vector<int> perm(size);
      for (unsigned int ind = 0; ind < perm.size(); ++ind) 
	perm[ind] = (int)ind;
      std::random_shuffle ( perm.begin(), perm.end() );
      A.assign_from_sparse(rowInd, colInd, values, rows, cols, perm, perm);
      
      std::vector<int> invPerm(size);
      for (unsigned int ind = 0; ind < perm.size(); ++ind) 
	invPerm[perm[ind]] = (int)ind;

      std::vector<int> rowIndConfirm;
      std::vector<int> colIndConfirm;
      std::vector<real> valuesConfirm;
      A.get_all_values(rowIndConfirm, colIndConfirm, valuesConfirm,
		       invPerm,invPerm);
      
#if 0
      // print out content:
      std::cout << "row, col, val:\n";
      for (unsigned int ind = 0;ind < values.size(); ++ind) 
	std::cout << rowInd[ind] << ", " 
		  << colInd[ind] << ", " 
		  << values[ind] << std::endl;

      std::cout << "row, col, val (conf):\n";
      for (unsigned int ind = 0;ind < valuesConfirm.size(); ++ind) 
	std::cout << rowIndConfirm[ind] << ", " 
		  << colIndConfirm[ind] << ", " 
		  << valuesConfirm[ind] << std::endl;

      //      A.resetSizesAndBlocks(rows,cols);
#endif 

    }

    /****** Test assign from full, general matrix */ {
      std::vector<real> fullMat(size*size);
      for (unsigned int ind = 0; ind < fullMat.size(); ++ind)
      	fullMat[ind] = (rand() / (real)RAND_MAX);
      normalMatrix A;
      A.resetSizesAndBlocks(rows,cols);
      A.assignFromFull(fullMat);
      std::vector<real> fullMatCopy(size*size);
      A.fullMatrix(fullMatCopy);
      assert(maxdiff(fullMat, fullMatCopy) < epsilon);
    } /* End: Test assign from full */

    /****** Test assign from full, symmetric matrix */ {
      std::vector<real> fullMat(size*size);
      real randomNumber;
      for (int col = 0; col < size; ++col)
	for (int row = 0; row < size; ++row) {
	  randomNumber = (rand() / (real)RAND_MAX);
	  fullMat[row + col * size] = randomNumber;
	  fullMat[col + row * size] = randomNumber;
	}
      symmMatrix syA;
      syA.resetSizesAndBlocks(rows,cols);
      syA.assignFromFull(fullMat);
      std::vector<real> fullMatCopy(size*size);
      syA.fullMatrix(fullMatCopy);
      assert(maxdiff(fullMat, fullMatCopy) < epsilon);
      
    } /* End: Test assign from full */

    /****** Test assign from full using permutation */{
      /* Get a random permutation */
      std::vector<int> permutationVec(size);
      for (int ind = 0; ind < size; ++ind )
	permutationVec[ind] = ind;
      std::random_shuffle ( permutationVec.begin(), permutationVec.end() );
      /* Get inverse permutation */
      std::vector<int> inversePermutationVec(size);
      for (int ind = 0; ind < size; ++ind) {
	inversePermutationVec[permutationVec[ind]] = ind;
      }
      
      std::vector<real> fullMat(size*size);
      real randomNumber;
      for (int col = 0; col < size; ++col)
	for (int row = 0; row < size; ++row) {
	  randomNumber = (rand() / (real)RAND_MAX);
	  fullMat[row + col * size] = randomNumber;
	  fullMat[col + row * size] = randomNumber;
	}

      symmMatrix syA;
      syA.resetSizesAndBlocks(rows,cols);
      syA.assignFromFull(fullMat,permutationVec, permutationVec);
      std::vector<real> fullMatCopy(size*size);
      syA.fullMatrix(fullMatCopy, inversePermutationVec, inversePermutationVec);
      
      std::cout<<" Max Absolute Error: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<maxdiff(fullMat, fullMatCopy);
      
      assert(maxdiff(fullMat, fullMatCopy) < epsilon);
      
    } /* End: test assign from full using permutation */
    
    


    /* Test copy and assignment between symmetric and general matrix types 
     * This tests the nosymToSym and symToNosym functions as well as 
     * transpose
     *
     */ {
      symmMatrix syA, syD;
      normalMatrix B;
      syA.resetSizesAndBlocks(rows,cols);
      syA.randomZeroStructure(0.7);
      //      syA.random();
      syA.add_identity(100);
      B = syA;
      normalMatrix C(syA); /* Test copy */
      syD = B;
      symmMatrix syE(B);
      
      std::vector<real> valuesA;
      syA.get_values(rowindexUpperTriangle, colindexUpperTriangle, valuesA);
      std::vector<real> valuesB;
      B.get_values(rowindexUpperTriangle, colindexUpperTriangle, valuesB);
      std::vector<real> valuesC;
      C.get_values(rowindexUpperTriangle, colindexUpperTriangle, valuesC);
      std::vector<real> valuesD;
      syD.get_values(rowindexUpperTriangle, colindexUpperTriangle, valuesD);
      std::vector<real> valuesE;
      syE.get_values(rowindexUpperTriangle, colindexUpperTriangle, valuesE);
      assert(valuesA.size() == valuesB.size());
      assert(valuesA.size() == valuesC.size());
      assert(valuesA.size() == valuesD.size());
      assert(valuesA.size() == valuesE.size());
      for (unsigned int ind = 0; ind < valuesA.size(); ++ind) {
	assert(valuesA[ind] == valuesB[ind]);
	assert(valuesA[ind] == valuesC[ind]);
	assert(valuesA[ind] == valuesD[ind]);
	assert(valuesA[ind] == valuesE[ind]);
      }
    } /* End of copy between symmetric and general matrix types. */

    /* Test copy and assignment from triangular to general matrix type 
     * Nothing should really go wrong here.
     */ {
      triangMatrix trA;
      trA.resetSizesAndBlocks(rows,cols);
      trA.random();
      normalMatrix B(trA);
      normalMatrix C;
      C = trA;
      std::vector<real> valuesA;
      trA.get_values(rowindexUpperTriangle, colindexUpperTriangle, valuesA);
      std::vector<real> valuesB;
      B.get_values(rowindexUpperTriangle, colindexUpperTriangle, valuesB);
      std::vector<real> valuesC;
      C.get_values(rowindexUpperTriangle, colindexUpperTriangle, valuesC);
      assert(valuesA.size() == valuesB.size());
      assert(valuesA.size() == valuesC.size());
      for (unsigned int ind = 0; ind < valuesA.size(); ++ind) {
	assert(valuesA[ind] == valuesB[ind]);
	assert(valuesA[ind] == valuesC[ind]);
      }
    } /* End of test of copy from triangular to general matrix type */
    
    
    /* *********  First gemm test: C = alpha * A * B     */ {
      
      const real zero = 0.0;
      const real alpha = 2.346;
      normalMatrix A, B, C;
      A.resetSizesAndBlocks(rows,cols);
      A.randomZeroStructure(0.7);
      //A.random();
      B.resetSizesAndBlocks(rows,cols);
      B.randomZeroStructure(0.7);
      //B.random();
      std::vector<real> valuesA;
      A.get_values(rowindex, colindex, valuesA);
      std::vector<real> valuesB;
      B.get_values(rowindex, colindex, valuesB);
      std::vector<real> valuesCblas(valuesA.size(),0);
      std::cout<<"\nMatrix multiply test:  \n";
      gemm("N", "N", &size, &size, &size, &alpha, &valuesA[0], &size, 
	   &valuesB[0], &size, &zero, &valuesCblas[0], &size);
      C = alpha * A * B;
      std::vector<real> valuesC;
      C.get_values(rowindex, colindex, valuesC);
      
      std::cout<<" Max Absolute Error: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<maxdiff(valuesC,valuesCblas);
      if (maxdiff(valuesC,valuesCblas) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      }
    } /* End of first gemm test */

    /* *********  Second gemm test: C = alpha * A' * B + beta * C     */ {
      
      const real alpha = 2.346;
      const real beta  = 3.348;      
      normalMatrix A, B, C;
      A.resetSizesAndBlocks(rows,cols);
      A.randomZeroStructure(0.7);
      //      A.random();
      B.resetSizesAndBlocks(rows,cols);
      B.randomZeroStructure(0.7);
      //      B.random();
      C.resetSizesAndBlocks(rows,cols);
      C.randomZeroStructure(0.7);
      //      C.random();
      std::vector<real> valuesA;
      A.get_values(rowindex, colindex, valuesA);
      std::vector<real> valuesB;
      B.get_values(rowindex, colindex, valuesB);
      std::vector<real> valuesCblas;
      C.get_values(rowindex, colindex, valuesCblas);

      std::cout<<"\nTransposed matrix multiply test:  \n";
      gemm("T", "N", &size, &size, &size, &alpha, &valuesA[0], &size,
	   &valuesB[0], &size, &beta, &valuesCblas[0], &size);
      C = alpha * transpose(A) * B + beta * C;
      std::vector<real> valuesC;
      C.get_values(rowindex, colindex, valuesC);

      std::cout<<" Max Absolute Error: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<maxdiff(valuesC,valuesCblas);
      if (maxdiff(valuesC,valuesCblas) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      }
    }  /* End of second gemm test */



    /* *********  Symm test: C = alpha * syA * B + beta * C     */ {

      const real alpha = 2.346;
      const real beta  = 3.348;      
      symmMatrix syA;
      normalMatrix B, C;
      syA.resetSizesAndBlocks(rows,cols);
      syA.randomZeroStructure(0.7);
      //      syA.random();
      B.resetSizesAndBlocks(rows,cols);
      B.randomZeroStructure(0.7);
      //      B.random();
      C.resetSizesAndBlocks(rows,cols);
      C.randomZeroStructure(0.7);
      //      C.random();
      std::vector<real> valuesAPacked;
      syA.get_values(rowindexUpperTriangle, colindexUpperTriangle, 
		     valuesAPacked);
      std::vector<real> valuesB;
      B.get_values(rowindex, colindex, valuesB);
      std::vector<real> valuesCblas;
      C.get_values(rowindex, colindex, valuesCblas);
      std::vector<real> valuesA(valuesB.size(),0);
      packedtofull(&valuesAPacked[0], &valuesA[0], size);
      
      std::cout<<"\nSymmetric matrix mul (symm) test: \n";
      symm("L", "U", &size, &size, &alpha, &valuesA[0], &size, &valuesB[0], 
	   &size, &beta, &valuesCblas[0], &size);
      C = alpha * syA * B + beta * C;
      std::vector<real> valuesC;
      C.get_values(rowindex, colindex, valuesC);
      
      std::cout<<" Max Absolute Error: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<maxdiff(valuesC,valuesCblas);
      if (maxdiff(valuesC,valuesCblas) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      }
      
    }  /* End of symm test */
    
    
    /* *********  Syrk test: syB = alpha * A * A' + beta * syB     */ {

      const real alpha = 2.346;
      const real beta  = 3.348;      
      normalMatrix A;
      symmMatrix syB;
      A.resetSizesAndBlocks(rows,cols);
      A.randomZeroStructure(0.7);
      //      A.random();
      syB.resetSizesAndBlocks(rows,cols);
      syB.randomZeroStructure(0.7);
      //      syB.random();
      std::vector<real> valuesA;
      A.get_values(rowindex, colindex, valuesA);
      std::vector<real> valuesBPacked;
      syB.get_values(rowindexUpperTriangle, colindexUpperTriangle, 
		     valuesBPacked);
      std::vector<real> valuesBblas(valuesA.size(),0);
      packedtofull(&valuesBPacked[0], &valuesBblas[0], size);
      
      std::cout<<"\nSymmetric rank-k update (syrk) test: \n";
      syrk("U", "N", &size, &size, &alpha, &valuesA[0], &size,
	   &beta, &valuesBblas[0], &size);
      syB = alpha * A * transpose(A) + beta * syB;
      syB.get_values(rowindexUpperTriangle, colindexUpperTriangle, 
		     valuesBPacked);
      std::vector<real> valuesB(valuesA.size(),0);
      packedtofull(&valuesBPacked[0], &valuesB[0], size);
      
      std::cout<<" Max Absolute Error: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<maxdiff_tri(&valuesB[0], &valuesBblas[0], size);
      if (maxdiff_tri(&valuesB[0], &valuesBblas[0], size) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      }
      
    }  /* End of syrk test */
    

    /* *********  Sysq test: syB = alpha * syA * syA + beta * syB     */ {
      const real alpha = 2.346;
      const real beta  = 3.348;      
      symmMatrix syA, syB;
      syA.resetSizesAndBlocks(rows,cols);
      syA.randomZeroStructure(0.7);
      //      syA.random();
      syB.resetSizesAndBlocks(rows,cols);
      syB.randomZeroStructure(0.7);
      //      syB.random();
      std::vector<real> valuesAPacked;
      syA.get_values(rowindexUpperTriangle, colindexUpperTriangle, 
		     valuesAPacked);
      std::vector<real> valuesBPacked;
      syB.get_values(rowindexUpperTriangle, colindexUpperTriangle, 
		     valuesBPacked);
      std::vector<real> valuesA(size * size,0);
      packedtofull(&valuesAPacked[0], &valuesA[0], size);
      std::vector<real> valuesBblas(size * size,0);
      packedtofull(&valuesBPacked[0], &valuesBblas[0], size);
      
      std::cout<<"\nSymmetric matrix square test: \n";
      gemm("N", "N", &size, &size, &size, &alpha, &valuesA[0], &size, 
	   &valuesA[0], &size, &beta, &valuesBblas[0], &size);
      syB = alpha * syA * syA + beta * syB;
      syB.get_values(rowindexUpperTriangle, colindexUpperTriangle, 
		     valuesBPacked);
      std::vector<real> valuesB(size * size, 0);
      packedtofull(&valuesBPacked[0], &valuesB[0], size);
      
      std::cout<<" Max Absolute Error: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<maxdiff(valuesB, valuesBblas);
      if (maxdiff(valuesB, valuesBblas) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      } 
    }  /* End of sysq test */
    
    /* *********  Ssmm test: C = alpha * syA * syB + beta * C     */ {
      const real alpha = 2.346;
      const real beta  = 3.348;      
      symmMatrix syA, syB;
      normalMatrix C;
      syA.resetSizesAndBlocks(rows,cols);
      syA.randomZeroStructure(0.7);
      //      syA.random();
      syB.resetSizesAndBlocks(rows,cols);
      syB.randomZeroStructure(0.7);
      //      syB.random();
      C.resetSizesAndBlocks(rows,cols);
      C.randomZeroStructure(0.7);
      //      C.random();
      std::vector<real> valuesAPacked;
      syA.get_values(rowindexUpperTriangle, colindexUpperTriangle, 
		     valuesAPacked);
      std::vector<real> valuesA(size * size,0);
      packedtofull(&valuesAPacked[0], &valuesA[0], size);
      std::vector<real> valuesBPacked;
      syB.get_values(rowindexUpperTriangle, colindexUpperTriangle, 
		     valuesBPacked);
      std::vector<real> valuesB(size * size,0);
      packedtofull(&valuesBPacked[0], &valuesB[0], size);
      std::vector<real> valuesCblas;
      C.get_values(rowindex, colindex, valuesCblas);
      
      std::cout<<"\nSymmetric-symmetric matrix multiply test: \n";
      gemm("N", "N", &size, &size, &size, &alpha, &valuesA[0], &size, 
	   &valuesB[0], &size, &beta, &valuesCblas[0], &size);
      C = alpha * syA * syB + beta * C;
      std::vector<real> valuesC;
      C.get_values(rowindex, colindex, valuesC);
      std::cout<<" Max Absolute Error: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<maxdiff(valuesC, valuesCblas);
      if (maxdiff(valuesC, valuesCblas) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      } 
    }  /* End of ssmm test */
    
    /* *********  Trmm test 1: C = alpha * transpose(trA) * C     */ {
      const real alpha = 2.346;
      triangMatrix trA;
      normalMatrix C;
      trA.resetSizesAndBlocks(rows, cols);
      trA.random();
      C.resetSizesAndBlocks(rows, cols);
      C.randomZeroStructure(0.7);
      //      C.random();
      std::vector<real> valuesAPacked;
      trA.get_values(rowindexUpperTriangle, colindexUpperTriangle, 
		     valuesAPacked);
      std::vector<real> valuesA(size * size,0);
      tripackedtofull(&valuesAPacked[0], &valuesA[0], size);
      std::vector<real> valuesCblas;
      C.get_values(rowindex, colindex, valuesCblas);
      
      std::cout<<"\nTriangular-general matrix multiply test 1: \n";
      trmm("L", "U", "T", "N", &size, &size, &alpha, &valuesA[0],
	   &size, &valuesCblas[0], &size);
      C = alpha * transpose(trA) * C;
      std::vector<real> valuesC;
      C.get_values(rowindex, colindex, valuesC);
      std::cout<<" Max Absolute Error: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<maxdiff(valuesC, valuesCblas);
      if (maxdiff(valuesC, valuesCblas) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      } 
    }  /* End of trmm test 1 */
    
    /* *********  Trmm test 2: C = alpha * C * trA     */ {
      const real alpha = 2.346;
      triangMatrix trA;
      normalMatrix C;
      trA.resetSizesAndBlocks(rows, cols);
      trA.random();
      C.resetSizesAndBlocks(rows, cols);
      C.randomZeroStructure(0.7);
      //      C.random();
      std::vector<real> valuesAPacked;
      trA.get_values(rowindexUpperTriangle, colindexUpperTriangle, 
		     valuesAPacked);
      std::vector<real> valuesA(size * size,0);
      tripackedtofull(&valuesAPacked[0], &valuesA[0], size);
      std::vector<real> valuesCblas;
      C.get_values(rowindex, colindex, valuesCblas);
      
      std::cout<<"\nTriangular-general matrix multiply test 2: \n";
      trmm("R", "U", "N", "N", &size, &size, &alpha, &valuesA[0],
	   &size, &valuesCblas[0], &size);
      C = alpha * C * trA;
      std::vector<real> valuesC;
      C.get_values(rowindex, colindex, valuesC);
      std::cout<<" Max Absolute Error: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<maxdiff(valuesC, valuesCblas);
      if (maxdiff(valuesC, valuesCblas) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      } 
    }  /* End of trmm test 2 */
    
    /* *********  Frobenius norm test: general matrix     */ {
      normalMatrix A;
      A.resetSizesAndBlocks(rows,cols);
      A.randomZeroStructure(0.7);
      //      A.random();
      std::vector<real> valuesA;
      A.get_values(rowindex, colindex, valuesA);

      std::cout<<"\nFrobenius norm (general matrix) test: \n";
      real frobNormRefValue = 0;
      for (unsigned int ind = 0; ind < valuesA.size(); ++ind) {
	frobNormRefValue += valuesA[ind] * valuesA[ind];
      }
      frobNormRefValue = template_blas_sqrt(frobNormRefValue);
      real frobNorm = A.frob();
      std::cout<<" Absolute Difference: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<template_blas_fabs(frobNormRefValue - frobNorm);
      if (template_blas_fabs(frobNormRefValue - frobNorm) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      } 
    } /* End of Frobenius norm test general matrix */
      
    /* *********  Frobenius norm test: symmetric matrix     */ {
      symmMatrix syA;
      syA.resetSizesAndBlocks(rows,cols);
      syA.randomZeroStructure(0.7);
      //      syA.random();
      std::vector<real> valuesAPacked;
      syA.get_values(rowindexUpperTriangle, colindexUpperTriangle, 
		     valuesAPacked);
      std::vector<real> valuesA(size * size,0);
      packedtofull(&valuesAPacked[0], &valuesA[0], size);

      std::cout<<"\nFrobenius norm (symmetric matrix) test: \n";
      real frobNormRefValue = 0;
      for (unsigned int ind = 0; ind < valuesA.size(); ++ind) {
	frobNormRefValue += valuesA[ind] * valuesA[ind];
      }
      frobNormRefValue = template_blas_sqrt(frobNormRefValue);
      real frobNorm = syA.frob();
      std::cout<<" Absolute Difference: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<template_blas_fabs(frobNormRefValue - frobNorm);
      if (template_blas_fabs(frobNormRefValue - frobNorm) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      } 
    } /* End of Frobenius norm test symmetric matrix */

    /* *********  Frobenius norm test: triangular matrix     */ {
      triangMatrix trA;
      trA.resetSizesAndBlocks(rows,cols);
      trA.random();
      std::vector<real> valuesAPacked;
      trA.get_values(rowindexUpperTriangle, colindexUpperTriangle, 
		     valuesAPacked);
      std::vector<real> valuesA(size * size,0);
      tripackedtofull(&valuesAPacked[0], &valuesA[0], size);

      std::cout<<"\nFrobenius norm (triangular matrix) test: \n";
      real frobNormRefValue = 0;
      for (unsigned int ind = 0; ind < valuesA.size(); ++ind) {
	frobNormRefValue += valuesA[ind] * valuesA[ind];
      }
      frobNormRefValue = template_blas_sqrt(frobNormRefValue);
      real frobNorm = trA.frob();
      std::cout<<" Absolute Difference: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<template_blas_fabs(frobNormRefValue - frobNorm);
      if (template_blas_fabs(frobNormRefValue - frobNorm) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      } 
    } /* End of Frobenius norm test triangular matrix */
    
    /* *********  Frobenius norm diff test: general matrix     */ {
      normalMatrix A, B;
      A.resetSizesAndBlocks(rows,cols);
      A.randomZeroStructure(0.7);
      //      A.random();
      B = A;
      B.add_identity(1);
      std::cout<<"\nFrobenius norm diff (general matrix) test: \n";
      real frobDiff = normalMatrix::frob_diff(A, B);
      real frobDiffRefValue = template_blas_sqrt((real)size);
      std::cout<<" Absolute Error: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<template_blas_fabs(frobDiffRefValue - frobDiff);
      if (template_blas_fabs(frobDiffRefValue - frobDiff) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      } 
    } /* End of Frobenius norm diff test general matrix */
    
    /* *********  Frobenius norm diff test: symmetric matrix     */ {
      symmMatrix A, B;
      A.resetSizesAndBlocks(rows,cols);
      A.randomZeroStructure(0.7);
      //      A.random();
      B = A;
      B.add_identity(1);
      std::cout<<"\nFrobenius norm diff (symmetric matrix) test: \n";
      real frobDiff = symmMatrix::frob_diff(A, B);
      real frobDiffRefValue = template_blas_sqrt((real)size);
      std::cout<<" Absolute Error: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<template_blas_fabs(frobDiffRefValue - frobDiff);
      if (template_blas_fabs(frobDiffRefValue - frobDiff) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      } 
    } /* End of Frobenius norm diff test symmetric matrix */
    
    /* *********  Write to file test: general matrix     */ {
      normalMatrix A, B;
      A.resetSizesAndBlocks(rows,cols);
      A.randomZeroStructure(0.7);
      //      A.random();
      B = A;
      /* Test write to file */
      A.writeToFile();
      A.readFromFile();
      std::vector<real> valuesA;
      A.get_values(rowindex, colindex, valuesA);
      std::vector<real> valuesB;
      B.get_values(rowindex, colindex, valuesB);
      std::cout<<"\nWrite to file (general matrix) test: \n";
      std::cout<<" Absolute Error: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<maxdiff(valuesA, valuesB);
      if (maxdiff(valuesA, valuesB) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      } 
    } /* End of write to file test general matrix */
    
    /* *********  Trace test: general matrix     */ {
      normalMatrix A;
      A.resetSizesAndBlocks(rows,cols);
      A.randomZeroStructure(0.7);
      //      A.random();
      std::vector<real> valuesA;
      A.get_values(rowindex, colindex, valuesA);

      std::cout<<"\nTrace (general matrix) test: \n";
      real traceRefValue = 0;
      for (int rc = 0; rc < size; ++rc) {
	traceRefValue += valuesA[rc * size + rc];
      }
      real traceValue = A.trace();
      std::cout<<" Absolute Difference: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<template_blas_fabs(traceRefValue - traceValue);
      if (template_blas_fabs(traceRefValue - traceValue) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      } 
    } /* End of trace test general matrix */
    
    /* *********  Addition test: general matrix     */ {
      normalMatrix A, B, C;
      A.resetSizesAndBlocks(rows,cols);
      A.randomZeroStructure(0.7);
      //      A.random();
      std::vector<real> valuesA;
      A.get_values(rowindex, colindex, valuesA);
      B.resetSizesAndBlocks(rows,cols);
      B.randomZeroStructure(0.7);
      //      B.random();
      std::vector<real> valuesB;
      B.get_values(rowindex, colindex, valuesB);
      std::cout<<"\nAddition (general matrix) test: \n";
      C = A + B;
      std::vector<real> valuesC;
      C.get_values(rowindex, colindex, valuesC);
      std::vector<real> valuesCRef(valuesA.size(),0);
      for (unsigned int ind = 0; ind < valuesA.size(); ++ind)
	valuesCRef[ind] = valuesA[ind] + valuesB[ind];
      std::cout<<" Absolute Difference: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<maxdiff(valuesCRef, valuesC);
      if (maxdiff(valuesCRef, valuesC) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      } 
    } /* End of addition test general matrix */

    /* *********  Addition test: symmetric matrix     */ {
      symmMatrix A, B, C;
      A.resetSizesAndBlocks(rows,cols);
      A.randomZeroStructure(0.7);
      //      A.random();
      std::vector<real> valuesA;
      A.get_values(rowindexUpperTriangle, colindexUpperTriangle, valuesA);
      B.resetSizesAndBlocks(rows,cols);
      B.randomZeroStructure(0.7);
      //      B.random();
      std::vector<real> valuesB;
      B.get_values(rowindexUpperTriangle, colindexUpperTriangle, valuesB);
      std::cout<<"\nAddition (symmetric matrix) test: \n";
      C = A + B;
      std::vector<real> valuesC;
      C.get_values(rowindexUpperTriangle, colindexUpperTriangle, valuesC);
      std::vector<real> valuesCRef(valuesA.size(),0);
      for (unsigned int ind = 0; ind < valuesA.size(); ++ind)
	valuesCRef[ind] = valuesA[ind] + valuesB[ind];
      std::cout<<" Absolute Difference: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<maxdiff(valuesCRef, valuesC);
      if (maxdiff(valuesCRef, valuesC) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      } 
    } /* End of addition test symmetric matrix */



    /* *********  trsytriplemm test 1 : A = Z' * A * Z
     *            where A is symmetric and Z upper triangular
     * This also tests the routines  gemm_upper_tr_only, 
     * trmm_upper_tr_only, and sytr_upper_tr_only which are called from 
     * trsytriplemm.
     *
     */ {
      real ONE = 1.0;
      symmMatrix syA;
      triangMatrix trZ;
      syA.resetSizesAndBlocks(rows,cols);
      syA.randomZeroStructure(0.7);
      //      syA.random();
      normalMatrix noA(syA);
      std::vector<real> valuesAblas;
      noA.get_values(rowindex, colindex, valuesAblas);
      trZ.resetSizesAndBlocks(rows,cols);
      trZ.random();
      std::vector<real> valuesZPacked;
      trZ.get_values(rowindexUpperTriangle, colindexUpperTriangle, 
		     valuesZPacked);
      std::vector<real> valuesZ(size * size,0);
      tripackedtofull(&valuesZPacked[0], &valuesZ[0], size);
      std::cout<<"\nTriple matrix multiply test 1: \n";
      trmm("L", "U", "T", "N", &size, &size, &ONE, &valuesZ[0],
	   &size, &valuesAblas[0], &size);
      trmm("R", "U", "N", "N", &size, &size, &ONE, &valuesZ[0],
	   &size, &valuesAblas[0], &size);
      syA = transpose(trZ) * syA * trZ;
      std::vector<real> valuesA;
      noA = syA;
      noA.get_values(rowindex, colindex, valuesA);
      std::cout<<" Max Absolute Error: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<maxdiff(valuesA, valuesAblas);
      if (maxdiff(valuesA, valuesAblas) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      } 
    } /* End of trsytriplemm test 1 */

    /* *********  trsytriplemm test 2 : A = Z * A * Z'    
     *            where A is symmetric and Z upper triangular
     * This also tests the routines  gemm_upper_tr_only, 
     * trmm_upper_tr_only, and sytr_upper_tr_only which are called from 
     * trsytriplemm.
     *
     */ {
      real ONE = 1.0;
      symmMatrix syA;
      triangMatrix trZ;
      syA.resetSizesAndBlocks(rows,cols);
      syA.randomZeroStructure(0.7);
      //      syA.random();
      normalMatrix noA(syA);
      std::vector<real> valuesAblas;
      noA.get_values(rowindex, colindex, valuesAblas);
      trZ.resetSizesAndBlocks(rows,cols);
      trZ.random();
      std::vector<real> valuesZPacked;
      trZ.get_values(rowindexUpperTriangle, colindexUpperTriangle, 
		     valuesZPacked);
      std::vector<real> valuesZ(size * size,0);
      tripackedtofull(&valuesZPacked[0], &valuesZ[0], size);
      std::cout<<"\nTriple matrix multiply test 2: \n";
      trmm("L", "U", "N", "N", &size, &size, &ONE, &valuesZ[0],
	   &size, &valuesAblas[0], &size);
      trmm("R", "U", "T", "N", &size, &size, &ONE, &valuesZ[0],
	   &size, &valuesAblas[0], &size);
      syA = trZ * syA * transpose(trZ);
      std::vector<real> valuesA;
      noA = syA;
      noA.get_values(rowindex, colindex, valuesA);
      std::cout<<" Max Absolute Error: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<maxdiff(valuesA, valuesAblas);
      if (maxdiff(valuesA, valuesAblas) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      } 
      
    } /* End of trsytriplemm test 2 */


    /********* Test Frobenius thresholding */
    {
      normalMatrix noA, noB;
      noA.resetSizesAndBlocks(rows,cols);
      noA.randomZeroStructure(0.1);
      //      noA.random();

      real thr;
      int ntests = 4;
      real frobE;
      real ONE = 1.0;
      for (int nr = ntests; nr > 0; --nr) {
	thr = pow((real)10,-nr+2);
	noB = noA;
	noB.frob_thresh(thr);
	noB += -ONE * noA;
	frobE = noB.frob();
	std::cout<<"\nFrobenius thresholding test "<<ntests-nr+1
		 <<": \n Requested threshold:"
		 <<std::setprecision(10)<<std::setw(15)
		 <<thr
		 <<"\n Received Error:     "
		 <<std::setprecision(10)<<std::setw(15)
		 <<frobE;
	if (frobE < thr)
	  std::cout<<"   OK" <<std::endl;
	else {
	  std::cout<<"   ERROR" <<std::endl;
	  std::exit(1);
	}
      }
    }
    
    
    /* *********  Inverse Cholesky test: symmetric matrix     */ 
    if (!realIsSingle<real>()) {
      
      symmMatrix syA;
      triangMatrix trZ;
      syA.resetSizesAndBlocks(rows,cols);
      syA.randomZeroStructure(0.7);
      //      syA.random();
      syA.add_identity(100.0);
      std::vector<real> valuesZPackedLapack;
      syA.get_values(rowindexUpperTriangle, 
		     colindexUpperTriangle, 
		     valuesZPackedLapack);
      //      for (int ind = 0; ind < valuesZPackedLapack.size();ind++)
      //	std::cout<< valuesZPackedLapack[ind] << "     ";
      std::cout<<"\nInverse Cholesky test: \n";
      pptrf("U",&size,&valuesZPackedLapack[0],&info); 
      if (info)
	{std::cout<<"Error in Lapack: pptrf info="<<info<<std::endl;
	  std::exit(1);}
      tptri("U","N",&size,&valuesZPackedLapack[0],&info);
      if (info)
	{std::cout<<"Error in Lapack: tptri info="<<info<<std::endl;
	  std::exit(1);}  

      //      tripackedtofull(overpacked,inchfull,size);
      trZ.inch(syA,0);
      std::vector<real> valuesZPacked;
      trZ.get_values(rowindexUpperTriangle, 
		     colindexUpperTriangle, 
		     valuesZPacked);
      std::cout<<" Max Absolute Error: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<maxdiff(valuesZPacked, valuesZPackedLapack);
      if (maxdiff(valuesZPacked, valuesZPackedLapack) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      } 

    } /* End of Inverse Cholesky test symmetric matrix */
    else
      std::cout<<"\nSkipping inverse Cholesky test since lapack \n"
	       <<"single precision Cholesky (spptrf) is unstable. \n";
    /* Turning off single precision test because lapack Cholesky (pptrf) 
     * seems to be broken for single precision. 
     */
      

    /* *********  Gershgorin test: symmetric matrix     */ {
      symmMatrix syA;
      syA.resetSizesAndBlocks(rows,cols);
      syA.randomZeroStructure(0.2);
      //      syA.random();
      normalMatrix noA(syA);
      std::vector<real> valuesA;
      noA.get_values(rowindex, colindex, valuesA);
      /* Compute reference */
      real refMin;
      real refMax;
      std::vector<real> colAbsSums(size);
      for (int col = 0; col < size; ++col) {
	colAbsSums[col] = 0;
	for (int row = 0; row < size; ++row)
	  colAbsSums[col] += template_blas_fabs(valuesA[col*size+row]);
      }
      std::vector<real> diag(size);
      for (int rc = 0; rc < size; ++rc) 
	diag[rc] = valuesA[rc*size+rc];
      real tmp1 = colAbsSums[0] - template_blas_fabs(diag[0]);
      real tmp2;
      refMin = diag[0] - tmp1;
      refMax = diag[0] + tmp1;
      for (int col = 1; col < size; col++) {
	tmp1 = colAbsSums[col] - template_blas_fabs(diag[col]);
	tmp2 = diag[col] - tmp1;
	refMin = (tmp2 < refMin ? tmp2 : refMin);
	tmp2 = diag[col] + tmp1;
	refMax = (tmp2 > refMax ? tmp2 : refMax);
      }
      std::cout<<"\nGershgorin eigenvalue bounds test: \n";
      real gersMin;
      real gersMax;
      syA.gershgorin(gersMin, gersMax);
      std::cout
	<<"Reference: "
	<<std::setprecision(10)<<std::setw(15)
	<<refMin<<", "<<refMax<<std::endl
	<<" Computed: "
	<<std::setprecision(10)<<std::setw(15)
	<<gersMin<<", "<<gersMax;	
      if (template_blas_fabs(refMin-gersMin) < 100*epsilon &&
	  template_blas_fabs(refMax-gersMax) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      } 
    } /* End of Gershgorin test symmetric matrix */
    
    /* *********  Memory usage     */ {
      symmMatrix syA;
      syA.resetSizesAndBlocks(rows,cols);
      syA.random();
      normalMatrix noA(syA);
      std::cout<<"\nMemory usage function test: \n";
      std::cout
	<<"General matrix memory usage  : "
	<<std::setprecision(10)<<std::setw(15)
	<<noA.memory_usage()<<" bytes";
      if (noA.memory_usage() == sizeof(real) * size * size)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      } 
      std::cout
	<<"Symmetric matrix memory usage: "
	<<std::setprecision(10)<<std::setw(15)
	<<syA.memory_usage()<<" bytes";
      if (syA.memory_usage() >= sizeof(real) * (size+1)*size/2 &&
	  syA.memory_usage() <= sizeof(real) * size * size)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      } 
    } /* End of memory usage test */


    /* *********  Accumulation symm     */ {
      symmMatrix syA;
      syA.resetSizesAndBlocks(rows,cols);
      syA.randomZeroStructure(0.7);
      //      syA.random();
      normalMatrix noA(syA);
      std::vector<real> valuesA;
      noA.get_values(rowindex, colindex, valuesA);
      testAccumulation(syA, size, &valuesA[0]);
    } /* End of accumulation test */
    

    /******* Test of vectors ***************************************/
    
    /**** Vectors: construct, copy, and assign. */ {
      std::cout<<"\nVector contructors, copy, and assignment test: ";
      normalVector a, b;
      a.resetSizesAndBlocks(rows);
      a.rand();
      b = a;             /* Test assignment.         */
      normalVector c(a); /* Test copy constructor.   */
      std::vector<real> full_a;
      a.fullvector(full_a);
      std::vector<real> full_b;
      b.fullvector(full_b);
      std::vector<real> full_c;
      c.fullvector(full_c);
      assert(full_a.size() == full_b.size());
      assert(full_a.size() == full_c.size());
      for (unsigned int ind = 0; ind < full_a.size(); ++ind) {
	assert(full_a[ind] == full_b[ind]);
	assert(full_a[ind] == full_c[ind]);
      }
      std::cout<<"OK" <<std::endl;
    }  /*** End: Vectors: construct, copy, and assign. */

    /*** Test dot product */ {
      normalVector a, b;
      a.resetSizesAndBlocks(rows);
      a.rand(); /* normalized rand */
      b = a;
      real res = transpose(a) * b;
      std::cout<<"\nDot product test: \n Error:              "
	       <<std::setprecision(10)<<std::setw(15)
	       <<template_blas_fabs(res - 1.0);
      if (template_blas_fabs(res - 1.0) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      }
    }  /*** End dot product test */
    
    /******* Test euclidean norm */ {
      /* FIXME: This test is kind of stupid since it uses rand which uses eucl.
       */
      normalVector a;
      a.resetSizesAndBlocks(rows);
      a.rand(); /* normalized rand */
      real res = a.eucl();
      std::cout<<"\nEuclidean norm (vector) test: \n Error:              "
	       <<std::setprecision(10)<<std::setw(15)
	       <<template_blas_fabs(res - 1.0);
      if (template_blas_fabs(res - 1.0) < 0.0001)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      }
    } /**** Test euclidean norm */

    /*** Test axpy */ {
      real alpha = 4.59659;
      normalVector a, b;
      a.resetSizesAndBlocks(rows);
      a.rand(); /* normalized rand */
      std::vector<real> valuesa;
      a.fullvector(valuesa);
      b.resetSizesAndBlocks(rows);
      b.rand(); /* normalized rand */
      std::vector<real> valuesbRef;
      b.fullvector(valuesbRef);
      std::vector<real> valuesb;
      b += alpha * a;
      b.fullvector(valuesb);
      for (unsigned int ind = 0; ind < valuesbRef.size(); ++ind) {
	valuesbRef[ind] += alpha * valuesa[ind];
      }
      
      std::cout<<"\naxpy test: \n Error:              "
	       <<std::setprecision(10)<<std::setw(15)
	       <<template_blas_fabs(maxdiff(valuesb,valuesbRef));
      if (template_blas_fabs(maxdiff(valuesb,valuesbRef)) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      }
    }  /*** End axpy test */
    


    /********* Test general matrix-vector product */ {
      int  const ONEint = 1;
      real const ONEreal = 1.0;
      //      int  const ZEROint = 0;
      real const ZEROreal = 0.0;
      real alpha = 2.3456;
      real beta = 6.3467;
      /* M-V-prod 1 */ {
	normalMatrix A;
	A.resetSizesAndBlocks(rows,cols);
	A.randomZeroStructure(0.7);
	//	A.random();
	std::vector<real> valuesA;
	A.get_values(rowindex, colindex, valuesA);

	normalVector b,c;
	b.resetSizesAndBlocks(rows);
	b.rand();
	std::vector<real> valuesb;
	b.fullvector(valuesb);
	std::vector<real> valuescBlas(valuesb.size());
	std::vector<real> valuesc;	
	gemv("N", &size, &size, &alpha, &valuesA[0], &size, &valuesb[0], 
	     &ONEint, &ZEROreal, &valuescBlas[0], &ONEint);
	c = alpha * A * b;
	c.fullvector(valuesc);
	std::cout<<"\nMatrix-vector product test 1: \n Max Absolute Error: "
		 <<std::setprecision(10)<<std::setw(15)
		 <<maxdiff(valuesc, valuescBlas);
	if (maxdiff(valuesc, valuescBlas) < 100*epsilon)
	  std::cout<<"   OK" <<std::endl;
	else {
	  std::cout<<"   ERROR" <<std::endl;
	  std::exit(1);
	}
      } /* End: M-V-prod 1 */

      /* M-V-prod 2 */ {
	normalMatrix A;
	A.resetSizesAndBlocks(rows,cols);
	A.randomZeroStructure(0.7);
	//	A.random();
	std::vector<real> valuesA;
	A.get_values(rowindex, colindex, valuesA);

	normalVector b,c;
	b.resetSizesAndBlocks(rows);
	b.rand();
	c.resetSizesAndBlocks(rows);
	c.rand();
	std::vector<real> valuesb;
	b.fullvector(valuesb);
	std::vector<real> valuescBlas;
	c.fullvector(valuescBlas);
	std::vector<real> valuesc;	
	gemv("T", &size, &size, &alpha, &valuesA[0], &size, &valuesb[0], 
	     &ONEint, &ONEreal, &valuescBlas[0], &ONEint);
	c += alpha * transpose(A) * b;
	c.fullvector(valuesc);
	
	std::cout<<"\nMatrix-vector product test 2: \n Max Absolute Error: "
		 <<std::setprecision(10)<<std::setw(15)
		 <<maxdiff(valuesc,valuescBlas);
	if (maxdiff(valuesc,valuescBlas) < 100*epsilon)
	  std::cout<<"   OK" <<std::endl;
	else {
	  std::cout<<"   ERROR" <<std::endl;
	  std::exit(1);
	}
      } /* End: M-V-prod 2 */

      /* M-V-prod 3 */ {
	normalMatrix A;
	A.resetSizesAndBlocks(rows,cols);
	A.randomZeroStructure(0.7);
	//	A.random();
	std::vector<real> valuesA;
	A.get_values(rowindex, colindex, valuesA);

	normalVector b,c;
	b.resetSizesAndBlocks(rows);
	b.rand();
	c.resetSizesAndBlocks(rows);
	c.rand();
	std::vector<real> valuesb;
	b.fullvector(valuesb);
	std::vector<real> valuescBlas;
	c.fullvector(valuescBlas);
	std::vector<real> valuesc;	
	gemv("T", &size, &size, &alpha, &valuesA[0], &size, &valuesb[0], 
	     &ONEint, &beta, &valuescBlas[0], &ONEint);
	c = alpha * transpose(A) * b + beta * c;
	c.fullvector(valuesc);

    	std::cout<<"\nMatrix-vector product test 3: \n Max Absolute Error: "
		 <<std::setprecision(10)<<std::setw(15)
		 <<maxdiff(valuesc,valuescBlas);
	if (maxdiff(valuesc,valuescBlas) < 100*epsilon)
	  std::cout<<"   OK" <<std::endl;
	else {
	  std::cout<<"   ERROR" <<std::endl;
	  std::exit(1);
	}
      } /* End: M-V-prod 2 */
    } /* End: Test general matrix-vector product */
    
    /******************* Test symmetric matrix-vector product */ {
      int  const ONEint = 1;
      real const ONEreal = 1.0;
      // int  const ZEROint = 0;
      real const ZEROreal = 0.0;
      real alpha = 2.3456;
      real beta = 6.3467;
      
      /* Symm M-V-prod 1 */ {
	symmMatrix syA;
	syA.resetSizesAndBlocks(rows,cols);
	syA.randomZeroStructure(0.7);
	//	syA.random();	
	normalMatrix noA(syA);
	std::vector<real> valuesA;
	noA.get_values(rowindex, colindex, valuesA);
	
	normalVector b,c;
	b.resetSizesAndBlocks(rows);
	b.rand();
	std::vector<real> valuesb;
	b.fullvector(valuesb);
	std::vector<real> valuescBlas(valuesb.size());
	std::vector<real> valuesc;	

	symv("U", &size, &alpha, &valuesA[0], &size, &valuesb[0], &ONEint, 
	     &ZEROreal, &valuescBlas[0], &ONEint);
	c = alpha * syA * b;
	c.fullvector(valuesc);
    	std::cout<<"\nSymm Matrix-vector product test 1: \n Max Absolute Error: "
		 <<std::setprecision(10)<<std::setw(15)
		 <<maxdiff(valuesc,valuescBlas);
	if (maxdiff(valuesc,valuescBlas) < 100*epsilon)
	  std::cout<<"   OK" <<std::endl;
	else {
	  std::cout<<"   ERROR" <<std::endl;
	  std::exit(1);
	}
      } /* End: Symm M-V-prod 1 */

      /* Symm M-V-prod 2 */ {
	symmMatrix syA;
	syA.resetSizesAndBlocks(rows,cols);
	syA.randomZeroStructure(0.7);
	//	syA.random();	
	normalMatrix noA(syA);
	std::vector<real> valuesA;
	noA.get_values(rowindex, colindex, valuesA);
	
	normalVector b,c;
	b.resetSizesAndBlocks(rows);
	b.rand();
	c.resetSizesAndBlocks(rows);
	c.rand();
	std::vector<real> valuesb;
	b.fullvector(valuesb);
	std::vector<real> valuescBlas;
	c.fullvector(valuescBlas);
	std::vector<real> valuesc;	
	symv("U", &size, &alpha, &valuesA[0], &size, &valuesb[0], &ONEint, 
	     &ONEreal, &valuescBlas[0], &ONEint);
	c += alpha * syA * b;
	c.fullvector(valuesc);
	
	std::cout<<"\nSymm Matrix-vector product test 2: \n Max Absolute Error: "
		 <<std::setprecision(10)<<std::setw(15)
		 <<maxdiff(valuesc,valuescBlas);
	if (maxdiff(valuesc,valuescBlas) < 100*epsilon)
	  std::cout<<"   OK" <<std::endl;
	else {
	  std::cout<<"   ERROR" <<std::endl;
	  std::exit(1);
	}
      } /* End: Symm M-V-prod 2 */

      /* Symm M-V-prod 3 */ {
	symmMatrix syA;
	syA.resetSizesAndBlocks(rows,cols);
	syA.randomZeroStructure(0.7);
	//	syA.random();	
	normalMatrix noA(syA);
	std::vector<real> valuesA;
	noA.get_values(rowindex, colindex, valuesA);
	
	normalVector b,c;
	b.resetSizesAndBlocks(rows);
	b.rand();
	c.resetSizesAndBlocks(rows);
	c.rand();
	std::vector<real> valuesb;
	b.fullvector(valuesb);
	std::vector<real> valuescBlas;
	c.fullvector(valuescBlas);
	std::vector<real> valuesc;	
	symv("U", &size, &alpha, &valuesA[0], &size, &valuesb[0], &ONEint, 
	     &beta, &valuescBlas[0], &ONEint);
	c = alpha * syA * b + beta * c;
	c.fullvector(valuesc);
	
	std::cout<<"\nSymm Matrix-vector product test 3: \n Max Absolute Error: "
		 <<std::setprecision(10)<<std::setw(15)
		 <<maxdiff(valuesc,valuescBlas);
	if (maxdiff(valuesc,valuescBlas) < 100*epsilon)
	  std::cout<<"   OK" <<std::endl;
	else {
	  std::cout<<"   ERROR" <<std::endl;
	  std::exit(1);
	}
      } /* End: Symm M-V-prod 3 */
    } /* End: Test symmetric matrix-vector product */

    /* Triangular matrix-vector product test 1 */ {
      int  const ONEint = 1;
      triangMatrix trA;
      trA.resetSizesAndBlocks(rows,cols);
      //trA.randomZeroStructure(0.7);
      trA.random();
      normalMatrix noA(trA);
      std::vector<real> valuesA;
      noA.get_values(rowindex, colindex, valuesA);
#if 0
      std::cout << "TRIANG EUCL NORM: " << trA.eucl(1e-4,400) <<std::endl;
      std::cout << "A=[\n";
      for (int rowi = 0; rowi < size; ++rowi) {
	for (int coli = 0;coli < size; coli++) 
	  std::cout << valuesA[rowi+coli*size] << "   ";
	std::cout << std::endl;
      }
      std::cout << "]; norm(A,2)\n";
#endif 
      normalVector x;
      x.resetSizesAndBlocks(rows);
      x.rand();
      std::vector<real> valuesx;
      std::vector<real> valuesxBlas;
      x.fullvector(valuesx);
      x.fullvector(valuesxBlas);
      trmv("U","N","N", &size, &valuesA[0], &size, &valuesxBlas[0], &ONEint);
      x = trA * x;
      x.fullvector(valuesx);
      std::cout<<"\nTriang Matrix-vector product test 1: \n Max Absolute Error: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<maxdiff(valuesx,valuesxBlas);
      if (maxdiff(valuesx,valuesxBlas) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      }
    } /* End: Triang M-V-prod 1 */

    /* Triangular matrix-vector product test 2 */ {
      int  const ONEint = 1;
      triangMatrix trA;
      trA.resetSizesAndBlocks(rows,cols);
      //trA.randomZeroStructure(0.7);
      trA.random();	
      normalMatrix noA(trA);
      std::vector<real> valuesA;
      noA.get_values(rowindex, colindex, valuesA);
      normalVector x;
      x.resetSizesAndBlocks(rows);
      x.rand();
      std::vector<real> valuesx;
      std::vector<real> valuesxBlas;
      x.fullvector(valuesx);
      x.fullvector(valuesxBlas);
      trmv("U","T","N", &size, &valuesA[0], &size, &valuesxBlas[0], &ONEint);
      x = transpose(trA) * x;
      x.fullvector(valuesx);
      std::cout<<"\nTriang Matrix-vector product test 2: \n Max Absolute Error: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<maxdiff(valuesx,valuesxBlas);
      if (maxdiff(valuesx,valuesxBlas) < 100*epsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      }
    } /* End: Triang M-V-prod 2 */


    if (!realIsSingle<real>()) { /********** Test Lanczos */ 
      real const ONEreal = 1.0;
      real lanEpsilon = epsilon*1e-1;
      symmMatrix syA;
      syA.resetSizesAndBlocks(rows,cols);
      syA.randomZeroStructure(0.3);
      //      syA.random();	
      normalVector x;
      x.resetSizesAndBlocks(rows);
      x.rand();
      int maxit = 400;
      myLanczosType lan(syA, x, maxit);
      lan.setAbsTol( lanEpsilon );
      lan.run();
      real eigenValue = 0;
      real accuracy = 0;
      lan.getLargestMagnitudeEig(eigenValue, accuracy);
      
      normalVector eigenVec;
      lan.getLargestMagnitudeEigPair(eigenValue, eigenVec, accuracy);
      normalVector resVec(eigenVec);
      resVec *= eigenValue;
      resVec += -ONEreal * syA * eigenVec;
      std::cout<<"\nLanczos largest magnitude test : \n" 
	       <<" Requested accuracy: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<lanEpsilon
	       <<"\n Indicated Error:    "
	       <<std::setprecision(10)<<std::setw(15)
	       <<accuracy
	       <<"\n Actual Error:       "
	       <<std::setprecision(10)<<std::setw(15)
	       <<resVec.eucl();
      if (accuracy < lanEpsilon && 
	  (resVec.eucl() < accuracy || 
	   resVec.eucl() < mat::getRelPrecision<real>() * 100))
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      }    
    } /***** End Lanczos test */
    else
      std::cout<<"\nSkipping Lanczos test for single precision. \n";
    /* Lanczos is difficult in single precision. 
     * Would probably work better with non-random matrix. 
     */


    if (!realIsSingle<real>()) { /**** Test Lanczos again for a
				       difficult matrix (this test
				       added by Elias 2010-03-20) */ 
      for(int kk = 0; kk < 20; kk++) {
	real const ONEreal = 1.0;
	symmMatrix syA;
	int sizeTmp = 5;
	int nlevelsTmp = 3;
	std::vector<int> blockSizesTmp(nlevelsTmp);
	blockSizesTmp[nlevelsTmp - 1] = 1;
	blockSizesTmp[nlevelsTmp - 2] = 1;
	blockSizesTmp[nlevelsTmp - 3] = 5;
	SizesAndBlocks rowsTmp(blockSizesTmp, sizeTmp);
	SizesAndBlocks colsTmp(blockSizesTmp, sizeTmp);
	syA.resetSizesAndBlocks(rowsTmp,colsTmp);
	// The 5 by 5 test matrix is zero except three values on the
	// diagonal. This particular matrix appeared as the X-X2 matrix
	// in purification in one of the ergo tests.
	std::vector<real> fullMat(5*5);
	for(int i = 0; i < 5*5; i++)
	  fullMat[i] = 0;
	fullMat[2*5+2] = 0.062498936885316151713;
	fullMat[3*5+3] = 0.062498936885437686439;
	fullMat[4*5+4] = 0.062498936885316151713;
	real correctNormValue = 0.062498936885437686439;
	syA.assignFromFull(fullMat);
	normalVector x;
	x.resetSizesAndBlocks(rowsTmp);
	x.rand();
	int maxit = 400;
	real requestedAccuracy = 1e-13; // The value 5.55e-14 was used in purification.
	myLanczosType lan(syA, x, maxit);
	lan.setAbsTol( requestedAccuracy );
      	lan.run();
	real eigenValue = 0;
	real accuracy = 0;
	lan.getLargestMagnitudeEig(eigenValue, accuracy);
	normalVector eigenVec;
	lan.getLargestMagnitudeEigPair(eigenValue, eigenVec, accuracy);
	normalVector resVec(eigenVec);
	resVec *= eigenValue;
	resVec += -ONEreal * syA * eigenVec;
	std::cout<<"\nLanczos largest magnitude test again (" << kk << ") : \n" 
		 <<  " Requested accuracy:           "
		 <<std::setprecision(10)<<std::setw(15)
		 <<requestedAccuracy
		 <<"\n Indicated Error:              "
		 <<std::setprecision(10)<<std::setw(15)
		 <<accuracy
		 <<"\n 'Actual Error' (resVec norm): "
		 <<std::setprecision(10)<<std::setw(15)
		 <<resVec.eucl()
		 <<"\n Actual Error (really):        "
		 <<std::setprecision(10)<<std::setw(15)
		 <<template_blas_fabs(correctNormValue-eigenValue);
	if (accuracy < requestedAccuracy && 
	    (resVec.eucl() < accuracy || 
	     resVec.eucl() < mat::getRelPrecision<real>() * 100) &&
	    template_blas_fabs(correctNormValue-eigenValue) < requestedAccuracy)
	  std::cout<<"   OK" <<std::endl;
	else {
	  std::cout<<"   ERROR" <<std::endl;
	  std::exit(1);
	}
      } // end for kk
    } /***** End Lanczos test */
    else
      std::cout<<"\nSkipping Lanczos test for single precision. \n";
    /* Lanczos is difficult in single precision. 
     * Would probably work better with non-random matrix. 
     */

    if (!realIsSingle<real>()) { /********** Test Lanczos for diffmatrix */ 
      //      real const ONEreal = 1.0;
      real lanEpsilon = epsilon*1e-1;
      symmMatrix syA;
      syA.resetSizesAndBlocks(rows,cols);
      syA.randomZeroStructure(0.3);
      symmMatrix syB;
      syB.resetSizesAndBlocks(rows,cols);
      std::vector<int> rowindex(1);    rowindex[0] = 0;
      std::vector<int> colindex(1);    colindex[0] = 0;
      std::vector<real> values_tmp(1); values_tmp[0] = 1;
      symmMatrix tmp;
      tmp.assign_from_sparse(rowindex, colindex, values_tmp, rows, cols);
      syB = syA + tmp;
      normalVector x;
      x.resetSizesAndBlocks(rows);
      x.rand();
      int maxit = 400;
      DiffMatrix<symmMatrix,real> Diff(syA,syB);
      arn::LanczosLargestMagnitudeEig
	<real, DiffMatrix<symmMatrix,real>, normalVector>
	lan(Diff, x, maxit);
      lan.setAbsTol( lanEpsilon );
      lan.run();
      real eigenValue = 0;
      real accuracy = 0;
      lan.getLargestMagnitudeEig(eigenValue, accuracy);
      
      std::cout<<"\nLanczos largest magnitude test (Diff matrix) : \n" 
	       <<" Requested accuracy: "
	       <<std::setprecision(10)<<std::setw(15)
	       <<lanEpsilon
	       <<"\n Indicated Error:    "
	       <<std::setprecision(10)<<std::setw(15)
	       <<accuracy
	       <<"\n Actual Error:       "
	       <<std::setprecision(10)<<std::setw(15)
	       <<eigenValue+1;
      if (accuracy < lanEpsilon)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      }    
    } /***** End Lanczos test */
    else
      std::cout<<"\nSkipping Lanczos (Diff) test for single precision. \n";
    


    /********* Test Euclidean thresholding */ {
    real thr = 0;
    int ntests = 4;
    real ONE = 1.0;
    symmMatrix syA;
    syA.resetSizesAndBlocks(rows,cols);

#if 1
    expRule<real> er;
    syA.setElementsByRule(er);
#else
    syA.randomZeroStructure(0.1);
    //    syA.random();	
#endif

    symmMatrix syB;

    for (int nr = ntests; nr > 0; --nr) {
      thr = pow((real)10,-nr+2);
      syB = syA;
      syB.eucl_thresh(thr);
      syB += -ONE * syA;
      std::cout<<"\nEuclidean  truncation test, symmetric matrices "<<ntests-nr+1
	       <<": \n Requested threshold:"
	       <<std::setprecision(10)<<std::setw(15)
	       <<thr
	       <<"\n Received Error:     "
	       <<std::setprecision(10)<<std::setw(15)
	       <<syB.eucl(thr * 1e-7);
      
      if (syB.eucl(thr * 1e-3) < thr)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      }
    }
    } /*** End Euclidean thresholding test */

    /********* Test mixed norm thresholding */ {
      real thr = 0;
      int ntests = 4;
      real ONE = 1.0;
      symmMatrix syA;
      syA.resetSizesAndBlocks(rows,cols);

#if 1
      expRule<real> er;
      syA.setElementsByRule(er);
#else
      syA.randomZeroStructure(0.1);
      //    syA.random();	
#endif

      symmMatrix syB;

      for (int nr = ntests; nr > 0; --nr) {
	thr = pow((real)10,-nr+2);
	syB = syA;
	std::cout <<"mixed norm res: " << syB.mixed_norm_thresh(thr) << std::endl;
	syB += -ONE * syA;
	std::cout<<"\nMixed norm truncation test, symmetric matrices "<<ntests-nr+1
		 <<": \n Requested threshold:"
		 <<std::setprecision(10)<<std::setw(15)
		 <<thr
		 <<"\n Received Error (Euclidean norm):     "
		 <<std::setprecision(10)<<std::setw(15)
		 <<syB.eucl(thr * 1e-7);
      
	if (syB.eucl(thr * 1e-3) < thr)
	  std::cout<<"   OK" <<std::endl;
	else {
	  std::cout<<"   ERROR" <<std::endl;
	  std::exit(1);
	}
      }
    } /*** End mixed norm thresholding test */

#endif
    
    /********* Test Euclidean thresholding, triangular matrices */ {
      real thr = 0;
      int ntests = 4;
      real ONE = 1.0;
      triangMatrix trA;
      trA.resetSizesAndBlocks(rows,cols);
      expRule<real> er;
      trA.setElementsByRule(er);      
      triangMatrix trB;
      
      for (int nr = ntests; nr > 0; --nr) {
	thr = pow((real)10,-nr+2);
	trB = trA;
	trB.eucl_thresh(thr);
	trB += -ONE * trA;
	std::cout<<"\nEuclidean truncation test, triangular matrices "<<ntests-nr+1
		 <<": \n Requested threshold:"
		 <<std::setprecision(10)<<std::setw(15)
		 <<thr
		 <<"\n Received Error:     "
		 <<std::setprecision(10)<<std::setw(15)
		 <<trB.eucl(1e-4);
      
	if (trB.eucl(thr * 1e-2) < thr)
	  std::cout<<"   OK" <<std::endl;
	else {
	  std::cout<<"   ERROR" <<std::endl;
	  std::exit(1);
	}
      }
    } /*** End Euclidean thresholding test, triangular matrices */

    /********* Test Euclidean thresholding, general matrices */ {
      real thr = 0;
      int ntests = 4;
      real ONE = 1.0;
      normalMatrix noA;
      noA.resetSizesAndBlocks(rows,cols);
      expRule<real> er;
      noA.setElementsByRule(er);      
      normalMatrix noB;
      
      for (int nr = ntests; nr > 0; --nr) {
	thr = pow((real)10,-nr+2);
	noB = noA;
	noB.eucl_thresh(thr);
	noB += -ONE * noA;
	std::cout<<"\nEuclidean truncation test, general matrices "<<ntests-nr+1
		 <<": \n Requested threshold:"
		 <<std::setprecision(10)<<std::setw(15)
		 <<thr
		 <<"\n Received Error:     "
		 <<std::setprecision(10)<<std::setw(15)
		 <<noB.eucl(1e-4);
      
	if (noB.eucl(thr * 1e-2) < thr)
	  std::cout<<"   OK" <<std::endl;
	else {
	  std::cout<<"   ERROR" <<std::endl;
	  std::exit(1);
	}
      }
    } /*** End Euclidean thresholding test, general matrices */

#if 1
    /* *********  Test truncation using congr trans measure
     *
     * Truncation: Zt = Z + E
     * Measure:     ||E^T * B * E + E^T * B * A + A^T * B * E||_2 < threshold
     *
     */ {
      
      real ONE = 1.0;
      symmMatrix syA;
      syA.resetSizesAndBlocks(rows,cols);
      syA.randomZeroStructure(0.7);
      symmMatrix syARef(syA);
      triangMatrix trZ;
      trZ.resetSizesAndBlocks(rows,cols);
      expRule<real> er;
      trZ.setElementsByRule(er);      
      triangMatrix trZCopy(trZ);

      std::cout<<"\nTest truncation with triple matrix measure: \n";
      syARef = transpose(trZ) * syARef * trZ;
      real thr = 0.01;
      real reported_error = trZ.eucl_thresh_congr_trans_measure(thr, syA);
      syA = transpose(trZ) * syA * trZ;
      syA += -ONE * syARef;
      real err = syA.eucl(1e-7);
      
      std::cout<<"\n Requested threshold:"
	       <<std::setprecision(10)<<std::setw(15)
	       <<thr
	       <<"\n Received Error:     "
	       <<std::setprecision(10)<<std::setw(15)
	       <<err
	       <<"\n Reported Error:     "
	       <<std::setprecision(10)<<std::setw(15)
	       <<reported_error;
      
      if (err < thr)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      }

    } /* End of truncation using triple matrix as measure */

    /* *********  Test truncation using congr trans measure again
     *
     * Truncation: Zt = Z + E
     * Measure:     ||E^T * B * E + E^T * B * A + A^T * B * E||_2 < threshold
     *
     */ {

      /********** Initialization of SizesAndBlocks                       */
      int size = 87; /* Use weird size to find more bugs. */
      int nlevels = 3;
      std::vector<int> blockSizes(nlevels);
      blockSizes[nlevels - 1] = 1;
      blockSizes[nlevels - 2] = 1;
      blockSizes[nlevels - 3] = 32;

      std::cout << "Running tests with blocksize vector: ";
      for (int ind = 0; ind < nlevels; ind++)
	std::cout << blockSizes[ind] << "  ";
      std::cout << std::endl;
      
      SizesAndBlocks rows2(blockSizes, size);
      SizesAndBlocks cols2(blockSizes, size);
      // Now we have rows2 and cols2, now we do the test again
      real ONE = 1.0;
      symmMatrix syA;
      syA.resetSizesAndBlocks(rows2,cols2);
      syA.randomZeroStructure(0.7);
      symmMatrix syARef(syA);
      triangMatrix trZ;
      trZ.resetSizesAndBlocks(rows2,cols2);
      expRule<real> er;
      trZ.setElementsByRule(er);      
      triangMatrix trZCopy(trZ);

      std::cout<<"\nTest truncation with triple matrix measure again: \n";
      syARef = transpose(trZ) * syARef * trZ;
      real thr = 0.01;
      real reported_error = trZ.eucl_thresh_congr_trans_measure(thr, syA);
      syA = transpose(trZ) * syA * trZ;
      syA += -ONE * syARef;
      real err = syA.eucl(1e-7);
      
      std::cout<<"\n Requested threshold:"
	       <<std::setprecision(10)<<std::setw(15)
	       <<thr
	       <<"\n Received Error:     "
	       <<std::setprecision(10)<<std::setw(15)
	       <<err
	       <<"\n Reported Error:     "
	       <<std::setprecision(10)<<std::setw(15)
	       <<reported_error;
      
      if (err < thr)
	std::cout<<"   OK" <<std::endl;
      else {
	std::cout<<"   ERROR" <<std::endl;
	std::exit(1);
      }

      

    } /* End of truncation using triple matrix as measure second time */

#endif


    /*** Test timers in FileWritable */    {
      //      double wallTimeWrite; 
      //      double wallTimeRead;
      //      double wallTimeCopyAndAssign;
      symmMatrix syA;
      syA.resetSizesAndBlocks(rows,cols);
      syA.randomZeroStructure(0.1);
      normalVector x;
      x.resetSizesAndBlocks(rows);
      x.rand();
      std::cout << std::endl;
      std::cout << "Test timers in FileWritable: " << std::endl;
      FileWritable::resetStats();
      std::cout << " Initial values         :  " << std::endl;
      std::cout << "  Write           : " 
		<< FileWritable::getStatsTimeWrite() << std::endl;
      std::cout << "  Read            : " 
		<< FileWritable::getStatsTimeRead() << std::endl;
      std::cout << "  Copy and assign : " 
		<< FileWritable::getStatsTimeCopyAndAssign() << std::endl;
      syA.writeToFile();
      std::cout << " After 1st writeToFile  :  "
		<< std::endl;
      std::cout << "  Write           : " 
		<< FileWritable::getStatsTimeWrite() << std::endl;
      std::cout << "  Read            : " 
		<< FileWritable::getStatsTimeRead() << std::endl;
      std::cout << "  Copy and assign : " 
		<< FileWritable::getStatsTimeCopyAndAssign() << std::endl;
      syA.readFromFile();
      std::cout << " After 1st readFromFile :  " 
		<< std::endl;
      std::cout << "  Write           : " 
		<< FileWritable::getStatsTimeWrite() << std::endl;
      std::cout << "  Read            : " 
		<< FileWritable::getStatsTimeRead() << std::endl;
      std::cout << "  Copy and assign : " 
		<< FileWritable::getStatsTimeCopyAndAssign() << std::endl;
      x.writeToFile();
      std::cout << " After 2nd writeToFile  :  " 
		<< std::endl;
      std::cout << "  Write           : " 
		<< FileWritable::getStatsTimeWrite() << std::endl;
      std::cout << "  Read            : " 
		<< FileWritable::getStatsTimeRead() << std::endl;
      std::cout << "  Copy and assign : " 
		<< FileWritable::getStatsTimeCopyAndAssign() << std::endl;
      x.readFromFile();
      std::cout << " After 2nd readFromFile :  " 
		<< std::endl;
      std::cout << "  Write           : " 
		<< FileWritable::getStatsTimeWrite() << std::endl;
      std::cout << "  Read            : " 
		<< FileWritable::getStatsTimeRead() << std::endl;
      std::cout << "  Copy and assign : " 
		<< FileWritable::getStatsTimeCopyAndAssign() << std::endl;
      
      std::cout << " Total counts :  " << std::endl;
      std::cout << "  Write           : " 
		<< FileWritable::getStatsCountWrite() << std::endl;
      std::cout << "  Read            : " 
		<< FileWritable::getStatsCountRead() << std::endl;
      std::cout << "  Copy and assign : " 
		<< FileWritable::getStatsCountCopyAndAssign() << std::endl;

    }
    
  }
  catch (Failure & e) {
    std::cout << "Failure caught: "<<e.what() << std::endl;
    std::exit(1);
  }
  catch (std::exception & e) {
    std::cout << "Exception caught: "<<e.what() << std::endl;
    std::exit(1);
  }
  return 0;
}

int main(int argc,char* argv[]){
  
  std::cout<<"=========================================================\n"
	   <<"Testing matrix library with double precision:\n"
	   <<"=========================================================\n";
  if (!mainFun<double>(argc,argv)) {
    std::cout
      <<"Matrix library tests with double precision completed successfully!" 
      <<std::endl<<std::endl;
  }
  
  std::cout<<"=========================================================\n"
	   <<"Testing matrix library with single precision:\n"
	   <<"=========================================================\n";
  if (!mainFun<float>(argc,argv))
    std::cout
      <<"Matrix library tests with single precision completed successfully!" 
      <<std::endl<<std::endl;
  std::exit(0);
}


template<class Treal>
static Treal maxdiff(std::vector<Treal> const & f1,
		     std::vector<Treal> const & f2) {
  Treal diff = 0;
  Treal tmpdiff;
  assert(f1.size() == f2.size());
  for(unsigned int i = 0; i < f1.size(); i++) {
    tmpdiff = template_blas_fabs(f1[i] - f2[i]);
    if (tmpdiff > 0) {
      diff = (diff > tmpdiff ? diff : tmpdiff);
    }
  }
  return diff;
}

template<class Treal>
static Treal maxdiff_tri(const Treal* f1,const Treal* f2,int size) {
  Treal diff = 0;
  Treal tmpdiff;
  for (int col = 0; col < size; col++)
    for (int row = 0; row < col + 1; row++) {
      tmpdiff = template_blas_fabs(f1[col * size + row] - f2[col * size + row]);
      diff = (diff > tmpdiff ? diff : tmpdiff);
    }
  return diff;
}


template<class Treal>
static Treal frobdiff(const Treal* f1,const Treal* f2,int size) {
  Treal diff = 0;
  Treal tmp;
  for(int i = 0; i < size * size; i++) {
    tmp = f1[i] - f2[i];
    diff += tmp * tmp;
  }
  return template_blas_sqrt(diff);
}

