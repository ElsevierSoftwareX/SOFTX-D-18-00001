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

/** @file MatrixTridiagSymmetric.h Class for tridiagonal symmetric matrices
 *
 * Copyright(c) Emanuel Rubensson 2006
 *
 * @author Emanuel Rubensson
 * @date December 2006
 *
 */
#ifndef MAT_MATRIXTRIDIAGSYMMETRIC
#define MAT_MATRIXTRIDIAGSYMMETRIC
#include "mat_gblas.h"
namespace mat { /* Matrix namespace */
  namespace arn { /* Arnoldi type methods namespace */
    /** Tridiagonal symmetric matrix class template
     *
     */
    template<typename Treal>
      class MatrixTridiagSymmetric {
    public:
      explicit MatrixTridiagSymmetric(int k = 100) 
	: alphaVec(new Treal[k]), betaVec(new Treal[k]),
        size(0), capacity(k) {}
	void increase(Treal const & alpha, Treal const & beta) {
	  if (size + 1 > capacity)
	    increaseCapacity(capacity * 2);
	  ++size;
	  alphaVec[size - 1] = alpha;
	  betaVec[size - 1] = beta;
	}
	virtual ~MatrixTridiagSymmetric() {
	  delete[] alphaVec;
	  delete[] betaVec;
	}
	

	void update_beta(Treal const & beta)
	{
	  betaVec[size-1] = beta;
	}

	void getEigsByInterval(Treal* eigVals,
			       Treal* eigVectors, 
			       Treal* acc, 
			       int & nEigsFound, 
			       Treal const lowBound,
			       Treal const uppBound,
			       Treal const abstol = 0) const;
	void getEigsByIndex(Treal* eigVals, 
			    Treal* eigVectors, 
			    Treal* acc, 
			    int const lowInd,
			    int const uppInd,
			    Treal const abstol = 0) const;
	inline void clear() {
	  size = 0;
	}
    protected:
	Treal* alphaVec;
	Treal* betaVec;
	int size;
	int capacity;
	void increaseCapacity(int const newCapacity);
    private:
    };

    template<typename Treal>
      void MatrixTridiagSymmetric<Treal>::
      getEigsByInterval(Treal* eigVals, /* length: >= nEigsFound */
			Treal* eigVectors, /* length: >= size * nEigsFound */
			Treal* acc, /* length: size */
			int & nEigsFound, /* The number of found eigenpairs. */
			Treal const lowBound,
			Treal const uppBound,
			Treal const absTol) const {
      Treal* eigArray  = new Treal[size];
      Treal* alphaCopy = new Treal[size];
      Treal* betaCopy  = new Treal[size];
      Treal* work = new Treal[5 * size];
      int* iwork = new int[5 * size];
      int* ifail = new int[size];
      for (int ind = 0; ind < size; ind++){
	alphaCopy[ind] = alphaVec[ind];
	betaCopy[ind] = betaVec[ind];
      }
      int dummy = -1;
      int info;
      /* Find eigenvalues */
      /* FIXME: change to stevr */
      mat::stevx("V", "V", &size, alphaCopy, betaCopy, 
		 &lowBound, &uppBound, &dummy, &dummy,
		 &absTol, 
		 &nEigsFound, eigArray, eigVectors, &size, 
		 work, iwork, ifail, 
		 &info);
      assert(info == 0);
      for (int ind = 0; ind < nEigsFound; ind++) {
	eigVals[ind] = eigArray[ind];
	acc[ind] = betaCopy[size - 1] * 
	  template_blas_fabs(eigVectors[(ind * size) + size - 1]) / 0.9; // FIXME: WHY IS THERE A FACTOR 0.9 HERE ?!?
      }
      delete[] eigArray;
      delete[] alphaCopy;
      delete[] betaCopy;
      delete[] work;
      delete[] iwork;
      delete[] ifail;
    }

    template<typename Treal>
      void MatrixTridiagSymmetric<Treal>::
      getEigsByIndex(Treal* eigVals, /* length: uppInd-lowInd+1 */
		     Treal* eigVectors, /* length: size*(uppInd-lowInd+1) */
		     Treal* acc, /* length: size */
		     int const lowInd,
		     int const uppInd,
		     Treal const abstol) const {
      Treal* eigArray  = new Treal[size];
      Treal* alphaCopy = new Treal[size];
      Treal* betaCopy  = new Treal[size];
      for (int ind = 0; ind < size; ind++){
	alphaCopy[ind] = alphaVec[ind];
	betaCopy[ind] = betaVec[ind];
      }
#if 1
      // Emanuel note 2010-03-14:
      // The following code uses stevr instead of stevx for two reasons: 
      // 1) Due to a bug in LAPACK we previously computed all
      // eigenvalues (see Elias' note below) which turned out to be
      // too time consuming in some cases.
      // 2) Contrary to stevx, stevr should never fail to compute the
      // desired eigenpairs unless there is a bug in the implementation
      // or erroneous input.
      int const lowIndNew(lowInd + 1);
      int const uppIndNew(uppInd + 1);
      int nEigsWanted = uppInd - lowInd + 1;
      int nEigsFound = 0;
      int* isuppz = new int[2*nEigsWanted];
      Treal* work;
      int lwork = -1;
      int* iwork;
      int liwork = -1;
      Treal dummy = -1.0;
      int info;
      // First do a workspace query:
      Treal work_query;
      int iwork_query;
      mat::stevr("V", "I", &size, alphaCopy, betaCopy, 
		 &dummy, &dummy, &lowIndNew, &uppIndNew,
		 &abstol,
		 &nEigsFound, eigArray, eigVectors, &size, 
		 isuppz, 
		 &work_query, &lwork, &iwork_query, &liwork, &info);
      lwork = int(work_query);
      liwork = iwork_query;
      work = new Treal[lwork];
      iwork = new int[liwork];
      mat::stevr("V", "I", &size, alphaCopy, betaCopy, 
		 &dummy, &dummy, &lowIndNew, &uppIndNew,
		 &abstol,
		 &nEigsFound, eigArray, eigVectors, &size, 
		 isuppz, 
		 work, &lwork, iwork, &liwork, &info);
      if (info)
	std::cout << "info = " << info <<std::endl;
      assert(info == 0);
      assert(nEigsFound == nEigsWanted);
      for (int ind = 0; ind < nEigsFound; ind++) {
	eigVals[ind] = eigArray[ind];
	acc[ind] = betaCopy[size - 1] * 
	  template_blas_fabs(eigVectors[(ind * size) + size - 1]) / 0.9; // FIXME: WHY IS THERE A FACTOR 0.9 HERE ?!?
      }
      delete[] eigArray;
      delete[] alphaCopy;
      delete[] betaCopy;
      delete[] isuppz;
      delete[] work;
      delete[] iwork;
      
#else
      Treal* work = new Treal[5 * size];
      int* iwork = new int[5 * size];
      int* ifail = new int[size];
      Treal dummy = -1.0;
      int info;
      int nEigsFound = 0;
      /* 
	 Elias note 2007-07-02:
	 There have been some problems with stevx returning with info=0
	 at the same time as nEigsFound != uppInd - lowInd + 1.
	 This is due to a bug in LAPACK which has been reported and confirmed.
	 To avoid this problem Elias changed the code so that ALL eigenvectors
	 are computed and then the desired ones are selected from the
	 complete list. 
      */
#if 1
      /* Original version of code calling stevx to get only the 
	 desired eigenvalues/vectors. */
      int const lowIndNew(lowInd + 1);
      int const uppIndNew(uppInd + 1);
      mat::stevx("V", "I", &size, alphaCopy, betaCopy, 
		 &dummy, &dummy, &lowIndNew, &uppIndNew,
		 &abstol, 
		 &nEigsFound, eigArray, eigVectors, &size, 
		 work, iwork, ifail, 
		 &info);
      assert(info == 0);
      assert(nEigsFound == uppInd - lowInd + 1);
      for (int ind = 0; ind < nEigsFound; ind++) {
	eigVals[ind] = eigArray[ind];
	acc[ind] = betaCopy[size - 1] * 
	  template_blas_fabs(eigVectors[(ind * size) + size - 1]) / 0.9; // FIXME: WHY IS THERE A FACTOR 0.9 HERE ?!?
      }
#else
      /* Modified version of code calling stevx to get ALL 
	 eigenvalues/vectors, and then picking out the desired ones. */
      Treal* eigVectorsTmp = new Treal[size*size];
      int const lowIndNew(1);
      int const uppIndNew(size);
      mat::stevx("V", "A", &size, alphaCopy, betaCopy, 
		 &dummy, &dummy, &lowIndNew, &uppIndNew, 
		 &abstol, 
		 &nEigsFound, eigArray, eigVectorsTmp, &size, 
		 work, iwork, ifail, 
		 &info);
      assert(info == 0);
      assert(nEigsFound == size);
      int nEigsWanted = uppInd - lowInd + 1;
      /* Copy desired eigenvectors from eigVectorsTmp to eigVectors. */
      for(int i = 0; i < nEigsWanted; i++)
	for(int j = 0; j < size; j++)
	  eigVectors[i*size+j] = eigVectorsTmp[(lowInd+i)*size+j];
      delete [] eigVectorsTmp;
      for (int ind = 0; ind < nEigsWanted; ind++) {
	eigVals[ind] = eigArray[lowInd+ind];
	acc[ind] = betaCopy[size - 1] * 
	  template_blas_fabs(eigVectors[(ind * size) + size - 1]) / 0.9; // FIXME: WHY IS THERE A FACTOR 0.9 HERE ?!?
      }
#endif
      delete[] eigArray;
      delete[] alphaCopy;
      delete[] betaCopy;
      delete[] work;
      delete[] iwork;
      delete[] ifail;
#endif
    }



    template<typename Treal>
      void MatrixTridiagSymmetric<Treal>::
      increaseCapacity(int const newCapacity) {
      capacity = newCapacity;
      Treal* alphaNew = new Treal[capacity];
      Treal* betaNew  = new Treal[capacity];
      for (int ind = 0; ind < size; ind++){
	alphaNew[ind] = alphaVec[ind];
	betaNew[ind] = betaVec[ind];
      }
      delete[] alphaVec;
      delete[] betaVec;
      alphaVec = alphaNew;
      betaVec = betaNew;
    }


  } /* end namespace arn */
} /* end namespace mat */
#endif
