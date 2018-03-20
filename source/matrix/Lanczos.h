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

/** @file Lanczos.h Class for building Krylov subspaces with the Lanczos method
 *
 * Copyright(c) Emanuel Rubensson 2006
 *
 * @author Emanuel Rubensson
 * @date December 2006
 *
 */
#ifndef MAT_LANCZOS
#define MAT_LANCZOS
#include "MatrixTridiagSymmetric.h"
#include "mat_utils.h"
namespace mat { /* Matrix namespace */
  namespace arn { /* Arnoldi type methods namespace */

    /** Class template for building Krylov subspaces with Lanczos
     *
     * Build up Krylov subspace for symmetric matrix with 
     * a Lanczos process.   
     * 
     * 
     * Treal: Type for real numbers
     *
     * Tmatrix: The matrix class
     *
     * Tvector: Vector class  
     * 
     */  
    template<typename Treal, typename Tmatrix, typename Tvector>
      class Lanczos {
    public:
      Lanczos(Tmatrix const & AA, Tvector const & startVec, 
	      int maxIt = 100, int cap = 100) 
	: A(AA), v(new Tvector[cap]), capacity(cap), j(0), maxIter(maxIt),    
	alpha(0), beta(0) {
	assert(cap > 1);
	Treal const ONE = 1.0;	  
	v[0] = startVec;
	if(v[0].eucl() < template_blas_sqrt(getRelPrecision<Treal>())) {
	  v[0].rand();
	}
	v[0] *= (ONE / v[0].eucl());
	r = v[0];
      }
	void restart(Tvector const & startVec) {
	  j = 0;
	  alpha = 0;
	  beta = 0;
	  delete[] v;
	  v = new Tvector[capacity];
	  Treal const ONE = 1.0;
	  v[0] = startVec;
	  v[0] *= (ONE / startVec.eucl());
	  r = startVec;	  
	  Tri.clear();
	}
	
	virtual void run() {
	  do {
	    step();
	    update();
	    if (j > maxIter)
	      throw AcceptableMaxIter("Lanczos::run() did not converge within maxIter");
	  }
	  while (!converged());
	}

	inline void copyTridiag(MatrixTridiagSymmetric<Treal> & Tricopy) {
	  Tricopy = Tri;
	}
	virtual ~Lanczos() {
	  delete[] v;
	}
    protected:
      Tmatrix const & A;
      Tvector* v; /** Vectors spanning Krylov subspace. 
		   *  In step j: Vectors 0 : j-2 is on file
		   *             Vectors j-1 : j is in memory
		   */

      Tvector r; /** Residual vector */
      MatrixTridiagSymmetric<Treal> Tri;
      int capacity;
      int j;     /** Current step */
      int maxIter;
      void increaseCapacity(int const newCapacity);
      void step();
      void getEigVector(Tvector& eigVec, Treal const * const eVecTri) const;
      virtual void update() = 0;
      virtual bool converged() const = 0;
    private:
      Treal alpha;
      Treal beta;
    }; /* end class definition Lanczos */

    template<typename Treal, typename Tmatrix, typename Tvector>
      void Lanczos<Treal, Tmatrix, Tvector>::step() {
      if (j + 1 >= capacity)
	increaseCapacity(capacity * 2);
      Treal const ONE = 1.0;
      A.matVecProd(r, v[j]);        // r = A * v[j] 
      alpha = transpose(v[j]) * r;
      r += (-alpha) * v[j];
      if (j) {
	r += (-beta) * v[j-1];
	v[j-1].writeToFile();
      }
      beta = r.eucl();
      v[j+1] = r;
      v[j+1] *= ONE / beta;
      Tri.increase(alpha, beta);
      ++j;
    }
    
    /*  FIXME: If several eigenvectors are needed it is more optimal to
     *  compute all of them at once since then the krylov subspace vectors
     *  only need to be read from memory once.
     */
    template<typename Treal, typename Tmatrix, typename Tvector>
      void Lanczos<Treal, Tmatrix, Tvector>::
      getEigVector(Tvector& eigVec, Treal const * const eVecTri) const {
      if (j <= 1) {
	eigVec = v[0];
      }	
      else {
	v[0].readFromFile();
	eigVec = v[0];
	v[0].writeToFile();
      }      
      eigVec *= eVecTri[0];
      for (int ind = 1; ind <= j - 2; ++ind) {
	v[ind].readFromFile();
     	eigVec += eVecTri[ind] * v[ind];
	v[ind].writeToFile();
      }
      eigVec += eVecTri[j-1] * v[j-1];
    }

    template<typename Treal, typename Tmatrix, typename Tvector>
      void Lanczos<Treal, Tmatrix, Tvector>::
      increaseCapacity(int const newCapacity) {
      assert(newCapacity > capacity);
      assert(j > 0);
      capacity = newCapacity;
      Tvector* new_v = new Tvector[capacity];
      /* FIXME: Fix so that file is copied when operator= is called in Vector
       * class
       */
      for (int ind = 0; ind <= j - 2; ind++){
	v[ind].readFromFile();
	new_v[ind] = v[ind];
	new_v[ind].writeToFile();
      }
      for (int ind = j - 1; ind <= j; ind++){
	new_v[ind] = v[ind];
      }
      delete[] v;
      v = new_v;
    }


  } /* end namespace arn */
} /* end namespace mat */
#endif
