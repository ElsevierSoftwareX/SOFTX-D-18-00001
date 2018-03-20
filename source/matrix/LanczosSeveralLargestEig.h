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

/** @file LanczosSeveralLargestEig.h Class for computing several largest 
 *  (note: not by magnitude) eigenvalues of a symmetric matrix with the Lanczos method.
 *
 * Copyright(c) Anastasia Kruchinina 2015
 *
 * @author Anastasia Kruchinina
 * @date   December 2015
 *
 */

#ifndef MAT_LANCZOSSEVERALLARGESTMAGNITUDEEIG
#define MAT_LANCZOSSEVERALLARGESTMAGNITUDEEIG

#include <limits>
#include <vector>

namespace mat { /* Matrix namespace */
  namespace arn { /* Arnoldi type methods namespace */

    template<typename Treal, typename Tmatrix, typename Tvector>
      class LanczosSeveralLargestEig 
      {
      public:
	// AA              - matrix
	// startVec        - starting guess vector
	// num_eigs        - number of eigenvalues to compute
	// maxIter(100)    - number of iterations
	// cap(100)        - estimated number of vectors in the Krylov subspace, will be increased if needed automatically
	// deflVec_(NULL)  - (deflation) vector corresponding to an uninteresting eigenvalue
	// sigma_(0)       - (deflation) shift of an uninteresting eigenvalue (to put it in the uninteresting part of the spectrum)
      LanczosSeveralLargestEig(Tmatrix const & AA, Tvector const & startVec, int num_eigs,
			       int maxit = 100, int cap = 100, Tvector * deflVec_ = NULL, Treal sigma_ = 0) 
	: A(AA),
	  v(new Tvector[cap]),
	  eigVectorTri(0),
	  capacity(cap),
	  j(0),
	  maxIter(maxit),
	  eValTmp(0),
	  accTmp(0),
	  number_of_eigenv(num_eigs),
	  alpha(0),
	  beta(0),
	  use_selective_orth(false),
	  use_full_orth(true),
	  counter_all(0),
	  counter_orth(0), 
	  deflVec(deflVec_)
	    {
	  assert(cap > 1);
	  Treal const ONE = 1.0;	  
	  v[0] = startVec;
	  if(v[0].eucl() < template_blas_sqrt(getRelPrecision<Treal>())) {
	    v[0].rand();
	  }
	  v[0] *= (ONE / v[0].eucl());
	  r = v[0];

	  if(number_of_eigenv == 1)
	    {unset_use_full_orth(); unset_use_selective_orth();}

    absTol = 1e-12;
    relTol = 1e-12; 
    sigma = sigma_;

	}

	// Absolute and relative tolerances
	// Absolute accuracy is measured by the residual ||Ax-lambda*x||
	// Realtive accuracy is measured by the relative residual ||Ax-lambda*x||/|lambda|
	void setRelTol(Treal const newTol) { relTol = newTol; }
	void setAbsTol(Treal const newTol) { absTol = newTol; }

	void set_use_selective_orth(){ use_selective_orth = true; }
	void set_use_full_orth(){ use_full_orth = true; }
	void unset_use_selective_orth(){ use_selective_orth = false; }
	void unset_use_full_orth(){ use_full_orth = false; }
	
	virtual void run() {
	  do {
	    if(j > 1 && use_selective_orth)
	      selective_orth();
	    step();
	    update();
	    if (j > maxIter)
	      throw AcceptableMaxIter("Lanczos::run() did not converge within maxIter");
	  }
	  while (!converged());
	}




	// i is a number of eigenvalue (1 is the largest, 2 is the second largest and so on)
	virtual void get_ith_eigenpair(int i, Treal& eigVal, Tvector& eigVec, Treal & acc)
	{
	  assert(i > 0);
	  assert(i <= size_accTmp);
	  eigVal = eValTmp[size_accTmp - i]; // array
	  assert(eigVectorTri);
	  getEigVector(eigVec, &eigVectorTri[j * (size_accTmp - i)]);
	  acc = accTmp[size_accTmp - i];
	}

	int get_num_iter() const{ return j;}
	
	virtual ~LanczosSeveralLargestEig() {

	  if(use_selective_orth)
	    printf("Orthogonalized %d of total possible %d, this is %lf %%\n", counter_orth, counter_all, (double)counter_orth/counter_all*100);

	  delete[] eigVectorTri;
	  delete[] eValTmp;
	  delete[] accTmp;
	  delete[] v;
	}


	inline void copyTridiag(MatrixTridiagSymmetric<Treal> & Tricopy) {
	  Tricopy = Tri;
	}


      protected:
	Tmatrix const & A;
	Tvector* v; /** Vectors spanning Krylov subspace. 
		     *  In step j: Vectors 0 : j-2 is on file
		     *             Vectors j-1 : j is in memory
		     */

	Tvector r; /** Residual vector */
	MatrixTridiagSymmetric<Treal> Tri;
	Treal* eigVectorTri; // Eigenvectors of the tridiagonal matrix 
	int capacity;
	int j;     /** Current step */
	int maxIter;
	void increaseCapacity(int const newCapacity);
	void getEigVector(Tvector& eigVec, Treal const * const eVecTri) const;
	Treal absTol;
	Treal relTol;
	virtual void step(); 
	virtual void computeEigenPairTri();
	virtual void update() {
	  computeEigenPairTri();
	}
	void selective_orth();
	virtual bool converged() const;
	virtual bool converged_ith(int i) const;
	Treal* eValTmp; // current computed eigenvalues (less or equal to number_of_eigenv) 
	Treal* accTmp;  // residuals
	int number_of_eigenv; // eigenvalues are saved in the decreasing order, thus the largest one has index 1
	int size_accTmp; // size of accTmp (number of computed eigenvalues of the matrix T)
      private:
	Treal alpha;
	Treal beta;

	bool use_selective_orth;
	bool use_full_orth;

	int counter_all;
	int counter_orth;

	// if deflation is used
	Tvector * deflVec;
	Treal sigma; 
      };

    template<typename Treal, typename Tmatrix, typename Tvector>
      void LanczosSeveralLargestEig<Treal, Tmatrix, Tvector>::
      selective_orth()
      {
	int j_curr = j-1;

	Treal coeff = 0, res;
	Treal normT = 0; // spectral norm of T (since norm of A is not available)
	// find largest by absolute value eigenvalue of T
	for(int i = 0; i <= j_curr; ++i)
	  if(template_blas_fabs(eValTmp[i]) > normT) normT = template_blas_fabs(eValTmp[i]);

	Treal epsilon = mat::getMachineEpsilon<Treal>();
	Tvector tmp;
	tmp = v[j_curr+1];
	tmp *= beta; // return non-normalized value

	for(int i = j_curr; i >= 0; --i)
	  {
	    counter_all++;
	    // get residual for this eigenpair
	    res = accTmp[i];
	    Treal tol = template_blas_sqrt(epsilon) * normT;
	     if(res <= tol) // b_{j} * |VT_i(j)| <= sqrt(eps) * norm(A), but we do not have norm(A)
	      {
		counter_orth++;
		Tvector eigVec;
		getEigVector(eigVec, &eigVectorTri[j_curr * i]);  // y = U*VT(:, i); % ith Ritz vector
		coeff = transpose(eigVec) * tmp;
		tmp += (-coeff) * (eigVec); // v = v - (y'*v)*y
	      }
	  }



	v[j_curr+1] = tmp;
	beta = v[j_curr+1].eucl(); // update beta
	Treal const ONE = 1.0;
	v[j_curr+1] *= ONE / beta; // normalized
	Tri.update_beta(beta);

	/* /\* // check orthogonality *\/ */
	/* try */
	/*   { */
	/*     for(int k = 0; k < j_curr; ++k) */
	/*       { */
	/* 	v[k].readFromFile(); */
	/* 	Treal val = transpose(v[k]) * v[j_curr+1]; */
	/* 	std::cout << val << ", "; */
	/* 	v[k].writeToFile(); */
	/*       } */
	/*     std::cout << std::endl; */
	/*   } */
	/* catch(const std::exception &e) */
	/*   { */
	/*     std::cout << "Exception: " << e.what() << std::endl; */
	/*   } */
      }




    template<typename Treal, typename Tmatrix, typename Tvector>
      void LanczosSeveralLargestEig<Treal, Tmatrix, Tvector>::
      step()
      {
	if (j + 1 >= capacity)
	  increaseCapacity(capacity * 2);
	Treal const ONE = 1.0;
	A.matVecProd(r, v[j]);        // r = A * v[j] 
	alpha = transpose(v[j]) * r;  // alpha = v[j]'*A*v[j] 

	/* 
	   If one wants to use deflation with vector
	   x_1:=deflVec (usually it is an eigenvector 
	   corresponding to an eigenvalue lambda_1 of A)
	   and thus compute eigenvalues of the matrix 
	   An = A-sigma*x_1*x_1'
	   Note: if lambga_i are eigenvalues of A corresponding to x_i, then
	   An will have eigenvalues (lambda_1-sigma, lambda_2, ..., lambda_N)
	   and unchanged eigenvectors x_i.
	 */

	if(deflVec != NULL)
	  {
	    /*
	      r = (A*vj - sigma*(x_1'*vj)*x_1) - alpha*vj - beta*v{j-1}
	      where 
	      alpha = vj'*An*vj = vj'*A*vj - sigma * (x_1'*vj)^2
	     */
	    Treal gamma = transpose(*deflVec) * v[j];  // dot product x' * v_j
	    alpha -= sigma*gamma*gamma;

	    r += (-sigma*gamma) * (*deflVec);
	  }

	r += (-alpha) * v[j];
	
	if (j) {
	  r += (-beta) * v[j-1];
	  v[j-1].writeToFile();
	}



	if(use_full_orth)
	  {
	    // full orthogonalization
	    // r = r - (q_1'*r)*q_1 - ... - (q_{j-1}'*r)*q_{j-1}
	    Treal gamma_i = 0;
	    Tvector tmp;
	    tmp = r;
	    for(int i = 0; i < j; ++i )
	      {
		v[i].readFromFile();
		gamma_i = transpose(tmp) * v[i]; // r'*v_i
		r += (-gamma_i) * v[i]; // (r'*vi) * v_i
		v[i].writeToFile();
	      }
	    tmp.clear();
	  }


	beta = r.eucl();
	v[j+1] = r;
	v[j+1] *= ONE / beta;
	Tri.increase(alpha, beta);


	/* /\* // check orthogonality *\/ */
	/* try */
	/*   { */
	/*     for(int k = 0; k < j; ++k) */
	/*       { */
	/* 	v[k].readFromFile(); */
	/* 	Treal val = transpose(v[k]) * v[j+1]; */
	/* 	std::cout << val << ", "; */
	/* 	v[k].writeToFile(); */
	/*       } */
	/*     std::cout << std::endl << "-----" << std::endl; */
	/*   } */
	/* catch(const std::exception &e) */
	/*   { */
	/*     std::cout << "Exception: " << e.what() << std::endl; */
	/*   } */


	j++;
      }


    /*
      Compute eigenvectors of the tridiagonal matrix
    */
    template<typename Treal, typename Tmatrix, typename Tvector>
      void LanczosSeveralLargestEig<Treal, Tmatrix, Tvector>::
      computeEigenPairTri() {
      if( eigVectorTri != NULL ) delete[] eigVectorTri;
      if( accTmp       != NULL ) delete[] accTmp;
      if( eValTmp      != NULL ) delete[] eValTmp;

      int num_compute_eigenvalues;
      if(use_selective_orth)      
	num_compute_eigenvalues = j; //  we need all eigenvectors of T
      else
	num_compute_eigenvalues = number_of_eigenv; // it is enough just number_of_eigenv of T

      /* Get largest eigenvalues */
      int const max_ind = j-1; // eigenvalue count starts with 0
      int const min_ind = std::max(j - num_compute_eigenvalues, 0);

      Treal* eigVectors = new Treal[j * num_compute_eigenvalues]; // every vector of size j
      Treal* eigVals = new Treal[num_compute_eigenvalues];
      Treal* accMax  = new Treal[num_compute_eigenvalues];
      assert(eigVectors != NULL);
      assert(eigVals != NULL);
      assert(accMax != NULL);

      Tri.getEigsByIndex(eigVals, eigVectors, accMax,  
			 min_ind,  max_ind);

      eValTmp = eigVals;


      eigVectorTri = eigVectors;
      accTmp = accMax;
      size_accTmp = num_compute_eigenvalues;

      // set unused pointers to NULL
      eigVectors = NULL;
      eigVals = NULL;
      accMax = NULL;
    }




    /*  FIXME: If several eigenvectors are needed it is more optimal to
     *  compute all of them at once since then the krylov subspace vectors
     *  only need to be read from memory once.
     */
    template<typename Treal, typename Tmatrix, typename Tvector>
      void LanczosSeveralLargestEig<Treal, Tmatrix, Tvector>::
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

      // normalized
      Treal norm_eigVec = eigVec.eucl(); 
      Treal const ONE = 1.0;
      eigVec *= ONE / norm_eigVec; 
    }


    // we want lowest eigenvalue to converge
    template<typename Treal, typename Tmatrix, typename Tvector>
      bool LanczosSeveralLargestEig<Treal, Tmatrix, Tvector>::
      converged() const {

      if(j < number_of_eigenv) return false;
      bool conv = converged_ith(number_of_eigenv); // if the last needed eigenvalue converged

      return conv;
    }

    // check convergence of ith eigenpair
    template<typename Treal, typename Tmatrix, typename Tvector>
      bool LanczosSeveralLargestEig<Treal, Tmatrix, Tvector>::
      converged_ith(int i) const {
      assert(size_accTmp >= i);

      bool conv = true;  //accTmp[size_accTmp - i] < absTol;                 /* Do not use absolute accuracy */ 
      if (template_blas_fabs(eValTmp[size_accTmp - i]) > 0) {
	conv = conv && 
	  accTmp[size_accTmp - i] / template_blas_fabs(eValTmp[size_accTmp - i]) < relTol; /* Relative acc.*/
      }
      return conv;
    }
    

    template<typename Treal, typename Tmatrix, typename Tvector>
      void LanczosSeveralLargestEig<Treal, Tmatrix, Tvector>::
      increaseCapacity(int const newCapacity) {
      assert(newCapacity > capacity);
      assert(j > 0);
      capacity = newCapacity;
      Tvector* new_v = new Tvector[capacity];
      assert(new_v != NULL);
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
