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

/** @file TC2.h Trace correcting purification class
 *
 * Copyright(c) Emanuel Rubensson 2006
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date May 2006
 *
 */
#ifndef MAT_TC2
#define MAT_TC2
#include <math.h>
#include "bisection.h"
namespace mat {
  /** Trace correcting purification.
   *  This template instantiates the trace correcting purification algorithm
   *  developed by Niklasson [Phys. Rev. B 66, 155115 (2002)] with 
   *  modifications by Rubensson and Rudberg [unpublished].
   *  The template can be used with any matrix class Tmatrix that has the 
   *  following member functions:
   *  - gershgorin(Treal&, Treal&) const
   *  - add_identity(Treal)
   *  - operator*=(Treal)
   *  - operator=(Tmatrix const &)
   *  - trace() const
   *  - frob_thresh(Treal)
   *
   *  The matrix class should also support the following syntax:
   *  - A = alpha * B * B 
   *  - A = alpha * B * B + beta * A 
   *
   *  where A and B are of type Tmatrix and alpha and beta are of type Treal.
   * 
   */  
  template<typename Treal, typename Tmatrix>
    class TC2 {
  public:    
    TC2(Tmatrix& F, /**< Fock/Kohn-Sham matrix (input/workspace) */ 
	Tmatrix& DM, /**< Density matrix (output) */
	const int size, /**< System size (Number of basis functions)*/
	const int noc, /**< Number of occupied orbitals. */
	const Treal trunc = 0,/**< Threshold for truncation in Frobenius norm.
			       */
	const int maxmm = 100 /**< Maximum aloud number of mm-multiplications.
			       */
	); 
    /**<  Constructor 
     * Initializes everything.
     */
    Treal fermi_level(Treal tol = 1e-15 /**< Fault-tolerance for result. */
		      ) const; 
    /**< Returns the Fermi level.
     * Run after call to purify().
     */
    Treal homo(Treal tol = 1e-15 /**< Fault-tolerance for result. */
	       ) const;
    /**< Returns upper bound of the HOMO eigenvalue.
     * Run after call to purify().
     */
    Treal lumo(Treal tol = 1e-15 /**< Fault-tolerance for result. */
	       ) const;
    /**< Returns lower bound of the LUMO eigenvalue.
     * Run after call to purify().
     */

    inline int n_multiplies() const {
      return nmul;
    }
    /**< Returns the number of used matrix matrix multiplications */
    void print_data(int const start, int const stop) const;
    virtual ~TC2() {
      delete[] idemerror;
      delete[] tracediff;
      delete[] polys;
    } /**< Destructor. */
  protected:
    Tmatrix& X; /**<  Fock / Kohn-Sham matrix at initialization. 
		 *    Then used as workspace by purify().
		 *    Empty after call to purify().
		 */
    Tmatrix& D; /**< Density matrix after purification. */
    const int n; /**< System size. */
    const int nocc;   /**< Number of occupied orbitals. */
    const Treal frob_trunc; /**< Threshold for the truncation. */
    const int maxmul; /**< Number of tolerated matrix multiplications. */
    Treal lmin; /**< Lower bound for eigenvalue spectrum. */
    Treal lmax; /**< Upper bound for eigenvalue spectrum. */
    int nmul; /**< Number of used matrix multiplications. */
    int nmul_firstpart; /**< Number of used matrix multiplications in 
			 *   the first part of the purification. 
			 */
    Treal* idemerror; /**< Upper bound of euclidean norm ||D-D^2||_2 before 
		       *   each step. 
		       *   This means: idemerror[i] = norm(D[i]-D[i]^2)
		       *   where D[0] is the initial matrix and D[i] is the
		       *   matrix after i steps in the purification.
		       *   This value is calculated after the step since
		       *   D[i]^2 or 2D[i] - D[i]^2 is needed.
		       *   Length: nmul 
		       */
    Treal* tracediff; /**< The difference between the trace of the matrix and 
		       *   the number of occupied orbitals before each step.
		       *	 Length: nmul + 1     
		       */
    int* polys; /**< Choices of polynomials 0 for x^2 and 1 for 2x-x^2 
		 *   Length: nmul
		 */ 
    void purify(); /**< Runs purification. 
		    * Run by constructor.
		    */
  private:
    class Fun;
  };
  /** Help class for bisection root finding calls.
   * @see fermi_level
   * @see homo
   * @see lumo
   */
  template<typename Treal, typename Tmatrix>
    class TC2<Treal, Tmatrix>::Fun {
  public:
    Fun(int const* const p, int const pl, Treal const s)
      :pol(p), pollength(pl), shift(s) {
      assert(shift <= 1 && shift >= 0);
      assert(pollength >= 0);
    }
      Treal eval(Treal const x) const {
	Treal y = x;
	for (int ind = 0; ind < pollength; ind++ )
	  y = 2 * pol[ind] * y + (1 - 2 * pol[ind]) * y * y;
	/*
	 * pol[ind] == 0 --> y = y * y
	 * pol[ind] == 1 --> y = 2 * y - y * y
	 */
	return y - shift;
      }
  protected:
  private:
      int const* const pol;
      int const pollength;
      Treal const shift;
  };
  

  template<typename Treal, typename Tmatrix>
    TC2<Treal, Tmatrix>::TC2(Tmatrix& F, Tmatrix& DM, const int size, 
			     const int noc, 
			     const Treal trunc, const int maxmm)
    :X(F), D(DM), n(size), nocc(noc), frob_trunc(trunc), maxmul(maxmm), 
    lmin(0), lmax(0), nmul(0), nmul_firstpart(0), 
    idemerror(0), tracediff(0), polys(0) {
    assert(frob_trunc >= 0);
    assert(nocc >= 0);
    assert(maxmul >= 0);
    X.gershgorin(lmin, lmax);    /* Find eigenvalue bounds              */
    X.add_identity(-lmax);      /* Scale to [0, 1] interval and negate */
    X *= ((Treal)1.0 / (lmin - lmax));
    D = X;
    idemerror = new Treal[maxmul];
    tracediff = new Treal[maxmul + 1];
    polys = new int[maxmul];
    tracediff[0] = X.trace() - nocc;
    purify(); /**< Run purification */
  } /**< Constructor */


  template<typename Treal, typename Tmatrix>
    void TC2<Treal, Tmatrix>::purify() {
    assert(nmul == 0);
    assert(nmul_firstpart == 0);
    Treal delta, beta, trD2;
    int ind;
    Treal const ONE = 1;
    Treal const TWO = 2;
    do {
      if (nmul >= maxmul) {
	print_data(0, nmul);
	throw AcceptableMaxIter("TC2<Treal, Tmatrix>::purify(): "
				"Purification reached maxmul"
				" without convergence", maxmul);
      }
      if (tracediff[nmul] > 0) {
	D = ONE * X * X;
	polys[nmul] = 0;
      }
      else {
	D = -ONE * X * X + TWO * D;
	polys[nmul] = 1;
      }
      D.frob_thresh(frob_trunc);
      idemerror[nmul] = Tmatrix::frob_diff(D, X);
      ++nmul;
      tracediff[nmul] = D.trace() - nocc;
      X = D;
      /* Setting up convergence criteria */
      beta = (3 - template_blas_sqrt(5)) / 2 - frob_trunc;
      if (idemerror[nmul - 1] < 1 / (Treal)4 && 
	  (1 - template_blas_sqrt(1 - 4 * idemerror[nmul - 1])) / 2 < beta)
	beta = (1 + template_blas_sqrt(1 - 4 * idemerror[nmul - 1])) / 2;
      trD2  = (tracediff[nmul] + nocc - 
	       2 * polys[nmul - 1] * (tracediff[nmul - 1] + nocc)) / 
	(1 - 2 * polys[nmul - 1]);
      delta = frob_trunc;
      ind = nmul - 1;
      while (ind > 0 && polys[ind] == polys[ind - 1]) {
	delta = delta + frob_trunc;
	ind--;
      }
      delta = delta < (template_blas_sqrt(1 + 4 * idemerror[nmul - 1]) - 1) / 2 ?
	delta : (template_blas_sqrt(1 + 4 * idemerror[nmul - 1]) - 1) / 2;
    } while((trD2 + beta * (1 + delta) * n - (1 + delta + beta) * 
	     (tracediff[nmul - 1] + nocc)) / 
	    ((1 + 2 * delta) * (delta + beta)) < n - nocc - 1 ||
	    (trD2 - delta * (1 - beta) * n - (1 - delta - beta) * 
	     (tracediff[nmul - 1] + nocc)) / 
	    ((1 + 2 * delta) * (delta + beta)) < nocc - 1);
    
    /* Note that: */
    /* tracediff[i] - tracediff[i - 1] = trace(D[i]) - trace(D[i - 1]) */
    /* i.e. the change of the trace. */
    
    /* Take one step to make sure the eigenvalues stays in */
    /* { [ 0 , 2 * epsilon [ , ] 1 - 2 * epsilon , 1] }    */
    if (tracediff[nmul - 1] > 0) { 
      /* The same tracediff as in the last step is used since we want to */
      /* take a step with the other direction (with the other polynomial).*/
      D = -ONE * X * X + TWO * D; /* This is correct!! */
      polys[nmul] = 1;
    }
    else {
      D = ONE * X * X; /* This is correct!! */
      polys[nmul] = 0;
    }
    D.frob_thresh(frob_trunc);
    idemerror[nmul] = Tmatrix::frob_diff(D, X);
    ++nmul;
    tracediff[nmul] = D.trace() - nocc;

    nmul_firstpart = nmul; /* First part of purification finished. At this */
    /* point the eigenvalues are separated but have not yet converged.     */
    /* Use second order convergence polynomials to converge completely:    */
    do { 
      if (nmul + 1 >= maxmul) {
	print_data(0, nmul);
	throw AcceptableMaxIter("TC2<Treal, Tmatrix>::purify(): "
				"Purification reached maxmul"
				" without convergence", maxmul);
      }
      if (tracediff[nmul] > 0) {
	X = ONE * D * D;
	idemerror[nmul] = Tmatrix::frob_diff(D, X);
	D = X;
	polys[nmul] = 0;
	++nmul;
	tracediff[nmul] = D.trace() - nocc;
	
	D = -ONE * X * X + TWO * D;
	idemerror[nmul] = Tmatrix::frob_diff(D, X);
	polys[nmul] = 1;
	++nmul;
	tracediff[nmul] = D.trace() - nocc;
      }
      else {
	X = D;
	X = -ONE * D * D + TWO * X;
	idemerror[nmul] = Tmatrix::frob_diff(D, X);
	polys[nmul] = 1;
	++nmul;
	tracediff[nmul] = X.trace() - nocc;
	
	D = ONE * X * X;
	idemerror[nmul] = Tmatrix::frob_diff(D, X);
	polys[nmul] = 0;
	++nmul;
	tracediff[nmul] = D.trace() - nocc;
      }
      D.frob_thresh(frob_trunc);
#if 0
    } while (idemerror[nmul - 1] > frob_trunc); /* FIXME Check conv. crit. */
#else
  } while ((1 - template_blas_sqrt(1 - 4 * idemerror[nmul - 1])) / 2 > frob_trunc); 
#endif
  X.clear();
}

  template<typename Treal, typename Tmatrix>
    Treal TC2<Treal, Tmatrix>::fermi_level(Treal tol) const {
    Fun const fermifun(polys, nmul, 0.5);
    Treal chempot = bisection(fermifun, (Treal)0, (Treal)1, tol);
    return (lmin - lmax) * chempot + lmax;
  }

  template<typename Treal, typename Tmatrix>
    Treal TC2<Treal, Tmatrix>::homo(Treal tol) const {
    Treal homo = 0;
    Treal tmp;
    for (int mul = nmul_firstpart; mul < nmul; mul++) {
      if (idemerror[mul] < 1.0 / 4) {
	Fun const homofun(polys, mul, (1 + template_blas_sqrt(1 - 4 * idemerror[mul])) / 2);
	tmp = bisection(homofun, (Treal)0, (Treal)1, tol);
	/*
	  std::cout <<  tmp << " ,   ";
	  std::cout << (lmin - lmax) * tmp + lmax << std::endl;
	*/
	homo = tmp > homo ? tmp : homo;
      }
    }
    return  (lmin - lmax) * homo + lmax;
  }
  
  template<typename Treal, typename Tmatrix>
    Treal  TC2<Treal, Tmatrix>::lumo(Treal tol) const {
    Treal lumo = 1;
    Treal tmp;
    for (int mul = nmul_firstpart; mul < nmul; mul++) {
      if (idemerror[mul] < 1.0 / 4) {
	Fun const lumofun(polys, mul, (1 - template_blas_sqrt(1 - 4 * idemerror[mul])) / 2);
	tmp = bisection(lumofun, (Treal)0, (Treal)1, tol);
	/*
	  std::cout <<  tmp << " ,   ";
	  std::cout << (lmin - lmax) * tmp + lmax << std::endl;
	*/
	lumo = tmp < lumo ? tmp : lumo;
      }
    }
    return  (lmin - lmax) * lumo + lmax;
  }
   
template<typename Treal, typename Tmatrix>
  void TC2<Treal, Tmatrix>::print_data(int const start, int const stop) const {
  for (int ind = start; ind < stop; ind ++) {
    std::cout << "Iteration: " << ind
	      << "   Idempotency error: " << idemerror[ind]
	      << "   Tracediff: " << tracediff[ind]
	      << "   Poly: " << polys[ind] 
	      << std::endl;
  } 
}


} /* end namespace mat */
#endif
