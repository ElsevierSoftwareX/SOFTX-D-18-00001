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

/** @file Interval.h Interval class
 *
 * Copyright(c) Emanuel Rubensson 2007
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date January 2007
 *
 */
#ifndef MAT_INTERVAL
#define MAT_INTERVAL
#include <math.h>
#include <iomanip>
namespace mat {

  /* FIXME represent interval as midpoint and length to get better precision */
  template<typename Treal>
    class Interval {
  public:
    explicit Interval(Treal low = 1, Treal upp = -1)
      : lowerBound(low), upperBound(upp) {
    }
      inline bool empty() const {return lowerBound > upperBound;}
      
      static Interval intersect(Interval const & A, Interval const & B) {
	if (A.empty()) 
	  return A;
	else if (B.empty()) 
	  return B;
	else
	  return Interval(A.lowerBound > B.lowerBound ? 
			  A.lowerBound : B.lowerBound,
			  A.upperBound < B.upperBound ? 
			  A.upperBound : B.upperBound);
      }
      inline void intersect(Interval const & other) {
	if (this->empty()) 
	  return;
	else if (other.empty()) {
	  *this = other;
	  return;
	}
	else {
	  this->lowerBound = this->lowerBound > other.lowerBound ? 
	    this->lowerBound : other.lowerBound;
	  this->upperBound = this->upperBound < other.upperBound ? 
	    this->upperBound : other.upperBound;
	  return;
	}
      }

      inline void intersect_always_non_empty(Interval const & other) {
	if (this->empty()) { 
	  *this = other;
	  return;
	}
	if (other.empty()) 
	  return;

	Interval intersection = intersect(*this, other);
	if ( !intersection.empty() ) {
	  *this = intersection;
	  return;
	}
	// Ok, the intersection is actually empty
	Treal tmp_val;
	if ( this->lowerBound > other.upperBound ) 
	  tmp_val = (this->lowerBound + other.upperBound) / 2;
	else if ( this->upperBound < other.lowerBound ) 
	  tmp_val = (this->upperBound + other.lowerBound) / 2;
	else
	  assert(0); // This should never happen!
	this->lowerBound = tmp_val;
	this->upperBound = tmp_val;
	return;
      }

      /** Returns the length of the interval.
       *  0 if empty. 
       */
      inline Treal length() const {
	if (empty())
	  return 0.0;
	else
	  return upperBound - lowerBound;
      }
      inline Treal midPoint() const {
	assert(!empty());
	return (upperBound + lowerBound) / 2;
      }
      inline bool cover(Treal const value) const {
	if (empty())
	  return false;
	else
	  return (value <= upperBound && value >= lowerBound);
      }
      inline bool overlap(Interval const & other) const {
	Interval tmp = intersect(*this, other);
	return !tmp.empty();
      }

      /** Increases interval with value in both directions.
       *  Useful for error control.
       */
      inline void increase(Treal const value) {
	if (empty())
	  throw Failure("Interval<Treal>::increase(Treal const) : "
			"Attempt to increase empty interval.");
	lowerBound -= value;
	upperBound += value;
      }
      inline void decrease(Treal const value) {
	lowerBound += value;
	upperBound -= value;
      }
      inline Treal low() const {return lowerBound;}
      inline Treal upp() const {return upperBound;}

#if 0
      inline Interval<Treal> operator*(Interval<Treal> const & other) const {
	assert(lowerBound >= 0);
	assert(other.lowerBound >= 0);
	return Interval<Treal>(lowerBound * other.lowerBound,
			       upperBound * other.upperBound);
      }
#endif

      inline Interval<Treal> operator*(Treal const & value) const {
	if (value < 0)
	  return Interval<Treal>(upperBound * value, lowerBound * value);
	else
	  return Interval<Treal>(lowerBound * value, upperBound * value);
      }

      /* Minus operator well defined? */
      inline Interval<Treal> operator-(Interval<Treal> const & other) const {
	return *this + (-1.0 * other);
	//	return Interval<Treal>(lowerBound - other.lowerBound,
	//		       upperBound - other.upperBound);
      }

      inline Interval<Treal> operator+(Interval<Treal> const & other) const {
	return Interval<Treal>(lowerBound + other.lowerBound,
			       upperBound + other.upperBound);
      }
      inline Interval<Treal> operator/(Treal const & value) const {
	if (value < 0)
	  return Interval<Treal>(upperBound / value, lowerBound / value);
	else
	  return Interval<Treal>(lowerBound / value, upperBound / value);
      }
      inline Interval<Treal> operator-(Treal const & value) const {
	return Interval<Treal>(lowerBound - value, upperBound - value);
      }
      inline Interval<Treal> operator+(Treal const & value) const {
	return Interval<Treal>(lowerBound + value, upperBound + value);
      }



      void puriStep(int poly);
      void invPuriStep(int poly);
      // The two functions above really not needed now, just special
      // case of functions below with alpha == 1
      void puriStep(int poly, Treal alpha);
      void invPuriStep(int poly, Treal alpha);
  protected:
      Treal lowerBound;
      Treal upperBound;
  private:
  };

#if 0
  /* Minus operator is not well defined for intervals */
  template<typename Treal>
    inline Interval<Treal> operator-(Treal const value, 
				     Interval<Treal> const & other) {
    return Interval<Treal>(value - other.lowerBound,
			   value - other.upperBound);
  }
  template<typename Treal>
    inline Interval<Treal> operator+(Treal const value, 
				     Interval<Treal> const & other) {
    return Interval<Treal>(value + other.lowerBound,
			   value + other.upperBound);
  }
#endif

#if 1
  template<typename Treal>
    Interval<Treal> sqrtInt(Interval<Treal> const & other) {
    if (other.empty())
      return other;
    else {
      assert(other.low() >= 0);
      return Interval<Treal>(template_blas_sqrt(other.low()), template_blas_sqrt(other.upp()));
    }
  }
#endif

  template<typename Treal>
    void Interval<Treal>::puriStep(int poly) {
    if (upperBound > 2.0 || lowerBound < -1.0)
      throw Failure("Interval<Treal>::puriStep(int) : It is assumed here "
		    "that the interval I is within [-1.0, 2.0]");
    Interval<Treal> oneInt(-1.0,1.0);
    Interval<Treal> zeroInt(0.0,2.0);
    /* Elias note 2010-09-12:
       Sometimes the polynomial 2*x-x*x makes a non-empty interval in [0,1] 
       become empty because of rounding errors. For some reason this problem 
       showed up when using Intel compiler version 11.1 but not when using 
       version 10.1. Many test cases failed on Kalkyl when using 
       the compiler "icpc (ICC) 11.1 20100414".
       Anyway, we know that the function f(x) = 2*x-x*x is growing 
       in the interval [0,1] so we know a non-empty interval inside [0,1] should 
       always stay non-empty. To avoid assertion failures in purification, 
       the code below now forces the interval to be non-empty if it became 
       empty due to rounding errors. */
    bool nonEmptyIntervalInZeroToOne = false;
    if(upperBound > lowerBound && lowerBound >= 0 && upperBound <= 1)
      nonEmptyIntervalInZeroToOne = true;
    if (poly) {
      /* 2*x - x*x */
      *this = Interval<Treal>::intersect(*this,oneInt);
      //      assert(upperBound <= 1.0);
      upperBound = 2 * upperBound - upperBound * upperBound;
      lowerBound = 2 * lowerBound - lowerBound * lowerBound;
      if(nonEmptyIntervalInZeroToOne && upperBound < lowerBound) {
        // Force interval to be non-empty.
        Treal midPoint = (lowerBound + upperBound) / 2;
        upperBound = lowerBound = midPoint;
      }
    }
    else { /* x*x */
      *this = Interval<Treal>::intersect(*this,zeroInt);
      //      assert(lowerBound >= 0.0);
      upperBound = upperBound * upperBound;
      lowerBound = lowerBound * lowerBound;
    }
  }

  template<typename Treal>
    void Interval<Treal>::invPuriStep(int poly) {
    if (upperBound > 2.0 || lowerBound < -1.0)
      throw Failure("Interval<Treal>::puriStep(int) : It is assumed here "
		    "that the interval I is within [-1.0, 1.0]");
    Interval<Treal> oneInt(-1.0,1.0);
    Interval<Treal> zeroInt(0.0,2.0);
    if (poly) {
      /* 1 - sqrt(1 - x) */
      *this = Interval<Treal>::intersect(*this,oneInt);
      //      assert(upperBound <= 1.0);
      upperBound = 1.0 - template_blas_sqrt(1.0 - upperBound);
      lowerBound = 1.0 - template_blas_sqrt(1.0 - lowerBound);
    }
    else { /* sqrt(x) */
      *this = Interval<Treal>::intersect(*this,zeroInt);
      //      assert(lowerBound >= 0.0);
      upperBound = template_blas_sqrt(upperBound);
      lowerBound = template_blas_sqrt(lowerBound);
    }
  }

#if 1
  template<typename Treal>
    void Interval<Treal>::puriStep(int poly, Treal alpha) {
    if (upperBound > 2.0 || lowerBound < -1.0)
      throw Failure("Interval<Treal>::puriStep(int, real) : It is assumed here "
		    "that the interval I is within [-1.0, 2.0]");
    Interval<Treal> oneInt(-1.0,1.0);
    Interval<Treal> zeroInt(0.0,2.0);
    /* Elias note 2010-09-12:
       Sometimes the polynomial 2*x-x*x makes a non-empty interval in [0,1] 
       become empty because of rounding errors. For some reason this problem 
       showed up when using Intel compiler version 11.1 but not when using 
       version 10.1. Many test cases failed on Kalkyl when using 
       the compiler "icpc (ICC) 11.1 20100414".
       Anyway, we know that the function f(x) = 2*x-x*x is growing 
       in the interval [0,1] so we know a non-empty interval inside [0,1] should 
       always stay non-empty. To avoid assertion failures in purification, 
       the code below now forces the interval to be non-empty if it became 
       empty due to rounding errors. */
    bool nonEmptyIntervalInZeroToOne = false;
    if(upperBound > lowerBound && lowerBound >= 0 && upperBound <= 1)
      nonEmptyIntervalInZeroToOne = true;
    if (poly) {
      /* 2*alpha*x - alpha*alpha*x*x */
      *this = Interval<Treal>::intersect(*this,oneInt);
      //      assert(upperBound <= 1.0);
      upperBound = 2 * alpha * upperBound - alpha * alpha * upperBound * upperBound;
      lowerBound = 2 * alpha * lowerBound - alpha * alpha * lowerBound * lowerBound;
      if(nonEmptyIntervalInZeroToOne && upperBound < lowerBound) {
        // Force interval to be non-empty.
        Treal midPoint = (lowerBound + upperBound) / 2;
        upperBound = lowerBound = midPoint;
      }
    }
    else { /* ( alpha*x + (1-alpha) )^2 */
      *this = Interval<Treal>::intersect(*this,zeroInt);
      //      assert(lowerBound >= 0.0);
      upperBound = ( alpha * upperBound + (1-alpha) ) * ( alpha * upperBound + (1-alpha) );
      lowerBound = ( alpha * lowerBound + (1-alpha) ) * ( alpha * lowerBound + (1-alpha) );
    }
  }

  template<typename Treal>
    void Interval<Treal>::invPuriStep(int poly, Treal alpha) {
    if (upperBound > 2.0 || lowerBound < -1.0)
      throw Failure("Interval<Treal>::invPuriStep(int, real) : It is assumed here "
		    "that the interval I is within [-1.0, 2.0]");
    Interval<Treal> oneInt(-1.0,1.0);
    Interval<Treal> zeroInt(0.0,2.0);
    if (poly) {
      /* ( 1 - sqrt(1 - x) ) / alpha */
      *this = Interval<Treal>::intersect(*this,oneInt);
      //      assert(upperBound <= 1.0);
      upperBound = ( 1.0 - template_blas_sqrt(1.0 - upperBound) ) / alpha;
      lowerBound = ( 1.0 - template_blas_sqrt(1.0 - lowerBound) ) / alpha;
    }
    else { /* ( sqrt(x) - (1-alpha) ) / alpha */
      *this = Interval<Treal>::intersect(*this,zeroInt);
      //      assert(lowerBound >= 0.0);
      upperBound = ( template_blas_sqrt(upperBound) - (1-alpha) ) / alpha;
      lowerBound = ( template_blas_sqrt(lowerBound) - (1-alpha) ) / alpha;
    }
  }

#endif

  template<typename Treal>
    std::ostream& operator<<(std::ostream& s, 
			     Interval<Treal> const & in) {
    if (in.empty())
      s<<"[empty]";
    else
      s<<"["<<in.low()<<", "<<in.upp()<<"]";
    return s;
  }

} /* end namespace mat */
#endif
