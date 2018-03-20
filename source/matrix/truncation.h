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

/** @file truncation.h Classes for truncation of small matrix elements.
 *
 * Copyright(c) Emanuel Rubensson 2010
 *
 * @author Emanuel Rubensson
 * @date April 2010
 * 
 * Most of this is essentially code that used to lie in
 * MatrixSymmetric.h somewhat rewritten for better structure and
 * reusability.
 *  
 *
 */
#ifndef MAT_TRUNCATION
#define MAT_TRUNCATION
#include <limits>
#include <stdexcept>
#include <cmath>
namespace mat { /* Matrix namespace */

  // Stuff for Euclidean norm based truncation
  
  template<typename Tmatrix, typename Treal>
    class EuclTruncationBase {
  public:
    explicit EuclTruncationBase( Tmatrix & A_ );
    Treal run(Treal const threshold);
    virtual ~EuclTruncationBase() {}
  protected:
    virtual void getFrobTruncBounds( Treal & lowTrunc, Treal & highTrunc, 
				     Treal const threshold ) = 0;
    virtual void getFrobSqNorms( std::vector<Treal> & frobsq_norms ) = 0;
    virtual void frobThreshLowLevel( Treal const threshold ) = 0;
    virtual Interval<Treal> euclIfSmall( Treal const absTol,
					 Treal const threshold ) = 0;
    Tmatrix & A; // Matrix to be truncated
    Tmatrix E; // Error matrix    
  };

  template<typename Tmatrix, typename Treal>
    EuclTruncationBase<Tmatrix, Treal>::EuclTruncationBase( Tmatrix & A_ ) 
    : A(A_) {
    SizesAndBlocks rows;
    SizesAndBlocks cols;
    A.getRows(rows);
    A.getCols(cols);
    E.resetSizesAndBlocks(rows, cols);      
  }

  template<typename Tmatrix, typename Treal>
    Treal EuclTruncationBase<Tmatrix, Treal>::run( Treal const threshold ) {
    assert(threshold >= (Treal)0.0);
    if (threshold == (Treal)0.0)
      return (Treal)0;
    std::vector<Treal> frobsq_norms;
    this->getFrobSqNorms( frobsq_norms ); /*=======*/
    std::sort(frobsq_norms.begin(), frobsq_norms.end());
    int low = -1;
    int high = frobsq_norms.size() - 1; 
    Treal lowFrobTrunc, highFrobTrunc;
    this->getFrobTruncBounds( lowFrobTrunc, highFrobTrunc, threshold ); /*=======*/
    Treal frobsqSum = 0;
    while( low < (int)frobsq_norms.size() - 1 && frobsqSum < lowFrobTrunc ) {
      ++low;
      frobsqSum += frobsq_norms[low];
    }
    high = low; /* Removing all tom high is to much */
    --low;
    while( high < (int)frobsq_norms.size() - 1 && frobsqSum < highFrobTrunc ) {
      ++high;
      frobsqSum += frobsq_norms[high];
    }
    // Now we have low and high
    int minStep   = int( 0.01 * frobsq_norms.size() ); // Consider elements in chunks of at least 1 percent of all elements at a time to not get too many iterations
    minStep = minStep > 0 ? minStep : 1; // step is at least one
    int testIndex = high;
    int previousTestIndex = high * 2;
    // Now, removing everything up to and including testIndex is too much
    Interval<Treal> euclEInt(0, threshold * 2);
    // We go from above (too many elements in the error matrix) and stop as soon as the error matrix is small enough
    while ( euclEInt.upp() > threshold ) {
      // Removing everything up to and including testIndex is too much, update high:
      high = testIndex;
      int stepSize = (int)((high - low) * 0.01); // We can accept that only 99% of elements possible to remove are removed
      // stepSize must be at least minStep:
      stepSize = stepSize >= minStep ? stepSize : minStep; 
      previousTestIndex = testIndex;
      testIndex -= stepSize;
      // testIndex cannot be smaller than low
      testIndex = testIndex > low ? testIndex : low; 
      /* Now we have decided the testIndex we would like to
	 use. However, we must be careful to handle the case when
	 there are several identical values in the frobsq_norms
	 list. In that case we use a modified value. */
      while(testIndex >= 0 && frobsq_norms[testIndex] == frobsq_norms[testIndex+1])
	testIndex--;
      /* Note that because of the above while loop, at this point it
	 is possible that testIndex < low. */
      if ( testIndex < 0 ) 
	// testIndex == -1, we have to break  
	break;
      assert( previousTestIndex != testIndex );
      Treal currentFrobTrunc = frobsq_norms[testIndex];
      frobThreshLowLevel( currentFrobTrunc ); /*=======*/
      euclEInt = euclIfSmall( Treal(threshold * 1e-2), threshold ); /*=======*/
      // Now we have an interval containing the Euclidean norm of E for the current testIndex
    } // end while
    Treal euclE; 
    if ( testIndex <= -1 ) {
      frobThreshLowLevel( (Treal)0.0 ); /*=======*/
      euclE = 0;
    }
    else {
      euclE = euclEInt.upp();
    }
    return euclE;    
  } // end run
  

  /** Truncation of symmetric matrices 
   *
   *
   */  
  template<typename Tmatrix, typename Treal>
    class EuclTruncationSymm : public EuclTruncationBase<Tmatrix, Treal> {
  public:
    explicit EuclTruncationSymm( Tmatrix & A_ ) 
      : EuclTruncationBase<Tmatrix, Treal>(A_) {}
  protected:
    virtual void getFrobTruncBounds( Treal & lowTrunc, Treal & highTrunc, 
				     Treal const threshold );
    virtual void getFrobSqNorms( std::vector<Treal> & frobsq_norms );
    virtual void frobThreshLowLevel( Treal const threshold );
    virtual Interval<Treal> euclIfSmall( Treal const absTol,
					 Treal const threshold );
  }; // end class EuclTruncationSymm
  
  template<typename Tmatrix, typename Treal>
    void EuclTruncationSymm<Tmatrix, Treal>::getFrobTruncBounds( Treal & lowTrunc, Treal & highTrunc, 
								 Treal const threshold ) {
    /* Divide by 2 because of symmetry */
    lowTrunc = (threshold * threshold) / 2;
    highTrunc = (threshold * threshold * this->A.get_nrows()) / 2;    
  }
  
  template<typename Tmatrix, typename Treal>
    void EuclTruncationSymm<Tmatrix, Treal>::getFrobSqNorms( std::vector<Treal> & frobsq_norms ) {
    this->A.getMatrix().getFrobSqLowestLevel(frobsq_norms);
  }
  
  template<typename Tmatrix, typename Treal>
    void EuclTruncationSymm<Tmatrix, Treal>::frobThreshLowLevel( Treal const threshold ) {
    this->A.getMatrix().frobThreshLowestLevel( threshold, &this->E.getMatrix() );
  }
  
  template<typename Tmatrix, typename Treal>
    Interval<Treal> EuclTruncationSymm<Tmatrix, Treal>::euclIfSmall( Treal const absTol,
								     Treal const threshold ) {
    Treal relTol = template_blas_sqrt(template_blas_sqrt(mat::getMachineEpsilon<Treal>()));
    Interval<Treal> tmpInterval =  mat::euclIfSmall(this->E, absTol, relTol, threshold); 
    if ( tmpInterval.length() < 2*absTol )
      return Interval<Treal>( tmpInterval.midPoint()-absTol, 
			      tmpInterval.midPoint()+absTol );
    return tmpInterval;
  }

  /** Truncation of symmetric matrices with Z 
   *
   * Truncation of a symmetric matrix A giving a truncated matrix B =
   * A + E such that the norm of the congruently transformed error
   * matrix ||Z^T * E * Z||_2 < threshold
   */  
  template<typename Tmatrix, typename TmatrixZ, typename Treal>
    class EuclTruncationSymmWithZ : public EuclTruncationSymm<Tmatrix, Treal> {
  public:
  EuclTruncationSymmWithZ( Tmatrix & A_, TmatrixZ const & Z_ ) 
    : EuclTruncationSymm<Tmatrix, Treal>(A_), Z(Z_) {}
  protected:
    virtual void getFrobTruncBounds( Treal & lowTrunc, Treal & highTrunc, 
				     Treal const threshold );
    // getFrobSqNorms(...)     from EuclTruncationSymm
    // frobThreshLowLevel(...) from EuclTruncationSymm
    virtual Interval<Treal> euclIfSmall( Treal const absTol,
					 Treal const threshold );
    TmatrixZ const & Z;
  }; // end class EuclTruncationSymmWithZ
  
  template<typename Tmatrix, typename TmatrixZ, typename Treal>
    void EuclTruncationSymmWithZ<Tmatrix, TmatrixZ, Treal>::getFrobTruncBounds( Treal & lowTrunc, Treal & highTrunc, 
										Treal const threshold ) {
    Treal Zfrob = Z.frob();
    Treal thresholdTakingZIntoAccount = threshold / (Zfrob * Zfrob);
    /* Divide by 2 because of symmetry */
    lowTrunc  = thresholdTakingZIntoAccount * thresholdTakingZIntoAccount / 2.0;
    highTrunc = template_blas_get_num_limit_max<Treal>();
  }
  
  template<typename Tmatrix, typename TmatrixZ, typename Treal>
    Interval<Treal> EuclTruncationSymmWithZ<Tmatrix, TmatrixZ, Treal>::euclIfSmall( Treal const absTol,
										    Treal const threshold ) {
    Treal relTol = template_blas_sqrt(template_blas_sqrt(mat::getMachineEpsilon<Treal>()));
    mat::TripleMatrix<Tmatrix, TmatrixZ, Treal> ErrMatTriple( this->E, Z);
    Interval<Treal> tmpInterval = mat::euclIfSmall(ErrMatTriple, absTol, relTol, threshold);      
    if ( tmpInterval.length() < 2*absTol )
      return Interval<Treal>( tmpInterval.midPoint()-absTol, 
			      tmpInterval.midPoint()+absTol );
    return tmpInterval;
  }
  
  /** Truncation of symmetric matrices at the element level (used for mixed norm truncation) 
   *
   * Works as EuclTruncationSymm but goes all the way to single matrix
   * elements. That is, it moves single matrix elements to and from
   * the error matrix.
   */
  template<typename Tmatrix, typename Treal>
    class EuclTruncationSymmElementLevel : public EuclTruncationSymm<Tmatrix, Treal> {
  public:
  explicit EuclTruncationSymmElementLevel( Tmatrix & A_ ) 
    : EuclTruncationSymm<Tmatrix, Treal>(A_) {}
  protected:
    // getFrobTruncBounds(...) from EuclTruncationSymm
    virtual void getFrobSqNorms( std::vector<Treal> & frobsq_norms );
    virtual void frobThreshLowLevel( Treal const threshold );
    // Interval<Treal> euclIfSmall(...) from EuclTruncationSymm    
  }; // end class EuclTruncationSymmElementLevel

  template<typename Tmatrix, typename Treal>
    void EuclTruncationSymmElementLevel<Tmatrix, Treal>::getFrobSqNorms( std::vector<Treal> & frobsq_norms ) {
    this->A.getMatrix().getFrobSqElementLevel(frobsq_norms);
  }

  template<typename Tmatrix, typename Treal>
    void EuclTruncationSymmElementLevel<Tmatrix, Treal>::frobThreshLowLevel( Treal const threshold ) {
    this->A.getMatrix().frobThreshElementLevel(threshold, &this->E.getMatrix() );
  }

  /** Truncation of general matrices 
   *
   *
   */  
  template<typename Tmatrix, typename Treal>
    class EuclTruncationGeneral : public EuclTruncationBase<Tmatrix, Treal> {
  public:
    explicit EuclTruncationGeneral( Tmatrix & A_ ) 
      : EuclTruncationBase<Tmatrix, Treal>(A_) {}
  protected:
    virtual void getFrobTruncBounds( Treal & lowTrunc, Treal & highTrunc, 
				     Treal const threshold );
    virtual void getFrobSqNorms( std::vector<Treal> & frobsq_norms );
    virtual void frobThreshLowLevel( Treal const threshold );
    virtual Interval<Treal> euclIfSmall( Treal const absTol,
					 Treal const threshold );
  }; // end class EuclTruncationGeneral
  
  template<typename Tmatrix, typename Treal>
    void EuclTruncationGeneral<Tmatrix, Treal>::getFrobTruncBounds( Treal & lowTrunc, Treal & highTrunc, 
								    Treal const threshold ) {
    // Try to improve bounds based on the Frobenius norm
    /*  ||E||_F^2 <= thres^2  -> 
     *  ||E||_F   <= thres    -> 
     *  ||E||_2   <= thresh 
     */
    lowTrunc = (threshold * threshold);
    /*  ||E||_F^2 >= thres^2 * n     -> 
     *  ||E||_F   >= thres * sqrt(n) -> 
     *  ||E||_2   >= thresh 
     */
    highTrunc = (threshold * threshold * this->A.get_nrows());    
  }
  
  template<typename Tmatrix, typename Treal>
    void EuclTruncationGeneral<Tmatrix, Treal>::getFrobSqNorms( std::vector<Treal> & frobsq_norms ) {
    this->A.getMatrix().getFrobSqLowestLevel(frobsq_norms);
  }
  
  template<typename Tmatrix, typename Treal>
    void EuclTruncationGeneral<Tmatrix, Treal>::frobThreshLowLevel( Treal const threshold ) {
    this->A.getMatrix().frobThreshLowestLevel( threshold, &this->E.getMatrix() );
  }
  
  template<typename Tmatrix, typename Treal>
    Interval<Treal> EuclTruncationGeneral<Tmatrix, Treal>::euclIfSmall( Treal const absTol,
									Treal const threshold ) {
    // FIXME: this should be changed (for all trunc classes) so that
    // some relative precision is always requested instead of the input
    // absTol which in the current case is not used(!)
    mat::ATAMatrix<Tmatrix, Treal> EtE(this->E);    
    Treal absTolDummy = template_blas_get_num_limit_max<Treal>(); // Treal(threshold * 1e-2)
    Treal relTol = 100 * mat::getMachineEpsilon<Treal>();
    Interval<Treal> tmpInterval = mat::euclIfSmall(EtE, absTolDummy, relTol, threshold);
    tmpInterval = Interval<Treal>( template_blas_sqrt(tmpInterval.low()), template_blas_sqrt(tmpInterval.upp()) );
    if ( tmpInterval.length() < 2*absTol )
      return Interval<Treal>( tmpInterval.midPoint()-absTol, 
			      tmpInterval.midPoint()+absTol );
    return tmpInterval;
  }




  /** Truncation of general matrices with impact on matrix triple multiply as error measure 
   *
   * Truncation of a matrix A giving a truncated matrix At =
   * A + E such that the norm of the congruently transformed error
   * matrix ||E^T * B * E + E^T * B * A + A^T * B * E||_2 < threshold
   */  
  template<typename Tmatrix, typename TmatrixB, typename Treal>
    class EuclTruncationCongrTransMeasure : public EuclTruncationGeneral<Tmatrix, Treal> {
  public:
  EuclTruncationCongrTransMeasure( Tmatrix & A_, TmatrixB const & B_ ) 
    : EuclTruncationGeneral<Tmatrix, Treal>(A_), B(B_) {}
  protected:
    virtual void getFrobTruncBounds( Treal & lowTrunc, Treal & highTrunc, 
				     Treal const threshold );
    // getFrobSqNorms(...)     from EuclTruncationGeneral
    // frobThreshLowLevel(...) from EuclTruncationGeneral
    virtual Interval<Treal> euclIfSmall( Treal const absTol,
					 Treal const threshold );
    TmatrixB const & B;
  }; // end class EuclTruncationCongrTransMeasure
  
  template<typename Tmatrix, typename TmatrixB, typename Treal>
    void EuclTruncationCongrTransMeasure<Tmatrix, TmatrixB, Treal>::getFrobTruncBounds( Treal & lowTrunc, 
											Treal & highTrunc, 
											Treal const threshold ) {
    Treal Afrob = this->A.frob();
    Treal Bfrob = B.frob();
    Treal tmp = -Afrob + template_blas_sqrt( Afrob*Afrob + threshold / Bfrob );
    lowTrunc = tmp*tmp; 
    highTrunc   = template_blas_get_num_limit_max<Treal>();
  }
  
  template<typename Tmatrix, typename TmatrixB, typename Treal>
    Interval<Treal> EuclTruncationCongrTransMeasure<Tmatrix, TmatrixB, Treal>::euclIfSmall( Treal const absTol,
											    Treal const threshold ) {
    Treal relTol = template_blas_sqrt(template_blas_sqrt(mat::getMachineEpsilon<Treal>()));
    mat::CongrTransErrorMatrix<TmatrixB, Tmatrix, Treal> ErrMatTriple( B, this->A, this->E );
    Interval<Treal> tmpInterval = mat::euclIfSmall(ErrMatTriple, absTol, relTol, threshold);      
    if ( tmpInterval.length() < 2*absTol ) {
      return Interval<Treal>( tmpInterval.midPoint()-absTol, 
			      tmpInterval.midPoint()+absTol );
    }
    return tmpInterval;
  }


} /* end namespace mat */
#endif
