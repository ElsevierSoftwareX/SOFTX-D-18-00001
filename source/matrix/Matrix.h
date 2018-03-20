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

/** @file Matrix.h The heart of the matrix library
 *
 * Copyright(c) Emanuel Rubensson 2005
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date September 2005
 *
 */


/** @section Matrix design principles

Heart of matrix library:
MatrixHierarchicBase is the base class both for higher levels in the 
hierarchic Matrix (first part of this file) and for the specialization 
at the lowest level (second part of this file).
Hence, functionality common to these class templates can be implemented in 
MatrixHierarchicBase.

The interface:
The interface is in MatrixBase.h MatrixSymmetric.h MatrixGeneral.h
MatrixTriangular.h.

Functions added to this file will not be visible. They have to be called 
from the interface classes. 

Example usage: Check out API_test.cc
*/


#ifndef MAT_MATRIX
#define MAT_MATRIX

#include <math.h>
#include <cstdlib>
#include <algorithm>

#include "MatrixHierarchicBase.h"
#include "matrix_proxy.h"
#include "mat_gblas.h"
#include "sort.h"
#include "allocate.h"
//#define MAT_NAIVE

namespace mat{
  enum side {left, right};
  enum inchversion {unstable, stable};
  template<class Treal, class Telement>
    class Vector;
  template<typename Tperm>
    struct AccessMap;

  class SingletonForTimings {
  private:
    double accumulatedTime;
  public:
    void reset() { accumulatedTime = 0; }
    double getAccumulatedTime() { return accumulatedTime; }
    void addTime(double timeToAdd) { 
#ifdef _OPENMP
      #pragma omp critical
#endif
      {
	accumulatedTime += timeToAdd;
      }
    }
    static SingletonForTimings & instance() {
      static SingletonForTimings theInstance;
      return theInstance;
    }
  private:
    SingletonForTimings(const SingletonForTimings & other);
    SingletonForTimings() : accumulatedTime(0) { }
  };


  /** Matrix class and heart of the matrix library
   *
   * This class is used to obtain the hierarchic data structure.
   *  
   * @see MatrixHierarchicBase
   * @see Permutation
   *
   */
  template<class Treal, class Telement = Treal>
    class Matrix: public MatrixHierarchicBase<Treal, Telement> {
    public:
    typedef Telement ElementType;
    typedef Vector<Treal, typename ElementType::VectorType> VectorType;

    friend class Vector<Treal, Telement>;
    Matrix():MatrixHierarchicBase<Treal, Telement>(){}
    

    void allocate() {
      assert(!this->is_empty());
      assert(this->is_zero());
      //#define MAT_USE_ALLOC_TIMER
#ifdef MAT_USE_ALLOC_TIMER
      mat::Time theTimer; theTimer.tic();
#endif
      this->elements = allocateElements<Telement>(this->nElements());
#ifdef MAT_USE_ALLOC_TIMER
      SingletonForTimings::instance().addTime(theTimer.toc());
#endif
      SizesAndBlocks colSAB;
      SizesAndBlocks rowSAB;
      for (int col = 0; col < this->cols.getNBlocks(); col++) {
	colSAB = this->cols.getSizesAndBlocksForLowerLevel(col);
	for (int row = 0; row < this->rows.getNBlocks(); row++) {
	  /* This could be improved - precompute rowSAB as well as colSAB */
	  rowSAB = this->rows.getSizesAndBlocksForLowerLevel(row);
	  (*this)(row,col).resetRows(rowSAB);
	  (*this)(row,col).resetCols(colSAB);
	}
      }
    }

    /* Full matrix assigns etc */
    void assignFromFull(std::vector<Treal> const & fullMat);
    void fullMatrix(std::vector<Treal> & fullMat) const; 
    void syFullMatrix(std::vector<Treal> & fullMat) const;
    void syUpTriFullMatrix(std::vector<Treal> & fullMat) const;
    
    /* Sparse matrix assigns etc */    
    void assignFromSparse(std::vector<int> const & rowind, 
			  std::vector<int> const & colind, 
			  std::vector<Treal> const & values);
    void assignFromSparse(std::vector<int> const & rowind, 
			  std::vector<int> const & colind, 
			  std::vector<Treal> const & values,
			  std::vector<int> const & indexes);
    /* Adds values (+=) to elements */
    void addValues(std::vector<int> const & rowind, 
		   std::vector<int> const & colind, 
		   std::vector<Treal> const & values);
    void addValues(std::vector<int> const & rowind, 
		   std::vector<int> const & colind, 
		   std::vector<Treal> const & values,
		   std::vector<int> const & indexes);

    void syAssignFromSparse(std::vector<int> const & rowind, 
			    std::vector<int> const & colind, 
			    std::vector<Treal> const & values);
    
    void syAddValues(std::vector<int> const & rowind, 
		     std::vector<int> const & colind, 
		     std::vector<Treal> const & values);
    
    void getValues(std::vector<int> const & rowind, 
		   std::vector<int> const & colind, 
		   std::vector<Treal> & values) const;
    void getValues(std::vector<int> const & rowind, 
		   std::vector<int> const & colind, 
		   std::vector<Treal> &,
		   std::vector<int> const & indexes) const;
    void syGetValues(std::vector<int> const & rowind, 
		     std::vector<int> const & colind, 
		     std::vector<Treal> & values) const;
    void getAllValues(std::vector<int> & rowind, 
		      std::vector<int> & colind, 
		      std::vector<Treal> &) const;
    void syGetAllValues(std::vector<int> & rowind, 
			std::vector<int> & colind, 
			std::vector<Treal> &) const;
    

    Matrix<Treal, Telement>& 
    operator=(const Matrix<Treal, Telement>& mat) {
      MatrixHierarchicBase<Treal, Telement>::operator=(mat);
      return *this;
    } 
    

    void clear(); 
    ~Matrix() {
      clear();
    }
    void writeToFile(std::ofstream & file) const;
    void readFromFile(std::ifstream & file);
    
    void random();
    void syRandom();
    /** Get a random zero structure with a specified probability that
     *  each submatrix is zero.
     */
    void randomZeroStructure(Treal probabilityBeingZero);
    void syRandomZeroStructure(Treal probabilityBeingZero);
    
    template<typename TRule>
    void setElementsByRule(TRule & rule);
    template<typename TRule>
    void sySetElementsByRule(TRule & rule);
    template<typename TRule>
    void trSetElementsByRule(TRule & rule) {
      // Setting elements for triangular is the same as for symmetric
      sySetElementsByRule(rule);
    }
    
    void addIdentity(Treal alpha);                  /* C += alpha * I     */

    static void transpose(Matrix<Treal, Telement> const & A,
			  Matrix<Treal, Telement> & AT);

    void symToNosym();
    void nosymToSym();

    /* Basic linear algebra routines */

    /* Set matrix to zero (k = 0) or identity (k = 1)                     */
    Matrix<Treal, Telement>& operator=(int const k);
    
    Matrix<Treal, Telement>& operator*=(const Treal alpha);

    static void gemm(const bool tA, const bool tB, const Treal alpha, 
		     const Matrix<Treal, Telement>& A, 
		     const Matrix<Treal, Telement>& B, 
		     const Treal beta, 
		     Matrix<Treal, Telement>& C);
    static void symm(const char side, const char uplo, const Treal alpha, 
		     const Matrix<Treal, Telement>& A, 
		     const Matrix<Treal, Telement>& B, 
		     const Treal beta, 
		     Matrix<Treal, Telement>& C);
    static void syrk(const char uplo, const bool tA, const Treal alpha, 
		     const Matrix<Treal, Telement>& A, 
		     const Treal beta, 
		     Matrix<Treal, Telement>& C);
    /* C = alpha * A * A + beta * C  where A and C are symmetric */
    static void sysq(const char uplo, const Treal alpha, 
		     const Matrix<Treal, Telement>& A, 
		     const Treal beta, 
		     Matrix<Treal, Telement>& C);
    /* C = alpha * A * B + beta * C where A and B are symmetric */
    static void ssmm(const Treal alpha, 
		     const Matrix<Treal, Telement>& A, 
		     const Matrix<Treal, Telement>& B, 
		     const Treal beta, 
		     Matrix<Treal, Telement>& C);
    /* C = alpha * A * B + beta * C where A and B are symmetric
     * and only the upper triangle of C is computed.
     */
    static void ssmm_upper_tr_only(const Treal alpha, 
				   const Matrix<Treal, Telement>& A, 
				   const Matrix<Treal, Telement>& B, 
				   const Treal beta, 
				   Matrix<Treal, Telement>& C);

    static void trmm(const char side, const char uplo, const bool tA, 
		     const Treal alpha, 
		     const Matrix<Treal, Telement>& A, 
		     Matrix<Treal, Telement>& B);

    /* Frobenius norms */

    /* Returns the Frobenius norm of the matrix.    */
    inline Treal frob() const { 
      return template_blas_sqrt(this->frobSquared());  
    } 
    /* Returns the squared Frobenius norm */
    Treal frobSquared() const;     
    inline Treal syFrob() const {
      return template_blas_sqrt(this->syFrobSquared());
    }
    Treal syFrobSquared() const;
    
    inline static Treal frobDiff
    (const Matrix<Treal, Telement>& A,
     const Matrix<Treal, Telement>& B) {
      return template_blas_sqrt(frobSquaredDiff(A, B));
    }
    static Treal frobSquaredDiff
    (const Matrix<Treal, Telement>& A,
     const Matrix<Treal, Telement>& B);
    
    inline static Treal syFrobDiff
    (const Matrix<Treal, Telement>& A,
     const Matrix<Treal, Telement>& B) {
      return template_blas_sqrt(syFrobSquaredDiff(A, B));	
    }
    static Treal syFrobSquaredDiff
    (const Matrix<Treal, Telement>& A,
     const Matrix<Treal, Telement>& B);      

    /* Traces */
    Treal trace() const;
    static Treal trace_ab(const Matrix<Treal, Telement>& A,
			  const Matrix<Treal, Telement>& B); 
    static Treal trace_aTb(const Matrix<Treal, Telement>& A,
			   const Matrix<Treal, Telement>& B); 
    static Treal sy_trace_ab(const Matrix<Treal, Telement>& A,
			     const Matrix<Treal, Telement>& B); 
    
    static void add(const Treal alpha,   /* B += alpha * A */
		    const Matrix<Treal, Telement>& A, 
		    Matrix<Treal, Telement>& B);
    void assign(Treal const  alpha,   /* *this = alpha * A */
		Matrix<Treal, Telement> const & A);


    /********** Help functions for thresholding */
    //  int getnnzBlocksLowestLevel() const;
    void getFrobSqLowestLevel(std::vector<Treal> & frobsq) const;
    void frobThreshLowestLevel
    (Treal const threshold, Matrix<Treal, Telement> * ErrorMatrix);

    void getFrobSqElementLevel(std::vector<Treal> & frobsq) const;
    void frobThreshElementLevel
    (Treal const threshold, Matrix<Treal, Telement> * ErrorMatrix);
    

#if 0
    inline void frobThreshLowestLevel
    (Treal const threshold, 
     Matrix<Treal, Telement> * ErrorMatrix) {
      bool a,b;
      frobThreshLowestLevel(threshold, ErrorMatrix, a, b);
    }
#endif

    /** Build a matrix with single matrix elements at the lowest level
	containing the Frobenius norms of the submatrices of A. */
    void assignFrobNormsLowestLevel
    ( Matrix<Treal, Matrix<Treal, Telement> > const & A ); 
    /** Version of assignFrobNormsLowestLevelToMatrix for symmetric matrices. */
    void syAssignFrobNormsLowestLevel
    ( Matrix<Treal, Matrix<Treal, Telement> > const & A ); 

    /** Same as assignFrobNormsLowestLevel except that the Frobenius
	norms of the differences between submatrices of A and B are
	assigned.  */
    void assignDiffFrobNormsLowestLevel
    ( Matrix<Treal, Matrix<Treal, Telement> > const & A,
      Matrix<Treal, Matrix<Treal, Telement> > const & B ); 
    /** Same as syAssignFrobNormsLowestLevel except that the Frobenius
	norms of the differences between submatrices of A and B are
	assigned.  */
    void syAssignDiffFrobNormsLowestLevel
    ( Matrix<Treal, Matrix<Treal, Telement> > const & A,
      Matrix<Treal, Matrix<Treal, Telement> > const & B ); 

    /** Truncate matrix A according to the sparsity pattern of the
     *  this matrix (frobNormMat).*/
    void truncateAccordingToSparsityPattern( Matrix<Treal, Matrix<Treal, Telement> > & A ) const;    


      /********** End of help functions for thresholding */

    static void gemm_upper_tr_only(const bool tA, const bool tB, 
				   const Treal alpha, 
				   const Matrix<Treal, Telement>& A, 
				   const Matrix<Treal, Telement>& B, 
				   const Treal beta, 
				   Matrix<Treal, Telement>& C);
    static void sytr_upper_tr_only(char const side, const Treal alpha,
				   Matrix<Treal, Telement>& A,
				   const Matrix<Treal, Telement>& Z);
    static void trmm_upper_tr_only(const char side, const char uplo, 
				   const bool tA, const Treal alpha, 
				   const Matrix<Treal, Telement>& A, 
				   Matrix<Treal, Telement>& B);
    static void trsytriplemm(char const side, 
			     const Matrix<Treal, Telement>& Z,
			     Matrix<Treal, Telement>& A);

    inline Treal frob_thresh
    (Treal const threshold,
     Matrix<Treal, Telement> * ErrorMatrix = 0) { 
      return template_blas_sqrt(frob_squared_thresh(threshold * threshold, 
						    ErrorMatrix)); 
    }                         
    /**< Removes small elements so that the introduced error is smaller than 
     * the threshold in the Frobenius norm
     * Returns the Frobenius norm of the introduced error.   
     */
    
    Treal frob_squared_thresh
    (Treal const threshold,
     Matrix<Treal, Telement> * ErrorMatrix = 0); 
    /**< Removes small elements so that the introduced error is smaller 
     * than threshold in the squared Frobenius norm, returns squared 
     * frobenius norm of the introduced error added to ErrorMatrix
     */
    
    static void syInch(const Matrix<Treal, Telement>& A,
		       Matrix<Treal, Telement>& Z,
		       const Treal threshold = 0, 
		       const side looking = left,
		       const inchversion version = unstable);

    void gershgorin(Treal& lmin, Treal& lmax) const; /* Computes bounds for*/
    /* real part of eigenvalues. The matrix must be quadratic (of course) */
    void sy_gershgorin(Treal& lmin, Treal& lmax) const {
      Matrix<Treal, Telement> tmp = (*this);
      tmp.symToNosym();
      tmp.gershgorin(lmin, lmax);
      return;
    }

    void add_abs_col_sums(Treal* abscolsums) const; /* Adds the absolute   */
    /* column sums to the abscolsums array.                                */
    /* abscolsums(i) += sum(abs(C(:,i))) for all i       (C == *this)      */
    /* Used by e.g. gershgorin eigenvalue inclusion                         */
    void get_diagonal(Treal* diag) const; /*Copy diagonal to the diag array*/
    
    size_t memory_usage() const; /* Returns memory usage in bytes */

    /* Note: we use size_t instead of int for nnz and nvalues to avoid
       integer overflow. */
    size_t nnz() const;    /**< Returns number of nonzeros in matrix. */
    size_t sy_nnz() const; /**< Returns number of nonzeros in matrix
			    *   including lower triangle elements.
			    */
    
    inline size_t nvalues() const {
      return nnz();
    } /**< Returns number of stored values in matrix.
       *   Returns same number as nnz()
       */
    size_t sy_nvalues() const; /**< Returns number of stored values in matrix.
				*   Lower triangle is not included.
				*   Lower triangle in diagonal submatrices
				*   is not included as well.
				*   Different from sy_nnz().
				*/
    
    template<typename Top> 
    Treal syAccumulateWith(Top & op) {
      Treal res = 0;
      if (!this->is_zero()) {
	for (int col = 0; col < this->ncols(); col++) {
	  for (int row = 0; row < col; row++) 
	    res += 2 * (*this)(row, col).geAccumulateWith(op);
	  res += (*this)(col, col).syAccumulateWith(op);
	}
      }
      return res;
    }
    /** Accumulation algorithm for general matrices */
    template<typename Top>
    Treal geAccumulateWith(Top & op) {
      Treal res = 0;
      if (!this->is_zero()) {
	for (int col = 0; col < this->ncols(); col++)
	  for (int row = 0; row < this->nrows(); row++)
	    res += (*this)(row, col).geAccumulateWith(op);
      }
      return res;
    }
    
    static inline unsigned int level() {return Telement::level() + 1;}

    Treal maxAbsValue() const {
      if (this->is_zero())
	return 0;
      else {
	Treal maxAbsGlobal = 0;
	Treal maxAbsLocal  = 0;
	for (int ind = 0; ind < this->nElements(); ++ind) {
	  maxAbsLocal = this->elements[ind].maxAbsValue();
	  maxAbsGlobal = maxAbsGlobal > maxAbsLocal ? 
	    maxAbsGlobal : maxAbsLocal; 
	} /* end for */
	return maxAbsGlobal;
      } 
    }
    
      protected:
      private:
  }; /* end class */
  

#if 1
  /* Full matrix assigns etc */
  template<class Treal, class Telement> 
    void Matrix<Treal, Telement>::
    assignFromFull(std::vector<Treal> const & fullMat) {
    int nTotalRows = this->rows.getNTotalScalars();
    int nTotalCols = this->cols.getNTotalScalars();
    assert((int)fullMat.size() == nTotalRows * nTotalCols);
    if (this->is_zero())
      allocate();
    for (int col = 0; col < this->ncols(); col++) 
      for (int row = 0; row < this->nrows(); row++) 
	(*this)(row, col).assignFromFull(fullMat);
  }
  
  template<class Treal, class Telement> 
    void Matrix<Treal, Telement>::
    fullMatrix(std::vector<Treal> & fullMat) const {
    int nTotalRows = this->rows.getNTotalScalars();
    int nTotalCols = this->cols.getNTotalScalars();
    fullMat.resize(nTotalRows * nTotalCols);
    if (this->is_zero()) {
      int rowOffset = this->rows.getOffset();
      int colOffset = this->cols.getOffset();
      for (int col = 0; col < this->nScalarsCols(); col++)
	for (int row = 0; row < this->nScalarsRows(); row++) 
	  fullMat[(colOffset + col) * nTotalRows + (rowOffset + row)] = 0;
    } 
    else {
      for (int col = 0; col < this->ncols(); col++) 
	for (int row = 0; row < this->nrows(); row++) 
	  (*this)(row, col).fullMatrix(fullMat);
    }
  }
  
  template<class Treal, class Telement> 
    void Matrix<Treal, Telement>::
    syFullMatrix(std::vector<Treal> & fullMat) const {
    int nTotalRows = this->rows.getNTotalScalars();
    int nTotalCols = this->cols.getNTotalScalars();
    fullMat.resize(nTotalRows * nTotalCols);
    if (this->is_zero()) {
      int rowOffset = this->rows.getOffset();
      int colOffset = this->cols.getOffset();
      for (int col = 0; col < this->nScalarsCols(); col++)
	for (int row = 0; row <= col; row++) {
	  fullMat[(colOffset + col) * nTotalRows + (rowOffset + row)] = 0;
	  fullMat[(rowOffset + row) * nTotalRows + (colOffset + col)] = 0;
	}
    } 
    else {
      for (int col = 0; col < this->ncols(); col++) {
	for (int row = 0; row < col; row++) 
	  (*this)(row, col).syUpTriFullMatrix(fullMat); 
	(*this)(col, col).syFullMatrix(fullMat);	
      }
    }
  }
  
  template<class Treal, class Telement> 
    void Matrix<Treal, Telement>::
    syUpTriFullMatrix(std::vector<Treal> & fullMat) const {
    int nTotalRows = this->rows.getNTotalScalars();
    int nTotalCols = this->cols.getNTotalScalars();
    fullMat.resize(nTotalRows * nTotalCols);
    if (this->is_zero()) {
      int rowOffset = this->rows.getOffset();
      int colOffset = this->cols.getOffset();
      for (int col = 0; col < this->nScalarsCols(); col++)
	for (int row = 0; row <= this->nScalarsRows(); row++) {
	  fullMat[(colOffset + col) * nTotalRows + (rowOffset + row)] = 0;
	  fullMat[(rowOffset + row) * nTotalRows + (colOffset + col)] = 0;
	}
    } 
    else {
      for (int col = 0; col < this->ncols(); col++) 
	for (int row = 0; row < this->nrows(); row++) 
	  (*this)(row, col).syUpTriFullMatrix(fullMat);
    }
  }
  
#endif


  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::
    assignFromSparse(std::vector<int> const & rowind, 
		     std::vector<int> const & colind, 
		     std::vector<Treal> const & values) {
    assert(rowind.size() == colind.size() &&
	   rowind.size() == values.size());
    std::vector<int> indexes(values.size());
    for (unsigned int ind = 0; ind < values.size(); ++ind) 
      indexes[ind] = ind;
    assignFromSparse(rowind, colind, values, indexes);
  }
  
  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::
    assignFromSparse(std::vector<int> const & rowind, 
		     std::vector<int> const & colind, 
		     std::vector<Treal> const & values,
		     std::vector<int> const & indexes) {
    if (indexes.empty()) {
      this->clear();
      return;
    }
    if (this->is_zero())
      allocate();
    
    std::vector<std::vector<int> > columnBuckets (this->cols.getNBlocks());
    std::vector<std::vector<int> > rowBuckets (this->rows.getNBlocks());    
    int currentInd = 0;
    

    std::vector<int>::const_iterator it;
    for ( it = indexes.begin(); it < indexes.end(); it++ )
      columnBuckets[ this->cols.whichBlock( colind[*it] ) ].push_back (*it);
    
    /* Go through all column buckets. */
    for (int col = 0; col < this->ncols(); col++) {
      /* Go through current column bucket and distribute to row buckets. */
      while (!columnBuckets[col].empty()) {
	currentInd = columnBuckets[col].back();
	columnBuckets[col].pop_back();
	rowBuckets[ this->rows.whichBlock
		    ( rowind[currentInd] ) ].push_back (currentInd);
      }
      /* Make calls to lower level for every row bucket. */
      for (int row = 0; row < this->nrows(); row++) {
	(*this)(row,col).assignFromSparse(rowind, colind, values, rowBuckets[row]);
	rowBuckets[row].clear();
      } /* end row loop */
    } /* end column loop */
  }  

  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::
    addValues(std::vector<int> const & rowind, 
	      std::vector<int> const & colind, 
	      std::vector<Treal> const & values) {
    assert(rowind.size() == colind.size() &&
	   rowind.size() == values.size());
    std::vector<int> indexes(values.size());
    for (unsigned int ind = 0; ind < values.size(); ++ind) 
      indexes[ind] = ind;
    addValues(rowind, colind, values, indexes);
  }
  
  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::
    addValues(std::vector<int> const & rowind, 
	      std::vector<int> const & colind, 
	      std::vector<Treal> const & values,
	      std::vector<int> const & indexes) {
    if (indexes.empty())
      return;
    if (this->is_zero())
      allocate();
    
    std::vector<std::vector<int> > columnBuckets (this->cols.getNBlocks());
    std::vector<std::vector<int> > rowBuckets (this->rows.getNBlocks());    
    int currentInd = 0;
    
    std::vector<int>::const_iterator it;
    for ( it = indexes.begin(); it < indexes.end(); it++ )
      columnBuckets[ this->cols.whichBlock( colind[*it] ) ].push_back (*it);
    
    /* Go through all column buckets. */
    for (int col = 0; col < this->ncols(); col++) {
      /* Go through current column bucket and distribute to row buckets. */
      while (!columnBuckets[col].empty()) {
	currentInd = columnBuckets[col].back();
	columnBuckets[col].pop_back();
	rowBuckets[ this->rows.whichBlock
		    ( rowind[currentInd] ) ].push_back (currentInd);
      }
      /* Make calls to lower level for every row bucket. */
      for (int row = 0; row < this->nrows(); row++) {
	(*this)(row,col).addValues(rowind, colind, values, rowBuckets[row]);
	rowBuckets[row].clear();
      } /* end row loop */
    } /* end column loop */
  }

  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::
    syAssignFromSparse(std::vector<int> const & rowind, 
		       std::vector<int> const & colind, 
		       std::vector<Treal> const & values) {
    assert(rowind.size() == colind.size() &&
	   rowind.size() == values.size());
    bool upperTriangleOnly = true;
    for (unsigned int ind = 0; ind < values.size(); ++ind) {
      upperTriangleOnly = 
	upperTriangleOnly && colind[ind] >= rowind[ind];
    }
    if (!upperTriangleOnly)
      throw Failure("Matrix<Treal, Telement>::"
		    "syAddValues(std::vector<int> const &, "
		    "std::vector<int> const &, "
		    "std::vector<Treal> const &, int const): "
		    "Only upper triangle can contain elements when assigning "
		    "symmetric or triangular matrix ");
    assignFromSparse(rowind, colind, values);
  }
  
  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::
    syAddValues(std::vector<int> const & rowind, 
		std::vector<int> const & colind, 
		std::vector<Treal> const & values) {
    assert(rowind.size() == colind.size() &&
	   rowind.size() == values.size());
    bool upperTriangleOnly = true;
    for (unsigned int ind = 0; ind < values.size(); ++ind) {
      upperTriangleOnly = 
	upperTriangleOnly && colind[ind] >= rowind[ind];
    }
    if (!upperTriangleOnly)
      throw Failure("Matrix<Treal, Telement>::"
		    "syAddValues(std::vector<int> const &, "
		    "std::vector<int> const &, "
		    "std::vector<Treal> const &, int const): "
		    "Only upper triangle can contain elements when assigning "
		    "symmetric or triangular matrix ");
    addValues(rowind, colind, values);
  }

  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::
    getValues(std::vector<int> const & rowind, 
	      std::vector<int> const & colind, 
	      std::vector<Treal> & values) const {
    assert(rowind.size() == colind.size());
    values.resize(rowind.size(),0);
    std::vector<int> indexes(rowind.size());
    for (unsigned int ind = 0; ind < rowind.size(); ++ind) 
      indexes[ind] = ind;
    getValues(rowind, colind, values, indexes);
  }
  
  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::
    getValues(std::vector<int> const & rowind, 
	      std::vector<int> const & colind, 
	      std::vector<Treal> & values,
	      std::vector<int> const & indexes) const {
    assert(!this->is_empty());
    if (indexes.empty())
      return;
    std::vector<int>::const_iterator it;
    if (this->is_zero()) {
      for ( it = indexes.begin(); it < indexes.end(); it++ )
	values[*it] = 0;
      return;
    }
    
    std::vector<std::vector<int> > columnBuckets (this->cols.getNBlocks());
    std::vector<std::vector<int> > rowBuckets (this->rows.getNBlocks());    
    int currentInd = 0;
    
    for ( it = indexes.begin(); it < indexes.end(); it++ )
      columnBuckets[ this->cols.whichBlock( colind[*it] ) ].push_back (*it);
    
    /* Go through all column buckets. */
    for (int col = 0; col < this->ncols(); col++) {
      /* Go through current column bucket and distribute to row buckets. */
      while (!columnBuckets[col].empty()) {
	currentInd = columnBuckets[col].back();
	columnBuckets[col].pop_back();
	rowBuckets[ this->rows.whichBlock
		    ( rowind[currentInd] ) ].push_back (currentInd);
      }
      /* Make calls to lower level for every row bucket. */
      for (int row = 0; row < this->nrows(); row++) {
	(*this)(row,col).getValues(rowind, colind, values, rowBuckets[row]);
	rowBuckets[row].clear();
      } /* end row loop */
    } /* end column loop */
  }
  
  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::
    syGetValues(std::vector<int> const & rowind, 
		std::vector<int> const & colind, 
		std::vector<Treal> & values) const {
    assert(rowind.size() == colind.size());
    bool upperTriangleOnly = true;
    for (int unsigned ind = 0; ind < rowind.size(); ++ind) {
      upperTriangleOnly = 
	upperTriangleOnly && colind[ind] >= rowind[ind];
    }
    if (!upperTriangleOnly)
      throw Failure("Matrix<Treal, Telement>::"
		    "syGetValues(std::vector<int> const &, "
		    "std::vector<int> const &, "
		    "std::vector<Treal> const &, int const): "
		    "Only upper triangle when retrieving elements from "
		    "symmetric or triangular matrix ");
    getValues(rowind, colind, values);
  }
  

  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::
    getAllValues(std::vector<int> & rowind, 
		 std::vector<int> & colind, 
		 std::vector<Treal> & values) const {
    assert(rowind.size() == colind.size() &&
	   rowind.size() == values.size());
    if (!this->is_zero()) 
      for (int ind = 0; ind < this->nElements(); ++ind)
	this->elements[ind].getAllValues(rowind, colind, values);
  }

  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::
    syGetAllValues(std::vector<int> & rowind, 
		   std::vector<int> & colind, 
		   std::vector<Treal> & values) const {
    assert(rowind.size() == colind.size() &&
	   rowind.size() == values.size());
    if (!this->is_zero()) 
      for (int col = 0; col < this->ncols(); ++col) {
	for (int row = 0; row < col; ++row)
	  (*this)(row, col).getAllValues(rowind, colind, values);
	(*this)(col, col).syGetAllValues(rowind, colind, values);
      }
  }
  
  
  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::clear() {
    if (!this->is_zero())
      for (int i = 0; i < this->nElements(); i++) 
	this->elements[i].clear();
    freeElements(this->elements);
    this->elements = 0;
  }
  
  template<class Treal, class Telement>
    void Matrix<Treal, Telement>:: 
    writeToFile(std::ofstream & file) const {
    int const ZERO = 0;
    int const ONE  = 1;
    if (this->is_zero()) {
      char * tmp = (char*)&ZERO;
      file.write(tmp,sizeof(int));
    }
    else {
      char * tmp = (char*)&ONE;
      file.write(tmp,sizeof(int));
      for (int i = 0; i < this->nElements(); i++)
	this->elements[i].writeToFile(file);
    }
  }
  template<class Treal, class Telement>
    void Matrix<Treal, Telement>:: 
    readFromFile(std::ifstream & file) {
    int const ZERO = 0;
    int const ONE  = 1;
    char tmp[sizeof(int)];
    file.read(tmp, (std::ifstream::pos_type)sizeof(int));
    switch ((int)*tmp) {
    case ZERO:
      this->clear();
      break;
    case ONE:
      if (this->is_zero()) 
	allocate();
      for (int i = 0; i < this->nElements(); i++)
	this->elements[i].readFromFile(file);
      break;
    default:
      throw Failure("Matrix<Treal, Telement>::" 
		    "readFromFile(std::ifstream & file):"
		    "File corruption int value not 0 or 1");
    }
  }

  template<class Treal, class Telement> 
    void Matrix<Treal, Telement>::random() {
    if (this->is_zero()) 
      allocate();
    for (int ind = 0; ind < this->nElements(); ind++)
      this->elements[ind].random();    	
  }

  template<class Treal, class Telement> 
    void Matrix<Treal, Telement>::syRandom() {
    if (this->is_zero()) 
      allocate();
    /* Above diagonal */
    for (int col = 1; col < this->ncols(); col++)
      for (int row = 0; row < col; row++)
	(*this)(row, col).random();
    /* Diagonal */
    for (int rc = 0; rc < this->nrows(); rc++)
      (*this)(rc,rc).syRandom();
  }

  template<class Treal, class Telement> 
    void Matrix<Treal, Telement>::
    randomZeroStructure(Treal probabilityBeingZero) {
    if (!this->highestLevel() && 
	probabilityBeingZero > rand() / (Treal)RAND_MAX)
      this->clear();
    else {
      if (this->is_zero()) 
	allocate();
      for (int ind = 0; ind < this->nElements(); ind++)
	this->elements[ind].randomZeroStructure(probabilityBeingZero);    	
    }      
  }

  template<class Treal, class Telement> 
    void  Matrix<Treal, Telement>::
    syRandomZeroStructure(Treal probabilityBeingZero) {
    if (!this->highestLevel() && 
	probabilityBeingZero > rand() / (Treal)RAND_MAX)
      this->clear();
    else {
      if (this->is_zero()) 
	allocate();
      /* Above diagonal */
      for (int col = 1; col < this->ncols(); col++)
	for (int row = 0; row < col; row++)
	  (*this)(row, col).randomZeroStructure(probabilityBeingZero);
      /* Diagonal */
      for (int rc = 0; rc < this->nrows(); rc++)
	(*this)(rc,rc).syRandomZeroStructure(probabilityBeingZero);
    }
  }

  template<class Treal, class Telement> 
    template<typename TRule>
    void Matrix<Treal, Telement>::
    setElementsByRule(TRule & rule) {
    if (this->is_zero()) 
      allocate();
    for (int ind = 0; ind < this->nElements(); ind++)
      this->elements[ind].setElementsByRule(rule);    	
  }
  template<class Treal, class Telement> 
    template<typename TRule>
    void Matrix<Treal, Telement>::
    sySetElementsByRule(TRule & rule) {
    if (this->is_zero()) 
      allocate();
    /* Above diagonal */
    for (int col = 1; col < this->ncols(); col++)
      for (int row = 0; row < col; row++)
	(*this)(row, col).setElementsByRule(rule);
    /* Diagonal */
    for (int rc = 0; rc < this->nrows(); rc++)
      (*this)(rc,rc).sySetElementsByRule(rule);
  }
  

  template<class Treal, class Telement> 
    void Matrix<Treal, Telement>::
    addIdentity(Treal alpha) {
    if (this->is_empty())
      throw Failure("Matrix<Treal, Telement>::addIdentity(Treal): "
		    "Cannot add identity to empty matrix.");
    if (this->ncols() != this->nrows()) 
      throw Failure("Matrix<Treal, Telement>::addIdentity(Treal): "
		    "Matrix must be square to add identity");
    if (this->is_zero()) 
      allocate();
    for (int ind = 0; ind < this->ncols(); ind++)
      (*this)(ind,ind).addIdentity(alpha);
  }
  
  template<class Treal, class Telement> 
    void Matrix<Treal, Telement>::
    transpose(Matrix<Treal, Telement> const & A,
	      Matrix<Treal, Telement> & AT) {
    if (A.is_zero()) { /* Condition also matches empty matrices. */
      AT.rows = A.cols;
      AT.cols = A.rows;
      freeElements(AT.elements);
      AT.elements = 0;
      return;
    }
    if (AT.is_zero() || (AT.nElements() != A.nElements())) {
      freeElements(AT.elements);
      AT.elements = allocateElements<Telement>(A.nElements());
    }
    AT.cols = A.rows;
    AT.rows = A.cols;
    for (int row = 0; row < AT.nrows(); row++)
      for (int col = 0; col < AT.ncols(); col++)
	Telement::transpose(A(col,row), AT(row,col));
  }

  template<class Treal, class Telement> 
    void Matrix<Treal, Telement>::
    symToNosym() {
    try {
      if (this->nrows() == this->ncols()) {
	if (!this->is_zero()) {
	  /* Fix the diagonal: */
	  for (int rc = 0; rc < this->ncols(); rc++)
	    (*this)(rc, rc).symToNosym();
	  /* Fix the lower triangle */
	  for (int row = 1; row < this->nrows(); row++)
	    for (int col = 0; col < row; col++)
	      Telement::transpose((*this)(col, row), (*this)(row, col));
	}
      }
      else
	throw Failure("Matrix<Treal, Telement>::symToNosym(): "
		      "Only quadratic matrices can be symmetric");      
    }
    catch(Failure& f) {
      std::cout<<"Failure caught:Matrix<Treal, Telement>::symToNosym()"
	       <<std::endl;
      throw f;
    }
  }

  template<class Treal, class Telement> 
    void Matrix<Treal, Telement>::
    nosymToSym() {
    if (this->nrows() == this->ncols()) {
      if (!this->is_zero()) {
	/* Fix the diagonal: */
	for (int rc = 0; rc < this->ncols(); rc++)
	  (*this)(rc, rc).nosymToSym();
	/* Fix the lower triangle */
	for (int col = 0; col < this->ncols() - 1; col++)
	  for (int row = col + 1; row < this->nrows(); row++)
	    (*this)(row, col).clear();
      }    
    }
    else
      throw Failure("Matrix<Treal, Telement>::nosymToSym(): "
		    "Only quadratic matrices can be symmetric");      
  }
  
  /* Basic linear algebra routines */

  template<class Treal, class Telement>
    Matrix<Treal, Telement>& 
    Matrix<Treal, Telement>::operator=(int const k) {
    switch (k) {
    case 0:
      this->clear();
      break;
    case 1:
      if (this->ncols() != this->nrows()) 
	throw Failure
	  ("Matrix::operator=(int k = 1): "
	   "Matrix must be quadratic to become identity matrix.");
      this->clear();
      this->allocate();
      for (int rc = 0; rc < this->ncols(); rc++) /*Set diagonal to identity*/
	(*this)(rc,rc) = 1;
      break;
    default:
      throw Failure("Matrix::operator=(int k) only "
		    "implemented for k = 0, k = 1");
    }
    return *this;
  }


  template<class Treal, class Telement>
    Matrix<Treal, Telement>& Matrix<Treal, Telement>:: 
    operator*=(const Treal alpha) {
    if (!this->is_zero() && alpha != 1) {
      for (int ind = 0; ind < this->nElements(); ind++)
	this->elements[ind] *= alpha;
    }
    return *this;
  }

  /* C = beta * C + alpha * A * B */
  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::
    gemm(const bool tA, const bool tB, const Treal alpha, 
	 const Matrix<Treal, Telement>& A, 
	 const Matrix<Treal, Telement>& B, 
	 const Treal beta, 
	 Matrix<Treal, Telement>& C) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (C.is_empty()) {
      assert(beta == 0);
      if (tA)
	C.resetRows(A.cols);
      else
	C.resetRows(A.rows);
      if (tB)
	C.resetCols(B.rows);
      else
	C.resetCols(B.cols);
    }      
    
    if ((A.is_zero() || B.is_zero() || alpha == 0) && 
	(C.is_zero() || beta == 0))
      C = 0;
    else {
      Treal beta_tmp = beta;
      if (C.is_zero()) {
	C.allocate();
	beta_tmp = 0;
      }
      if (!A.is_zero() && !B.is_zero() && alpha != 0) {
	MAT_OMP_INIT;
	if (!tA && !tB) 
	  if (A.ncols() == B.nrows() && 
	      A.nrows() == C.nrows() &&
	      B.ncols() == C.ncols()) {	
	    int C_ncols = C.ncols();
#ifdef _OPENMP
#pragma omp parallel for if(A.level() == Params::getMatrixParallelLevel()) num_threads(Params::getNProcs()) default(shared) schedule(dynamic) 
#endif
	    for (int col = 0; col < C_ncols; col++) {
	      MAT_OMP_START;
	      for (int row = 0; row < C.nrows(); row++) {
		Telement::gemm(tA, tB, alpha, 
			       A(row, 0), B(0, col), 
			       beta_tmp, 
			       C(row, col));
		for(int cola = 1; cola < A.ncols(); cola++)
		  Telement::gemm(tA, tB, alpha, 
				 A(row, cola), B(cola, col), 
				 1.0, 
				 C(row, col));
	      }
	      MAT_OMP_END;
	    } /* end omp for */
	  }
	  else 
	    throw Failure("Matrix<class Treal, class Telement>::"
			  "gemm(bool, bool, Treal, "
			  "Matrix<Treal, Telement>, "
			  "Matrix<Treal, Telement>, Treal, "
			  "Matrix<Treal, Telement>): "
			  "Incorrect matrixdimensions for multiplication");
	else if (tA && !tB)
	  if (A.nrows() == B.nrows() &&
	      A.ncols() == C.nrows() &&
	      B.ncols() == C.ncols()) {
	    int C_ncols = C.ncols();
#ifdef _OPENMP
#pragma omp parallel for if(A.level() == Params::getMatrixParallelLevel()) num_threads(Params::getNProcs()) default(shared) schedule(dynamic) 
#endif
	    for (int col = 0; col < C_ncols; col++) {
	      MAT_OMP_START;
	      for (int row = 0; row < C.nrows(); row++) {
		Telement::gemm(tA, tB, alpha,
			       A(0,row), B(0,col),
			       beta_tmp,
			       C(row,col));
		for(int cola = 1; cola < A.nrows(); cola++)
		  Telement::gemm(tA, tB, alpha, 
				 A(cola, row), B(cola, col), 
				 1.0, 
				 C(row,col));
	      }
	      MAT_OMP_END;
	    } /* end omp for */
	  }
	  else
	    throw Failure("Matrix<class Treal, class Telement>::"
			  "gemm(bool, bool, Treal, "
			  "Matrix<Treal, Telement>, "
			  "Matrix<Treal, Telement>, Treal, "
			  "Matrix<Treal, Telement>): "
			  "Incorrect matrixdimensions for multiplication");
	else if (!tA && tB)
	  if (A.ncols() == B.ncols() && 
	      A.nrows() == C.nrows() &&
	      B.nrows() == C.ncols()) {
	    int C_ncols = C.ncols();
#ifdef _OPENMP
#pragma omp parallel for if(A.level() == Params::getMatrixParallelLevel()) num_threads(Params::getNProcs()) default(shared) schedule(dynamic) 
#endif
	    for (int col = 0; col < C_ncols; col++) {
	      MAT_OMP_START;
	      for (int row = 0; row < C.nrows(); row++) {
		Telement::gemm(tA, tB, alpha, 
			       A(row, 0), B(col, 0), 
			       beta_tmp, 
			       C(row,col));
		for(int cola = 1; cola < A.ncols(); cola++)
		  Telement::gemm(tA, tB, alpha, 
				 A(row, cola), B(col, cola), 
				 1.0, 
				 C(row,col));
	      }
	      MAT_OMP_END;
	    } /* end omp for */
	  }
      	  else 
	    throw Failure("Matrix<class Treal, class Telement>::"
			  "gemm(bool, bool, Treal, "
			  "Matrix<Treal, Telement>, "
			  "Matrix<Treal, Telement>, Treal, "
			  "Matrix<Treal, Telement>): "
			  "Incorrect matrixdimensions for multiplication");
	else if (tA && tB)
	  if (A.nrows() == B.ncols() && 
	      A.ncols() == C.nrows() &&
              B.nrows() == C.ncols()) {
	    int C_ncols = C.ncols();
#ifdef _OPENMP
#pragma omp parallel for if(A.level() == Params::getMatrixParallelLevel()) num_threads(Params::getNProcs()) default(shared) schedule(dynamic) 
#endif
	    for (int col = 0; col < C_ncols; col++) {
	      MAT_OMP_START;
	      for (int row = 0; row < C.nrows(); row++) {
		Telement::gemm(tA, tB, alpha, 
			       A(0, row), B(col, 0), 
			       beta_tmp, 
			       C(row,col));
		for(int cola = 1; cola < A.nrows(); cola++)
		  Telement::gemm(tA, tB, alpha, 
				 A(cola, row), B(col, cola), 
				 1.0, 
				 C(row,col));
	      }
	      MAT_OMP_END;
	    } /* end omp for */
	  }
	  else 
	    throw Failure("Matrix<class Treal, class Telement>::"
			  "gemm(bool, bool, Treal, "
			  "Matrix<Treal, Telement>, "
			  "Matrix<Treal, Telement>, "
			  "Treal, Matrix<Treal, Telement>): "
			  "Incorrect matrixdimensions for multiplication");
	else throw Failure("Matrix<class Treal, class Telement>::"
			   "gemm(bool, bool, Treal, "
			   "Matrix<Treal, Telement>, "
			   "Matrix<Treal, Telement>, Treal, "
			   "Matrix<Treal, Telement>):"
			   "Very strange error!!");
	MAT_OMP_FINALIZE;
      }    
      else
	C *= beta;
    }
  }

  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::
    symm(const char side, const char uplo, const Treal alpha, 
	 const Matrix<Treal, Telement>& A, 
	 const Matrix<Treal, Telement>& B, 
	 const Treal beta, 
	 Matrix<Treal, Telement>& C) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    assert(A.nrows() == A.ncols());
    int dimA = A.nrows();
    if (C.is_empty()) {
      assert(beta == 0);
      if (side =='L') {
	C.resetRows(A.rows);
	C.resetCols(B.cols);
      } 
      else {
	assert(side == 'R');
	C.resetRows(B.rows);
	C.resetCols(A.cols);
      }
    }      
    if ((A.is_zero() || B.is_zero() || alpha == 0) && 
	(C.is_zero() || beta == 0))
      C = 0;
    else {
      if (uplo == 'U') {
	Treal beta_tmp = beta;
	if (C.is_zero()) {
	  C.allocate();
	  beta_tmp = 0;
	}
	if (!A.is_zero() && !B.is_zero() && alpha != 0) {
	  MAT_OMP_INIT;
	  if (side =='L') 
	    if (dimA == B.nrows() && 
		dimA == C.nrows() &&
		B.ncols() == C.ncols()) {
	      int C_ncols = C.ncols();
#ifdef _OPENMP
#pragma omp parallel for if(A.level() == Params::getMatrixParallelLevel()) num_threads(Params::getNProcs()) default(shared) schedule(dynamic) 
#endif
	      for (int col = 0; col < C_ncols; col++) {
		MAT_OMP_START;
		for (int row = 0; row < C.nrows(); row++) {
		  /* Diagonal element in matrix A */
		  Telement::symm(side, uplo, alpha, A(row, row),
				 B(row, col), beta_tmp, C(row, col));
		  /* Elements below diagonal in A */
		  for(int ind = 0; ind < row; ind++)
		    Telement::gemm(true, false, alpha, 
				   A(ind, row), B(ind, col), 
				   1.0, C(row,col));
		  /* Elements above diagonal in A */
		  for(int ind = row + 1; ind < dimA; ind++)
		    Telement::gemm(false, false, alpha, 
				   A(row, ind), B(ind, col), 
				   1.0, C(row,col));
		}
		MAT_OMP_END;
	      } /* end omp for */
	    }
	    else 
	      throw Failure("Matrix<class Treal, class Telement>"
			    "::symm(bool, bool, Treal, Matrix<Treal, "
			    "Telement>, Matrix<Treal, Telement>, "
			    "Treal, Matrix<Treal, Telement>): "
			    "Incorrect matrixdimensions for multiplication");
	  else { /* side == 'R' */
	    assert(side == 'R');
	    if (B.ncols() == dimA && 
		B.nrows() == C.nrows() &&
		dimA == C.ncols()) {
	      int C_ncols = C.ncols();
#ifdef _OPENMP
#pragma omp parallel for if(A.level() == Params::getMatrixParallelLevel()) num_threads(Params::getNProcs()) default(shared) schedule(dynamic) 
#endif
	      for (int col = 0; col < C_ncols; col++) {
		MAT_OMP_START;
		for (int row = 0; row < C.nrows(); row++) {
		  /* Diagonal element in matrix A */
		  Telement::symm(side, uplo, alpha, A(col, col),
				 B(row, col), beta_tmp, C(row, col));
		  /* Elements below diagonal in A */
		  for(int ind = col + 1; ind < dimA; ind++)
		    Telement::gemm(false, true, alpha, 
				   B(row, ind), A(col, ind),  
				   1.0, C(row,col));
		  /* Elements above diagonal in A */
		  for(int ind = 0; ind < col; ind++)
		    Telement::gemm(false, false, alpha, 
				   B(row, ind), A(ind, col),  
				   1.0, C(row,col));
		}
		MAT_OMP_END;
	      } /* end omp for */
	    }
	    else 
	      throw Failure("Matrix<class Treal, class Telement>"
			    "::symm(bool, bool, Treal, Matrix<Treal, "
			    "Telement>, Matrix<Treal, Telement>, "
			    "Treal, Matrix<Treal, Telement>): "
			    "Incorrect matrixdimensions for multiplication");
	  }
	  MAT_OMP_FINALIZE;
	}
	else
	  C *= beta;
      }
      else
	throw Failure("Matrix<class Treal, class Telement>::"
		      "symm only implemented for symmetric matrices in "
		      "upper triangular storage");
    }    
  }

  template<class Treal, class Telement>
    void Matrix<Treal, Telement>:: 
    syrk(const char uplo, const bool tA, const Treal alpha, 
	 const Matrix<Treal, Telement>& A, 
	 const Treal beta, 
	 Matrix<Treal, Telement>& C) {
    assert(!A.is_empty());
    if (C.is_empty()) {
      assert(beta == 0);
      if (tA) {
	C.resetRows(A.cols);
	C.resetCols(A.cols);
      } 
      else {
	C.resetRows(A.rows);
	C.resetCols(A.rows);
      }
    }
    
    if (C.nrows() == C.ncols() && 
	((!tA && A.nrows() == C.nrows()) || (tA && A.ncols() == C.nrows())))
      if (alpha != 0 && !A.is_zero()) {
	Treal beta_tmp = beta;
	if (C.is_zero()) {
	  C.allocate();
	  beta_tmp = 0;
	}
	MAT_OMP_INIT;
	if (!tA && uplo == 'U') { /* C = alpha * A * A' + beta * C */
#ifdef _OPENMP
#pragma omp parallel if(A.level() == Params::getMatrixParallelLevel()) num_threads(Params::getNProcs()) default(shared)
#endif
	  {
	    int C_ncols = C.ncols();
#ifdef _OPENMP
#pragma omp for  schedule(dynamic) nowait
#endif
	    for (int rc = 0; rc < C_ncols; rc++) {
	      MAT_OMP_START;
	      Telement::syrk(uplo, tA, alpha, A(rc, 0), beta_tmp, C(rc, rc));
	      for (int cola = 1; cola < A.ncols(); cola++)
		Telement::syrk(uplo, tA, alpha, A(rc, cola), 1.0, C(rc, rc));
	      MAT_OMP_END;
	    }
	    int C_nrows = C.nrows();
#ifdef _OPENMP
#pragma omp for  schedule(dynamic) nowait
#endif
	    for (int row = 0; row < C_nrows - 1; row++) {
	      MAT_OMP_START;
	      for (int col = row + 1; col < C.ncols(); col++) {
		Telement::gemm(tA, !tA, alpha, 
			       A(row, 0), A(col,0), beta_tmp, C(row,col));
		for (int ind = 1; ind < A.ncols(); ind++)
		  Telement::gemm(tA, !tA, alpha, 
				 A(row, ind), A(col,ind), 1.0, C(row,col));
	      }
	      MAT_OMP_END;
	    }
	  } /* end omp parallel */
	} /* end if (!tA) */
	else if (tA && uplo == 'U') { /* C = alpha * A' * A + beta * C */
#ifdef _OPENMP
#pragma omp parallel if(A.level() == Params::getMatrixParallelLevel()) num_threads(Params::getNProcs()) default(shared)
#endif
	  {
	    int C_ncols = C.ncols();
#ifdef _OPENMP
#pragma omp for  schedule(dynamic) nowait
#endif
	    for (int rc = 0; rc < C_ncols; rc++) {
	      MAT_OMP_START;
	      Telement::syrk(uplo, tA, alpha, A(0, rc), beta_tmp, C(rc, rc));
	      for (int rowa = 1; rowa < A.nrows(); rowa++)
		Telement::syrk(uplo, tA, alpha, A(rowa, rc), 1.0, C(rc, rc));
	      MAT_OMP_END;
	    }
	    int C_nrows = C.nrows();
#ifdef _OPENMP
#pragma omp for  schedule(dynamic) nowait
#endif
	    for (int row = 0; row < C_nrows - 1; row++) {
	      MAT_OMP_START;
	      for (int col = row + 1; col < C.ncols(); col++) {
		Telement::gemm(tA, !tA, alpha, 
			       A(0, row), A(0, col), beta_tmp, C(row,col));
		for (int ind = 1; ind < A.nrows(); ind++)
		  Telement::gemm(tA, !tA, alpha, 
				 A(ind, row), A(ind, col), 1.0, C(row,col));
	      }
	      MAT_OMP_END;
	    }
	  } /* end omp parallel */
	} /* end if (tA) */  
	else
	  throw Failure("Matrix<class Treal, class Telement>::"
			"syrk not implemented for lower triangular storage");
	MAT_OMP_FINALIZE;
      }
      else {
	C *= beta;
      }
    else
      throw Failure("Matrix<class Treal, class Telement>::syrk: "
		    "Incorrect matrix dimensions for symmetric rank-k update");
  }

  template<class Treal, class Telement>
    void Matrix<Treal, Telement>:: 
    sysq(const char uplo, const Treal alpha, 
	 const Matrix<Treal, Telement>& A, 
	 const Treal beta, 
	 Matrix<Treal, Telement>& C) {
    assert(!A.is_empty());
    if (C.is_empty()) {
      assert(beta == 0);
      C.resetRows(A.rows);
      C.resetCols(A.cols);
    }
    if (C.nrows() == C.ncols() && 
        A.nrows() == C.nrows() && A.nrows() == A.ncols())
      if (alpha != 0 && !A.is_zero()) {
	if (uplo == 'U') {
	  Treal beta_tmp = beta;
	  if (C.is_zero()) {
	    C.allocate();
	    beta_tmp = 0;
	  }
	  MAT_OMP_INIT;
#ifdef _OPENMP
#pragma omp parallel if(A.level() == Params::getMatrixParallelLevel()) num_threads(Params::getNProcs()) default(shared)
#endif
	  {
	    int C_ncols = C.ncols();
#ifdef _OPENMP
#pragma omp for  schedule(dynamic) nowait
#endif
	    for (int rc = 0; rc < C_ncols; rc++) {
	      MAT_OMP_START;
	      Telement::sysq(uplo, alpha, A(rc, rc), beta_tmp, C(rc, rc));
	      for (int cola = 0; cola < rc; cola++)
		Telement::syrk(uplo, true, alpha, A(cola, rc), 1.0, C(rc, rc));
	      for (int cola = rc + 1; cola < A.ncols(); cola++)
		Telement::syrk(uplo, false, alpha, A(rc, cola), 1.0, C(rc, rc));
	      MAT_OMP_END;
	    }
	    /* Maste anvanda symm? */
	    int C_nrows = C.nrows();
#ifdef _OPENMP
#pragma omp for schedule(dynamic) nowait
#endif
	    for (int row = 0; row < C_nrows - 1; row++) {
	      MAT_OMP_START;
	      for (int col = row + 1; col < C.ncols(); col++) {
		/* First the two operations involving diagonal elements */
		Telement::symm('L', 'U', alpha, A(row, row), A(row, col), 
			       beta_tmp, C(row, col));
		Telement::symm('R', 'U', alpha, A(col, col), A(row, col), 
			       1.0, C(row, col));
		/* First element in product is below  the diagonal */
		for (int ind = 0; ind < row; ind++)
		  Telement::gemm(true, false, alpha, 
				 A(ind, row), A(ind, col), 1.0, C(row, col));
		/* None of the elements are below the diagonal     */
		for (int ind = row + 1; ind < col; ind++)
		  Telement::gemm(false, false, alpha, 
				 A(row, ind), A(ind, col), 1.0, C(row, col));
		/* Second element is below the diagonal            */
		for (int ind = col + 1; ind < A.ncols(); ind++)
		  Telement::gemm(false, true, alpha, 
				 A(row, ind), A(col, ind), 1.0, C(row, col));
	      }
	      MAT_OMP_END;
	    } /* end omp for */
	  } /* end omp parallel */
	  MAT_OMP_FINALIZE;
	}
	else
	  throw Failure("Matrix<class Treal, class Telement>::"
			"sysq only implemented for symmetric matrices in "
			"upper triangular storage");;
      }
      else {
	C *= beta;
      }
    else
      throw Failure("Matrix<class Treal, class Telement>::sysq: "
		    "Incorrect matrix dimensions for symmetric square "
		    "operation");
  }

  /* C = alpha * A * B + beta * C where A and B are symmetric */
  template<class Treal, class Telement>
    void Matrix<Treal, Telement>:: 
    ssmm(const Treal alpha, 
	 const Matrix<Treal, Telement>& A, 
	 const Matrix<Treal, Telement>& B, 
	 const Treal beta, 
	 Matrix<Treal, Telement>& C) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (C.is_empty()) {
      assert(beta == 0);
      C.resetRows(A.rows);
      C.resetCols(B.cols);
    }
    if (A.ncols() != B.nrows() || 
	A.nrows() != C.nrows() ||
	B.ncols() != C.ncols() ||
	A.nrows() != A.ncols() ||
	B.nrows() != B.ncols()) {
      throw Failure("Matrix<class Treal, class Telement>::"
		    "ssmm(Treal, "
		    "Matrix<Treal, Telement>, "
		    "Matrix<Treal, Telement>, Treal, "
		    "Matrix<Treal, Telement>): "
		    "Incorrect matrixdimensions for ssmm multiplication");
    }
    if ((A.is_zero() || B.is_zero() || alpha == 0) && 
	(C.is_zero() || beta == 0)) {
      C = 0;
      return;
    }
    Treal beta_tmp = beta;
    if (C.is_zero()) {
      C.allocate();
      beta_tmp = 0;
    }
    if (A.is_zero() || B.is_zero() || alpha == 0) {
      C *= beta;
      return;
    }
    MAT_OMP_INIT;
    int C_ncols = C.ncols();
#ifdef _OPENMP
#pragma omp parallel for if(A.level() == Params::getMatrixParallelLevel()) num_threads(Params::getNProcs()) default(shared) schedule(dynamic) 
#endif
    for (int col = 0; col < C_ncols; col++) {
      MAT_OMP_START;
      /* Diagonal */
      /* C(col, col) = alpha * A(col, col) * B(col, col) + beta * C(col, col)*/
      Telement::ssmm(alpha, A(col,col), B(col, col),
		     beta_tmp, C(col,col));
      for (int ind = 0; ind < col; ind++)
	/* C(col, col) += alpha * A(ind, col)' * B(ind, col) */
	Telement::gemm(true, false, 
		       alpha, A(ind,col), B(ind, col),
		       1.0, C(col,col));
      for (int ind = col + 1; ind < A.ncols(); ind++)
	/* C(col, col) += alpha * A(col, ind) * B(col, ind)' */
	Telement::gemm(false, true,  
		       alpha, A(col, ind), B(col, ind),
		       1.0, C(col,col));
      /* Upper half */
      for (int row = 0; row < col; row++) {
	/* C(row, col) = alpha * A(row, row) * B(row, col) + 
	 *               beta * C(row, col)
	 */
	Telement::symm('L', 'U', 
		       alpha, A(row, row), B(row, col), 
		       beta_tmp, C(row, col));
	/* C(row, col) += alpha * A(row, col) * B(col, col)	 */
	Telement::symm('R', 'U', 
		       alpha, B(col, col), A(row, col), 
		       1.0, C(row, col));
	for (int ind = 0; ind < row; ind++)
	  /* C(row, col) += alpha * A(ind, row)' * B(ind, col)	 */
	  Telement::gemm(true, false,  
			 alpha, A(ind, row), B(ind, col),
			 1.0, C(row,col));
	for (int ind = row + 1; ind < col; ind++)
	  /* C(row, col) += alpha * A(row, ind) * B(ind, col)	 */
	  Telement::gemm(false, false,  
			 alpha, A(row, ind), B(ind, col),
			 1.0, C(row,col));
	for (int ind = col + 1; ind < A.ncols(); ind++)
	  /* C(row, col) += alpha * A(row, ind) * B(col, ind)'	 */
	  Telement::gemm(false, true,  
			 alpha, A(row, ind), B(col, ind),
			 1.0, C(row,col));
      }
      /* Lower half */
      Telement tmpSubMat;
      for (int row = col + 1; row < C.nrows(); row++) {
	Telement::transpose(C(row, col), tmpSubMat);
	/* tmpSubMat = alpha * B(col, col) * A(col, row) + 
	 *             beta * tmpSubMat
	 */
	Telement::symm('L', 'U', 
		       alpha, B(col, col), A(col, row), 
		       beta_tmp, tmpSubMat);
	/* tmpSubMat += alpha * B(col, row) * A(row, row)	 */
	Telement::symm('R', 'U', 
		       alpha, A(row, row), B(col, row), 
		       1.0, tmpSubMat);
	for (int ind = 0; ind < col; ind++)
	  /* tmpSubMat += alpha * B(ind, col)' * A(ind, row)	 */
	  Telement::gemm(true, false,  
			 alpha, B(ind, col), A(ind, row),
			 1.0, tmpSubMat);
	for (int ind = col + 1; ind < row; ind++)
	  /* tmpSubMat += alpha * B(col, ind) * A(ind, row)	 */
	  Telement::gemm(false, false,  
			 alpha, B(col, ind), A(ind, row),
			 1.0, tmpSubMat);
	for (int ind = row + 1; ind < B.nrows(); ind++)
	  /* tmpSubMat += alpha * B(col, ind) * A(row, ind)'	 */
	  Telement::gemm(false, true,  
			 alpha, B(col, ind), A(row, ind),
			 1.0, tmpSubMat);
	Telement::transpose(tmpSubMat, C(row, col));
      }
      MAT_OMP_END;
    }
    MAT_OMP_FINALIZE;
  } /* end ssmm */
  

    /* C = alpha * A * B + beta * C where A and B are symmetric
     * and only the upper triangle of C is computed.
     */
  template<class Treal, class Telement>
    void Matrix<Treal, Telement>:: 
    ssmm_upper_tr_only(const Treal alpha, 
		       const Matrix<Treal, Telement>& A, 
		       const Matrix<Treal, Telement>& B, 
		       const Treal beta, 
		       Matrix<Treal, Telement>& C) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (C.is_empty()) {
      assert(beta == 0);
      C.resetRows(A.rows);
      C.resetCols(B.cols);
    }
    if (A.ncols() != B.nrows() || 
	A.nrows() != C.nrows() ||
	B.ncols() != C.ncols() ||
	A.nrows() != A.ncols() ||
	B.nrows() != B.ncols()) {
      throw Failure("Matrix<class Treal, class Telement>::"
		    "ssmm_upper_tr_only(Treal, "
		    "Matrix<Treal, Telement>, "
		    "Matrix<Treal, Telement>, Treal, "
		    "Matrix<Treal, Telement>): "
		    "Incorrect matrixdimensions for ssmm_upper_tr_only "
		    "multiplication");
    }
    if ((A.is_zero() || B.is_zero() || alpha == 0) && 
	(C.is_zero() || beta == 0)) {
      C = 0;
      return;
    }
    Treal beta_tmp = beta;
    if (C.is_zero()) {
      C.allocate();
      beta_tmp = 0;
    }
    if (A.is_zero() || B.is_zero() || alpha == 0) {
      C *= beta;
      return;
    }
    MAT_OMP_INIT;
    int C_ncols = C.ncols();
#ifdef _OPENMP
#pragma omp parallel for if(A.level() == Params::getMatrixParallelLevel()) num_threads(Params::getNProcs()) default(shared) schedule(dynamic) 
#endif
    for (int col = 0; col < C_ncols; col++) {
      MAT_OMP_START;
      /* Diagonal */
      /* C(col, col) = alpha * A(col, col) * B(col, col) + beta * C(col, col)*/
      Telement::ssmm_upper_tr_only(alpha, A(col,col), B(col, col),
				   beta_tmp, C(col,col));
      for (int ind = 0; ind < col; ind++)
	/* C(col, col) += alpha * A(ind, col)' * B(ind, col) */
	Telement::gemm_upper_tr_only(true, false, 
				     alpha, A(ind,col), B(ind, col),
				     1.0, C(col,col));
      for (int ind = col + 1; ind < A.ncols(); ind++)
	/* C(col, col) += alpha * A(col, ind) * B(col, ind)' */
	Telement::gemm_upper_tr_only(false, true,  
				     alpha, A(col, ind), B(col, ind),
				     1.0, C(col,col));
      /* Upper half */
      for (int row = 0; row < col; row++) {
	/* C(row, col) = alpha * A(row, row) * B(row, col) + 
	 *               beta * C(row, col)
	 */
	Telement::symm('L', 'U', 
		       alpha, A(row, row), B(row, col), 
		       beta_tmp, C(row, col));
	/* C(row, col) += alpha * A(row, col) * B(col, col)	 */
	Telement::symm('R', 'U', 
		       alpha, B(col, col), A(row, col), 
		       1.0, C(row, col));
	for (int ind = 0; ind < row; ind++)
	  /* C(row, col) += alpha * A(ind, row)' * B(ind, col)	 */
	  Telement::gemm(true, false,  
			 alpha, A(ind, row), B(ind, col),
			 1.0, C(row,col));
	for (int ind = row + 1; ind < col; ind++)
	  /* C(row, col) += alpha * A(row, ind) * B(ind, col)	 */
	  Telement::gemm(false, false,  
			 alpha, A(row, ind), B(ind, col),
			 1.0, C(row,col));
	for (int ind = col + 1; ind < A.ncols(); ind++)
	  /* C(row, col) += alpha * A(row, ind) * B(col, ind)'	 */
	  Telement::gemm(false, true,  
			 alpha, A(row, ind), B(col, ind),
			 1.0, C(row,col));
      }
      MAT_OMP_END;
    }
    MAT_OMP_FINALIZE;
  } /* end ssmm_upper_tr_only */
  


  template<class Treal, class Telement>
    void Matrix<Treal, Telement>:: 
    trmm(const char side, const char uplo, const bool tA, const Treal alpha, 
	 const Matrix<Treal, Telement>& A, 
	 Matrix<Treal, Telement>& B) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (alpha != 0 && !A.is_zero() && !B.is_zero())
      if (((side == 'R' && B.ncols() == A.nrows()) || 
	   (side == 'L' && A.ncols() == B.nrows())) &&
	  A.nrows() == A.ncols())
	if (uplo == 'U')
	  if (!tA) { 
	    if (side == 'R') {
	      /* Last column must be calculated first                        */
	      for (int col = B.ncols() - 1; col >= 0; col--)
		for (int row = 0; row < B.nrows(); row++) {
		  /* Use first the diagonal element in A                     */
		  /* Otherwise the faster call to trmm cannot be utilized    */
		  Telement::trmm(side, uplo, tA, alpha, 
				 A(col, col), B(row,col));
		  /* And the rest:                                           */
		  for (int ind = 0; ind < col; ind++)
		    Telement::gemm(false, tA, alpha, 
				   B(row,ind), A(ind, col), 
				   1.0, B(row,col));
		}
	    } /* end if (side == 'R')*/
	    else {
	      assert(side == 'L');
	      /* First row must be calculated first                          */
	      for (int row = 0; row < B.nrows(); row++ )
		for (int col = 0; col < B.ncols(); col++) {
		  Telement::trmm(side, uplo, tA, alpha, 
				 A(row,row), B(row,col));
		  for (int ind = row + 1 ; ind < B.nrows(); ind++)
		    Telement::gemm(tA, false, alpha,
				   A(row,ind), B(ind,col),
				   1.0, B(row,col));
		}
	    } /* end else (side == 'L')*/
	  } /* end if (tA == false) */
	  else {
	    assert(tA == true);
	    if (side == 'R') {
	      /* First column must be calculated first                       */
	      for (int col = 0; col < B.ncols(); col++)
		for (int row = 0; row < B.nrows(); row++) {
		  Telement::trmm(side, uplo, tA, alpha,
				 A(col,col), B(row,col));
		  for (int ind = col + 1; ind < A.ncols(); ind++)
		    Telement::gemm(false, tA, alpha,
				   B(row,ind), A(col,ind),
				   1.0, B(row,col));
		}
	    } /* end if (side == 'R')*/
	    else {
	      assert(side == 'L');
	      /* Last row must be calculated first                           */
	      for (int row = B.nrows() - 1; row >= 0; row--)
		for (int col = 0; col < B.ncols(); col++) {
		  Telement::trmm(side, uplo, tA, alpha,
				 A(row,row), B(row,col));
		  for (int ind = 0; ind < row; ind++)
		    Telement::gemm(tA, false, alpha,
				   A(ind,row), B(ind,col),
				   1.0, B(row,col));
		}
	    } /* end else (side == 'L')*/
	  } /* end else (tA == true)*/
	else
	  throw Failure("Matrix<class Treal, class Telement>::"
			"trmm not implemented for lower triangular matrices");
      else
	throw Failure("Matrix<class Treal, class Telement>::trmm"
		      ": Incorrect matrix dimensions for multiplication");
    else
      B = 0;
  }


  template<class Treal, class Telement>
    Treal Matrix<Treal, Telement>::frobSquared() const {
    assert(!this->is_empty());
    if (this->is_zero()) 
      return 0;
    else {
      Treal sum(0);
      for (int i = 0; i < this->nElements(); i++)
	sum += this->elements[i].frobSquared();
      return sum;
    }
  } 
  
  template<class Treal, class Telement>
    Treal Matrix<Treal, Telement>::
    syFrobSquared() const {
    assert(!this->is_empty());
    if (this->nrows() != this->ncols())
      throw Failure("Matrix<Treal, Telement>::syFrobSquared(): "
		    "Matrix must be have equal number of rows "
		    "and cols to be symmetric");
    Treal sum(0);
    if (!this->is_zero()) {
      for (int col = 1; col < this->ncols(); col++)
	for (int row = 0; row < col; row++)
	  sum += 2 * (*this)(row, col).frobSquared();
      for (int rc = 0; rc < this->ncols(); rc++)
	sum += (*this)(rc, rc).syFrobSquared();
    }
    return sum;
  }
  
  template<class Treal, class Telement>
    Treal Matrix<Treal, Telement>::
    frobSquaredDiff
    (const Matrix<Treal, Telement>& A,
     const Matrix<Treal, Telement>& B) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (A.nrows() != B.nrows() || A.ncols() != B.ncols()) 
      throw Failure("Matrix<Treal, Telement>::"
		    "frobSquaredDiff: Incorrect matrix dimensions");
    Treal sum(0);
    if (!A.is_zero() && !B.is_zero()) 
      for (int i = 0; i < A.nElements(); i++)
	sum += Telement::frobSquaredDiff(A.elements[i],B.elements[i]);
    else if (A.is_zero() && !B.is_zero()) 
      sum = B.frobSquared();
    else if (!A.is_zero() && B.is_zero())
      sum = A.frobSquared();
    /* sum is already zero if A.is_zero() && B.is_zero() */
    return sum;
  }  
  
  template<class Treal, class Telement>
    Treal Matrix<Treal, Telement>::
    syFrobSquaredDiff
    (const Matrix<Treal, Telement>& A,
     const Matrix<Treal, Telement>& B) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (A.nrows() != B.nrows() ||
	A.ncols() != B.ncols() || 
	A.nrows() != A.ncols()) 
      throw Failure("Matrix<Treal, Telement>::syFrobSquaredDiff: "
		    "Incorrect matrix dimensions");
    Treal sum(0);
    if (!A.is_zero() && !B.is_zero()) {
      for (int col = 1; col < A.ncols(); col++)
	for (int row = 0; row < col; row++)
	  sum += 2 * Telement::frobSquaredDiff(A(row, col), B(row, col));
      for (int rc = 0; rc < A.ncols(); rc++)
	sum += Telement::syFrobSquaredDiff(A(rc, rc),B(rc, rc));
    }
    else if (A.is_zero() && !B.is_zero()) 
      sum = B.syFrobSquared();
    else if (!A.is_zero() && B.is_zero())
      sum = A.syFrobSquared();
    /* sum is already zero if A.is_zero() && B.is_zero() */
    return sum;
  }
  
  template<class Treal, class Telement> 
    Treal Matrix<Treal, Telement>::
    trace() const {
    assert(!this->is_empty());
    if (this->nrows() != this->ncols())
      throw Failure("Matrix<Treal, Telement>::trace(): "
		    "Matrix must be quadratic"); 	  
    if (this->is_zero())
      return 0;
    else {
      Treal sum = 0;
      for (int rc = 0; rc < this->ncols(); rc++)
	sum += (*this)(rc,rc).trace();
      return sum;
    }
  }
  
  template<class Treal, class Telement> 
    Treal Matrix<Treal, Telement>::
    trace_ab(const Matrix<Treal, Telement>& A,
	     const Matrix<Treal, Telement>& B) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (A.nrows() != B.ncols() || A.ncols() != B.nrows()) 
      throw Failure("Matrix<Treal, Telement>::trace_ab: "
		    "Wrong matrix dimensions for trace(A * B)"); 
    Treal tr = 0;
    if (!A.is_zero() && !B.is_zero())
      for (int rc = 0; rc < A.nrows(); rc++)
	for (int colA = 0; colA < A.ncols(); colA++)
	  tr += Telement::trace_ab(A(rc,colA), B(colA, rc));
    return tr;
  } 
  
  template<class Treal, class Telement> 
    Treal  Matrix<Treal, Telement>::
    trace_aTb(const Matrix<Treal, Telement>& A,
	      const Matrix<Treal, Telement>& B) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (A.ncols() != B.ncols() || A.nrows() != B.nrows())
      throw Failure("Matrix<Treal, Telement>::trace_aTb: "
		    "Wrong matrix dimensions for trace(A' * B)"); 
    Treal tr = 0;
    if (!A.is_zero() && !B.is_zero())
      for (int rc = 0; rc < A.ncols(); rc++)
	for (int rowB = 0; rowB < B.nrows(); rowB++)
	  tr += Telement::trace_aTb(A(rowB,rc), B(rowB, rc));
    return tr;
  }
  
  template<class Treal, class Telement> 
    Treal Matrix<Treal, Telement>::
    sy_trace_ab(const Matrix<Treal, Telement>& A,
		const Matrix<Treal, Telement>& B) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (A.nrows() != B.ncols() || A.ncols() != B.nrows() || 
	A.nrows() != A.ncols()) 
      throw Failure("Matrix<Treal, Telement>::sy_trace_ab: "
		    "Wrong matrix dimensions for symmetric trace(A * B)");
    Treal tr = 0;
    if (!A.is_zero() && !B.is_zero()) {
      /* Diagonal first */
      for (int rc = 0; rc < A.nrows(); rc++)
	tr += Telement::sy_trace_ab(A(rc,rc), B(rc, rc));
      /* Using that trace of transpose is equal to that without transpose: */
      for (int colA = 1; colA < A.ncols(); colA++)
	for (int rowA = 0; rowA < colA; rowA++)
	  tr += 2 * Telement::trace_aTb(A(rowA, colA), B(rowA, colA));
    }
    return tr;
  }

  template<class Treal, class Telement>
    void Matrix<Treal, Telement>:: 
    add(const Treal alpha,   /* B += alpha * A */
	const Matrix<Treal, Telement>& A, 
	Matrix<Treal, Telement>& B) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (A.nrows() != B.nrows() || A.ncols() != B.ncols()) 
      throw Failure("Matrix<Treal, Telement>::add: "
		    "Wrong matrix dimensions for addition");
    if (!A.is_zero() && alpha != 0) {
      if (B.is_zero())
	B.allocate();
      for (int ind = 0; ind < A.nElements(); ind++)
	Telement::add(alpha, A.elements[ind], B.elements[ind]);
    }
  }

  template<class Treal, class Telement>
    void Matrix<Treal, Telement>:: 
    assign(Treal const  alpha,   /* *this = alpha * A */
	   Matrix<Treal, Telement> const & A) {
    assert(!A.is_empty());
    if (this->is_empty()) {
      this->resetRows(A.rows);
      this->resetCols(A.cols);
    }
    *this = 0;
    Matrix<Treal, Telement>:: 
      add(alpha, A, *this);
  }


  /********** Help functions for thresholding */

  
  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::
    getFrobSqLowestLevel(std::vector<Treal> & frobsq) const {
    if (!this->is_zero())
      for (int ind = 0; ind < this->nElements(); ind++)
	this->elements[ind].getFrobSqLowestLevel(frobsq);
  }

  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::
    frobThreshLowestLevel
    (Treal const threshold, Matrix<Treal, Telement> * ErrorMatrix) {
    assert(!this->is_empty());
    bool thisMatIsZero = true;
    if (ErrorMatrix) {
      assert(!ErrorMatrix->is_empty());
      bool EMatIsZero = true;
      if (!ErrorMatrix->is_zero() || !this->is_zero()) {
	if (ErrorMatrix->is_zero())
	  ErrorMatrix->allocate();
	if (this->is_zero())
	  this->allocate();
	for (int ind = 0; ind < this->nElements(); ind++) {
	  this->elements[ind].
	    frobThreshLowestLevel(threshold, &ErrorMatrix->elements[ind]);
	  thisMatIsZero = thisMatIsZero && this->elements[ind].is_zero();
	  EMatIsZero    = EMatIsZero && ErrorMatrix->elements[ind].is_zero();
	}
	if (thisMatIsZero)
	  this->clear();
	if (EMatIsZero)
	  ErrorMatrix->clear();
      }
    }      
    else
      if (!this->is_zero()) {
	for (int ind = 0; ind < this->nElements(); ind++) {
	  this->elements[ind].frobThreshLowestLevel(threshold, 0);
	  thisMatIsZero = thisMatIsZero && this->elements[ind].is_zero();
	}
	if (thisMatIsZero)
	  this->clear();
      }
  }

  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::
    getFrobSqElementLevel(std::vector<Treal> & frobsq) const {
    if (!this->is_zero())
      for (int ind = 0; ind < this->nElements(); ind++)
	this->elements[ind].getFrobSqElementLevel(frobsq);
  }

  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::
    frobThreshElementLevel
    (Treal const threshold, Matrix<Treal, Telement> * ErrorMatrix) {
    assert(!this->is_empty());
    bool thisMatIsZero = true;
    if (ErrorMatrix) {
      assert(!ErrorMatrix->is_empty());
      bool EMatIsZero = true;
      if (!ErrorMatrix->is_zero() || !this->is_zero()) {
	if (ErrorMatrix->is_zero())
	  ErrorMatrix->allocate();
	if (this->is_zero())
	  this->allocate();
	for (int ind = 0; ind < this->nElements(); ind++) {
	  this->elements[ind].
	    frobThreshElementLevel(threshold, &ErrorMatrix->elements[ind]);
	  thisMatIsZero = thisMatIsZero && this->elements[ind].is_zero();
	  EMatIsZero    = EMatIsZero && ErrorMatrix->elements[ind].is_zero();
	}
	if (thisMatIsZero)
	  this->clear();
	if (EMatIsZero)
	  ErrorMatrix->clear();
      }
    }      
    else
      if (!this->is_zero()) {
	for (int ind = 0; ind < this->nElements(); ind++) {
	  this->elements[ind].frobThreshElementLevel(threshold, 0);
	  thisMatIsZero = thisMatIsZero && this->elements[ind].is_zero();
	}
	if (thisMatIsZero)
	  this->clear();
      }    
  }
  


  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::assignFrobNormsLowestLevel
    ( Matrix<Treal, Matrix<Treal, Telement> > const & A ) {
    if (!A.is_zero()) {
      if ( this->is_zero() )
	this->allocate();
      assert( this->nElements() == A.nElements() );
      for (int ind = 0; ind < this->nElements(); ind++)
	this->elements[ind].assignFrobNormsLowestLevel(A[ind]);
    }
    else
      this->clear();
  } 

  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::syAssignFrobNormsLowestLevel
    ( Matrix<Treal, Matrix<Treal, Telement> > const & A ) {
    assert(!A.is_empty());
    if (A.nrows() != A.ncols())
      throw Failure("Matrix<Treal, Telement>::syAssignFrobNormsLowestLevel(...): "
		    "Matrix must be have equal number of rows "
		    "and cols to be symmetric");
    if (!A.is_zero()) {
      if ( this->is_zero() )
	this->allocate();
      assert( this->nElements() == A.nElements() );
      for (int col = 1; col < this->ncols(); col++)
	for (int row = 0; row < col; row++)
	  (*this)(row, col).assignFrobNormsLowestLevel( A(row,col) );
      for (int rc = 0; rc < this->ncols(); rc++)
	(*this)(rc, rc).syAssignFrobNormsLowestLevel( A(rc,rc) );
    }
    else 
      this->clear();    
  }

  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::assignDiffFrobNormsLowestLevel
    ( Matrix<Treal, Matrix<Treal, Telement> > const & A,
      Matrix<Treal, Matrix<Treal, Telement> > const & B ) {
    if ( A.is_zero() && B.is_zero() ) {
      // Both A and B are zero
      this->clear();
      return;
    }
    // At least one of A and B is nonzero
    if ( this->is_zero() )
      this->allocate();
    if ( !A.is_zero() && !B.is_zero() ) {
      // Both are nonzero
      assert( this->nElements() == A.nElements() );
      assert( this->nElements() == B.nElements() );
      for (int ind = 0; ind < this->nElements(); ind++)
	this->elements[ind].assignDiffFrobNormsLowestLevel( A[ind], B[ind] );
      return;
    }
    if ( !A.is_zero() ) {
      // A is nonzero
      this->assignFrobNormsLowestLevel( A );
      return;
    }
    if ( !B.is_zero() ) {
      // B is nonzero
      this->assignFrobNormsLowestLevel( B );
      return;
    }
  }
  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::syAssignDiffFrobNormsLowestLevel
    ( Matrix<Treal, Matrix<Treal, Telement> > const & A,
      Matrix<Treal, Matrix<Treal, Telement> > const & B ) {
    if ( A.is_zero() && B.is_zero() ) {
      // Both A and B are zero
      this->clear();
      return;
    }
    // At least one of A and B is nonzero
    if ( this->is_zero() )
      this->allocate();
    if ( !A.is_zero() && !B.is_zero() ) {
      // Both are nonzero
      assert( this->nElements() == A.nElements() );
      assert( this->nElements() == B.nElements() );
      for (int col = 1; col < this->ncols(); col++)
	for (int row = 0; row < col; row++)
	  (*this)(row, col).assignDiffFrobNormsLowestLevel( A(row,col), B(row,col) );
      for (int rc = 0; rc < this->ncols(); rc++)
	(*this)(rc, rc).syAssignDiffFrobNormsLowestLevel( A(rc,rc), B(rc,rc) );
      return;
    }
    if ( !A.is_zero() ) {
      // A is nonzero
      this->syAssignFrobNormsLowestLevel( A );
      return;
    }
    if ( !B.is_zero() ) {
      // B is nonzero
      this->syAssignFrobNormsLowestLevel( B );
      return;
    }
  } 



  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::
    truncateAccordingToSparsityPattern( Matrix<Treal, Matrix<Treal, Telement> > & A ) const {
    if ( this->is_zero() ) {
      A.clear();
    }
    else {
      assert( !A.is_zero() );
      assert( this->nElements() == A.nElements() );
      for (int ind = 0; ind < this->nElements(); ind++)
	this->elements[ind].truncateAccordingToSparsityPattern( A[ind] );
    }
  }  
  


  /********** End of help functions for thresholding */
  
  /* C = beta * C + alpha * A * B where only the upper triangle of C is */
  /* referenced and updated */
  template<class Treal, class Telement>
    void Matrix<Treal, Telement>::
    gemm_upper_tr_only(const bool tA, const bool tB, const Treal alpha, 
		       const Matrix<Treal, Telement>& A, 
		       const Matrix<Treal, Telement>& B, 
		       const Treal beta, 
		       Matrix<Treal, Telement>& C) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (C.is_empty()) {
      assert(beta == 0);
      if (tA)
	C.resetRows(A.cols);
      else
	C.resetRows(A.rows);
      if (tB)
	C.resetCols(B.rows);
      else
	C.resetCols(B.cols);
    }      
    if ((A.is_zero() || B.is_zero() || alpha == 0) && 
	(C.is_zero() || beta == 0))
      C = 0;
    else {
      Treal beta_tmp = beta;
      if (C.is_zero()) {
	C.allocate();
	beta_tmp = 0;
      }
      if (!A.is_zero() && !B.is_zero() && alpha != 0) {
	if (!tA && !tB) 
	  if (A.ncols() == B.nrows() && 
	      A.nrows() == C.nrows() &&
	      B.ncols() == C.ncols()) {
	    for (int col = 0; col < C.ncols(); col++) {
	      Telement::gemm_upper_tr_only(tA, tB, alpha, 
					   A(col, 0), B(0, col), 
					   beta_tmp, 
					   C(col, col));
	      for(int cola = 1; cola < A.ncols(); cola++)
		Telement::gemm_upper_tr_only(tA, tB, alpha, 
					     A(col, cola), B(cola, col), 
					     1.0, 
					     C(col,col));
	      for (int row = 0; row < col; row++) {
		Telement::gemm(tA, tB, alpha, 
			       A(row, 0), B(0, col), 
			       beta_tmp, 
			       C(row,col));
		for(int cola = 1; cola < A.ncols(); cola++)
		  Telement::gemm(tA, tB, alpha, 
				 A(row, cola), B(cola, col), 
				 1.0, 
				 C(row,col));
	      }
	    }
	  }
	  else 
	    throw Failure("Matrix<class Treal, class Telement>::"
			  "gemm_upper_tr_only(bool, bool, Treal, "
			  "Matrix<Treal, Telement>, "
			  "Matrix<Treal, Telement>, "
			  "Treal, Matrix<Treal, Telement>): "
			  "Incorrect matrixdimensions for multiplication");
	else if (tA && !tB)
	  if (A.nrows() == B.nrows() &&
	      A.ncols() == C.nrows() &&
	      B.ncols() == C.ncols()) {
	    for (int col = 0; col < C.ncols(); col++) {
	      Telement::gemm_upper_tr_only(tA, tB, alpha,
					   A(0,col), B(0,col),
					   beta_tmp,
					   C(col,col));
	      for(int cola = 1; cola < A.nrows(); cola++)
		Telement::gemm_upper_tr_only(tA, tB, alpha, 
					     A(cola, col), B(cola, col), 
					     1.0, 
					     C(col, col));
	      for (int row = 0; row < col; row++) {
		Telement::gemm(tA, tB, alpha,
			       A(0,row), B(0,col),
			       beta_tmp,
			       C(row,col));
		for(int cola = 1; cola < A.nrows(); cola++)
		  Telement::gemm(tA, tB, alpha, 
				 A(cola, row), B(cola, col), 
				 1.0, 
				 C(row,col));
	      }
	    }
	  }
	  else
	    throw Failure("Matrix<class Treal, class Telement>::"
			  "gemm_upper_tr_only(bool, bool, Treal, "
			  "Matrix<Treal, Telement>, "
			  "Matrix<Treal, Telement>, "
			  "Treal, Matrix<Treal, Telement>): "
			  "Incorrect matrixdimensions for multiplication");
	else if (!tA && tB)
	  if (A.ncols() == B.ncols() && 
	      A.nrows() == C.nrows() &&
	      B.nrows() == C.ncols()) {
	    for (int col = 0; col < C.ncols(); col++) {
	      Telement::gemm_upper_tr_only(tA, tB, alpha, 
					   A(col, 0), B(col, 0), 
					   beta_tmp, 
					   C(col,col));
	      for(int cola = 1; cola < A.ncols(); cola++)
		Telement::gemm_upper_tr_only(tA, tB, alpha, 
					     A(col, cola), B(col, cola), 
					     1.0, 
					     C(col,col));
	      for (int row = 0; row < col; row++) {
		Telement::gemm(tA, tB, alpha, 
			       A(row, 0), B(col, 0), 
			       beta_tmp, 
			       C(row,col));
		for(int cola = 1; cola < A.ncols(); cola++)
		  Telement::gemm(tA, tB, alpha, 
				 A(row, cola), B(col, cola), 
				 1.0, 
				 C(row,col));
	      }
	    }
	  }
	  else 
	    throw Failure("Matrix<class Treal, class Telement>::"
			  "gemm_upper_tr_only(bool, bool, Treal, "
			  "Matrix<Treal, Telement>, "
			  "Matrix<Treal, Telement>, "
			  "Treal, Matrix<Treal, Telement>): "
			  "Incorrect matrixdimensions for multiplication");
	else if (tA && tB)
	  if (A.nrows() == B.ncols() && 
	      A.ncols() == C.nrows() &&
	      B.nrows() == C.ncols()) {
	    for (int col = 0; col < C.ncols(); col++) {
	      Telement::gemm_upper_tr_only(tA, tB, alpha, 
					   A(0, col), B(col, 0), 
					   beta_tmp, 
					   C(col,col));
	      for(int cola = 1; cola < A.nrows(); cola++)
		Telement::gemm_upper_tr_only(tA, tB, alpha, 
					     A(cola, col), B(col, cola), 
					     1.0, 
					     C(col,col));
	      for (int row = 0; row < col; row++) {
		Telement::gemm(tA, tB, alpha, 
			       A(0, row), B(col, 0), 
			       beta_tmp, 
			       C(row,col));
		for(int cola = 1; cola < A.nrows(); cola++)
		  Telement::gemm(tA, tB, alpha, 
				 A(cola, row), B(col, cola), 
				 1.0, 
				 C(row,col));
	      }
	    }
	  }
	  else 
	    throw Failure("Matrix<class Treal, class Telement>::"
			  "gemm_upper_tr_only(bool, bool, Treal, "
			  "Matrix<Treal, Telement>, "
			  "Matrix<Treal, Telement>, Treal, "
			  "Matrix<Treal, Telement>): "
			  "Incorrect matrixdimensions for multiplication");
	else throw Failure("Matrix<class Treal, class Telement>::"
			   "gemm_upper_tr_only(bool, bool, Treal, "
			   "Matrix<Treal, Telement>, "
			   "Matrix<Treal, Telement>, Treal, "
			   "Matrix<Treal, Telement>):"
			   "Very strange error!!");
      }    
      else
	C *= beta;
    }
  }

  /* A = alpha * A * Z or A = alpha * Z * A where A is symmetric, */
  /* Z is upper triangular and */
  /* only the upper triangle of the result is calculated */
  /* side == 'R' for A = alpha * A * Z */
  /* side == 'L' for A = alpha * Z * A */
  template<class Treal, class Telement>
    void Matrix<Treal, Telement>:: 
    sytr_upper_tr_only(char const side, const Treal alpha,
		       Matrix<Treal, Telement>& A,
		       const Matrix<Treal, Telement>& Z) {
    assert(!A.is_empty());
    assert(!Z.is_empty());
    if (alpha != 0 && !A.is_zero() && !Z.is_zero()) {
      if (A.nrows() == A.ncols() && 
	  Z.nrows() == Z.ncols() && 
	  A.nrows() == Z.nrows()) {
	if (side == 'R') {
	  /* Last column must be calculated first                       */
    	  for (int col = A.ncols() - 1; col >= 0; col--) {
	    // A(col, col) = alpha * A(col, col) * Z(col, col)
	    Telement::sytr_upper_tr_only(side, alpha, 
					 A(col, col), Z(col, col));
    	    for (int ind = 0; ind < col; ind++) {
	      // A(col, col) += alpha * A(ind, col)' * Z(ind, col)
	      Telement::gemm_upper_tr_only(true, false, alpha, A(ind, col), 
					   Z(ind, col), 1.0, A(col, col));
	    }
	    /* Last row must be calculated first                       */
	    for (int row = col - 1; row >= 0; row--) {
	      // A(row, col) = alpha * A(row, col) * Z(col, col);
	      Telement::trmm(side, 'U', false, 
			     alpha, Z(col, col), A(row, col));
	      // A(row, col) += alpha * A(row, row) * Z(row, col);
	      Telement::symm('L', 'U', alpha, A(row, row), Z(row, col), 
			     1.0, A(row, col));
	      for (int ind = 0; ind < row; ind++) {
		// A(row, col) += alpha * A(ind, row)' * Z(ind, col);
		Telement::gemm(true, false, alpha, A(ind, row), Z(ind, col),
			       1.0, A(row, col));
	      }
	      for (int ind = row + 1; ind < col; ind++) {
		// A(row, col) += alpha * A(row, ind) * Z(ind, col);
		Telement::gemm(false, false, alpha, A(row, ind), Z(ind, col),
			       1.0, A(row, col));
	      }
	    }
	  }
	}
	else { /* side == 'L' */
	  assert(side == 'L');
	  /* First column must be calculated first                       */
	  for (int col = 0; col < A.ncols(); col++) {
	    /* First row must be calculated first                       */
	    for (int row = 0; row < col; row++) {
	      // A(row, col) = alpha * Z(row, row) * A(row, col)
	      Telement::trmm(side, 'U', false, alpha, 
			     Z(row, row), A(row, col));
	      // A(row, col) += alpha * Z(row, col) * A(col, col)
	      Telement::symm('R', 'U', alpha, A(col, col), Z(row, col), 
			     1.0, A(row, col));
	      for (int ind = row + 1; ind < col; ind++)
		// A(row, col) += alpha * Z(row, ind) * A(ind, col)
		Telement::gemm(false, false, alpha, Z(row, ind), A(ind, col),
			       1.0, A(row, col));
	      for (int ind = col + 1; ind < A.nrows(); ind++)
		// A(row, col) += alpha * Z(row, ind) * A(col, ind)'
		Telement::gemm(false, true, alpha, Z(row, ind), A(col, ind),
			       1.0, A(row, col));
	    }
	    Telement::sytr_upper_tr_only(side, alpha, 
					 A(col, col), Z(col, col));
	    for (int ind = col + 1; ind < A.ncols(); ind++)
	      Telement::gemm_upper_tr_only(false, true, alpha, Z(col, ind), 
					   A(col, ind), 1.0, A(col, col));
	  }
	}
      }
      else 	
	throw Failure("Matrix<class Treal, class Telement>::"
		      "sytr_upper_tr_only: Incorrect matrix dimensions "
		      "for symmetric-triangular multiplication");
    }     
    else
      A = 0;
  }

  /* The result B is assumed to be symmetric, i.e. only upper triangle */
  /* calculated and hence only upper triangle of input B is used. */
  template<class Treal, class Telement>
    void Matrix<Treal, Telement>:: 
    trmm_upper_tr_only(const char side, const char uplo, 
		       const bool tA, const Treal alpha, 
		       const Matrix<Treal, Telement>& A, 
		       Matrix<Treal, Telement>& B) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (alpha != 0 && !A.is_zero() && !B.is_zero())
      if (((side == 'R' && B.ncols() == A.nrows()) || 
	   (side == 'L' && A.ncols() == B.nrows())) &&
	  A.nrows() == A.ncols())
	if (uplo == 'U')
	  if (!tA) { 
	    throw Failure("Matrix<Treal, class Telement>::"
			  "trmm_upper_tr_only : "
			  "not possible for upper triangular not transposed "
			  "matrices A since lower triangle of B is needed");
	  } /* end if (tA == false) */
	  else {
	    assert(tA == true);
	    if (side == 'R') {
	      /* First column must be calculated first                       */
	      for (int col = 0; col < B.ncols(); col++) {
		Telement::trmm_upper_tr_only(side, uplo, tA, alpha,
					     A(col,col), B(col,col));
		for (int ind = col + 1; ind < A.ncols(); ind++)
		  Telement::gemm_upper_tr_only(false, tA, alpha,
					       B(col,ind), A(col,ind),
					       1.0, B(col,col));
		for (int row = 0; row < col; row++) {
		  Telement::trmm(side, uplo, tA, alpha,
				 A(col,col), B(row,col));
		  for (int ind = col + 1; ind < A.ncols(); ind++)
		    Telement::gemm(false, tA, alpha,
				   B(row,ind), A(col,ind),
				   1.0, B(row,col));
		}
	      }
	    } /* end if (side == 'R')*/
	    else {
	      assert(side == 'L');
	      /* Last row must be calculated first                           */
	      for (int row = B.nrows() - 1; row >= 0; row--) {
		Telement::trmm_upper_tr_only(side, uplo, tA, alpha,
					     A(row,row), B(row,row));
		for (int ind = 0; ind < row; ind++)
		  Telement::gemm_upper_tr_only(tA, false, alpha,
					       A(ind,row), B(ind,row),
					       1.0, B(row,row));
		for (int col = row + 1; col < B.ncols(); col++) {
		  Telement::trmm(side, uplo, tA, alpha,
				 A(row,row), B(row,col));
		  for (int ind = 0; ind < row; ind++)
		    Telement::gemm(tA, false, alpha,
				   A(ind,row), B(ind,col),
				   1.0, B(row,col));
		}
	      }
	    } /* end else (side == 'L')*/
	  } /* end else (tA == true)*/
	else
	  throw Failure("Matrix<class Treal, class Telement>::"
			"trmm_upper_tr_only not implemented for lower "
			"triangular matrices");
      else
	throw Failure("Matrix<class Treal, class Telement>::"
		      "trmm_upper_tr_only: Incorrect matrix dimensions "
		      "for multiplication");
    else
      B = 0;
  }

  /* A = Z' * A * Z or A = Z * A * Z' */
  /* where A is symmetric and Z is (nonunit) upper triangular */
  /* side == 'R' for A = Z' * A * Z */
  /* side == 'L' for A = Z * A * Z' */
  template<class Treal, class Telement>
    void Matrix<Treal, Telement>:: 
    trsytriplemm(char const side, 
		 const Matrix<Treal, Telement>& Z,
		 Matrix<Treal, Telement>& A) {
    if (side == 'R') {
      Matrix<Treal, Telement>::
	sytr_upper_tr_only('R', 1.0, A, Z);
      Matrix<Treal, Telement>::
	trmm_upper_tr_only('L', 'U', true, 1.0, Z, A);
    }
    else {
      assert(side == 'L');
      Matrix<Treal, Telement>::
	sytr_upper_tr_only('L', 1.0, A, Z);
      Matrix<Treal, Telement>::
	trmm_upper_tr_only('R', 'U', true, 1.0, Z, A);
    }
  }

  template<class Treal, class Telement>
    Treal Matrix<Treal, Telement>::
    frob_squared_thresh(Treal const threshold,
			Matrix<Treal, Telement> * ErrorMatrix) {
    assert(!this->is_empty());
    if (ErrorMatrix && ErrorMatrix->is_empty()) {
      ErrorMatrix->resetRows(this->rows);
      ErrorMatrix->resetCols(this->cols);
    }
    assert(threshold >= (Treal)0.0);
    if (threshold == (Treal)0.0)
      return 0;
    
    std::vector<Treal> frobsq_norms;
    getFrobSqLowestLevel(frobsq_norms);
    /* Sort in ascending order */
    std::sort(frobsq_norms.begin(), frobsq_norms.end());
    int lastRemoved = -1;
    Treal frobsqSum = 0;
    int nnzBlocks = frobsq_norms.size();
    while(lastRemoved < nnzBlocks - 1 && frobsqSum < threshold) {
      ++lastRemoved;
      frobsqSum += frobsq_norms[lastRemoved];
    }

    /* Check if entire matrix will be removed */
    if (lastRemoved == nnzBlocks - 1 && frobsqSum < threshold) {
      if (ErrorMatrix) 
	Matrix<Treal, Telement>::swap(*this, *ErrorMatrix);
      else
	*this = 0;
    }
    else {
      frobsqSum -= frobsq_norms[lastRemoved];
      --lastRemoved;
      while(lastRemoved > -1 && frobsq_norms[lastRemoved] == 
	    frobsq_norms[lastRemoved + 1]) {
	frobsqSum -= frobsq_norms[lastRemoved];
	--lastRemoved;
      }
      if (lastRemoved > -1) {
	Treal threshLowestLevel = 
	  (frobsq_norms[lastRemoved + 1] + 
	   frobsq_norms[lastRemoved]) / 2;
	this->frobThreshLowestLevel(threshLowestLevel, ErrorMatrix);
      }
    }
    return frobsqSum;
  }

  template<class Treal, class Telement> 
    void Matrix<Treal, Telement>::
    syInch(const Matrix<Treal, Telement>& A,
	   Matrix<Treal, Telement>& Z,
	   const Treal threshold, const side looking,
	   const inchversion version) {
    assert(!A.is_empty());
    if (A.nrows() != A.ncols()) 
      throw Failure("Matrix<Treal, Telement>::sy_inch: "
		    "Matrix must be quadratic!");
    if (A.is_zero())
      throw Failure("Matrix<Treal>::sy_inch: Matrix is zero! "
		    "Not possible to compute inverse cholesky.");
    if (version == stable) /* STABILIZED */
      throw Failure("Matrix<Treal>::sy_inch: Only unstable "
		    "version of sy_inch implemented.");
    Treal myThresh = 0;
    if (threshold != 0)
      myThresh = threshold / (Z.ncols() * Z.nrows());
    Z.resetRows(A.rows);
    Z.resetCols(A.cols);
    Z.allocate();
    
    if (looking == left) {  /* LEFT-LOOKING INCH */
      if (threshold != 0)
	throw Failure("Matrix<Treal, Telement>::syInch: "
		      "Thresholding not implemented for left-looking inch.");
      /* Left and unstable */
      Telement::syInch(A(0,0), Z(0,0), threshold, looking, version);
      Telement Ptmp;//, tmp;
      for (int i = 1; i < Z.ncols(); i++) {
	for (int j = 0; j < i; j++) {	
	  /* Z(k,i) is nonzero for k = 0, 1, ...,j - 2, j - 1, i         */
	  /* and Z(i,i) = I                       (yes it is i ^)        */
	  Ptmp = A(j,i);                    /* (Z(i,i) == I)             */
	  for (int k = 0; k < j; k++)       /* Ptmp = A(j,:) * Z(:,i)  */
	    Telement::gemm(true, false, 1.0, /* SYMMETRY USED */
			   A(k,j), Z(k,i), 1.0, Ptmp);
	  Telement::trmm('L', 'U', true, 1.0, Z(j,j), Ptmp);
	
	  for (int k = 0; k < j; k++)       /* Z(:,i) -= Z(:,j) * Ptmp   */
	    Telement::gemm(false, false, -1.0,
			   Z(k,j), Ptmp, 1.0, Z(k,i));
	  /* Z(j,j) is triangular: */
	  Telement::trmm('L', 'U', false, -1.0, Z(j,j), Ptmp);
	  Telement::add(1.0, Ptmp, Z(j,i));
	}
	Ptmp = A(i,i); /* Z(i,i) == I */
	for (int k = 0; k < i; k++) /* SYMMETRY USED */
	  Telement::gemm_upper_tr_only(true, false, 1.0, 
				       A(k,i), Z(k,i), 1.0, Ptmp);
	/* Z(i,i) == I !!                                                */
	/* Z(:,i) *= INCH(Ptmp) */
	Telement::syInch(Ptmp, Z(i,i), threshold, looking, version);
	for (int k = 0; k < i; k++) {         
	  Telement::trmm('R', 'U', false, 1.0, Z(i,i), Z(k,i));
	  /* INCH(Ptmp) is upper triangular*/
	} 
      }
    } /* end if left-looking inch */
    else { /* RIGHT-LOOKING INCH */
      assert(looking == right);  /* right and unstable */
      Telement Ptmp;
      Treal newThresh = 0;
      if (myThresh != 0 && Z.ncols() > 1)
	newThresh = myThresh * 0.0001;
      else
	newThresh = myThresh;

      for (int i = 0; i < Z.ncols(); i++) {
	/* Ptmp =  A(i,:) * Z(:,i)   */
	Ptmp = A(i,i); /* Z(i,i) == I */
	for (int k = 0; k < i; k++)
	  Telement::gemm_upper_tr_only(true, false, 1.0, /* SYMMETRY USED */ 
				       A(k,i), Z(k,i), 1.0, Ptmp);
	
	/* Z(:,i) *= INCH(Ptmp)  */
	Telement::syInch(Ptmp, Z(i,i), newThresh, looking, version);
	for (int k = 0; k < i; k++) {         
	  Telement::trmm('R', 'U', false, 1.0, Z(i,i), Z(k,i));
	  /* INCH(Ptmp) is upper triangular */
	}
	  
	for (int j = i + 1; j < Z.ncols(); j++) {
	  /* Compute Ptmp = Z(i,i)^T * A(i,:) * Z(:,j)  */
	  /* Z(k,j) is nonzero for k = 0, 1, ...,i - 2, i - 1, j         */
	  /* and Z(j,j) = I                                              */
	  Ptmp = A(i,j);                  /* (Z(j,j) == I)             */
	  for (int k = 0; k < i; k++)       /* Ptmp = A(i,:) * Z(:,j)  */
	    Telement::gemm(true, false, 1.0, /* SYMMETRY USED */
			   A(k,i), Z(k,j), 1.0, Ptmp);
	  Telement::trmm('L', 'U', true, 1.0, Z(i,i), Ptmp);
	    
	  for (int k = 0; k < i; k++)       /* Z(:,j) -= Z(:,i) * Ptmp   */
	    Telement::gemm(false, false, -1.0,
			   Z(k,i), Ptmp, 1.0, Z(k,j));
	  /* Z(i,i) is triangular: */
	  Telement::trmm('L', 'U', false, -1.0, Z(i,i), Ptmp);
	  Telement::add(1.0, Ptmp, Z(i,j));
	} /* end for j */
	
	/* Truncation starts here */
	if (threshold != 0) {
	  for (int k = 0; k < i; k++) 
	    Z(k,i).frob_thresh(myThresh);
	}
      } /* end for i */
    } /* end else right-looking inch */
  }

  template<class Treal, class Telement>
    void Matrix<Treal, Telement>:: 
    gershgorin(Treal& lmin, Treal& lmax) const {
    assert(!this->is_empty());
    if (this->nScalarsRows() == this->nScalarsCols()) {
      int size = this->nScalarsCols();
      Treal* colsums = new Treal[size];
      Treal* diag    = new Treal[size];
      for (int ind = 0; ind < size; ind++)
	colsums[ind] = 0;
      this->add_abs_col_sums(colsums);
      this->get_diagonal(diag);
      Treal tmp1 = colsums[0] - template_blas_fabs(diag[0]);
      Treal tmp2;
      lmin = diag[0] - tmp1;
      lmax = diag[0] + tmp1;
      for (int col = 1; col < size; col++) {
	tmp1 = colsums[col] - template_blas_fabs(diag[col]);
	tmp2 = diag[col] - tmp1;
	lmin = (tmp2 < lmin ? tmp2 : lmin);
	tmp2 = diag[col] + tmp1;
	lmax = (tmp2 > lmax ? tmp2 : lmax);
      }
      delete[] diag;
      delete[] colsums;
    }
    else
      throw Failure("Matrix<Treal, Telement, int>::gershgorin(Treal, Treal): "
		    "Matrix must be quadratic");
  }


  template<class Treal, class Telement> 
    void Matrix<Treal, Telement>::
    add_abs_col_sums(Treal* sums) const {
    assert(sums);
    if (!this->is_zero()) {
      int offset = 0;
      for (int col = 0; col < this->ncols(); col++) {   
	for (int row = 0; row < this->nrows(); row++) {
	  (*this)(row,col).add_abs_col_sums(&sums[offset]);
	}
	offset += (*this)(0,col).nScalarsCols();
      }
    }  
  }

  template<class Treal, class Telement> 
    void Matrix<Treal, Telement>::
    get_diagonal(Treal* diag) const {
    assert(diag);
    assert(this->nrows() == this->ncols());
    if (this->is_zero()) 
      for (int rc = 0; rc < this->nScalarsCols(); rc++)
	diag[rc] = 0;
    else {
      int offset = 0;
      for (int rc = 0; rc < this->ncols(); rc++) {
	(*this)(rc,rc).get_diagonal(&diag[offset]);
	offset += (*this)(rc,rc).nScalarsCols();
      }
    }
  }

  template<class Treal, class Telement>
    size_t Matrix<Treal, Telement>::memory_usage() const {
    assert(!this->is_empty());
    size_t sum = 0; 
    if (this->is_zero())
      return (size_t)0;
    for (int ind = 0; ind < this->nElements(); ind++)
      sum += this->elements[ind].memory_usage();
    return sum;
  }

  template<class Treal, class Telement>
    size_t Matrix<Treal, Telement>::nnz() const {
    size_t sum = 0;
    if (!this->is_zero()) {
      for (int ind = 0; ind < this->nElements(); ind++)
	sum += this->elements[ind].nnz();
    }
    return sum;
  }
  template<class Treal, class Telement>
    size_t Matrix<Treal, Telement>::sy_nnz() const {
    size_t sum = 0;    
    if (!this->is_zero()) {
      /* Above diagonal */
      for (int col = 1; col < this->ncols(); col++)
	for (int row = 0; row < col; row++)
	  sum += (*this)(row, col).nnz();
      /* Below diagonal */
      sum *= 2;
      /* Diagonal */
      for (int rc = 0; rc < this->nrows(); rc++)
	sum += (*this)(rc,rc).sy_nnz();
    }
    return sum;
  }

  template<class Treal, class Telement>
    size_t Matrix<Treal, Telement>::sy_nvalues() const {
    size_t sum = 0;
    if (!this->is_zero()) {
      /* Above diagonal */
      for (int col = 1; col < this->ncols(); col++)
	for (int row = 0; row < col; row++)
	  sum += (*this)(row, col).nvalues();
      /* Diagonal */
      for (int rc = 0; rc < this->nrows(); rc++)
	sum += (*this)(rc,rc).sy_nvalues();
    }
    return sum;
  }





  /***************************************************************************/
  /***************************************************************************/
  /*           Specialization for Telement = Treal                           */
  /***************************************************************************/
  /***************************************************************************/

  template<class Treal>
    class Matrix<Treal>: public MatrixHierarchicBase<Treal> {
   
    public:
    typedef Vector<Treal, Treal> VectorType;
    friend class Vector<Treal, Treal>;
    
    Matrix()
      :MatrixHierarchicBase<Treal>() {
    }
      /*    Matrix(const Atomblock<Treal>& row_atoms, 
	    const Atomblock<Treal>& col_atoms,
	    bool z = true, int nr = 0, int nc = 0, char tr = 'N')
	    :MatrixHierarchicBase<Treal>(row_atoms, col_atoms, z, nr, nc,tr) {}
      */
      
    void allocate() {
      assert(!this->is_empty());
      assert(this->is_zero());
      this->elements = allocateElements<Treal>(this->nElements());
      for (int ind = 0; ind < this->nElements(); ++ind)
	this->elements[ind] = 0;
    }
    
    /* Full matrix assigns etc */
    void assignFromFull(std::vector<Treal> const & fullMat);
    void fullMatrix(std::vector<Treal> & fullMat) const; 
    void syFullMatrix(std::vector<Treal> & fullMat) const;
    void syUpTriFullMatrix(std::vector<Treal> & fullMat) const;
    
    /* Sparse matrix assigns etc */
    void assignFromSparse(std::vector<int> const & rowind, 
			  std::vector<int> const & colind, 
			  std::vector<Treal> const & values);
    void assignFromSparse(std::vector<int> const & rowind, 
			  std::vector<int> const & colind, 
			  std::vector<Treal> const & values,
			  std::vector<int> const & indexes);
    
    /* Adds values (+=) to elements */
    void addValues(std::vector<int> const & rowind, 
		   std::vector<int> const & colind, 
		   std::vector<Treal> const & values);
    void addValues(std::vector<int> const & rowind, 
		   std::vector<int> const & colind, 
		   std::vector<Treal> const & values,
		   std::vector<int> const & indexes);
    
    void syAssignFromSparse(std::vector<int> const & rowind, 
			    std::vector<int> const & colind, 
			    std::vector<Treal> const & values);
    
    void syAddValues(std::vector<int> const & rowind, 
		     std::vector<int> const & colind, 
		     std::vector<Treal> const & values);
    
    void getValues(std::vector<int> const & rowind, 
		   std::vector<int> const & colind, 
		   std::vector<Treal> & values) const;
    void getValues(std::vector<int> const & rowind, 
		   std::vector<int> const & colind, 
		   std::vector<Treal> & values,
		   std::vector<int> const & indexes) const;
    void syGetValues(std::vector<int> const & rowind, 
		     std::vector<int> const & colind, 
		     std::vector<Treal> & values) const;
    
    void getAllValues(std::vector<int> & rowind, 
		     std::vector<int> & colind, 
		     std::vector<Treal> & values) const;
    void syGetAllValues(std::vector<int> & rowind, 
			std::vector<int> & colind, 
			std::vector<Treal> & values) const;

    Matrix<Treal>& 
      operator=(const Matrix<Treal>& mat) {
      MatrixHierarchicBase<Treal>::operator=(mat);
      return *this;
    } 

    void clear(); /**< Set matrix to zero and delete all arrays */
    ~Matrix() {
      clear();
    }

    void writeToFile(std::ofstream & file) const;
    void readFromFile(std::ifstream & file);

    void random();
    void syRandom();
    void randomZeroStructure(Treal probabilityBeingZero);
    void syRandomZeroStructure(Treal probabilityBeingZero);
    template<typename TRule>
      void setElementsByRule(TRule & rule);
    template<typename TRule>
      void sySetElementsByRule(TRule & rule);
    

    void addIdentity(Treal alpha);        /* C += alpha * I */

    static void transpose(Matrix<Treal> const & A,
			  Matrix<Treal> & AT);

    void symToNosym();
    void nosymToSym();
    
    /* Set matrix to zero (k = 0) or identity (k = 1)                      */
    Matrix<Treal>& operator=(int const k);

    Matrix<Treal>& operator*=(const Treal alpha);
    
    static void gemm(const bool tA, const bool tB, const Treal alpha, 
		     const Matrix<Treal>& A, 
		     const Matrix<Treal>& B, 
		     const Treal beta, 
		     Matrix<Treal>& C);
    static void symm(const char side, const char uplo, const Treal alpha, 
		     const Matrix<Treal>& A, 
		     const Matrix<Treal>& B, 
		     const Treal beta, 
		     Matrix<Treal>& C);
    static void syrk(const char uplo, const bool tA, const Treal alpha, 
		     const Matrix<Treal>& A, 
		     const Treal beta, 
		     Matrix<Treal>& C);
    /* C = beta * C + alpha * A * A where A and C are symmetric */
    static void sysq(const char uplo, const Treal alpha, 
		     const Matrix<Treal>& A, 
		     const Treal beta, 
		     Matrix<Treal>& C);
    /* C = alpha * A * B + beta * C where A and B are symmetric */
    static void ssmm(const Treal alpha, 
		     const Matrix<Treal>& A, 
		     const Matrix<Treal>& B, 
		     const Treal beta, 
		     Matrix<Treal>& C);
    /* C = alpha * A * B + beta * C where A and B are symmetric
     * and only the upper triangle of C is computed.
     */
    static void ssmm_upper_tr_only(const Treal alpha, 
				   const Matrix<Treal>& A, 
				   const Matrix<Treal>& B, 
				   const Treal beta, 
				   Matrix<Treal>& C);

    static void trmm(const char side, const char uplo, const bool tA, 
		     const Treal alpha, 
		     const Matrix<Treal>& A, 
		     Matrix<Treal>& B);

    /* Frobenius norms */
    inline Treal frob() const {return template_blas_sqrt(frobSquared());} 
    Treal frobSquared() const;
    inline Treal syFrob() const {
      return template_blas_sqrt(this->syFrobSquared());
    }
    Treal syFrobSquared() const;
    
    inline static Treal frobDiff
      (const Matrix<Treal>& A,
       const Matrix<Treal>& B) {
      return template_blas_sqrt(frobSquaredDiff(A, B));
    }
    static Treal frobSquaredDiff
      (const Matrix<Treal>& A,
       const Matrix<Treal>& B);
    
    inline static Treal syFrobDiff
      (const Matrix<Treal>& A,
       const Matrix<Treal>& B) {
      return template_blas_sqrt(syFrobSquaredDiff(A, B));	
    }
    static Treal syFrobSquaredDiff
      (const Matrix<Treal>& A,
       const Matrix<Treal>& B);      

    Treal trace() const;
    static Treal trace_ab(const Matrix<Treal>& A,
			  const Matrix<Treal>& B); 
    static Treal trace_aTb(const Matrix<Treal>& A,
			   const Matrix<Treal>& B); 
    static Treal sy_trace_ab(const Matrix<Treal>& A,
			     const Matrix<Treal>& B); 

    static void add(const Treal alpha,   /* B += alpha * A */
		    const Matrix<Treal>& A, 
		    Matrix<Treal>& B);
    void assign(Treal const  alpha,   /* *this = alpha * A */
		Matrix<Treal> const & A);
    

    /********** Help functions for thresholding */
    //    int getnnzBlocksLowestLevel() const;
    void getFrobSqLowestLevel(std::vector<Treal> & frobsq) const;
    void frobThreshLowestLevel
      (Treal const threshold, Matrix<Treal> * ErrorMatrix);

    void getFrobSqElementLevel(std::vector<Treal> & frobsq) const;
    void frobThreshElementLevel
      (Treal const threshold, Matrix<Treal> * ErrorMatrix);
    

#if 0
    inline void frobThreshLowestLevel
      (Treal const threshold, 
       Matrix<Treal> * ErrorMatrix) {
      bool a,b;
      frobThreshLowestLevel(threshold, ErrorMatrix, a, b);
    }
#endif

    void assignFrobNormsLowestLevel
      ( Matrix<Treal, Matrix<Treal> > const & A ) {
      if (!A.is_zero()) {
	if ( this->is_zero() )
	  this->allocate();
	assert( this->nElements() == A.nElements() );
	for (int ind = 0; ind < this->nElements(); ind++)
	  this->elements[ind] = A[ind].frob();
      }
      else
	this->clear();
    } 

    void syAssignFrobNormsLowestLevel( Matrix<Treal, Matrix<Treal> > const & A ) {
      if (!A.is_zero()) {
	if ( this->is_zero() )
	  this->allocate();
	assert( this->nElements() == A.nElements() );
	for (int col = 1; col < this->ncols(); col++)
	  for (int row = 0; row < col; row++)
	    (*this)(row,col) = A(row,col).frob();
	for (int rc = 0; rc < this->nrows(); rc++)
	  (*this)(rc,rc) = A(rc,rc).syFrob();
      }
      else
	this->clear();
    }

    void assignDiffFrobNormsLowestLevel( Matrix<Treal, Matrix<Treal> > const & A,
					 Matrix<Treal, Matrix<Treal> > const & B ) {
      if ( A.is_zero() && B.is_zero() ) {
	// Both A and B are zero
	this->clear();
	return;
      }
      // At least one of A and B is nonzero
      if ( this->is_zero() )
	this->allocate();
      if ( !A.is_zero() && !B.is_zero() ) {
	// Both are nonzero
	assert( this->nElements() == A.nElements() );
	assert( this->nElements() == B.nElements() );
	for (int ind = 0; ind < this->nElements(); ind++) {
	  Matrix<Treal> Diff(A[ind]);
	  Matrix<Treal>::add( -1.0, B[ind], Diff );
	  this->elements[ind] = Diff.frob();
	}
	return;
      }
      if ( !A.is_zero() ) {
	// A is nonzero
	this->assignFrobNormsLowestLevel( A );
	return;
      }
      if ( !B.is_zero() ) {
	// B is nonzero
	this->assignFrobNormsLowestLevel( B );
	return;
      }
    } 
    void syAssignDiffFrobNormsLowestLevel( Matrix<Treal, Matrix<Treal> > const & A,
					   Matrix<Treal, Matrix<Treal> > const & B ) {
      if ( A.is_zero() && B.is_zero() ) {
	// Both A and B are zero
	this->clear();
	return;
      }
      // At least one of A and B is nonzero
      if ( this->is_zero() )
	this->allocate();
      if ( !A.is_zero() && !B.is_zero() ) {
	// Both are nonzero
	assert( this->nElements() == A.nElements() );
	assert( this->nElements() == B.nElements() );
	for (int col = 1; col < this->ncols(); col++)
	  for (int row = 0; row < col; row++) {
	    Matrix<Treal> Diff(A(row,col));
	    Matrix<Treal>::add( -1.0, B(row,col), Diff );
	    (*this)(row, col) = Diff.frob();
	  }
	for (int rc = 0; rc < this->ncols(); rc++) {
	  Matrix<Treal> Diff( A(rc,rc) );
	  Matrix<Treal>::add( -1.0, B(rc,rc), Diff );
	  (*this)(rc, rc) = Diff.syFrob();
	}
	return;
      }
      if ( !A.is_zero() ) {
	// A is nonzero
	this->syAssignFrobNormsLowestLevel( A );
	return;
      }
      if ( !B.is_zero() ) {
	// B is nonzero
	this->syAssignFrobNormsLowestLevel( B );
	return;
      }
    }


    void truncateAccordingToSparsityPattern( Matrix<Treal, Matrix<Treal> > & A ) const {
      if ( this->is_zero() ) 
	A.clear();
      else {
	assert( !A.is_zero() );
	assert( this->nElements() == A.nElements() );
	for (int ind = 0; ind < this->nElements(); ind++) 
	  if (this->elements[ind] == 0)
	    A[ind].clear();
      }	
    }
    

    /********** End of help functions for thresholding */

    static void gemm_upper_tr_only(const bool tA, const bool tB, 
				   const Treal alpha, 
				   const Matrix<Treal>& A, 
				   const Matrix<Treal>& B, 
				   const Treal beta, 
				   Matrix<Treal>& C);
    static void sytr_upper_tr_only(char const side, const Treal alpha,
				   Matrix<Treal>& A,
				   const Matrix<Treal>& Z);
    static void trmm_upper_tr_only(const char side, const char uplo, 
				   const bool tA, const Treal alpha, 
				   const Matrix<Treal>& A, 
				   Matrix<Treal>& B);
    static void trsytriplemm(char const side, 
			     const Matrix<Treal>& Z,
			     Matrix<Treal>& A);

    inline Treal frob_thresh(Treal const threshold,
			     Matrix<Treal> * ErrorMatrix = 0) {     
      return template_blas_sqrt
	(frob_squared_thresh(threshold * threshold, ErrorMatrix));
    }   
    /* Returns the Frobenius norm of the introduced error */
    
    Treal frob_squared_thresh(Treal const threshold,
			      Matrix<Treal> * ErrorMatrix = 0);
    

    static void inch(const Matrix<Treal>& A,
		     Matrix<Treal>& Z,
		     const Treal threshold = 0,
		     const side looking = left,
		     const inchversion version = unstable);
    static void syInch(const Matrix<Treal>& A,
			Matrix<Treal>& Z,
			const Treal threshold = 0, 
			const side looking = left,
			const inchversion version = unstable) {
      inch(A, Z, threshold, looking, version);
    }

    void gershgorin(Treal& lmin, Treal& lmax) const;
    void sy_gershgorin(Treal& lmin, Treal& lmax) const {
      Matrix<Treal> tmp = (*this);
      tmp.symToNosym();
      tmp.gershgorin(lmin, lmax);
      return;
    }
    void add_abs_col_sums(Treal* abscolsums) const;
    void get_diagonal(Treal* diag) const; /* Copy diagonal to the diag array */
    
    inline size_t memory_usage() const { /* Returns memory usage in bytes */
      assert(!this->is_empty());
      if (this->is_zero())
	return (size_t)0;
      else
	return (size_t)this->nElements() * sizeof(Treal);
    }

    inline size_t nnz() const {
      if (this->is_zero())
	return 0;
      else
	return this->nElements();
    } /**< Returns number of nonzeros in matrix. */
    inline size_t sy_nnz() const {
      if (this->is_zero())
	return 0;
      else
	return this->nElements();
    } /**< Returns number of nonzeros in matrix
       *   including lower triangle elements.
       */

    inline size_t nvalues() const {
      return nnz();
    } /**< Returns number of stored values in matrix.
       *   Returns same number as nnz()
       */
    size_t sy_nvalues() const {
      assert(this->nScalarsRows() == this->nScalarsCols());
      int n = this->nrows();
      return this->is_zero() ? 0 : n * (n + 1) / 2; 
    } /**< Returns number of stored values in matrix.
       *   Lower triangle is not included.
       *   Different from sy_nnz().
       */
    
    template<class Top> 
      Treal syAccumulateWith(Top & op) {
      Treal res = 0;
      if (!this->is_zero()) {
	int rowOffset = this->rows.getOffset();
	int colOffset = this->cols.getOffset();
	for (int col = 0; col < this->ncols(); col++) {
	  for (int row = 0; row < col; row++) {
	    res += 2*op.accumulate((*this)(row, col),
				   rowOffset + row,
				   colOffset + col);
	  }
	  res += op.accumulate((*this)(col, col),
			       rowOffset + col,
			       colOffset + col);
	}
      }
      return res;
    }
    template<class Top>
      Treal geAccumulateWith(Top & op) {
      Treal res = 0;
      if (!this->is_zero()) {
	int rowOffset = this->rows.getOffset();
	int colOffset = this->cols.getOffset();
	for (int col = 0; col < this->ncols(); col++)
	  for (int row = 0; row < this->nrows(); row++)
	    res += op.accumulate((*this)(row, col),
				 rowOffset + row,
				 colOffset + col);
      }
      return res;
    }
    
    static inline unsigned int level() {return 0;}
    
    Treal maxAbsValue() const {
      if (this->is_zero())
	return 0;
      else {
	Treal maxAbsGlobal = 0;
	Treal maxAbsLocal  = 0;
	for (int ind = 0; ind < this->nElements(); ++ind) {
	  maxAbsLocal = template_blas_fabs(this->elements[ind]);
	  maxAbsGlobal = maxAbsGlobal > maxAbsLocal ? 
	    maxAbsGlobal : maxAbsLocal; 
	} /* end for */
	return maxAbsGlobal;
      } 
    }
    
    /* New routines above */    
    
#if 0 /* OLD ROUTINES */   


#if 0
    inline Matrix<Treal>& operator=(const Matrix<Treal>& mat) {
      this->MatrixHierarchicBase<Treal>::operator=(mat);
      std::cout<<"operator="<<std::endl;
    }
#endif




    





    void trim_memory_usage();
#if 1
    void write_to_buffer_count(int& zb_length, int& vb_length) const;
    void write_to_buffer(int* zerobuf, const int zb_length, 
			 Treal* valuebuf, const int vb_length,
			 int& zb_index, int& vb_index) const;
    void read_from_buffer(int* zerobuf, const int zb_length, 
			  Treal* valuebuf, const int vb_length,
			  int& zb_index, int& vb_index);
#endif




    
    
    /* continue here */

    

    



    
    inline bool lowestLevel() const {return true;}
    //    inline unsigned int level() const {return 0;}

#endif     /* END OF OLD ROUTINES */   
  protected:
  private:
    static const Treal ZERO;
    static const Treal ONE;
  }; /* end class specialization Matrix<Treal> */

  template<class Treal>
    const Treal Matrix<Treal>::ZERO = 0;
  template<class Treal>
    const Treal Matrix<Treal>::ONE = 1;

#if 1
  /* Full matrix assigns etc */
  template<class Treal> 
    void Matrix<Treal>::
    assignFromFull(std::vector<Treal> const & fullMat) {
    int nTotalRows = this->rows.getNTotalScalars();
    int nTotalCols = this->cols.getNTotalScalars();
    assert((int)fullMat.size() == nTotalRows * nTotalCols);
    int rowOffset = this->rows.getOffset();
    int colOffset = this->cols.getOffset();
    if (this->is_zero())
      allocate();
    for (int col = 0; col < this->ncols(); col++) 
      for (int row = 0; row < this->nrows(); row++) 
	(*this)(row, col) = 
	  fullMat[(colOffset + col) * nTotalRows + (rowOffset + row)];
  }

  template<class Treal> 
    void Matrix<Treal>::
    fullMatrix(std::vector<Treal> & fullMat) const {
    int nTotalRows = this->rows.getNTotalScalars();
    int nTotalCols = this->cols.getNTotalScalars();
    fullMat.resize(nTotalRows * nTotalCols);
    int rowOffset = this->rows.getOffset();
    int colOffset = this->cols.getOffset();
    if (this->is_zero()) {
      for (int col = 0; col < this->nScalarsCols(); col++)
	for (int row = 0; row < this->nScalarsRows(); row++) 
	  fullMat[(colOffset + col) * nTotalRows + (rowOffset + row)] = 0;
    } 
    else {
      for (int col = 0; col < this->ncols(); col++) 
	for (int row = 0; row < this->nrows(); row++) 
	  fullMat[(colOffset + col) * nTotalRows + (rowOffset + row)] = 
	    (*this)(row, col);
    }
  }

  template<class Treal> 
    void Matrix<Treal>::
    syFullMatrix(std::vector<Treal> & fullMat) const {
    int nTotalRows = this->rows.getNTotalScalars();
    int nTotalCols = this->cols.getNTotalScalars();
    fullMat.resize(nTotalRows * nTotalCols);
    int rowOffset = this->rows.getOffset();
    int colOffset = this->cols.getOffset();
    if (this->is_zero()) {
      for (int col = 0; col < this->nScalarsCols(); col++)
	for (int row = 0; row <= col; row++) {
	  fullMat[(colOffset + col) * nTotalRows + (rowOffset + row)] = 0;
	  fullMat[(rowOffset + row) * nTotalRows + (colOffset + col)] = 0;
	}
    } 
    else {
      for (int col = 0; col < this->ncols(); col++) 
	for (int row = 0; row <= col; row++) {
	  fullMat[(colOffset + col) * nTotalRows + (rowOffset + row)] = 
	    (*this)(row, col);
	  fullMat[(rowOffset + row) * nTotalRows + (colOffset + col)] = 
	    (*this)(row, col);
	}
    }
  }
  
  template<class Treal> 
    void Matrix<Treal>::
    syUpTriFullMatrix(std::vector<Treal> & fullMat) const {
    int nTotalRows = this->rows.getNTotalScalars();
    int nTotalCols = this->cols.getNTotalScalars();
    fullMat.resize(nTotalRows * nTotalCols);
    int rowOffset = this->rows.getOffset();
    int colOffset = this->cols.getOffset();
    if (this->is_zero()) {
      for (int col = 0; col < this->nScalarsCols(); col++)
	for (int row = 0; row <= this->nScalarsRows(); row++) {
	  fullMat[(colOffset + col) * nTotalRows + (rowOffset + row)] = 0;
	  fullMat[(rowOffset + row) * nTotalRows + (colOffset + col)] = 0;
	}
    } 
    else {
      for (int col = 0; col < this->ncols(); col++) 
	for (int row = 0; row < this->nrows(); row++) {
	  fullMat[(colOffset + col) * nTotalRows + (rowOffset + row)] = 
	    (*this)(row, col);
	  fullMat[(rowOffset + row) * nTotalRows + (colOffset + col)] = 
	    (*this)(row, col);
	}
    }
  }
  
#endif
  
  template<class Treal>
    void Matrix<Treal>::
    assignFromSparse(std::vector<int> const & rowind, 
		     std::vector<int> const & colind, 
		     std::vector<Treal> const & values) {
    assert(rowind.size() == colind.size() &&
	   rowind.size() == values.size());
    std::vector<int> indexes(values.size());
    for (int ind = 0; ind < values.size(); ++ind) {
      indexes[ind] = ind;
    }
    assignFromSparse(rowind, colind, values, indexes);
  }

  template<class Treal>
    void Matrix<Treal>::
    assignFromSparse(std::vector<int> const & rowind, 
		     std::vector<int> const & colind, 
		     std::vector<Treal> const & values,
		     std::vector<int> const & indexes) {
    if (indexes.empty()) {
      this->clear();
      return;
    }
    if (this->is_zero())
      allocate();
    for (int ind = 0; ind < this->nElements(); ++ind)
      this->elements[ind] = 0;
    std::vector<int>::const_iterator it;
    for ( it = indexes.begin(); it < indexes.end(); it++ ) {
      (*this)(this->rows.whichBlock( rowind[*it] ),
	      this->cols.whichBlock( colind[*it] )) = values[*it];
    }
  }

  
  template<class Treal>
    void Matrix<Treal>::
    addValues(std::vector<int> const & rowind, 
	      std::vector<int> const & colind, 
	      std::vector<Treal> const & values) {
    assert(rowind.size() == colind.size() &&
	   rowind.size() == values.size());
    std::vector<int> indexes(values.size());
    for (int ind = 0; ind < values.size(); ++ind) {
      indexes[ind] = ind;
    }
    addValues(rowind, colind, values, indexes);
  }
  
  template<class Treal>
    void Matrix<Treal>::
    addValues(std::vector<int> const & rowind, 
	      std::vector<int> const & colind, 
	      std::vector<Treal> const & values,
	      std::vector<int> const & indexes) {
    if (indexes.empty())
      return;
    if (this->is_zero())
      allocate();
    std::vector<int>::const_iterator it;
    for ( it = indexes.begin(); it < indexes.end(); it++ ) {
      (*this)(this->rows.whichBlock( rowind[*it] ),
	      this->cols.whichBlock( colind[*it] )) += values[*it];
    }
  }
    
  template<class Treal>
    void Matrix<Treal>::
    syAssignFromSparse(std::vector<int> const & rowind, 
		       std::vector<int> const & colind, 
		       std::vector<Treal> const & values) {
    assert(rowind.size() == colind.size() &&
	   rowind.size() == values.size());
    bool upperTriangleOnly = true;
    for (int ind = 0; ind < values.size(); ++ind) {
      upperTriangleOnly = 
	upperTriangleOnly && colind[ind] >= rowind[ind];
    }
    if (!upperTriangleOnly)
      throw Failure("Matrix<Treal>::"
		    "syAddValues(std::vector<int> const &, "
		    "std::vector<int> const &, "
		    "std::vector<Treal> const &, int const): "
		    "Only upper triangle can contain elements when assigning "
		    "symmetric or triangular matrix ");
    assignFromSparse(rowind, colind, values);
  }
    
  template<class Treal>
    void Matrix<Treal>::
    syAddValues(std::vector<int> const & rowind, 
		std::vector<int> const & colind, 
		std::vector<Treal> const & values) {
    assert(rowind.size() == colind.size() &&
	   rowind.size() == values.size());
    bool upperTriangleOnly = true;
    for (int ind = 0; ind < values.size(); ++ind) {
      upperTriangleOnly = 
	upperTriangleOnly && colind[ind] >= rowind[ind];
    }
    if (!upperTriangleOnly)
      throw Failure("Matrix<Treal>::"
		    "syAddValues(std::vector<int> const &, "
		    "std::vector<int> const &, "
		    "std::vector<Treal> const &, int const): "
		    "Only upper triangle can contain elements when assigning "
		    "symmetric or triangular matrix ");
    addValues(rowind, colind, values);
  }
    
  template<class Treal>
    void Matrix<Treal>::
    getValues(std::vector<int> const & rowind, 
	      std::vector<int> const & colind, 
	      std::vector<Treal> & values) const {
    assert(rowind.size() == colind.size());
    values.resize(rowind.size(),0);
    std::vector<int> indexes(rowind.size());
    for (int ind = 0; ind < rowind.size(); ++ind) {
      indexes[ind] = ind;
    }
    getValues(rowind, colind, values, indexes);
  }
  
  template<class Treal>
    void Matrix<Treal>::
    getValues(std::vector<int> const & rowind, 
	      std::vector<int> const & colind, 
	      std::vector<Treal> & values,
	      std::vector<int> const & indexes) const {
    assert(!this->is_empty());
    if (indexes.empty())
      return;
    std::vector<int>::const_iterator it;
    if (this->is_zero()) {
      for ( it = indexes.begin(); it < indexes.end(); it++ )
	values[*it] = 0;
      return;
    }
    for ( it = indexes.begin(); it < indexes.end(); it++ ) 
      values[*it] = (*this)(this->rows.whichBlock( rowind[*it] ),
			    this->cols.whichBlock( colind[*it] ));
  }
  

  template<class Treal>
    void Matrix<Treal>::
    syGetValues(std::vector<int> const & rowind, 
		std::vector<int> const & colind, 
		std::vector<Treal> & values) const {
    assert(rowind.size() == colind.size());
    bool upperTriangleOnly = true;
    for (int ind = 0; ind < rowind.size(); ++ind) {
      upperTriangleOnly = 
	upperTriangleOnly && colind[ind] >= rowind[ind];
    }
    if (!upperTriangleOnly)
      throw Failure("Matrix<Treal>::"
		    "syGetValues(std::vector<int> const &, "
		    "std::vector<int> const &, "
		    "std::vector<Treal> const &, int const): "
		    "Only upper triangle when retrieving elements from "
		    "symmetric or triangular matrix ");
    getValues(rowind, colind, values);
  }

  template<class Treal>
    void Matrix<Treal>::
    getAllValues(std::vector<int> & rowind, 
		 std::vector<int> & colind, 
		 std::vector<Treal> & values) const {
    assert(rowind.size() == colind.size() &&
	   rowind.size() == values.size());
    if (!this->is_zero()) 
      for (int col = 0; col < this->ncols(); col++)
	for (int row = 0; row < this->nrows(); row++) {
 	  rowind.push_back(this->rows.getOffset() + row);
	  colind.push_back(this->cols.getOffset() + col);
	  values.push_back((*this)(row, col));
	}
  }
  
  template<class Treal>
    void Matrix<Treal>::
    syGetAllValues(std::vector<int> & rowind, 
		   std::vector<int> & colind, 
		   std::vector<Treal> & values) const {
    assert(rowind.size() == colind.size() &&
	   rowind.size() == values.size());
    if (!this->is_zero()) 
      for (int col = 0; col < this->ncols(); ++col) 
	for (int row = 0; row <= col; ++row) {
	  rowind.push_back(this->rows.getOffset() + row);
	  colind.push_back(this->cols.getOffset() + col);
	  values.push_back((*this)(row, col));
	}
  }

    
  template<class Treal>
    void Matrix<Treal>::clear() {
    freeElements(this->elements);
    this->elements = 0;
  }

  template<class Treal>
    void Matrix<Treal>:: 
    writeToFile(std::ofstream & file) const {
    int const ZERO = 0;
    int const ONE  = 1;
    if (this->is_zero()) {
      char * tmp = (char*)&ZERO;
      file.write(tmp,sizeof(int));
    }
    else {
      char * tmp = (char*)&ONE;
      file.write(tmp,sizeof(int));
      char * tmpel = (char*)this->elements;
      file.write(tmpel,sizeof(Treal) * this->nElements());
    }
  }

  template<class Treal>
    void Matrix<Treal>:: 
    readFromFile(std::ifstream & file) {
    int const ZERO = 0;
    int const ONE  = 1;
    char tmp[sizeof(int)];
    file.read(tmp, (std::ifstream::pos_type)sizeof(int));
    switch ((int)*tmp) {
    case ZERO:
      this->clear();
      break;
    case ONE:
      if (this->is_zero())
	allocate();
      file.read((char*)this->elements, sizeof(Treal) * this->nElements());
      break;
    default:
      throw Failure("Matrix<Treal>::" 
		    "readFromFile(std::ifstream & file):"
		    "File corruption, int value not 0 or 1");
    }
  }

  template<class Treal> 
    void Matrix<Treal>::random() {
    if (this->is_zero())
      allocate();
    for (int ind = 0; ind < this->nElements(); ind++)
      this->elements[ind] = rand() / (Treal)RAND_MAX;
  }
  
  template<class Treal> 
    void Matrix<Treal>::syRandom() {
    if (this->is_zero())
      allocate();
    /* Above diagonal */
    for (int col = 1; col < this->ncols(); col++)
      for (int row = 0; row < col; row++)
	(*this)(row, col) = rand() / (Treal)RAND_MAX;;
    /* Diagonal */
    for (int rc = 0; rc < this->nrows(); rc++)
      (*this)(rc,rc) = rand() / (Treal)RAND_MAX;;
  }
  
  template<class Treal> 
    void Matrix<Treal>::
    randomZeroStructure(Treal probabilityBeingZero) {
    if (!this->highestLevel() && 
	probabilityBeingZero > rand() / (Treal)RAND_MAX)
      this->clear();
    else 
      this->random();
  }

  template<class Treal> 
    void  Matrix<Treal>::
    syRandomZeroStructure(Treal probabilityBeingZero) {
    if (!this->highestLevel() && 
	probabilityBeingZero > rand() / (Treal)RAND_MAX)
      this->clear();
    else 
      this->syRandom();
  }

  template<class Treal> 
    template<typename TRule>
    void Matrix<Treal>::
    setElementsByRule(TRule & rule) {
    if (this->is_zero()) 
      allocate();
    for (int col = 0; col < this->ncols(); col++)
      for (int row = 0; row < this->nrows(); row++)
	(*this)(row,col) = rule.set(this->rows.getOffset() + row,
				      this->cols.getOffset() + col);
  }
  template<class Treal> 
    template<typename TRule>
    void Matrix<Treal>::
    sySetElementsByRule(TRule & rule) {
    if (this->is_zero()) 
      allocate(); 
    /* Upper triangle */
    for (int col = 0; col < this->ncols(); col++)
      for (int row = 0; row <= col; row++)
	(*this)(row,col) = rule.set(this->rows.getOffset() + row,
				      this->cols.getOffset() + col);
  }


  template<class Treal> 
    void Matrix<Treal>::
    addIdentity(Treal alpha) {
    if (this->is_empty())
      throw Failure("Matrix<Treal>::addIdentity(Treal): "
		    "Cannot add identity to empty matrix.");
    if (this->ncols() != this->nrows()) 
      throw Failure("Matrix<Treal, Telement>::addIdentity(Treal): "
		    "Matrix must be square to add identity");
    if (this->is_zero()) 
      allocate();
    for (int ind = 0; ind < this->ncols(); ind++)
      (*this)(ind,ind) += alpha;
  }

  template<class Treal> 
    void Matrix<Treal>::
    transpose(Matrix<Treal> const & A, Matrix<Treal> & AT) {
    if (A.is_zero()) { /* Condition also matches empty matrices. */
      AT.rows = A.cols;
      AT.cols = A.rows;
      freeElements(AT.elements);
      AT.elements = 0;
      return;
    }
    if (AT.is_zero() || (AT.nElements() != A.nElements())) {
      freeElements(AT.elements);
      AT.elements = allocateElements<Treal>(A.nElements());
    }
    AT.cols = A.rows;
    AT.rows = A.cols;
    for (int row = 0; row < AT.nrows(); row++)
      for (int col = 0; col < AT.ncols(); col++)
	AT(row,col) = A(col,row);
  }

  template<class Treal> 
    void Matrix<Treal>::
    symToNosym() {
    if (this->nrows() == this->ncols()) {
      if (!this->is_zero()) {
	/* Diagonal should be fine */
	/* Fix the lower triangle */
	for (int row = 1; row < this->nrows(); row++)
	  for (int col = 0; col < row; col++)
	    (*this)(row, col) = (*this)(col, row);
      }
    }
    else
      throw Failure("Matrix<Treal>::symToNosym(): "
		    "Only quadratic matrices can be symmetric");
  }

  template<class Treal> 
    void Matrix<Treal>::
    nosymToSym() {
    if (this->nrows() == this->ncols()) {
      if (!this->is_zero()) {
	/* Diagonal should be fine */
	/* Fix the lower triangle */
	for (int col = 0; col < this->ncols() - 1; col++)
	  for (int row = col + 1; row < this->nrows(); row++)
	    (*this)(row, col) = 0;
      }
    }
    else
      throw Failure("Matrix<Treal>::nosymToSym(): "
		    "Only quadratic matrices can be symmetric");
  }

  template<class Treal>
    Matrix<Treal>& 
    Matrix<Treal>::operator=(int const k) {
    switch (k) {
    case 0:
      this->clear();
      break;
    case 1:
      if (this->ncols() != this->nrows())
	throw Failure("Matrix<Treal>::operator=(int k = 1): "
		      "Matrix must be quadratic to "
		      "become identity matrix.");
      this->clear();
      this->allocate();
      for (int rc = 0; rc < this->ncols(); rc++) /*Set diagonal to identity*/
	(*this)(rc,rc) = 1;
      break;
    default:
      throw Failure
	("Matrix<Treal>::operator=(int k) only implemented for k = 0, k = 1");
    }
    return *this;
  }

  template<class Treal>
    Matrix<Treal>& Matrix<Treal>:: 
    operator*=(const Treal alpha) {
    if (!this->is_zero() && alpha != 1) {
      for (int ind = 0; ind < this->nElements(); ind++)
	this->elements[ind] *= alpha;
    }
    return *this;
  }

  template<class Treal> 
    void Matrix<Treal>::
    gemm(const bool tA, const bool tB, const Treal alpha, 
	 const Matrix<Treal>& A, 
	 const Matrix<Treal>& B, 
	 const Treal beta, 
	 Matrix<Treal>& C) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (C.is_empty()) {
      assert(beta == 0);
      if (tA)
	C.resetRows(A.cols);
      else
	C.resetRows(A.rows);
      if (tB)
	C.resetCols(B.rows);
      else
	C.resetCols(B.cols);
    }      
    
    if ((A.is_zero() || B.is_zero() || alpha == 0) && 
	(C.is_zero() || beta == 0))
      C = 0;
    else {
      Treal beta_tmp = beta;
      if (C.is_zero()) {
	C.allocate();
	beta_tmp = 0;
      }
	
      if (!A.is_zero() && !B.is_zero() && alpha != 0) {
	if (!tA && !tB) 
	  if (A.ncols() == B.nrows() && 
	      A.nrows() == C.nrows() &&
	      B.ncols() == C.ncols()) 
	    mat::gemm("N", "N", &A.nrows(), &B.ncols(), &A.ncols(), &alpha,
		      A.elements, &A.nrows(), B.elements, &B.nrows(),
		      &beta_tmp, C.elements, &C.nrows());
	  else 
	    throw Failure("Matrix<Treal>::"
			  "gemm(bool, bool, Treal, Matrix<Treal>, "
			  "Matrix<Treal>, Treal, "
			  "Matrix<Treal>): "
			  "Incorrect matrixdimensions for multiplication");
	else if (tA && !tB)
	  if (A.nrows() == B.nrows() &&
	      A.ncols() == C.nrows() &&
	      B.ncols() == C.ncols()) 
	    mat::gemm("T", "N", &A.ncols(), &B.ncols(), &A.nrows(), &alpha,
		      A.elements, &A.nrows(), B.elements, &B.nrows(),
		      &beta_tmp, C.elements, &C.nrows());
	  else
	    throw Failure("Matrix<Treal>::"
			  "gemm(bool, bool, Treal, Matrix<Treal>, "
			  "Matrix<Treal>, Treal, "
			  "Matrix<Treal>): "
			  "Incorrect matrixdimensions for multiplication");
	else if (!tA && tB)
	  if (A.ncols() == B.ncols() && 
	      A.nrows() == C.nrows() &&
	      B.nrows() == C.ncols()) 
	    mat::gemm("N", "T", &A.nrows(), &B.nrows(), &A.ncols(), &alpha,
		      A.elements, &A.nrows(), B.elements, &B.nrows(),
		      &beta_tmp, C.elements, &C.nrows());
	  else 
	    throw Failure("Matrix<Treal>::"
			  "gemm(bool, bool, Treal, Matrix<Treal>, "
			  "Matrix<Treal>, Treal, "
			  "Matrix<Treal>): "
			  "Incorrect matrixdimensions for multiplication");
	else if (tA && tB)
	  if (A.nrows() == B.ncols() && 
	      A.ncols() == C.nrows() &&
	      B.nrows() == C.ncols()) 
	    mat::gemm("T", "T", &A.ncols(), &B.nrows(), &A.nrows(), &alpha,
		      A.elements, &A.nrows(), B.elements, &B.nrows(),
		      &beta_tmp, C.elements, &C.nrows());
	  else 
	    throw Failure("Matrix<Treal>::"
			  "gemm(bool, bool, Treal, Matrix<Treal>, "
			  "Matrix<Treal>, Treal, "
			  "Matrix<Treal>): "
			  "Incorrect matrixdimensions for multiplication");
	else throw Failure("Matrix<Treal>::"
			   "gemm(bool, bool, Treal, Matrix<Treal>, "
			   "Matrix<Treal>, Treal, "
			   "Matrix<Treal>):Very strange error!!");
      }
      else 
	C *= beta;
    }
  }


  template<class Treal> 
    void Matrix<Treal>::
    symm(const char side, const char uplo, const Treal alpha, 
	 const Matrix<Treal>& A, 
	 const Matrix<Treal>& B, 
	 const Treal beta, 
	 Matrix<Treal>& C) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    assert(A.nrows() == A.ncols());
    //    int dimA = A.nrows();
    if (C.is_empty()) {
      assert(beta == 0);
      if (side =='L') {
	C.resetRows(A.rows);
	C.resetCols(B.cols);
      } 
      else {
	assert(side == 'R');
	C.resetRows(B.rows);
	C.resetCols(A.cols);
      }
    }      
    
    if ((A.is_zero() || B.is_zero() || alpha == 0) && 
	(C.is_zero() || beta == 0))
      C = 0;
    else {
      Treal beta_tmp = beta;
      if (C.is_zero()) {
	C.allocate();
	beta_tmp = 0;
      }
      if (!A.is_zero() && !B.is_zero() && alpha != 0) {
	if (A.nrows() == A.ncols() && C.nrows() == B.nrows() && 
	    C.ncols() == B.ncols() && 
	    ((side == 'L' && A.ncols() == B.nrows()) ||
	     (side == 'R' && B.ncols() == A.nrows())))
	  mat::symm(&side, &uplo, &C.nrows(), &C.ncols(), &alpha, 
		    A.elements, &A.nrows(), B.elements, &B.nrows(),
		    &beta_tmp, C.elements, &C.nrows());
	else
	  throw Failure("Matrix<Treal>::symm: Incorrect matrix "
			"dimensions for symmetric multiply");
      } /* end if (!A.is_zero() && !B.is_zero() && alpha != 0) */
      else 
	C *= beta;
    }
  }

  template<class Treal> 
    void Matrix<Treal>::
    syrk(const char uplo, const bool tA, const Treal alpha, 
	 const Matrix<Treal>& A, 
	 const Treal beta, 
	 Matrix<Treal>& C) {
    assert(!A.is_empty());
    if (C.is_empty()) {
      assert(beta == 0);
      if (tA) {
	C.resetRows(A.cols);
	C.resetCols(A.cols);
      } 
      else {
	C.resetRows(A.rows);
	C.resetCols(A.rows);
      }
    }      
    if (C.nrows() == C.ncols() && 
	((!tA && A.nrows() == C.nrows()) || (tA && A.ncols() == C.nrows())))
      if (alpha != 0 && !A.is_zero()) {
	Treal beta_tmp = beta;
	if (C.is_zero()) {
	  C.allocate();
	  beta_tmp = 0;
	}	
	if (!tA) {
	  mat::syrk(&uplo, "N", &C.nrows(), &A.ncols(), &alpha, 
		    A.elements, &A.nrows(), &beta_tmp, 
		    C.elements, &C.nrows());
	} /* end if (!tA) */
	else
	  mat::syrk(&uplo, "T", &C.nrows(), &A.nrows(), &alpha, 
		    A.elements, &A.nrows(), &beta_tmp, 
		    C.elements, &C.nrows());
      } 
      else
	C *= beta;
    else
      throw Failure("Matrix<Treal>::syrk: Incorrect matrix "
		    "dimensions for symmetric rank-k update");
  }
  
  template<class Treal> 
    void Matrix<Treal>::
    sysq(const char uplo, const Treal alpha, 
	 const Matrix<Treal>& A, 
	 const Treal beta, 
	 Matrix<Treal>& C) {
    assert(!A.is_empty());
    if (C.is_empty()) {
      assert(beta == 0);
      C.resetRows(A.rows);
      C.resetCols(A.cols);
    }
    /* FIXME: slow copy */
    Matrix<Treal> tmpA = A; 
    tmpA.symToNosym();
    Matrix<Treal>::syrk(uplo, false, alpha, tmpA, beta, C);
  }

  /* C = alpha * A * B + beta * C where A and B are symmetric */
  template<class Treal> 
    void Matrix<Treal>::
    ssmm(const Treal alpha, 
	 const Matrix<Treal>& A, 
	 const Matrix<Treal>& B, 
	 const Treal beta, 
	 Matrix<Treal>& C) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (C.is_empty()) {
      assert(beta == 0);
      C.resetRows(A.rows);
      C.resetCols(B.cols);
    }
    /* FIXME: slow copy */
    Matrix<Treal> tmpB = B; 
    tmpB.symToNosym();
    Matrix<Treal>::symm('L', 'U', alpha, A, tmpB, beta, C);
  }
  
  /* C = alpha * A * B + beta * C where A and B are symmetric
   * and only the upper triangle of C is computed.
   */
  template<class Treal> 
    void Matrix<Treal>::
    ssmm_upper_tr_only(const Treal alpha, 
		       const Matrix<Treal>& A, 
		       const Matrix<Treal>& B, 
		       const Treal beta, 
		       Matrix<Treal>& C) {
    /* FIXME: Symmetry is not utilized. */
    ssmm(alpha, A, B, beta, C);
    C.nosymToSym();   
  }


  template<class Treal> 
    void Matrix<Treal>::
    trmm(const char side, const char uplo, const bool tA, 
	 const Treal alpha, 
	 const Matrix<Treal>& A, 
	 Matrix<Treal>& B) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (alpha != 0 && !A.is_zero() && !B.is_zero())
      if (((side == 'R' && B.ncols() == A.nrows()) || 
	   (side == 'L' && A.ncols() == B.nrows())) &&
	  A.nrows() == A.ncols())
	if (!tA)
	  mat::trmm(&side, &uplo, "N", "N", &B.nrows(), &B.ncols(),
		    &alpha, A.elements, &A.nrows(), B.elements, &B.nrows());
	else
	  mat::trmm(&side, &uplo, "T", "N", &B.nrows(), &B.ncols(),
		    &alpha, A.elements, &A.nrows(), B.elements, &B.nrows());
      else
	throw Failure("Matrix<Treal>::trmm: "
		      "Incorrect matrix dimensions for multiplication");
    else
      B = 0;
  }

  template<class Treal>
    Treal Matrix<Treal>::frobSquared() const {
    assert(!this->is_empty());
    if (this->is_zero()) 
      return 0;
    else {
      Treal sum(0);
      for (int i = 0; i < this->nElements(); i++)
	sum += this->elements[i] * this->elements[i];
      return sum;
    }
  }

  template<class Treal>
    Treal Matrix<Treal>::
    syFrobSquared() const {
    assert(!this->is_empty());
    if (this->nrows() != this->ncols()) 
      throw Failure("Matrix<Treal>::syFrobSquared(): "
		    "Matrix must be have equal number of rows "
		    "and cols to be symmetric");
    Treal sum(0);
    if (!this->is_zero()) {
      for (int col = 1; col < this->ncols(); col++)
	for (int row = 0; row < col; row++)
	  sum += 2 * (*this)(row, col) * (*this)(row, col);
      for (int rc = 0; rc < this->ncols(); rc++)
	sum += (*this)(rc, rc) * (*this)(rc, rc);
    }
    return sum;
  }

  template<class Treal>
    Treal Matrix<Treal>::
    frobSquaredDiff
    (const Matrix<Treal>& A,
     const Matrix<Treal>& B) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (A.nrows() != B.nrows() || A.ncols() != B.ncols()) 
      throw Failure("Matrix<Treal>::frobSquaredDiff: "
		    "Incorrect matrix dimensions");
    Treal sum(0);
    if (!A.is_zero() && !B.is_zero()) {
      Treal diff;
      for (int i = 0; i < A.nElements(); i++) {
	diff = A.elements[i] - B.elements[i];
	sum += diff * diff;
      }
    }
    else if (A.is_zero() && !B.is_zero()) 
      sum = B.frobSquared();
    else if (!A.is_zero() && B.is_zero())
      sum = A.frobSquared();
    /* sum is already zero if A.is_zero() && B.is_zero() */
    return sum;
  }

  
  template<class Treal>
    Treal Matrix<Treal>::
    syFrobSquaredDiff
    (const Matrix<Treal>& A,
     const Matrix<Treal>& B) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (A.nrows() != B.nrows() || 
	A.ncols() != B.ncols() || 
	A.nrows() != A.ncols()) 
      throw Failure("Matrix<Treal>::syFrobSquaredDiff: "
		    "Incorrect matrix dimensions");
    Treal sum(0);
    if (!A.is_zero() && !B.is_zero()) {
      Treal diff;
      for (int col = 1; col < A.ncols(); col++)
	for (int row = 0; row < col; row++) {
	  diff = A(row, col) - B(row, col);
	  sum += 2 * diff * diff;
	}
      for (int rc = 0; rc < A.ncols(); rc++) {
	diff = A(rc, rc) - B(rc, rc);
	sum += diff * diff;
      }
    }
    else if (A.is_zero() && !B.is_zero()) 
      sum = B.syFrobSquared();
    else if (!A.is_zero() && B.is_zero())
      sum = A.syFrobSquared();
    /* sum is already zero if A.is_zero() && B.is_zero() */
    return sum;
  }

  template<class Treal> 
    Treal Matrix<Treal>::
    trace() const {
    assert(!this->is_empty());
    if (this->nrows() != this->ncols())
      throw Failure("Matrix<Treal>::trace(): Matrix must be quadratic");  
    if (this->is_zero())
      return 0;
    else {
      Treal sum = 0;
      for (int rc = 0; rc < this->ncols(); rc++)
	sum += (*this)(rc,rc);
      return sum;
    }
  }

  template<class Treal> 
    Treal Matrix<Treal>::
    trace_ab(const Matrix<Treal>& A,
	     const Matrix<Treal>& B) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (A.nrows() != B.ncols() || A.ncols() != B.nrows()) 
      throw Failure("Matrix<Treal>::trace_ab: "
		    "Wrong matrix dimensions for trace(A * B)"); 
    Treal tr = 0;      
    if (!A.is_zero() && !B.is_zero())
      for (int rc = 0; rc < A.nrows(); rc++)
	for (int colA = 0; colA < A.ncols(); colA++)
	  tr += A(rc,colA) * B(colA, rc);
    return tr;
  } 
  
  template<class Treal> 
    Treal  Matrix<Treal>::
    trace_aTb(const Matrix<Treal>& A,
	      const Matrix<Treal>& B) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (A.ncols() != B.ncols() || A.nrows() != B.nrows())
      throw Failure("Matrix<Treal>::trace_aTb: "
		    "Wrong matrix dimensions for trace(A' * B)"); 
    Treal tr = 0;
    if (!A.is_zero() && !B.is_zero())
      for (int rc = 0; rc < A.ncols(); rc++)
	for (int rowB = 0; rowB < B.nrows(); rowB++)
	  tr += A(rowB,rc) * B(rowB, rc);
    return tr;
  }

  template<class Treal> 
    Treal Matrix<Treal>::
    sy_trace_ab(const Matrix<Treal>& A,
		const Matrix<Treal>& B) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (A.nrows() != B.ncols() || A.ncols() != B.nrows() || 
	A.nrows() != A.ncols()) 
      throw Failure("Matrix<Treal>::sy_trace_ab: "
		    "Wrong matrix dimensions for symmetric trace(A * B)");
    if (A.is_zero() || B.is_zero())
      return 0;
    /* Now we know both A and B are nonzero */
    Treal tr = 0;
    /* Diagonal first */
    for (int rc = 0; rc < A.nrows(); rc++)
      tr += A(rc,rc) * B(rc, rc);
    /* Using that trace of transpose is equal to that without transpose: */
    for (int colA = 1; colA < A.ncols(); colA++)
      for (int rowA = 0; rowA < colA; rowA++)
	tr += 2 * A(rowA, colA) * B(rowA, colA);
    return tr;
  }

  template<class Treal>
    void Matrix<Treal>:: 
    add(const Treal alpha,   /* B += alpha * A */
	const Matrix<Treal>& A, 
	Matrix<Treal>& B) {
    assert(!A.is_empty());
    assert(!B.is_empty());
    if (A.nrows() != B.nrows() || A.ncols() != B.ncols())
      throw Failure("Matrix<Treal>::add: "
		    "Wrong matrix dimensions for addition");
    if (!A.is_zero() && alpha != 0) {
      if (B.is_zero())
	B.allocate();
      for (int ind = 0; ind < A.nElements(); ind++)
	B.elements[ind] += alpha * A.elements[ind];
    }
  }

  template<class Treal>
    void Matrix<Treal>:: 
    assign(Treal const  alpha,   /* *this = alpha * A */
	   Matrix<Treal> const & A) {
    assert(!A.is_empty());
    if (this->is_empty()) {
      this->resetRows(A.rows);
      this->resetCols(A.cols);
    }
    Matrix<Treal>::add(alpha, A, *this);
  }
  

  /********** Help functions for thresholding */
  
  template<class Treal>
    void Matrix<Treal>::
    getFrobSqLowestLevel(std::vector<Treal> & frobsq) const {
    if (!this->is_zero()) 
      frobsq.push_back(this->frobSquared());
  }

  template<class Treal>
    void Matrix<Treal>::
    getFrobSqElementLevel(std::vector<Treal> & frobsq) const {
    if (!this->is_zero()) 
      for (int ind = 0; ind < this->nElements(); ind++) 
	if ( this->elements[ind] != 0 ) // Add nonzero elements only
	  frobsq.push_back(this->elements[ind] * this->elements[ind]);    
  }


  template<class Treal>
    void Matrix<Treal>::
    frobThreshLowestLevel
    (Treal const threshold, Matrix<Treal> * ErrorMatrix) {
    if (ErrorMatrix) {
      if ((!this->is_zero() && this->frobSquared() <= threshold) ||
	  (!ErrorMatrix->is_zero() && 
	   ErrorMatrix->frobSquared() > threshold)) 
	Matrix<Treal>::swap(*this,*ErrorMatrix);
    }      
    else if (!this->is_zero() && this->frobSquared() <= threshold) 
      this->clear();
  }

  template<class Treal>
    void Matrix<Treal>::
    frobThreshElementLevel
    (Treal const threshold, Matrix<Treal> * ErrorMatrix) {
    assert(!this->is_empty());
    bool thisMatIsZero = true;
    if (ErrorMatrix) {
      assert(!ErrorMatrix->is_empty());
      bool EMatIsZero = true;
      if (!ErrorMatrix->is_zero() || !this->is_zero()) {
	if (ErrorMatrix->is_zero())
	  ErrorMatrix->allocate();
	if (this->is_zero())
	  this->allocate();
	for (int ind = 0; ind < this->nElements(); ind++) {
	  if ( this->elements[ind] != 0 ) {
	    assert(ErrorMatrix->elements[ind] == 0);
	    // ok, let's check if we want to move the element to the error matrix
	    if ( this->elements[ind] * this->elements[ind] <= threshold ) {
	      // we want to move the element
	      ErrorMatrix->elements[ind] = this->elements[ind];
	      this->elements[ind] = 0;
	      EMatIsZero = false; // at least one element is nonzero
	    }
	    else
	      thisMatIsZero = false; // at least one element is nonzero 
	    continue;
	  }
	  if ( ErrorMatrix->elements[ind] != 0 ) {
	    assert(this->elements[ind] == 0);
	    // ok, let's check if we want to move the element from the error matrix
	    if ( ErrorMatrix->elements[ind] * ErrorMatrix->elements[ind] > threshold ) {
	      // we want to move the element
	      this->elements[ind] = ErrorMatrix->elements[ind];
	      ErrorMatrix->elements[ind] = 0;
	      thisMatIsZero = false; // at least one element is nonzero 
	    }
	    else
	      EMatIsZero = false; // at least one element is nonzero	    
	  }
	}
	if (thisMatIsZero) {
#if 0
	  for (int ind = 0; ind < this->nElements(); ind++) 
	    assert( this->elements[ind] == 0);
#endif
	  this->clear();
	}
	if (EMatIsZero) {
#if 0
	  for (int ind = 0; ind < this->nElements(); ind++) 
	    assert( ErrorMatrix->elements[ind] == 0);
#endif
	  ErrorMatrix->clear();
	}
      }
    }      
    else
      if (!this->is_zero()) {
	for (int ind = 0; ind < this->nElements(); ind++) {
	  if ( this->elements[ind] * this->elements[ind] <= threshold ) 
	    /* FIXME BUG? EMANUEL LOOK AT THIS! */
	    // this->elements[ind] == 0; OLD CODE. Compiler complained that "statement has no effect".
	    this->elements[ind] = 0;  /* New code. Changed from == to =. */
	  else
	    thisMatIsZero = false;
	}
	if (thisMatIsZero)
	  this->clear();
      }    
  }
  

  
  /********** End of help functions for thresholding */

  /* C = beta * C + alpha * A * B where only the upper triangle of C is */
  /* referenced and updated */
  template<class Treal>
    void Matrix<Treal>:: 
    gemm_upper_tr_only(const bool tA, const bool tB, 
		       const Treal alpha, 
		       const Matrix<Treal>& A, 
		       const Matrix<Treal>& B, 
		       const Treal beta, 
		       Matrix<Treal>& C) {
    /* FIXME: Symmetry is not utilized. */
    Matrix<Treal>::gemm(tA, tB, alpha, A, B, beta, C);
    C.nosymToSym();
  }

  /* A = alpha * A * Z or A = alpha * Z * A where A is symmetric, */
  /* Z is upper triangular and */
  /* only the upper triangle of the result is calculated */
  /* side == 'R' for A = alpha * A * Z */
  /* side == 'L' for A = alpha * Z * A */
  template<class Treal>
    void Matrix<Treal>:: 
    sytr_upper_tr_only(char const side, const Treal alpha,
		       Matrix<Treal>& A,
		       const Matrix<Treal>& Z) {
    /* FIXME: Symmetry is not utilized. */
    A.symToNosym(); 
    Matrix<Treal>::trmm(side, 'U', false, alpha, Z, A);
    A.nosymToSym();
  }
  
  /* The result B is assumed to be symmetric, i.e. only upper triangle */
  /* calculated and hence only upper triangle of input B is used. */
  template<class Treal>
    void Matrix<Treal>:: 
    trmm_upper_tr_only(const char side, const char uplo, 
		       const bool tA, const Treal alpha, 
		       const Matrix<Treal>& A, 
		       Matrix<Treal>& B) {
    /* FIXME: Symmetry is not utilized. */
    assert(tA);
    B.symToNosym();
    Matrix<Treal>::trmm(side, uplo, tA, alpha, A, B);
    B.nosymToSym();
  }

  /* A = Z' * A * Z or A = Z * A * Z' */
  /* where A is symmetric and Z is (nonunit) upper triangular */
  /* side == 'R' for A = Z' * A * Z */
  /* side == 'L' for A = Z * A * Z' */
  template<class Treal>
    void Matrix<Treal>:: 
    trsytriplemm(char const side, 
		 const Matrix<Treal>& Z,
		 Matrix<Treal>& A) {
    if (side == 'R') {
      Matrix<Treal>::
	sytr_upper_tr_only('R', 1.0, A, Z);
      Matrix<Treal>::
	trmm_upper_tr_only('L', 'U', true, 1.0, Z, A);
    }
    else {
      assert(side == 'L');
      Matrix<Treal>::
	sytr_upper_tr_only('L', 1.0, A, Z);
      Matrix<Treal>::
	trmm_upper_tr_only('R', 'U', true, 1.0, Z, A);
    }
  }


  template<class Treal>
    Treal Matrix<Treal>::frob_squared_thresh
    (Treal const threshold, Matrix<Treal> * ErrorMatrix) {
    assert(!this->is_empty());
    if (ErrorMatrix && ErrorMatrix->is_empty()) {
      ErrorMatrix->resetRows(this->rows);
      ErrorMatrix->resetCols(this->cols);
    }
    Treal fs = frobSquared();
    if (fs < threshold) {
      if (ErrorMatrix)
	Matrix<Treal>::swap(*this, *ErrorMatrix);
      return fs;
    }
    else
      return 0;
  }


  template<class Treal> 
    void Matrix<Treal>::
    inch(const Matrix<Treal>& A,
	 Matrix<Treal>& Z,
	 const Treal threshold, const side looking,
	 const inchversion version) {
    assert(!A.is_empty());
    if (A.nrows() != A.ncols()) 
      throw Failure("Matrix<Treal>::inch: Matrix must be quadratic!");
    if (A.is_zero())
      throw Failure("Matrix<Treal>::inch: Matrix is zero! "
		    "Not possible to compute inverse cholesky.");
    Z = A;
    int info;
    potrf("U", &A.nrows(), Z.elements, &A.nrows(), &info);
    if (info > 0)
      throw Failure("Matrix<Treal>::inch: potrf  returned info > 0. The matrix is not positive definite.");
    if (info < 0)
      throw Failure("Matrix<Treal>::inch: potrf  returned info < 0");
      
    trtri("U", "N", &A.nrows(), Z.elements, &A.nrows(), &info);
    if (info > 0)
      throw Failure("Matrix<Treal>::inch: trtri  returned info > 0. The matrix is not positive definite.");
    if (info < 0)
      throw Failure("Matrix<Treal>::inch: trtri  returned info < 0");
    /* Fill lower triangle with zeroes just to be safe */
    trifulltofull(Z.elements, A.nrows()); 
  }

  template<class Treal>
    void Matrix<Treal>:: 
    gershgorin(Treal& lmin, Treal& lmax) const {
    assert(!this->is_empty());
    if (this->nScalarsRows() == this->nScalarsCols()) {
      int size = this->nScalarsCols();
      Treal* colsums = new Treal[size];
      Treal* diag    = new Treal[size];
      for (int ind = 0; ind < size; ind++)
	colsums[ind] = 0;
      this->add_abs_col_sums(colsums);
      this->get_diagonal(diag);
      Treal tmp1 = colsums[0] - template_blas_fabs(diag[0]);
      Treal tmp2;
      lmin = diag[0] - tmp1;
      lmax = diag[0] + tmp1;
      for (int col = 1; col < size; col++) {
	tmp1 = colsums[col] - template_blas_fabs(diag[col]);
	tmp2 = diag[col] - tmp1;
	lmin = (tmp2 < lmin ? tmp2 : lmin);
	tmp2 = diag[col] + tmp1;
	lmax = (tmp2 > lmax ? tmp2 : lmax);
      }
      delete[] diag;
      delete[] colsums;
    }
    else
      throw Failure("Matrix<Treal>::gershgorin(Treal, Treal): Matrix must be quadratic");
  }


  template<class Treal> 
    void Matrix<Treal>::
    add_abs_col_sums(Treal* sums) const {
    assert(sums);
    if (!this->is_zero()) 
      for (int col = 0; col < this->ncols(); col++) 
	for (int row = 0; row < this->nrows(); row++) 
	  sums[col] += template_blas_fabs((*this)(row,col));
  }
 
  template<class Treal> 
    void Matrix<Treal>::
    get_diagonal(Treal* diag) const {
    assert(diag);
    assert(this->nrows() == this->ncols());
    if (this->is_zero()) 
      for (int rc = 0; rc < this->nScalarsCols(); rc++)
	diag[rc] = 0;
    else 
      for (int rc = 0; rc < this->ncols(); rc++) 
	diag[rc] = (*this)(rc,rc);
  }


    /* New routines above */

#if 0 /* OLD ROUTINES */   


  



 
   

  


  





  

  template<class Treal>
    void Matrix<Treal>::trim_memory_usage() {
    if (this->is_zero() && this->cap > 0) {
      freeElements(this->elements);
      this->elements = NULL;
      this->cap = 0;
      this->nel = 0;
    }
    else if (this->cap > this->nel) {
      Treal* tmp = new Treal[this->nel];
      for (int i = 0; i < this->nel; i++) 
	tmp[i] = this->elements[i];
      freeElements(this->elements);
      this->cap = this->nel;
      this->elements = tmp;
    }
    assert(this->cap == this->nel);
  }



#if 1

  template<class Treal>
    void Matrix<Treal>::
    write_to_buffer_count(int& zb_length, int& vb_length) const {
    if (this->is_zero()) {
      ++zb_length;	
    }
    else {
      ++zb_length;	
      vb_length += this->nel;
    }
  }

  template<class Treal>
    void Matrix<Treal>::
    write_to_buffer(int* zerobuf, const int zb_length, 
		    Treal* valuebuf, const int vb_length,
		    int& zb_index, int& vb_index) const {
    if (zb_index < zb_length) {
      if (this->is_zero()) {
	zerobuf[zb_index] = 0;
	++zb_index;	
      }
      else {
	if (vb_index + this->nel < vb_length + 1) {
	  zerobuf[zb_index] = 1;
	  ++zb_index;	
	  for (int i = 0; i < this->nel; i++)
	    valuebuf[vb_index + i] = this->elements[i];
	  vb_index += this->nel;
	}
	else
	  throw Failure("Matrix<Treal, Telement>::write_to_buffer: "
			"Insufficient space in buffers");   
      }
    }
    else
      throw Failure("Matrix<Treal, Telement>::write_to_buffer: "
		    "Insufficient space in buffers");   
  }

  template<class Treal>
    void Matrix<Treal>::
    read_from_buffer(int* zerobuf, const int zb_length, 
		     Treal* valuebuf, const int vb_length,
		     int& zb_index, int& vb_index) {
    if (zb_index < zb_length) {
      if (zerobuf[zb_index] == 0) {
	(*this) = 0;
	++zb_index;	
      }
      else {
	this->content = ful;
	this->nel = this->nrows() * this->ncols();
	this->assert_alloc();	
	if (vb_index + this->nel < vb_length + 1) {
	  assert(zerobuf[zb_index] == 1);
	  ++zb_index;	
	  for (int i = 0; i < this->nel; i++)
	    this->elements[i] = valuebuf[vb_index + i];
	  vb_index += this->nel;
	}
	else
	  throw Failure("Matrix<Treal, Telement>::read_from_buffer: "
			"Mismatch, buffers does not match matrix");   
      }
    }
    else
      throw Failure("Matrix<Treal, Telement>::read_from_buffer: "
		    "Mismatch, buffers does not match matrix");   
  }
  
#endif





  


  /* continue here */


    

    









#if 1

  
  
#endif

#endif     /* END OF OLD ROUTINES */   


} /* end namespace mat */

#endif
