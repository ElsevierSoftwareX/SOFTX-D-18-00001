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

/** @file MatrixSymmetric.h Symmetric matrix class
 *
 * Copyright(c) Emanuel Rubensson 2006
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date January 2006
 *
 */
#ifndef MAT_MatrixSymmetric
#define MAT_MatrixSymmetric

#include <algorithm>

#include "MatrixBase.h"
#include "Interval.h"
#include "LanczosLargestMagnitudeEig.h"
#include "mat_utils.h"
#include "truncation.h"

namespace mat {
  /** Symmetric matrix 
   *
   *
   * This class belongs to the matrix API
   *
   * The matrix is stored in the upper triangle.  
   *
   * Treal: Type for real numbers
   *
   * Tmatrix: The matrix class
   *
   * @see MatrixBase
   * @see MatrixGeneral
   * @see MatrixTriangular
   * 
   * 
   */  
  template<typename Treal, typename Tmatrix>
    class MatrixSymmetric : public MatrixBase<Treal, Tmatrix> {
  public:
    typedef VectorGeneral<Treal, typename Tmatrix::VectorType> VectorType;
    typedef Treal real;

    MatrixSymmetric()
      :MatrixBase<Treal, Tmatrix>() {} /**< Default constructor  */
    explicit MatrixSymmetric(const MatrixSymmetric<Treal, Tmatrix>& symm)
      :MatrixBase<Treal, Tmatrix>(symm) {} /**< Copy constructor  */
    explicit MatrixSymmetric(const XY<Treal, MatrixSymmetric<Treal, Tmatrix> >& sm)
      :MatrixBase<Treal, Tmatrix>() { *this = sm.A * sm.B; }
    explicit MatrixSymmetric(const MatrixGeneral<Treal, Tmatrix>& matr)
      :MatrixBase<Treal, Tmatrix>(matr) {
      this->matrixPtr->nosymToSym();
    } /**< 'Copy from normal matrix' - constructor  */
    

#if 0
    template<typename Tfull>
      inline void assign_from_full
      (Tfull const* const fullmatrix, 
       int const totnrows, int const totncols) {
       assert(totnrows == totncols);
      this->matrixPtr->assign_from_full(fullmatrix, totnrows, totncols);      
      this->matrixPtr->nosym_to_sym();
    }    
    inline void assign_from_full
      (Treal const* const fullmatrix, 
       int const totnrows, int const totncols) {
      assert(totnrows == totncols);
      this->matrixPtr->assign_from_full(fullmatrix, totnrows, totncols);      
      this->matrixPtr->nosym_to_sym();
    }
#endif

    inline void assignFromFull
      (std::vector<Treal> const & fullMat) {
      assert((int)fullMat.size() == this->get_nrows() * this->get_ncols());
      this->matrixPtr->assignFromFull(fullMat);      
      this->matrixPtr->nosymToSym();
    }
    
    inline void assignFromFull
      (std::vector<Treal> const & fullMat,
       std::vector<int> const & rowPermutation, 
       std::vector<int> const & colPermutation) {
      assert((int)fullMat.size() == this->get_nrows() * this->get_ncols());
      std::vector<int> rowind(fullMat.size());
      std::vector<int> colind(fullMat.size());
      int ind = 0;
      for (int col = 0; col < this->get_ncols(); ++col)
	for (int row = 0; row < this->get_nrows(); ++row) {
	  rowind[ind] = row;
	  colind[ind] = col;
	  ++ind;
	}
      this->assign_from_sparse(rowind, 
			       colind, 
			       fullMat, 
			       rowPermutation, 
			       colPermutation);
      this->matrixPtr->nosymToSym();
    }

    inline void fullMatrix(std::vector<Treal> & fullMat) const {
      this->matrixPtr->syFullMatrix(fullMat);
    }
    /**< Save matrix as full matrix.
     * Whole matrix is written in columnwise order.
     * Both lower and upper triangle.
     * NOTE that no permutation is used in this operation.
     */
    
    inline void fullMatrix
      (std::vector<Treal> & fullMat,
       std::vector<int> const & rowInversePermutation, 
       std::vector<int> const & colInversePermutation) const {
      std::vector<int> rowind;
      std::vector<int> colind; 
      std::vector<Treal> values;
      get_all_values(rowind, colind, values, 
		     rowInversePermutation, 
		     colInversePermutation);
      fullMat.assign(this->get_nrows()*this->get_ncols(),0);
      assert(rowind.size() == values.size());
      assert(rowind.size() == colind.size());
      for (unsigned int ind = 0; ind < values.size(); ++ind) {
	assert(rowind[ind] + this->get_nrows() * colind[ind] < 
	       this->get_nrows()*this->get_ncols());
	fullMat[rowind[ind] + this->get_nrows() * colind[ind]] =
	  values[ind];
	fullMat[colind[ind] + this->get_nrows() * rowind[ind]] =
	  values[ind];
      }
    }
    /**< Save matrix as full matrix.
     * Whole matrix is written in columnwise order.
     * Both lower and upper triangle.
     * Permutation is used.
     */

    inline void assign_from_sparse
      (std::vector<int> const & rowind, 
       std::vector<int> const & colind, 
       std::vector<Treal> const & values) {
      this->matrixPtr->syAssignFromSparse(rowind, colind, values);
    }
    /**< Assign from sparse matrix given by three vectors. 
     * The vectors contain row indices, column indices and values.
     * The indices start at zero.
     * The elements to be added must be given in upper triangluar storage.
     * Information about sizes and blocks for rows as well as columns 
     * must also be given.
     * Assumes that sizes and blocks are already set.
     * @warning All indexing start at zero.
     */


    inline void assign_from_sparse
      (std::vector<int> const & rowind, 
       std::vector<int> const & colind, 
       std::vector<Treal> const & values, 
       SizesAndBlocks const & newRows,
       SizesAndBlocks const & newCols) {
      this->resetSizesAndBlocks(newRows, newCols);
      this->matrixPtr->syAssignFromSparse(rowind, colind, values);
    }
    /**< Assign from sparse matrix given by three vectors. 
     * The vectors contain row indices, column indices and values.
     * The indices start at zero.
     * The elements to be added must be given in upper triangluar storage.
     * Information about sizes and blocks for rows as well as columns 
     * must also be given.
     * @warning All indexing start at zero.
     */
 
    inline void assign_from_sparse
      (std::vector<int> const & rowind, 
       std::vector<int> const & colind, 
       std::vector<Treal> const & values, 
       std::vector<int> const & rowPermutation, 
       std::vector<int> const & colPermutation) {
      std::vector<int> newRowind;
      std::vector<int> newColind; 
      this->getPermutedAndSymmetrized(rowind, rowPermutation, newRowind,
				      colind, colPermutation, newColind);
      
      this->matrixPtr->syAssignFromSparse(newRowind, newColind, values);
    }
    /**< Same as above, except taking two additional arguments 
     *   specifying the permutation of rows and columns.
     *   Also, assuming that sizes and blocks are already set.
     */

    inline void assign_from_sparse
      (std::vector<int> const & rowind, 
       std::vector<int> const & colind, 
       std::vector<Treal> const & values, 
       SizesAndBlocks const & newRows,
       SizesAndBlocks const & newCols,
       std::vector<int> const & rowPermutation, 
       std::vector<int> const & colPermutation) {
      this->resetSizesAndBlocks(newRows, newCols);
      assign_from_sparse(rowind, colind, values, 
			 rowPermutation, colPermutation);
    }
    /**< Same as above, except taking sizes and blocks arguments.
     */

    
    /** Add given set of values to the matrix. 
     *  The values should be given in upper triangular storage.
     */
    inline void add_values
      (std::vector<int> const & rowind, 
       std::vector<int> const & colind, 
       std::vector<Treal> const & values) {
      this->matrixPtr->syAddValues(rowind, colind, values);
    }

    /** Same as above, except taking two additional arguments 
     *  specifying the permutation of rows and columns.
     */
    inline void add_values
      (std::vector<int> const & rowind, 
       std::vector<int> const & colind, 
       std::vector<Treal> const & values,
       std::vector<int> const & rowPermutation, 
       std::vector<int> const & colPermutation) {
      std::vector<int> newRowind;
      std::vector<int> newColind; 
      this->getPermutedAndSymmetrized(rowind, rowPermutation, newRowind,
				      colind, colPermutation, newColind);
      this->matrixPtr->syAddValues(newRowind, newColind, values);
    }



    inline void get_values
      (std::vector<int> const & rowind, 
       std::vector<int> const & colind, 
       std::vector<Treal> & values) const {
      this->matrixPtr->syGetValues(rowind, colind, values);
    }
    /**< Get values given by row and column index lists.
     * Input arrays contain row and column indices.
     * The wanted elements must be given in upper triangluar storage.
     * The output array contains values for the given indices.
     * @warning All indexing start at zero.
     */    

    inline void get_values
      (std::vector<int> const & rowind, 
       std::vector<int> const & colind, 
       std::vector<Treal> & values,
       std::vector<int> const & rowPermutation, 
       std::vector<int> const & colPermutation) const {
      std::vector<int> newRowind;
      std::vector<int> newColind; 
      this->getPermutedAndSymmetrized(rowind, rowPermutation, newRowind,
				      colind, colPermutation, newColind);
      this->matrixPtr->syGetValues(newRowind, newColind, values);
    }
    /**< Same as above, except taking two additional arguments 
     *   specifying the permutation of rows and columns.
     */    
    

    inline void get_all_values
      (std::vector<int> & rowind, 
       std::vector<int> & colind, 
       std::vector<Treal> & values) const {
      rowind.resize(0);
      colind.resize(0);
      values.resize(0);
      rowind.reserve(nnz());
      colind.reserve(nnz());
      values.reserve(nnz());
      this->matrixPtr->syGetAllValues(rowind, colind, values);
    }
    /**< Get all values and corresponding row and column index lists,
     * in matrix. Only upper triangle values are returned. 
     * @warning All indexing start at zero.
     */    

    inline void get_all_values
      (std::vector<int> & rowind, 
       std::vector<int> & colind, 
       std::vector<Treal> & values,
       std::vector<int> const & rowInversePermutation, 
       std::vector<int> const & colInversePermutation) const {
      std::vector<int> tmpRowind;
      std::vector<int> tmpColind; 
      tmpRowind.reserve(rowind.capacity());
      tmpColind.reserve(colind.capacity());
      values.resize(0);
      this->matrixPtr->syGetAllValues(tmpRowind, tmpColind, values);
      this->getPermutedAndSymmetrized(tmpRowind, rowInversePermutation, rowind,
				      tmpColind, colInversePermutation, colind);
    }
    /**< Same as above, except taking two additional arguments 
     *   specifying the permutation of rows and columns.
     *   Note, however, that this permutation is the inverse 
     *   permutation compared to the permutations provided in the
     *   functions "assign_from_sparse", "add_values", and "get_values"
     *   @warning permutation is inverse compared to other functions
     */    



    MatrixSymmetric<Treal, Tmatrix>& 
      operator=(const MatrixSymmetric<Treal, Tmatrix>& symm) {
      MatrixBase<Treal, Tmatrix>::operator=(symm);
      return *this;
    } 
    MatrixSymmetric<Treal, Tmatrix>& 
      operator=(const MatrixGeneral<Treal, Tmatrix>& matr) {
      MatrixBase<Treal, Tmatrix>::operator=(matr);
      this->matrixPtr->nosymToSym();
      return *this;
    } 
    inline MatrixSymmetric<Treal, Tmatrix>& operator=(int const k) {
      *this->matrixPtr = k;
      return *this;
    }

    inline Treal frob() const {
      return this->matrixPtr->syFrob();
    }
    Treal mixed_norm(Treal const requestedAccuracy,
		     int maxIter = -1) const;
    Treal eucl(Treal const requestedAccuracy,
	       int maxIter = -1) const;

    void quickEuclBounds(Treal & euclLowerBound, 
			 Treal & euclUpperBound) const {
      Treal frobTmp = frob();
      euclLowerBound = frobTmp  / template_blas_sqrt( (Treal)this->get_nrows() );
      euclUpperBound = frobTmp;
    }


    /** Returns interval containing the Euclidean norm of A - B 
     * ( || A - B ||_2 )
     * @see eucl_diff
     * @see frob_diff 
     */
      static Interval<Treal> 
	diff(const MatrixSymmetric<Treal, Tmatrix>& A,
	     const MatrixSymmetric<Treal, Tmatrix>& B,
	     normType const norm, 
	     Treal const requestedAccuracy);
    /** Returns interval containing the Euclidean norm of A - B 
     *  ( || A - B ||_2 ) based on the chosen norm.
     *  BUT, in the case of Euclidean norm, the norm is only computed with
     *  the requested accuracy if it is smaller than 'maxAbsVal'.
     *  @see euclDiffIfSmall
     *  @see frob_diff
     */
      static Interval<Treal> 
      diffIfSmall(const MatrixSymmetric<Treal, Tmatrix>& A,
		  const MatrixSymmetric<Treal, Tmatrix>& B,
		  normType const norm, 
		  Treal const requestedAccuracy,
		  Treal const maxAbsVal);
    /** Returns the Frobenius norm of A - B 
     *  ( || A - B ||_F )
     */
    static inline Treal frob_diff
      (const MatrixSymmetric<Treal, Tmatrix>& A,
       const MatrixSymmetric<Treal, Tmatrix>& B) {
      return Tmatrix::syFrobDiff(*A.matrixPtr, *B.matrixPtr);
    }

    /** Returns the Euclidean norm of A - B 
     *  ( || A - B ||_2 )
     */
    static Treal eucl_diff
      (const MatrixSymmetric<Treal, Tmatrix>& A,
       const MatrixSymmetric<Treal, Tmatrix>& B,
       Treal const requestedAccuracy,
       int maxIter = -1);
    
    /** Returns the 'mixed' norm of A - B 
     *  ( || A - B ||_mixed )
     */
    static Treal mixed_diff
      (const MatrixSymmetric<Treal, Tmatrix>& A,
       const MatrixSymmetric<Treal, Tmatrix>& B,
       Treal const requestedAccuracy);
    
    /** Returns interval containing the Euclidean norm of A - B 
     *  ( || A - B ||_2 ).
     *  BUT, the norm is only computed with
     *  the requested accuracy if it is smaller than 'maxAbsVal'.
     *  Otherwise, the Frobenius norm is used to get the bounds.
     */
    static Interval<Treal> euclDiffIfSmall
      (const MatrixSymmetric<Treal, Tmatrix>& A,
       const MatrixSymmetric<Treal, Tmatrix>& B,
       Treal const requestedAccuracy,
       Treal const maxAbsVal,
       VectorType * const eVecPtr = 0);
      

    /** Does thresholding so that the error in the chosen norm is below
     *  the given threshold. Returns the actual introduced error.
     *  In case of the Frobenius norm the return value may be an upper bound.
     *  In case of the Euclidean norm the return value is sometimes an 
     *  upper bound as well but it can only happen if the whole matrix 
     *  is removed. 
     *     
     *  @see frob_thresh(Treal)
     *  @see eucl_thresh(Treal const)
     */
    Treal thresh(Treal const threshold,
		 normType const norm);
    
    /** Does thresholding so that the error in the Frobenius norm
     *  is below the given threshold.
     *  Returns an upper bound of the introduced error.
     *  If no elements on the block diagonal are removed the return value
     *  is equal to the introduced error.
     */
    inline Treal frob_thresh(Treal const threshold) {
      return 2.0 * this->matrixPtr->frob_thresh(threshold / 2);
    } 
    
    Treal eucl_thresh(Treal const threshold,
		      MatrixTriangular<Treal, Tmatrix> const * const Zptr = NULL); 
    
    Treal eucl_element_level_thresh(Treal const threshold); 
    
    void getSizesAndBlocksForFrobNormMat
      ( SizesAndBlocks & rows_new, SizesAndBlocks & cols_new ) const;

    Treal mixed_norm_thresh(Treal const threshold);
    
    void simple_blockwise_frob_thresh(Treal const threshold) {
      this->matrixPtr->frobThreshLowestLevel(threshold*threshold, 0);
    }
    
    inline void gershgorin(Treal& lmin, Treal& lmax) const {
      this->matrixPtr->sy_gershgorin(lmin, lmax);
    }
    static inline Treal trace_ab
      (const MatrixSymmetric<Treal, Tmatrix>& A,
       const MatrixSymmetric<Treal, Tmatrix>& B) {
      return Tmatrix::sy_trace_ab(*A.matrixPtr, *B.matrixPtr);
    }
    inline size_t nnz() const { /* Note: size_t instead of int here to avoid integer overflow. */
      return this->matrixPtr->sy_nnz();
    }
    inline size_t nvalues() const { /* Note: size_t instead of int here to avoid integer overflow. */
      return this->matrixPtr->sy_nvalues();
    }
    inline void write_to_buffer(void* buffer, const int n_bytes) const {
      this->write_to_buffer_base(buffer, n_bytes, matrix_symm);
    }
    inline void read_from_buffer(void* buffer, const int n_bytes) {
      this->read_from_buffer_base(buffer, n_bytes, matrix_symm);
    }


    /** B = alpha * A   : A and B are symmetric*/
    MatrixSymmetric<Treal, Tmatrix>& operator= 
      (XY<Treal, MatrixSymmetric<Treal, Tmatrix> > const & sm);
    /** C = alpha * A * A + beta * C          : A and C are symmetric    */
    MatrixSymmetric<Treal, Tmatrix>& operator=
      (const XYZpUV<Treal, 
       MatrixSymmetric<Treal, Tmatrix>, 
       MatrixSymmetric<Treal, Tmatrix>,
       Treal, 
       MatrixSymmetric<Treal, Tmatrix> >& sm2psm);
    /** C = alpha * A * A                     : A and C are symmetric    */
    MatrixSymmetric<Treal, Tmatrix>& operator=
      (const XYZ<Treal, 
       MatrixSymmetric<Treal, Tmatrix>, 
       MatrixSymmetric<Treal, Tmatrix> >& sm2);
    /** C += alpha * A * A                    : A and C are symmetric    */
    MatrixSymmetric<Treal, Tmatrix>& operator+=
      (const XYZ<Treal, 
       MatrixSymmetric<Treal, Tmatrix>, 
       MatrixSymmetric<Treal, Tmatrix> >& sm2);
    /** C = alpha * A * transpose(A) + beta * C        : C is symmetric  */
    MatrixSymmetric<Treal, Tmatrix>& operator= 
      (const XYZpUV<Treal, 
       MatrixGeneral<Treal, Tmatrix>,
       MatrixGeneral<Treal, Tmatrix>,
       Treal,
       MatrixSymmetric<Treal, Tmatrix> >& smmpsm);
    /** C = alpha * A * transpose(A)                   : C is symmetric  */
    MatrixSymmetric<Treal, Tmatrix>& operator= 
      (const XYZ<Treal, 
       MatrixGeneral<Treal, Tmatrix>, 
       MatrixGeneral<Treal, Tmatrix> >& smm);
    /** C += alpha * A * transpose(A)                   : C is symmetric  */
    MatrixSymmetric<Treal, Tmatrix>& operator+= 
      (const XYZ<Treal, 
       MatrixGeneral<Treal, Tmatrix>, 
       MatrixGeneral<Treal, Tmatrix> >& smm);
    
    /** A = Z * A * transpose(Z) : Z is upper triangular and A is symmetric;
     *  A = transpose(Z) * A * Z : Z is upper triangular and A is symmetric 
     */
    MatrixSymmetric<Treal, Tmatrix>& operator= 
      (const XYZ<MatrixTriangular<Treal, Tmatrix>, 
       MatrixSymmetric<Treal, Tmatrix>, 
       MatrixTriangular<Treal, Tmatrix> >& zaz);
    
    /** C = alpha * A * B + beta * C where A and B are symmetric
     *  and only the upper triangle of C is computed, 
     *  C is enforced to be symmetric!
     */
    static void ssmmUpperTriangleOnly(const Treal alpha, 
				      const MatrixSymmetric<Treal, Tmatrix>& A, 
				      const MatrixSymmetric<Treal, Tmatrix>& B, 
				      const Treal beta, 
				      MatrixSymmetric<Treal, Tmatrix>& C);
    
    
    /* Addition */    
    /** C =  A + B       : A, B, and C are symmetric  */
    MatrixSymmetric<Treal, Tmatrix>& operator= 
      (XpY<MatrixSymmetric<Treal, Tmatrix>,
       MatrixSymmetric<Treal, Tmatrix> > const & mpm);
    /** C =  A - B       : A, B, and C are symmetric  */
    MatrixSymmetric<Treal, Tmatrix>& operator= 
      (XmY<MatrixSymmetric<Treal, Tmatrix>,
       MatrixSymmetric<Treal, Tmatrix> > const & mm);
    /** B += A           : A and B are symmetric */
    MatrixSymmetric<Treal, Tmatrix>& operator+= 
      (MatrixSymmetric<Treal, Tmatrix> const & A);
    /** B -= A           : A and B are symmetric */
    MatrixSymmetric<Treal, Tmatrix>& operator-=
      (MatrixSymmetric<Treal, Tmatrix> const & A);

    /** B += alpha * A   : A and B are symmetric*/
    MatrixSymmetric<Treal, Tmatrix>& operator+= 
      (XY<Treal, MatrixSymmetric<Treal, Tmatrix> > const & sm);

    /** B -= alpha * A   : A and B are symmetric*/
    MatrixSymmetric<Treal, Tmatrix>& operator-= 
      (XY<Treal, MatrixSymmetric<Treal, Tmatrix> > const & sm);
    
    template<typename Top>
      Treal accumulateWith(Top & op) {
      return this->matrixPtr->syAccumulateWith(op);
    }

    void random() {
      this->matrixPtr->syRandom();
    }
    
    void randomZeroStructure(Treal probabilityBeingZero) {
      this->matrixPtr->syRandomZeroStructure(probabilityBeingZero);
    }
    
    /** Uses rule depending on the row and column indexes to set matrix elements
     * The Trule class should have the function "Treal = set(int row,int col)"
     * which is used to set the elements.
     */
    template<typename TRule>
      void setElementsByRule(TRule & rule) {
      this->matrixPtr->sySetElementsByRule(rule);
      return;
    }

    /** Transfer this matrix to dest, clearing previous content of
	dest if any. */ 
    void transfer( MatrixSymmetric<Treal, Tmatrix> & dest ) {
      ValidPtr<Tmatrix>::swap( this->matrixPtr, dest.matrixPtr );
      // *this now contains previous content of dest
      this->clear(); 
    }

    template<typename Tvector>
      void matVecProd(Tvector & y, Tvector const & x) const {
      Treal const ONE = 1.0;
      y = (ONE * (*this) * x);
    }
    
    
    std::string obj_type_id() const {return "MatrixSymmetric";}
  protected:
    inline void writeToFileProt(std::ofstream & file) const {
      this->writeToFileBase(file, matrix_symm);
    }
    inline void readFromFileProt(std::ifstream & file) {
      this->readFromFileBase(file, matrix_symm);
    }
    
    /** This function permutes row and column indices according to the 
     *  specified permutation and gives the indices as upper triangle 
     *  in the new permutation.
     *  @warning Duplicate indices are kept.
     */
    static void getPermutedAndSymmetrized
      (std::vector<int> const & rowind, 
       std::vector<int> const & rowPermutation,
       std::vector<int> & newRowind,
       std::vector<int> const & colind, 
       std::vector<int> const & colPermutation,
       std::vector<int> & newColind) {
      MatrixBase<Treal, Tmatrix>::
	getPermutedIndexes(rowind, rowPermutation, newRowind);
      MatrixBase<Treal, Tmatrix>::
	getPermutedIndexes(colind, colPermutation, newColind);
      int tmp;
      for (unsigned int i = 0; i < newRowind.size(); ++i) {
	if (newRowind[i] > newColind[i]) {
	  tmp = newRowind[i];
	  newRowind[i] = newColind[i];
	  newColind[i] = tmp;
	}
      }
    }
  private:
  };
  
  template<typename Treal, typename Tmatrix>
    Treal MatrixSymmetric<Treal, Tmatrix>::
    mixed_norm(Treal const requestedAccuracy,
	       int maxIter) const {
    // Construct SizesAndBlocks for frobNormMat
    SizesAndBlocks rows_new;
    SizesAndBlocks cols_new;
    this->getSizesAndBlocksForFrobNormMat( rows_new, cols_new );
    // Now we can construct an empty matrix where the Frobenius norms
    // of lowest level nonzero submatrices will be stored
    MatrixSymmetric<Treal, typename Tmatrix::ElementType> frobNormMat;
    frobNormMat.resetSizesAndBlocks(rows_new, cols_new);
    frobNormMat.getMatrix().syAssignFrobNormsLowestLevel(this->getMatrix());
    return frobNormMat.eucl(requestedAccuracy, maxIter);
  }


  template<typename Treal, typename Tmatrix>
    Treal MatrixSymmetric<Treal, Tmatrix>::
    eucl(Treal const requestedAccuracy,
	 int maxIter) const {
    assert(requestedAccuracy >= 0);
    /* Check if norm is really small, in that case quick return */
    Treal frobTmp = this->frob(); 
    if (frobTmp < requestedAccuracy)
      return (Treal)0.0;
    if (maxIter < 0)
      maxIter = this->get_nrows() * 100;
    VectorType guess;
    SizesAndBlocks cols;
    this->getCols(cols);
    guess.resetSizesAndBlocks(cols);
    guess.rand();
    // Elias note 2010-03-26: changed this back from "new code" to "old code" to reduce memory usage. 
#if 0 // "new code"
    MatrixSymmetric<Treal, Tmatrix> Copy(*this);
    Copy.frob_thresh(requestedAccuracy / 2.0);
    arn::LanczosLargestMagnitudeEig
      <Treal, MatrixSymmetric<Treal, Tmatrix>, VectorType>
      lan(Copy, guess, maxIter);
    lan.setAbsTol( requestedAccuracy / 2.0 );    
#else // "old code"
    arn::LanczosLargestMagnitudeEig
      <Treal, MatrixSymmetric<Treal, Tmatrix>, VectorType>
      lan(*this, guess, maxIter);
    lan.setAbsTol( requestedAccuracy );    
#endif
    lan.run();
    Treal eVal = 0;
    Treal acc = 0;
    lan.getLargestMagnitudeEig(eVal, acc);
    return template_blas_fabs(eVal);
  }

  template<typename Treal, typename Tmatrix>
    Interval<Treal> MatrixSymmetric<Treal, Tmatrix>::
    diff(const MatrixSymmetric<Treal, Tmatrix>& A,
	 const MatrixSymmetric<Treal, Tmatrix>& B,
	 normType const norm, Treal const requestedAccuracy) {
    Treal diff;
    Treal eNMin;
    switch (norm) {
    case frobNorm:
      diff = frob_diff(A, B);
      return Interval<Treal>(diff / template_blas_sqrt((Treal)A.get_nrows()), diff);
      break;
    case euclNorm:
      diff = eucl_diff(A, B, requestedAccuracy);
      eNMin = diff - requestedAccuracy;
      eNMin = eNMin >= 0 ? eNMin : 0; 
      return Interval<Treal>(eNMin, diff + requestedAccuracy);
      break;
    default:
      throw Failure("MatrixSymmetric<Treal, Tmatrix>::"
		    "diff(const MatrixSymmetric<Treal, Tmatrix>&, "
		    "const MatrixSymmetric<Treal, Tmatrix>&, "
		    "normType const, Treal): "
		    "Diff not implemented for selected norm");
    }    
  }
  
#if 1
  template<typename Treal, typename Tmatrix>
    Interval<Treal> MatrixSymmetric<Treal, Tmatrix>::
    diffIfSmall(const MatrixSymmetric<Treal, Tmatrix>& A,
		const MatrixSymmetric<Treal, Tmatrix>& B,
		normType const norm, 
		Treal const requestedAccuracy,
		Treal const maxAbsVal) {
    Treal diff;
    switch (norm) {
    case frobNorm:
      {
	diff = frob_diff(A, B);
	return Interval<Treal>(diff / template_blas_sqrt((Treal)A.get_nrows()), diff);
      }
      break;
    case euclNorm:
      return euclDiffIfSmall(A, B, requestedAccuracy, maxAbsVal);
      break;
    default:
      throw Failure("MatrixSymmetric<Treal, Tmatrix>::"
		    "diffIfSmall"
		    "(const MatrixSymmetric<Treal, Tmatrix>&, "
		    "const MatrixSymmetric<Treal, Tmatrix>&, "
		    "normType const, Treal const, Treal const): "
		    "Diff not implemented for selected norm");
    }    
  }

#endif


  template<typename Treal, typename Tmatrix>
    Treal MatrixSymmetric<Treal, Tmatrix>::eucl_diff
    (const MatrixSymmetric<Treal, Tmatrix>& A,
     const MatrixSymmetric<Treal, Tmatrix>& B,
     Treal const requestedAccuracy,
     int maxIter) {
    // DiffMatrix is a lightweight proxy object:
    mat::DiffMatrix< MatrixSymmetric<Treal, Tmatrix>, Treal> Diff(A, B);
    Treal maxAbsVal = 2 * frob_diff(A,B); 
    // Now, maxAbsVal should be larger than the Eucl norm
    // Note that mat::euclIfSmall lies outside this class 
    Treal relTol = template_blas_sqrt(template_blas_sqrt(mat::getMachineEpsilon<Treal>()));
    VectorType * const eVecPtrNotUsed = 0;
    Interval<Treal> euclInt = 
      mat::euclIfSmall(Diff, requestedAccuracy, relTol, maxAbsVal, eVecPtrNotUsed, maxIter);
    return euclInt.midPoint();
  }

  template<typename Treal, typename Tmatrix>
    Treal MatrixSymmetric<Treal, Tmatrix>::mixed_diff
    (const MatrixSymmetric<Treal, Tmatrix>& A,
     const MatrixSymmetric<Treal, Tmatrix>& B,
     Treal const requestedAccuracy) {
    MatrixSymmetric<Treal, typename Tmatrix::ElementType> frobNormMat;
    {
      SizesAndBlocks rows_new;
      SizesAndBlocks cols_new;
      A.getSizesAndBlocksForFrobNormMat( rows_new, cols_new );
      frobNormMat.resetSizesAndBlocks(rows_new, cols_new);
      frobNormMat.getMatrix().syAssignDiffFrobNormsLowestLevel(A.getMatrix(),B.getMatrix());
    }
    return frobNormMat.eucl(requestedAccuracy);
  }

#if 1
  template<typename Treal, typename Tmatrix>
    Interval<Treal> MatrixSymmetric<Treal, Tmatrix>::euclDiffIfSmall
    (const MatrixSymmetric<Treal, Tmatrix>& A,
     const MatrixSymmetric<Treal, Tmatrix>& B,
     Treal const requestedAccuracy,
     Treal const maxAbsVal,
     VectorType * const eVecPtr) {
    // DiffMatrix is a lightweight proxy object:
    mat::DiffMatrix< MatrixSymmetric<Treal, Tmatrix>, Treal> Diff(A, B);
    // Note that this function lies outside this class 
    Treal relTol = template_blas_sqrt(template_blas_sqrt(mat::getMachineEpsilon<Treal>()));
    Interval<Treal> tmpInterval = mat::euclIfSmall(Diff, requestedAccuracy, relTol, maxAbsVal, eVecPtr);
    // Emanuel note: Ugly fix to make certain tests pass, we expand
    // the interval up to the requested accuracy. Note that larger
    // intervals may occur if the norm is not 'small'. It happens that
    // Lanczos misconverges to for example the second largest
    // eigenvalue. This happens in particular when the first and second
    // eigenvalues are very close (of the order of the requested
    // accuracy). Expanding the interval makes the largest eigenvalue
    // (at least for certain cases) end up inside the interval even
    // though Lanczos has misconverged.
    if ( tmpInterval.length() < 2*requestedAccuracy )
      return Interval<Treal>( tmpInterval.midPoint()-requestedAccuracy, tmpInterval.midPoint()+requestedAccuracy );
    return tmpInterval;
  }
  
#endif


  template<typename Treal, typename Tmatrix>
    Treal MatrixSymmetric<Treal, Tmatrix>::
    thresh(Treal const threshold,
	   normType const norm) {
    switch (norm) {
    case frobNorm:
      return this->frob_thresh(threshold);
      break;
    case euclNorm:
      return this->eucl_thresh(threshold);
      break;
    case mixedNorm:
      return this->mixed_norm_thresh(threshold);
      break;
    default:
      throw Failure("MatrixSymmetric<Treal, Tmatrix>::"
		    "thresh(Treal const, "
		    "normType const): "
		    "Thresholding not implemented for selected norm");
    }
  }
  
#if 1

  template<typename Treal, typename Tmatrix>
    Treal MatrixSymmetric<Treal, Tmatrix>::
    eucl_thresh(Treal const threshold,
		MatrixTriangular<Treal, Tmatrix> const * const Zptr) {    
    if ( Zptr == NULL ) {
      EuclTruncationSymm<MatrixSymmetric<Treal, Tmatrix>, Treal> TruncObj(*this);
      return TruncObj.run( threshold );
    }
    EuclTruncationSymmWithZ<MatrixSymmetric<Treal, Tmatrix>, MatrixTriangular<Treal, Tmatrix>, Treal> TruncObj(*this, *Zptr);
    return TruncObj.run( threshold );
  }
  
#endif  


  template<typename Treal, typename Tmatrix>
    void MatrixSymmetric<Treal, Tmatrix>::getSizesAndBlocksForFrobNormMat
    ( SizesAndBlocks & rows_new, SizesAndBlocks & cols_new ) const {
    std::vector<int> rows_block_sizes;
    std::vector<int> cols_block_sizes;
    int n_rows;
    int n_cols;
    {
      SizesAndBlocks rows;
      SizesAndBlocks cols;
      this->getRows(rows);
      this->getCols(cols);
      rows.getBlockSizeVector( rows_block_sizes );
      cols.getBlockSizeVector( cols_block_sizes );
      rows_block_sizes.pop_back(); // Remove the '1' at the end
      cols_block_sizes.pop_back(); // Remove the '1' at the end
      n_rows = rows.getNTotalScalars();
      n_cols = cols.getNTotalScalars();
      int factor_rows = rows_block_sizes[rows_block_sizes.size()-1]; 
      int factor_cols = cols_block_sizes[cols_block_sizes.size()-1]; 
      for (unsigned int ind = 0; ind < rows_block_sizes.size(); ++ind) 
	rows_block_sizes[ind] = rows_block_sizes[ind] / factor_rows;
      for (unsigned int ind = 0; ind < cols_block_sizes.size(); ++ind) 
	cols_block_sizes[ind] = cols_block_sizes[ind] / factor_cols;
      // Now set the number of (scalar) rows and cols, should be equal
      // to the number of blocks at the lowest level of the original
      // matrix
      if (n_rows % factor_rows) 
	n_rows = n_rows / factor_rows + 1; 
      else
	n_rows = n_rows / factor_rows;
      if (n_cols % factor_cols) 
	n_cols = n_cols / factor_cols + 1; 
      else
	n_cols = n_cols / factor_cols;
    }
    rows_new = SizesAndBlocks( rows_block_sizes, n_rows );
    cols_new = SizesAndBlocks( cols_block_sizes, n_cols );
  }


  template<typename Treal, typename Tmatrix>
    Treal MatrixSymmetric<Treal, Tmatrix>::
    mixed_norm_thresh(Treal const threshold) {
    assert(threshold >= (Treal)0.0);
    if (threshold == (Treal)0.0)
      return (Treal)0;
    // Construct SizesAndBlocks for frobNormMat
    SizesAndBlocks rows_new;
    SizesAndBlocks cols_new;
    this->getSizesAndBlocksForFrobNormMat( rows_new, cols_new );
    // Now we can construct an empty matrix where the Frobenius norms
    // of lowest level nonzero submatrices will be stored
    MatrixSymmetric<Treal, typename Tmatrix::ElementType> frobNormMat;
    frobNormMat.resetSizesAndBlocks(rows_new, cols_new);
    // We want the following step to dominate the mixed_norm truncation (this
    // is where Frobenius norms of submatrices are computed, i.e. it
    // is here we loop over all matrix elements.)
    frobNormMat.getMatrix().syAssignFrobNormsLowestLevel(this->getMatrix());
    EuclTruncationSymmElementLevel<MatrixSymmetric<Treal, typename Tmatrix::ElementType>, Treal> TruncObj( frobNormMat );
    Treal mixed_norm_result = TruncObj.run( threshold );
    frobNormMat.getMatrix().truncateAccordingToSparsityPattern(this->getMatrix());
    return mixed_norm_result;
  }


  /* B = alpha * A */
  template<typename Treal, typename Tmatrix>
    inline MatrixSymmetric<Treal, Tmatrix>& 
    MatrixSymmetric<Treal, Tmatrix>::operator= 
    (XY<Treal, MatrixSymmetric<Treal, Tmatrix> > const & sm) {
  
  if(this == &sm.B) // A = B
      {
        *this *= sm.A; // B *= alpha
        return *this;
      }
    assert(!sm.tB);
    /* Data structure set by assign - therefore set haveDataStructure to true */
    this->matrixPtr.haveDataStructureSet(true); 
    this->matrixPtr->assign(sm.A, *sm.B.matrixPtr);
    return *this;
  }
  /* C = alpha * A * A + beta * C          : A and C are symmetric    */
  template<typename Treal, typename Tmatrix>
    MatrixSymmetric<Treal, Tmatrix>& 
    MatrixSymmetric<Treal, Tmatrix>::operator=
    (const XYZpUV<Treal, 
     MatrixSymmetric<Treal, Tmatrix>, 
     MatrixSymmetric<Treal, Tmatrix>,
     Treal, 
     MatrixSymmetric<Treal, Tmatrix> >& sm2psm) {
    assert(this != &sm2psm.B);
    if (this == &sm2psm.E && &sm2psm.B == &sm2psm.C) { 
      /* Operation is C = alpha * A * A + beta * C */
      Tmatrix::sysq('U', 
		    sm2psm.A, *sm2psm.B.matrixPtr, 
		    sm2psm.D, *this->matrixPtr);
      return *this;
    }
    else /* this != &sm2psm.C || &sm2psm.B != &sm2psm.C */
      throw Failure("MatrixSymmetric<Treal, Tmatrix>::operator="
		    "(const XYZpUV<Treal, MatrixSymmetric"
		    "<Treal, Tmatrix> >& sm2psm) :  "
		    "D = alpha * A * B + beta * C not supported for C != D"
		    " and for symmetric matrices not for A != B since this "
		    "generally will result in a nonsymmetric matrix");
  }

  /* C = alpha * A * A                     : A and C are symmetric    */
  template<typename Treal, typename Tmatrix>
    MatrixSymmetric<Treal, Tmatrix>& 
    MatrixSymmetric<Treal, Tmatrix>::operator=
    (const XYZ<Treal, 
     MatrixSymmetric<Treal, Tmatrix>, 
     MatrixSymmetric<Treal, Tmatrix> >& sm2) {
    assert(this != &sm2.B);
    if (&sm2.B == &sm2.C) { 
      this->matrixPtr.haveDataStructureSet(true); 
      Tmatrix::sysq('U', sm2.A, *sm2.B.matrixPtr, 0, *this->matrixPtr);
      return *this;
    }
    else
      throw Failure("MatrixSymmetric<Treal, Tmatrix>::operator="
		    "(const XYZ<Treal, MatrixSymmetric<Treal, Tmatrix>,"
		    " MatrixSymmetric<Treal, Tmatrix> >& sm2) :  "
		    "Operation C = alpha * A * B with only symmetric "
		    "matrices not supported for A != B");
  }
  
  /* C += alpha * A * A                    : A and C are symmetric    */
  template<typename Treal, typename Tmatrix>
    MatrixSymmetric<Treal, Tmatrix>& 
    MatrixSymmetric<Treal, Tmatrix>::operator+=
    (const XYZ<Treal, 
     MatrixSymmetric<Treal, Tmatrix>, 
     MatrixSymmetric<Treal, Tmatrix> >& sm2) {
    assert(this != &sm2.B);
    if (&sm2.B == &sm2.C) { 
      Tmatrix::sysq('U', sm2.A, *sm2.B.matrixPtr, 1, *this->matrixPtr);
      return *this;
    }
    else
      throw Failure("MatrixSymmetric<Treal, Tmatrix>::operator+="
		    "(const XYZ<Treal, MatrixSymmetric<Treal, Tmatrix>,"
		    " MatrixSymmetric<Treal, Tmatrix> >& sm2) :  "
		    "Operation C += alpha * A * B with only symmetric "
		    "matrices not supported for A != B");
  }

  /* C = alpha * A * transpose(A) + beta * C        : C is symmetric  */
  template<typename Treal, typename Tmatrix>
    MatrixSymmetric<Treal, Tmatrix>& 
    MatrixSymmetric<Treal, Tmatrix>::operator= 
    (const XYZpUV<Treal, 
     MatrixGeneral<Treal, Tmatrix>, 
     MatrixGeneral<Treal, Tmatrix>, 
     Treal,
     MatrixSymmetric<Treal, Tmatrix> >& smmpsm) {
    if (this == &smmpsm.E)
      if (&smmpsm.B == &smmpsm.C)
	if (!smmpsm.tB && smmpsm.tC) {
	  Tmatrix::syrk('U', false, 
			smmpsm.A, *smmpsm.B.matrixPtr, 
			smmpsm.D, *this->matrixPtr);
	}
	else 
	  throw Failure("MatrixSymmetric<Treal, Tmatrix>::operator="
			"(const XYZpUV<Treal, MatrixGeneral"
			"<Treal, Tmatrix>, "
			"MatrixGeneral<Treal, Tmatrix>, Treal, "
			"MatrixSymmetric<Treal, Tmatrix> >&) : "
			"C = alpha * A' * A + beta * C, not implemented"
			" only C = alpha * A * A' + beta * C");
      else
	throw Failure("MatrixSymmetric<Treal, Tmatrix>::operator="
		      "(const XYZpUV<"
		      "Treal, MatrixGeneral<Treal, Tmatrix>, "
		      "MatrixGeneral<Treal, Tmatrix>, Treal, "
		      "MatrixSymmetric<Treal, Tmatrix> >&) : "
		      "You are trying to call C = alpha * A * A' + beta * C"
		      " with  A and A' being different objects");
    else
      throw Failure("MatrixSymmetric<Treal, Tmatrix>::operator="
		    "(const XYZpUV<"
		    "Treal, MatrixGeneral<Treal, Tmatrix>, "
		    "MatrixGeneral<Treal, Tmatrix>, Treal, "
		    "MatrixSymmetric<Treal, Tmatrix> >&) :  "
		    "D = alpha * A * A' + beta * C not supported for C != D");
    return *this;
  }

  /* C = alpha * A * transpose(A)                   : C is symmetric  */
  template<typename Treal, typename Tmatrix>
    MatrixSymmetric<Treal, Tmatrix>& 
    MatrixSymmetric<Treal, Tmatrix>::operator= 
    (const XYZ<Treal, 
     MatrixGeneral<Treal, Tmatrix>, 
     MatrixGeneral<Treal, Tmatrix> >& smm) {
    if (&smm.B == &smm.C)
      if (!smm.tB && smm.tC) {
	Tmatrix::syrk('U', false, 
		      smm.A, *smm.B.matrixPtr, 
		      0, *this->matrixPtr);
      }
      else 
	throw Failure("MatrixSymmetric<Treal, Tmatrix>::operator="
		      "(const XYZ<"
		      "Treal, MatrixGeneral<Treal, Tmatrix>, "
		      "MatrixGeneral<Treal, Tmatrix> >&) : "
		      "C = alpha * A' * A, not implemented "
		      "only C = alpha * A * A'");
    else
      throw Failure("MatrixSymmetric<Treal, Tmatrix>::operator="
		    "(const XYZ<"
		    "Treal, MatrixGeneral<Treal, Tmatrix>, "
		    "MatrixGeneral<Treal, Tmatrix> >&) : "
		    "You are trying to call C = alpha * A * A' "
		    "with A and A' being different objects");
    return *this;
  }
  
  /* C += alpha * A * transpose(A)                   : C is symmetric  */
  template<typename Treal, typename Tmatrix>
    MatrixSymmetric<Treal, Tmatrix>& 
    MatrixSymmetric<Treal, Tmatrix>::operator+= 
    (const XYZ<Treal, 
     MatrixGeneral<Treal, Tmatrix>, 
     MatrixGeneral<Treal, Tmatrix> >& smm) {
    if (&smm.B == &smm.C)
      if (!smm.tB && smm.tC) {
	Tmatrix::syrk('U', false, 
		      smm.A, *smm.B.matrixPtr, 
		      1, *this->matrixPtr);
      }
      else 
	throw Failure("MatrixSymmetric<Treal, Tmatrix>::operator+="
		      "(const XYZ<"
		      "Treal, MatrixGeneral<Treal, Tmatrix>, "
		      "MatrixGeneral<Treal, Tmatrix> >&) : "
		      "C += alpha * A' * A, not implemented "
		      "only C += alpha * A * A'");
    else
      throw Failure("MatrixSymmetric<Treal, Tmatrix>::operator+="
		    "(const XYZ<"
		    "Treal, MatrixGeneral<Treal, Tmatrix>, "
		    "MatrixGeneral<Treal, Tmatrix> >&) : "
		    "You are trying to call C += alpha * A * A' "
		    "with  A and A' being different objects");
    return *this;
  }
  
#if 1
  /* A = op1(Z) * A * op2(Z)   : Z is upper triangular and A is symmetric */
  /* Either op1() or op2() is the transpose operation. */
  template<typename Treal, typename Tmatrix>
    MatrixSymmetric<Treal, Tmatrix>& 
    MatrixSymmetric<Treal, Tmatrix>::operator= 
    (const XYZ<MatrixTriangular<Treal, Tmatrix>, 
     MatrixSymmetric<Treal, Tmatrix>, 
     MatrixTriangular<Treal, Tmatrix> >& zaz) {
    if (this == &zaz.B) {
      if (&zaz.A == &zaz.C) {
	if (zaz.tA && !zaz.tC) {
	  Tmatrix::trsytriplemm('R', *zaz.A.matrixPtr, *this->matrixPtr); 
	}
	else if (!zaz.tA && zaz.tC) {
	  Tmatrix::trsytriplemm('L', *zaz.A.matrixPtr, *this->matrixPtr); 
	}
	else
	  throw Failure("MatrixSymmetric<Treal, Tmatrix>::operator=" 
			"(const XYZ<MatrixTriangular<Treal, Tmatrix>," 
			"MatrixSymmetric<Treal, Tmatrix>," 
			"MatrixTriangular<Treal, Tmatrix> >&) : "
			"A = op1(Z) * A * op2(Z) : Either op1 xor op2 must be "
			"the transpose operation.");
      }
      else
	throw Failure("MatrixSymmetric<Treal, Tmatrix>::operator=" 
		      "(const XYZ<MatrixTriangular<Treal, Tmatrix>," 
		      "MatrixSymmetric<Treal, Tmatrix>," 
		      "MatrixTriangular<Treal, Tmatrix> >&) : "
		      "A = op1(Z1) * A * op2(Z2) : Z1 and Z2 must be the same "
		      "object");
    }
    else
      throw Failure("MatrixSymmetric<Treal, Tmatrix>::operator=" 
		    "(const XYZ<MatrixTriangular<Treal, Tmatrix>," 
		    "MatrixSymmetric<Treal, Tmatrix>," 
		    "MatrixTriangular<Treal, Tmatrix> >&) : "
		    "C = op1(Z) * A * op2(Z) : A and C must be the same "
		    "object");
    return *this;
  }

#endif


  /** C = alpha * A * B + beta * C where A and B are symmetric
   *  and only the upper triangle of C is computed, 
   *  C is enforced to be symmetric!
   */
  template<typename Treal, typename Tmatrix>
    void MatrixSymmetric<Treal, Tmatrix>::
    ssmmUpperTriangleOnly(const Treal alpha, 
			  const MatrixSymmetric<Treal, Tmatrix>& A, 
			  const MatrixSymmetric<Treal, Tmatrix>& B, 
			  const Treal beta, 
			  MatrixSymmetric<Treal, Tmatrix>& C) {
    Tmatrix::ssmm_upper_tr_only(alpha, *A.matrixPtr, *B.matrixPtr,
				beta, *C.matrixPtr);
  }

  



  /* Addition */  
  /* C =  A + B   */
  template<typename Treal, typename Tmatrix>
    inline MatrixSymmetric<Treal, Tmatrix>& 
    MatrixSymmetric<Treal, Tmatrix>::operator= 
    (const XpY<MatrixSymmetric<Treal, Tmatrix>,
     MatrixSymmetric<Treal, Tmatrix> >& mpm) {
    assert(this != &mpm.A);
    (*this) = mpm.B;
    Tmatrix::add(1.0, *mpm.A.matrixPtr, *this->matrixPtr);
    return *this;
  }
  /* C =  A - B   */
  template<typename Treal, typename Tmatrix>
    inline MatrixSymmetric<Treal, Tmatrix>& 
    MatrixSymmetric<Treal, Tmatrix>::operator= 
    (const XmY<MatrixSymmetric<Treal, Tmatrix>,
     MatrixSymmetric<Treal, Tmatrix> >& mmm) {
    assert(this != &mmm.B);
    (*this) = mmm.A;
    Tmatrix::add(-1.0, *mmm.B.matrixPtr, *this->matrixPtr);
    return *this;
  }
  /* B += A */
  template<typename Treal, typename Tmatrix>
    inline MatrixSymmetric<Treal, Tmatrix>& 
    MatrixSymmetric<Treal, Tmatrix>::operator+= 
    (MatrixSymmetric<Treal, Tmatrix> const & A) {
    Tmatrix::add(1.0, *A.matrixPtr, *this->matrixPtr);
    return *this;
  }
  /* B -= A */
  template<typename Treal, typename Tmatrix>
    inline MatrixSymmetric<Treal, Tmatrix>& 
    MatrixSymmetric<Treal, Tmatrix>::operator-=
    (MatrixSymmetric<Treal, Tmatrix> const & A) {
    Tmatrix::add(-1.0, *A.matrixPtr, *this->matrixPtr);
    return *this;
  }
  /* B += alpha * A */
  template<typename Treal, typename Tmatrix>
    inline MatrixSymmetric<Treal, Tmatrix>& 
    MatrixSymmetric<Treal, Tmatrix>::operator+= 
    (XY<Treal, MatrixSymmetric<Treal, Tmatrix> > const & sm) {
    assert(!sm.tB);
    Tmatrix::add(sm.A, *sm.B.matrixPtr, *this->matrixPtr);
    return *this;
  }

  /* B -= alpha * A */
  template<typename Treal, typename Tmatrix>
    inline MatrixSymmetric<Treal, Tmatrix>& 
    MatrixSymmetric<Treal, Tmatrix>::operator-= 
    (XY<Treal, MatrixSymmetric<Treal, Tmatrix> > const & sm) {
    assert(!sm.tB);
    Tmatrix::add(-sm.A, *sm.B.matrixPtr, *this->matrixPtr);
    return *this;
  }

  /** Performs operation specified in 'op' on all nonzero matrix elements 
   *  and sums up the result and returns it.
   * 
   */
  template<typename Treal, typename Tmatrix, typename Top>
    Treal accumulate(MatrixSymmetric<Treal, Tmatrix> & A, 
		     Top & op) {
    return A.accumulateWith(op);
  }


  
} /* end namespace mat */
#endif

