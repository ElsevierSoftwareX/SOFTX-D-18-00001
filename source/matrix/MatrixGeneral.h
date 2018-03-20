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

/** @file MatrixGeneral.h General matrix class
 *
 * Copyright(c) Emanuel Rubensson 2006
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date January 2006
 *
 */
#ifndef MAT_MATRIXGENERAL
#define MAT_MATRIXGENERAL
#include "MatrixBase.h"
namespace mat {
  /** Normal matrix 
   *
   * This class belongs to the matrix API
   *
   *  
   * Treal: Type for real numbers
   *
   * Tmatrix: The matrix class
   *
   * Tperm: Permutation used in the matrix class  
   *
   * @see MatrixBase
   * @see MatrixSymmetric
   * @see MatrixTriangular
   * 
   */  
  template<typename Treal, typename Tmatrix>
    class MatrixGeneral : public MatrixBase<Treal, Tmatrix> {
  public:
    typedef VectorGeneral<Treal, typename Tmatrix::VectorType> VectorType;

    MatrixGeneral()
      :MatrixBase<Treal, Tmatrix>() {} /**< Default constructor  */
    explicit MatrixGeneral(const MatrixGeneral<Treal, Tmatrix>& matr)
      :MatrixBase<Treal, Tmatrix>(matr) {} /**< Copy constructor  */
    explicit MatrixGeneral(const MatrixSymmetric<Treal, Tmatrix>& symm)
      :MatrixBase<Treal, Tmatrix>(symm) { 
      this->matrixPtr->symToNosym();
    }/**< Copy from symmetric matrix constructor  */
    explicit MatrixGeneral(const MatrixTriangular<Treal, Tmatrix>& triang)
      :MatrixBase<Treal, Tmatrix>(triang) {}
    /**< Copy from triangular matrix constructor  */

#if 0
    template<typename Tfull>
      inline void assign_from_full
      (Tfull const* const fullmatrix, 
       const int totnrows, const int totncols) {
      this->matrixPtr->assign_from_full(fullmatrix, totnrows, totncols);      
    }    
    inline void assign_from_full
      (Treal const* const fullmatrix, 
       const int totnrows, const int totncols) {
      this->matrixPtr->assign_from_full(fullmatrix, totnrows, totncols);      
    }
#endif

    inline void assignFromFull
      (std::vector<Treal> const & fullMat) {
      assert((int)fullMat.size() == this->get_nrows() * this->get_ncols());
      this->matrixPtr->assignFromFull(fullMat);      
    }

    inline void fullMatrix(std::vector<Treal> & fullMat) const {
      this->matrixPtr->fullMatrix(fullMat);
    }    

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
      for (unsigned int ind = 0; ind < values.size(); ++ind) 
	fullMat[rowind[ind] + this->get_nrows() * colind[ind]] =
	  values[ind];
    }
    /**< Save matrix as full matrix.
     * Whole matrix is written in columnwise order.
     * Both lower and upper triangle.
     * Permutation is used.
     */


    inline void assign_from_sparse
      (std::vector<int> const & rowind, 
       std::vector<int> const & colind, 
       std::vector<Treal> const & values, 
       SizesAndBlocks const & newRows,
       SizesAndBlocks const & newCols) {
      this->resetSizesAndBlocks(newRows, newCols);
      this->matrixPtr->assignFromSparse(rowind, colind, values);
    }
    /**< Assign from sparse matrix given by three arrays. 
     * The arrays contain row indices, column indices and values.
     * The indices start at zero.
     * nval is the length of the three arrays.
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
      MatrixBase<Treal, Tmatrix>::
	getPermutedIndexes(rowind, rowPermutation, newRowind);
      MatrixBase<Treal, Tmatrix>::
	getPermutedIndexes(colind, colPermutation, newColind);
      this->matrixPtr->assignFromSparse(newRowind, newColind, values);
    }
    /**< Same as above, except taking two additional arguments 
     *   specifying the permutation of rows and columns.
     *   Also assuming that sizes and blocks are already known
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
    /**< Same as above, except not assuming that sizes and blocks are set.
     */


    inline void get_values
      (std::vector<int> const & rowind, 
       std::vector<int> const & colind, 
       std::vector<Treal> & values) const {
      this->matrixPtr->getValues(rowind, colind, values);
    }
    /**< Get values given by row and column index lists.
     * Input arrays contain row and column indices.
     * The output array contains values for the given indices.
     * nval is the length of the three arrays.
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
      MatrixBase<Treal, Tmatrix>::
	getPermutedIndexes(rowind, rowPermutation, newRowind);
      MatrixBase<Treal, Tmatrix>::
	getPermutedIndexes(colind, colPermutation, newColind);
      this->matrixPtr->getValues(newRowind, newColind, values);
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
      this->matrixPtr->getAllValues(rowind, colind, values);
    }
    /**< Get all values and corresponding row and column index lists,
     * in matrix.
     * nval is the length of the three arrays and is preferably 
     * computed with nvalues() before hand.
     * Returns the number of values.
     * @see nvalues()
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
      this->matrixPtr->getAllValues(tmpRowind, tmpColind, values);
      MatrixBase<Treal, Tmatrix>::
	getPermutedIndexes(tmpRowind, rowInversePermutation, rowind);
      MatrixBase<Treal, Tmatrix>::
	getPermutedIndexes(tmpColind, colInversePermutation, colind);
      
    }
    /**< Same as above, except taking two additional arguments 
     *   specifying the permutation of rows and columns.
     *   Note, however, that this permutation is the inverse 
     *   permutation compared to the permutations provided in the
     *   functions "assign_from_sparse", "add_values", and "get_values"
     *   @warning permutation is inverse compared to other functions
     */    


#if 0
    inline void fullmatrix(Treal* const full, 
			   const int totnrows, const int totncols) const {
      this->matrixPtr->fullmatrix(full, totnrows, totncols);
    }
#endif
    MatrixGeneral<Treal, Tmatrix>& 
      operator=(const MatrixGeneral<Treal, Tmatrix>& mat) {
      MatrixBase<Treal, Tmatrix>::operator=(mat);
      return *this;
    } 
    inline MatrixGeneral<Treal, Tmatrix>& 
      operator=(const Xtrans<MatrixGeneral<Treal, Tmatrix> >& mt) {
      if (mt.tA)
	MatrixBase<Treal, Tmatrix>::operator=(transpose(mt.A));
      else
	MatrixBase<Treal, Tmatrix>::operator=(mt.A);
      return *this;
    }

    MatrixGeneral<Treal, Tmatrix>& 
      operator=(const MatrixSymmetric<Treal, Tmatrix>& symm) {
      MatrixBase<Treal, Tmatrix>::operator=(symm);
      this->matrixPtr->symToNosym();
      return *this;
    } 
    MatrixGeneral<Treal, Tmatrix>& 
      operator=(const MatrixTriangular<Treal, Tmatrix>& triang) {
      MatrixBase<Treal, Tmatrix>::operator=(triang);
      return *this;
    } 

    inline MatrixGeneral<Treal, Tmatrix>& operator=(int const k) {
      *this->matrixPtr = k;
      return *this;
    }
    inline Treal frob() const {
      return this->matrixPtr->frob();
    }
    static inline Treal frob_diff
      (const MatrixGeneral<Treal, Tmatrix>& A,
       const MatrixGeneral<Treal, Tmatrix>& B) {
      return Tmatrix::frobDiff(*A.matrixPtr, *B.matrixPtr);
    }
    Treal eucl(Treal const requestedAccuracy,
	       int maxIter = -1) const;


    void thresh(Treal const threshold,
		normType const norm);

    inline void frob_thresh(Treal threshold) {
      this->matrixPtr->frob_thresh(threshold);
    }
    Treal eucl_thresh(Treal const threshold);
    
    inline void gershgorin(Treal& lmin, Treal& lmax) {
      this->matrixPtr->gershgorin(lmin, lmax);
    }    
    static inline Treal trace_ab
      (const MatrixGeneral<Treal, Tmatrix>& A,
       const MatrixGeneral<Treal, Tmatrix>& B) {
      return Tmatrix::trace_ab(*A.matrixPtr, *B.matrixPtr);
    }
    static inline Treal trace_aTb
      (const MatrixGeneral<Treal, Tmatrix>& A,
       const MatrixGeneral<Treal, Tmatrix>& B) {
      return Tmatrix::trace_aTb(*A.matrixPtr, *B.matrixPtr);
    }
    inline size_t nnz() const {  /* Note: size_t instead of int here to avoid integer overflow. */
      return this->matrixPtr->nnz();
    }
    inline size_t nvalues() const { /* Note: size_t instead of int here to avoid integer overflow. */
      return this->matrixPtr->nvalues();
    }

    inline void write_to_buffer(void* buffer, const int n_bytes) const {
      this->write_to_buffer_base(buffer, n_bytes, matrix_matr);
    }
    inline void read_from_buffer(void* buffer, const int n_bytes) {
      this->read_from_buffer_base(buffer, n_bytes, matrix_matr);
    }
    
    /* OPERATIONS ONLY INVOLVING ORDINARY MATRICES */
    /** C = alpha * op(A) * op(B)  */
    MatrixGeneral<Treal, Tmatrix>& operator= 
      (const XYZ<Treal, 
       MatrixGeneral<Treal, Tmatrix>, 
       MatrixGeneral<Treal, Tmatrix> >& smm);

    /** C = op(A) * op(B)  */
    MatrixGeneral<Treal, Tmatrix>& operator= 
      (const XY<MatrixGeneral<Treal, Tmatrix>, 
       MatrixGeneral<Treal, Tmatrix> >& mm);

    /** C += alpha * op(A) * op(B) */
    MatrixGeneral<Treal, Tmatrix>& operator+=   
      (const XYZ<Treal, 
       MatrixGeneral<Treal, Tmatrix>, 
       MatrixGeneral<Treal, Tmatrix> >& smm);

    /** C = alpha * op(A) * op(B) + beta * C  */
    MatrixGeneral<Treal, Tmatrix>& operator= 
      (const XYZpUV<Treal, 
       MatrixGeneral<Treal, Tmatrix>,
       MatrixGeneral<Treal, Tmatrix>,
       Treal, 
       MatrixGeneral<Treal, Tmatrix> >& smmpsm);

    /** C =  A + B */
    MatrixGeneral<Treal, Tmatrix>& operator= 
      (XpY<MatrixGeneral<Treal, Tmatrix>,
       MatrixGeneral<Treal, Tmatrix> > const & mpm);
    /** B += A */
    MatrixGeneral<Treal, Tmatrix>& operator+= 
      (MatrixGeneral<Treal, Tmatrix> const & A);
    MatrixGeneral<Treal, Tmatrix>& operator-= 
      (MatrixGeneral<Treal, Tmatrix> const & A);
    /** B += alpha * A */
    MatrixGeneral<Treal, Tmatrix>& operator+= 
      (XY<Treal, MatrixGeneral<Treal, Tmatrix> > const & sm);
        

    /* OPERATIONS INVOLVING SYMMETRIC MATRICES */
    /** C = alpha * A * B                      : A is symmetric */
    MatrixGeneral<Treal, Tmatrix>& operator= 
      (const XYZ<Treal, 
       MatrixSymmetric<Treal, Tmatrix>, 
       MatrixGeneral<Treal, Tmatrix> >& smm);
    /** C = A * B                      : A is symmetric */
    MatrixGeneral<Treal, Tmatrix>& operator= 
      (const XY<MatrixSymmetric<Treal, Tmatrix>, 
       MatrixGeneral<Treal, Tmatrix> >& mm);
    /** C += alpha * A * B                     : A is symmetric */
    MatrixGeneral<Treal, Tmatrix>& operator+=   
      (const XYZ<Treal, 
       MatrixSymmetric<Treal, Tmatrix>, 
       MatrixGeneral<Treal, Tmatrix> >& smm);
    /** C = alpha * A * B + beta * C           : A is symmetric */
    MatrixGeneral<Treal, Tmatrix>& operator= 
      (const XYZpUV<Treal, 
       MatrixSymmetric<Treal, Tmatrix>, 
       MatrixGeneral<Treal, Tmatrix>, 
       Treal,
       MatrixGeneral<Treal, Tmatrix> >& smmpsm);
    /** C = alpha * B * A                      : A is symmetric */
    MatrixGeneral<Treal, Tmatrix>& operator= 
      (const XYZ<Treal, 
       MatrixGeneral<Treal, Tmatrix>,
       MatrixSymmetric<Treal, Tmatrix> >& smm);
    /** C = B * A                      : A is symmetric */
    MatrixGeneral<Treal, Tmatrix>& operator= 
      (const XY<MatrixGeneral<Treal, Tmatrix>,
       MatrixSymmetric<Treal, Tmatrix> >& mm);
    /** C += alpha * B * A                     : A is symmetric */
    MatrixGeneral<Treal, Tmatrix>& operator+=   
      (const XYZ<Treal, 
       MatrixGeneral<Treal, Tmatrix>, 
       MatrixSymmetric<Treal, Tmatrix> >& smm);
    /** C = alpha * B * A + beta * C           : A is symmetric */
    MatrixGeneral<Treal, Tmatrix>& operator= 
      (const XYZpUV<Treal, 
       MatrixGeneral<Treal, Tmatrix>, 
       MatrixSymmetric<Treal, Tmatrix>, 
       Treal,
       MatrixGeneral<Treal, Tmatrix> >& smmpsm);
    /** C = alpha * A * B                      : A and B are symmetric */
    MatrixGeneral<Treal, Tmatrix>& operator= 
      (const XYZ<Treal, 
       MatrixSymmetric<Treal, Tmatrix>, 
       MatrixSymmetric<Treal, Tmatrix> >& smm);
    /** C = A * B                      : A and B are symmetric */
    MatrixGeneral<Treal, Tmatrix>& operator= 
      (const XY<MatrixSymmetric<Treal, Tmatrix>, 
       MatrixSymmetric<Treal, Tmatrix> >& mm);
    /** C += alpha * A * B                     : A and B are symmetric */
    MatrixGeneral<Treal, Tmatrix>& operator+= 
      (const XYZ<Treal, 
       MatrixSymmetric<Treal, Tmatrix>, 
       MatrixSymmetric<Treal, Tmatrix> >& smm);
    /** C = alpha * A * B + beta * C           : A and B are symmetric */
    MatrixGeneral<Treal, Tmatrix>& operator= 
      (const XYZpUV<Treal, 
       MatrixSymmetric<Treal, Tmatrix>, 
       MatrixSymmetric<Treal, Tmatrix>, 
       Treal,
       MatrixGeneral<Treal, Tmatrix> >& smmpsm);
    
    /* OPERATIONS INVOLVING UPPER TRIANGULAR MATRICES */    
    /** B = alpha * op(A) * B                  : A is upper triangular */
    MatrixGeneral<Treal, Tmatrix>& operator= 
      (const XYZ<Treal, 
       MatrixTriangular<Treal, Tmatrix>, 
       MatrixGeneral<Treal, Tmatrix> >& smm);
    /** B = alpha * B * op(A)                  : A is upper triangular */
    MatrixGeneral<Treal, Tmatrix>& operator= 
      (const XYZ<Treal, 
       MatrixGeneral<Treal, Tmatrix>, 
       MatrixTriangular<Treal, Tmatrix> >& smm);
    
    void random() {
      this->matrixPtr->random();
    }

    void randomZeroStructure(Treal probabilityBeingZero) {
      this->matrixPtr->randomZeroStructure(probabilityBeingZero);
    }

    template<typename TRule>
      void setElementsByRule(TRule & rule) {
      this->matrixPtr->setElementsByRule(rule);
      return;
    }

    std::string obj_type_id() const {return "MatrixGeneral";}
    protected:
    inline void writeToFileProt(std::ofstream & file) const {
      this->writeToFileBase(file, matrix_matr);
    }
    inline void readFromFileProt(std::ifstream & file) {
      this->readFromFileBase(file, matrix_matr);
    }
    private:

  };

  template<typename Treal, typename Tmatrix>
    Treal MatrixGeneral<Treal, Tmatrix>::
    eucl(Treal const requestedAccuracy,
	 int maxIter) const {
    VectorType guess;
    SizesAndBlocks cols;
    this->getCols(cols);
    guess.resetSizesAndBlocks(cols);
    guess.rand();
    mat::ATAMatrix<MatrixGeneral<Treal, Tmatrix>, Treal> ata(*this);    
    if (maxIter < 0)
      maxIter = this->get_nrows() * 100;
    arn::LanczosLargestMagnitudeEig
      <Treal, ATAMatrix<MatrixGeneral<Treal, Tmatrix>, Treal>, VectorType>
      lan(ata, guess, maxIter);
    lan.setRelTol( 100 * mat::getMachineEpsilon<Treal>() );
    lan.run();
    Treal eVal = 0;
    Treal acc = 0;
    lan.getLargestMagnitudeEig(eVal, acc);
    Interval<Treal> euclInt( sqrt(eVal-acc),
			     sqrt(eVal+acc) );    
    if ( euclInt.low() < 0 )
      euclInt = Interval<Treal>( 0, sqrt(eVal+acc) );
    if ( euclInt.length() / 2.0 > requestedAccuracy ) {
      std::cout << "req: " << requestedAccuracy
		<< "  obt: " << euclInt.length() / 2.0 << std::endl;
      throw std::runtime_error("Desired accuracy not obtained in Lanczos.");
    }
    return euclInt.midPoint();
  }


  template<typename Treal, typename Tmatrix>
    void MatrixGeneral<Treal, Tmatrix>::
    thresh(Treal const threshold,
	   normType const norm) {
    switch (norm) {
    case frobNorm:
      this->frob_thresh(threshold);
      break;
    default:
      throw Failure("MatrixGeneral<Treal, Tmatrix>::"
		    "thresh(Treal const, "
		    "normType const): "
		    "Thresholding not imlpemented for selected norm");
    }
  }

  template<typename Treal, typename Tmatrix>
    Treal MatrixGeneral<Treal, Tmatrix>::
    eucl_thresh(Treal const threshold) {
    EuclTruncationGeneral<MatrixGeneral<Treal, Tmatrix>, Treal> TruncObj( *this );
    return TruncObj.run( threshold );    
  }


  /* OPERATIONS ONLY INVOLVING ORDINARY MATRICES */
    /* C = alpha * op(A) * op(B)  */
  template<typename Treal, typename Tmatrix>
    MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator= 
    (const XYZ<Treal, 
     MatrixGeneral<Treal, Tmatrix>, 
     MatrixGeneral<Treal, Tmatrix> >& smm) {
    assert(this != &smm.B && this != &smm.C);
    this->matrixPtr.haveDataStructureSet(true);
    Tmatrix::gemm(smm.tB, smm.tC, smm.A, 
		  *smm.B.matrixPtr, *smm.C.matrixPtr, 
		  0, *this->matrixPtr);
    return *this;
  }

    /* C = op(A) * op(B)  */
  template<typename Treal, typename Tmatrix>
    MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator= 
    (const XY<MatrixGeneral<Treal, Tmatrix>, 
     MatrixGeneral<Treal, Tmatrix> >& mm) {
    assert(this != &mm.A && this != &mm.B);
    Tmatrix::gemm(mm.tA, mm.tB, 1.0, 
		  *mm.A.matrixPtr, *mm.B.matrixPtr, 
		  0, *this->matrixPtr);
    return *this;
  }
   
  /* C += alpha * op(A) * op(B) */
  template<typename Treal, typename Tmatrix>
    MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator+=   
    (const XYZ<Treal, 
     MatrixGeneral<Treal, Tmatrix>, 
     MatrixGeneral<Treal, Tmatrix> >& smm) {
    assert(this != &smm.B && this != &smm.C);
    Tmatrix::gemm(smm.tB, smm.tC, smm.A, 
		  *smm.B.matrixPtr, *smm.C.matrixPtr, 
		  1, *this->matrixPtr);
    return *this;
  }
  
  /* C = alpha * op(A) * op(B) + beta * C  */
  template<typename Treal, typename Tmatrix>
    MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator= 
    (const XYZpUV<Treal, 
     MatrixGeneral<Treal, Tmatrix>, 
     MatrixGeneral<Treal, Tmatrix>,
     Treal,
     MatrixGeneral<Treal, Tmatrix> >& smmpsm) {
    assert(this != &smmpsm.B && this != &smmpsm.C);
    assert(!smmpsm.tE);
    if (this == &smmpsm.E)
      Tmatrix::gemm(smmpsm.tB, smmpsm.tC, smmpsm.A, 
		    *smmpsm.B.matrixPtr, *smmpsm.C.matrixPtr, 
		    smmpsm.D, *this->matrixPtr);
    else
      throw Failure("MatrixGeneral<Treal, Tmatrix>::operator="
		    "(const XYZpUV<Treal, MatrixGeneral"
		    "<Treal, Tmatrix> >&) :  D = alpha "
		    "* op(A) * op(B) + beta * C not supported for C != D");
    return *this;
  }


  /* C =  A + B   */
  template<typename Treal, typename Tmatrix>
    inline MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator= 
    (const XpY<MatrixGeneral<Treal, Tmatrix>,
     MatrixGeneral<Treal, Tmatrix> >& mpm) {
    assert(this != &mpm.A);
    (*this) = mpm.B;
    Tmatrix::add(1.0, *mpm.A.matrixPtr, *this->matrixPtr);
    return *this;
  }
  /* B += A */
  template<typename Treal, typename Tmatrix>
    inline MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator+= 
    (MatrixGeneral<Treal, Tmatrix> const & A) {
    Tmatrix::add(1.0, *A.matrixPtr, *this->matrixPtr);
    return *this;
  }
  /* B -= A */
  template<typename Treal, typename Tmatrix>
    inline MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator-= 
    (MatrixGeneral<Treal, Tmatrix> const & A) {
    Tmatrix::add(-1.0, *A.matrixPtr, *this->matrixPtr);
    return *this;
  }
  /* B += alpha * A */
  template<typename Treal, typename Tmatrix>
    inline MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator+= 
    (XY<Treal, MatrixGeneral<Treal, Tmatrix> > const & sm) {
    assert(!sm.tB);
    Tmatrix::add(sm.A, *sm.B.matrixPtr, *this->matrixPtr);
    return *this;
  }

    
  /* OPERATIONS INVOLVING SYMMETRIC MATRICES */
  /* C = alpha * A * B                      : A is symmetric */
  template<typename Treal, typename Tmatrix>
    MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator= 
    (const XYZ<Treal, 
     MatrixSymmetric<Treal, Tmatrix>, 
     MatrixGeneral<Treal, Tmatrix> >& smm) {
    assert(this != &smm.C);
    assert(!smm.tB && !smm.tC);
    this->matrixPtr.haveDataStructureSet(true);
    Tmatrix::symm('L', 'U', smm.A, 
		  *smm.B.matrixPtr, *smm.C.matrixPtr,
		  0, *this->matrixPtr);
    return *this;
  }

  /* C = A * B                      : A is symmetric */
  template<typename Treal, typename Tmatrix>
    MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator= 
    (const XY<MatrixSymmetric<Treal, Tmatrix>, 
     MatrixGeneral<Treal, Tmatrix> >& mm) {
    assert(this != &mm.B);
    assert(!mm.tB);
    Tmatrix::symm('L', 'U', 1.0, 
		  *mm.A.matrixPtr, *mm.B.matrixPtr,
		  0, *this->matrixPtr);
    return *this;
  }
  
  /* C += alpha * A * B                     : A is symmetric */
  template<typename Treal, typename Tmatrix>
    MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator+=   
    (const XYZ<Treal, 
     MatrixSymmetric<Treal, Tmatrix>, 
     MatrixGeneral<Treal, Tmatrix> >& smm) {
    assert(this != &smm.C);
    assert(!smm.tB && !smm.tC);
    Tmatrix::symm('L', 'U', smm.A, 
		  *smm.B.matrixPtr, *smm.C.matrixPtr,
		  1, *this->matrixPtr);
    return *this;
  }
    
  /* C = alpha * A * B + beta * C           : A is symmetric */
  template<typename Treal, typename Tmatrix>
    MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator= 
    (const XYZpUV<Treal, 
     MatrixSymmetric<Treal, Tmatrix>, 
     MatrixGeneral<Treal, Tmatrix>, 
     Treal,
     MatrixGeneral<Treal, Tmatrix> >& smmpsm) {
    assert(this != &smmpsm.C);
    assert(!smmpsm.tB && !smmpsm.tC && !smmpsm.tE);
    if (this == &smmpsm.E)
      Tmatrix::symm('L', 'U', smmpsm.A, 
		    *smmpsm.B.matrixPtr, *smmpsm.C.matrixPtr,
		    smmpsm.D, *this->matrixPtr);
    else
      throw Failure("MatrixGeneral<Treal, Tmatrix>::operator="
		    "(const XYZpUV<Treal, MatrixGeneral"
		    "<Treal, Tmatrix>, MatrixSymmetric<Treal, "
		    "Tmatrix>, Treal, MatrixGeneral"
		    "<Treal, Tmatrix> >&) "
		    ":  D = alpha * A * B + beta * C (with A symmetric)"
		    " not supported for C != D");
    return *this;
  }

  /* C = alpha * B * A                      : A is symmetric */
  template<typename Treal, typename Tmatrix>
  MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator= 
  (const XYZ<Treal, 
   MatrixGeneral<Treal, Tmatrix>, 
   MatrixSymmetric<Treal, Tmatrix> >& smm) {
    assert(this != &smm.B);
    assert(!smm.tB && !smm.tC);
    this->matrixPtr.haveDataStructureSet(true);
    Tmatrix::symm('R', 'U', smm.A, 
		  *smm.C.matrixPtr, *smm.B.matrixPtr,
		  0, *this->matrixPtr);
    return *this;
  }

  /* C = B * A                      : A is symmetric */
  template<typename Treal, typename Tmatrix>
  MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator= 
  (const XY<MatrixGeneral<Treal, Tmatrix>, 
   MatrixSymmetric<Treal, Tmatrix> >& mm) {
    assert(this != &mm.A);
    assert(!mm.tA && !mm.tB);
    Tmatrix::symm('R', 'U', 1.0, 
		  *mm.B.matrixPtr, *mm.A.matrixPtr,
		  0, *this->matrixPtr);
    return *this;
  }
  
  /* C += alpha * B * A                      : A is symmetric */
  template<typename Treal, typename Tmatrix>
    MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator+=   
    (const XYZ<Treal, 
     MatrixGeneral<Treal, Tmatrix>, 
     MatrixSymmetric<Treal, Tmatrix> >& smm) {
    assert(this != &smm.B);
    assert(!smm.tB && !smm.tC);
    Tmatrix::symm('R', 'U', smm.A, 
		  *smm.C.matrixPtr, *smm.B.matrixPtr,
		  1, *this->matrixPtr);
    return *this;
  }
    
  /* C = alpha * B * A + beta * C           : A is symmetric */
  template<typename Treal, typename Tmatrix>
    MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator= 
    (const XYZpUV<Treal, 
     MatrixGeneral<Treal, Tmatrix>, 
     MatrixSymmetric<Treal, Tmatrix>, 
     Treal,
     MatrixGeneral<Treal, Tmatrix> >& smmpsm) {
    assert(this != &smmpsm.B);
    assert(!smmpsm.tB && !smmpsm.tC && !smmpsm.tE);
    if (this == &smmpsm.E)
      Tmatrix::symm('R', 'U', smmpsm.A, 
		    *smmpsm.C.matrixPtr, *smmpsm.B.matrixPtr,
		    smmpsm.D, *this->matrixPtr);
    else
      throw Failure("MatrixGeneral<Treal, Tmatrix>::operator="
		    "(const XYZpUV<Treal, MatrixSymmetric"
		    "<Treal, Tmatrix>, MatrixGeneral<Treal, "
		    "Tmatrix>, Treal, MatrixGeneral"
		    "<Treal, Tmatrix> >&) "
		    ":  D = alpha * B * A + beta * C (with A symmetric)"
		    " not supported for C != D");
    return *this;
  }


  /** C = alpha * A * B           : A and B are symmetric */
  template<typename Treal, typename Tmatrix>
    MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator= 
    (const XYZ<Treal, 
     MatrixSymmetric<Treal, Tmatrix>, 
     MatrixSymmetric<Treal, Tmatrix> >& smm) {
    assert(!smm.tB && !smm.tC);
    this->matrixPtr.haveDataStructureSet(true);
    Tmatrix::ssmm(smm.A, 
		  *smm.B.matrixPtr, 
		  *smm.C.matrixPtr,
		  0, 
		  *this->matrixPtr);
    return *this;    
  }

  /** C = A * B           : A and B are symmetric */
  template<typename Treal, typename Tmatrix>
    MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator= 
    (const XY<MatrixSymmetric<Treal, Tmatrix>, 
     MatrixSymmetric<Treal, Tmatrix> >& mm) {
    assert(!mm.tB);
    Tmatrix::ssmm(1.0, 
		  *mm.A.matrixPtr, 
		  *mm.B.matrixPtr,
		  0, 
		  *this->matrixPtr);
    return *this;    
  }

  /** C += alpha * A * B           : A and B are symmetric */
  template<typename Treal, typename Tmatrix>
    MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator+= 
    (const XYZ<Treal, 
     MatrixSymmetric<Treal, Tmatrix>, 
     MatrixSymmetric<Treal, Tmatrix> >& smm) {
    assert(!smm.tB && !smm.tC);
    Tmatrix::ssmm(smm.A, 
		  *smm.B.matrixPtr, 
		  *smm.C.matrixPtr,
		  1, 
		  *this->matrixPtr);
    return *this;    
  }


  /** C = alpha * A * B + beta * C           : A and B are symmetric */
  template<typename Treal, typename Tmatrix>
    MatrixGeneral<Treal, Tmatrix>& 
    MatrixGeneral<Treal, Tmatrix>::operator= 
    (const XYZpUV<Treal, 
     MatrixSymmetric<Treal, Tmatrix>, 
     MatrixSymmetric<Treal, Tmatrix>, 
     Treal,
     MatrixGeneral<Treal, Tmatrix> >& smmpsm) {
    assert(!smmpsm.tB && !smmpsm.tC && !smmpsm.tE);
    if (this == &smmpsm.E)
      Tmatrix::ssmm(smmpsm.A, 
		    *smmpsm.B.matrixPtr, 
		    *smmpsm.C.matrixPtr,
		    smmpsm.D, 
		    *this->matrixPtr);
    else
      throw Failure("MatrixGeneral<Treal, Tmatrix>::"
		    "operator=(const XYZpUV<"
		    "Treal, MatrixSymmetric<Treal, Tmatrix>, "
		    "MatrixSymmetric<Treal, Tmatrix>, Treal, "
		    "MatrixGeneral<Treal, Tmatrix> >&) "
		    ":  D = alpha * A * B + beta * C (with A and B symmetric)"
		    " not supported for C != D");
    return *this;    
  }
   

    
    /* OPERATIONS INVOLVING UPPER TRIANGULAR MATRICES */
    
    /* B = alpha * op(A) * B                  : A is upper triangular */
    template<typename Treal, typename Tmatrix>
      MatrixGeneral<Treal, Tmatrix>& 
      MatrixGeneral<Treal, Tmatrix>::operator= 
      (const XYZ<Treal, 
       MatrixTriangular<Treal, Tmatrix>, 
       MatrixGeneral<Treal, Tmatrix> >& smm) {
      assert(!smm.tC);
      if (this == &smm.C)
	Tmatrix::trmm('L', 'U', smm.tB, smm.A, 
		      *smm.B.matrixPtr, *this->matrixPtr);
      else
	throw Failure("MatrixGeneral<Treal, Tmatrix>::operator="
		      "(const XYZ<Treal, MatrixTriangular"
		      "<Treal, Tmatrix>, MatrixGeneral<Treal,"
		      " Tmatrix> >& : D = alpha * op(A) * B (with"
		      " A upper triangular) not supported for B != D");
      return *this;
    }


    /* A = alpha * A * op(B)                  : B is upper triangular */
    template<typename Treal, typename Tmatrix>
      MatrixGeneral<Treal, Tmatrix>& 
      MatrixGeneral<Treal, Tmatrix>::operator= 
      (const XYZ<Treal, 
       MatrixGeneral<Treal, Tmatrix>, 
       MatrixTriangular<Treal, Tmatrix> >& smm) {
      assert(!smm.tB);
      if (this == &smm.B)
	Tmatrix::trmm('R', 'U', smm.tC, smm.A, 
		      *smm.C.matrixPtr, *this->matrixPtr);
      else
	throw Failure("MatrixGeneral<Treal, Tmatrix>::operator="
		      "(const XYZ<Treal, MatrixGeneral"
		      "<Treal, Tmatrix>, MatrixTriangular<Treal,"
		      " Tmatrix> >& : D = alpha * A * op(B) (with"
		      " B upper triangular) not supported for A != D");
      return *this;
    }
    

    /******* FUNCTIONS DECLARED OUTSIDE CLASS */
    template<typename Treal, typename Tmatrix>
      Treal trace(const XYZ<Treal, 
		  MatrixGeneral<Treal, Tmatrix>, 
		  MatrixGeneral<Treal, Tmatrix> >& smm) {
      if ((!smm.tB && !smm.tC) || (smm.tB && smm.tC))
	return smm.A * MatrixGeneral<Treal, Tmatrix>::
	  trace_ab(smm.B, smm.C);
      else if (smm.tB) 
	return smm.A * MatrixGeneral<Treal, Tmatrix>::
	  trace_aTb(smm.B, smm.C);
      else
	return smm.A * MatrixGeneral<Treal, Tmatrix>::
	  trace_aTb(smm.C, smm.B);
    }


} /* end namespace mat */
#endif


