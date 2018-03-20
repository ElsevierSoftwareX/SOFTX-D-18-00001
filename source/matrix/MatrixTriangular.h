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

/** @file MatrixTriangular.h Triangular matrix class
 *
 * Copyright(c) Emanuel Rubensson 2006
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date January 2006
 *
 */
#ifndef MAT_MATRIXTRIANGULAR
#define MAT_MATRIXTRIANGULAR
#include <stdexcept>
#include "MatrixBase.h"
namespace mat {
  /** Upper non-unit triangular matrix 
   *
   *
   * This class belongs to the matrix API
   *
   * Treal: Type for real numbers
   *
   * Tmatrix: The matrix class
   *
   * @see MatrixBase
   * @see MatrixGeneral
   * @see MatrixSymmetric
   * 
   * 
   */  
  template<typename Treal, typename Tmatrix>
    class MatrixTriangular : public MatrixBase<Treal, Tmatrix> {
  public:
    typedef VectorGeneral<Treal, typename Tmatrix::VectorType> VectorType;
    typedef Treal real;

    MatrixTriangular()
      :MatrixBase<Treal, Tmatrix>() {} /**< Default constructor  */
    explicit MatrixTriangular(const MatrixTriangular<Treal, Tmatrix>& tri)
      :MatrixBase<Treal, Tmatrix>(tri) {} /**< Copy constructor  */
      
    MatrixTriangular<Treal, Tmatrix>& 
      operator=(const MatrixTriangular<Treal, Tmatrix>& tri) {
      MatrixBase<Treal, Tmatrix>::operator=(tri);
      return *this;
    } 
    
    inline MatrixTriangular<Treal, Tmatrix>& operator=(int const k) {
      *this->matrixPtr = k;
      return *this;
    } /**< Set matrix to zero or identity: A = 0 or A = 1 
       *
       * Only zero and one are valid arguments.
       *
       * 
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
     * The elements to be added must be in upper triangle.
     * Information about sizes and blocks for rows as well as columns 
     * must also be given.
     * @warning All indexing start at zero.
     */

    /** Add given set of values to the matrix (+=). 
     *  The values should be in upper triangle.
     */
    inline void add_values
      (std::vector<int> const & rowind, 
       std::vector<int> const & colind, 
       std::vector<Treal> const & values) {
      this->matrixPtr->syAddValues(rowind, colind, values);
    }


    inline void get_values
      (std::vector<int> const & rowind, 
       std::vector<int> const & colind, 
       std::vector<Treal> & values) const {
      this->matrixPtr->syGetValues(rowind, colind, values);
    }
    /**< Get values given by row and column index lists.
     * Input arrays contain row and column indices.
     * The wanted elements must in upper triangle.
     * The output array contains values for the given indices.
     * @warning All indexing start at zero.
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

#if 0
    inline void fullmatrix(Treal* const full, int const size) const {
      this->matrixPtr->fullmatrix(full, size, size);
    } /* FIXME? Should triangular matrix always have zeros below diagonal? */
#endif

    inline void inch(const MatrixGeneral<Treal, Tmatrix>& SPD, 
		     const Treal threshold,
		     const side looking = left,
		     const inchversion version = unstable) {
      Tmatrix::inch(*SPD.matrixPtr, *this->matrixPtr, 
		    threshold, looking, version);
    }
    inline void inch(const MatrixSymmetric<Treal, Tmatrix>& SPD, 
		     const Treal threshold,
		     const side looking = left,
		     const inchversion version = unstable) {
      this->matrixPtr.haveDataStructureSet(true);
      Tmatrix::syInch(*SPD.matrixPtr, *this->matrixPtr, 
		      threshold, looking, version);
    }
    
    void thresh(Treal const threshold,
		normType const norm);

    inline Treal frob() const {
      return this->matrixPtr->frob();
    }
    Treal eucl(Treal const requestedAccuracy,
	       int maxIter = -1) const;

    Treal eucl_thresh(Treal const threshold);
    Treal eucl_thresh_congr_trans_measure(Treal const threshold,
					  MatrixSymmetric<Treal, Tmatrix> & trA);
    inline void frob_thresh(Treal threshold) {
      this->matrixPtr->frob_thresh(threshold);
    }    
    inline size_t nnz() const { /* Note: size_t instead of int here to avoid integer overflow. */
      return this->matrixPtr->nnz();
    }
    inline size_t nvalues() const { /* Note: size_t instead of int here to avoid integer overflow. */
      return this->matrixPtr->nvalues();
    }


    inline void write_to_buffer(void* buffer, const int n_bytes) const {
      this->write_to_buffer_base(buffer, n_bytes, matrix_triang);
    }
    inline void read_from_buffer(void* buffer, const int n_bytes) {
      this->read_from_buffer_base(buffer, n_bytes, matrix_triang);
    }

    void random() {
      this->matrixPtr->syRandom();
    }

    /** Uses rule depending on the row and column indexes to set matrix elements
     * The Trule class should have the function "Treal = set(int row,int col)"
     * which is used to set the elements.
     */
    template<typename TRule>
      void setElementsByRule(TRule & rule) {
      this->matrixPtr->trSetElementsByRule(rule);
      return;
    }

    /** B += alpha * A */
    MatrixTriangular<Treal, Tmatrix>& operator+= 
      (XY<Treal, MatrixTriangular<Treal, Tmatrix> > const & sm);
 

    std::string obj_type_id() const {return "MatrixTriangular";}
    protected:
    inline void writeToFileProt(std::ofstream & file) const {
      this->writeToFileBase(file, matrix_triang);
    }
    inline void readFromFileProt(std::ifstream & file) {
      this->readFromFileBase(file, matrix_triang);
    }

    private:

  };

  /* B += alpha * A */
  template<typename Treal, typename Tmatrix>
    inline MatrixTriangular<Treal, Tmatrix>& 
    MatrixTriangular<Treal, Tmatrix>::operator+= 
    (XY<Treal, MatrixTriangular<Treal, Tmatrix> > const & sm) {
    assert(!sm.tB);
    Tmatrix::add(sm.A, *sm.B.matrixPtr, *this->matrixPtr);
    return *this;
  }


  template<typename Treal, typename Tmatrix>
    void MatrixTriangular<Treal, Tmatrix>::
    thresh(Treal const threshold,
	   normType const norm) {
    switch (norm) {
    case frobNorm:
      this->frob_thresh(threshold);
      break;
    default:
      throw Failure("MatrixTriangular<Treal, Tmatrix>::"
		    "thresh(Treal const, "
		    "normType const): "
		    "Thresholding not imlpemented for selected norm");
    }
  }

  template<typename Treal, typename Tmatrix>
    Treal MatrixTriangular<Treal, Tmatrix>::
    eucl(Treal const requestedAccuracy,
	 int maxIter) const {
    VectorType guess;
    SizesAndBlocks cols;
    this->getCols(cols);
    guess.resetSizesAndBlocks(cols);
    guess.rand();
    mat::ATAMatrix<MatrixTriangular<Treal, Tmatrix>, Treal> ztz(*this);    
    if (maxIter < 0)
      maxIter = this->get_nrows() * 100;
    arn::LanczosLargestMagnitudeEig
      <Treal, ATAMatrix<MatrixTriangular<Treal, Tmatrix>, Treal>, VectorType>
      lan(ztz, guess, maxIter);
    lan.setRelTol( 100 * mat::getMachineEpsilon<Treal>() );
    lan.run();
    Treal eVal = 0;
    Treal acc = 0;
    lan.getLargestMagnitudeEig(eVal, acc);
    Interval<Treal> euclInt( template_blas_sqrt(eVal-acc),
			     template_blas_sqrt(eVal+acc) );    
    if ( euclInt.low() < 0 )
      euclInt = Interval<Treal>( 0, template_blas_sqrt(eVal+acc) );
    if ( euclInt.length() / 2.0 > requestedAccuracy ) {
      std::cout << "req: " << (double)requestedAccuracy
		<< "  obt: " << (double)(euclInt.length() / 2.0) << std::endl;
      throw std::runtime_error("Desired accuracy not obtained in Lanczos.");
    }
    return euclInt.midPoint();
  }

#if 1

  template<typename Treal, typename Tmatrix>
    Treal MatrixTriangular<Treal, Tmatrix>::
    eucl_thresh(Treal const threshold) {
    EuclTruncationGeneral<MatrixTriangular<Treal, Tmatrix>, Treal> TruncObj( *this );
    return TruncObj.run( threshold );    
  }
  
#endif  

  template<typename Treal, typename Tmatrix>
    Treal MatrixTriangular<Treal, Tmatrix>::
    eucl_thresh_congr_trans_measure(Treal const threshold,
			       MatrixSymmetric<Treal, Tmatrix> & trA) {
    EuclTruncationCongrTransMeasure<MatrixTriangular<Treal, Tmatrix>, MatrixSymmetric<Treal, Tmatrix>, Treal> TruncObj(*this, trA);
    return TruncObj.run( threshold );    
  }
  

} /* end namespace mat */
#endif


