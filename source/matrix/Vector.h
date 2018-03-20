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

/** @file Vector.h 
 *
 * Copyright(c) Emanuel Rubensson 2006
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date October 2006
 *
 */

#ifndef MAT_VECTOR
#define MAT_VECTOR
#include "VectorHierarchicBase.h"
namespace mat{
  template<class Treal, class Telement>
    class Matrix;
  /** Vector class
   *
   * This class is used to obtain the hierarchic vector data structure.
   *  
   * @see VectorHierarchicBase
   * @see Permutation
   *
   */
  template<class Treal, class Telement = Treal>
    class Vector: public VectorHierarchicBase<Treal, Telement> {
    public:
    typedef Telement ElementType;
    //    template<typename TmatrixElement>
    //    friend class Matrix<Treal, TmatrixElement>;
    Vector():VectorHierarchicBase<Treal, Telement>(){}

    void allocate() {
      assert(!this->is_empty());
      assert(this->is_zero());
      this->elements = allocateElements<Telement>(this->n());
      SizesAndBlocks rowSAB;
      for (int row = 0; row < this->rows.getNBlocks(); row++) {
	rowSAB = this->rows.getSizesAndBlocksForLowerLevel(row);
	(*this)(row).resetRows(rowSAB);
      }
    }
    
    void assignFromFull(std::vector<Treal> const & fullVector);
    
    void addFromFull(std::vector<Treal> const & fullVector);

    void fullVector(std::vector<Treal> & fullVector) const;

    Vector<Treal, Telement>& 
    operator=(const Vector<Treal, Telement>& vec) {
      VectorHierarchicBase<Treal, Telement>::operator=(vec);
      return *this;
    } 

    void clear();

    void writeToFile(std::ofstream & file) const;
    void readFromFile(std::ifstream & file);
    Vector<Treal, Telement>& operator=(int const k);

    inline void randomNormalized() {
      this->random();
      (*this) *= (1.0 / this->eucl());
    }
    void random();

    inline Treal eucl() const {
      return template_blas_sqrt(dot(*this,*this));
    }

    /* LEVEL 1 operations */
    Vector<Treal, Telement>& operator*=(const Treal alpha);
    static Treal dot(Vector<Treal, Telement> const & x, 
		     Vector<Treal, Telement> const & y);

    /* y += alpha * x */
    static void axpy(Treal const & alpha, 
		     Vector<Treal, Telement> const & x, 
		     Vector<Treal, Telement> & y);


    /* LEVEL 2 operations */
    /** gemv:
     * y = alpha * A * x + beta * y,   or   
     * y = alpha * transpose(A) * x + beta * y 
     */
    template<typename TmatrixElement>
    static void gemv(bool const tA, Treal const alpha, 
		     Matrix<Treal, TmatrixElement> const & A,
		     Vector<Treal, Telement> const & x,
		     Treal const beta, 
		     Vector<Treal, Telement>& y);

    /** symv:
     * y = alpha * A * x + beta * y, where A is symmetric
     */
    template<typename TmatrixElement>
    static void symv(char const uplo, Treal const alpha, 
		     Matrix<Treal, TmatrixElement> const & A,
		     Vector<Treal, Telement> const & x,
		     Treal const beta, 
		     Vector<Treal, Telement>& y);
    /** trmv:
     * x = A * x,    or 
     * x = transpose(A) * x,  where A is triangular
     */
    template<typename TmatrixElement>
    static void trmv(char const uplo, const bool tA,  
		     Matrix<Treal, TmatrixElement> const & A,
		     Vector<Treal, Telement> & x);


#if 0    /* OLD ROUTINES */
    void assign_from_full(Treal const * const fullvector, const int totn);
    /* Convert to full vector */
    void fullvector(Treal * const full, const int totn) const; 
    
    
    

    




#endif /* END OLD ROUTINES */
  }; /* end class Vector */


  template<class Treal, class Telement>
    void Vector<Treal, Telement>::
    assignFromFull(std::vector<Treal> const & fullVector) {
    addFromFull(fullVector);
  }

  template<class Treal, class Telement>
    void Vector<Treal, Telement>::
    addFromFull(std::vector<Treal> const & fullVector) {
    if (this->is_zero())
      allocate();
    for (int ind = 0; ind < this->n(); ind++)
      (*this)(ind).addFromFull(fullVector);
  }

  template<class Treal, class Telement> 
    void Vector<Treal, Telement>::
    fullVector(std::vector<Treal> & fullVec) const {
    if (this->is_zero()) {
      fullVec.resize(this->rows.getNTotalScalars());
      for (int row = 0; row < this->nScalars(); ++row )
	fullVec[this->rows.getOffset()+row] = 0;
    }
    else
      for (int ind = 0; ind < this->n(); ind++)
	(*this)(ind).fullVector(fullVec);
  }
  

  template<class Treal, class Telement>
    void Vector<Treal, Telement>::clear() {
    freeElements(this->elements);
    this->elements = 0;
  }

  template<class Treal, class Telement>
    void Vector<Treal, Telement>:: 
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
      for (int i = 0; i < this->n(); i++)
	this->elements[i].writeToFile(file);
    }
  }
  template<class Treal, class Telement>
    void Vector<Treal, Telement>:: 
    readFromFile(std::ifstream & file) {
    int const ZERO = 0;
    int const ONE  = 1;
    char tmp[sizeof(int)];
    file.read(tmp, (std::ifstream::pos_type)sizeof(int));
    switch ((int)*tmp) {
    case ZERO:
      (*this) = 0;
      break;
    case ONE:
      if (this->is_zero()) 
	allocate();
      for (int i = 0; i < this->n(); i++)
	this->elements[i].readFromFile(file);
      break;
    default:
      throw Failure("Vector<Treal, Telement>::" 
		    "readFromFile(std::ifstream & file):"
		    "File corruption int value not 0 or 1");
    }
  }

  template<class Treal, class Telement>
    Vector<Treal, Telement>& Vector<Treal, Telement>:: 
    operator=(int const k) {
    if (k == 0) 
      this->clear();
    else
      throw Failure("Vector::operator=(int k) only "
		    "implemented for k = 0");
    return *this;
  }

  template<class Treal, class Telement>
    void Vector<Treal, Telement>::random() {
    if (this->is_zero()) 
      allocate();
    for (int ind = 0; ind < this->n(); ind++)
      (*this)(ind).random();    
  }

  /* LEVEL 1 operations */
  
  template<class Treal, class Telement>
    Vector<Treal, Telement>& Vector<Treal, Telement>:: 
    operator*=(const Treal alpha) {
    if (!this->is_zero() && alpha != 1) {
      for (int ind = 0; ind < this->n(); ind++)
	(*this)(ind) *= alpha;
    }
    return *this;
  }

  template<class Treal, class Telement>
    Treal Vector<Treal, Telement>:: 
    dot(Vector<Treal, Telement> const & x, 
	Vector<Treal, Telement> const & y) {
    assert(x.n() == y.n());
    if (x.is_zero() || y.is_zero())
      return 0;
    Treal dotProduct = 0;
    for (int ind = 0; ind < x.n(); ind++)
      dotProduct += Telement::dot(x(ind), y(ind));    
    return dotProduct;
  }

  /* y += alpha * x */
  template<class Treal, class Telement>
    void Vector<Treal, Telement>:: 
    axpy(Treal const & alpha, 
	 Vector<Treal, Telement> const & x, 
	 Vector<Treal, Telement> & y) {
    assert(x.n() == y.n());
    if (x.is_zero())
      return;
    if (y.is_zero()) {
      y.allocate();
    }
    for (int ind = 0; ind < x.n(); ind++)
      Telement::axpy(alpha, x(ind), y(ind));
  }

  /* LEVEL 2 operations */
  
  /** gemv:
   * y = alpha * A * x + beta * y,   or   
   * y = alpha * transpose(A) * x + beta * y 
   */
  template<class Treal, class Telement>
    template<typename TmatrixElement>
    void Vector<Treal, Telement>::
    gemv(bool const tA, Treal const alpha, 
	 Matrix<Treal, TmatrixElement> const & A,
	 Vector<Treal, Telement> const & x,
	 Treal const beta, 
	 Vector<Treal, Telement>& y) {
    if (y.is_empty()) {
      assert(beta == 0);
      y.resetRows(x.rows);
    }
    if ((A.is_zero() || x.is_zero() || alpha == 0) && 
	(y.is_zero() || beta == 0))
      y = 0;
    else {
      Treal beta_tmp = beta;
      if (y.is_zero()) {
	y.allocate();
	beta_tmp = 0;
      }
      if (A.is_zero() || x.is_zero() || alpha == 0)
	y *= beta_tmp;
      else {
	MAT_OMP_INIT;
	if (!tA) {
	  if (A.ncols() != x.n() || A.nrows() != y.n())
	    throw Failure("Vector<Treal, Telement>::"
			  "gemv(bool const, Treal const, "
			  "const Matrix<Treal, Telement>&, "
			  "const Vector<Treal, Telement>&, "
			  "Treal const, const Vector<Treal, "
			  "Telement>&): "
			  "Incorrect dimensions for matrix-vector product");
	  else {
	    int A_nrows = A.nrows();
#ifdef _OPENMP
#pragma omp parallel for if(A.level() == Params::getMatrixParallelLevel()) num_threads(Params::getNProcs()) default(shared) schedule(dynamic) 
#endif
	    for (int row = 0; row < A_nrows; row++) {
	      MAT_OMP_START;
	      Telement::gemv(tA, alpha, A(row, 0), x(0), beta_tmp, y(row));
	      for (int col = 1; col < A.ncols(); col++) 
		Telement::gemv(tA, alpha, A(row, col), x(col), 1.0, y(row));
	      MAT_OMP_END;
	    }
	  } /* end else */
	} /* end if (!tA) */
	else {
	  assert(tA);
	  if (A.nrows() != x.n() || A.ncols() != y.n())
	    throw Failure("Vector<Treal, Telement>::"
			  "gemv(bool const, Treal const, "
			  "const Matrix<Treal, Telement>&, "
			  "const Vector<Treal, Telement>&, "
			  "Treal const, const Vector<Treal, "
			  "Telement>&): "
			  "Incorrect dimensions for matrix-vector product");
	  else {
	    int A_ncols = A.ncols();
#ifdef _OPENMP
#pragma omp parallel for if(A.level() == Params::getMatrixParallelLevel()) num_threads(Params::getNProcs()) default(shared) schedule(dynamic) 
#endif
	    for (int col = 0; col < A_ncols; col++) {
	      MAT_OMP_START;
	      Telement::gemv(tA, alpha, A(0, col), x(0), beta_tmp, y(col));
	      for (int row = 1; row < A.nrows(); row++)
		Telement::gemv(tA, alpha, A(row, col), x(row), 1.0, y(col));
	      MAT_OMP_END;
	    }
	  } /* end else */
	} /* end else */
	MAT_OMP_FINALIZE;
      } /* end else */
    } /* end else */
  }

  /** symv:
   * y = alpha * A * x + beta * y, where A is symmetric
   */
  template<class Treal, class Telement>
    template<typename TmatrixElement>
    void Vector<Treal, Telement>::
    symv(char const uplo, Treal const alpha, 
	 Matrix<Treal, TmatrixElement> const & A,
	 Vector<Treal, Telement> const & x,
	 Treal const beta, 
	 Vector<Treal, Telement>& y) {
    if (y.is_empty()) {
      assert(beta == 0);
      y.resetRows(x.rows);
    }
    if (x.n() != y.n() || A.nrows() != A.ncols() || A.ncols() != x.n())
      throw Failure("Vector<Treal, Telement>::"
		    "symv(char const uplo, Treal const, " 
		    "const Matrix<Treal, Telement>&, "
		    "const Vector<Treal, Telement>&, "
		    "Treal const, const Vector<Treal, Telement>&):"
		    "Incorrect dimensions for symmetric "
		    "matrix-vector product");
    if (uplo != 'U') 
      throw Failure("Vector<class Treal, class Telement>::"
		    "symv only implemented for symmetric matrices in "
		    "upper triangular storage");
    if ((A.is_zero() || x.is_zero() || alpha == 0) && 
	(y.is_zero() || beta == 0))
      y = 0;
    else {
      Treal beta_tmp = beta;
      if (y.is_zero()) {
	y.allocate();
	beta_tmp = 0;
      }
      if (A.is_zero() || x.is_zero() || alpha == 0)
	y *= beta_tmp;
      else {
	MAT_OMP_INIT;
#ifdef _OPENMP
#pragma omp parallel if(A.level() == Params::getMatrixParallelLevel()) num_threads(Params::getNProcs()) default(shared)
#endif
	{
	  /* Diagonal */
	  int A_ncols = A.ncols();
#ifdef _OPENMP
#pragma omp for  schedule(dynamic)
#endif
	  for (int rc = 0; rc < A_ncols; rc++) {
	    MAT_OMP_START;
	    Telement::symv(uplo, alpha, A(rc,rc), x(rc), beta_tmp, y(rc));
	    MAT_OMP_END;
	  }
	  /* Upper triangle */
	  int A_nrows = A.nrows();
#ifdef _OPENMP
#pragma omp for  schedule(dynamic)
#endif
	  for (int row = 0; row < A_nrows - 1; row++) {
	    MAT_OMP_START;
	    for (int col = row + 1; col < A.ncols(); col++)
	      Telement::gemv(false, alpha, A(row, col), x(col), 1.0, y(row));
	    MAT_OMP_END;
	  }
	  /* Lower triangle */
#ifdef _OPENMP
#pragma omp for  schedule(dynamic)
#endif
	  for (int row = 1; row < A_nrows; row++) {
	    MAT_OMP_START;
	    for (int col = 0; col < row; col++)
	      Telement::gemv(true, alpha, A(col, row), x(col), 1.0, y(row));
	    MAT_OMP_END;
	  }
	} /* end omp parallel*/
	MAT_OMP_FINALIZE;
      } /* end else */
    } /* end else */
  }
  
  template<class Treal, class Telement>
    template<typename TmatrixElement>
    void Vector<Treal, Telement>::
    trmv(char const uplo, const bool tA,  
	 Matrix<Treal, TmatrixElement> const & A,
	 Vector<Treal, Telement> & x) {
    if (A.nrows() != A.ncols() || A.ncols() != x.n())
      throw Failure("Vector<Treal, Telement>::"
		    "trmv(...):"
		    "Incorrect dimensions for triangular "
		    "matrix-vector product");
    if (uplo != 'U') 
      throw Failure("Vector<class Treal, class Telement>::"
		    "trmv only implemented for upper triangular matrices");
    if ( ( A.is_zero() || x.is_zero() ) ) {
      x = 0;
      return;
    }
    if (!tA) {
      // not transposed
      for (int row = 0; row < A.nrows(); row++) {
	Telement::trmv(uplo, tA, A(row,row), x(row));
	for (int col = row + 1; col < A.ncols(); col++)
	  Telement::gemv(tA, (Treal)1.0, A(row, col), x(col), 1.0, x(row));
      }
      return;
    }
    // transposed
    for (int col = A.ncols() - 1; col >= 0; col--) {
      Telement::trmv(uplo, tA, A(col,col), x(col));
      for (int row = 0; row < col; row++) 
	Telement::gemv(tA, (Treal)1.0, A(row, col), x(row), 1.0, x(col));
    }
  }
  

  
  
  /***************************************************************************/
  /***************************************************************************/
  /*           Specialization for Telement = Treal                           */
  /***************************************************************************/
  /***************************************************************************/
  
  template<class Treal>
    class Vector<Treal>: public VectorHierarchicBase<Treal> {
  public:
    friend class Matrix<Treal>;
    Vector()
      :VectorHierarchicBase<Treal>(){}
    
    void allocate() {
      assert(!this->is_empty());
      assert(this->is_zero());
      this->elements = allocateElements<Treal>(this->n());
      for (int ind = 0; ind < this->n(); ind++) 
	this->elements[ind] = 0;
    }
    
    void assignFromFull(std::vector<Treal> const & fullVector);
    
    void addFromFull(std::vector<Treal> const & fullVector);
    
    void fullVector(std::vector<Treal> & fullVector) const;
    
    
    Vector<Treal>& 
      operator=(const Vector<Treal>& vec) {
      VectorHierarchicBase<Treal>::operator=(vec);
      return *this;
    } 

    void clear(); /**< Set vector to zero and delete all arrays */
    
    void writeToFile(std::ofstream & file) const;
    void readFromFile(std::ifstream & file);
    
    Vector<Treal>& operator=(int const k);
    

    inline void randomNormalized() {
      this->random();
      (*this) *= 1 / this->eucl();
    }
    void random();

    inline Treal eucl() const {
      return template_blas_sqrt(dot(*this,*this));
    }

    /* LEVEL 1 operations */
    Vector<Treal>& operator*=(const Treal alpha);

    static Treal dot(Vector<Treal> const & x, 
		     Vector<Treal> const & y);

    
    /* y += alpha * x */
    static void axpy(Treal const & alpha, 
		     Vector<Treal> const & x, 
		     Vector<Treal> & y);

    /* LEVEL 2 operations */
    /** gemv:
     * y = alpha * A * x + beta * y,   or   
     * y = alpha * transpose(A) * x + beta * y 
     */
     static void gemv(bool const tA, Treal const alpha, 
		     Matrix<Treal> const & A,
		     Vector<Treal> const & x,
		     Treal const beta, 
		     Vector<Treal>& y);

    /** symv:
     * y = alpha * A * x + beta * y, where A is symmetric
     */
    static void symv(char const uplo, Treal const alpha, 
		     Matrix<Treal> const & A,
		     Vector<Treal> const & x,
		     Treal const beta, 
		     Vector<Treal>& y);

    /** trmv:
     * x = A * x,    or 
     * x = transpose(A) * x,  where A is triangular
     */
    static void trmv(char const uplo, const bool tA,  
		     Matrix<Treal> const & A,
		     Vector<Treal> & x);
      
  }; /* end class Vector specialization */
  
  
  template<class Treal>
    void Vector<Treal>::
    assignFromFull(std::vector<Treal> const & fullVector) {
    addFromFull(fullVector);
  }
  
  template<class Treal>
    void Vector<Treal>::
    addFromFull(std::vector<Treal> const & fullVector) {
    if (this->is_zero())
      allocate();
    assert((unsigned)this->rows.getNTotalScalars() == fullVector.size()); 
    /*  Assertion AFTER empty check done 
     *  by allocate()
     */
    for (int row = 0; row < this->n(); ++row )
      (*this)(row) += fullVector[this->rows.getOffset()+row];
  }

  template<class Treal> 
    void Vector<Treal>::
    fullVector(std::vector<Treal> & fullVec) const {
    fullVec.resize(this->rows.getNTotalScalars());
    if (this->is_zero()) 
      for (int row = 0; row < this->nScalars(); ++row )
	fullVec[this->rows.getOffset()+row] = 0;
    else
      for (int row = 0; row < this->n(); ++row )
	fullVec[this->rows.getOffset()+row] = (*this)(row);
  }
  

  template<class Treal>
    void Vector<Treal>::clear() {
    freeElements(this->elements);
    this->elements = 0;
  }


  template<class Treal>
    void Vector<Treal>:: 
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
      file.write(tmpel,sizeof(Treal) * this->n());
    }
  }

  template<class Treal>
    void Vector<Treal>:: 
    readFromFile(std::ifstream & file) {
    int const ZERO = 0;
    int const ONE  = 1;
    char tmp[sizeof(int)];
    file.read(tmp, (std::ifstream::pos_type)sizeof(int));
    switch ((int)*tmp) {
    case ZERO:
      (*this) = 0;
      break;
    case ONE:
      if (this->is_zero())
	allocate();
      file.read((char*)this->elements, sizeof(Treal) * this->n());
      break;
    default:
      throw Failure("Vector<Treal>::" 
		    "readFromFile(std::ifstream & file):"
		    "File corruption, int value not 0 or 1");
    }
  }

  template<class Treal>
    Vector<Treal>& Vector<Treal>:: 
    operator=(int const k) {
    if (k == 0) 
      this->clear();
    else
      throw Failure("Vector::operator=(int k) only implemented for k = 0");
    return *this;
  }

  template<class Treal>
    void Vector<Treal>::random() {
    if (this->is_zero())
      allocate();
    for (int ind = 0; ind < this->n(); ind++)
      (*this)(ind) = rand() / (Treal)RAND_MAX;    
  }

  /* LEVEL 1 operations */
  template<class Treal>
    Vector<Treal>& Vector<Treal>:: 
    operator*=(const Treal alpha) {
    if (!this->is_zero() && alpha != 1) {
      int const ONE = 1;
      mat::scal(&this->n(),&alpha,this->elements,&ONE);
    }
    return *this;
  }

  template<class Treal>
    Treal Vector<Treal>:: 
    dot(Vector<Treal> const & x, 
	Vector<Treal> const & y) {
    assert(x.n() == y.n());
    if (x.is_zero() || y.is_zero())
      return 0;
    else {
      int const ONE = 1;
      return mat::dot(&x.n(), x.elements, &ONE, y.elements, &ONE);
    }
  }

  /* y += alpha * x */
  template<class Treal>
    void Vector<Treal>:: 
    axpy(Treal const & alpha, 
	 Vector<Treal> const & x, 
	 Vector<Treal> & y) {
    assert(x.n() == y.n());
    if (x.is_zero())
      return;
    if (y.is_zero()) {
      y.allocate();
      for (int ind = 0; ind < y.n(); ind++)
	y.elements[ind] = 0; /* fill with zeros */
    }
    int const ONE = 1;
    mat::axpy(&x.n(), &alpha, x.elements, &ONE, y.elements, &ONE);
  }


  /* LEVEL 2 operations */
  /** gemv:
   * y = alpha * A * x + beta * y,   or   
   * y = alpha * transpose(A) * x + beta * y 
   */
  template<class Treal>
    void Vector<Treal>:: 
    gemv(bool const tA, Treal const alpha, 
	 Matrix<Treal> const & A,
	 Vector<Treal> const & x,
	 Treal const beta, 
	 Vector<Treal>& y) {
    if (y.is_empty()) {
      assert(beta == 0);
      y.resetRows(x.rows);
    }
    if ((A.is_zero() || x.is_zero() || alpha == 0) && 
	(y.is_zero() || beta == 0))
      y = 0;
    else {
      Treal beta_tmp = beta;
      if (y.is_zero()) {
	y.allocate();
	beta_tmp = 0;
      }
      if (A.is_zero() || x.is_zero() || alpha == 0)
	y *= beta_tmp;
      else {
	int const ONE = 1;
	if (!tA) {
	  if (A.ncols() != x.n() || A.nrows() != y.n())
	    throw Failure("Vector<Treal, Telement>::"
			  "gemv(bool const, Treal const, "
			  "const Matrix<Treal, Telement>&, "
			  "const Vector<Treal, Telement>&, "
			  "Treal const, const Vector<Treal, "
			  "Telement>&): "
			  "Incorrect dimensions for matrix-vector product");
	  else {
	    mat::gemv("N", &A.nrows(), &A.ncols(), &alpha, A.elements,
		      &A.nrows(),x.elements,&ONE,&beta_tmp,y.elements,&ONE);
	  } /* end else */
	} /* end if (!tA) */
	else {
	  assert(tA);
	  if (A.nrows() != x.n()  || A.ncols() != y.n())
	    throw Failure("Vector<Treal, Telement>::"
			  "gemv(bool const, Treal const, "
			  "const Matrix<Treal, Telement>&, "
			  "const Vector<Treal, Telement>&, "
			  "Treal const, const Vector<Treal, "
			  "Telement>&): "
			  "Incorrect dimensions for matrix-vector product");
	  else {
	    mat::gemv("T", &A.nrows(), &A.ncols(), &alpha, A.elements,
		      &A.nrows(),x.elements,&ONE,&beta_tmp,y.elements,&ONE);
	  } /* end else */
	} /* end else */
      } /* end else */
    } /* end else */
  }

  /** symv:
   * y = alpha * A * x + beta * y, where A is symmetric
   */
  template<class Treal>
    void Vector<Treal>:: 
    symv(char const uplo, Treal const alpha, 
	 Matrix<Treal> const & A,
	 Vector<Treal> const & x,
	 Treal const beta, 
	 Vector<Treal>& y) {
    if (y.is_empty()) {
      assert(beta == 0);
      y.resetRows(x.rows);
    }    
    if (x.n() != y.n() || A.nrows() != A.ncols() || A.ncols() != x.n())
      throw Failure("Vector<Treal>::"
		    "symv(char const uplo, Treal const, " 
		    "const Matrix<Treal>&, "
		    "const Vector<Treal>&, "
		    "Treal const, const Vector<Treal>&):"
		    "Incorrect dimensions for symmetric "
		    "matrix-vector product");
    if ((A.is_zero() || x.is_zero() || alpha == 0) && 
	(y.is_zero() || beta == 0))
      y = 0;
    else {
      Treal beta_tmp = beta;
      if (y.is_zero()) {
	y.allocate();
	beta_tmp = 0;
      }
      if (A.is_zero() || x.is_zero() || alpha == 0)
	y *= beta_tmp;
      else {
	int const ONE = 1;
	mat::symv(&uplo, &x.n(), &alpha, A.elements, &A.nrows(), 
		  x.elements, &ONE, &beta, y.elements, &ONE);
      } /* end else */
    } /* end else */
  }

  template<class Treal>
    void Vector<Treal>:: 
    trmv(char const uplo, const bool tA,  
	 Matrix<Treal> const & A,
	 Vector<Treal> & x) {
    if (A.nrows() != A.ncols() || A.ncols() != x.n())
      throw Failure("Vector<Treal>::"
		    "trmv(...): Incorrect dimensions for triangular "
		    "matrix-vector product");
    if (uplo != 'U') 
      throw Failure("Vector<class Treal>::"
		    "trmv only implemented for upper triangular matrices");
    if ( ( A.is_zero() || x.is_zero() ) ) {
      x = 0;
      return;
    }
    int const ONE = 1;
    if (!tA)
      mat::trmv(&uplo, "N", "N",  &x.n(), A.elements, &A.nrows(), 
		x.elements, &ONE);
    else
      mat::trmv(&uplo, "T", "N",  &x.n(), A.elements, &A.nrows(), 
		x.elements, &ONE);
  }
  



} /* end namespace mat */
#endif
