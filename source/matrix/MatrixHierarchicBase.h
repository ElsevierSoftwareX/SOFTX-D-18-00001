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

/** @file MatrixHierarchicBase.h Base class for Matrix
 *
 * Copyright(c) Emanuel Rubensson 2005
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date September 2005
 *
 */


#ifndef MAT_MATRIXHIERARCHICBASE
#define MAT_MATRIXHIERARCHICBASE
#include "matInclude.h"
#include "allocate.h"
namespace mat{
  /** Base class for Matrix and Matrix specialization
   *
   * @see Matrix
   * @see Permutation
   *
   */
  template<class Treal, class Telement = Treal>
    class MatrixHierarchicBase {
    public:  
    /* No public constructors (!) */
    inline bool operator==(int k) const {
      if (k == 0)
	return this->is_zero();
      else 
	throw Failure("Matrix::operator== only implemented for k == 0");
    }
    /* Check if matrix is zero (k = 0)                                     */
#if 1
    inline const int& nScalarsRows() const
    {return rows.getNScalars();}
    inline const int& nScalarsCols() const
    {return cols.getNScalars();}
#endif

    inline const int& nrows() const   /* Number of rows in Telement matrix */
    {return rows.getNBlocks();} 
    inline const int& ncols() const   /* Number of cols in Telement matrix */
    {return cols.getNBlocks();} 

    inline Telement& operator() /* Returns the element A(row, col)        */
    (int row, int col) {
      assert(elements);
      assert(row >= 0);
      assert(col >= 0);
      assert(row < nrows());
      assert(col < ncols());
      return elements[row + col * nrows()];
    }   
    inline const Telement& operator() /*Write protected reference returned*/
    (int row, int col) const {
      assert(elements);
      assert(row >= 0);
      assert(col >= 0);
      assert(row < nrows());
      assert(col < ncols());
      return elements[row + col * nrows()];
    }

    inline Telement& operator[] 
    (int index) {
      assert(elements);
      assert(index >= 0);
      assert(index < nElements());
      return elements[index];
    }   
    inline Telement const & operator[] 
    (int index) const {
      assert(elements);
      assert(index >= 0);
      assert(index < nElements());
      return elements[index];
    }   

    inline bool is_zero() const {return !elements;}

    inline int nElements() const {
      return rows.getNBlocks() * cols.getNBlocks();
    }

    inline void resetRows(SizesAndBlocks const & newRows) {
      freeElements(elements);
      elements = 0;
      rows = newRows;
    }
    inline void resetCols(SizesAndBlocks const & newCols) {
      freeElements(elements);
      elements = 0;
      cols = newCols;
    }

    inline void getRows(SizesAndBlocks & rowsCopy) const {
      rowsCopy = rows;
    }
    inline void getCols(SizesAndBlocks & colsCopy) const {
      colsCopy = cols;
    }


    inline bool highestLevel() const {
      return (rows.getNTotalScalars() == rows.getNScalars() &&
	      cols.getNTotalScalars() == cols.getNScalars());
    }

    /** Check if matrix is empty
     *  Empty is different from zero, a zero matrix contains information
     *  about blocksizes etc.  
     *
     */
    inline bool is_empty() const {
      return rows.is_empty() || cols.is_empty();
    }
    protected:
    
    
    MatrixHierarchicBase() 
    : elements(0) {}
    MatrixHierarchicBase(SizesAndBlocks const & rowsInp, 
			 SizesAndBlocks const & colsInp) 
    : rows(rowsInp), cols(colsInp), elements(0) {}
    MatrixHierarchicBase(const MatrixHierarchicBase<Treal, Telement>& mat);
    
    
    MatrixHierarchicBase<Treal, Telement>& 
    operator=(const MatrixHierarchicBase<Treal, Telement>& mat);
    
    static void swap(MatrixHierarchicBase<Treal, Telement>& A,
		     MatrixHierarchicBase<Treal, Telement>& B);      
    
    virtual ~MatrixHierarchicBase();      
    SizesAndBlocks rows;
    SizesAndBlocks cols;
    Telement* elements; /* Length is nRows * nCols unless 0 */
    
    
    
#if 0
    inline void assert_alloc() {
      if (this->cap < this->nel) {
	freeElements(this->elements);
	this->cap = this->nel;      
	this->elements = allocateElements<Telement>(this->cap);
	for (int ind = 0; ind < this->cap; ind++)
	  this->elements[ind] = 0;
      }
    } 
#endif
    private:
    
  }; /* end class */



  template<class Treal, class Telement> /* Copy constructor     */
    MatrixHierarchicBase<Treal, Telement>::
    MatrixHierarchicBase(const MatrixHierarchicBase<Treal, Telement>& mat) 
    : rows(mat.rows), cols(mat.cols), elements(0) {
    if (!mat.is_zero()) {
      elements = allocateElements<Telement>(nElements());
      for (int i = 0; i < nElements(); i++) 
	elements[i] = mat.elements[i];
    }
  }
    

    template<class Treal, class Telement> /*Assignment operator*/
      MatrixHierarchicBase<Treal, Telement>& 
      MatrixHierarchicBase<Treal, Telement>::
      operator=(const MatrixHierarchicBase<Treal, Telement>& mat) {
      if (mat.is_zero()) { /* Condition also matches empty matrices. */
	rows = mat.rows;
	cols = mat.cols;
	freeElements(elements);
	elements = 0;
	return *this;
      }
      if (is_zero() || (nElements() != mat.nElements())) {
	freeElements(elements);
	elements = allocateElements<Telement>(mat.nElements());
      }
      rows = mat.rows;
      cols = mat.cols;
      for (int i = 0; i < nElements(); i++)
	elements[i] = mat.elements[i];
      return *this;
    }

    template<class Treal, class Telement>
      void MatrixHierarchicBase<Treal, Telement>::
      swap(MatrixHierarchicBase<Treal, Telement>& A,
	   MatrixHierarchicBase<Treal, Telement>& B) {
      assert(A.rows == B.rows && A.cols == B.cols);
      Telement* elementsTmp = A.elements;
      A.elements = B.elements;
      B.elements = elementsTmp;
    }

  
    template<class Treal, class Telement>
      MatrixHierarchicBase<Treal, Telement>::~MatrixHierarchicBase() {
      freeElements(elements);
      elements = 0;
    }
  
     
} /* end namespace mat */

#endif
