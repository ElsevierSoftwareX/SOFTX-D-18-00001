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

/** @file VectorHierarchicBase.h Base class for Vector
 *
 * Copyright(c) Emanuel Rubensson 2006
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date October 2006
 *
 */


#ifndef MAT_VECTORHIERARCHICBASE
#define MAT_VECTORHIERARCHICBASE
#include "matInclude.h"
namespace mat{
  /** Base class for Vector and Vector specialization
   *
   * @see Vector
   * @see Permutation
   *
   */
  template<class Treal, class Telement = Treal>
    class VectorHierarchicBase {
    public:

#if 1
    inline const int& nScalars() const
    {return rows.getNScalars();}
#endif
    inline const int& n() const   /* Number of elements in Telement matrix */
    {return rows.getNBlocks();} 

    inline Telement& operator() /* Returns the element v(ind)        */
    (int ind) {
      assert(elements);
      assert(ind >= 0);
      assert(ind < n());
      return elements[ind];
    }   
    inline const Telement& operator() /*Write protected reference returned*/
    (int ind) const {
      assert(elements);
      assert(ind >= 0);
      assert(ind < n());
      return elements[ind];
    }
    inline bool is_zero() const {return !elements;}
    
    inline void resetRows(SizesAndBlocks const & newRows) {
      freeElements(elements);
      elements = 0;
      rows = newRows;
    }

    protected:
    /** Check if vector is empty
     *  Empty is different from zero, a zero matrix contains information
     *  about blocksizes etc.  
     *
     */
    inline bool is_empty() const {
      return rows.is_empty();
    }

    VectorHierarchicBase()
    : elements(0) {}

    explicit VectorHierarchicBase(SizesAndBlocks const & rowsInp)
    :elements(0) {}
    VectorHierarchicBase
    (const VectorHierarchicBase<Treal, Telement>& vec);
    VectorHierarchicBase<Treal, Telement>& 
    operator=(const VectorHierarchicBase<Treal, Telement>& vec);
    virtual ~VectorHierarchicBase();      
    
    
    SizesAndBlocks rows;

    //    const Tperm* perm;
    Telement* elements;
    //    int cap;             /* The length of the elements array         */
    //    int nel;             /* Number of USED elements in the elements  */
    /*                 array (can be positive even though content == zero) */

    //    property content; /* content can be one of the properties listed */
    /* in the enum type property (matInclude.h) for example:               */
    /* zero: An all zero matrix                                            */
    /* ful : An ordinary matrix with the values in the elements array      */

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

  }; /* end class VectorHierarchicBase */

  template<class Treal, class Telement> /* Copy constructor     */
    VectorHierarchicBase<Treal, Telement>::
    VectorHierarchicBase
    (const VectorHierarchicBase<Treal, Telement>& vec) 
    : rows(vec.rows) {
      if (!vec.is_zero()) {
	elements = allocateElements<Telement>(n());
	for (int i = 0; i < n(); i++) 
	  elements[i] = vec.elements[i];
      }
    }   
  

    template<class Treal, class Telement> /* Assignment operator*/
      VectorHierarchicBase<Treal, Telement>& 
      VectorHierarchicBase<Treal, Telement>::
      operator=(const VectorHierarchicBase<Treal, Telement>& vec) {
      if (vec.is_zero()) { /* Condition also matches empty matrices. */
	rows = vec.rows;
	freeElements(elements);
	elements = 0;
	return *this;
      }
      if (is_zero() || (n() != vec.n())) {
	freeElements(elements);
	elements = allocateElements<Telement>(vec.n());
      }
      rows = vec.rows;
      for (int i = 0; i < n(); i++)
	elements[i] = vec.elements[i];
      return *this;
    }

  
    template<class Treal, class Telement>
      VectorHierarchicBase<Treal, Telement>::
      ~VectorHierarchicBase() {
      freeElements(elements);
    }
    
} /* end namespace mat */
#endif
