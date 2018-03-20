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

/** @file VectorGeneral.h General vector class
 *
 * Copyright(c) Emanuel Rubensson 2006
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date January 2006
 *
 */
#ifndef MAT_VECTORGENERAL
#define MAT_VECTORGENERAL
#include <iostream>
#include <fstream>
#include <ios>
#include "FileWritable.h"
#include "matrix_proxy.h"
#include "ValidPtr.h"
namespace mat {
  template<typename Treal, typename Tvector>
    class VectorGeneral : public FileWritable {
  public:

    inline void resetSizesAndBlocks(SizesAndBlocks const & newRows) {
      vectorPtr.haveDataStructureSet(true);
      vectorPtr->resetRows(newRows);
    }

    inline bool is_empty() const {
      return !vectorPtr.haveDataStructureGet();
    }

    inline void clear_structure(){
      vectorPtr.haveDataStructureSet(false);
    }

    VectorGeneral(SizesAndBlocks const & newRows):vectorPtr(new Tvector) {
      resetSizesAndBlocks(newRows);      
    }
    VectorGeneral():vectorPtr(new Tvector) {}
    explicit VectorGeneral(const VectorGeneral<Treal, Tvector>& other)
      :FileWritable(other), vectorPtr(new Tvector) {
      if (other.vectorPtr.haveDataStructureGet()) {
	vectorPtr.haveDataStructureSet(true);
      }
      *vectorPtr = *other.vectorPtr;
    }
      
    inline void assign_from_full
      (std::vector<Treal> const & fullVector, 
       SizesAndBlocks const & newRows) {
      resetSizesAndBlocks(newRows);
      this->vectorPtr->assignFromFull(fullVector);
    }
    inline void fullvector(std::vector<Treal> & fullVector) const {
      this->vectorPtr->fullVector(fullVector);
    }
    VectorGeneral<Treal, Tvector>& 
      operator=(const VectorGeneral<Treal, Tvector>& other) {
      if (other.vectorPtr.haveDataStructureGet()) {
	vectorPtr.haveDataStructureSet(true);
      }
      *this->vectorPtr = *other.vectorPtr; 
      return *this;
    } 

    inline void clear() {
      if (is_empty())
	// This means that the object's data structure has not been set
	// There is nothing to clear and the vectorPtr is not valid either
	return;
      vectorPtr->clear();
    }

    inline void rand() {
      vectorPtr->randomNormalized();
    }
    
    /* LEVEL 2 operations */
    /* OPERATIONS INVOLVING ORDINARY MATRICES */
    /** y = alpha * op(A) * x */
    template<typename Tmatrix>
      VectorGeneral<Treal, Tvector>& operator=
      (const XYZ<Treal,
       MatrixGeneral<Treal, Tmatrix>,
       VectorGeneral<Treal, Tvector> >& smv);
    
    /** y += alpha * op(A) * x */
    template<typename Tmatrix>
      VectorGeneral<Treal, Tvector>& operator+=
      (const XYZ<Treal,
       MatrixGeneral<Treal, Tmatrix>,
       VectorGeneral<Treal, Tvector> >& smv);
    /** y = alpha * op(A) * x + beta * y */
    template<typename Tmatrix>
      VectorGeneral<Treal, Tvector>& operator=
      (const XYZpUV<Treal,
       MatrixGeneral<Treal, Tmatrix>,
       VectorGeneral<Treal, Tvector>,
       Treal,
       VectorGeneral<Treal, Tvector> >& smvpsv);

    /** y = op(A) * x                      : A is general */
    template<typename Tmatrix>
      VectorGeneral<Treal, Tvector>& operator=
      (const XY<MatrixGeneral<Treal, Tmatrix>,
       VectorGeneral<Treal, Tvector> >& mv) {
      Treal ONE = 1.0;
      return this->operator=(XYZ<Treal, MatrixGeneral<Treal, Tmatrix>,
			     VectorGeneral<Treal, Tvector> >(ONE, mv.A, mv.B,
							     false, mv.tA, mv.tB));      
    }
    
    /* OPERATIONS INVOLVING SYMMETRIC MATRICES */
    /** y = alpha * A * x                      : A is symmetric */
    template<typename Tmatrix>
      VectorGeneral<Treal, Tvector>& operator=
      (const XYZ<Treal,
       MatrixSymmetric<Treal, Tmatrix>,
       VectorGeneral<Treal, Tvector> >& smv);
    /** y += alpha * A * x                     : A is symmetric */
    template<typename Tmatrix>
      VectorGeneral<Treal, Tvector>& operator+=
      (const XYZ<Treal,
       MatrixSymmetric<Treal, Tmatrix>,
       VectorGeneral<Treal, Tvector> >& smv);
    /** y = alpha * A * x + beta * y           : A is symmetric */
    template<typename Tmatrix>
      VectorGeneral<Treal, Tvector>& operator=
      (const XYZpUV<Treal,
       MatrixSymmetric<Treal, Tmatrix>,
       VectorGeneral<Treal, Tvector>,
       Treal,
       VectorGeneral<Treal, Tvector> >& smvpsv);

    /* OPERATIONS INVOLVING TRIANGULAR MATRICES */
    /** y = op(A) * x                      : A is triangular */
    template<typename Tmatrix>
      VectorGeneral<Treal, Tvector>& operator=
      (const XY<MatrixTriangular<Treal, Tmatrix>,
       VectorGeneral<Treal, Tvector> >& mv);


    /* LEVEL 1 operations */
    inline Treal eucl() const {
      return vectorPtr->eucl();
    }

    inline VectorGeneral<Treal, Tvector>& 
      operator*=(Treal const alpha) {
      *vectorPtr *= alpha;
      return *this;
    }

    inline VectorGeneral<Treal, Tvector>& 
      operator=(int const k) {
      *vectorPtr = k;
      return *this;
    }

    /** y += alpha * x */
    VectorGeneral<Treal, Tvector>& operator+=
      (const XY<Treal, VectorGeneral<Treal, Tvector> >& sv);
    

    inline Tvector const & getVector() const {return *vectorPtr;}
    
    std::string obj_type_id() const {return "VectorGeneral";}
  protected:
    ValidPtr<Tvector> vectorPtr;

    inline void writeToFileProt(std::ofstream & file) const {
      if (is_empty())
	// This means that the object's data structure has not been set
	return;
      vectorPtr->writeToFile(file);
    }
    inline void readFromFileProt(std::ifstream & file) {
      if (is_empty())
	// This means that the object's data structure has not been set
	return;
      vectorPtr->readFromFile(file);
    }

    inline void inMemorySet(bool inMem) {
      vectorPtr.inMemorySet(inMem);
    }

  private:
    
  }; /* end class VectorGeneral */


    /* LEVEL 2 operations */
    /* OPERATIONS INVOLVING ORDINARY MATRICES */
    /** y = alpha * op(A) * x */
  template<typename Treal, typename Tvector>
    template<typename Tmatrix>
    VectorGeneral<Treal, Tvector>& 
    VectorGeneral<Treal, Tvector>::operator=
    (const XYZ<Treal,
     MatrixGeneral<Treal, Tmatrix>,
     VectorGeneral<Treal, Tvector> >& smv) {
    assert(!smv.tC);
    vectorPtr.haveDataStructureSet(true);
    if ( this == &smv.C ) {
      // We need a copy of the smv.C vector since it is the same as *this
      VectorGeneral<Treal, Tvector> tmp(smv.C);
      Tvector::gemv(smv.tB, smv.A, smv.B.getMatrix(),
		    *tmp.vectorPtr, 0, *this->vectorPtr);
    }
    else       
      Tvector::gemv(smv.tB, smv.A, smv.B.getMatrix(),
		    *smv.C.vectorPtr, 0, *this->vectorPtr);
    return *this;
  }
  
  /** y += alpha * op(A) * x */
  template<typename Treal, typename Tvector>
    template<typename Tmatrix>
    VectorGeneral<Treal, Tvector>& 
    VectorGeneral<Treal, Tvector>::operator+=
    (const XYZ<Treal,
     MatrixGeneral<Treal, Tmatrix>,
     VectorGeneral<Treal, Tvector> >& smv) {
    assert(!smv.tC);
    assert(this != &smv.C);
    Tvector::gemv(smv.tB, smv.A, smv.B.getMatrix(),
		  *smv.C.vectorPtr, 1, *this->vectorPtr);
    return *this;
  }
  
  
  /** y = alpha * op(A) * x + beta * y */
  template<typename Treal, typename Tvector>
    template<typename Tmatrix>
    VectorGeneral<Treal, Tvector>& 
    VectorGeneral<Treal, Tvector>::operator=
    (const XYZpUV<Treal,
     MatrixGeneral<Treal, Tmatrix>,
     VectorGeneral<Treal, Tvector>,
     Treal,
     VectorGeneral<Treal, Tvector> >& smvpsv) {
    assert(!smvpsv.tC && !smvpsv.tE);
    assert(this != &smvpsv.C);
    if (this == &smvpsv.E)
      Tvector::gemv(smvpsv.tB, smvpsv.A, smvpsv.B.getMatrix(),
		    *smvpsv.C.vectorPtr, smvpsv.D, *this->vectorPtr);
    else
      throw Failure("VectorGeneral<Treal, Tvector>::operator="
		    "(const XYZpUV<Treal, "
		    "MatrixGeneral<Treal, Tmatrix>, "
		    "VectorGeneral<Treal, Tvector>, "
		    "Treal, "
		    "VectorGeneral<Treal, Tvector> >&) : "
		    "y = alpha * op(A) * x + beta * z "
		    "not supported for z != y");
    return *this;
  }

  
  /* OPERATIONS INVOLVING SYMMETRIC MATRICES */
  /** y = alpha * A * x                      : A is symmetric */
  template<typename Treal, typename Tvector>
    template<typename Tmatrix>
    VectorGeneral<Treal, Tvector>& 
    VectorGeneral<Treal, Tvector>::operator=
    (const XYZ<Treal,
     MatrixSymmetric<Treal, Tmatrix>,
     VectorGeneral<Treal, Tvector> >& smv) {
    assert(!smv.tC);
    assert(this != &smv.C);
    vectorPtr.haveDataStructureSet(true);
    Tvector::symv('U', smv.A, smv.B.getMatrix(),
    		  *smv.C.vectorPtr, 0, *this->vectorPtr);
    return *this;    
  }
  

  /** y += alpha * A * x                     : A is symmetric */
  template<typename Treal, typename Tvector>
    template<typename Tmatrix>
    VectorGeneral<Treal, Tvector>& 
    VectorGeneral<Treal, Tvector>::operator+=
    (const XYZ<Treal,
     MatrixSymmetric<Treal, Tmatrix>,
     VectorGeneral<Treal, Tvector> >& smv) {
    assert(!smv.tC);
    assert(this != &smv.C);
    Tvector::symv('U', smv.A, smv.B.getMatrix(),
    		  *smv.C.vectorPtr, 1, *this->vectorPtr);
    return *this;    
  }

  /** y = alpha * A * x + beta * y           : A is symmetric */
  template<typename Treal, typename Tvector>
    template<typename Tmatrix>
    VectorGeneral<Treal, Tvector>& 
    VectorGeneral<Treal, Tvector>::operator=
    (const XYZpUV<Treal,
     MatrixSymmetric<Treal, Tmatrix>,
     VectorGeneral<Treal, Tvector>,
     Treal,
     VectorGeneral<Treal, Tvector> >& smvpsv) {
    assert(!smvpsv.tC && !smvpsv.tE);
    assert(this != &smvpsv.C);
    if (this == &smvpsv.E)
      Tvector::symv('U', smvpsv.A, smvpsv.B.getMatrix(),
		    *smvpsv.C.vectorPtr, smvpsv.D, *this->vectorPtr);
    else
      throw Failure("VectorGeneral<Treal, Tvector>::operator="
		    "(const XYZpUV<Treal, "
		    "MatrixSymmetric<Treal, Tmatrix>, "
		    "VectorGeneral<Treal, Tvector>, "
		    "Treal, "
		    "VectorGeneral<Treal, Tvector> >&) : "
		    "y = alpha * A * x + beta * z "
		    "not supported for z != y");
    return *this;    
  }

  /* OPERATIONS INVOLVING TRIANGULAR MATRICES */
  /** x = op(A) * x                      : A is triangular */

  template<typename Treal, typename Tvector>
    template<typename Tmatrix>
    VectorGeneral<Treal, Tvector>& 
    VectorGeneral<Treal, Tvector>::operator=
    (const XY<MatrixTriangular<Treal, Tmatrix>,
     VectorGeneral<Treal, Tvector> >& mv) {
    assert(!mv.tB);
    if (this != &mv.B)
      throw Failure("y = A * x not supported for y != x ");
    Tvector::trmv('U', mv.tA, 
		  mv.A.getMatrix(), 
		  *this->vectorPtr);
    return *this;
  }
  
  /* LEVEL 1 operations */

  /** y += alpha * x */
  template<typename Treal, typename Tvector>
    VectorGeneral<Treal, Tvector>& 
    VectorGeneral<Treal, Tvector>::operator+=
    (const XY<Treal, VectorGeneral<Treal, Tvector> >& sv) {
    assert(!sv.tB);
    assert(this != &sv.B);
    Tvector::axpy(sv.A, *sv.B.vectorPtr, *this->vectorPtr);
    return *this;
  }



  /* Defined outside class */
  /** transpose(x) * y 
   * Scalar (dot) product of two vectors 
   */
  template<typename Treal, typename Tvector>
    Treal operator*(Xtrans<VectorGeneral<Treal, Tvector> > const & xT,
		    VectorGeneral<Treal, Tvector> const & y) {
    if (xT.tA == false)
      throw Failure("operator*("
		    "Xtrans<VectorGeneral<Treal, Tvector> > const &,"
		    " VectorGeneral<Treal, Tvector> const &): "
		    "Dimension mismatch in vector operation");
    return Tvector::dot(xT.A.getVector(), y.getVector());
  }
  
  

}  /* end namespace mat */
#endif

