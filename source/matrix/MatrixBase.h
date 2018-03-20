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

/** @file MatrixBase.h Base class for matrix API
 *
 * Copyright(c) Emanuel Rubensson 2006
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date January 2006
 *
 */
#ifndef MAT_MATRIXBASE
#define MAT_MATRIXBASE
#include <iostream>
#include <fstream>
#include <ios>
#include "FileWritable.h"
#include "matrix_proxy.h"
#include "ValidPtr.h"
#include "SizesAndBlocks.h"
namespace mat {
  template<typename Treal, typename Tmatrix>
    class MatrixGeneral; 
  template<typename Treal, typename Tmatrix>
    class MatrixSymmetric;  
  template<typename Treal, typename Tmatrix>
    class MatrixTriangular; 
  template<typename Treal, typename Tvector>
    class VectorGeneral;
  enum matrix_type {matrix_matr, matrix_symm, matrix_triang};
  
  /** Base class for matrix API
   *
   * This class provides a base for an API to a matrix library built 
   * up from three types which are also the template arguments to this class.
   * 
   * Treal: Type for real numbers
   *
   * Tmatrix: The matrix class
   *
   */  
  template<typename Treal, typename Tmatrix>
    class MatrixBase : public FileWritable {
  public:
    friend class MatrixGeneral<Treal, Tmatrix>;
    friend class MatrixSymmetric<Treal, Tmatrix>;
    friend class MatrixTriangular<Treal, Tmatrix>;
    

    inline void resetSizesAndBlocks(SizesAndBlocks const & newRows,
                                    SizesAndBlocks const & newCols) {
      matrixPtr.haveDataStructureSet(true);
      matrixPtr->resetRows(newRows);
      matrixPtr->resetCols(newCols);
    }
    inline void getRows(SizesAndBlocks & rowsCopy) const {
      matrixPtr->getRows(rowsCopy);
    }
    inline void getCols(SizesAndBlocks & colsCopy) const {
      matrixPtr->getCols(colsCopy);
    }

    /** Check if matrix is empty.
	Being empty is not the same as being zero.
	A matrix being empty means that the data structure has not been set. 
     */
    inline bool is_empty() const {
      return !matrixPtr.haveDataStructureGet();
    }

    inline Treal trace() const {
      return matrixPtr->trace();
    }
    
    inline void add_identity(Treal alpha) {
      matrixPtr->addIdentity(alpha);
    }
    inline MatrixBase<Treal, Tmatrix>& operator*=(Treal const alpha) {
      *matrixPtr *= alpha;
      return *this;
    }

    inline bool operator==(int k) const {
      if (k == 0)
	return *matrixPtr == 0;
      else 
	throw Failure("MatrixBase::operator== only implemented for k == 0");
    }

    
    
    inline void clear() {
      if (is_empty())
	// This means that the object's data structure has not been set
	// There is nothing to clear and the matrixPtr is not valid either
	return;
      matrixPtr->clear();
    }
    
    inline size_t memory_usage() const {
      return matrixPtr->memory_usage();
    }

    inline void write_to_buffer_count(int& n_bytes) const {
      int ib_length = 3;
      int vb_length = 0;
      this->matrixPtr->write_to_buffer_count(ib_length, vb_length);
      n_bytes = vb_length * sizeof(Treal) + ib_length * sizeof(int);
    }

#if 1
    inline int get_nrows() const {
      return matrixPtr->nScalarsRows();
    }
    inline int get_ncols() const {
      return matrixPtr->nScalarsCols();
    }
#endif
    
    inline Tmatrix const & getMatrix() const {return *matrixPtr;}
    inline Tmatrix & getMatrix() {return *matrixPtr;}
 
    /** Get largest absolute value of matrix element in the matrix. */
    inline Treal maxAbsValue() const {return matrixPtr->maxAbsValue();}

  protected:
    ValidPtr<Tmatrix> matrixPtr;
    
    MatrixBase():matrixPtr(new Tmatrix) {}
    MatrixBase(const MatrixBase<Treal, Tmatrix>& other)
      :FileWritable(other), matrixPtr(new Tmatrix) {
      matrixPtr.haveDataStructureSet(other.matrixPtr.haveDataStructureGet());
      /* getConstRefForCopying() is used here to make sure it works
	 also in the case when the matrix is written to file. */
      *matrixPtr = other.matrixPtr.getConstRefForCopying();
      matrixPtr.inMemorySet(other.matrixPtr.inMemoryGet());
    }
      
    MatrixBase<Treal, Tmatrix>& 
      operator=(const MatrixBase<Treal, Tmatrix>& other) {
      FileWritable::operator=(other); /* Allows us to copy mat on file */
      matrixPtr.haveDataStructureSet(other.matrixPtr.haveDataStructureGet());
      /* getConstRefForCopying() is used here to make sure it works
	 also in the case when the matrix is written to file. */
      *matrixPtr = other.matrixPtr.getConstRefForCopying();
      matrixPtr.inMemorySet(other.matrixPtr.inMemoryGet());
      return *this;
    } 
    
    MatrixBase<Treal, Tmatrix>& 
      operator=(const Xtrans<MatrixGeneral<Treal, Tmatrix> >& mt) {
      if (mt.A.matrixPtr.haveDataStructureGet()) {
	matrixPtr.haveDataStructureSet(true);
      }
      if (mt.tA)
	Tmatrix::transpose(*mt.A.matrixPtr, *this->matrixPtr);
      else
	*this->matrixPtr = *mt.A.matrixPtr; 
      return *this;
      // FileWritable::operator=(other);/*Could be used to copy mat on file*/
    } 
    

    void write_to_buffer_base(void* buffer, const int n_bytes,
			 const matrix_type mattype) const;     
    void read_from_buffer_base(void* buffer, const int n_bytes,
			  const matrix_type mattype);

    void writeToFileBase(std::ofstream & file, 
			 matrix_type const mattype) const;
    void readFromFileBase(std::ifstream & file, 
			  matrix_type const mattype);

    std::string obj_type_id() const {return "MatrixBase";}
    inline void inMemorySet(bool inMem) {
      matrixPtr.inMemorySet(inMem);
    }

    static void getPermutedIndexes(std::vector<int> const & index, 
				   std::vector<int> const & permutation,
				   std::vector<int> & newIndex) {
      newIndex.resize(index.size());
      for (unsigned int i = 0; i < index.size(); ++i) 
	newIndex[i] = permutation[index[i]];
    }


  private:
    
  };

   
  template<typename Treal, typename Tmatrix>
    void MatrixBase<Treal, Tmatrix>::
    writeToFileBase(std::ofstream & file, 
		matrix_type const mattype) const {
    int type = (int)mattype;
    file.write((char*)&type,sizeof(int));
    
    if (is_empty())
      // This means that the object's data structure has not been set
      // The ValidPtr prevents setting the data structure between
      // calls to writeToFile and readFromFile 
      return;
    matrixPtr->writeToFile(file);
  }
  
  template<typename Treal, typename Tmatrix>
    void MatrixBase<Treal, Tmatrix>::
    readFromFileBase(std::ifstream & file, 
		 matrix_type const mattype) {
    char type[sizeof(int)];
    file.read(type, sizeof(int));
    if (((int)*type) != mattype)
      throw Failure("MatrixBase<Treal, Tmatrix>::"
		    "readFromFile(std::ifstream &, " 
		    "matrix_type const): Wrong matrix type");
    if (is_empty())
      // This means that the object's data structure has not been set
      return;
    matrixPtr->readFromFile(file);
  }



  template<typename Treal, typename Tmatrix>
    void MatrixBase<Treal, Tmatrix>::
    write_to_buffer_base(void* buffer, const int n_bytes,
			 const matrix_type mattype) const {
    int ib_length = 3; /* Length of integer buffer, at least 3: matrix_type, */
    /*                    ib_length and vb_length                            */
    int vb_length = 0; /* Length of value buffer                             */
    this->matrixPtr->write_to_buffer_count(ib_length, vb_length);
    if (n_bytes >= 
	(int)(vb_length * sizeof(Treal) + ib_length * sizeof(int))) {
      int* int_buf = (int*)buffer;
      int_buf[0] = mattype;
      int_buf[1] = ib_length;
      int_buf[2] = vb_length;
      Treal* value_buf = (Treal*)&(int_buf[ib_length]); /* Value buffer      */
      /* begins after integer buffer end                                     */
      int ib_index = 0; 
      int vb_index = 0;
      this->matrixPtr->write_to_buffer(&int_buf[3], ib_length - 3, 
				       value_buf, vb_length,
				       ib_index, vb_index);
    }
    else {
      throw Failure("MatrixBase::write_to_buffer: Buffer is too small");
    }
  }
    
  template<typename Treal, typename Tmatrix>
    void MatrixBase<Treal, Tmatrix>::
    read_from_buffer_base(void* buffer, const int n_bytes,
		     const matrix_type mattype) {
    int* int_buf = (int*)buffer;
    if(int_buf[0] == mattype) {
      int ib_length = int_buf[1];
      int vb_length = int_buf[2];
      int ib_index = 0; 
      int vb_index = 0;
      Treal* value_buf = (Treal*)&(int_buf[ib_length]);
      this->matrixPtr->read_from_buffer(&int_buf[3], ib_length - 3,
					value_buf, vb_length,
					ib_index, vb_index);
    }
    else {
      throw Failure("MatrixBase::read_from_buffer: Wrong matrix type");
    }
  }


} /* end namespace mat */
#endif


