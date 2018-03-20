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

/** @file mat_utils.h

    \brief Utilities used by other parts of the hierarchical matrix
    library.
*/

#ifndef MAT_UTILS_HEADER
#define MAT_UTILS_HEADER
#include "Interval.h"
#include "matrix_proxy.h"
namespace mat {

  template<typename Tmatrix, typename Treal>
    struct DiffMatrix {
      typedef typename Tmatrix::VectorType VectorType;
      void getCols(SizesAndBlocks & colsCopy) const {
	A.getCols(colsCopy);
      }
      int get_nrows() const { 
	assert( A.get_nrows() == B.get_nrows() );
	return A.get_nrows(); 
      }
      Treal frob() const {
	return Tmatrix::frob_diff(A, B);
      }
      void quickEuclBounds(Treal & euclLowerBound, 
			   Treal & euclUpperBound) const {
	Treal frobTmp = frob();
	euclLowerBound = frobTmp  / template_blas_sqrt( (Treal)get_nrows() );
	euclUpperBound = frobTmp;
      }

      Tmatrix const & A;
      Tmatrix const & B;
      DiffMatrix(Tmatrix const & A_, Tmatrix const & B_)
      : A(A_), B(B_) {}
      template<typename Tvector>
      void matVecProd(Tvector & y, Tvector const & x) const {
	Tvector tmp(y);
	tmp = (Treal)-1.0 * B * x;   // -B * x
	y   = (Treal)1.0 * A * x;    // A * x
	y  += (Treal)1.0 * tmp;        // A * x - B * x  => (A - B) * x
      }
    };


  // ATAMatrix AT*A 
  template<typename Tmatrix, typename Treal>
    struct ATAMatrix {
      typedef typename Tmatrix::VectorType VectorType;
      Tmatrix const & A;
      explicit ATAMatrix(Tmatrix const & A_)
      : A(A_) {}
      void getCols(SizesAndBlocks & colsCopy) const {
	A.getRows(colsCopy);
      }
      void quickEuclBounds(Treal & euclLowerBound, 
			   Treal & euclUpperBound) const {
	Treal frobA = A.frob();
	euclLowerBound = 0;
	euclUpperBound = frobA * frobA;
      }
      
      // y = AT*A*x
      template<typename Tvector>
      void matVecProd(Tvector & y, Tvector const & x) const {
	y = x;
	y = A * y;
	y = transpose(A) * y;
      }
      // Number of rows of A^T * A is the number of columns of A 
      int get_nrows() const { return A.get_ncols(); }       
    };


  template<typename Tmatrix, typename Tmatrix2, typename Treal>
    struct TripleMatrix {
      typedef typename Tmatrix::VectorType VectorType;
      void getCols(SizesAndBlocks & colsCopy) const {
	A.getCols(colsCopy);
      }
      int get_nrows() const { 
	assert( A.get_nrows() == Z.get_nrows() );
	return A.get_nrows(); 
      }
      void quickEuclBounds(Treal & euclLowerBound, 
			   Treal & euclUpperBound) const {
	Treal frobA = A.frob();
	Treal frobZ = Z.frob();
	euclLowerBound = 0;
	euclUpperBound = frobA * frobZ * frobZ;
      }
      
      Tmatrix  const & A;
      Tmatrix2 const & Z;
      TripleMatrix(Tmatrix const & A_, Tmatrix2 const & Z_)
      : A(A_), Z(Z_) {}
      void matVecProd(VectorType & y, VectorType const & x) const {
	VectorType tmp(x);
	tmp = Z * tmp;            // Z * x
	y = (Treal)1.0 * A * tmp; // A * Z * x
	y = transpose(Z) * y;     // Z^T * A * Z * x
      }
    };


  template<typename Tmatrix, typename Tmatrix2, typename Treal>
    struct CongrTransErrorMatrix {
      typedef typename Tmatrix::VectorType VectorType;
      void getCols(SizesAndBlocks & colsCopy) const {
	E.getRows(colsCopy);
      }
      int get_nrows() const { 
	return E.get_ncols(); 
      }
      void quickEuclBounds(Treal & euclLowerBound, 
			   Treal & euclUpperBound) const {
	Treal frobA = A.frob();
	Treal frobZ = Zt.frob();
	Treal frobE = E.frob();
	euclLowerBound = 0;
	euclUpperBound = frobA * frobE * frobE + 2 * frobA * frobE * frobZ;
      }
      
      Tmatrix  const & A;
      Tmatrix2 const & Zt;
      Tmatrix2 const & E;
      
      CongrTransErrorMatrix(Tmatrix const & A_, 
			    Tmatrix2 const & Z_,
			    Tmatrix2 const & E_)
      : A(A_), Zt(Z_), E(E_) {}
      void matVecProd(VectorType & y, VectorType const & x) const {
	
	VectorType tmp(x);
	tmp = E * tmp;               // E * x
	y = (Treal)-1.0 * A * tmp;   // -A * E * x
	y = transpose(E) * y;        // -E^T * A * E * x
	
	VectorType tmp1;
	tmp = x;
	tmp = Zt * tmp;              // Zt * x
	tmp1 = (Treal)1.0 * A * tmp; // A * Zt * x
	tmp1 = transpose(E) * tmp1;  // E^T * A * Zt * x
	y += (Treal)1.0 * tmp1;

	tmp = x;
	tmp = E * tmp;               // E * x
	tmp1 = (Treal)1.0 * A * tmp; // A * E * x
	tmp1 = transpose(Zt) * tmp1; // Zt^T * A * E * x
	y += (Treal)1.0 * tmp1;	
      }
    };



}  /* end namespace mat */
#endif
