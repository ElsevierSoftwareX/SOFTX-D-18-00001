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

#if !defined(_SPARSE_MATRIX_H_)
#define _SPARSE_MATRIX_H_ 1

/** @file sparse_matrix.h Declares a sparse matrix optimized for the
    XC code.  The object provides methods for fast preallocation of
    the matrix elements, and some matrix elements iterators as needed
    for the numerical matrix element integration schemes. */

#include <stdio.h>
#include <vector>
#include <algorithm>


#include "realtype.h"
#include "matrix_typedefs.h"
#include "basisinfo.h"
#include "sparse_pattern.h"

#if !defined(BEGIN_NAMESPACE)
#define BEGIN_NAMESPACE(x) namespace x {
#define END_NAMESPACE(x)   } /* x */
#endif

BEGIN_NAMESPACE(Dft)

/** Sparse matrix structure optimized for XC data access pattern. */
class SparseMatrix {
  class Exception : public std::exception {
    const char *msg;
  public:
  explicit Exception(const char *msg_) : msg(msg_) {}
    virtual const char *what() const throw() { return msg; }
  };

  const SparsePattern& pattern;
  ergo_real **columns;
  int **offsets; /**< for accelerated at() and add() methods. */
  int **his;     /**< for accelerated at() and add() methods. */
  int *cnt;      /**< for accelerated at() and add() methods. */
  int n;
 /** Fills in offsets and his based on pattern. */
  void createOffsets(const SparsePattern& pattern);
 public:
  /** Constructs a square matrix and preallocate according to the
      specified pattern. */ 
  explicit SparseMatrix(const SparsePattern& pattern_);
  SparseMatrix(const SparsePattern& pattern_,
               const symmMatrix& m, const int *aoMap,
               std::vector<int> const & permutationHML);

  ~SparseMatrix() {
    for(int i=0; i<n; i++) {
      delete [](columns[i]);
      delete [](offsets[i]);
      delete [](his[i]);
    }
    delete []columns;
    delete []offsets;
    delete []his;
    delete []cnt;
  }

  void print(const char *title) const;

  /** Assigns itself to a given hierarchic matrix. */
  void addSymmetrizedTo(symmMatrix& sMat,
                        const int *aoMap,
                        std::vector<int> const & permutationHML) const;

  /** Adds given value to an element in given row and column.
      Checking against intervals.end() is *terribly* expensive!!!
      Luckily, we do not have to do it.*/
  void add(int row, int col, ergo_real val) {
    ergo_real *columnData = columns[col];
    const int *hi = his[col];
    int idx;
    for(idx = 0; idx < cnt[col] && row >hi[idx]; ++idx);
    //int idx = std::upper_bound(hi, hi+cnt[col], row)-hi;
    if(idx >= cnt[col])
      throw Exception("SparseMatrix::add called with incorrect args");
    int offset = offsets[col][idx];
    /* Add it... */
    //assert(row-offset>=0);
    //assert(row-offset<pattern.getColumnSize(col));
    columnData[row-offset] += val;
  }

  /* This operator[] syntax glue that we could in principle use can be
     expensive performance-wise, do it the old-fashioned way.
     Checking against intervals.end() is *terribly* expensive!!!
  */
  ergo_real at(int row, int col) const {
    const ergo_real *columnData = columns[col];
    const int *hi = his[col]; 
    int idx; for(idx = 0; idx < cnt[col] && row >hi[idx]; ++idx);
    if(idx >= cnt[col])
      throw Exception("SparseMatrix::at called with incorrect args");
    //int idx = std::upper_bound(hi, hi+cnt[col], row)-hi;
    int offset = offsets[col][idx];
    /* return it... */
    //assert(row-offset>=0);
    //assert(row-offset<pattern.getColumnSize(col));
    return columnData[row-offset];
  }
};


END_NAMESPACE(Dft)

void
getrho_blocked_lda(int nbast, const Dft::SparseMatrix& dmat,
                   const ergo_real * gao,
                   const int* nblocks, const int (*iblocks)[2],
                   int ldaib, ergo_real *tmp, int nvclen, ergo_real *rho);
void
getrho_blocked_gga(int nbast, const Dft::SparseMatrix& dmat,
                   const ergo_real * gao,
                   const int* nblocks, const int (*iblocks)[2],
                   int ldaib, ergo_real *tmp, int nvclen,
                   ergo_real *rho, ergo_real (*grad)[3]);

#endif /* _SPARSE_MATRIX_H_ */
