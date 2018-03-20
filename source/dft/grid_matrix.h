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

/** @file grid_matrix.h

    @brief Generic matrix interface. It is not optimized for speed.

    @author: Pawel Salek <em>responsible</em>
*/

#if !defined(_GRID_MATRIX_H_)
#define _GRID_MATRIX_H_ 1

#include "sparse_matrix.h"

namespace Dft {

  class Matrix {
  public:
    virtual ergo_real at(int row, int col) const = 0;
    virtual bool isSparse() const = 0;
    virtual const SparseMatrix* asSparse() const = 0;
    virtual const ergo_real* asFull() const = 0;
    virtual ~Matrix() {}
  };

  class FullMatrix {
  public:
    ergo_real* mat;
    int nbast;
    bool owned;
  explicit FullMatrix(int nbast_)
    : mat(new ergo_real[nbast_*nbast_]), nbast(nbast_), owned(true)
      {
        for(int i= nbast*nbast-1; i >=0; --i) mat[i] = 0.0;
      }
  FullMatrix(ergo_real *m, int nbast_)
    : mat(m), nbast(nbast_), owned(false)
      {
      }
    /** ugly-hack constructor. Remove it! */
  FullMatrix(const ergo_real *m, int nbast_)
    : mat( (ergo_real*)(m)), nbast(nbast_), owned(false)
      {
      }

    ~FullMatrix() { if (owned && mat) delete []mat; }
    void add(int row, int col, ergo_real val)
    { 
      mat[row + col*nbast] += val;
    }
    ergo_real at(int row, int col) const
    {
      return mat[row + col*nbast];
    }
  };

}

#endif /* _GRID_MATRIX_H_ */
