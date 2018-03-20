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

/** @file matrix_typedefs_chtml.h

    @brief Header file with typedefs for matrix types, using either
    the hierarchical matrix library (HML) or the Chunks and Tasks
    matrix library (CHTML).

    @author Anastasia Kruchinina <em>responsible</em>
*/
#ifndef MATRIX_TYPEDEFS_CHTML_HEADER
#define MATRIX_TYPEDEFS_CHTML_HEADER

#include "matrix_typedefs.h"

#ifdef USE_CHUNKS_AND_TASKS
//#define USE_SYMMETRIC // TODO
#include "CHTGeneralMatrix.h"
#include "CHTSymmMatrix.h"
#include "CHTTriangMatrix.h"
#ifdef USE_CHUNKS_AND_TASKS_BSM
#include "block_sparse_matrix_lib.h"
typedef bsm::BlockSparseMatrix<ergo_real>         LeafMatrixType;
#else
#include "basic_matrix_lib.h"
typedef bml::FullMatrix<ergo_real>                LeafMatrixType;
#endif
typedef chtml::CHTMatrixParamsType<ergo_real>            ParamsType;
#ifdef USE_SYMMETRIC
typedef chtml::CHTSymmMatrix<ergo_real, ParamsType>      symmMatrixWrap;
#else
typedef chtml::CHTGeneralMatrix<ergo_real, ParamsType>   symmMatrixWrap;
#endif
typedef chtml::CHTGeneralMatrix<ergo_real, ParamsType>   normalMatrixWrap;
typedef chtml::CHTTriangMatrix<ergo_real, ParamsType>    triangMatrixWrap;
#else // not CHT
class MatrixParamsType
{
public:
  MatrixParamsType(){}
};
typedef MatrixParamsType                          ParamsType;
typedef symmMatrix                                symmMatrixWrap;
typedef normalMatrix                              normalMatrixWrap;
typedef triangMatrix                              triangMatrixWrap;
#endif




#endif
