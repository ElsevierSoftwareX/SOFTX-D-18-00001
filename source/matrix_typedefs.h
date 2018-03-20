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

/** @file matrix_typedefs.h

\brief Header file with typedefs for matrix and vector types. The levels of hierarchic matrices are defined here.

    @author: Elias Rudberg <em>responsible</em>. 
*/
#ifndef MATRIX_TYPEDEFS_HEADER
#define MATRIX_TYPEDEFS_HEADER


#include "realtype.h"

#include "Matrix.h"
#include "Vector.h"
#include "MatrixSymmetric.h"
#include "MatrixTriangular.h"
#include "MatrixGeneral.h"
#include "VectorGeneral.h"
#include "mat_gblas.h"


// Matrix typedefs for different levels
typedef mat::Matrix<ergo_real, ergo_real> Mat_1;
typedef mat::Matrix<ergo_real, Mat_1    > Mat_2;
typedef mat::Matrix<ergo_real, Mat_2    > Mat_3;
typedef mat::Matrix<ergo_real, Mat_3    > Mat_4;
typedef mat::Matrix<ergo_real, Mat_4    > Mat_5;

// Vector typedefs for different levels
typedef mat::Vector<ergo_real, ergo_real> Vec_1;
typedef mat::Vector<ergo_real, Vec_1    > Vec_2;
typedef mat::Vector<ergo_real, Vec_2    > Vec_3;
typedef mat::Vector<ergo_real, Vec_3    > Vec_4;
typedef mat::Vector<ergo_real, Vec_4    > Vec_5;

// perm and matri are the types actually used.
// This is the point where we choose how many levels to use.
#define MATLEVEL 5
typedef Mat_5   Matri;
typedef Vec_5   Vectorrr;

// The typedefs for symmMatrix, triangMatrix, normalMatrix 
// are the same regardless of how many levels are used.
typedef mat::MatrixSymmetric<ergo_real, Matri> symmMatrix;
typedef mat::MatrixTriangular<ergo_real, Matri> triangMatrix;
typedef mat::MatrixGeneral<ergo_real, Matri> normalMatrix;
typedef mat::VectorGeneral<ergo_real, Vectorrr> generalVector;

typedef mat::Interval<ergo_real> intervalType;



#endif
