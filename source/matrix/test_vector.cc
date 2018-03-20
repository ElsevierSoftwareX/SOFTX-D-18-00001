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

/** @file test_vector.cc

    \brief Test code testing creation of vector objects using the
    hierarchical vector data structure in the hierarchical matrix
    library.
*/

#include <fstream>  /* For ifstream */
#include <iomanip> /* For setprecision in fstream */
#include <iostream>
#include <cmath>
#include <stdio.h> /* For FILE */
#include "SizesAndBlocks.h"
#include "Matrix.h"
#include "Vector.h"
#include "MatrixSymmetric.h"
#include "MatrixTriangular.h"
#include "MatrixGeneral.h"
#include "VectorGeneral.h"
#include "mat_gblas.h"
#include "Lanczos.h"

using namespace mat;

typedef double real;
typedef Vector<real, real> Vec_3;
//typedef Vector<real, Vec_1> Vec_2;
//typedef Vector<real, Vec_2> Vec_3;

typedef Vec_3 vect;
typedef VectorGeneral<real, vect> normalVector;

int main()
{
  int size=87;
  int nlevels=1;
  std::vector<int> blockSizes(nlevels);
  blockSizes[nlevels - 1] = 1;  
  SizesAndBlocks rows(blockSizes, size);
  std::vector<normalVector> V(2,normalVector(rows));
  return 0;
}
