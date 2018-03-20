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

/** @file source/test/registration.cc

    @brief Registration of task types and chunk types, relevant only
    if configuring ergo to use Chunks and Tasks.
*/

#include "registration.h"

#ifdef USE_CHUNKS_AND_TASKS


/* CHTTL registration stuff */
CHTTL_REGISTER_ALL_FOR_FUND_ARITHMETIC_TYPES;


/* CHTML registration stuff */
CHTML_REGISTER_CHUNK_TYPE((chtml::Matrix<LeafMatrixType >));
CHTML_REGISTER_CHUNK_TYPE((chtml::MatrixParams<LeafMatrixType >));
// 
CHTML_REGISTER_TASK_TYPE((chtml::MatrixGetElements<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixMultiply<LeafMatrixType, false, true >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixMultiply<LeafMatrixType, false, false >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixMultiply<LeafMatrixType, true, true >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixMultiply<LeafMatrixType, true, false >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixSymmMultiply<LeafMatrixType, true, false>));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixSymmMultiply<LeafMatrixType, false, true>));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixSymmSquare<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixSymmRK<LeafMatrixType, true >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixSymmRK<LeafMatrixType, false >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixRescale<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixNNZ<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixMaxAbsElement<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixAssignFromChunkIDs<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixCombineElements<double>));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixAdd<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixAddScaledIdentity<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixAssignToScaledIdentity<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixAssignFromSparse<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixSquaredFrobOfErrorMatrix<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixNormFrobenius<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixFrobTruncLowestLevel<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixInvChol<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixInvCholWithTrunc<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixTranspose<LeafMatrixType >));

CHTML_REGISTER_CHUNK_TYPE((chtml::ChunkMatrixElements<real>));
CHTML_REGISTER_TASK_TYPE((chtml::ChunkMatrixElementsAppend<real>));
CHTML_REGISTER_TASK_TYPE((chtml::GetSubmatrixElements<LeafMatrixType>));
CHTML_REGISTER_TASK_TYPE((chtml::SetSubmatrix<LeafMatrixType>));
CHTML_REGISTER_TASK_TYPE((chtml::SetSubmatrix_B_null<LeafMatrixType>));
CHTML_REGISTER_TASK_TYPE((chtml::SetSubmatrix_A_null<LeafMatrixType>));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixTrace<LeafMatrixType >));

CHTML_REGISTER_TASK_TYPE((chtml::MatrixColsSumsPart<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixRowsSumsPart<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixGetDiagPart<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::CombineVectors<chttl::ChunkVector<real> >));
CHTML_REGISTER_TASK_TYPE((chtml::CombineVectors_a1_null<chttl::ChunkVector<real> >));
CHTML_REGISTER_TASK_TYPE((chtml::CombineVectors_a2_null<chttl::ChunkVector<real> >));
CHTML_REGISTER_TASK_TYPE((chtml::getGershgorinBoundsGenTask<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::getGershgorinBoundsSymmTask<LeafMatrixType >));
CHTML_REGISTER_TASK_TYPE((chtml::GetEigMinMaxPart<real>));
CHTML_REGISTER_TASK_TYPE((chtml::CompareEigMinMaxParts<real>));
CHTML_REGISTER_TASK_TYPE((chtml::GetEigMinMaxPartZeroDiag<real>));


CHTML_REGISTER_TASK_TYPE((chtml::RegisterMatrixByBlocks<LeafMatrixType>));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixNNZDiagLeafs<LeafMatrixType>));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixSqFrobNormDiffGen<LeafMatrixType>));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixSquaredFrobOfErrorMatrixSymm<LeafMatrixType>));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixFrobTruncLowestLevelSymm<LeafMatrixType>));
CHTML_REGISTER_TASK_TYPE((chtml::GeneralToSymmetric<LeafMatrixType>));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixSymmMultiplyUpperTr<LeafMatrixType, false, true>));
CHTML_REGISTER_TASK_TYPE((chtml::MatrixSymmMultiplyUpperTr<LeafMatrixType, true, false>));


#endif
