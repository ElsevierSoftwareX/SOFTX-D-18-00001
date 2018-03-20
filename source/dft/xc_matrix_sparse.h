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

/** @file xc_matrix_sparse.h

    @brief The sparse XC matrix evaluator.

    @author: Pawel Salek <em>responsible</em>
*/

#if !defined(_XC_MATRIX_SPARSE_H_)
#define _XC_MATRIX_SPARSE_H_ 1

#define BEGIN_NAMESPACE(x) namespace x {
#define END_NAMESPACE(x)   } /* x */

#include "basisinfo.h"
#include "matrix_typedefs.h"
#include "realtype.h"


typedef ergo_real real;

BEGIN_NAMESPACE(Dft)

real getXC_seq(const BasisInfoStruct& bis, const IntegralInfo& integralInfo,
               const Molecule& mol,  const Dft::GridParams& gss,
               int nelectrons, const symmMatrix& dmat,
               symmMatrix& ksm, real* edfty, 
               std::vector<int> const & permutationHML);

real getXC_mt(const BasisInfoStruct& bis, const IntegralInfo& integralInfo,
              const Molecule& mol,  const Dft::GridParams& gss,
              int nElectrons, const symmMatrix& dens,
              symmMatrix& xcm, real* xcEnergy,
              std::vector<int> const & permutationHML);

real getUXC_seq(const BasisInfoStruct& bis, const IntegralInfo& integralInfo,
                const Molecule& mol, const Dft::GridParams& gss, int nElectrons,
                const symmMatrix& densA, const symmMatrix& densB,
                symmMatrix& xcA, symmMatrix& xcB, real* xcEnergy,
                std::vector<int> const & permutationHML);

real getUXC_mt(const BasisInfoStruct& bis, const IntegralInfo& integralInfo,
               const Molecule& mol, const Dft::GridParams& gss, int nElectrons,
               const symmMatrix& densA, const symmMatrix& densB,
               symmMatrix& xcA, symmMatrix& xcB, real* xcEnergy,
               std::vector<int> const & permutationHML);

END_NAMESPACE(Dft)

#endif /* _XC_MATRIX_SPARSE_H_ */
