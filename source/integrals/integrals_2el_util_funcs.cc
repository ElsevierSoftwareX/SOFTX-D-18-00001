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

/** @file integrals_2el_util_funcs.cc

    \brief Code for utility functions used by 2-electron integral
    computation (i.e. computation of J and K matrices).

    @author: Elias Rudberg <em>responsible</em>.
*/

#include "integrals_2el_util_funcs.h"

/* ELIAS NOTE 2014-05-29: The do_summedIntegralList_contribs_std
   routine defined in this file is responsible for a large part of the
   computational effort for both J and K matrix
   construction. Therefore, this routine would be a good candidate for
   further optimization attempts.  */


/* This is the simple implementation, without unrolling.  It turned
   out that this could be optimized significantly by unrolling the
   outer loop.  */
/*
ELIAS NOTE 2014-07-13: Commented out this unused routine to silence compiler warning.
static void do_summedIntegralList_contribs_std_simple(const i_j_val_struct* conv_mat_1_sp, int conv_mat_1_sp_nnz,
						      const i_j_val_struct* conv_mat_2_sp, int conv_mat_2_sp_nnz,
						      int noOfMonomials_1, int noOfMonomials_2,
						      const ergo_real* primitiveIntegralList,
						      int noOfBasisFuncPairs_1, int noOfBasisFuncPairs_2,
						      ergo_real* summedIntegralList) {
  for(int idx_i = 0; idx_i < conv_mat_1_sp_nnz; idx_i++) {
    int idx_1 = conv_mat_1_sp[idx_i].i;
    int ii = conv_mat_1_sp[idx_i].j;
    ergo_real value_i = conv_mat_1_sp[idx_i].value;
    const ergo_real* primitiveIntegralListPtr = &primitiveIntegralList[ii*noOfMonomials_2];
    ergo_real* summedIntegralListPtr = &summedIntegralList[idx_1*noOfBasisFuncPairs_2];
    int idx_j = 0;
    while(idx_j < conv_mat_2_sp_nnz) {
      int nn = conv_mat_2_sp[idx_j].same_i_count;
      ergo_real sum = 0;
      for(int kk = 0; kk < nn; kk++) {
	int jj = conv_mat_2_sp[idx_j+kk].j;
	sum += primitiveIntegralListPtr[jj] * conv_mat_2_sp[idx_j+kk].value;
      }
      int idx_2 = conv_mat_2_sp[idx_j].i;
      summedIntegralListPtr[idx_2] += sum * value_i;
      idx_j += nn;
    }
  }
}
*/


void do_summedIntegralList_contribs_std(const i_j_val_struct* conv_mat_1_sp, int conv_mat_1_sp_nnz,
					const i_j_val_struct* conv_mat_2_sp, int conv_mat_2_sp_nnz,
					int noOfMonomials_1, int noOfMonomials_2,
					const ergo_real* primitiveIntegralList,
					int noOfBasisFuncPairs_1, int noOfBasisFuncPairs_2,
					ergo_real* summedIntegralList) {
  int idx_i = 0;
  while(idx_i < conv_mat_1_sp_nnz) {
    if(idx_i == conv_mat_1_sp_nnz-1) {
      // Only one left; treat in the old way
      int idx_1 = conv_mat_1_sp[idx_i].i;
      int ii = conv_mat_1_sp[idx_i].j;
      ergo_real value_i = conv_mat_1_sp[idx_i].value;
      const ergo_real* primitiveIntegralListPtr = &primitiveIntegralList[ii*noOfMonomials_2];
      ergo_real* summedIntegralListPtr = &summedIntegralList[idx_1*noOfBasisFuncPairs_2];
      int idx_j = 0;
      while(idx_j < conv_mat_2_sp_nnz) {
	int nn = conv_mat_2_sp[idx_j].same_i_count;
	ergo_real sum = 0;
	for(int kk = 0; kk < nn; kk++) {
	  int jj = conv_mat_2_sp[idx_j+kk].j;
	  sum += primitiveIntegralListPtr[jj] * conv_mat_2_sp[idx_j+kk].value;
	}
	int idx_2 = conv_mat_2_sp[idx_j].i;
	summedIntegralListPtr[idx_2] += sum * value_i;
	idx_j += nn;
      }
      idx_i++;
    }
    else if(idx_i == conv_mat_1_sp_nnz-2) {
      // Unroll by 2
      int idx_1_A = conv_mat_1_sp[idx_i].i;
      int ii_A = conv_mat_1_sp[idx_i].j;
      ergo_real value_i_A = conv_mat_1_sp[idx_i].value;
      const ergo_real* primitiveIntegralListPtr_A = &primitiveIntegralList[ii_A*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_A = &summedIntegralList[idx_1_A*noOfBasisFuncPairs_2];
      int idx_1_B = conv_mat_1_sp[idx_i+1].i;
      int ii_B = conv_mat_1_sp[idx_i+1].j;
      ergo_real value_i_B = conv_mat_1_sp[idx_i+1].value;
      const ergo_real* primitiveIntegralListPtr_B = &primitiveIntegralList[ii_B*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_B = &summedIntegralList[idx_1_B*noOfBasisFuncPairs_2];
      int idx_j = 0;
      while(idx_j < conv_mat_2_sp_nnz) {
	int nn = conv_mat_2_sp[idx_j].same_i_count;
	ergo_real sum_A = 0;
	ergo_real sum_B = 0;
	for(int kk = 0; kk < nn; kk++) {
	  int jj = conv_mat_2_sp[idx_j+kk].j;
	  sum_A += primitiveIntegralListPtr_A[jj] * conv_mat_2_sp[idx_j+kk].value;
	  sum_B += primitiveIntegralListPtr_B[jj] * conv_mat_2_sp[idx_j+kk].value;
	}
	int idx_2 = conv_mat_2_sp[idx_j].i;
	summedIntegralListPtr_A[idx_2] += sum_A * value_i_A;
	summedIntegralListPtr_B[idx_2] += sum_B * value_i_B;
	idx_j += nn;
      }
      idx_i += 2;
    }
    else if(idx_i == conv_mat_1_sp_nnz-3) {
      // Unroll by 3
      int idx_1_A = conv_mat_1_sp[idx_i].i;
      int ii_A = conv_mat_1_sp[idx_i].j;
      ergo_real value_i_A = conv_mat_1_sp[idx_i].value;
      const ergo_real* primitiveIntegralListPtr_A = &primitiveIntegralList[ii_A*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_A = &summedIntegralList[idx_1_A*noOfBasisFuncPairs_2];
      int idx_1_B = conv_mat_1_sp[idx_i+1].i;
      int ii_B = conv_mat_1_sp[idx_i+1].j;
      ergo_real value_i_B = conv_mat_1_sp[idx_i+1].value;
      const ergo_real* primitiveIntegralListPtr_B = &primitiveIntegralList[ii_B*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_B = &summedIntegralList[idx_1_B*noOfBasisFuncPairs_2];
      int idx_1_C = conv_mat_1_sp[idx_i+2].i;
      int ii_C = conv_mat_1_sp[idx_i+2].j;
      ergo_real value_i_C = conv_mat_1_sp[idx_i+2].value;
      const ergo_real* primitiveIntegralListPtr_C = &primitiveIntegralList[ii_C*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_C = &summedIntegralList[idx_1_C*noOfBasisFuncPairs_2];
      int idx_j = 0;
      while(idx_j < conv_mat_2_sp_nnz) {
	int nn = conv_mat_2_sp[idx_j].same_i_count;
	ergo_real sum_A = 0;
	ergo_real sum_B = 0;
	ergo_real sum_C = 0;
	for(int kk = 0; kk < nn; kk++) {
	  int jj = conv_mat_2_sp[idx_j+kk].j;
	  sum_A += primitiveIntegralListPtr_A[jj] * conv_mat_2_sp[idx_j+kk].value;
	  sum_B += primitiveIntegralListPtr_B[jj] * conv_mat_2_sp[idx_j+kk].value;
	  sum_C += primitiveIntegralListPtr_C[jj] * conv_mat_2_sp[idx_j+kk].value;
	}
	int idx_2 = conv_mat_2_sp[idx_j].i;
	summedIntegralListPtr_A[idx_2] += sum_A * value_i_A;
	summedIntegralListPtr_B[idx_2] += sum_B * value_i_B;
	summedIntegralListPtr_C[idx_2] += sum_C * value_i_C;
	idx_j += nn;
      }
      idx_i += 3;
    }
    else if(idx_i == conv_mat_1_sp_nnz-4) {
      // Unroll by 4
      int idx_1_A = conv_mat_1_sp[idx_i].i;
      int ii_A = conv_mat_1_sp[idx_i].j;
      ergo_real value_i_A = conv_mat_1_sp[idx_i].value;
      const ergo_real* primitiveIntegralListPtr_A = &primitiveIntegralList[ii_A*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_A = &summedIntegralList[idx_1_A*noOfBasisFuncPairs_2];
      int idx_1_B = conv_mat_1_sp[idx_i+1].i;
      int ii_B = conv_mat_1_sp[idx_i+1].j;
      ergo_real value_i_B = conv_mat_1_sp[idx_i+1].value;
      const ergo_real* primitiveIntegralListPtr_B = &primitiveIntegralList[ii_B*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_B = &summedIntegralList[idx_1_B*noOfBasisFuncPairs_2];
      int idx_1_C = conv_mat_1_sp[idx_i+2].i;
      int ii_C = conv_mat_1_sp[idx_i+2].j;
      ergo_real value_i_C = conv_mat_1_sp[idx_i+2].value;
      const ergo_real* primitiveIntegralListPtr_C = &primitiveIntegralList[ii_C*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_C = &summedIntegralList[idx_1_C*noOfBasisFuncPairs_2];
      int idx_1_D = conv_mat_1_sp[idx_i+3].i;
      int ii_D = conv_mat_1_sp[idx_i+3].j;
      ergo_real value_i_D = conv_mat_1_sp[idx_i+3].value;
      const ergo_real* primitiveIntegralListPtr_D = &primitiveIntegralList[ii_D*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_D = &summedIntegralList[idx_1_D*noOfBasisFuncPairs_2];
      int idx_j = 0;
      while(idx_j < conv_mat_2_sp_nnz) {
	int nn = conv_mat_2_sp[idx_j].same_i_count;
	ergo_real sum_A = 0;
	ergo_real sum_B = 0;
	ergo_real sum_C = 0;
	ergo_real sum_D = 0;
	for(int kk = 0; kk < nn; kk++) {
	  int jj = conv_mat_2_sp[idx_j+kk].j;
	  sum_A += primitiveIntegralListPtr_A[jj] * conv_mat_2_sp[idx_j+kk].value;
	  sum_B += primitiveIntegralListPtr_B[jj] * conv_mat_2_sp[idx_j+kk].value;
	  sum_C += primitiveIntegralListPtr_C[jj] * conv_mat_2_sp[idx_j+kk].value;
	  sum_D += primitiveIntegralListPtr_D[jj] * conv_mat_2_sp[idx_j+kk].value;
	}
	int idx_2 = conv_mat_2_sp[idx_j].i;
	summedIntegralListPtr_A[idx_2] += sum_A * value_i_A;
	summedIntegralListPtr_B[idx_2] += sum_B * value_i_B;
	summedIntegralListPtr_C[idx_2] += sum_C * value_i_C;
	summedIntegralListPtr_D[idx_2] += sum_D * value_i_D;
	idx_j += nn;
      }
      idx_i += 4;
    }
    else if(idx_i == conv_mat_1_sp_nnz-5) {
      // Unroll by 5
      int idx_1_A = conv_mat_1_sp[idx_i].i;
      int ii_A = conv_mat_1_sp[idx_i].j;
      ergo_real value_i_A = conv_mat_1_sp[idx_i].value;
      const ergo_real* primitiveIntegralListPtr_A = &primitiveIntegralList[ii_A*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_A = &summedIntegralList[idx_1_A*noOfBasisFuncPairs_2];
      int idx_1_B = conv_mat_1_sp[idx_i+1].i;
      int ii_B = conv_mat_1_sp[idx_i+1].j;
      ergo_real value_i_B = conv_mat_1_sp[idx_i+1].value;
      const ergo_real* primitiveIntegralListPtr_B = &primitiveIntegralList[ii_B*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_B = &summedIntegralList[idx_1_B*noOfBasisFuncPairs_2];
      int idx_1_C = conv_mat_1_sp[idx_i+2].i;
      int ii_C = conv_mat_1_sp[idx_i+2].j;
      ergo_real value_i_C = conv_mat_1_sp[idx_i+2].value;
      const ergo_real* primitiveIntegralListPtr_C = &primitiveIntegralList[ii_C*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_C = &summedIntegralList[idx_1_C*noOfBasisFuncPairs_2];
      int idx_1_D = conv_mat_1_sp[idx_i+3].i;
      int ii_D = conv_mat_1_sp[idx_i+3].j;
      ergo_real value_i_D = conv_mat_1_sp[idx_i+3].value;
      const ergo_real* primitiveIntegralListPtr_D = &primitiveIntegralList[ii_D*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_D = &summedIntegralList[idx_1_D*noOfBasisFuncPairs_2];
      int idx_1_E = conv_mat_1_sp[idx_i+4].i;
      int ii_E = conv_mat_1_sp[idx_i+4].j;
      ergo_real value_i_E = conv_mat_1_sp[idx_i+4].value;
      const ergo_real* primitiveIntegralListPtr_E = &primitiveIntegralList[ii_E*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_E = &summedIntegralList[idx_1_E*noOfBasisFuncPairs_2];
      int idx_j = 0;
      while(idx_j < conv_mat_2_sp_nnz) {
	int nn = conv_mat_2_sp[idx_j].same_i_count;
	ergo_real sum_A = 0;
	ergo_real sum_B = 0;
	ergo_real sum_C = 0;
	ergo_real sum_D = 0;
	ergo_real sum_E = 0;
	for(int kk = 0; kk < nn; kk++) {
	  int jj = conv_mat_2_sp[idx_j+kk].j;
	  sum_A += primitiveIntegralListPtr_A[jj] * conv_mat_2_sp[idx_j+kk].value;
	  sum_B += primitiveIntegralListPtr_B[jj] * conv_mat_2_sp[idx_j+kk].value;
	  sum_C += primitiveIntegralListPtr_C[jj] * conv_mat_2_sp[idx_j+kk].value;
	  sum_D += primitiveIntegralListPtr_D[jj] * conv_mat_2_sp[idx_j+kk].value;
	  sum_E += primitiveIntegralListPtr_E[jj] * conv_mat_2_sp[idx_j+kk].value;
	}
	int idx_2 = conv_mat_2_sp[idx_j].i;
	summedIntegralListPtr_A[idx_2] += sum_A * value_i_A;
	summedIntegralListPtr_B[idx_2] += sum_B * value_i_B;
	summedIntegralListPtr_C[idx_2] += sum_C * value_i_C;
	summedIntegralListPtr_D[idx_2] += sum_D * value_i_D;
	summedIntegralListPtr_E[idx_2] += sum_E * value_i_E;
	idx_j += nn;
      }
      idx_i += 5;
    }
    else {
      // Unroll by 6
      int idx_1_A = conv_mat_1_sp[idx_i].i;
      int ii_A = conv_mat_1_sp[idx_i].j;
      ergo_real value_i_A = conv_mat_1_sp[idx_i].value;
      const ergo_real* primitiveIntegralListPtr_A = &primitiveIntegralList[ii_A*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_A = &summedIntegralList[idx_1_A*noOfBasisFuncPairs_2];
      int idx_1_B = conv_mat_1_sp[idx_i+1].i;
      int ii_B = conv_mat_1_sp[idx_i+1].j;
      ergo_real value_i_B = conv_mat_1_sp[idx_i+1].value;
      const ergo_real* primitiveIntegralListPtr_B = &primitiveIntegralList[ii_B*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_B = &summedIntegralList[idx_1_B*noOfBasisFuncPairs_2];
      int idx_1_C = conv_mat_1_sp[idx_i+2].i;
      int ii_C = conv_mat_1_sp[idx_i+2].j;
      ergo_real value_i_C = conv_mat_1_sp[idx_i+2].value;
      const ergo_real* primitiveIntegralListPtr_C = &primitiveIntegralList[ii_C*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_C = &summedIntegralList[idx_1_C*noOfBasisFuncPairs_2];
      int idx_1_D = conv_mat_1_sp[idx_i+3].i;
      int ii_D = conv_mat_1_sp[idx_i+3].j;
      ergo_real value_i_D = conv_mat_1_sp[idx_i+3].value;
      const ergo_real* primitiveIntegralListPtr_D = &primitiveIntegralList[ii_D*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_D = &summedIntegralList[idx_1_D*noOfBasisFuncPairs_2];
      int idx_1_E = conv_mat_1_sp[idx_i+4].i;
      int ii_E = conv_mat_1_sp[idx_i+4].j;
      ergo_real value_i_E = conv_mat_1_sp[idx_i+4].value;
      const ergo_real* primitiveIntegralListPtr_E = &primitiveIntegralList[ii_E*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_E = &summedIntegralList[idx_1_E*noOfBasisFuncPairs_2];
      int idx_1_F = conv_mat_1_sp[idx_i+5].i;
      int ii_F = conv_mat_1_sp[idx_i+5].j;
      ergo_real value_i_F = conv_mat_1_sp[idx_i+5].value;
      const ergo_real* primitiveIntegralListPtr_F = &primitiveIntegralList[ii_F*noOfMonomials_2];
      ergo_real* summedIntegralListPtr_F = &summedIntegralList[idx_1_F*noOfBasisFuncPairs_2];
      int idx_j = 0;
      while(idx_j < conv_mat_2_sp_nnz) {
	int nn = conv_mat_2_sp[idx_j].same_i_count;
	ergo_real sum_A = 0;
	ergo_real sum_B = 0;
	ergo_real sum_C = 0;
	ergo_real sum_D = 0;
	ergo_real sum_E = 0;
	ergo_real sum_F = 0;
	for(int kk = 0; kk < nn; kk++) {
	  int jj = conv_mat_2_sp[idx_j+kk].j;
	  sum_A += primitiveIntegralListPtr_A[jj] * conv_mat_2_sp[idx_j+kk].value;
	  sum_B += primitiveIntegralListPtr_B[jj] * conv_mat_2_sp[idx_j+kk].value;
	  sum_C += primitiveIntegralListPtr_C[jj] * conv_mat_2_sp[idx_j+kk].value;
	  sum_D += primitiveIntegralListPtr_D[jj] * conv_mat_2_sp[idx_j+kk].value;
	  sum_E += primitiveIntegralListPtr_E[jj] * conv_mat_2_sp[idx_j+kk].value;
	  sum_F += primitiveIntegralListPtr_F[jj] * conv_mat_2_sp[idx_j+kk].value;
	}
	int idx_2 = conv_mat_2_sp[idx_j].i;
	summedIntegralListPtr_A[idx_2] += sum_A * value_i_A;
	summedIntegralListPtr_B[idx_2] += sum_B * value_i_B;
	summedIntegralListPtr_C[idx_2] += sum_C * value_i_C;
	summedIntegralListPtr_D[idx_2] += sum_D * value_i_D;
	summedIntegralListPtr_E[idx_2] += sum_E * value_i_E;
	summedIntegralListPtr_F[idx_2] += sum_F * value_i_F;
	idx_j += nn;
      }
      idx_i += 6;
    }
  }
}


void do_summedIntegralList_contribs_self(const i_j_val_struct* conv_mat_1_sp, int conv_mat_1_sp_nnz,
					 const i_j_val_struct* conv_mat_2_sp, int conv_mat_2_sp_nnz,
					 int noOfMonomials_1, int noOfMonomials_2,
					 const ergo_real* primitiveIntegralList,
					 int noOfBasisFuncPairs_1, int noOfBasisFuncPairs_2,
					 ergo_real* summedIntegralList) {
  // Special interactionWithSelf case
  for(int idx_i = 0; idx_i < conv_mat_1_sp_nnz; idx_i++) {
    int idx_1 = conv_mat_1_sp[idx_i].i;
    int ii = conv_mat_1_sp[idx_i].j;
    ergo_real value_i = conv_mat_1_sp[idx_i].value;
    const ergo_real* primitiveIntegralListPtr = &primitiveIntegralList[ii*noOfMonomials_2];
    ergo_real* summedIntegralListPtr = &summedIntegralList[idx_1*noOfBasisFuncPairs_2];
    int idx_j = 0;
    while(idx_j < conv_mat_2_sp_nnz) {
      int nn = conv_mat_2_sp[idx_j].same_i_count;
      ergo_real sum = 0;
      for(int kk = 0; kk < nn; kk++) {
	int jj = conv_mat_2_sp[idx_j+kk].j;
	sum += primitiveIntegralListPtr[jj] * conv_mat_2_sp[idx_j+kk].value;
      }
      int idx_2 = conv_mat_2_sp[idx_j].i;
      if(idx_1 == idx_2)
	summedIntegralListPtr[idx_2] += sum * value_i * 0.5;
      else if(idx_1 > idx_2)
	summedIntegralListPtr[idx_2] += sum * value_i;
      idx_j += nn;
    }
  }
}
