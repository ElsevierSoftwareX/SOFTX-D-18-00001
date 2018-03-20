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

/** @file gemm_sse.h
 *
 * @brief Generalized matrix matrix multiplication using SSE intrinsics.
 *
 * @author Emanuel H. Rubensson
 * @date 2009
 *
 */


#ifndef GEMM_SSE_H
#define GEMM_SSE_H
#include <stdexcept>
#include "mm_kernel_inner_sse2_A.h"
#include "mm_kernel_outer_A.h"


template<typename real, typename regType, 
  int m_kernel, int n_kernel, int k_kernel,
  int m_block, int n_block>
  static void gemm_sse(real const * const A,
		       real const * const B,
		       real *             C,
		       size_t const m,
		       size_t const n,
		       size_t const k,
		       real * A_packed,
		       real * B_packed,
		       real * C_packed,
		       size_t const ap_size,
		       size_t const bp_size,
		       size_t const cp_size) {
  // typedef double real; typedef __m128d regType;
  // typedef float real; typedef __m128 regType;
  typedef MM_kernel_inner_sse2_A<real, regType, m_kernel, n_kernel, k_kernel> MM_inner;
  typedef MM_kernel_outer_A<MM_inner, m_block, n_block> MM_outer;
  if (m != m_kernel*m_block) 
    throw std::runtime_error("Error in gemm_sse(...): m != m_kernel*m_block");
  if (n != n_kernel*n_block) 
    throw std::runtime_error("Error in gemm_sse(...): n != n_kernel*n_block");
  if (k != k_kernel) 
    throw std::runtime_error("Error in gemm_sse(...): k != k_kernel");
  if (ap_size < MM_outer::Pack_type_A::size_packed) 
    throw std::runtime_error("Error in gemm_sse(...): "
			     "ap_size < MM_outer::Pack_type_A::size_packed");
  if (bp_size < MM_outer::Pack_type_B::size_packed) 
    throw std::runtime_error("Error in gemm_sse(...): "
			     "bp_size < MM_outer::Pack_type_B::size_packed");
  if (cp_size < MM_outer::Pack_type_C::size_packed) 
    throw std::runtime_error("Error in gemm_sse(...): "
			     "cp_size < MM_outer::Pack_type_C::size_packed");
  MM_outer::Pack_type_C::template pack<Ordering_col_wise>( C, C_packed, m, n);
  MM_outer::Pack_type_A::template pack<Ordering_col_wise>( A, A_packed, m, k);
  MM_outer::Pack_type_B::template pack<Ordering_col_wise>( B, B_packed, k, n);
  MM_outer::exec(&A_packed, &B_packed, C_packed);
  MM_outer::Pack_type_C::template unpack<Ordering_col_wise>(C, C_packed, m, n);
}

template<typename real>
static void gemm_sse(real const * const A,
		     real const * const B,
		     real *             C,
		     size_t const m,
		     size_t const n,
		     size_t const k,
		     real * A_packed,
		     real * B_packed,
		     real * C_packed,
		     size_t const ap_size,
		     size_t const bp_size,
		     size_t const cp_size) {
  throw std::runtime_error("gemm_sse not implemented for chosen real type.");
}

template<>
void gemm_sse(double const * const A,
	      double const * const B,
	      double *             C,
	      size_t const m,
	      size_t const n,
	      size_t const k,
	      double * A_packed,
	      double * B_packed,
	      double * C_packed,
	      size_t const ap_size,
	      size_t const bp_size,
	      size_t const cp_size) {
  gemm_sse<double, __m128d, 4, 4, 32, 8, 8>
    (A, B, C, m, n, k, 
     A_packed, B_packed, C_packed, ap_size, bp_size, cp_size);
}

template<>
void gemm_sse(float const * const A,
	      float const * const B,
	      float *             C,
	      size_t const m,
	      size_t const n,
	      size_t const k,
	      float * A_packed,
	      float * B_packed,
	      float * C_packed,
	      size_t const ap_size,
	      size_t const bp_size,
	      size_t const cp_size) {
  gemm_sse<float, __m128, 8, 4, 32, 4, 8>
    (A, B, C, m, n, k, 
     A_packed, B_packed, C_packed, ap_size, bp_size, cp_size);
}

#endif
