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

/** @file mm_kernel_outer_A.h
 *
 * @brief Templates for efficient gemm kernels. Outer kernel that
 * wraps around given inner kernel.
 *
 * @author Emanuel H. Rubensson
 * @date 2009
 *
 */

#ifndef MM_KERNEL_OUTER_A_H
#define MM_KERNEL_OUTER_A_H
#include "common.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/** Template for matrix matrix multiplication that wraps around a kernel given as template argument. 
 *
 *  The idea is that the inner kernel should be fully unrolled and block for registers.
 *  
 */
template<typename T_gemm_kernel, int T_M_block, int T_N_block>
  class MM_kernel_outer_A {
  template<int T_rows_block, int T_cols_block, typename T_ordering_block, typename T_pack_type_kernel>
  class Pack;
 public:
  static int const M_kernel = T_gemm_kernel::M;        /**< Number of rows of A and C kernels.            */
  static int const N_kernel = T_gemm_kernel::N;        /**< Number of columns of B and C kernels.         */
  static int const K_kernel = T_gemm_kernel::K;        /**< Number of columns of A kernels and rows of B kernels. */
  static int const M_block  = T_M_block;               /**< Number of rows of A and C (blocks).            */
  static int const N_block  = T_N_block;               /**< Number of columns of B and C (blocks).         */
  static int const K_block  = 1;                       /**< Number of columns of A and rows of B (blocks). */
  static int const M = M_kernel * M_block;             /**< Number of rows of A and C.            */
  static int const N = N_kernel * N_block;             /**< Number of columns of B and C.         */
  static int const K = K_kernel * K_block;             /**< Number of columns of A and rows of B. */
  typedef typename T_gemm_kernel::real real;           /**< Real number type (usually float or double) */
  
  typedef Ordering_col_wise Ordering_block_A;
  typedef Ordering_col_wise Ordering_block_B;
  typedef Ordering_col_wise Ordering_block_C;
  
  typedef Pack< M_block, K_block, Ordering_block_A, typename T_gemm_kernel::Pack_type_A > Pack_type_A;
  typedef Pack< K_block, N_block, Ordering_block_B, typename T_gemm_kernel::Pack_type_B > Pack_type_B;
  typedef Pack< M_block, N_block, Ordering_block_C, typename T_gemm_kernel::Pack_type_C > Pack_type_C;
  /** Executes the matrix-matrix multiply C += A B with the three matrices A, B, and C 
   *  stored using the packing types of this class. 
   */
  static void exec( real const * const * const A, 
		    real const * const * const B, 
		    real * const C,
		    int const i = 1);
  
};

template<typename T_gemm_kernel, int T_M_block, int T_N_block>
  void MM_kernel_outer_A<T_gemm_kernel, T_M_block, T_N_block>::exec( real const * const * const A, 
								     real const * const * const B, 
								     real * const C,
								     int const n_mul ) {
#if 1
  for ( int n = 0; n < N_block; ++n )
    for ( int m = 0; m < M_block; ++m ) {
      T_gemm_kernel::exec( A, B, C, n_mul,  
			   Ordering_block_A::get( m, 0, M_block, K_block ) * T_gemm_kernel::Pack_type_A::size_packed,
			   Ordering_block_B::get( 0, n, K_block, N_block ) * T_gemm_kernel::Pack_type_B::size_packed,
			   Ordering_block_C::get( m, n, M_block, N_block ) * T_gemm_kernel::Pack_type_C::size_packed );      
    }
  
#else
#if 1
  // FIXME: This is faster since the offsets are known at compile time, TODO: unroll for loops...
  T_gemm_kernel::template exec<Ordering_block_A::template Get<0, 0, M_block, K_block>::index * T_gemm_kernel::Pack_type_A::size_packed,
    Ordering_block_B::template Get<0, 0, K_block, N_block>::index * T_gemm_kernel::Pack_type_B::size_packed,
    Ordering_block_C::template Get<0, 0, M_block, N_block>::index * T_gemm_kernel::Pack_type_C::size_packed>( A, B, C, n_mul );      
  T_gemm_kernel::template exec<Ordering_block_A::template Get<1, 0, M_block, K_block>::index * T_gemm_kernel::Pack_type_A::size_packed,
    Ordering_block_B::template Get<0, 0, K_block, N_block>::index * T_gemm_kernel::Pack_type_B::size_packed,
    Ordering_block_C::template Get<1, 0, M_block, N_block>::index * T_gemm_kernel::Pack_type_C::size_packed>( A, B, C, n_mul );      
#else
  T_gemm_kernel::exec( A, B, C, n_mul,  
		       Ordering_block_A::get( 0, 0, M_block, K_block ) * T_gemm_kernel::Pack_type_A::size_packed,
		       Ordering_block_B::get( 0, 0, K_block, N_block ) * T_gemm_kernel::Pack_type_B::size_packed,
		       Ordering_block_C::get( 0, 0, M_block, N_block ) * T_gemm_kernel::Pack_type_C::size_packed );      
#endif
#endif
}


/** Template for for translations between unpacked and packed matrix storage. 
 *
 *  Template arguments:
 *    - T_rows_block       : number of rows (blocks)
 *    - T_cols_block       : number of columns (blocks)
 *    - T_ordering_block   : Type that specifies how the matrix blocks are stored
 *    - T_pack_type_kernel : Type specifying how each matrix block should be packed 
 *
 */
template<typename T_gemm_kernel, int T_M_block, int T_N_block>
  template<int T_rows_block, int T_cols_block, typename T_ordering_block, typename T_pack_type_kernel>
  class MM_kernel_outer_A<T_gemm_kernel, T_M_block, T_N_block>::Pack {
  static int const rows_kernel    = T_pack_type_kernel::rows;
  static int const cols_kernel    = T_pack_type_kernel::cols;
 public:
  static int const rows     = rows_kernel * T_rows_block;  /**< Number of rows in the matrix. */
  static int const cols     = cols_kernel * T_cols_block;  /**< Number of columns in the matrix. */
  /** Memory needed to store the matrix in packed form. 
   *  (Can be used in the allocation of memory for the packed matrix.) 
   */
  //  static int const size_packed    = rows * cols * T_pack_type_kernel::size_packed;
  static unsigned int const size_packed    = T_rows_block * T_cols_block * T_pack_type_kernel::size_packed;
  //  typedef Packed<T_real, T_rows_block, T_cols_block, T_rows_kernel, T_cols_kernel, T_kernel_index, T_repetitions> ThisType;
  
  template<typename T_ordering_matrix> 
  struct Assign_to_packed : public T_pack_type_kernel::template Assign_to_packed<T_ordering_matrix> {
    typedef T_ordering_matrix Ordering_matrix;
  };
  template<typename T_ordering_matrix> 
  struct Extract_from_packed : public T_pack_type_kernel::template Extract_from_packed<T_ordering_matrix> {
    typedef T_ordering_matrix Ordering_matrix;
  };
  
  
  /** Elaborate function that can be called either to assign to or extract from packed format. 
   *
   *  Use of types in T_assign automatically sets the const qualifier on the desired argument. 
   */
  template<template<typename T_ordering> class T_assign, typename T_ordering_matrix>
    static void exec(typename T_assign<T_ordering_matrix>::PtrType X, typename T_assign<T_ordering_matrix>::PtrTypePacked X_packed,
		     int const rows_total_matrix, int const cols_total_matrix) {
    // Loop over column blocks of new packed matrix
    for ( int col_b = 0; col_b < T_cols_block; ++col_b ) {
      // Loop over row blocks of new packed matrix
      for ( int row_b = 0; row_b < T_rows_block; ++row_b ) {
	T_pack_type_kernel::template exec< T_assign, T_ordering_matrix >
	  ( &X[ T_assign<T_ordering_matrix>::Ordering_matrix::get( row_b * rows_kernel, col_b * cols_kernel, 
								   rows_total_matrix, cols_total_matrix ) ], 
	    &X_packed[ T_ordering_block::get( row_b, col_b, T_rows_block, T_cols_block ) *
		       T_pack_type_kernel::size_packed ],
	    rows_total_matrix, cols_total_matrix );
	// Indexes of original matrix             : ( row_b * rows_kernel, col_b * cols_kernel )
	// Block indexes (packed matrix)          : ( row_b, col_b )
	// Number of reals needed for each kernel : T_pack_type_kernel::size_packed
      }
    }
  } // end exec()
  
  /** Convenience function for assignments to packed matrix.
   *  The template argument specifies how the original (unpacked) matrix
   *  is stored (e.g. column or row wise)
   */
  template<typename T_ordering_matrix>
    inline static void pack(real const * const X, real * X_packed,
			    int const rows_total_matrix, int const cols_total_matrix) {
    exec< Assign_to_packed, T_ordering_matrix >(X, X_packed, rows_total_matrix, cols_total_matrix);
  }
  /** Convenience function for extracting matrix from packed matrix.
   *  The template argument specifies how the unpacked matrix
   *  is stored (e.g. column or row wise)
   */
  template<typename T_ordering_matrix>
    inline static void unpack(real * X, real const * const X_packed,
			      int const rows_total_matrix, int const cols_total_matrix) {
    exec< Extract_from_packed, T_ordering_matrix >(X, X_packed, rows_total_matrix, cols_total_matrix);
  }
  
  //  real * values;
};
#endif
