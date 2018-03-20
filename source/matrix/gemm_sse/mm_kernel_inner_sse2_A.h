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

/** @file mm_kernel_inner_sse2_A.h
 *
 * @brief Templates for efficient gemm kernels.  For architectures
 * with SSE2 or higher.
 *
 * @author Emanuel H. Rubensson
 * @date 2009
 *
 */
#ifndef MM_KERNEL_INNER_SSE2_A_H
#define MM_KERNEL_INNER_SSE2_A_H
#include "common.h"
#include "vector_intrin.h"

/** Matrix multiplication template for architectures with SSE2 or higher 
 *  and compilers that support C++ intrinsics for access to SSE 
 *  instructions.  
 *
 *  Choice of template parameters:
 *   - T_M and T_N should be chosen so that the T_M x T_N matrix C  
 *     fits in registers. For example T_M == T_N == 4
 *   - T_K should be chosen so that the generated code fits in L1 instruction cache.
 *     For example T_K == 128.
 *   - T_real and T_reg must go together.
 *       Example:
 *       - <T_real, T_reg> == <double, __m128d> 
 *       - <T_real, T_reg> == <float, __m128> 
 *   
 *  The public typedefs and static members specify how the matrices must be 
 *  stored. 
 *
 */
template<typename T_real, typename T_reg, int T_M, int T_N, int T_K>
  class MM_kernel_inner_sse2_A {
 public:
  typedef T_real real;              /**< Real number type (usually float or double) */
  static int const M = T_M;         /**< Number of rows of A and C.                 */
  static int const N = T_N;         /**< Number of columns of B and C.              */
  static int const K = T_K;         /**< Number of columns of A and rows of B.      */
 protected: 
  static int const floats_per_register = ( sizeof(T_reg) / sizeof(real) );
  /**< Number of real numbers that fit in one register. */
  
 private:
  /** Template for packing of matrix elements. */
  template<int T_ROWS_kernel, int T_COLS_kernel, typename T_ordering_kernel, int T_repetitions>
    class Pack;
 public:
  typedef Pack< M, K, Ordering_col_wise, 1                   > Pack_type_A; /**< Type that can (should) be used to pack A. */
  typedef Pack< K, N, Ordering_row_wise, floats_per_register > Pack_type_B; /**< Type that can (should) be used to pack B. */
  typedef Pack< M, N, Ordering_col_wise, 1                   > Pack_type_C; /**< Type that can (should) be used to pack C. */
 
 
 
  // Consider passing derived class from std::vector as arguments
  // that have compile time constant length that can be checked at 
  // compile time using static asserts
  /** Executes the matrix-matrix multiply C += A B with the three matrices A, B, and C 
   *  stored according to the static members and typedefs of this class. 
   */
  static void exec( real const * const * const A, 
		    real const * const * const B, 
		    real * const C,
		    int const i = 1,
		    int const offset_A = 0,
		    int const offset_B = 0,
		    int const offset_C = 0 );

  template<int T_offset_A, int T_offset_B, int T_offset_C>
    static void exec( real const * const * const A, 
		      real const * const * const B, 
		      real * const C,
		      int const i = 1 );
  
 
 
 protected:
  template<int T_loop_index, int T_end>
    struct Loop {
      static inline void ALWAYS_INLINE set_to_zero( Vector_intrin<real, T_reg> * X_reg ) {
	X_reg[T_loop_index].set_to_zero();
	Loop<T_loop_index+1, T_end>::set_to_zero( X_reg );
      }
      static inline void ALWAYS_INLINE inner( int const row_A_reg, // == row_C_reg
					      int const row_B,     
					      Vector_intrin<real, T_reg> const & A_reg,
					      Vector_intrin<real, T_reg> * C_reg,
					      real  const * B_packed ) {
	Vector_intrin<real, T_reg> B_reg;
	B_reg.load_p( &B_packed[row_B * T_N * floats_per_register + 
				T_loop_index * floats_per_register] );
	B_reg *= A_reg;
	C_reg[row_A_reg + T_loop_index * T_M / floats_per_register] += B_reg;
	Loop<T_loop_index+1, T_end>::inner( row_A_reg, row_B, 
					    A_reg, C_reg, 
					    B_packed );     
      }
      static inline void ALWAYS_INLINE middle( int const col_A, // == row_B
					       Vector_intrin<real, T_reg> * C_reg,
					       real  const * A,
					       real  const * B_packed ) {
	Vector_intrin<real, T_reg> A_reg;
	A_reg.load_p( &A[col_A * T_M + T_loop_index * floats_per_register] );
	// Loop over cols of B
	Loop<0, T_N>::inner( T_loop_index, // == row_A_reg == row_C_reg
			     col_A,        // == row_B
			     A_reg,
			     C_reg,
			     B_packed );
	Loop<T_loop_index+1, T_end>::middle( col_A, C_reg, A, B_packed );     
      }
      static inline void ALWAYS_INLINE outer( int const  start_i,
					      Vector_intrin<real, T_reg> * C_reg,
					      real const * A,
					      real const * B_packed ) {  
	// Loop over (register) rows of A and C 
	Loop<0, T_M/floats_per_register>::middle( start_i + T_loop_index, 
						  C_reg, 
						  A,
						  B_packed );
	Loop<T_loop_index+1, T_end>::outer( start_i, C_reg, A, B_packed );     
      }
      static inline void ALWAYS_INLINE add( Vector_intrin<real, T_reg> * X_reg,
					    real const * X ) {
	X_reg[T_loop_index] += &X[T_loop_index * floats_per_register];
	Loop<T_loop_index+1, T_end>::add( X_reg, X );
      }   
      static inline void ALWAYS_INLINE store( Vector_intrin<real, T_reg> const * X_reg,
					      real * X ) {
	X_reg[T_loop_index].store_p( &X[T_loop_index * floats_per_register] );
	Loop<T_loop_index+1, T_end>::store( X_reg, X );
      }

      static inline void ALWAYS_INLINE multiple_loop( Vector_intrin<real, T_reg> * C_reg,
						      real const * const * const A, 
						      real const * const * const B ) {
	// Loop over columns of A and rows of B^T
	Loop<0, T_K>::outer( 0, C_reg, A[T_loop_index], B[T_loop_index] );
	Loop<T_loop_index+1, T_end>::multiple_loop( C_reg, A, B );
      }
    };
 
  template<int T_end>
    struct Loop<T_end, T_end> {
    static inline void ALWAYS_INLINE set_to_zero( Vector_intrin<real, T_reg>  * X_reg ) {}
    static inline void ALWAYS_INLINE inner( int const row_A_reg, // == row_C_reg
					    int const row_B,     
					    Vector_intrin<real, T_reg> const & A_reg,
					    Vector_intrin<real, T_reg> * C_reg,
					    real  const * B_packed ) {}
    static inline void ALWAYS_INLINE middle( int const col_A, // == row_B
					     Vector_intrin<real, T_reg> * C_reg,
					     real  const * A,
					     real  const * B_packed ) {}
    static inline void ALWAYS_INLINE outer( int const  start_i,
					    Vector_intrin<real, T_reg> * C_reg,
					    real const * A,
					    real const * B_packed ) {}
    static inline void ALWAYS_INLINE add( Vector_intrin<real, T_reg> * X_reg,
					  real const * X ) {}
    static inline void ALWAYS_INLINE store( Vector_intrin<real, T_reg> const * X_reg,
					    real * X ) {}
    static inline void ALWAYS_INLINE multiple_loop( Vector_intrin<real, T_reg> * C_reg,
						    real const * const * const A, 
						    real const * const * const B ) {}
  };
};

// Doesn't matter if inlined or not... same performance for 4, 4, 5
// IT DOES MATTER WHEN OUTER CONTROL STRUCTURE IS ADDED ON TOP!!!
// THEN, INLINING IS BAD 
template<typename real, typename T_reg, int T_M, int T_N, int T_K>
  void MM_kernel_inner_sse2_A<real, T_reg, T_M, T_N, T_K>::exec(real const * const * const A, 
								real const * const * const B, 
								real * const C,
								int const i,
								int const offset_A,
								int const offset_B,
								int const offset_C) {
  STATIC_ASSERT_DEBUG(!(T_M%floats_per_register), TEMPLATE_ARGUMENT_T_M_MUST_BE_MULTIPLE_OF_floats_per_register);
  Vector_intrin<real, T_reg> C_reg[T_M * T_N / floats_per_register];
  MM_kernel_inner_sse2_A<real, T_reg, T_M, T_N, T_K>::template Loop<0, T_M * T_N / floats_per_register>::set_to_zero( C_reg );
#if 1 // I loose a bit performance because of the offsets
  for (int ind = 0; ind < i; ++ind)
    MM_kernel_inner_sse2_A<real, T_reg, T_M, T_N, T_K>::template Loop<0, T_K>::outer( 0, C_reg, A[ind] + offset_A, B[ind] + offset_B );
  //// Loop over columns of A and rows of B^T
  //  MM_kernel_inner_sse2_A<real, T_reg, T_M, T_N, T_K>::template Loop<0, T_K>::outer( 0, C_reg, A, B );
  MM_kernel_inner_sse2_A<real, T_reg, T_M, T_N, T_K>::template Loop<0, T_M * T_N / floats_per_register>::add( C_reg, C + offset_C);
  MM_kernel_inner_sse2_A<real, T_reg, T_M, T_N, T_K>::template Loop<0, T_M * T_N / floats_per_register>::store( C_reg, C + offset_C);
#else
  for (int ind = 0; ind < i; ++ind)
    MM_kernel_inner_sse2_A<real, T_reg, T_M, T_N, T_K>::template Loop<0, T_K>::outer( 0, C_reg, A[ind], B[ind] );
  //// Loop over columns of A and rows of B^T
  //  MM_kernel_inner_sse2_A<real, T_reg, T_M, T_N, T_K>::template Loop<0, T_K>::outer( 0, C_reg, A, B );
  MM_kernel_inner_sse2_A<real, T_reg, T_M, T_N, T_K>::template Loop<0, T_M * T_N / floats_per_register>::add( C_reg, C);
  MM_kernel_inner_sse2_A<real, T_reg, T_M, T_N, T_K>::template Loop<0, T_M * T_N / floats_per_register>::store( C_reg, C);
#endif
} // end exec

template<typename real, typename T_reg, int T_M, int T_N, int T_K>
  template<int T_offset_A, int T_offset_B, int T_offset_C>
  void MM_kernel_inner_sse2_A<real, T_reg, T_M, T_N, T_K>::exec( real const * const * const A, 
								 real const * const * const B, 
								 real * const C,
								 int const i ) {
  STATIC_ASSERT_DEBUG(!(T_M%floats_per_register), TEMPLATE_ARGUMENT_T_M_MUST_BE_MULTIPLE_OF_floats_per_register);
  Vector_intrin<real, T_reg> C_reg[T_M * T_N / floats_per_register];
  MM_kernel_inner_sse2_A<real, T_reg, T_M, T_N, T_K>::template Loop<0, T_M * T_N / floats_per_register>::set_to_zero( C_reg );
  for (int ind = 0; ind < i; ++ind)
    MM_kernel_inner_sse2_A<real, T_reg, T_M, T_N, T_K>::template Loop<0, T_K>::outer( 0, C_reg, A[ind] + T_offset_A, B[ind] + T_offset_B );
  //// Loop over columns of A and rows of B^T
  //  MM_kernel_inner_sse2_A<real, T_reg, T_M, T_N, T_K>::template Loop<0, T_K>::outer( 0, C_reg, A, B );
  MM_kernel_inner_sse2_A<real, T_reg, T_M, T_N, T_K>::template Loop<0, T_M * T_N / floats_per_register>::add( C_reg, C + T_offset_C);
  MM_kernel_inner_sse2_A<real, T_reg, T_M, T_N, T_K>::template Loop<0, T_M * T_N / floats_per_register>::store( C_reg, C + T_offset_C);
} // end exec template






/** Class template for packing of matrix elements prior to matrix-matrix multiply. */
template<typename real, typename T_reg, int T_M, int T_N, int T_K>
  template<int T_rows, int T_cols, typename T_ordering_kernel, int T_repetitions>
  class MM_kernel_inner_sse2_A<real, T_reg, T_M, T_N, T_K>::Pack {
 public:
  static int const size_packed = T_rows * T_cols * T_repetitions;
  static int const rows = T_rows;
  static int const cols = T_cols;
  
  template<typename T_ordering_matrix>
    struct Assign_to_packed {
      typedef real * PtrTypePacked;        /**< Type of packed pointer - note the absence of const qualifiers. */
      typedef real const * const PtrType;  /**< Type of matrix pointer - note the presence of const qualifiers. */
      inline static void exec( PtrType X, PtrTypePacked X_packed,
			       int const row_k,
			       int const col_k,
			       int const rows_total_matrix,
			       int const cols_total_matrix ) {
	for ( int ir = 0; ir < T_repetitions; ++ir)
	  X_packed[ T_ordering_kernel::get( row_k, col_k, T_rows, T_cols ) * T_repetitions + ir ]
	    = X[ T_ordering_matrix::get(row_k, col_k, rows_total_matrix, cols_total_matrix) ];
      }      
    };
  
  template<typename T_ordering_matrix>
    struct Extract_from_packed {
      typedef real const * const PtrTypePacked; /**< Type of packed pointer - note the presence of const qualifiers. */
      typedef real * PtrType;                   /**< Type of matrix pointer - note the absence of const qualifiers. */
      inline static void exec( PtrType X, PtrTypePacked X_packed,
			       int const row_k,
			       int const col_k,
			       int const rows_total_matrix,
			       int const cols_total_matrix ) {
	for ( int ir = 0; ir < T_repetitions; ++ir)
	  X[ T_ordering_matrix::get(row_k, col_k, rows_total_matrix, cols_total_matrix) ] = 
	    X_packed[ T_ordering_kernel::get( row_k, col_k, T_rows, T_cols ) * T_repetitions + ir ];
      }      
    };
  
  template<template<typename T_ordering> class T_assign, typename T_ordering_matrix>
    static void exec(typename T_assign<T_ordering_matrix>::PtrType X, 
		     typename T_assign<T_ordering_matrix>::PtrTypePacked X_packed,
		     int const rows_total_matrix, int const cols_total_matrix) {
    // Loop over columns of kernel	
    for ( int col_k = 0; col_k < T_cols; ++col_k ) {
      // Loop over rows of kernel	
      for ( int row_k = 0; row_k < T_rows; ++row_k ) {
	T_assign<T_ordering_matrix>::exec( X, X_packed, 
					   row_k, col_k,
					   rows_total_matrix, cols_total_matrix );
      }
    }
  } // end exec()
  
  /** Convenience function for assignments to packed matrix.
   *  The template argument specifies how the original (unpacked) matrix
   *  is stored (e.g. column- or rowwise)
   */
  template<typename T_ordering_matrix>
    inline static void pack(real const * const X, real * X_packed,
			    int const rows_total_matrix, int const cols_total_matrix) {
    exec< Assign_to_packed, T_ordering_matrix >(X, X_packed, rows_total_matrix, cols_total_matrix);
  }
  /** Convenience function for extracting matrix from packed matrix.
   *  The template argument specifies how the unpacked matrix
   *  is stored (e.g. column- or rowwise)
   */
  template<typename T_ordering_matrix>
    inline static void unpack(real * X, real const * const X_packed,
			      int const rows_total_matrix, int const cols_total_matrix) {
    exec< Extract_from_packed, T_ordering_matrix >(X, X_packed, rows_total_matrix, cols_total_matrix);
  }
};
#endif
