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

/** @file common.h
 *
 * @brief Macros for inlining and static assertions and structs for
 * access to matrix elements specifying the layout of matrix entries
 * in memory.
 *
 * @author Emanuel H. Rubensson
 * @date 2009
 *
 */

#ifndef COMMON_H
#define COMMON_H
#include <cassert>

#define ALWAYS_INLINE  __attribute__((__always_inline__))
//#define ALWAYS_INLINE

//#define DEBUG_ON 

/** Class template for use in static asserts. 
 *
 */
template<bool>
struct CompileTimeChecker {
  CompileTimeChecker(...){}
};
/** Specialization of class template for use in static asserts. 
 *
 */
template<>
struct CompileTimeChecker<false> {
};
#define STATIC_ASSERT_ALWAYS(expr, msg)					\
  {									\
    class ERROR_##msg {};						\
    (CompileTimeChecker<(expr) != 0>(ERROR_##msg()));		\
  }

#ifdef DEBUG_ON
#define STATIC_ASSERT_DEBUG(expr, msg) STATIC_ASSERT_ALWAYS(expr, msg)
#else
#define STATIC_ASSERT_DEBUG(expr, msg)
#endif

//    (void)sizeof(CompileTimeChecker<(expr) != 0>((ERROR_##msg())));	\


// Store leading dimension (template argument) as static const
// Then one can either use "get" function (ROWS, COLS args not needed?) or 
// specialize templates depending on the type (transposed or regular).


/** Struct for access to matrix elements stored in row wise order.
 *  This struct used to specify how and in which order matrix elements are 
 *  stored. At the moment, only regular row or column wise ordering is 
 *  supported, but one could imagine symmetric or triangular storage.  
 * @see Ordering_col_wise
 */
struct Ordering_row_wise {
  inline static int get( int const row, int const col, 
			 int const rows, int const cols ) {
    return row * cols + col;
  }
  template<int T_row, int T_col, int T_rows, int T_cols>
    struct Get {
      static int const index = T_row * T_cols + T_col;
    };

};

/** Struct for access to matrix elements stored in column wise order.
 * @see Ordering_row_wise
 */
struct Ordering_col_wise {
  inline static int get( int const row, int const col, 
			 int const rows, int const cols ) {
    return row + col * rows;
  }
  template<int T_row, int T_col, int T_rows, int T_cols>
    struct Get {
      static int const index = T_row + T_col * T_rows;
    };
};


#endif
