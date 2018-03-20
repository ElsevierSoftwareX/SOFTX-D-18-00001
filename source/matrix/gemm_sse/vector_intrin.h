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

/** @file vector_intrin.h
 *
 * @brief Vector template for convenient access to SIMD operations.
 *
 * @author Emanuel H. Rubensson
 * @date 2009
 *
 */
/*
 * Templates for efficient gemm kernels.
 * For architectures with SSE2 or higher.
 *
 * Copyright Emanuel Rubensson, 2009 
 *
 *
 *
 *
 */

#ifndef VECTOR_INTRIN
#define VECTOR_INTRIN
#include "common.h"
#include "g_intrin.h"

/** Vector class template for access to SIMD operations.
 *
 *  Currently supports a limited set of double and single precision SSE operations.
 *
 */
template<typename Treal, typename Treg>
  class Vector_intrin {
 public:
  inline void ALWAYS_INLINE load_p(Treal const * ptr) {
    values = _mm_load_p(ptr);
  }
  inline void ALWAYS_INLINE load1_p(Treal const * ptr) {
    values = _mm_load1_p(ptr);
  }
  inline void ALWAYS_INLINE store_p(Treal * ptr) const {
    _mm_store_p ( ptr, values ); 
  }
  inline Vector_intrin<Treal, Treg>& ALWAYS_INLINE operator*= ( Vector_intrin<Treal, Treg> const & other ) {
    values = _mm_mul_p ( other.values, values );
    return *this;
  }
  inline Vector_intrin<Treal, Treg>& ALWAYS_INLINE operator+= ( Vector_intrin<Treal, Treg> const & other ) {
    values = _mm_add_p ( other.values, values );
    return *this;
  }
  inline Vector_intrin<Treal, Treg>& ALWAYS_INLINE operator+= ( Treal const * ptr ) {
    Treg tmp;
    tmp = _mm_load_p(ptr);
    values = _mm_add_p ( tmp, values );
    return *this;
  }
#if 0
  inline void ALWAYS_INLINE xor_p( Vector_intrin<Treal, Treg> const & other ) {
    values = _mm_xor_p ( other.values, values );
  }
#endif
  inline void ALWAYS_INLINE set_to_zero() {
    values = _mm_xor_p ( values, values );
  }
 protected:
  Treg values;
 private:
};

template<typename Treal>
class Vector_intrin<Treal, Treal> {
 public:
  inline void ALWAYS_INLINE load_p(Treal const * ptr) {
    values = *ptr;
  }
  inline void ALWAYS_INLINE load1_p(Treal const * ptr) {
    values = *ptr;
  }
  inline void ALWAYS_INLINE store_p(Treal * ptr) const {
    *ptr = values; 
  }
  inline Vector_intrin<Treal, Treal>& ALWAYS_INLINE operator*= ( Vector_intrin<Treal, Treal> const & other ) {
    values *= other.values;
    return *this;
  }
  inline Vector_intrin<Treal, Treal>& ALWAYS_INLINE operator+= ( Vector_intrin<Treal, Treal> const & other ) {
    values += other.values;
    return *this;
  }
  inline Vector_intrin<Treal, Treal>& ALWAYS_INLINE operator+= ( Treal const * ptr ) {
    values += *ptr;
    return *this;
  }
#if 0
  inline void ALWAYS_INLINE xor_p( Vector_intrin<Treal, Treg> const & other ) {
    values = _mm_xor_p ( other.values, values );
  }
#endif
  inline void ALWAYS_INLINE set_to_zero() {
    values = 0;
  }
 protected:
  Treal values;
 private:
};






#if 0
template<>
class Vector_intrin<double, double> {
 public:
  inline void ALWAYS_INLINE load_p(double const * ptr) {
    values[0] = *ptr;
    values[1] = ptr[1];
  }
  inline void ALWAYS_INLINE load1_p(double const * ptr) {
    values[0] = *ptr;
    values[1] = *ptr;
  }
  inline void ALWAYS_INLINE store_p(double * ptr) const {
    ptr[0] = values[0]; 
    ptr[1] = values[1]; 
  }
  inline Vector_intrin<double, double>& ALWAYS_INLINE operator*= ( Vector_intrin<double, double> const & other ) {
    values[0] *= other.values[0];
    values[1] *= other.values[1];
  }
  inline Vector_intrin<double, double>& ALWAYS_INLINE operator+= ( Vector_intrin<double, double> const & other ) {
    values[0] += other.values[0];
    values[1] += other.values[1];
  }
  inline Vector_intrin<double, double>& ALWAYS_INLINE operator+= ( double const * ptr ) {
    values[0] += *ptr;
    values[1] += ptr[1];
  }
#if 0
  inline void ALWAYS_INLINE xor_p( Vector_intrin<double, Treg> const & other ) {
    values = _mm_xor_p ( other.values, values );
  }
#endif
  inline void ALWAYS_INLINE set_to_zero() {
    values[0] = 0;
    values[1] = 0;
  }
 protected:
  double values[2];
 private:
};
#endif

#endif // VECTOR_INTRIN
