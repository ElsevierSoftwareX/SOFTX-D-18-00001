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

/** @file g_intrin.h
 *
 * @brief Templates for convenient access to intrinsic instructions.
 *
 * @author Emanuel H. Rubensson
 * @date 2009
 *
 */

#ifndef G_INTRIN
#define G_INTRIN
#include <emmintrin.h>
#ifdef __SSE3__
#include <pmmintrin.h>
#endif

/* Interface to load functions. */

/* load_p */
template<typename Treal, typename Treg>
  inline static Treg _mm_load_p (Treal const * ptr);

inline static __m128 _mm_load_p (float const * ptr) {
  return _mm_load_ps (ptr);
}

inline static __m128d _mm_load_p (double const * ptr) {
  return _mm_load_pd (ptr);
}

/* load1_p */
template<typename Treal, typename Treg>
  inline static Treg _mm_load1_p (Treal const * ptr);

inline static __m128 _mm_load1_p (float const * ptr) {
  return _mm_load1_ps (ptr);
}


inline static __m128d _mm_load1_p (double const * ptr) {
  return _mm_load1_pd (ptr);
}

/* set1_p */
template<typename Treal, typename Treg>
  inline static Treg _mm_set1_p (Treal const val);

inline static __m128 _mm_set1_p (float const val) {
  return _mm_set1_ps (val);
}


inline static __m128d _mm_set1_p (double const val) {
  return _mm_set1_pd (val);
}


/* Interface to store functions. */
template<typename Treal, typename Treg>
  inline static void _mm_store_p (Treal * ptr, Treg A);

inline static void  _mm_store_p (float * ptr, __m128 A) {
  _mm_store_ps (ptr, A);
}

inline static void  _mm_store_p (double * ptr, __m128d A) {
  _mm_store_pd (ptr, A);
}


/* Interface to add functions. */

template<typename Treg>
inline static Treg _mm_add_p (Treg A, Treg B);

inline static __m128 _mm_add_p (__m128 A, __m128 B) {
  return _mm_add_ps(A, B);
}

inline static __m128d _mm_add_p (__m128d A, __m128d B) {
  return _mm_add_pd(A, B);
}


/* Interface to mul functions. */

template<typename Treg>
inline static Treg _mm_mul_p (Treg A, Treg B);

inline static __m128 _mm_mul_p (__m128 A, __m128 B) {
  return _mm_mul_ps(A, B);
}

inline static __m128d _mm_mul_p (__m128d A, __m128d B) {
  return _mm_mul_pd(A, B);
}

/* pxor */

template<typename Treg>
inline static Treg _mm_xor_p (Treg A, Treg B);

inline static __m128 _mm_xor_p (__m128 A, __m128 B) {
  return _mm_xor_ps(A, B);
}

inline static __m128d _mm_xor_p (__m128d A, __m128d B) {
  return _mm_xor_pd(A, B);
}

#endif
