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

#if !defined(_ERGO_BARRIER_H_)
#define _ERGO_BARRIER_H_
/** @file barrier.h declares a pthread-compatible barrier. This is
    to be used with older pthread implementations that do not provide
    barriers. */

#if !defined(HAS_PTHREAD_BARRIER)

#ifdef __cplusplus
#define EXTERN_C extern "C"
#else
#define EXTERN_C extern
#endif

#if !defined(PTHREAD_BARRIER_SERIAL_THREAD)
#define PTHREAD_BARRIER_SERIAL_THREAD -1
#endif

typedef struct ergo_barrier {
  pthread_mutex_t barrierMutex;
  pthread_cond_t  conditionVar;
  int initCount;
  int currCount;
  int cycle;
} ergo_barrier_t;

EXTERN_C int ergo_barrier_init(ergo_barrier_t * barrier,
                               const void * attr_ignored,
                               unsigned int count);

EXTERN_C int ergo_barrier_destroy (ergo_barrier_t *__barrier);
EXTERN_C int ergo_barrier_wait (ergo_barrier_t *__barrier);

#else /* HAS_PTHREAD_BARRIER */
#define ergo_barrier_t       pthread_barrier_t
#define ergo_barrier_init    pthread_barrier_init
#define ergo_barrier_destroy pthread_barrier_destroy
#define ergo_barrier_wait    pthread_barrier_wait
#endif /* HAS_PTHREAD_BARRIER */

#endif /* _ERGO_BARRIER_H_ */
