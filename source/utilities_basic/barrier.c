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

/** @file barrier.c implements a pthread-compatible barrier. This is
    to be used with older pthread implementations that do not provide
    barriers. This implementation is applicable only in simple
    cases. Check section 7.1.1 of "Programming with POSIX threads" for
    a full-blown implementation.  In particular, this implementation
    does not check for some error conditions, like destroying the
    barrier when some threads wait on it. */

#include <errno.h>
#include <pthread.h>

#include "barrier.h"

int
ergo_barrier_init(ergo_barrier_t * barrier,
                  const void * attr_ignored,
                  unsigned int count)
{
  int r;

  /* This may be called only once for given barrier! */
  if( (r=pthread_mutex_init(&barrier->barrierMutex, NULL)) ) return r;
  if( (r=pthread_cond_init(&barrier->conditionVar, NULL)) )  return r;
  barrier->currCount = barrier->initCount = count;
  return 0;
}

int
ergo_barrier_destroy(ergo_barrier_t *barrier)
{
  int r;
  if( (r=pthread_mutex_destroy(&barrier->barrierMutex)) ) return r;
  if( (r=pthread_cond_destroy(&barrier->conditionVar)) )  return r;

  return 0;
}
     
int
ergo_barrier_wait(ergo_barrier_t *barrier)
{
  int r, ret;

  if( (r=pthread_mutex_lock(&barrier->barrierMutex)) ) return r;
  
  if (barrier->currCount == 0)
    return EINVAL;

  --barrier->currCount;

  if(barrier->currCount == 0) { /* I am the retarded one! :) */
    ++barrier->cycle;
    barrier->currCount = barrier->initCount;
    if( (r=pthread_cond_broadcast(&barrier->conditionVar)) ) return r;
    ret = PTHREAD_BARRIER_SERIAL_THREAD;
  } else {
    int curCycle = barrier->cycle;
    while(barrier->cycle == curCycle) {
      if( (r=pthread_cond_wait(&barrier->conditionVar,
                               &barrier->barrierMutex)) ) return r;
    }
    ret = 0;
  }
  if( (r=pthread_mutex_unlock(&barrier->barrierMutex)) ) return r;
  return ret;
}

#ifdef TEST
static void* thr(void *p)
{
  ergo_barrier_t *b = (ergo_barrier_t*)p;
  pthread_t a;
  int r;

  a  = pthread_self();
  r = 1 + (((int)a) >> 12) & 0x03; /* this will hopefully be a bit of
                                      a random number */
  printf("%x sleeps for %i s\n", (int)a, r);
  sleep(r);
  printf("%x waits on barrier\n", (int)a);
  r = ergo_barrier_wait(b);
  printf("%x continues with r=%d\n", (int)a, r);
  return NULL;
}
int main(int argc, char *argv)
{
#define THREAD_COUNT 10
  pthread_t pid[THREAD_COUNT];
  ergo_barrier_t barrier;
  int i;

  ergo_barrier_init(&barrier, NULL, THREAD_COUNT);
  for(i=0; i<THREAD_COUNT; i++)
    pthread_create(&pid[i], NULL, thr, &barrier);

  for(i=0; i<THREAD_COUNT; i++) {
    pthread_join(pid[i], NULL);
  }
  ergo_barrier_destroy(&barrier);

  return 0;
}
#endif /* TEST */
