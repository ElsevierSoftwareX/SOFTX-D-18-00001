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
 
 /* This file belongs to the template_lapack part of the Ergo source 
  * code. The source files in the template_lapack directory are modified
  * versions of files originally distributed as CLAPACK, see the
  * Copyright/license notice in the file template_lapack/COPYING.
  */
 


#include "template_blas_common.h"

#include <stdexcept>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <pthread.h>

typedef double realtype;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;


static void*
thread_func(void* arg)
{
  for(int loop_k = 0; loop_k < 22222; loop_k++)
    {
      realtype A[4];
      A[0] = 1.3;  A[1] = 2.4;
      A[2] = 1.1;  A[3] = 2.7;
  
      realtype B[4];
      B[0] = 0.3;  B[1] = 0.4;
      B[2] = 0.1;  B[3] = 0.7;
  
      realtype C[4];
      C[0] = 0;  C[1] = 0;
      C[2] = 0;  C[3] = 0;
  
      int m = 2;
      int n = 2;
      int k = 2;
      realtype alpha = 1.0;
      realtype beta = 0.0;
      template_blas_gemm("N", "N", &m, &n, &k, &alpha, A, &m, B, &m, &beta, C, &m);

      pthread_mutex_lock(&mutex);
      //  std::cout << "template_blas_gemm finished." << std::endl;
      pthread_mutex_unlock(&mutex);

#if 0
      for(int i = 0; i < n; i++)
	{
	  for(int j = 0; j < n; j++)
	    std::cout << "  " << C[i*n+j];
	  std::cout << std::endl;
	}
#endif

      realtype D[4];
      for(int i = 0; i < n; i++)
	for(int j = 0; j < n; j++)
	  {
	    realtype sum = 0;
	    for(int k = 0; k < n; k++)
	      sum += B[i*n+k] * A[k*n+j];
	    D[i*n+j] = sum;
	  }

#if 0
      for(int i = 0; i < n; i++)
	{
	  for(int j = 0; j < n; j++)
	    std::cout << "  " << D[i*n+j];
	  std::cout << std::endl;
	}
#endif
      
      for(int i = 0; i < n; i++)
	for(int j = 0; j < n; j++)
	  {
	    if(fabs(C[i*n+j] - D[i*n+j]) > 1e-11)
	      throw std::runtime_error("ERROR! wrong result from gemm!");
	  }
    }

  pthread_mutex_lock(&mutex);
  std::cout << "Thread " << pthread_self() << " finished OK." << std::endl;
  pthread_mutex_unlock(&mutex);

  return NULL;
}


int main()
{
  const int nThreads = 8;
  pthread_t threads[nThreads];
  for(int i = 0; i < nThreads; i++) {
    if(pthread_create(&threads[i], NULL, thread_func, NULL) != 0)
      throw std::runtime_error("Error in pthread_create.");
  }
  for(int i = 0; i < nThreads; i++) {
    if(pthread_join(threads[i], NULL) != 0)
      throw std::runtime_error("Error in pthread_join.");
  }

  realtype piTest = template_blas_compute_pi_BBP((realtype)0);
  printf("pi = %22.15f\n", (double)piTest);

  printf("template_blas_test_threaded test ended OK.\n");

  return 0;
}

