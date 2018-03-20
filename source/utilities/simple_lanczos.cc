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

/** @file simple_lanczos.cc

    @brief Simple implementation of the Lanczos method.

    @author: Elias Rudberg <em>responsible</em>
*/

#include "../matrix/mat_gblas.h"
#include "simple_lanczos.h"
#include "output.h"
#include <cmath>
#include <cstring>

namespace simple_lanczos {

  ergo_real simple_lanczos_get_vector_norm(int n, const ergo_real* v) {
    ergo_real sqSum = 0;
    for(int i = 0; i < n; i++)
      sqSum += v[i] * v[i];
    return template_blas_sqrt(sqSum);
  }

  void simple_lanczos_normalize_vector(int n, ergo_real* v) {
    ergo_real factor = 1.0 / simple_lanczos_get_vector_norm(n, v);
    for(int i = 0; i < n; i++)
      v[i] *= factor;
  }

  void simple_lanczos_get_eigs(int n, ergo_real* M, 
			       ergo_real & currEig_lo, ergo_real* bestVector_lo, 
			       ergo_real & currEig_hi, ergo_real* bestVector_hi, 
			       ergo_real* eigValListResult)
  {
    int lwork = 3*n*n;
    ergo_real* work = new ergo_real[lwork];
    ergo_real* eigvalList = new ergo_real[n];
    ergo_real* A = new ergo_real[n*n];
    memcpy(A, M, n*n*sizeof(ergo_real));
    int info = 0;
    mat::syev("V", "U", &n, A,
	      &n, eigvalList, work, &lwork, 
	      &info);
    if(info != 0)
      {
	do_output(LOG_CAT_ERROR, LOG_AREA_CI, "error in syev, info = %i", info);
	throw("error in syev");
      }

    if(bestVector_lo) {
      for(int i = 0; i < n; i++)
	bestVector_lo[i] = A[0*n+i];
    }
    if(bestVector_hi) {
      for(int i = 0; i < n; i++)
	bestVector_hi[i] = A[(n-1)*n+i];
    }

    if(eigValListResult)
      {
	for(int i = 0; i < n; i++)
	  eigValListResult[i] = eigvalList[i];
      }

    currEig_lo = eigvalList[0];
    currEig_hi = eigvalList[n-1];

    delete [] work;
    delete [] eigvalList;
    delete [] A;
  }

} // end namespace simple_lanczos
