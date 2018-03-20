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

/** @file simple_lanczos.h

    @brief Simple implementation of the Lanczos method.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef SIMPLE_LANCZOS_HEADER
#define SIMPLE_LANCZOS_HEADER

#include "realtype.h"
#include <vector>
#include <cstdio>

namespace simple_lanczos {

  ergo_real simple_lanczos_get_vector_norm(int n, const ergo_real* v);
  void simple_lanczos_normalize_vector(int n, ergo_real* v);
  void simple_lanczos_get_eigs(int n, ergo_real* M, 
			       ergo_real & currEig_lo, ergo_real* bestVector_lo, 
			       ergo_real & currEig_hi, ergo_real* bestVector_hi, 
			       ergo_real* eigValListResult);

  template<typename Tmatvecmul>
    void do_lanczos_method(int n,
			   const ergo_real* guessVector,
			   ergo_real & resultEig_lo,
			   ergo_real* resultVec_lo,
			   ergo_real & resultEig_hi,
			   ergo_real* resultVec_hi,
			   const Tmatvecmul & matvecmul,
			   int maxIterations_in,
			   ergo_real shift,
			   ergo_real extraEnergy) {
    if(n == 1) {
      // Special case for n=1, in this case we need only one "matrix-vector" (really scalar) multiplication to get all info we need.
      ergo_real tmpVec1[1];
      tmpVec1[0] = 1;
      ergo_real tmpVec2[1];
      tmpVec2[0] = 0;
      matvecmul.do_mat_vec_mul(n, tmpVec1, tmpVec2);
      ergo_real eigenValue = tmpVec2[0];
      resultEig_lo = eigenValue;
      resultEig_hi = eigenValue;
      resultVec_lo[0] = 1;
      resultVec_hi[0] = 1;
    }
    typedef ergo_real* ergo_real_ptr;
    int maxIterations = maxIterations_in;
    if(maxIterations > n)
      maxIterations = n;
    ergo_real** q = new ergo_real_ptr[n+1];
    q[0] = new ergo_real[n];
    for(int i = 0; i < n; i++)
      q[0][i] = 0;
    q[1] = new ergo_real[n];
    for(int i = 0; i < n; i++)
      q[1][i] = guessVector[i];
    simple_lanczos_normalize_vector(n, q[1]);
    std::vector<ergo_real> z(n);
    std::vector<ergo_real> alpha(n+1);
    std::vector<ergo_real> beta(n+1);
    beta[0] = 0;
    std::vector<ergo_real> bestVector_lo(maxIterations+1);
    std::vector<ergo_real> bestVector_hi(maxIterations+1);
    ergo_real currEig_lo = 0;
    ergo_real currEig_hi = 0;
    ergo_real curr_E_lo = 0;
    ergo_real curr_E_hi = 0;
    for(int j = 1; j <= maxIterations; j++) {
      // Do matrix-vector multiplication
      matvecmul.do_mat_vec_mul(n, q[j], &z[0]);
      // OK, matrix-vector multiplication done
      alpha[j] = 0;
      for(int i = 0; i < n; i++)
	alpha[j] += q[j][i] * z[i];
      for(int i = 0; i < n; i++)
	z[i] = z[i] - alpha[j] * q[j][i] - beta[j-1] * q[j-1][i];
      beta[j] = simple_lanczos_get_vector_norm(n, &z[0]);
      ergo_real* T = new ergo_real[j*j];
      for(int i = 0; i < j*j; i++)
	T[i] = 0;
      for(int i = 0; i < j; i++)
	T[i*j+i] = alpha[i+1];
      for(int i = 0; i < j-1; i++) {
	T[i*j+(i+1)] = beta[i+1];
	T[(i+1)*j+i] = beta[i+1];
      }
      simple_lanczos_get_eigs(j, T, currEig_lo, &bestVector_lo[0], currEig_hi, &bestVector_hi[0], NULL);
      // Set resultVec_lo
      for(int k = 0; k < n; k++) {
	ergo_real sum = 0;
	for(int i = 1; i <= j; i++)
	  sum += bestVector_lo[i-1] * q[i][k];
	resultVec_lo[k] = sum;
      }
      // Set resultVec_hi
      for(int k = 0; k < n; k++) {
	ergo_real sum = 0;
	for(int i = 1; i <= j; i++)
	  sum += bestVector_hi[i-1] * q[i][k];
	resultVec_hi[k] = sum;
      }
      curr_E_lo = currEig_lo + extraEnergy + shift;
      curr_E_hi = currEig_hi + extraEnergy + shift;
      if(beta[j] < 1e-11 && j > 1)
	break;
      if(j == maxIterations)
	break;
      q[j+1] = new ergo_real[n];
      for(int i = 0; i < n; i++)
	q[j+1][i] = z[i] / beta[j];
    } // end for j
    resultEig_lo = curr_E_lo;
    resultEig_hi = curr_E_hi;
    simple_lanczos_normalize_vector(n, resultVec_lo);
    simple_lanczos_normalize_vector(n, resultVec_hi);
  }

} // end namespace simple_lanczos

#endif
