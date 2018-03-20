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

/** @file densfromf_full.cc

    @brief Routine get_dens_from_fock_full() for getting density
    matrix from a given Fock matrix using diagonalization.

    @author: Elias Rudberg <em>responsible</em>
*/

#include "densfromf_full.h"
#include "output.h"
#include <memory.h>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "memorymanag.h"
#include "machine_epsilon.h"
#include "utilities.h"
#include "matrix_algebra.h"
#include "units.h"

#include "mat_gblas.h"
#include "matInclude.h"



/** get_f_orbs: use diagonalization to find the molecular orbitals
 * corresponding to given Fock matrix f.
 */
int
get_F_orbs(int n, 
	   const ergo_real* F, 
	   const ergo_real* ovl, 
	   ergo_real* cmo, 
	   ergo_real* eigv)
{
  Util::TimeMeter timeMeter;
  static int ITYPE=1;
  int lwork = 10*n;
  ergo_real* work = (ergo_real*)ergo_malloc(lwork*sizeof(ergo_real));
  ergo_real* s    = (ergo_real*)ergo_malloc(n*n*sizeof(ergo_real));
  
  memcpy(cmo, F,   n*n*sizeof(ergo_real));
  memcpy(s,   ovl, n*n*sizeof(ergo_real));
  /* solve f*cmo = ovl*cmo*eigv; note that it destroys s! */
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "calling LAPACK routine sygv, n = %i", n);
  int info = -1;
  mat::sygv(&ITYPE, "V", "L", &n, cmo, &n, s, &n, 
	    eigv, work, &lwork, &info);
  timeMeter.print(LOG_AREA_DENSFROMF, "sygv diagonalization");
  if(info != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_DENSFROMF, "error in sygv");
      return -1;
    }

  ergo_free(work);
  ergo_free(s);

  return 0;
}


static void 
get_dens_from_cmo_zeroT(int n, 
			const ergo_real* cmo,
			const ergo_real* eigv,
			int noOfOccupiedOrbs, 
			ergo_real* dens,
			ergo_real & resultHomoLumoGap)
{
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "entering get_dens_from_cmo_zeroT, n = %i", n);
  Util::TimeMeter timeMeter;
  multiply_matrices_general_T_1(noOfOccupiedOrbs, n, noOfOccupiedOrbs, n, cmo, cmo, dens);
  for(int i = 0; i < n*n; i++)
    dens[i] *= 2;
  /*  compute bandgap */
  ergo_real E_HOMO = eigv[noOfOccupiedOrbs-1];
  ergo_real E_LUMO = eigv[noOfOccupiedOrbs-0];
  ergo_real E_min  = eigv[0];
  ergo_real E_max  = eigv[n-1];
  ergo_real bandGap = E_LUMO - E_HOMO;
  resultHomoLumoGap = bandGap;
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "get_dens_from_cmo finished, E_LUMO-E_HOMO = %12.8f Hartree = %12.8f eV", 
	    (double)bandGap,
	    (double)bandGap / UNIT_one_eV);
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "E_HOMO = %22.11f", (double)E_HOMO);
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "E_LUMO = %22.11f", (double)E_LUMO);
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "E_min  = %22.11f", (double)E_min );
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "E_max  = %22.11f", (double)E_max );
  timeMeter.print(LOG_AREA_DENSFROMF, "get_dens_from_cmo_zeroT");
}


static ergo_real x_times_ln_x(ergo_real x) {
  ergo_real eps = mat::getMachineEpsilon<ergo_real>();
  if(x < template_blas_sqrt(eps))
    return 0;
  return x * template_blas_log(x);
}

static void 
get_dens_from_cmo_FermiDiracDistr(int n, 
				  const ergo_real* cmo,
				  const ergo_real* eigv,
				  int noOfOccupiedOrbs, 
				  ergo_real* dens,
				  ergo_real electronicTemperature,
				  ergo_real & resultEntropyTerm)
{
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "entering get_dens_from_cmo_FermiDiracDistr, n = %i, nocc = %i, T = %12.6f", n, noOfOccupiedOrbs, electronicTemperature);
  Util::TimeMeter timeMeter;
  // Set occupation numbers according to Fermi-Dirac distribution.
  ergo_real T = electronicTemperature;
  ergo_real k_B = 1; // we use units such that k_B = 1
  ergo_real mu_min  = eigv[0];
  ergo_real mu_max  = eigv[n-1];
  ergo_real trace = 0;
  ergo_real mu;
  std::vector<ergo_real> occupationNumbers(n);
  int maxiter = 1000;
  int niter = 0;
  while(1) {
    ergo_real mu_hi = mu_max;
    ergo_real mu_lo = mu_min;
    while(1) {
      niter++;
      mu = (mu_hi + mu_lo) / 2;
      for(int i = 0; i < n; i++)
	occupationNumbers[i] = (ergo_real)1 / (template_blas_exp((eigv[i] - mu) / (k_B * T)) + 1);
      trace = 0;
      for(int i = 0; i < n; i++)
	trace += occupationNumbers[i];
      if(trace < noOfOccupiedOrbs)
	mu_lo = mu;
      else
	mu_hi = mu;
      // Emanuel comment: use some kind of relative error here,
      // otherwise we do not converge if mu is small or large
      if(mu_hi - mu_lo <= (template_blas_fabs(mu_hi)+template_blas_fabs(mu_lo))*100*mat::getMachineEpsilon<ergo_real>())
	break;
      if (niter >= maxiter)
	throw "Reached maxiter in Fermi function.";
    }
    // Check relative occupation error, if large we try to expand the
    // search interval
    if(template_blas_fabs(noOfOccupiedOrbs - trace)  < 
       std::abs(noOfOccupiedOrbs)*10000*mat::getMachineEpsilon<ergo_real>())
      break;
    ergo_real mu_min_max_diff = mu_max - mu_min;
    mu_min = mu_min - mu_min_max_diff;
    mu_max = mu_max + mu_min_max_diff;
  }
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "In get_dens_from_cmo_FermiDiracDistr, final trace = %12.6f, chemical potential = %12.6f", trace, mu);
  for(int j = 0; j < n; j++)
    for(int k = 0; k < n; k++)
      dens[j*n+k] = 0;
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++)
      for(int k = 0; k < n; k++)
	dens[j*n+k] += occupationNumbers[i] * cmo[i*n+j] * cmo[i*n+k];
  }
  for(int i = 0; i < n*n; i++)
    dens[i] *= 2;
  // Now compute electronic entropy term.
  resultEntropyTerm = 0;
  for(int i = 0; i < n; i++) {
    ergo_real lambda_i = occupationNumbers[i];
    resultEntropyTerm += (k_B * T) * ( x_times_ln_x(lambda_i) + x_times_ln_x(1-lambda_i) );
  }
  resultEntropyTerm *= 2;
  timeMeter.print(LOG_AREA_DENSFROMF, "get_dens_from_cmo_FermiDiracDistr");
}


int 
get_dens_from_fock_full(int n, 
			int noOfOccupiedOrbs, 
			ergo_real* result_P, 
			const ergo_real* F, 
			const ergo_real* ovl, 
			ergo_real factor,
			ergo_real electronicTemperature,
			ergo_real & resultEntropyTerm,
			ergo_real & resultHomoLumoGap,
			int store_all_eigenvalues_to_file, 
			ergo_real * const lumoVec,
			ergo_real * const homoVec)
{
  if(noOfOccupiedOrbs == 0) {
    memset(result_P, 0, n*n*sizeof(ergo_real));
    return 0;
  }
  Util::TimeMeter timeMeter;
  ergo_real* cmo  = (ergo_real*)ergo_malloc(n*n*sizeof(ergo_real));
  ergo_real* eigv = (ergo_real*)ergo_malloc(n*sizeof(ergo_real));
  if(get_F_orbs(n, F, ovl, cmo, eigv) != 0)
    return -1;
  if(electronicTemperature == 0)
    get_dens_from_cmo_zeroT(n, cmo, eigv, noOfOccupiedOrbs, result_P, resultHomoLumoGap);
  else {
    get_dens_from_cmo_FermiDiracDistr(n, cmo, eigv, noOfOccupiedOrbs, result_P, electronicTemperature, resultEntropyTerm);
    resultHomoLumoGap = -1; // Gap not defined in FermiDiracDistr case
  }

  // Save all eigenvalues to the file.
  // Functionality is used for plotting eigenspectrum of the matrix F.
  if(store_all_eigenvalues_to_file)
    {	
      do_output(LOG_CAT_WARNING, LOG_AREA_SCF, "Storing eigenvalues of the matrix F into the file.");

      std::ofstream file;
      file.open ("eigv_F.out");
      if (!file.is_open())
	throw "Error opening file for storing eigenvalues.";

      for (int ind = 0; ind < n; ind++) 
	file << std::fixed << std::setprecision(12) << (double)eigv[ind] << std::endl;
      file.close();  
    } // end of store_all_eigenvalues_to_file

  if ( lumoVec ) 
    for (int i = 0; i < n; i++) 
      lumoVec[i] = cmo[n*(noOfOccupiedOrbs)+i];
  if ( homoVec ) 
    for (int i = 0; i < n; i++) 
      homoVec[i] = cmo[n*(noOfOccupiedOrbs-1)+i];
  ergo_free(cmo);
  ergo_free(eigv);  

  // Take factor into account (factor is 2 for restricted case, 1 for unrestricted case).
  for(int i = 0; i < n*n; i++)
    result_P[i] *= factor / 2;
  resultEntropyTerm *= factor / 2;

  timeMeter.print(LOG_AREA_DENSFROMF, "get_dens_from_fock_full");
  return 0;
}


