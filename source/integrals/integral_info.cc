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

/** @file integral_info.cc

    @brief Defines IntegralInfo object, providing the coefficients
    needed for integral evaluation.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <memory.h>
#include <time.h>
#include <stdarg.h>
#include <stdexcept>
#include "integral_info.h"
#include "boysfunction.h"
#include "output.h"
#include "memorymanag.h"
#include "mat_gblas.h"



#define NBIN 28
static int BinCoeffs[NBIN*NBIN];

static void
setup_bin_coeffs()
{
  for(int i = 0; i < NBIN; i++) {
    BinCoeffs[i*NBIN+0] = 1;
    BinCoeffs[i*NBIN+i] = 1;
    for(int j = 1; j < i; j++) {
      BinCoeffs[i*NBIN+j] = 
	BinCoeffs[(i-1)*NBIN+j-1]
	+ BinCoeffs[(i-1)*NBIN+j];
    }
  }
}

static int getBinCoeff(int i, int j) {
  if(i >= NBIN || j >= NBIN)
    throw "Error in integral_info getBinCoeff: (i >= NBIN || j >= NBIN).";
  return BinCoeffs[i*NBIN+j];
}


/* Earlier, the factorial() routine here had return type int, but that
   gave problems with integer overflow. */
static ergo_real
factorial(int n) {
  if(n == 0)
    return 1;
  return n * factorial(n-1);
}

static int
get_real_solid_harmonic_poly(int l, int m,
			     basis_func_poly_struct* result) {
  setup_bin_coeffs();  
  
  ergo_real denominator;
  if(m == 0)
    denominator = 2;
  else
    denominator = 1;
  ergo_real NSlm = ((ergo_real)1 / (template_blas_pow((ergo_real)2, (ergo_real)abs(m))*factorial(l))) 
    * template_blas_sqrt(2*factorial(l+abs(m))*factorial(l-abs(m)) / denominator);  

  const int MAX_DEGREE = 10;
  ergo_real terms[MAX_DEGREE][MAX_DEGREE][MAX_DEGREE];
  memset(terms, 0, sizeof(terms));  

  for(int t = 0; t <= (l-abs(m))/2; t++)
    for(int u = 0; u <= t; u++) {
      ergo_real vm;
      int n; // n is the number of terms in the innermost sum
      if(m >= 0) {
	// m >= 0  ==>  vm = 0
	vm = 0;
	n = abs(m) / 2 + 1;
      }
      else
	{
	  // m < 0  ==>  vm = 0.5
	  vm = 0.5;
	  n = (abs(m) - 1) / 2 + 1;
	}
      for(int v_index = 0; v_index < n; v_index++) {
	ergo_real v = v_index + vm;
	int twov = (int)(2 * v);
	int ix = (int)(2*t + abs(m) - 2 * (u + v));
	int iy = (int)(2*(u+v));
	int iz = (int)(l - 2*t - abs(m));
	if(ix+iy+iz != l) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, 
		    "error in get_real_solid_harmonic_poly: ix iy iz = %i %i %i", ix, iy, iz);
	  return -1;
	}
	if(ix >= MAX_DEGREE || iy >= MAX_DEGREE || iz >= MAX_DEGREE) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, 
		    "error in get_real_solid_harmonic_poly: "
		    "(ix >= MAX_DEGREE || iy >= MAX_DEGREE || iz >= MAX_DEGREE), ix iy iz = %i %i %i", ix, iy, iz);
	  return -1;
	}
	if(ix < 0 || iy < 0 || iz < 0) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: (ix < 0 || iy < 0 || iz < 0)");
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "ix iy iz = %i %i %i", ix, iy, iz);
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "t = %i, m = %i, u = %i, v = %f", t, m, u, (double)v);
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "n = %i, v_index = %i", n, v_index);
	  return -1;
	}

	int power_t_plus_v_minus_vm = (int)(t + v - vm);
	ergo_real Clmtuv = 
			  template_blas_pow((ergo_real)-1, (ergo_real)power_t_plus_v_minus_vm) * 
	  template_blas_pow((ergo_real)0.25, (ergo_real)t) *
	  getBinCoeff(l, t)   * getBinCoeff(l-t, abs(m)+t)     * getBinCoeff(t, u)   * getBinCoeff(abs(m), twov);

	terms[ix][iy][iz] += NSlm * Clmtuv;
      } // END FOR v_index
    } // END FOR t u

  int termCount = 0;
  for(int ix = 0; ix < MAX_DEGREE; ix++)
    for(int iy = 0; iy < MAX_DEGREE; iy++)
      for(int iz = 0; iz < MAX_DEGREE; iz++) {
	if(terms[ix][iy][iz] != 0) {
	  result->termList[termCount].monomialInts[0] = ix;
	  result->termList[termCount].monomialInts[1] = iy;
	  result->termList[termCount].monomialInts[2] = iz;
	  result->termList[termCount].coeff = terms[ix][iy][iz];
	  termCount++;
	  if(termCount >= MAX_NO_OF_TERMS_IN_BASIS_FUNC_POLY) {
	    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, 
		      "error in get_real_solid_harmonic_poly: (termCount >= MAX_NO_OF_TERMS_IN_BASIS_FUNC_POLY)");
	    return -1;
	  }
	}
      } // END FOR ix iy iz

  result->noOfTerms = termCount;
    
  return 0;
}




int
setup_basis_func_polys(IntegralInfo* b)
{
  basis_func_poly_struct* curr = NULL;
  int count = 0;

  const int MAX_L_QUANTUM_NUMBER = BASIS_FUNC_POLY_MAX_DEGREE;
  ergo_real scaleFactorList[MAX_L_QUANTUM_NUMBER+1];
  scaleFactorList[0] = 1;
  scaleFactorList[1] = 1;
  scaleFactorList[2] = template_blas_sqrt(3.0);
  scaleFactorList[3] = template_blas_sqrt(15.0);
  // Set all other factors to same value. FIXME: find out if/how/why this matters.
  for(int ii = 4; ii <= MAX_L_QUANTUM_NUMBER; ii++)
    scaleFactorList[ii] = template_blas_sqrt(15.0); 
  
  for(int l = 0; l <= MAX_L_QUANTUM_NUMBER; l++) {
    for(int m = -l; m <= l; m++) {
      curr = &b->basis_func_poly_list[count];
      if(get_real_solid_harmonic_poly(l, m, curr) != 0) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_real_solid_harmonic_poly");
	return -1;
      }

      // Now curr contains the Solid Harmonic polynamial as given in table 6.3 in the book by Helgaker at al.
      // compute scaledSolidHarmonicPrefactor
      if(m == 0)
	curr->scaledSolidHarmonicPrefactor = 1;
      else
	curr->scaledSolidHarmonicPrefactor = template_blas_pow((ergo_real)-1, (ergo_real)m) / template_blas_sqrt((ergo_real)2);

      // use scalefactor
      for(int i = 0; i < curr->noOfTerms; i++)
	curr->termList[i].coeff /= scaleFactorList[l];

      // update curr->scaledSolidHarmonicPrefactor to compensate for that each term has been divided by scaleFactor
      curr->scaledSolidHarmonicPrefactor *= scaleFactorList[l];

      count++;
      if(count >= MAX_NO_OF_BASIS_FUNC_POLYS) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, 
		  "error in setup_basis_func_polys: (count >= MAX_NO_OF_BASIS_FUNC_POLYS)");	      
	return -1;
      }
    } // END FOR m
  } // END FOR l

  /* set monomialID for each term */
  for(int j = 0; j < count; j++) {
    curr = &b->basis_func_poly_list[j];
    for(int i = 0; i < curr->noOfTerms; i++) {
      basis_func_term_struct* currTerm = &curr->termList[i];
      int i0 = currTerm->monomialInts[0];
      int i1 = currTerm->monomialInts[1];
      int i2 = currTerm->monomialInts[2];
      currTerm->monomialID = b->monomial_info.monomial_index_list[i0][i1][i2];
    }
  }

  b->no_of_basis_func_polys = count;

  return 0;
}

ergo_real IntegralInfo::BoysFunction(int n, ergo_real x) const {
  if(!initialized)
    throw std::runtime_error("Error in IntegralInfo::BoysFunction: not initialized.");
  return boysFunctionManager.BoysFunction(n, x);
}

ergo_real IntegralInfo::BoysFunction_expensive(int n, ergo_real x, int noOfIntegrationIntervals) const {
  return boysFunctionManager.BoysFunction_expensive(n, x, noOfIntegrationIntervals);
}

int IntegralInfo::multiply_by_hermite_conversion_matrix_from_right(int n1max,        
								   int n2max,        
								   ergo_real a,      
								   ergo_real* A,     
								   ergo_real* result) const {
  return hermite_conversion_info.multiply_by_hermite_conversion_matrix_from_right(monomial_info, n1max, n2max, a, A, result);
}

int IntegralInfo::multiply_by_hermite_conversion_matrix_from_left(int n1max,        
								  int n2max,        
								  ergo_real a,      
								  ergo_real* A,     
								  ergo_real* result) const {
  return hermite_conversion_info.multiply_by_hermite_conversion_matrix_from_left(monomial_info, n1max, n2max, a, A, result);
}

int IntegralInfo::get_hermite_conversion_matrix_right(int nmax,
						      ergo_real a,
						      ergo_real* result) const {
  return hermite_conversion_info.get_hermite_conversion_matrix_right(monomial_info, nmax, a, result);
}

int IntegralInfo::get_hermite_conversion_matrix_left(int nmax,
						     ergo_real a,
						     ergo_real* result) const {
  return hermite_conversion_info.get_hermite_conversion_matrix_left(monomial_info, nmax, a, result);
}

int IntegralInfo::get_hermite_conversion_matrix_right_sparse(int nmax,
							     ergo_real a,
							     i_j_val_struct* result) const {
  return hermite_conversion_info.get_hermite_conversion_matrix_right_sparse(monomial_info, nmax, a, result);  
}

void IntegralInfo::init() {
  if(initialized)
    return;
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "IntegralInfo::init() calling boysFunctionManager.init().");
  boysFunctionManager.init();
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "IntegralInfo::init() calling multipolePrep.init().");
  multipolePrep.init();
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "IntegralInfo::init() calling mmLimitTable.init().");
  mmLimitTable.inittt(multipolePrep);
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "IntegralInfo::init() calling monomial_info.init().");
  monomial_info.init();
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "IntegralInfo::init() calling hermite_conversion_info.init().");
  hermite_conversion_info.init(monomial_info);
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "IntegralInfo::init() calling setup_basis_func_polys");
  if(setup_basis_func_polys(this) != 0)
    throw std::runtime_error("Error in IntegralInfo::init(), in setup_basis_func_polys()..");    
  initialized = true;
}

IntegralInfo::IntegralInfo(bool initialize) : initialized(false)
{
  if(initialize)
    init();
}

IntegralInfo::~IntegralInfo()
{
  /* Nothing is dynamically allocated, nothing needs to be released. */
}

/** Function needed for Chunks&Tasks usage. */
IntegralInfo::IntegralInfo(const IntegralInfo & ii) 
  : boysFunctionManager(ii.boysFunctionManager), 
    multipolePrep(ii.multipolePrep),
    mmLimitTable(ii.mmLimitTable),
    hermite_conversion_info(ii.hermite_conversion_info), 
    initialized(ii.initialized), 
    no_of_basis_func_polys(ii.no_of_basis_func_polys), 
    monomial_info(ii.monomial_info)
{
  memcpy(basis_func_poly_list, ii.basis_func_poly_list, sizeof(basis_func_poly_list));
}

/** Function needed for Chunks&Tasks usage. */
void IntegralInfo::write_to_buffer ( char * dataBuffer, size_t const bufferSize ) const {
  if(!initialized)
    throw std::runtime_error("Error: IntegralInfo::write_to_buffer called when not initialized.");
  char* p = dataBuffer;
  if(bufferSize < get_size())
    throw std::runtime_error("Error in IntegralInfo::write_to_buffer: bufferSize too small.");
  // boysFunctionManager
  boysFunctionManager.write_to_buffer(p, bufferSize - (p - dataBuffer));
  p += boysFunctionManager.get_size();
  // multipolePrep
  multipolePrep.write_to_buffer(p, bufferSize - (p - dataBuffer));
  p += multipolePrep.get_size();
  // mmLimitTable
  mmLimitTable.write_to_buffer(p, bufferSize - (p - dataBuffer));
  p += mmLimitTable.get_size();
  // hermite_conversion_info
  hermite_conversion_info.write_to_buffer(p, bufferSize - (p - dataBuffer));
  p += hermite_conversion_info.get_size();
  // initialized
  memcpy(p, &initialized, sizeof(bool));
  p += sizeof(bool);
  // basis_func_poly_list
  memcpy(p, basis_func_poly_list, MAX_NO_OF_BASIS_FUNC_POLYS*sizeof(basis_func_poly_struct));
  p += MAX_NO_OF_BASIS_FUNC_POLYS*sizeof(basis_func_poly_struct);
  // no_of_basis_func_polys
  memcpy(p, &no_of_basis_func_polys, sizeof(int));
  p += sizeof(int);
  // monomial_info
  monomial_info.write_to_buffer(p, bufferSize - (p - dataBuffer));
  p += monomial_info.get_size();
  // DONE!
}

/** Function needed for Chunks&Tasks usage. */
size_t IntegralInfo::get_size() const {
  return boysFunctionManager.get_size()
    + multipolePrep.get_size()
    + mmLimitTable.get_size()
    + hermite_conversion_info.get_size()
    + sizeof(bool)
    + MAX_NO_OF_BASIS_FUNC_POLYS*sizeof(basis_func_poly_struct)
    + sizeof(int)
    + monomial_info.get_size();
}

/** Function needed for Chunks&Tasks usage. */
void IntegralInfo::assign_from_buffer ( char const * dataBuffer, size_t const bufferSize) {
  const char* p = dataBuffer;
  // boysFunctionManager
  boysFunctionManager.assign_from_buffer(p, bufferSize - (p - dataBuffer));
  p += boysFunctionManager.get_size();
  // multipolePrep
  multipolePrep.assign_from_buffer(p, bufferSize - (p - dataBuffer));
  p += multipolePrep.get_size();
  // mmLimitTable
  mmLimitTable.assign_from_buffer(p, bufferSize - (p - dataBuffer));
  p += mmLimitTable.get_size();
  // hermite_conversion_info
  hermite_conversion_info.assign_from_buffer(p, bufferSize - (p - dataBuffer));
  p += hermite_conversion_info.get_size();
  // initialized
  memcpy(&initialized, p, sizeof(bool));
  p += sizeof(bool);
  // basis_func_poly_list
  memcpy(basis_func_poly_list, p, MAX_NO_OF_BASIS_FUNC_POLYS*sizeof(basis_func_poly_struct));
  p += MAX_NO_OF_BASIS_FUNC_POLYS*sizeof(basis_func_poly_struct);
  // no_of_basis_func_polys
  memcpy(&no_of_basis_func_polys, p, sizeof(int));
  p += sizeof(int);
  // monomial_info
  monomial_info.assign_from_buffer(p, bufferSize - (p - dataBuffer));
  p += monomial_info.get_size();
  // DONE!
  if(static_cast<size_t>(p-dataBuffer) > bufferSize)
    throw std::runtime_error("Error: (p > bufferSize).");  
}

