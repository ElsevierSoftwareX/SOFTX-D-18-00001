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

/** @file integral_info.h

    @brief Defines IntegralInfo object, providing the coefficients
    needed for integral evaluation.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef BASISINFO_BASIC_HEADER
#define BASISINFO_BASIC_HEADER


#include "realtype.h"
#include "monomial_info.h"
#include "hermite_conversion_prep.h"
#include "boysfunction.h"
#include "multipole_prep.h"
#include "mm_limit_table.h"

#ifndef BASIS_FUNC_POLY_MAX_DEGREE
#error The constant BASIS_FUNC_POLY_MAX_DEGREE must be defined.
#endif
#if BASIS_FUNC_POLY_MAX_DEGREE<6
#define MAX_NO_OF_TERMS_IN_BASIS_FUNC_POLY 12
#define MAX_NO_OF_POLY_12_TERMS 180
#define MAX_NO_OF_BASIS_FUNC_POLYS 50
#else
#define MAX_NO_OF_TERMS_IN_BASIS_FUNC_POLY 16
#define MAX_NO_OF_POLY_12_TERMS 360
#define MAX_NO_OF_BASIS_FUNC_POLYS 100
#endif

typedef struct
{
  ergo_real coeff;
  char monomialInts[4];  /* nx, ny, nz    */
  int monomialID;
} basis_func_term_struct;

typedef struct
{
  int noOfTerms;
  basis_func_term_struct termList[MAX_NO_OF_TERMS_IN_BASIS_FUNC_POLY];
  ergo_real scaledSolidHarmonicPrefactor;
} basis_func_poly_struct;

typedef struct
{
  int id_1;
  int id_2;
  ergo_real coeff;
} poly_12_term_struct;

typedef struct
{
  int noOfTerms;
  poly_12_term_struct termList[MAX_NO_OF_POLY_12_TERMS];
} poly_12_struct;

/** Contains coefficients needed for quick integral evaluation. This
    object is quite large and should always be allocated with
    new. Placing it on stack is a bad idea. */

class IntegralInfo
{
 private:
  BoysFunctionManager boysFunctionManager;
  MultipolePrepManager multipolePrep;
  MMLimitTable mmLimitTable;
  hermite_conversion_info_struct hermite_conversion_info;
  bool initialized;
  IntegralInfo(); // This is to make it forbidden to create it without argument.
 public:
  basis_func_poly_struct basis_func_poly_list[MAX_NO_OF_BASIS_FUNC_POLYS];
  int no_of_basis_func_polys;
  monomial_info_struct monomial_info;
  void init();
  ergo_real BoysFunction(int n, ergo_real x) const;
  ergo_real BoysFunction_expensive(int n, ergo_real x, int noOfIntegrationIntervals) const;
  const MultipolePrepManager & GetMultipolePrep() const { return multipolePrep; }
  const MMLimitTable & GetMMLimitTable() const { return mmLimitTable; }
  int multiply_by_hermite_conversion_matrix_from_right(int n1max,        
						       int n2max,        
						       ergo_real a,      
						       ergo_real* A,     
						       ergo_real* result) const;
  int multiply_by_hermite_conversion_matrix_from_left(int n1max,        
						      int n2max,        
						      ergo_real a,      
						      ergo_real* A,     
						      ergo_real* result) const;
  int get_hermite_conversion_matrix_right(int nmax,
					  ergo_real a,
					  ergo_real* result) const;
  int get_hermite_conversion_matrix_left(int nmax,
					 ergo_real a,
					 ergo_real* result) const;
  int get_hermite_conversion_matrix_right_sparse(int nmax,
                                                 ergo_real a,
                                                 i_j_val_struct* result) const;

  IntegralInfo(bool initialize);
  ~IntegralInfo();

  // Stuff needed for Chunks&Tasks usage
  IntegralInfo(const IntegralInfo & ii);
  void write_to_buffer ( char * dataBuffer, size_t const bufferSize ) const;
  size_t get_size() const;
  void assign_from_buffer ( char const * dataBuffer, size_t const bufferSize);
};


namespace JK {
/* Struct ExchWeights holds parameters for CAM-style range-separated HF
   exchange.  We use the following short-long range split
   (alpha+beta*erf(mu*r))*HF_exchange.
 */
struct ExchWeights
{
  ergo_real alpha;
  ergo_real beta;
  ergo_real mu;
  int computeRangeSeparatedExchange; /**< shortcut for |beta| != 0 */

ExchWeights() :
  alpha(0),
    beta(0),
    mu(0),
    computeRangeSeparatedExchange(0)
  {}
  
};

}


int get_poly_info_from_shell_type(int* polyid_start, int* poly_count, int shellType);

int get_no_of_basis_func_polys_used_from_no_of_shell_types(int no_of_shell_types);

int get_shell_type_from_basis_func_poly_id(int basfuncpolyid);


#endif
