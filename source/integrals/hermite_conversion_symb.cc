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

/** @file hermite_conversion_symb.cc

    @brief Code for conversion between integrals computed for Hermite
    Gaussians and Cartesian Gaussians, using a symbolic conversion
    matrix.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include "hermite_conversion_symb.h"


typedef struct
{
  int ix; // power of x
  int ia; // power of a
  ergo_real coeff;
} poly_1d_term_struct_symb;

#define MAX_NO_OF_1D_TERMS 888

typedef struct
{
  int noOfTerms;
  poly_1d_term_struct_symb termList[MAX_NO_OF_1D_TERMS];
} poly_1d_struct_symb;

typedef struct
{
  int monomialInts[3];
  int ia; // power of a
  ergo_real coeff;
} poly_3d_term_struct_symb;

#define MAX_NO_OF_3D_TERMS 888

typedef struct
{
  int noOfTerms;
  poly_3d_term_struct_symb termList[MAX_NO_OF_3D_TERMS];
} poly_3d_struct_symb;


static int
get_1d_hermite_poly_symb(poly_1d_struct_symb* result, int n)
{
  switch(n)
    {
    case 0:
      result->noOfTerms = 1;
      result->termList[0].ix = 0;
      result->termList[0].ia = 0;
      result->termList[0].coeff = 1;
      break;
    case 1:
      result->noOfTerms = 1;
      result->termList[0].ix = 1;
      result->termList[0].ia = 1;
      result->termList[0].coeff = 2;
      break;
    default:
      {
	// Create polys for n-1 and n-2
	poly_1d_struct_symb poly_n_m_1;
	poly_1d_struct_symb poly_n_m_2;
	get_1d_hermite_poly_symb(&poly_n_m_1, n - 1);
	get_1d_hermite_poly_symb(&poly_n_m_2, n - 2);
	assert(poly_n_m_1.noOfTerms + poly_n_m_2.noOfTerms < MAX_NO_OF_1D_TERMS);
	// Now the result is 2*a*x*poly_n_m_1 - (n-1)*2*a*poly_n_m_2
	for(int i = 0; i < poly_n_m_1.noOfTerms; i++)
	  {
	    result->termList[i] = poly_n_m_1.termList[i];
	    result->termList[i].ix++;
	    result->termList[i].ia++;
	    result->termList[i].coeff *= 2;
	  }
	int nn = poly_n_m_1.noOfTerms;
	for(int i = 0; i < poly_n_m_2.noOfTerms; i++)
	  {
	    result->termList[nn+i] = poly_n_m_2.termList[i];
	    result->termList[nn+i].ia++;
	    result->termList[nn+i].coeff *= -2 * (n-1);
	  }
	result->noOfTerms = poly_n_m_1.noOfTerms + poly_n_m_2.noOfTerms;
      }
    }
  return 0;
}

static int
get_1d_hermite_poly_inv_symb(poly_1d_struct_symb* result, int n)
{
  switch(n)
    {
    case 0:
      result->noOfTerms = 1;
      result->termList[0].ix = 0;
      result->termList[0].ia = 0;
      result->termList[0].coeff = 1;
      break;
    case 1:
      result->noOfTerms = 1;
      result->termList[0].ix = 1;
      result->termList[0].ia = -1;
      result->termList[0].coeff = 0.5;
      break;
    default:
      {
	// Create polys for n-1 and n-2
	poly_1d_struct_symb poly_n_m_1;
	poly_1d_struct_symb poly_n_m_2;
	get_1d_hermite_poly_inv_symb(&poly_n_m_1, n - 1);
	get_1d_hermite_poly_inv_symb(&poly_n_m_2, n - 2);
	assert(poly_n_m_1.noOfTerms + poly_n_m_2.noOfTerms < MAX_NO_OF_1D_TERMS);
	// Now the result is 0.5*(1/a)*x*poly_n_m_1 + (n-1)*(1/a)*0.5*poly_n_m_2
	for(int i = 0; i < poly_n_m_1.noOfTerms; i++)
	  {
	    result->termList[i] = poly_n_m_1.termList[i];
	    result->termList[i].ix++;
	    result->termList[i].ia--;
	    result->termList[i].coeff *= 0.5;
	  }
	int nn = poly_n_m_1.noOfTerms;
	for(int i = 0; i < poly_n_m_2.noOfTerms; i++)
	  {
	    result->termList[nn+i] = poly_n_m_2.termList[i];
	    result->termList[nn+i].ia--;
	    result->termList[nn+i].coeff *= (n-1) * 0.5;
	  }
	result->noOfTerms = poly_n_m_1.noOfTerms + poly_n_m_2.noOfTerms;
      }
    }
  return 0;
}


static int
create_3d_poly_from_1d_poly_symb(poly_3d_struct_symb* poly_3d, 
				 poly_1d_struct_symb* poly_1d, 
				 int coordIndex)
{
  memset(poly_3d, 0, sizeof(poly_3d_struct_symb));
  for(int i = 0; i < poly_1d->noOfTerms; i++)
    {
      poly_3d->termList[i].coeff = poly_1d->termList[i].coeff;
      poly_3d->termList[i].monomialInts[coordIndex] = poly_1d->termList[i].ix;
      poly_3d->termList[i].ia = poly_1d->termList[i].ia;
    }
  poly_3d->noOfTerms = poly_1d->noOfTerms;
  return 0;
}


static int
compute_product_of_3d_polys_symb(poly_3d_struct_symb* result, 
				 poly_3d_struct_symb* poly_1, 
				 poly_3d_struct_symb* poly_2)
{
  int termCount = 0;
  int termidx_1, termidx_2;
  for(termidx_1 = 0; termidx_1 < poly_1->noOfTerms; termidx_1++)
    for(termidx_2 = 0; termidx_2 < poly_2->noOfTerms; termidx_2++)
      {
	poly_3d_term_struct_symb* term_1 = &poly_1->termList[termidx_1];
	poly_3d_term_struct_symb* term_2 = &poly_2->termList[termidx_2];
	
	// Create product term
	poly_3d_term_struct_symb newTerm;
	newTerm.coeff = term_1->coeff * term_2->coeff;
	for(int k = 0; k < 3; k++)
	  newTerm.monomialInts[k] = 
	    term_1->monomialInts[k] + term_2->monomialInts[k];
	newTerm.ia = term_1->ia + term_2->ia;

	result->termList[termCount] = newTerm;
	termCount++;
	assert(termCount < MAX_NO_OF_3D_TERMS);
      } // END FOR termidx_1 termidx_2
  result->noOfTerms = termCount;
  return 0;
}


int 
get_hermite_conversion_matrix_symb(const monomial_info_struct* monomial_info,
				   int nmax, 
				   int inverseFlag,
				   symb_matrix_element* result)
{
  int noOfMonomials = monomial_info->no_of_monomials_list[nmax];
  memset(result, 0, noOfMonomials*noOfMonomials*sizeof(symb_matrix_element));

  int monomialIndex;
  for(monomialIndex = 0; monomialIndex < noOfMonomials; monomialIndex++)
    {

      // get monomialInts
      int ix = monomial_info->monomial_list[monomialIndex].ix;
      int iy = monomial_info->monomial_list[monomialIndex].iy;
      int iz = monomial_info->monomial_list[monomialIndex].iz;

      // Get x y z 1-d Hermite polynomials
      poly_1d_struct_symb hermitePoly_1d_x;
      poly_1d_struct_symb hermitePoly_1d_y;
      poly_1d_struct_symb hermitePoly_1d_z;
      if(inverseFlag == 1)
	{
	  get_1d_hermite_poly_inv_symb(&hermitePoly_1d_x, ix);
	  get_1d_hermite_poly_inv_symb(&hermitePoly_1d_y, iy);
	  get_1d_hermite_poly_inv_symb(&hermitePoly_1d_z, iz);
	}
      else
	{
	  get_1d_hermite_poly_symb(&hermitePoly_1d_x, ix);
	  get_1d_hermite_poly_symb(&hermitePoly_1d_y, iy);
	  get_1d_hermite_poly_symb(&hermitePoly_1d_z, iz);
	}

      // Store x y z Hermite polys as 3-d polys
      poly_3d_struct_symb hermitePoly_3d_x;
      poly_3d_struct_symb hermitePoly_3d_y;
      poly_3d_struct_symb hermitePoly_3d_z;
      create_3d_poly_from_1d_poly_symb(&hermitePoly_3d_x, &hermitePoly_1d_x, 0);
      create_3d_poly_from_1d_poly_symb(&hermitePoly_3d_y, &hermitePoly_1d_y, 1);
      create_3d_poly_from_1d_poly_symb(&hermitePoly_3d_z, &hermitePoly_1d_z, 2);

      // Compute product
      poly_3d_struct_symb hermitePoly_3d_xy;
      poly_3d_struct_symb hermitePoly_3d_xyz;
      compute_product_of_3d_polys_symb(&hermitePoly_3d_xy, 
				       &hermitePoly_3d_x, 
				       &hermitePoly_3d_y);
      compute_product_of_3d_polys_symb(&hermitePoly_3d_xyz, 
				       &hermitePoly_3d_xy, 
				       &hermitePoly_3d_z);
      
      // Go through result product poly, for each term get monomialIndex and add 
      // coeff to final result at position given by monomialIndex.
      for(int i = 0; i < hermitePoly_3d_xyz.noOfTerms; i++)
	{
	  poly_3d_term_struct_symb* currTerm = &hermitePoly_3d_xyz.termList[i];
	  // Get monomialIndex2
	  int ix = currTerm->monomialInts[0];
	  int iy = currTerm->monomialInts[1];
	  int iz = currTerm->monomialInts[2];
	  int monomialIndex2 = monomial_info->monomial_index_list[ix][iy][iz];
	  result[monomialIndex * noOfMonomials + monomialIndex2].coeff += currTerm->coeff;
	  if(result[monomialIndex * noOfMonomials + monomialIndex2].ia != 0)
	    assert(result[monomialIndex * noOfMonomials + monomialIndex2].ia == currTerm->ia);
	  result[monomialIndex * noOfMonomials + monomialIndex2].ia = currTerm->ia;
	} // END FOR i
      
    } // END FOR monomialIndex

  return 0;
}
