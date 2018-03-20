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

/** @file ci.cc

    @brief Configuration Interaction (CI) code.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <cstdlib>
#include <cstdio>
#include <memory.h>
#include <assert.h>
#include <vector>

#include "ci.h"
#include "output.h"
#include "utilities.h"
#include "integrals_2el_explicit.h"
#include "densfromf_full.h"
#include "simple_lanczos.h"
#include "operator_matrix.h"
#include "dipole_moment.h"

#include "../matrix/mat_gblas.h"


const int MAX_AOS = 30;
const int MAX_SOS = 2 * MAX_AOS;
const int MAX_ELECTRONS = 40;

const int SPIN_A = 1;
const int SPIN_B = 2;

typedef struct
{
  ergo_real x[MAX_AOS][MAX_AOS][MAX_AOS][MAX_AOS];
} four_idx_AO_struct;

typedef struct
{
  ergo_real x[MAX_SOS][MAX_SOS][MAX_SOS][MAX_SOS];
} four_idx_SO_struct;

typedef struct
{
  ergo_real x[MAX_SOS][MAX_SOS];
} two_idx_SO_struct;

typedef struct
{
  ergo_real coeffs[MAX_AOS];
  int spin;
} SO_struct;

typedef struct
{
  char SO_list[MAX_ELECTRONS];
  int startIndex; // at next level
  int count; // at next level
} SlaterDet_struct;

typedef struct
{
  int a;
  int b;
  int nDiff;
  char SOs_a[2];
  char SOs_b[2];
  char SOs_a_pos[2];
  char SOs_b_pos[2];
} SlaterDet_pair_struct;


static ergo_real get_vector_norm(int n, const ergo_real* v)
{
  ergo_real sqSum = 0;
  for(int i = 0; i < n; i++)
    sqSum += v[i] * v[i];
  return template_blas_sqrt(sqSum);
}

static void normalize_vector(int n, ergo_real* v)
{
  ergo_real factor = 1.0 / get_vector_norm(n, v);
  for(int i = 0; i < n; i++)
    v[i] *= factor;
}


void get_1el_energy_and_gradient(int nSOs,
				 int nEl,
				 int nSlaterDets, 
				 const SlaterDet_struct* SlaterDetList, 
				 int nSlaterDetPairs,
				 const SlaterDet_pair_struct* SlaterDet_pair_list,
				 const int* pairCountList,
				 int noOfTrialVectors,
				 ergo_real* energy_list,
				 ergo_real** coeffListList,
				 const two_idx_SO_struct* h_SO, 
				 ergo_real** resultGradient_list)
{
  for(int k = 0; k < noOfTrialVectors; k++)
    {
      energy_list[k] = 0;
      for(int i = 0; i < nSlaterDets; i++)
	resultGradient_list[k][i] = 0;
    }
  for(int pairIdx = 0; pairIdx < nSlaterDetPairs; pairIdx++)
    {      

      if(SlaterDet_pair_list[pairIdx].nDiff == 2)
	{
	  // Do nothing here!
	}
      else if(SlaterDet_pair_list[pairIdx].nDiff == 1)
	{

	  int a = SlaterDet_pair_list[pairIdx].a;
	  int b = SlaterDet_pair_list[pairIdx].b;

	  int p = SlaterDet_pair_list[pairIdx].SOs_a[0];

	  int count, savedCount;
	  int signa = 1;
	  // annihilate p from a
	  char SO_list_mod_a_1[MAX_ELECTRONS];
	  count = 0;
	  savedCount = -1;
	  for(int i = 0; i < nEl; i++)
	    {
	      if(SlaterDetList[a].SO_list[i] == p)
		savedCount = count;
	      else
		{
		  SO_list_mod_a_1[count] = SlaterDetList[a].SO_list[i];
		  count++;
		  if(savedCount < 0)
		    signa *= -1;
		}
	    }
	  if(savedCount < 0)
	    continue;

	  int q = SlaterDet_pair_list[pairIdx].SOs_b[0];

	  int signb = signa;
	  // annihilate q from b
	  char SO_list_mod_b_1[MAX_ELECTRONS];
	  count = 0;
	  savedCount = -1;
	  for(int i = 0; i < nEl; i++)
	    {
	      if(SlaterDetList[b].SO_list[i] == q)
		savedCount = count;
	      else
		{
		  SO_list_mod_b_1[count] = SlaterDetList[b].SO_list[i];
		  count++;
		  if(savedCount < 0)
		    signb *= -1;
		}
	    }
	  if(savedCount < 0)
	    continue;
		  
	  int equal = 1;
	  for(int i = 0; i < nEl-1; i++)
	    {
	      if(SO_list_mod_a_1[i] != SO_list_mod_b_1[i])
		equal = 0;
	    }
		  
	  if(equal)
	    {
	      // OK, we have a contribution
	      ergo_real x = signb;
	      x *= 2; // To account for a > b case
	      ergo_real h_pq = h_SO->x[p][q];
	      for(int k = 0; k < noOfTrialVectors; k++)
		{
		  energy_list[k] += x * h_pq * coeffListList[k][a] * coeffListList[k][b];
		  resultGradient_list[k][a] += x * h_pq * coeffListList[k][b];
		  resultGradient_list[k][b] += x * h_pq * coeffListList[k][a];
		}
	    }
	}
      else
	{
	  
	  int a = SlaterDet_pair_list[pairIdx].a;
	  int b = SlaterDet_pair_list[pairIdx].b;

	  for(int p = 0; p < nSOs; p++)
	    {
	      int count, savedCount;
	      int signa = 1;
	      // annihilate p from a
	      char SO_list_mod_a_1[MAX_ELECTRONS];
	      count = 0;
	      savedCount = -1;
	      for(int i = 0; i < nEl; i++)
		{
		  if(SlaterDetList[a].SO_list[i] == p)
		    savedCount = count;
		  else
		    {
		      SO_list_mod_a_1[count] = SlaterDetList[a].SO_list[i];
		      count++;
		      if(savedCount < 0)
			signa *= -1;
		    }
		}
	      if(savedCount < 0)
		continue;
	      for(int q = 0; q < nSOs; q++)
		{
		  int signb = signa;
		  // annihilate q from b
		  char SO_list_mod_b_1[MAX_ELECTRONS];
		  count = 0;
		  savedCount = -1;
		  for(int i = 0; i < nEl; i++)
		    {
		      if(SlaterDetList[b].SO_list[i] == q)
			savedCount = count;
		      else
			{
			  SO_list_mod_b_1[count] = SlaterDetList[b].SO_list[i];
			  count++;
			  if(savedCount < 0)
			    signb *= -1;
			}
		    }
		  if(savedCount < 0)
		    continue;
		  
		  int equal = 1;
		  for(int i = 0; i < nEl-1; i++)
		    {
		      if(SO_list_mod_a_1[i] != SO_list_mod_b_1[i])
			equal = 0;
		    }
		  
		  if(equal)
		    {
		      // OK, we have a contribution
		      ergo_real x = signb;
		      ergo_real h_pq = h_SO->x[p][q];
		      for(int k = 0; k < noOfTrialVectors; k++)
			{
			  energy_list[k] += x * h_pq * coeffListList[k][a] * coeffListList[k][b];
			  resultGradient_list[k][a] += x * h_pq * coeffListList[k][b];
			  resultGradient_list[k][b] += x * h_pq * coeffListList[k][a];
			}
		    }
		} // END FOR q
	    } // END FOR p
	} // END ELSE
    } // END FOR pairIdx
}




void get_1el_contribs_to_mult_or_dmat(int nSOs,
				      int nEl,
				      int nSlaterDets, 
				      const SlaterDet_struct* SlaterDetList, 
				      const SlaterDet_pair_struct* SlaterDetPair,
				      const two_idx_SO_struct* h_SO, 
				      const ergo_real* sourceVector,
				      ergo_real* resultVector, // if result of matrix-vector mult is requested
				      two_idx_SO_struct* resultdmat // if dmat is requested
				      )
{
  if(SlaterDetPair->nDiff == 2)
    {
      // Do nothing here!
    }
  else if(SlaterDetPair->nDiff == 1)
    {

      int a = SlaterDetPair->a;
      int b = SlaterDetPair->b;

      int p = SlaterDetPair->SOs_a[0];

      int count, savedCount;
      int signa = 1;
      // annihilate p from a
      char SO_list_mod_a_1[MAX_ELECTRONS];
      count = 0;
      savedCount = -1;
      for(int i = 0; i < nEl; i++)
	{
	  if(SlaterDetList[a].SO_list[i] == p)
	    savedCount = count;
	  else
	    {
	      SO_list_mod_a_1[count] = SlaterDetList[a].SO_list[i];
	      count++;
	      if(savedCount < 0)
		signa *= -1;
	    }
	}
      if(savedCount < 0)
	return;

      int q = SlaterDetPair->SOs_b[0];

      int signb = signa;
      // annihilate q from b
      char SO_list_mod_b_1[MAX_ELECTRONS];
      count = 0;
      savedCount = -1;
      for(int i = 0; i < nEl; i++)
	{
	  if(SlaterDetList[b].SO_list[i] == q)
	    savedCount = count;
	  else
	    {
	      SO_list_mod_b_1[count] = SlaterDetList[b].SO_list[i];
	      count++;
	      if(savedCount < 0)
		signb *= -1;
	    }
	}
      if(savedCount < 0)
	return;
		  
      int equal = 1;
      for(int i = 0; i < nEl-1; i++)
	{
	  if(SO_list_mod_a_1[i] != SO_list_mod_b_1[i])
	    equal = 0;
	}
		  
      if(equal)
	{
	  // OK, we have a contribution
	  ergo_real x = signb;
	  ergo_real h_pq = h_SO->x[p][q];
	  
	  ergo_real resultMatrix_ab = x * h_pq;
	  ergo_real resultMatrix_ba = x * h_pq;

	  if(resultVector)
	    {
	      resultVector[a] += resultMatrix_ab * sourceVector[b];
	      resultVector[b] += resultMatrix_ba * sourceVector[a];
	    }
	  else
	    {
	      resultdmat->x[p][q] += x * sourceVector[a] * sourceVector[b];
	      resultdmat->x[q][p] += x * sourceVector[a] * sourceVector[b];
	    }
	}
    }
  else
    {
	  
      int a = SlaterDetPair->a;
      int b = SlaterDetPair->b;

      for(int p = 0; p < nSOs; p++)
	{
	  int count, savedCount;
	  int signa = 1;
	  // annihilate p from a
	  char SO_list_mod_a_1[MAX_ELECTRONS];
	  count = 0;
	  savedCount = -1;
	  for(int i = 0; i < nEl; i++)
	    {
	      if(SlaterDetList[a].SO_list[i] == p)
		savedCount = count;
	      else
		{
		  SO_list_mod_a_1[count] = SlaterDetList[a].SO_list[i];
		  count++;
		  if(savedCount < 0)
		    signa *= -1;
		}
	    }
	  if(savedCount < 0)
	    continue;
	  for(int q = 0; q < nSOs; q++)
	    {
	      int signb = signa;
	      // annihilate q from b
	      char SO_list_mod_b_1[MAX_ELECTRONS];
	      count = 0;
	      savedCount = -1;
	      for(int i = 0; i < nEl; i++)
		{
		  if(SlaterDetList[b].SO_list[i] == q)
		    savedCount = count;
		  else
		    {
		      SO_list_mod_b_1[count] = SlaterDetList[b].SO_list[i];
		      count++;
		      if(savedCount < 0)
			signb *= -1;
		    }
		}
	      if(savedCount < 0)
		continue;
		  
	      int equal = 1;
	      for(int i = 0; i < nEl-1; i++)
		{
		  if(SO_list_mod_a_1[i] != SO_list_mod_b_1[i])
		    equal = 0;
		}
		  
	      if(equal)
		{
		  // OK, we have a contribution
		  ergo_real x = signb;
		  ergo_real h_pq = h_SO->x[p][q];

		  ergo_real resultMatrix_ab = x * h_pq;

		  if(resultVector)
		    resultVector[a] += resultMatrix_ab * sourceVector[b];
		  else
		    resultdmat->x[p][q] += x * sourceVector[a] * sourceVector[b];
		}
	    } // END FOR q
	} // END FOR p
    } // END ELSE
}






void get_1el_contribs(int nSOs,
		      int nEl,
		      int nSlaterDets, 
		      const SlaterDet_struct* SlaterDetList, 
		      int nSlaterDetPairs,
		      const SlaterDet_pair_struct* SlaterDet_pair_list,
		      const int* pairCountList,
		      const two_idx_SO_struct* h_SO, 
		      ergo_real* resultMatrix)
{
  for(int pairIdx = 0; pairIdx < nSlaterDetPairs; pairIdx++)
    {            
      if(SlaterDet_pair_list[pairIdx].nDiff == 2)
	{
	  // Do nothing here!
	}
      else if(SlaterDet_pair_list[pairIdx].nDiff == 1)
	{

	  int a = SlaterDet_pair_list[pairIdx].a;
	  int b = SlaterDet_pair_list[pairIdx].b;

	  int p = SlaterDet_pair_list[pairIdx].SOs_a[0];

	  int count, savedCount;
	  int signa = 1;
	  // annihilate p from a
	  char SO_list_mod_a_1[MAX_ELECTRONS];
	  count = 0;
	  savedCount = -1;
	  for(int i = 0; i < nEl; i++)
	    {
	      if(SlaterDetList[a].SO_list[i] == p)
		savedCount = count;
	      else
		{
		  SO_list_mod_a_1[count] = SlaterDetList[a].SO_list[i];
		  count++;
		  if(savedCount < 0)
		    signa *= -1;
		}
	    }
	  if(savedCount < 0)
	    continue;

	  int q = SlaterDet_pair_list[pairIdx].SOs_b[0];

	  int signb = signa;
	  // annihilate q from b
	  char SO_list_mod_b_1[MAX_ELECTRONS];
	  count = 0;
	  savedCount = -1;
	  for(int i = 0; i < nEl; i++)
	    {
	      if(SlaterDetList[b].SO_list[i] == q)
		savedCount = count;
	      else
		{
		  SO_list_mod_b_1[count] = SlaterDetList[b].SO_list[i];
		  count++;
		  if(savedCount < 0)
		    signb *= -1;
		}
	    }
	  if(savedCount < 0)
	    continue;
		  
	  int equal = 1;
	  for(int i = 0; i < nEl-1; i++)
	    {
	      if(SO_list_mod_a_1[i] != SO_list_mod_b_1[i])
		equal = 0;
	    }
		  
	  if(equal)
	    {
	      // OK, we have a contribution
	      ergo_real x = signb;
	      //x *= 2; // To account for a > b case
	      ergo_real h_pq = h_SO->x[p][q];

	      //printf("1el contrib to a b = %i %i : %22.11f\n", a, b, x * h_pq);

	      resultMatrix[a*nSlaterDets+b] += x * h_pq;
	      resultMatrix[b*nSlaterDets+a] += x * h_pq;
	    }
	}
      else
	{
	  
	  int a = SlaterDet_pair_list[pairIdx].a;
	  int b = SlaterDet_pair_list[pairIdx].b;

	  for(int p = 0; p < nSOs; p++)
	    {
	      int count, savedCount;
	      int signa = 1;
	      // annihilate p from a
	      char SO_list_mod_a_1[MAX_ELECTRONS];
	      count = 0;
	      savedCount = -1;
	      for(int i = 0; i < nEl; i++)
		{
		  if(SlaterDetList[a].SO_list[i] == p)
		    savedCount = count;
		  else
		    {
		      SO_list_mod_a_1[count] = SlaterDetList[a].SO_list[i];
		      count++;
		      if(savedCount < 0)
			signa *= -1;
		    }
		}
	      if(savedCount < 0)
		continue;
	      for(int q = 0; q < nSOs; q++)
		{
		  int signb = signa;
		  // annihilate q from b
		  char SO_list_mod_b_1[MAX_ELECTRONS];
		  count = 0;
		  savedCount = -1;
		  for(int i = 0; i < nEl; i++)
		    {
		      if(SlaterDetList[b].SO_list[i] == q)
			savedCount = count;
		      else
			{
			  SO_list_mod_b_1[count] = SlaterDetList[b].SO_list[i];
			  count++;
			  if(savedCount < 0)
			    signb *= -1;
			}
		    }
		  if(savedCount < 0)
		    continue;
		  
		  int equal = 1;
		  for(int i = 0; i < nEl-1; i++)
		    {
		      if(SO_list_mod_a_1[i] != SO_list_mod_b_1[i])
			equal = 0;
		    }
		  
		  if(equal)
		    {
		      // OK, we have a contribution
		      ergo_real x = signb;
		      ergo_real h_pq = h_SO->x[p][q];

		      //printf("1el contrib to a b = %i %i : %22.11f\n", a, b, x * h_pq);

		      resultMatrix[a*nSlaterDets+b] += x * h_pq;
		    }
		} // END FOR q
	    } // END FOR p
	} // END ELSE
    } // END FOR pairIdx
}






typedef struct
{
  int p;
  int q;
  int r;
  int s;
  int sign;
} contrib_debug_struct;



void get_2el_energy_and_gradient(int nSOs,
				 int nEl,
				 int nSlaterDets, 
				 const SlaterDet_struct* SlaterDetList, 
				 int nSlaterDetPairs,
				 const SlaterDet_pair_struct* SlaterDet_pair_list,
				 const int* pairCountList,
				 int noOfTrialVectors,
				 ergo_real* energy_list,
				 ergo_real** coeffListList,
				 const four_idx_SO_struct* g_SO, 
				 ergo_real** resultGradient_list)
{
  for(int k = 0; k < noOfTrialVectors; k++)
    {
      energy_list[k] = 0;
      for(int i = 0; i < nSlaterDets; i++)
	resultGradient_list[k][i] = 0;
    }
  for(int pairIdx = 0; pairIdx < nSlaterDetPairs; pairIdx++)
    {      
      int a = SlaterDet_pair_list[pairIdx].a;
      int b = SlaterDet_pair_list[pairIdx].b;

      if(SlaterDet_pair_list[pairIdx].nDiff == 2)
	{
	  int p = SlaterDet_pair_list[pairIdx].SOs_a[0];
	  int q = SlaterDet_pair_list[pairIdx].SOs_b[0];
	  int pos_p = SlaterDet_pair_list[pairIdx].SOs_a_pos[0];
	  int pos_q = SlaterDet_pair_list[pairIdx].SOs_b_pos[0];
	  
	  int r = SlaterDet_pair_list[pairIdx].SOs_a[1];
	  int s = SlaterDet_pair_list[pairIdx].SOs_b[1];
	  int pos_r = SlaterDet_pair_list[pairIdx].SOs_a_pos[1];
	  int pos_s = SlaterDet_pair_list[pairIdx].SOs_b_pos[1];

	  int sign_a = 1;
	  if(pos_p%2 == 1)
	    sign_a *= -1;
	  if(pos_r%2 == 1)
	    sign_a *= -1;
	  if(r > p)
	    sign_a *= -1;
	  
	  int sign_b = 1;
	  if(pos_q%2 == 1)
	    sign_b *= -1;
	  if(pos_s%2 == 1)
	    sign_b *= -1;
	  if(s > q)
	    sign_a *= -1;

	  // OK, we have a contribution
	  ergo_real x = sign_a * sign_b;
	  x *= 2; // To account for a > b case

	  for(int k = 0; k < noOfTrialVectors; k++)
	    {
	      // pqrs
	      ergo_real g_pqrs = g_SO->x[p][q][r][s];
	      energy_list[k]            += 0.5 * x * g_pqrs * coeffListList[k][a] * coeffListList[k][b];
	      resultGradient_list[k][a] += 0.5 * x * g_pqrs * coeffListList[k][b];
	      resultGradient_list[k][b] += 0.5 * x * g_pqrs * coeffListList[k][a];

	      // rspq
	      ergo_real g_rspq = g_SO->x[r][s][p][q];
	      energy_list[k]            += 0.5 * x * g_rspq * coeffListList[k][a] * coeffListList[k][b];
	      resultGradient_list[k][a] += 0.5 * x * g_rspq * coeffListList[k][b];
	      resultGradient_list[k][b] += 0.5 * x * g_rspq * coeffListList[k][a];

	      // rqps
	      ergo_real g_rqps = g_SO->x[r][q][p][s];
	      energy_list[k]            -= 0.5 * x * g_rqps * coeffListList[k][a] * coeffListList[k][b];
	      resultGradient_list[k][a] -= 0.5 * x * g_rqps * coeffListList[k][b];
	      resultGradient_list[k][b] -= 0.5 * x * g_rqps * coeffListList[k][a];

	      // psrq
	      ergo_real g_psrq = g_SO->x[p][s][r][q];
	      energy_list[k]            -= 0.5 * x * g_psrq * coeffListList[k][a] * coeffListList[k][b];
	      resultGradient_list[k][a] -= 0.5 * x * g_psrq * coeffListList[k][b];
	      resultGradient_list[k][b] -= 0.5 * x * g_psrq * coeffListList[k][a];
	    } // END FOR k
	}
      else if(SlaterDet_pair_list[pairIdx].nDiff == 1)
	{
	  int p = SlaterDet_pair_list[pairIdx].SOs_a[0];
	  int q = SlaterDet_pair_list[pairIdx].SOs_b[0];
	  int pos_p = SlaterDet_pair_list[pairIdx].SOs_a_pos[0];
	  int pos_q = SlaterDet_pair_list[pairIdx].SOs_b_pos[0];

	  for(int rr = 0; rr < nEl; rr++)
	    {
	      int r = SlaterDetList[a].SO_list[rr];
	      if(r == p)
		continue;

	      int sign_a = 1;
	      if(pos_p%2 == 1)
		sign_a *= -1;		
	      if(rr%2 == 1)
		sign_a *= -1;
	      if(r > p)
		sign_a *= -1;
	      
	      for(int ss = 0; ss < nEl; ss++)
		{
		  int s = SlaterDetList[b].SO_list[ss];
		  if(s == q)
		    continue;
		  if(s != r)
		    continue;

		  int sign_b = 1;
		  if(pos_q%2 == 1)
		    sign_b *= -1;
		  if(ss%2 == 1)
		    sign_b *= -1;
		  if(s > q)
		    sign_b *= -1;

		  // OK, we have a contribution
		  ergo_real x = sign_a * sign_b;
		  x *= 2; // To account for a > b case

		  for(int k = 0; k < noOfTrialVectors; k++)
		    {
		      // pqrs
		      ergo_real g_pqrs = g_SO->x[p][q][r][s];
		      energy_list[k]            += 0.5 * x * g_pqrs * coeffListList[k][a] * coeffListList[k][b];
		      resultGradient_list[k][a] += 0.5 * x * g_pqrs * coeffListList[k][b];
		      resultGradient_list[k][b] += 0.5 * x * g_pqrs * coeffListList[k][a];

		      // rspq
		      ergo_real g_rspq = g_SO->x[r][s][p][q];
		      energy_list[k]            += 0.5 * x * g_rspq * coeffListList[k][a] * coeffListList[k][b];
		      resultGradient_list[k][a] += 0.5 * x * g_rspq * coeffListList[k][b];
		      resultGradient_list[k][b] += 0.5 * x * g_rspq * coeffListList[k][a];

		      // rqps
		      ergo_real g_rqps = g_SO->x[r][q][p][s];
		      energy_list[k]            -= 0.5 * x * g_rqps * coeffListList[k][a] * coeffListList[k][b];
		      resultGradient_list[k][a] -= 0.5 * x * g_rqps * coeffListList[k][b];
		      resultGradient_list[k][b] -= 0.5 * x * g_rqps * coeffListList[k][a];

		      // psrq
		      ergo_real g_psrq = g_SO->x[p][s][r][q];
		      energy_list[k]            -= 0.5 * x * g_psrq * coeffListList[k][a] * coeffListList[k][b];
		      resultGradient_list[k][a] -= 0.5 * x * g_psrq * coeffListList[k][b];
		      resultGradient_list[k][b] -= 0.5 * x * g_psrq * coeffListList[k][a];
		    } // END FOR k

		} // END FOR ss
	    } // END FOR rr
	}
      else
	{
	  for(int pp = 0; pp < nEl; pp++)
	    {
	      int p = SlaterDetList[a].SO_list[pp];
	      for(int rr = 0; rr < nEl; rr++)
		{
		  int r = SlaterDetList[a].SO_list[rr];
		  int sign_a = 1;
		  if(pp%2 == 1)
		    sign_a *= -1;
		  if(rr%2 == 1)
		    sign_a *= -1;
		  if(pp < rr)
		    sign_a *= -1;		    
		  int qlist[2];
		  qlist[0] = pp;
		  qlist[1] = rr;
		  for(int qqq = 0; qqq < 2; qqq++)
		    {
		      int qq = qlist[qqq];
		      int q = SlaterDetList[b].SO_list[qq];
		      int slist[2];
		      slist[0] = pp;
		      slist[1] = rr;
		      for(int sss = 0; sss < 2; sss++)
			{
			  int ss = slist[sss];
			  int s = SlaterDetList[b].SO_list[ss];
			  if(s == q)
			    continue;
			  int sign_b = 1;
			  if(qq%2 == 1)
			    sign_b *= -1;
			  if(ss%2 == 1)
			    sign_b *= -1;
			  if(qq < ss)
			    sign_b *= -1;

			  // OK, we have a contribution
			  ergo_real x = sign_a * sign_b;
			  ergo_real g_pqrs = g_SO->x[p][q][r][s];

			  for(int k = 0; k < noOfTrialVectors; k++)
			    {
			      energy_list[k]            += 0.5 * x * g_pqrs * coeffListList[k][a] * coeffListList[k][b];
			      resultGradient_list[k][a] += 0.5 * x * g_pqrs * coeffListList[k][b];
			      resultGradient_list[k][b] += 0.5 * x * g_pqrs * coeffListList[k][a];
			    }
			  
			} // END FOR sss
		    } // END FOR qqq
		} // END FOR rr
	    } // END FOR pp

	} // END ELSE

    } // END FOR pair
  
}





void get_2el_contribs_to_mult_or_dmat(int nSOs,
				      int nEl,
				      int nSlaterDets, 
				      const SlaterDet_struct* SlaterDetList, 
				      const SlaterDet_pair_struct* SlaterDetPair,
				      const four_idx_SO_struct* g_SO, 
				      const ergo_real* sourceVector,
				      ergo_real* resultVector, // if result of matrix-vector mult is requested
				      four_idx_SO_struct* result_dmat_2el // if dmat is requested
				      )
{
  int a = SlaterDetPair->a;
  int b = SlaterDetPair->b;

  if(SlaterDetPair->nDiff == 2)
    {
      int p = SlaterDetPair->SOs_a[0];
      int q = SlaterDetPair->SOs_b[0];
      int pos_p = SlaterDetPair->SOs_a_pos[0];
      int pos_q = SlaterDetPair->SOs_b_pos[0];
      
      int r = SlaterDetPair->SOs_a[1];
      int s = SlaterDetPair->SOs_b[1];
      int pos_r = SlaterDetPair->SOs_a_pos[1];
      int pos_s = SlaterDetPair->SOs_b_pos[1];
      
      int sign_a = 1;
      if(pos_p%2 == 1)
	sign_a *= -1;
      if(pos_r%2 == 1)
	sign_a *= -1;
      if(r > p)
	sign_a *= -1;
      
      int sign_b = 1;
      if(pos_q%2 == 1)
	sign_b *= -1;
      if(pos_s%2 == 1)
	sign_b *= -1;
      if(s > q)
	sign_a *= -1;

      // OK, we have a contribution
      ergo_real x = sign_a * sign_b;

      ergo_real g_pqrs = g_SO->x[p][q][r][s];
      ergo_real g_rspq = g_SO->x[r][s][p][q];
      ergo_real g_rqps = g_SO->x[r][q][p][s];
      ergo_real g_psrq = g_SO->x[p][s][r][q];

      ergo_real resultMatrix_ab = 0;
      resultMatrix_ab += 0.5 * x * g_pqrs;
      resultMatrix_ab += 0.5 * x * g_rspq;
      resultMatrix_ab -= 0.5 * x * g_rqps;
      resultMatrix_ab -= 0.5 * x * g_psrq;
      ergo_real resultMatrix_ba = resultMatrix_ab;


      if(resultVector)
	{
	  resultVector[a] += resultMatrix_ab * sourceVector[b];
	  resultVector[b] += resultMatrix_ba * sourceVector[a];
	}
      else
	{
	  ergo_real value = x * sourceVector[a] * sourceVector[b];
	  result_dmat_2el->x[p][q][r][s] += value;
	  result_dmat_2el->x[q][p][s][r] += value;
	  result_dmat_2el->x[r][s][p][q] += value;
	  result_dmat_2el->x[s][r][q][p] += value;
	  result_dmat_2el->x[r][q][p][s] -= value;
	  result_dmat_2el->x[q][r][s][p] -= value;
	  result_dmat_2el->x[p][s][r][q] -= value;
	  result_dmat_2el->x[s][p][q][r] -= value;
	}
    }
  else if(SlaterDetPair->nDiff == 1)
    {
      int p = SlaterDetPair->SOs_a[0];
      int q = SlaterDetPair->SOs_b[0];
      int pos_p = SlaterDetPair->SOs_a_pos[0];
      int pos_q = SlaterDetPair->SOs_b_pos[0];
      
      for(int rr = 0; rr < nEl; rr++)
	{
	  int r = SlaterDetList[a].SO_list[rr];
	  if(r == p)
	    continue;

	  int sign_a = 1;
	  if(pos_p%2 == 1)
	    sign_a *= -1;		
	  if(rr%2 == 1)
	    sign_a *= -1;
	  if(r > p)
	    sign_a *= -1;
	  
	  for(int ss = 0; ss < nEl; ss++)
	    {
	      int s = SlaterDetList[b].SO_list[ss];
	      if(s == q)
		continue;
	      if(s != r)
		continue;

	      int sign_b = 1;
	      if(pos_q%2 == 1)
		sign_b *= -1;
	      if(ss%2 == 1)
		sign_b *= -1;
	      if(s > q)
		sign_b *= -1;

	      // OK, we have a contribution
	      ergo_real x = sign_a * sign_b;

	      ergo_real g_pqrs = g_SO->x[p][q][r][s];
	      ergo_real g_rspq = g_SO->x[r][s][p][q];
	      ergo_real g_rqps = g_SO->x[r][q][p][s];
	      ergo_real g_psrq = g_SO->x[p][s][r][q];

	      ergo_real resultMatrix_ab = 0;
	      resultMatrix_ab += 0.5 * x * g_pqrs;
	      resultMatrix_ab += 0.5 * x * g_rspq;
	      resultMatrix_ab -= 0.5 * x * g_rqps;
	      resultMatrix_ab -= 0.5 * x * g_psrq;
	      ergo_real resultMatrix_ba = resultMatrix_ab;

	      if(resultVector)
		{
		  resultVector[a] += resultMatrix_ab * sourceVector[b];
		  resultVector[b] += resultMatrix_ba * sourceVector[a];
		}
	      else
		{
		  ergo_real value = x * sourceVector[a] * sourceVector[b];
		  result_dmat_2el->x[p][q][r][s] += value;
		  result_dmat_2el->x[q][p][r][s] += value;
		  result_dmat_2el->x[r][s][p][q] += value;
		  result_dmat_2el->x[r][s][q][p] += value;
		  result_dmat_2el->x[r][q][p][s] -= value;
		  result_dmat_2el->x[q][r][s][p] -= value;
		  result_dmat_2el->x[p][s][r][q] -= value;
		  result_dmat_2el->x[s][p][q][r] -= value;
		}
	    } // END FOR ss
	} // END FOR rr
    }
  else
    {
      for(int pp = 0; pp < nEl; pp++)
	{
	  int p = SlaterDetList[a].SO_list[pp];
	  for(int rr = 0; rr < nEl; rr++)
	    {
	      int r = SlaterDetList[a].SO_list[rr];
	      int sign_a = 1;
	      if(pp%2 == 1)
		sign_a *= -1;
	      if(rr%2 == 1)
		sign_a *= -1;
	      if(pp < rr)
		sign_a *= -1;		    
	      int qlist[2];
	      qlist[0] = pp;
	      qlist[1] = rr;
	      for(int qqq = 0; qqq < 2; qqq++)
		{
		  int qq = qlist[qqq];
		  int q = SlaterDetList[b].SO_list[qq];
		  int slist[2];
		  slist[0] = pp;
		  slist[1] = rr;
		  for(int sss = 0; sss < 2; sss++)
		    {
		      int ss = slist[sss];
		      int s = SlaterDetList[b].SO_list[ss];
		      if(s == q)
			continue;
		      int sign_b = 1;
		      if(qq%2 == 1)
			sign_b *= -1;
		      if(ss%2 == 1)
			sign_b *= -1;
		      if(qq < ss)
			sign_b *= -1;

		      // OK, we have a contribution
		      ergo_real x = sign_a * sign_b;
		      ergo_real g_pqrs = g_SO->x[p][q][r][s];

		      ergo_real resultMatrix_ab = 0;
		      resultMatrix_ab += 0.5 * x * g_pqrs;

		      if(resultVector)
			{
			  resultVector[a] += resultMatrix_ab * sourceVector[b];
			}
		      else
			{
			  ergo_real value = x * sourceVector[a] * sourceVector[b];
			  result_dmat_2el->x[p][q][r][s] += value;
			}
		    } // END FOR sss
		} // END FOR qqq
	    } // END FOR rr
	} // END FOR pp
      
    } // END ELSE

}






void get_2el_contribs(int nSOs,
		      int nEl,
		      int nSlaterDets, 
		      const SlaterDet_struct* SlaterDetList, 
		      int nSlaterDetPairs,
		      const SlaterDet_pair_struct* SlaterDet_pair_list,
		      const int* pairCountList,
		      const four_idx_SO_struct* g_SO, 
		      ergo_real* resultMatrix)
{
  for(int pairIdx = 0; pairIdx < nSlaterDetPairs; pairIdx++)
    {      
      int a = SlaterDet_pair_list[pairIdx].a;
      int b = SlaterDet_pair_list[pairIdx].b;

      if(SlaterDet_pair_list[pairIdx].nDiff == 2)
	{
	  int p = SlaterDet_pair_list[pairIdx].SOs_a[0];
	  int q = SlaterDet_pair_list[pairIdx].SOs_b[0];
	  int pos_p = SlaterDet_pair_list[pairIdx].SOs_a_pos[0];
	  int pos_q = SlaterDet_pair_list[pairIdx].SOs_b_pos[0];
	  
	  int r = SlaterDet_pair_list[pairIdx].SOs_a[1];
	  int s = SlaterDet_pair_list[pairIdx].SOs_b[1];
	  int pos_r = SlaterDet_pair_list[pairIdx].SOs_a_pos[1];
	  int pos_s = SlaterDet_pair_list[pairIdx].SOs_b_pos[1];

	  int sign_a = 1;
	  if(pos_p%2 == 1)
	    sign_a *= -1;
	  if(pos_r%2 == 1)
	    sign_a *= -1;
	  if(r > p)
	    sign_a *= -1;
	  
	  int sign_b = 1;
	  if(pos_q%2 == 1)
	    sign_b *= -1;
	  if(pos_s%2 == 1)
	    sign_b *= -1;
	  if(s > q)
	    sign_a *= -1;

	  // OK, we have a contribution
	  ergo_real x = sign_a * sign_b;
	  //x *= 2; // To account for a > b case

	  ergo_real g_pqrs = g_SO->x[p][q][r][s];
	  ergo_real g_rspq = g_SO->x[r][s][p][q];
	  ergo_real g_rqps = g_SO->x[r][q][p][s];
	  ergo_real g_psrq = g_SO->x[p][s][r][q];

	  //printf("2el contrib to a b = %i %i : %22.11f\n", a, b, 0.5 * x * g_pqrs);

	  resultMatrix[a*nSlaterDets+b] += 0.5 * x * g_pqrs;
	  resultMatrix[a*nSlaterDets+b] += 0.5 * x * g_rspq;
	  resultMatrix[a*nSlaterDets+b] -= 0.5 * x * g_rqps;
	  resultMatrix[a*nSlaterDets+b] -= 0.5 * x * g_psrq;

	  resultMatrix[b*nSlaterDets+a] += 0.5 * x * g_pqrs;
	  resultMatrix[b*nSlaterDets+a] += 0.5 * x * g_rspq;
	  resultMatrix[b*nSlaterDets+a] -= 0.5 * x * g_rqps;
	  resultMatrix[b*nSlaterDets+a] -= 0.5 * x * g_psrq;

	}
      else if(SlaterDet_pair_list[pairIdx].nDiff == 1)
	{
	  int p = SlaterDet_pair_list[pairIdx].SOs_a[0];
	  int q = SlaterDet_pair_list[pairIdx].SOs_b[0];
	  int pos_p = SlaterDet_pair_list[pairIdx].SOs_a_pos[0];
	  int pos_q = SlaterDet_pair_list[pairIdx].SOs_b_pos[0];

	  for(int rr = 0; rr < nEl; rr++)
	    {
	      int r = SlaterDetList[a].SO_list[rr];
	      if(r == p)
		continue;

	      int sign_a = 1;
	      if(pos_p%2 == 1)
		sign_a *= -1;		
	      if(rr%2 == 1)
		sign_a *= -1;
	      if(r > p)
		sign_a *= -1;
	      
	      for(int ss = 0; ss < nEl; ss++)
		{
		  int s = SlaterDetList[b].SO_list[ss];
		  if(s == q)
		    continue;
		  if(s != r)
		    continue;

		  int sign_b = 1;
		  if(pos_q%2 == 1)
		    sign_b *= -1;
		  if(ss%2 == 1)
		    sign_b *= -1;
		  if(s > q)
		    sign_b *= -1;

		  // OK, we have a contribution
		  ergo_real x = sign_a * sign_b;
		  //x *= 2; // To account for a > b case

		  ergo_real g_pqrs = g_SO->x[p][q][r][s];
		  ergo_real g_rspq = g_SO->x[r][s][p][q];
		  ergo_real g_rqps = g_SO->x[r][q][p][s];
		  ergo_real g_psrq = g_SO->x[p][s][r][q];

		  //printf("2el contrib to a b = %i %i : %22.11f\n", a, b, 0.5 * x * g_pqrs);

		  resultMatrix[a*nSlaterDets+b] += 0.5 * x * g_pqrs;
		  resultMatrix[a*nSlaterDets+b] += 0.5 * x * g_rspq;
		  resultMatrix[a*nSlaterDets+b] -= 0.5 * x * g_rqps;
		  resultMatrix[a*nSlaterDets+b] -= 0.5 * x * g_psrq;

		  resultMatrix[b*nSlaterDets+a] += 0.5 * x * g_pqrs;
		  resultMatrix[b*nSlaterDets+a] += 0.5 * x * g_rspq;
		  resultMatrix[b*nSlaterDets+a] -= 0.5 * x * g_rqps;
		  resultMatrix[b*nSlaterDets+a] -= 0.5 * x * g_psrq;

		} // END FOR ss
	    } // END FOR rr
	}
      else
	{
	  for(int pp = 0; pp < nEl; pp++)
	    {
	      int p = SlaterDetList[a].SO_list[pp];
	      for(int rr = 0; rr < nEl; rr++)
		{
		  int r = SlaterDetList[a].SO_list[rr];
		  int sign_a = 1;
		  if(pp%2 == 1)
		    sign_a *= -1;
		  if(rr%2 == 1)
		    sign_a *= -1;
		  if(pp < rr)
		    sign_a *= -1;		    
		  int qlist[2];
		  qlist[0] = pp;
		  qlist[1] = rr;
		  for(int qqq = 0; qqq < 2; qqq++)
		    {
		      int qq = qlist[qqq];
		      int q = SlaterDetList[b].SO_list[qq];
		      int slist[2];
		      slist[0] = pp;
		      slist[1] = rr;
		      for(int sss = 0; sss < 2; sss++)
			{
			  int ss = slist[sss];
			  int s = SlaterDetList[b].SO_list[ss];
			  if(s == q)
			    continue;
			  int sign_b = 1;
			  if(qq%2 == 1)
			    sign_b *= -1;
			  if(ss%2 == 1)
			    sign_b *= -1;
			  if(qq < ss)
			    sign_b *= -1;

			  // OK, we have a contribution
			  ergo_real x = sign_a * sign_b;
			  ergo_real g_pqrs = g_SO->x[p][q][r][s];

			  //printf("2el contrib to a b = %i %i : %22.11f\n", a, b, 0.5 * x * g_pqrs);

			  resultMatrix[a*nSlaterDets+b] += 0.5 * x * g_pqrs;
			  //resultMatrix[b*nSlaterDets+a] += 0.5 * x * g_pqrs;

			} // END FOR sss
		    } // END FOR qqq
		} // END FOR rr
	    } // END FOR pp

	} // END ELSE

    } // END FOR pair
  
}








int get_1e_density_matrix(int nSOs,
			  int nEl,
			  two_idx_SO_struct* D,
			  int nSlaterDets, 
			  const SlaterDet_struct* SlaterDetList, 
			  const ergo_real* coeffList)
{
  int p, q;
  for(p = 0; p < nSOs; p++)
    for(q = 0; q < nSOs; q++)
      {
	// Compute element pq
	ergo_real sum = 0;
	int a, b;
	for(a = 0; a < nSlaterDets; a++)
	  for(b = 0; b < nSlaterDets; b++)
	    {
	      int i, count, savedCount;
	      int sign = 1;
	      // annihilate p from a
	      char SO_list_mod_a_1[MAX_ELECTRONS];
	      count = 0;
	      savedCount = -1;
	      for(i = 0; i < nEl; i++)
		{
		  if(SlaterDetList[a].SO_list[i] == p)
		    savedCount = count;
		  else
		    {
		      SO_list_mod_a_1[count] = SlaterDetList[a].SO_list[i];
		      count++;
		      if(savedCount < 0)
			sign *= -1;
		    }
		}
	      if(savedCount < 0)
		continue;

	      // annihilate q from b
	      char SO_list_mod_b_1[MAX_ELECTRONS];
	      count = 0;
	      savedCount = -1;
	      for(i = 0; i < nEl; i++)
		{
		  if(SlaterDetList[b].SO_list[i] == q)
		    savedCount = count;
		  else
		    {
		      SO_list_mod_b_1[count] = SlaterDetList[b].SO_list[i];
		      count++;
		      if(savedCount < 0)
			sign *= -1;
		    }
		}
	      if(savedCount < 0)
		continue;

	      int equal = 1;
	      for(i = 0; i < nEl-1; i++)
		{
		  if(SO_list_mod_a_1[i] != SO_list_mod_b_1[i])
		    equal = 0;
		}
	      
	      if(equal)
		sum += sign * coeffList[a] * coeffList[b];
	    } // END FOR a b
	D->x[p][q] = sum;
      } // END FOR p q

  // Verify symmetry
  for(p = 0; p < nSOs; p++)
    for(q = 0; q < nSOs; q++)
      {
	if(template_blas_fabs(D->x[p][q] - D->x[q][p]) > 1e-9)
	  printf("ERROR 11!\n");
      }
  ergo_real sum = 0;
  for(p = 0; p < nSOs; p++)
    sum += D->x[p][p];

  return 0;
}


int get_2e_density_matrix(int nSOs,
			  int nEl,
			  four_idx_SO_struct* d,
			  int nSlaterDets, 
			  const SlaterDet_struct* SlaterDetList, 
			  const ergo_real* coeffList)
{
  int p, q, r, s;
  int a, b;
  for(p = 0; p < nSOs; p++)
    for(r = 0; r < nSOs; r++)
      for(a = 0; a < nSlaterDets; a++)
	{
	  int i, count, savedCount;
	  int sign1 = 1;
	  // annihilate p from a
	  char SO_list_mod_a_1[MAX_ELECTRONS];
	  count = 0;
	  savedCount = -1;
	  for(i = 0; i < nEl; i++)
	    {
	      if(SlaterDetList[a].SO_list[i] == p)
		savedCount = count;
	      else
		{
		  SO_list_mod_a_1[count] = SlaterDetList[a].SO_list[i];
		  count++;
		  if(savedCount < 0)
		    sign1 *= -1;
		}
	    }
	  if(savedCount < 0)
	    continue;
	  // annihilate r from a
	  char SO_list_mod_a_2[MAX_ELECTRONS];
	  count = 0;
	  savedCount = -1;
	  for(i = 0; i < nEl-1; i++)
	    {
	      if(SO_list_mod_a_1[i] == r)
		savedCount = count;
	      else
		{
		  SO_list_mod_a_2[count] = SO_list_mod_a_1[i];
		  count++;
		  if(savedCount < 0)
		    sign1 *= -1;
		}
	    }
	  if(savedCount < 0)
	    continue;
	  
	  for(q = 0; q < nSOs; q++)
	    for(s = 0; s < nSOs; s++)
	      for(b = 0; b < nSlaterDets; b++)
		{
		  int sign = sign1;
		  
		  // annihilate q from b
		  char SO_list_mod_b_1[MAX_ELECTRONS];
		  count = 0;
		  savedCount = -1;
		  for(i = 0; i < nEl; i++)
		    {
		      if(SlaterDetList[b].SO_list[i] == q)
			savedCount = count;
		      else
			{
			  SO_list_mod_b_1[count] = SlaterDetList[b].SO_list[i];
			  count++;
			  if(savedCount < 0)
			    sign *= -1;
			}
		    }
		  if(savedCount < 0)
		    continue;
		  // annihilate s from b
		  char SO_list_mod_b_2[MAX_ELECTRONS];
		  count = 0;
		  savedCount = -1;
		  for(i = 0; i < nEl-1; i++)
		    {
		      if(SO_list_mod_b_1[i] == s)
			savedCount = count;
		      else
			{
			  SO_list_mod_b_2[count] = SO_list_mod_b_1[i];
			  count++;
			  if(savedCount < 0)
			    sign *= -1;
			}
		    }
		  if(savedCount < 0)
		    continue;
		  
		  int equal = 1;
		  for(i = 0; i < nEl-2; i++)
		    {
		      if(SO_list_mod_a_2[i] != SO_list_mod_b_2[i])
			equal = 0;
		    }
		  
		  if(equal)
		    d->x[p][q][r][s] += sign * coeffList[a] * coeffList[b];
		} // END FOR q s b
	} // END FOR p r a

  //printf("d done\n");

  // Verify symmetry
  for(p = 0; p < nSOs; p++)
    for(q = 0; q < nSOs; q++)
      for(r = 0; r < nSOs; r++)
	for(s = 0; s < nSOs; s++)
	  {
	    if(template_blas_fabs(d->x[p][q][r][s] - d->x[r][s][p][q]) > 1e-9)
	      printf("ERROR 1!");
	    if(template_blas_fabs(d->x[p][q][r][s] + d->x[r][q][p][s]) > 1e-9)
	      printf("ERROR 2!");
	    if(template_blas_fabs(d->x[p][q][r][s] + d->x[p][s][r][q]) > 1e-9)
	      printf("ERROR 3!");
	    if(template_blas_fabs(d->x[p][q][p][s]) > 1e-9)
	      printf("ERROR 4!");
	    if(template_blas_fabs(d->x[p][q][r][q]) > 1e-9)
	      printf("ERROR 5!");
	    if(template_blas_fabs(d->x[p][q][p][q]) > 1e-9)
	      printf("ERROR 6!");
	  }
  
  return 0;
}






ergo_real get_CI_energy(int nSOs,
			int nEl,
			const four_idx_SO_struct* g_SO, 
			const two_idx_SO_struct* h_SO, 
			int nSlaterDets, 
			const SlaterDet_struct* SlaterDetList, 
			const ergo_real* coeffList,
			ergo_real nuclearEnergy)
{
  two_idx_SO_struct* D = new two_idx_SO_struct;
  four_idx_SO_struct* d = new four_idx_SO_struct;

  get_1e_density_matrix(nSOs,
			nEl,
			D, 
			nSlaterDets, 
			SlaterDetList, 
			coeffList);
  
  get_2e_density_matrix(nSOs,
			nEl,
			d, 
			nSlaterDets, 
			SlaterDetList, 
			coeffList);

  int p, q, r, s;
  
  ergo_real energy_1el = 0;
  for(p = 0; p < nSOs; p++)
    for(q = 0; q < nSOs; q++)
      energy_1el += D->x[p][q] * h_SO->x[p][q];
  
  ergo_real energy_2el = 0;
  for(p = 0; p < nSOs; p++)
    for(q = 0; q < nSOs; q++)
      for(r = 0; r < nSOs; r++)
	for(s = 0; s < nSOs; s++)
	  energy_2el += 0.5 * d->x[p][q][r][s] * g_SO->x[p][q][r][s];

  delete D;
  delete d;
  
  return energy_1el + energy_2el + nuclearEnergy;
}



void get_CI_energy_and_gradient(int nSOs,
				int nEl,
				const four_idx_SO_struct* g_SO, 
				const two_idx_SO_struct* h_SO, 
				int nSlaterDets, 
				const SlaterDet_struct* SlaterDetList, 
				int nSlaterDetPairs,
				const SlaterDet_pair_struct* SlaterDet_pair_list,
				const int* pairCountList,
				int noOfTrialVectors,
				ergo_real* energyList,
				ergo_real** coeffListList,
				ergo_real** gradientList,
				ergo_real nuclearEnergy)
{
  ergo_real* gradient_1el_list[88];
  for(int i = 0; i < noOfTrialVectors; i++)
    gradient_1el_list[i] = new ergo_real[nSlaterDets];
  ergo_real energy_1el_list[88];
  get_1el_energy_and_gradient(nSOs,
			      nEl,
			      nSlaterDets, 
			      SlaterDetList, 
			      nSlaterDetPairs,
			      SlaterDet_pair_list,
			      pairCountList,
			      noOfTrialVectors,
			      energy_1el_list,
			      coeffListList,
			      h_SO,
			      gradient_1el_list);
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "1-electron part done");

  //printf("energy_1el_list[0] = %22.11f\n", energy_1el_list[0]);

  ergo_real* gradient_2el_list[88];
  for(int i = 0; i < noOfTrialVectors; i++)
    gradient_2el_list[i] = new ergo_real[nSlaterDets];
  ergo_real energy_2el_list[88];
  get_2el_energy_and_gradient(nSOs,
			      nEl,
			      nSlaterDets,
			      SlaterDetList,
			      nSlaterDetPairs,
			      SlaterDet_pair_list,
			      pairCountList,
			      noOfTrialVectors,
			      energy_2el_list,
			      coeffListList,
			      g_SO,
			      gradient_2el_list);
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "2-electron part done");
  
  // Take into account that we want the derivative when the sum of squares is fixed.
  for(int k = 0; k < noOfTrialVectors; k++)
    for(int i = 0; i < nSlaterDets; i++)
      gradientList[k][i] = gradient_1el_list[k][i] + gradient_2el_list[k][i] - 2 * (energy_1el_list[k] + energy_2el_list[k]) * coeffListList[k][i];
    
  for(int k = 0; k < noOfTrialVectors; k++)
    {
      delete gradient_1el_list[k];
      delete gradient_2el_list[k];
    }
  
  for(int k = 0; k < noOfTrialVectors; k++)
    energyList[k] = energy_1el_list[k] + energy_2el_list[k] + nuclearEnergy;
}



void mult_by_CI_matrix(int nSOs,
		       int nEl,
		       const four_idx_SO_struct* g_SO, 
		       const two_idx_SO_struct* h_SO, 
		       int nSlaterDets, 
		       const SlaterDet_struct* SlaterDetList, 
		       int nSlaterDetPairs,
		       const SlaterDet_pair_struct* SlaterDet_pair_list,
		       const int* pairCountList,
		       const ergo_real* sourceVector,
		       ergo_real* resultVector,
		       ergo_real shift)
{
  memset(resultVector, 0, nSlaterDets*sizeof(ergo_real));
  for(int i = 0; i < nSlaterDetPairs; i++)
    {
      get_1el_contribs_to_mult_or_dmat(nSOs,
				       nEl,
				       nSlaterDets, 
				       SlaterDetList, 
				       &SlaterDet_pair_list[i],
				       h_SO,
				       sourceVector,
				       resultVector,
				       NULL);
      get_2el_contribs_to_mult_or_dmat(nSOs,
				       nEl,
				       nSlaterDets,
				       SlaterDetList,
				       &SlaterDet_pair_list[i],
				       g_SO,
				       sourceVector,
				       resultVector,
				       NULL);
    } // END FOR i
  // Apply shift
  for(int i = 0; i < nSlaterDets; i++)
    resultVector[i] -= shift * sourceVector[i];
}



void get_CI_matrix(int nSOs,
		   int nEl,
		   const four_idx_SO_struct* g_SO, 
		   const two_idx_SO_struct* h_SO, 
		   int nSlaterDets, 
		   const SlaterDet_struct* SlaterDetList, 
		   int nSlaterDetPairs,
		   const SlaterDet_pair_struct* SlaterDet_pair_list,
		   const int* pairCountList,
		   ergo_real* resultMatrix)
{
  memset(resultMatrix, 0, nSlaterDets*nSlaterDets*sizeof(ergo_real));
  get_1el_contribs(nSOs,
		   nEl,
		   nSlaterDets, 
		   SlaterDetList, 
		   nSlaterDetPairs,
		   SlaterDet_pair_list,
		   pairCountList,
		   h_SO,
		   resultMatrix);
  get_2el_contribs(nSOs,
		   nEl,
		   nSlaterDets,
		   SlaterDetList,
		   nSlaterDetPairs,
		   SlaterDet_pair_list,
		   pairCountList,
		   g_SO,
		   resultMatrix);
}



int get_combinations(SlaterDet_struct* SlaterDetList, 
		     int nEl, 
		     int nSOs)
{
  if(nEl == 0)
    return 1;

  int list[MAX_SOS];
  memset(list, 0, MAX_SOS * sizeof(int));

  int i;
  for(i = 0; i < nEl; i++)
    list[i] = i;

  int finishedFlag = 0;
  int count = 0;
  while(finishedFlag == 0)
    {
      if(SlaterDetList)
	{
	  // Create det for current configuration
	  for(i = 0; i < nEl; i++)
	    SlaterDetList[count].SO_list[i] = list[i];
	}
      count++;

      // Move to next configuration
      int breakFlag = 0;
      int level = nEl - 1;      
      while(breakFlag == 0)
	{
	  // Check if value at current level can be increased
	  if(list[level] < nSOs-(nEl-level))
	    {
	      // OK, increase value
	      list[level]++;
	      // Set all values at higher levels as low as possible
	      for(i = level+1; i < nEl; i++)
		list[i] = list[i-1] + 1;
	      breakFlag = 1;
	    }
	  else
	    {
	      // Go to next level
	      level--;
	    }
	  if(level < 0)
	    {
	      breakFlag = 1;
	      finishedFlag = 1;
	    }
	}
    }

  return count;  
}


int get_FCI_Slater_dets_alpha_beta(SlaterDet_struct* SlaterDetList, 
				   int nEl_a, 
				   int nEl_b, 
				   int nSOs)
{
  int nCombs_a = get_combinations(NULL, nEl_a, nSOs/2);  
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "get_FCI_Slater_dets_alpha_beta, nCombs_a = %9i", nCombs_a);
  SlaterDet_struct* SlaterDetList_a = new SlaterDet_struct[nCombs_a];
  get_combinations(SlaterDetList_a, nEl_a, nSOs/2);

  int nCombs_b = get_combinations(NULL, nEl_b, nSOs/2);
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "get_FCI_Slater_dets_alpha_beta, nCombs_b = %9i", nCombs_b);
  SlaterDet_struct* SlaterDetList_b = new SlaterDet_struct[nCombs_b];
  get_combinations(SlaterDetList_b, nEl_b, nSOs/2);

  int count = 0;
  int ia, ib;
  for(ia = 0; ia < nCombs_a; ia++)
    for(ib = 0; ib < nCombs_b; ib++)
      {
	if(SlaterDetList)
	  {
	    for(int i = 0; i < nEl_a; i++)
	      SlaterDetList[count].SO_list[i] = SlaterDetList_a[ia].SO_list[i];
	    for(int i = 0; i < nEl_b; i++)
	      SlaterDetList[count].SO_list[nEl_a+i] = nSOs / 2 + SlaterDetList_b[ib].SO_list[i];
	    SlaterDetList[count].startIndex = -1;
	    SlaterDetList[count].count = 0;
	  }

	count++;
      }

  delete [] SlaterDetList_a;
  delete [] SlaterDetList_b;

  return count;
}



int get_FCI_Slater_dets_all(SlaterDet_struct* SlaterDetList, 
			    int nElTot, 
			    int nSOs)
{
  return get_combinations(SlaterDetList, nElTot, nSOs);
}



static ergo_real rand_m1_to_1()
{
  ergo_real randomNumber = (ergo_real)rand() / RAND_MAX;
  // Now randomNumber is between 0 and 1
  randomNumber *= 2;
  // Now randomNumber is between 0 and 2
  randomNumber -= 1;
  // Now randomNumber is between -1 and 1
  return randomNumber;
}



typedef struct
{
  int aDiffList[2];
  int bDiffList[2];
  int aDiffPosList[2];
  int bDiffPosList[2];
} pair_status_struct;





void get_relevant_SlaterDet_pairs_recursive_do_mult_or_dmat(SlaterDet_struct** groupList, 
							    int nEl,
							    int level,
							    int groupIdx1,
							    int groupIdx2,
							    int ia,
							    int ib,
							    int aDiffCount,
							    int bDiffCount,
							    pair_status_struct* status,
							    int nSOs,
							    int nSlaterDets,
							    const SlaterDet_struct* SlaterDetList,
							    const four_idx_SO_struct* g_SO, 
							    const two_idx_SO_struct* h_SO, 
							    const ergo_real* sourceVector,
							    ergo_real* resultVector,
							    two_idx_SO_struct* result_dmat_1el,
							    four_idx_SO_struct* result_dmat_2el
							    )
{
  if(groupIdx1 > groupIdx2)
    return;

  if(ia != level-1 && ib != level-1 && level != 0)
    throw "Error in CI: (ia != level && ib != level)";

  // ia or ib or both point to an electron at this level.
  const SlaterDet_struct* group_1 = &groupList[level][groupIdx1];
  const SlaterDet_struct* group_2 = &groupList[level][groupIdx2];

  if(level == nEl)
    {
      // final level
      while(ia < level || ib < level)
	{
	  if(ia == nEl)
	    {
	      // Only b left, must be diff
	      if(bDiffCount == 2)
		return;
	      status->bDiffList[bDiffCount] = group_2->SO_list[ib];
	      status->bDiffPosList[bDiffCount] = ib;
	      bDiffCount++;
	      ib++;
	      continue;
	    }
	  if(ib == nEl)
	    {
	      // Only a left, must be diff
	      if(aDiffCount == 2)
		return;
	      status->aDiffList[aDiffCount] = group_1->SO_list[ia];
	      status->aDiffPosList[aDiffCount] = ia;
	      aDiffCount++;
	      ia++;
	      continue;
	    }
	  while(ia < level && ib < level)
	    {
	      if(group_1->SO_list[ia] == group_2->SO_list[ib])
		{
		  ia++;
		  ib++;
		}
	      else if(group_1->SO_list[ia] > group_2->SO_list[ib])
		{
		  if(bDiffCount == 2)
		    return;
		  status->bDiffList[bDiffCount] = group_2->SO_list[ib];
		  status->bDiffPosList[bDiffCount] = ib;
		  bDiffCount++;
		  ib++;
		}
	      else
		{	
		  if(aDiffCount == 2)
		    return;
		  status->aDiffList[aDiffCount] = group_1->SO_list[ia];
		  status->aDiffPosList[aDiffCount] = ia;
		  aDiffCount++;
		  ia++;
		}
	    } // END WHILE
	} // END WHILE

      // Do contribution to matrix-vector multiplication.
      SlaterDet_pair_struct SlaterDetPair;
      SlaterDetPair.a = groupIdx1;
      SlaterDetPair.b = groupIdx2;
      SlaterDetPair.nDiff = aDiffCount;
      for(int i = 0; i < aDiffCount; i++)
	{
	  SlaterDetPair.SOs_a[i] = status->aDiffList[i];
	  SlaterDetPair.SOs_b[i] = status->bDiffList[i];
	  SlaterDetPair.SOs_a_pos[i] = status->aDiffPosList[i];
	  SlaterDetPair.SOs_b_pos[i] = status->bDiffPosList[i];
	}
      get_1el_contribs_to_mult_or_dmat(nSOs,
				       nEl,
				       nSlaterDets, 
				       SlaterDetList, 
				       &SlaterDetPair,
				       h_SO,
				       sourceVector,
				       resultVector,
				       result_dmat_1el);
      get_2el_contribs_to_mult_or_dmat(nSOs,
				       nEl,
				       nSlaterDets,
				       SlaterDetList,
				       &SlaterDetPair,
				       g_SO,
				       sourceVector,
				       resultVector,
				       result_dmat_2el);
      
      return;
    } // END IF final level

  if(level > 0)
    {
      while(ia < level && ib < level)
	{
	  if(group_1->SO_list[ia] == group_2->SO_list[ib])
	    {
	      ia++;
	      ib++;
	    }
	  else if(group_1->SO_list[ia] > group_2->SO_list[ib])
	    {
	      if(bDiffCount == 2)
		return;
	      status->bDiffList[bDiffCount] = group_2->SO_list[ib];
	      status->bDiffPosList[bDiffCount] = ib;
	      bDiffCount++;
	      ib++;
	    }
	  else
	    {
	      if(aDiffCount == 2)
		return;
	      status->aDiffList[aDiffCount] = group_1->SO_list[ia];
	      status->aDiffPosList[aDiffCount] = ia;
	      aDiffCount++;
	      ia++;
	    }
	} // END WHILE
    } // END IF level > 0
 
  // No, we could not skip. Go to next level.
  for(int i = 0; i < group_1->count; i++)
    for(int j = 0; j < group_2->count; j++)
      {
	get_relevant_SlaterDet_pairs_recursive_do_mult_or_dmat(groupList, 
							       nEl,
							       level+1,
							       group_1->startIndex + i,
							       group_2->startIndex + j,
							       ia,
							       ib,
							       aDiffCount,
							       bDiffCount,
							       status,
							       nSOs,
							       nSlaterDets,
							       SlaterDetList,
							       g_SO, 
							       h_SO, 
							       sourceVector,
							       resultVector,
							       result_dmat_1el,
							       result_dmat_2el
							       );
      } // END FOR i j
  return;
}





void mult_by_CI_matrix_direct(int nSOs,
			      int nEl,
			      const four_idx_SO_struct* g_SO, 
			      const two_idx_SO_struct* h_SO, 
			      int nSlaterDets, 
			      const SlaterDet_struct* SlaterDetList, 
			      SlaterDet_struct** groupList, 
			      const ergo_real* sourceVector,
			      ergo_real* resultVector,
			      ergo_real shift)
{
  memset(resultVector, 0, nSlaterDets*sizeof(ergo_real));

  pair_status_struct status;
  memset(&status, 0, sizeof(pair_status_struct));
  get_relevant_SlaterDet_pairs_recursive_do_mult_or_dmat(groupList, 
							 nEl,
							 0,
							 0,
							 0,
							 0, 0, 0, 0,
							 &status,
							 nSOs,
							 nSlaterDets,
							 SlaterDetList,
							 g_SO, 
							 h_SO, 
							 sourceVector,
							 resultVector,
							 NULL,
							 NULL
							 );
  
  // Apply shift
  for(int i = 0; i < nSlaterDets; i++)
    resultVector[i] -= shift * sourceVector[i];
}



static ergo_real get_SlaterDet_energy(int nSOs,
				      int nEl,
				      const four_idx_SO_struct* g_SO, 
				      const two_idx_SO_struct* h_SO, 
				      const SlaterDet_struct* SlaterDet)
{
  // Create "Slater det pair" consisting of this determinant paired with itself.
  SlaterDet_pair_struct SlaterDetPair;
  SlaterDetPair.a = 0;
  SlaterDetPair.b = 0;
  SlaterDetPair.nDiff = 0;
  const int nSlaterDets = 1;
  ergo_real sourceVector[1];
  ergo_real resultVector[1];
  sourceVector[0] = 1.0;
  resultVector[0] = 0.0;
  get_1el_contribs_to_mult_or_dmat(nSOs,
				   nEl,
				   nSlaterDets, 
				   SlaterDet, 
				   &SlaterDetPair,
				   h_SO,
				   sourceVector,
				   resultVector,
				   NULL);
  get_2el_contribs_to_mult_or_dmat(nSOs,
				   nEl,
				   nSlaterDets,
				   SlaterDet,
				   &SlaterDetPair,
				   g_SO,
				   sourceVector,
				   resultVector,
				   NULL);
  return resultVector[0];
}



int get_relevant_SlaterDet_pairs_recursive_2(SlaterDet_struct** groupList, 
					     int nEl,
					     SlaterDet_pair_struct* resultList,
					     int level,
					     int groupIdx1,
					     int groupIdx2,
					     int ia,
					     int ib,
					     int aDiffCount,
					     int bDiffCount,
					     pair_status_struct* status)
{
  if(groupIdx1 > groupIdx2)
    return 0;

  if(ia != level-1 && ib != level-1 && level != 0)
    throw "Error in CI: (ia != level && ib != level)";

  // ia or ib or both point to an electron at this level.
  SlaterDet_struct* group_1 = &groupList[level][groupIdx1];
  SlaterDet_struct* group_2 = &groupList[level][groupIdx2];

  if(level == nEl)
    {
      // final level
      while(ia < level || ib < level)
	{
	  if(ia == nEl)
	    {
	      // Only b left, must be diff
	      if(bDiffCount == 2)
		return 0;
	      status->bDiffList[bDiffCount] = group_2->SO_list[ib];
	      status->bDiffPosList[bDiffCount] = ib;
	      bDiffCount++;
	      ib++;
	      continue;
	    }
	  if(ib == nEl)
	    {
	      // Only a left, must be diff
	      if(aDiffCount == 2)
		return 0;
	      status->aDiffList[aDiffCount] = group_1->SO_list[ia];
	      status->aDiffPosList[aDiffCount] = ia;
	      aDiffCount++;
	      ia++;
	      continue;
	    }
	  while(ia < level && ib < level)
	    {
	      if(group_1->SO_list[ia] == group_2->SO_list[ib])
		{
		  ia++;
		  ib++;
		}
	      else if(group_1->SO_list[ia] > group_2->SO_list[ib])
		{
		  if(bDiffCount == 2)
		    return 0;
		  status->bDiffList[bDiffCount] = group_2->SO_list[ib];
		  status->bDiffPosList[bDiffCount] = ib;
		  bDiffCount++;
		  ib++;
		}
	      else
		{	
		  if(aDiffCount == 2)
		    return 0;
		  status->aDiffList[aDiffCount] = group_1->SO_list[ia];
		  status->aDiffPosList[aDiffCount] = ia;
		  aDiffCount++;
		  ia++;
		}
	    } // END WHILE
	} // END WHILE
      
      if(resultList)
	{
	  resultList[0].a = groupIdx1;
	  resultList[0].b = groupIdx2;
	  resultList[0].nDiff = aDiffCount;
	  for(int i = 0; i < aDiffCount; i++)
	    {
	      resultList[0].SOs_a[i] = status->aDiffList[i];
	      resultList[0].SOs_b[i] = status->bDiffList[i];
	      resultList[0].SOs_a_pos[i] = status->aDiffPosList[i];
	      resultList[0].SOs_b_pos[i] = status->bDiffPosList[i];
	    }
	}
      return 1;
    } // END IF final level

  if(level > 0)
    {
      while(ia < level && ib < level)
	{
	  if(group_1->SO_list[ia] == group_2->SO_list[ib])
	    {
	      ia++;
	      ib++;
	    }
	  else if(group_1->SO_list[ia] > group_2->SO_list[ib])
	    {
	      if(bDiffCount == 2)
		return 0;
	      status->bDiffList[bDiffCount] = group_2->SO_list[ib];
	      status->bDiffPosList[bDiffCount] = ib;
	      bDiffCount++;
	      ib++;
	    }
	  else
	    {
	      if(aDiffCount == 2)
		return 0;
	      status->aDiffList[aDiffCount] = group_1->SO_list[ia];
	      status->aDiffPosList[aDiffCount] = ia;
	      aDiffCount++;
	      ia++;
	    }
	} // END WHILE
    } // END IF level > 0
 
  // No, we could not skip. Go to next level.
  int count = 0;
  for(int i = 0; i < group_1->count; i++)
    for(int j = 0; j < group_2->count; j++)
      {
	SlaterDet_pair_struct* resultListPtr = NULL;
	if(resultList)
	  resultListPtr = &resultList[count];
	int currCount = get_relevant_SlaterDet_pairs_recursive_2(groupList, 
								 nEl,
								 resultListPtr,
								 level+1,
								 group_1->startIndex + i,
								 group_2->startIndex + j,
								 ia,
								 ib,
								 aDiffCount,
								 bDiffCount,
								 status);
	count += currCount;
      } // END FOR i j
  return count;
}






int get_relevant_SlaterDet_pairs_recursive(int nSlaterDets, 
					   SlaterDet_struct* SlaterDetList, 
					   SlaterDet_struct** groupList, 
					   int nEl,
					   SlaterDet_pair_struct* resultList,
					   int level,
					   int groupIdx1,
					   int groupIdx2)
{
  //printf("get_relevant_SlaterDet_pairs_recursive, level = %i\n", level);

  // Check if this pair of groups can be skipped.
  SlaterDet_struct* group_1 = &groupList[level][groupIdx1];
  SlaterDet_struct* group_2 = &groupList[level][groupIdx2];
  if(level < nEl)
    {
      int nSame = 0;
      int ia = 0;
      int ib = 0;
      int aDiffCount = 0;
      int bDiffCount = 0;
      //int aDiffList[MAX_SOS];
      //int bDiffList[MAX_SOS];
      //int aDiffPosList[MAX_SOS];
      //int bDiffPosList[MAX_SOS];
      while(ia < level && ib < level)
	{
	  if(group_1->SO_list[ia] == group_2->SO_list[ib])
	    {
	      nSame++;
	      ia++;
	      ib++;
	    }
	  else
	    {
	      if(group_1->SO_list[ia] > group_2->SO_list[ib])
		{
		  //bDiffList[bDiffCount] = group_2->SO_list[ib];
		  //bDiffPosList[bDiffCount] = ib;
		  bDiffCount++;
		  if(bDiffCount > 2)
		    break;
		  ib++;
		}
	      else
		{
		  //aDiffList[aDiffCount] = group_1->SO_list[ia];
		  //aDiffPosList[aDiffCount] = ia;
		  aDiffCount++;
		  if(aDiffCount > 2)
		    break;
		  ia++;
		}
	    }
	}
      if(aDiffCount > 2 || bDiffCount > 2)
	{
	  //printf("aDiffCount = %i  bDiffCount = %i  level = %i\n", aDiffCount, bDiffCount, level);
	  return 0;
	}
 
      // No, we could not skip. Go to next level.
      int count = 0;
      for(int i = 0; i < group_1->count; i++)
	for(int j = 0; j < group_2->count; j++)
	  {
	    SlaterDet_pair_struct* resultListPtr = NULL;
	    if(resultList)
	      resultListPtr = &resultList[count];
	    int currCount = get_relevant_SlaterDet_pairs_recursive(nSlaterDets, 
								   SlaterDetList, 
								   groupList, 
								   nEl,
								   resultListPtr,
								   level+1,
								   group_1->startIndex + i,
								   group_2->startIndex + j);
	    count += currCount;
	  } // END FOR i j
      return count;
    } // END IF level < nEl
  
  // We are at lowest level.

  //printf("lowest level!\n");

  int a = groupIdx1;
  int b = groupIdx2;

  if(b < a)
    return 0;
  
  int nSame = 0;
  int ia = 0;
  int ib = 0;
  int aDiffCount = 0;
  int bDiffCount = 0;
  int aDiffList[MAX_SOS];
  int bDiffList[MAX_SOS];
  int aDiffPosList[MAX_SOS];
  int bDiffPosList[MAX_SOS];
  //printf("ia ib loop starting.\n");
  while(ia < nEl || ib < nEl)
    {
      if(ia == nEl)
	{
	  // Only b left, must be diff
	  bDiffList[bDiffCount] = SlaterDetList[b].SO_list[ib];
	  bDiffPosList[bDiffCount] = ib;
	  bDiffCount++;
	  ib++;
	  continue;
	}
      if(ib == nEl)
	{
	  // Only a left, must be diff
	  aDiffList[aDiffCount] = SlaterDetList[a].SO_list[ia];
	  aDiffPosList[aDiffCount] = ia;
	  aDiffCount++;
	  ia++;
	  continue;
	}
      if(SlaterDetList[a].SO_list[ia] == SlaterDetList[b].SO_list[ib])
	{
	  nSame++;
	  ia++;
	  ib++;
	}
      else
	{
	  if(SlaterDetList[a].SO_list[ia] > SlaterDetList[b].SO_list[ib])
	    {
	      bDiffList[bDiffCount] = SlaterDetList[b].SO_list[ib];
	      bDiffPosList[bDiffCount] = ib;
	      bDiffCount++;
	      if(bDiffCount > 2)
		break;
	      ib++;
	    }
	  else
	    {
	      aDiffList[aDiffCount] = SlaterDetList[a].SO_list[ia];
	      aDiffPosList[aDiffCount] = ia;
	      aDiffCount++;
	      if(aDiffCount > 2)
		break;
	      ia++;
	    }
	}
    }
  if(nSame >= nEl - 2)
    {
      if(resultList)
	{
	  int pairCount = 0;
	  resultList[pairCount].a = a;
	  resultList[pairCount].b = b;
	  resultList[pairCount].nDiff = aDiffCount;
	  for(int i = 0; i < 2; i++)
	    {
	      resultList[pairCount].SOs_a[i] = -1;
	      resultList[pairCount].SOs_b[i] = -1;
	      resultList[pairCount].SOs_a_pos[i] = -1;
	      resultList[pairCount].SOs_b_pos[i] = -1;
	    }
	  for(int i = 0; i < aDiffCount; i++)
	    {
	      resultList[pairCount].SOs_a[i] = aDiffList[i];
	      resultList[pairCount].SOs_b[i] = bDiffList[i];
	      resultList[pairCount].SOs_a_pos[i] = aDiffPosList[i];
	      resultList[pairCount].SOs_b_pos[i] = bDiffPosList[i];
	    }
	}
      return 1;
    }
  else
    {
      return 0;
    }
}



int get_relevant_SlaterDet_pairs(int nSlaterDets, 
				 SlaterDet_struct* SlaterDetList, 
				 SlaterDet_struct** groupList, 
				 int nEl,
				 SlaterDet_pair_struct* resultList)
{
  pair_status_struct status;
  memset(&status, 0, sizeof(pair_status_struct));
  return get_relevant_SlaterDet_pairs_recursive_2(groupList, 
						  nEl,
						  resultList,
						  0,
						  0,
						  0,
						  0, 0, 0, 0,
						  &status);
}




ergo_real get_eigs(int n, ergo_real* M, ergo_real* bestVector, ergo_real* eigValListResult)
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
      exit(0);
    }

  if(bestVector)
    {
      for(int i = 0; i < n; i++)
	bestVector[i] = A[0*n+i];
    }

  if(eigValListResult)
    {
      for(int i = 0; i < n; i++)
	eigValListResult[i] = eigvalList[i];
    }

  ergo_real result = eigvalList[0];

  delete [] work;
  delete [] eigvalList;
  delete [] A;

  return result;
}






int get_Lowdin_orbitals(int n, const ergo_real* S, ergo_real* MOs)
{
  int lwork = 3*n*n;
  ergo_real* work = new ergo_real[lwork];
  ergo_real* eigvalList = new ergo_real[n];
  ergo_real* A = new ergo_real[n*n];
  memcpy(A, S, n*n*sizeof(ergo_real));
  int info = 0;

  mat::syev("V", "U", &n, A,
	    &n, eigvalList, work, &lwork, 
	    &info);
  if(info != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_CI, "error in syev, info = %i", info);
      exit(0);
    }
  
  for(int i = 0; i < n; i++)
    {
      assert(eigvalList[i] > 0);
    }

  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      {
	ergo_real sum = 0;
	for(int k = 0; k < n; k++)
	  sum += A[k*n+i] * A[k*n+j] * (1.0 / template_blas_sqrt(eigvalList[k]));
	MOs[i*n+j] = sum;
      }

  return 0;
}







struct MatVecMul {
  int nEl;
  int nSOs;
  int nSlaterDets;
  const four_idx_SO_struct* g_SO;
  const two_idx_SO_struct* h_SO;
  const SlaterDet_struct* SlaterDetList;
  SlaterDet_struct** groupList;
  ergo_real shift;
  void do_mat_vec_mul(int n, const ergo_real* sourceVector, ergo_real* resultVector) const {
    mult_by_CI_matrix_direct(nSOs,
			     nEl,
			     g_SO, 
			     h_SO, 
			     nSlaterDets, 
			     SlaterDetList, 
			     groupList, 
			     sourceVector,
			     resultVector,
			     shift);
  }
};



int do_power_method(int n,
		    ergo_real* A,
		    ergo_real* v)
{
  ergo_real* v2 = new ergo_real[n];
  int iter = 0;
  while(1)
    {
      iter++;
      normalize_vector(n, v);
      for(int i = 0; i < n; i++)
	{
	  ergo_real sum = 0;
	  for(int j = 0; j < n; j++)
	    sum += A[i*n+j] * v[j];
	  v2[i] = sum;
	} // END FOR i
      ergo_real norm = get_vector_norm(n, v2);
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "power method, iter = %5i, norm = %22.11f", iter, (double)norm);
      for(int i = 0; i < n; i++)
	v[i] = v2[i];
    } // END WHILE
  return 0;
}


struct pqrs_struct
{
  int p;
  int q;
  int r;
  int s;
  void assign(int pp, int qq, int rr, int ss)
  {
    p = pp;
    q = qq;
    r = rr;
    s = ss;
  }
  int compare(int pp, int qq, int rr, int ss)
  {
    if(pp == p && qq == q && rr == r && ss == s)
      return 0;
    else
      return -1;
  }
};


int do_CI(
	  const BasisInfoStruct & basisInfo, 
	  const IntegralInfo & integralInfo,
	  const CI::Options& options,
	  const Molecule & molecule,
	  const ergo_real* S,
	  const ergo_real* h_AO,
	  const ergo_real* F_a,
	  const ergo_real* F_b,
	  int n_el_a,
	  int n_el_b,
	  ergo_real nuclearEnergy,
	  ergo_real HF_energy
	  )
{
  int n = basisInfo.noOfBasisFuncs;
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "do_CI, n = %3i, HF_energy = %22.11f", n, (double)HF_energy);

  if(n > MAX_AOS)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_CI, "Error in do_CI: (n > MAX_AOS).");
      return -1;
    }

  output_current_memory_usage(LOG_AREA_CI, "beginning of do_CI");

  // Use basisInfo to get g_AO
  four_idx_AO_struct* g_AO = new four_idx_AO_struct;
  int p, q, r, s;
  for(p = 0; p < n; p++)
    for(q = 0; q < n; q++)
      for(r = 0; r < n; r++)
	for(s = 0; s < n; s++)
	  {
	    g_AO->x[p][q][r][s] = do_2e_integral(p, q, r, s, 
						 basisInfo, 
						 integralInfo);
	  }
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "g_AO done");
  
  ergo_real* MOs_a = new ergo_real[n*n];
  ergo_real* MOs_b = new ergo_real[n*n];
  ergo_real* eigv_a = new ergo_real[n];
  ergo_real* eigv_b = new ergo_real[n];


  if(options.use_random_orbitals == 1)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "Using random orbitals, different for alpha and beta.");
      ergo_real* F_random = new ergo_real[n*n];
      for(p = 0; p < n*n; p++)
	F_random[p] = rand_m1_to_1();
      get_F_orbs(n, F_random, S, MOs_a, eigv_a);
      for(p = 0; p < n*n; p++)
	F_random[p] = rand_m1_to_1();
      get_F_orbs(n, F_random, S, MOs_b, eigv_b);
    }
  else if(options.use_lowdin_orbitals == 1)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "Using Lowdin orbitals, same for alpha and beta.");
      get_Lowdin_orbitals(n, S, MOs_a);
      get_Lowdin_orbitals(n, S, MOs_b);
    }
  else
    {
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "Using HF orbitals obtained from alpha and beta Fock matrices.");
      get_F_orbs(n, F_a, S, MOs_a, eigv_a);  
      get_F_orbs(n, F_b, S, MOs_b, eigv_b);
    }
  
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "MOs_a MOs_b done.");

  int nSOs = 2 * n;
  
  SO_struct* SOs = new SO_struct[nSOs];
  int count = 0;
  int i;
  // ALPHA
  for(i = 0; i < n; i++)
    {
      SOs[count].spin = SPIN_A;
      int j;
      for(j = 0; j < n; j++)
	SOs[count].coeffs[j] = MOs_a[i*n+j];
      count++;
    }
  // BETA
  for(i = 0; i < n; i++)
    {
      SOs[count].spin = SPIN_B;
      int j;
      for(j = 0; j < n; j++)
	SOs[count].coeffs[j] = MOs_b[i*n+j];
      count++;
    }
  
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "SOs done.");

  four_idx_SO_struct* g_SO = new four_idx_SO_struct;
  for(p = 0; p < nSOs; p++)
    for(q = 0; q < nSOs; q++)
      for(r = 0; r < nSOs; r++)
	for(s = 0; s < nSOs; s++)
	  {
	    ergo_real value = 0;
	    if(SOs[p].spin == SOs[q].spin && SOs[r].spin == SOs[s].spin)
	      {
		int a, b, c, d;
		for(a = 0; a < n; a++)
		  for(b = 0; b < n; b++)
		    for(c = 0; c < n; c++)
		      for(d = 0; d < n; d++)
			value += SOs[p].coeffs[a] * SOs[q].coeffs[b] * SOs[r].coeffs[c] * SOs[s].coeffs[d] * g_AO->x[a][b][c][d];
	      }
	    g_SO->x[p][q][r][s] = value;
	  }
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "g_SO done");

  two_idx_SO_struct* h_SO = new two_idx_SO_struct;
  for(p = 0; p < nSOs; p++)
    for(q = 0; q < nSOs; q++)
      {
	ergo_real value = 0;
	if(SOs[p].spin == SOs[q].spin)
	  {
	    int a, b;
	    for(a = 0; a < n; a++)
	      for(b = 0; b < n; b++)
		value += SOs[p].coeffs[a] * SOs[q].coeffs[b] * h_AO[a*n+b];
	  }
	h_SO->x[p][q] = value;
      }
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "h_SO done.");

  do_output(LOG_CAT_INFO, LOG_AREA_CI, "n_el_a = %2i, n_el_b = %2i", n_el_a, n_el_b);

  int nElTot = n_el_a + n_el_b;

  if(nElTot > MAX_ELECTRONS)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_CI, "Error in do_CI: (nElTot > MAX_ELECTRONS).");
      return -1;
    }

  int nSlaterDets1 = get_FCI_Slater_dets_alpha_beta(NULL, n_el_a, n_el_b, nSOs);
  SlaterDet_struct* SlaterDetList1 = new SlaterDet_struct[nSlaterDets1];
  get_FCI_Slater_dets_alpha_beta(SlaterDetList1, n_el_a, n_el_b, nSOs);

  do_output(LOG_CAT_INFO, LOG_AREA_CI, "Slater determinants done, nSlaterDets = %5i", nSlaterDets1);

  do_output(LOG_CAT_INFO, LOG_AREA_CI, "Computing energy for each Slater determinant, nSlaterDets = %9i.", nSlaterDets1);
  std::vector<ergo_real> SlaterDetEnergyList(nSlaterDets1);
  for(int i = 0; i < nSlaterDets1; i++)
    {
      SlaterDetEnergyList[i] = get_SlaterDet_energy(nSOs,
						    nElTot,
						    g_SO, 
						    h_SO, 
						    &SlaterDetList1[i]);
    }

  SlaterDet_struct* SlaterDetList = new SlaterDet_struct[nSlaterDets1];
  int nSlaterDets;
  if(options.use_energy_diff_limit == 1)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "Selecting Slater dets using energy_diff_limit = %9.5f.", (double)options.energy_diff_limit);
      int detCount = 0;
      for(int i = 0; i < nSlaterDets1; i++)
	{
	  if(SlaterDetEnergyList[i] + nuclearEnergy - HF_energy < options.energy_diff_limit)
	    {
	      SlaterDetList[detCount] = SlaterDetList1[i];
	      detCount++;
	    }
	}
      nSlaterDets = detCount;
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "Selected Slater dets low enough in energy, nSlaterDets = %9i.", nSlaterDets);
    }
  else
    {
      for(int i = 0; i < nSlaterDets1; i++)
	SlaterDetList[i] = SlaterDetList1[i];
      nSlaterDets = nSlaterDets1;
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "Using all Slater dets, nSlaterDets = %9i.", nSlaterDets);
    }

  if(nSlaterDets < 1)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_CI, "Error: nSlaterDets = %i", nSlaterDets);
      return -1;
    }
  

  // Define groups of Slater determinants, one level at a time
  SlaterDet_struct* groupList[MAX_ELECTRONS+1];
  int groupCountList[MAX_ELECTRONS+1];
  groupList[nElTot] = SlaterDetList;
  groupCountList[nElTot] = nSlaterDets;
  std::vector<SlaterDet_struct> list(nSlaterDets);
  for(int i = 0; i < nElTot; i++)
    {
      // Define group with nElTot-1-i electrons.
      int nElCurrLevel = nElTot-1-i;
      int prevGroupIndex = nElTot-i;
      // Go through prev level to find groups at curr level.
      SlaterDet_struct* prevGroupList = groupList[prevGroupIndex];
      int prevCount = groupCountList[prevGroupIndex];
      int j = 0;
      int nGroupsCurrLevel = 0;
      while(j < prevCount)
	{
	  // Now j points to beginning of a new group
	  SlaterDet_struct newGroup;
	  memset(&newGroup, 0, sizeof(SlaterDet_struct));
	  for(int k = 0; k < nElCurrLevel; k++)
	    newGroup.SO_list[k] = prevGroupList[j].SO_list[k];
	  newGroup.startIndex = j;
	  int count = 0;
	  while(j < prevCount)
	    {
	      int ok = 1;
	      for(int k = 0; k < nElCurrLevel; k++)
		{
		  if(prevGroupList[j].SO_list[k] != newGroup.SO_list[k])
		    ok = 0;
		}
	      if(ok)
		{
		  j++;
		  count++;
		}
	      else
		break;
	    }
	  // Now j points to the beginning of next group or to the end of the list.
	  newGroup.count = count;
	  list[nGroupsCurrLevel] = newGroup;
	  nGroupsCurrLevel++;
	} // END WHILE
      int groupIndex = nElTot-1-i;
      groupList[groupIndex] = new SlaterDet_struct[nGroupsCurrLevel];
      memcpy(groupList[groupIndex], &list[0], nGroupsCurrLevel*sizeof(SlaterDet_struct));
      groupCountList[groupIndex] = nGroupsCurrLevel;
    }
  
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "Slater determinant groups done:");
  for(int i = 0; i < nElTot+1; i++)
    do_output(LOG_CAT_INFO, LOG_AREA_CI, "groupCountList[%2i] = %9i", i, groupCountList[i]);

  ergo_real* coeffList_guess = new ergo_real[nSlaterDets];
  if(options.use_random_starting_guess == 1)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "Using random CI starting guess.");
      for(i = 0; i < nSlaterDets; i++)
	coeffList_guess[i] = rand_m1_to_1();
    }
  else
    {
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "Using vector [ 1 0 0 ... 0 ] as CI starting guess.");
      do_output(LOG_CAT_INFO, LOG_AREA_CI, "If HF orbitals are used, this should correspond to the HF state.");
      for(i = 0; i < nSlaterDets; i++)
	coeffList_guess[i] = 0;
      coeffList_guess[0] = 1;
    }
  normalize_vector(nSlaterDets, coeffList_guess);

  ergo_real* coeffList_hi = new ergo_real[nSlaterDets];
  ergo_real* coeffList_lo = new ergo_real[nSlaterDets];

  ergo_real resultEig_lo = 0;
  ergo_real resultEig_hi = 0;

  MatVecMul matVecMul;
  matVecMul.nEl = nElTot;
  matVecMul.nSOs = nSOs;
  matVecMul.nSlaterDets = nSlaterDets;
  matVecMul.g_SO = g_SO;
  matVecMul.h_SO = h_SO;
  matVecMul.SlaterDetList = SlaterDetList;
  matVecMul.groupList = groupList;
  matVecMul.shift = options.shift;
  simple_lanczos::do_lanczos_method(nSlaterDets,
				    coeffList_guess,
				    resultEig_lo,
				    coeffList_lo,
				    resultEig_hi,
				    coeffList_hi,
				    matVecMul,
				    options.max_no_of_iterations,
				    options.shift,
				    nuclearEnergy);
  ergo_real CI_energy = resultEig_lo;

  do_output(LOG_CAT_INFO, LOG_AREA_CI, "FINAL CI ENERGY = %22.11f", (double)CI_energy);
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "FINAL CI CORRELATION ENERGY = %22.11f", (double)(CI_energy - HF_energy));

  // Now we have the final wavefunction vector in coeffList_lo. Use it to get the final CI density matrix.
  ergo_real* FinalDensityMatrix = new ergo_real[n*n];
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      FinalDensityMatrix[i*n+j] = 0;

  // FIXME: take care of memory leaks here (and elsewhere in CI code also)
  four_idx_SO_struct* dmat2e = new four_idx_SO_struct;
  memset(dmat2e, 0, sizeof(four_idx_SO_struct));
  two_idx_SO_struct* dmat1e = new two_idx_SO_struct;
  memset(dmat1e, 0, sizeof(two_idx_SO_struct));
  pair_status_struct status;
  memset(&status, 0, sizeof(pair_status_struct));
  get_relevant_SlaterDet_pairs_recursive_do_mult_or_dmat(groupList,
							 nElTot,
							 0,
							 0,
							 0,
							 0, 0, 0, 0,
							 &status,
							 nSOs,
							 nSlaterDets,
							 SlaterDetList,
							 g_SO,
							 h_SO,
							 coeffList_lo,
							 NULL,
							 dmat1e,
							 dmat2e
							 );
  for(p = 0; p < nSOs; p++)
    for(q = 0; q < nSOs; q++) {
      SO_struct & currSO_p = SOs[p];
      SO_struct & currSO_q = SOs[q];
      ergo_real* orbitalCoeffs_p = currSO_p.coeffs;
      ergo_real* orbitalCoeffs_q = currSO_q.coeffs;
      for(int a = 0; a < n; a++)
	for(int b = 0; b < n; b++)
	  FinalDensityMatrix[a*n+b] += dmat1e->x[p][q] * orbitalCoeffs_p[a] * orbitalCoeffs_q[b];
    }

  // Check FinalDensityMatrix by computing trace(D*S) which should be the number of electrons
  ergo_real sum = 0;
  for(int a = 0; a < n; a++)
    for(int b = 0; b < n; b++)
      sum += FinalDensityMatrix[a*n+b] * S[a*n+b];
  do_output(LOG_CAT_INFO, LOG_AREA_CI, "CI: trace(D*S) for FinalDensityMatrix: %22.11f", (double)sum);

  // Get dipole moment
  get_dipole_moment_fullmat(n, FinalDensityMatrix, basisInfo, molecule, LOG_AREA_CI, "CI");

  return 0;
}

