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

/** @file integrals_2el_boxed.cc

    @brief Code for 2-electron integrals, computation of Coulomb (J)
    and HF exchange (K) matrices using a single box.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <string.h>

#include "integrals_2el_boxed.h"
#include "integrals_2el_utils.h"
#include "organize_distrs.h"
#include "pi.h"
#include "utilities.h"


static const int HUGE_INTEGER_NUMBER = 2000000000;



typedef struct
{
  int a, b, c, d;
  int poly_ab_index;
  int poly_cd_index;
  int idx1;
  int idx2;
  ergo_real densValue;
} abcd_struct;

#define set_abcd_list_item_macro(i,A,B,C,D,v,i1,i2)			\
  list[i].a = A; list[i].b = B; list[i].c = C; list[i].d = D; list[i].densValue = v; list[i].idx1 = i1; list[i].idx2 = i2; 



static int 
get_JK_contribs_from_2_interacting_boxes(const BasisInfoStruct & basisInfo,
					 const IntegralInfo & integralInfo,
					 int maxNoOfMonomials,
					 ergo_real* J,
					 ergo_real* K,
					 const ergo_real* dens,
					 const minimal_distr_struct* minimalDistrList_1,
					 int noOfGroups_1,
					 const distr_group_struct* groupList_1,
					 const minimal_distr_struct* minimalDistrList_2,
					 int noOfGroups_2,
					 const distr_group_struct* groupList_2,
					 const cluster_struct* clusterList_1,
					 int nClusters_1,
					 const cluster_struct* clusterList_2,
					 int nClusters_2,
					 const batch_struct* batchList_1,
					 int nBatchs_1,
					 const batch_struct* batchList_2,
					 int nBatchs_2,
					 const basis_func_pair_struct* basisFuncPairList_1,
					 const basis_func_pair_struct* basisFuncPairList_2,
					 int interactionWithSelf,
					 ergo_real threshold,
					 JK_contribs_buffer_struct* bufferStructPtr)
{
  int n = basisInfo.noOfBasisFuncs;

  const JK::ExchWeights CAM_params_not_used;

  const ergo_real twoTimesPiToPow5half = 2 * pitopow52;// = 2 * pow(pi, 2.5);
  ergo_real* summedIntegralList = bufferStructPtr->summedIntegralList;
  ergo_real* primitiveIntegralList = bufferStructPtr->primitiveIntegralList;
  ergo_real* primitiveIntegralList_work = bufferStructPtr->primitiveIntegralList_work;

  for(int batch_i = 0; batch_i < nBatchs_1; batch_i++)
    {
      int batch_j_start = 0;
      if(interactionWithSelf == 1)
	batch_j_start = batch_i;
      for(int batch_j = batch_j_start; batch_j < nBatchs_2; batch_j++)
	{
	  int noOfBasisFuncPairs_1 = batchList_1[batch_i].noOfBasisFuncPairs;
	  int noOfBasisFuncPairs_2 = batchList_2[batch_j].noOfBasisFuncPairs;
	  // set integral list to zero
	  memset(summedIntegralList, 0, noOfBasisFuncPairs_1*noOfBasisFuncPairs_2*sizeof(ergo_real));

	  // get largest dmat element
	  ergo_real maxabsdmatelement = 0;
	  for(int i = 0; i < noOfBasisFuncPairs_1; i++)
	    for(int j = 0; j < noOfBasisFuncPairs_2; j++)
	      {
		int a = basisFuncPairList_1[batchList_1[batch_i].basisFuncPairListIndex+i].index_1;
		int b = basisFuncPairList_1[batchList_1[batch_i].basisFuncPairListIndex+i].index_2;
		int c = basisFuncPairList_2[batchList_2[batch_j].basisFuncPairListIndex+j].index_1;
		int d = basisFuncPairList_2[batchList_2[batch_j].basisFuncPairListIndex+j].index_2;
		ergo_real absval;
		if(J != NULL)
		  {
		    absval = template_blas_fabs(dens[a*n+b]);
		    if(absval > maxabsdmatelement)
		      maxabsdmatelement = absval;
		    absval = template_blas_fabs(dens[c*n+d]);
		    if(absval > maxabsdmatelement)
		      maxabsdmatelement = absval;
		  }
		if(K != NULL)
		  {
		    absval = template_blas_fabs(dens[a*n+c]);
		    if(absval > maxabsdmatelement)
		      maxabsdmatelement = absval;
		    absval = template_blas_fabs(dens[a*n+d]);
		    if(absval > maxabsdmatelement)
		      maxabsdmatelement = absval;
		    absval = template_blas_fabs(dens[b*n+c]);
		    if(absval > maxabsdmatelement)
		      maxabsdmatelement = absval;
		    absval = template_blas_fabs(dens[b*n+d]);
		    if(absval > maxabsdmatelement)
		      maxabsdmatelement = absval;
		  }
	      } // END FOR i j get largest dmat element

	  int cluster_i_start = batchList_1[batch_i].clusterStartIndex;
	  int clusterCount1 = batchList_1[batch_i].noOfClusters;
	  for(int cluster_i = cluster_i_start; cluster_i < cluster_i_start + clusterCount1; cluster_i++)
	    {
	      int cluster_j_start = batchList_2[batch_j].clusterStartIndex;
	      int clusterCount2 = batchList_2[batch_j].noOfClusters;
	      int cluterIndexEnd2 = cluster_j_start + clusterCount2;
	      if(interactionWithSelf == 1 && batch_i == batch_j)
		cluster_j_start = cluster_i;
	      for(int cluster_j = cluster_j_start; cluster_j < cluterIndexEnd2; cluster_j++)
		{
		  // check if we can skip this combination of clusters
		  if(clusterList_1[cluster_i].maxLimitingFactorForCluster * clusterList_2[cluster_j].maxLimitingFactorForCluster * maxabsdmatelement < threshold)
		    continue;

		  int group_i_start = clusterList_1[cluster_i].groupStartIndex;
		  int group_i_end = group_i_start + clusterList_1[cluster_i].noOfGroups;
		  int group_j_start = clusterList_2[cluster_j].groupStartIndex;
		  int group_j_end = group_j_start + clusterList_2[cluster_j].noOfGroups;

		  int n1max = clusterList_1[cluster_i].nmax;
		  int n2max = clusterList_2[cluster_j].nmax;

		  // Now we can precompute things that depend only on exponents
		  ergo_real alpha_1 = groupList_1[group_i_start].exponent;
		  ergo_real alpha_2 = groupList_2[group_j_start].exponent;
		  ergo_real alphasum = alpha_1 + alpha_2;
		  ergo_real alphaproduct = alpha_1 * alpha_2;
		  ergo_real alpha_0 = alphaproduct / alphasum;

		  ergo_real resultPreFactor = twoTimesPiToPow5half / (alphaproduct*template_blas_sqrt(alphasum));

		  for(int group_i = group_i_start; group_i < group_i_end; group_i++)
		    {
		      if(interactionWithSelf == 1 && batch_i == batch_j && cluster_i == cluster_j)
			group_j_start = group_i;
		      for(int group_j = group_j_start; group_j < group_j_end; group_j++)
			{
			  if(K == NULL)
			    {
			      // Only J is considered; we can use maxAbsDmatElementGroup
			      ergo_real maxabs_1 = groupList_1[group_i].maxAbsDmatElementGroup;
			      ergo_real maxabs_2 = groupList_2[group_j].maxAbsDmatElementGroup;
			      if((groupList_1[group_i].maxLimitingFactorGroup * groupList_2[group_j].maxLimitingFactorGroup * maxabs_1 < threshold) && 
				 (groupList_1[group_i].maxLimitingFactorGroup * groupList_2[group_j].maxLimitingFactorGroup * maxabs_2 < threshold))
				continue;
			    }
			  else
			    {
			      if(groupList_1[group_i].maxLimitingFactorGroup * groupList_2[group_j].maxLimitingFactorGroup * maxabsdmatelement < threshold)
				continue;
			    }

			  // now we can do all integrals needed for this pair of groups
			  ergo_real dx = groupList_2[group_j].centerCoords[0] - groupList_1[group_i].centerCoords[0];
			  ergo_real dy = groupList_2[group_j].centerCoords[1] - groupList_1[group_i].centerCoords[1];
			  ergo_real dz = groupList_2[group_j].centerCoords[2] - groupList_1[group_i].centerCoords[2];

			  // now we have dx dy dz alpha0 alpha1 n1max n2max. Get all integrals for this case.
			  int noOfMonomials_1 = integralInfo.monomial_info.no_of_monomials_list[n1max];
			  int noOfMonomials_2 = integralInfo.monomial_info.no_of_monomials_list[n2max];

			  if(get_related_integrals_h(integralInfo, 
						     CAM_params_not_used,
						     n1max, noOfMonomials_1,
						     n2max, noOfMonomials_2,
						     dx, dy, dz, alpha_1, alpha_2, alpha_0,
						     primitiveIntegralList,
						     primitiveIntegralList_work,
						     resultPreFactor
						     ) != 0)
			    {
			      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_related_integrals");
			      return -1;
			    }

			  int i_start = groupList_1[group_i].startIndex;
			  int i_end = i_start + groupList_1[group_i].distrCount;
			  for(int i = i_start; i < i_end; i++)
			    {
			      int idx_1 = minimalDistrList_1[i].basisFuncPairIndex;
			      int monomialIndex_1 = minimalDistrList_1[i].monomialIndex;
			      int j_start = groupList_2[group_j].startIndex;
			      int j_end = j_start + groupList_2[group_j].distrCount;
			      if(interactionWithSelf == 1 && group_j == group_i && batch_i == batch_j && cluster_i == cluster_j)
				{
				  // take care of case i = j separately
				  ergo_real integralValue = primitiveIntegralList[monomialIndex_1*noOfMonomials_2+monomialIndex_1];
				  ergo_real integralValueCurr = minimalDistrList_1[i].coeff * minimalDistrList_1[i].coeff * integralValue;
				  integralValueCurr *= 0.5;				  
				  summedIntegralList[idx_1*noOfBasisFuncPairs_2 + idx_1] += integralValueCurr;
				  j_start = i+1;
				}
			      for(int j = j_start; j < j_end; j++)
				{
				  int idx_2 = minimalDistrList_2[j].basisFuncPairIndex;
				  int monomialIndex_2 = minimalDistrList_2[j].monomialIndex;
				  ergo_real integralValue = primitiveIntegralList[monomialIndex_1*noOfMonomials_2+monomialIndex_2];
				  ergo_real integralValueCurr = minimalDistrList_1[i].coeff * minimalDistrList_2[j].coeff * integralValue;
				  summedIntegralList[idx_1*noOfBasisFuncPairs_2 + idx_2] += integralValueCurr;
				} // END FOR j 		      			      
			    } // END FOR i

			} // END FOR group_j
		    } // END FOR group_i
		} // END FOR cluster_j
	    } // END FOR cluster_i

	  for(int idx_1 = 0; idx_1 < noOfBasisFuncPairs_1; idx_1++)
	    for(int idx_2 = 0; idx_2 < noOfBasisFuncPairs_2; idx_2++)
	      {
		int a = basisFuncPairList_1[batchList_1[batch_i].basisFuncPairListIndex+idx_1].index_1;
		int b = basisFuncPairList_1[batchList_1[batch_i].basisFuncPairListIndex+idx_1].index_2;
		int c = basisFuncPairList_2[batchList_2[batch_j].basisFuncPairListIndex+idx_2].index_1;
		int d = basisFuncPairList_2[batchList_2[batch_j].basisFuncPairListIndex+idx_2].index_2;
		ergo_real integralValueCurr = summedIntegralList[idx_1*noOfBasisFuncPairs_2 + idx_2];
		
		if(a == c && b == d)
		  integralValueCurr *= 2;

		if(template_blas_fabs(integralValueCurr)*maxabsdmatelement < threshold)
		  continue;

		if(a != b && c != d && a != c && a != d && b != c && b != d)
		  {
		    if(J != NULL)
		      {
			J[a*n+b] += 2 * dens[c*n+d] * integralValueCurr;
			J[c*n+d] += 2 * dens[a*n+b] * integralValueCurr;
		      }
		    if(K != NULL)
		      {
			if(d >= a)
			  K[a*n+d] += -0.5 * dens[b*n+c] * integralValueCurr;
			else
			  K[d*n+a] += -0.5 * dens[b*n+c] * integralValueCurr;
			if(c >= a)
			  K[a*n+c] += -0.5 * dens[b*n+d] * integralValueCurr;
			else
			  K[c*n+a] += -0.5 * dens[b*n+d] * integralValueCurr;
			if(c >= b)
			  K[b*n+c] += -0.5 * dens[a*n+d] * integralValueCurr;
			else
			  K[c*n+b] += -0.5 * dens[a*n+d] * integralValueCurr;
			if(d >= b)
			  K[b*n+d] += -0.5 * dens[c*n+a] * integralValueCurr;
			else
			  K[d*n+b] += -0.5 * dens[c*n+a] * integralValueCurr;
		      }
		  }
		else
		  {
		    abcd_struct list[8];
      
		    /* determine unique configurations */
		    set_abcd_list_item_macro(0, a, b, c, d, 0, 0, 0);
		    set_abcd_list_item_macro(1, a, b, d, c, 0, 0, 0);
		    set_abcd_list_item_macro(2, b, a, c, d, 0, 0, 0);
		    set_abcd_list_item_macro(3, b, a, d, c, 0, 0, 0);
		    set_abcd_list_item_macro(4, c, d, a, b, 0, 0, 0);
		    set_abcd_list_item_macro(5, d, c, a, b, 0, 0, 0);
		    set_abcd_list_item_macro(6, c, d, b, a, 0, 0, 0);
		    set_abcd_list_item_macro(7, d, c, b, a, 0, 0, 0);

		    int ccc = 0;
	  
		    for(int ii = 0; ii < 8; ii++)
		      {
			abcd_struct* abcd = &list[ii];
			int aa, bb, cc, dd;

			/* check if this is a new unique configuration */
			int unique = 1;
			for(int jj = 0; jj < ii; jj++)
			  {
			    if(abcd->a == list[jj].a && 
			       abcd->b == list[jj].b && 
			       abcd->c == list[jj].c && 
			       abcd->d == list[jj].d)
			      unique = 0;
			  }
			if(unique == 0)
			  continue;
			/* now we know that this configuration is unique. */
			aa = abcd->a;
			bb = abcd->b;
			cc = abcd->c;
			dd = abcd->d;

			ccc++;

			/* add contribution to coulomb matrix */
			if(bb >= aa && J != NULL)
			  J[aa*n+bb] += dens[cc*n+dd] * integralValueCurr;

			if(dd >= aa && K != NULL)
			  K[aa*n+dd] += -0.5 * dens[bb*n+cc] * integralValueCurr;
	      
		      } /* END FOR ii go through 8 configurations */
		  }

	      } // END FOR idx_1 idx_2
	} // END FOR batch_j
    } // END FOR batch_i

  return 0;
}





typedef struct
{
  int id;
  ergo_real x[3];
} point_3d_struct;





int
compute_JK_single_box(const BasisInfoStruct & basisInfo,
		      const IntegralInfo & integralInfo,
		      ergo_real* J,
		      ergo_real* K,
		      const ergo_real* dens,
		      ergo_real threshold)
{
  Util::TimeMeter timeMeterTot;

  Util::TimeMeter timeMeterDistrList;

  int n = basisInfo.noOfBasisFuncs;
  ergo_real maxDensityMatrixElement = get_max_abs_vector_element(n*n, dens);

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "entering compute_JK_single_box, no of basis funcs = %5i, threshold = %7.3g", 
	    n, (double)threshold);

  // Require that threshold value is positive.
  if(threshold <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_JK_single_box: (threshold <= 0)");
      return -1;
    }


  // get largest limiting factor
  Util::TimeMeter timeMeterTmp1;
  ergo_real maxLimitingFactor = 0;
  if(get_list_of_labeled_distrs_maxLimitingFactor(basisInfo,
						  integralInfo,
						  threshold,
						  &maxLimitingFactor,
						  maxDensityMatrixElement) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_list_of_labeled_distrs_maxLimitingFactor");
      return -1;
    }
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "get_list_of_labeled_distrs_maxLimitingFactor done, maxLimitingFactor = %22.11f", 
	    (double)maxLimitingFactor);
  timeMeterTmp1.print(LOG_AREA_INTEGRALS, "get_list_of_labeled_distrs_maxLimitingFactor");
  
  // Get number of distributions
  Util::TimeMeter timeMeterTmp2;
  int distrCount = get_list_of_labeled_distrs(basisInfo,
					      integralInfo,
					      threshold,
					      NULL,
					      0,
					      maxLimitingFactor,
					      dens,
					      maxDensityMatrixElement);
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "getting distrCount done, distrCount = %12i", distrCount);
  timeMeterTmp2.print(LOG_AREA_INTEGRALS, "get_list_of_labeled_distrs for getting distrCount");
  if(distrCount == 0)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "compute_JK_single_box: (distrCount == 0), skipping.");
      memset(J, 0, n*n*sizeof(ergo_real));
      memset(K, 0, n*n*sizeof(ergo_real));
      return 0;
    }
  if(distrCount <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_JK_single_box: (distrCount <= 0)");
      return -1;
    }

  std::vector<DistributionSpecStructLabeled> distrList(distrCount);

  // create list of product primitives, with labels
  Util::TimeMeter timeMeterTmp3;
  int distrCountTemp = get_list_of_labeled_distrs(basisInfo,
						  integralInfo,
						  threshold,
						  &distrList[0],
						  distrCount,
						  maxLimitingFactor,
						  dens,
						  maxDensityMatrixElement);
  if(distrCountTemp != distrCount)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_JK_single_box:(distrCountTemp != distrCount)");
      return -1;
    }
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "get_list_of_labeled_distrs done, distrCount = %12i", distrCount);
  timeMeterTmp3.print(LOG_AREA_INTEGRALS, "get_list_of_labeled_distrs");


  // compute extent for all distrs
  Util::TimeMeter timeMeterComputeExtentForAllDistrs;
  compute_extent_for_list_of_distributions(distrCount, 
					   &distrList[0], 
					   threshold,
					   maxLimitingFactor,
					   maxDensityMatrixElement);
  timeMeterComputeExtentForAllDistrs.print(LOG_AREA_INTEGRALS, "Compute extent for all distrs");

  // get maximum number of monomials
  int maxNoOfMonomials = 0;
  for(int i = 0; i < distrCount; i++)
    {
      int degree = 0;
      for(int j = 0; j < 3; j++)
	degree += distrList[i].distr.monomialInts[j];
      int noOfMonomials = integralInfo.monomial_info.no_of_monomials_list[degree];
      if(noOfMonomials > maxNoOfMonomials)
	maxNoOfMonomials = noOfMonomials;
    } // END FOR ABcount

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "distrCount = %i", distrCount);

  std::vector<DistributionSpecStructLabeled> distrList2(distrCount);
  int jcounter = 0;
  for(int i = 0; i < distrCount; i++)
    {
      distrList2[jcounter] = distrList[i];
      jcounter++;
    }
  distrCount = jcounter;
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "distrCount = %i (after removing negligible products)", distrCount);

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "Creating list of distributions done, distrCount = %9i", distrCount);
  timeMeterDistrList.print(LOG_AREA_INTEGRALS, "Creating list of distributions");

#define NUMBER_OF_PARTS 1

  int n_list[NUMBER_OF_PARTS];

  distr_list_description_struct distr_list_description_list[NUMBER_OF_PARTS];

  for(int i = 0; i < NUMBER_OF_PARTS; i++)
    n_list[i] = 0;

  for(int i = 0; i < distrCount; i++)
    n_list[i % NUMBER_OF_PARTS]++;

  ergo_real centerCoords[3];
  memset(centerCoords, 0, 3*sizeof(ergo_real));
  int currIndex = 0;
  for(int i = 0; i < NUMBER_OF_PARTS; i++)
    {
      if(organize_distributions(integralInfo, 
				&distrList2[currIndex], 
				n_list[i], 
				&distr_list_description_list[i].org,
				centerCoords,
				HUGE_INTEGER_NUMBER) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in organize_distributions");
	  return -1;
	}
      currIndex += n_list[i];
    }


  // Set J to zero
  memset(J, 0, n*n*sizeof(ergo_real));

  // Set K to zero
  memset(K, 0, n*n*sizeof(ergo_real));


  // Allocate buffers needed by integral code
  JK_contribs_buffer_struct bufferStruct;
  allocate_buffers_needed_by_integral_code(integralInfo, maxNoOfMonomials, 0, &bufferStruct);

  for(int i = 0; i < NUMBER_OF_PARTS; i++)
    for(int j = i; j < NUMBER_OF_PARTS; j++)
      {
	int self = 0;
	if(i == j)
	  self = 1;
	Util::TimeMeter timeMeterJKcontribs;
	if(get_JK_contribs_from_2_interacting_boxes(basisInfo,
						    integralInfo,
						    maxNoOfMonomials,
						    J,
						    K,
						    dens,

						    &distr_list_description_list[i].org.minimalDistrList[0],
						    distr_list_description_list[i].org.groupList.size(),
						    &distr_list_description_list[i].org.groupList[0],

						    &distr_list_description_list[j].org.minimalDistrList[0],
						    distr_list_description_list[j].org.groupList.size(),
						    &distr_list_description_list[j].org.groupList[0],

						    &distr_list_description_list[i].org.clusterList[0],
						    distr_list_description_list[i].org.clusterList.size(),

						    &distr_list_description_list[j].org.clusterList[0],
						    distr_list_description_list[j].org.clusterList.size(),

						    &distr_list_description_list[i].org.batchList[0],
						    distr_list_description_list[i].org.batchList.size(),

						    &distr_list_description_list[j].org.batchList[0],
						    distr_list_description_list[j].org.batchList.size(),

						    &distr_list_description_list[i].org.basisFuncPairList[0],
						    &distr_list_description_list[j].org.basisFuncPairList[0],

						    self,
						    threshold,
						    &bufferStruct) != 0)
	  {
	    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_JK_contribs_from_2_interacting_boxes");
	    return -1;
	  }
	timeMeterJKcontribs.print(LOG_AREA_INTEGRALS, "get_JK_contribs_from_2_interacting_boxes for both J and K together");
      } // END FOR i j

  // Fill the other triangle of K
  for(int i = 0; i < n; i++)
    for(int j = 0; j < i; j++)
      K[i*n+j] = K[j*n+i];

  // Fill the other triangle of J
  for(int i = 0; i < n; i++)
    for(int j = 0; j < i; j++)
      J[i*n+j] = J[j*n+i];

  free_buffers_needed_by_integral_code(&bufferStruct);

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "compute_JK_single_box ending OK.");
  timeMeterTot.print(LOG_AREA_INTEGRALS, "compute_JK_single_box");

  return 0;
}
