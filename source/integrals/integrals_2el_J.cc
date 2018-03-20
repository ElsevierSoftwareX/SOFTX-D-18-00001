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

/** @file integrals_2el_J.cc

    \brief Code for computing the Coulomb matrix J.

    @author: Elias Rudberg <em>responsible</em>.
*/

#include <cstring>
#include <cstdio>
#include <cassert>

#include <pthread.h>

#include "integrals_2el_J_kernel.h"
#include "integrals_2el_J.h"
#include "integrals_2el_utils.h"
#include "mm_limit_table.h"
#include "basis_func_pair_list.h"
#include "integrals_2el_repeating.h"
#include "integrals_hermite.h"
#include "integrals_general.h"
#include "utilities.h"
#include "pi.h"
#include "integrals_2el_util_funcs.h"
#include "integrals_2el_J_mm_utils.h"
#include "integrals_2el_J_mm_kernel.h"


static const int HUGE_INTEGER_NUMBER = 2000000000;


typedef struct
{
  int boxIndex_1;
  int boxIndex_2;
  int branchIndex_1;
  int branchIndex_2;
} job_list_standard_entry_J_struct;

typedef struct
{
  int boxIndex;
  int multipoleBoxIndex;
  short int branchIndex;
  short int multipoleBranchIndex;
} job_list_multipole_entry_J_struct;


static int
add_multipole_jobs_for_2_boxes_branches_recursive(int multipoleBoxIndex,
						  int multipoleBranchIndex,
						  const box_struct* boxList,
						  int boxIndex,
						  int branchIndex,
						  int numberOfLevels,
						  int currLevel,
						  job_list_multipole_entry_J_struct* jobList_multipole,
						  int maxNoOfJobs_multipole
						  )
{
  if(currLevel == numberOfLevels - 1)
    {
      if(maxNoOfJobs_multipole <= 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in add_multipole_jobs_for_2_boxes_branches_recursive: (maxNoOfJobs_multipole <= 0)");
	  return -1;
	}
      if(jobList_multipole != NULL)
	{
	  jobList_multipole[0].boxIndex = boxIndex;
	  jobList_multipole[0].branchIndex = branchIndex;
	  jobList_multipole[0].multipoleBoxIndex = multipoleBoxIndex;
	  jobList_multipole[0].multipoleBranchIndex = multipoleBranchIndex;
	}
      return 1;
    }
  // go through children
  int noOfChildren = boxList[boxIndex].basicBox.noOfChildBoxes;
  int noOfNewJobs = 0;
  for(int i = 0; i < noOfChildren; i++)
    {
      int childIndex = boxList[boxIndex].basicBox.firstChildBoxIndex + i;
      job_list_multipole_entry_J_struct* jobListPtr = NULL;      
      if(jobList_multipole != NULL)
	jobListPtr = &jobList_multipole[noOfNewJobs];
      int nJobs = add_multipole_jobs_for_2_boxes_branches_recursive(multipoleBoxIndex,
								    multipoleBranchIndex,
								    boxList,
								    childIndex,
								    branchIndex,
								    numberOfLevels,
								    currLevel + 1,
								    jobListPtr,
								    maxNoOfJobs_multipole - noOfNewJobs
								    );
	if(nJobs < 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in add_multipole_jobs_for_2_boxes_branches_recursive");
	  return -1;
	}
      noOfNewJobs += nJobs;
    }  
  return noOfNewJobs;
}




static int
get_joblists_J_for_two_boxes_recursive(const IntegralInfo & integralInfo,
				       ergo_real threshold,
				       const box_struct* boxList,
				       int numberOfLevels,
				       int currLevel,
				       int boxIndex_1,
				       int boxIndex_2,
				       int branchIndex_1,
				       int branchIndex_2,

				       job_list_standard_entry_J_struct* jobList_standard,
				       int maxNoOfJobs_standard,
				       int* noOfNewJobs_standard,

				       job_list_multipole_entry_J_struct* jobList_multipole,
				       int maxNoOfJobs_multipole,
				       int* noOfNewJobs_multipole
				       )
{
  // Both boxes must have same width
  assert(template_blas_fabs(boxList[boxIndex_1].basicBox.width - boxList[boxIndex_2].basicBox.width) < 1e-4);

  int useMultipoleDescription = check_if_multipoles_can_be_used(integralInfo,
								threshold,
								boxList[boxIndex_1].basicBox.centerCoords,
								boxList[boxIndex_2].basicBox.centerCoords,
								boxList[boxIndex_1].basicBox.width,
								boxList[boxIndex_1].branchListForJ[branchIndex_1].org,
								boxList[boxIndex_1].branchListForJ[branchIndex_1].org_mm,
								boxList[boxIndex_2].branchListForJ[branchIndex_2].org,
								boxList[boxIndex_2].branchListForJ[branchIndex_2].org_mm);
  if(useMultipoleDescription < 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_joblists_J_for_two_boxes_recursive: (useMultipoleDescription < 0).");
    return -1;
  }

  if(useMultipoleDescription == 1)
    {
      // Use multipole description
      int noOfNewJobs_1 = add_multipole_jobs_for_2_boxes_branches_recursive(boxIndex_2,
									    branchIndex_2,
									    boxList,
									    boxIndex_1,
									    branchIndex_1,
									    numberOfLevels,
									    currLevel,
									    jobList_multipole,
									    maxNoOfJobs_multipole
									    );
      job_list_multipole_entry_J_struct* secondPtr = NULL;
      if(jobList_multipole != NULL)
	secondPtr = &jobList_multipole[noOfNewJobs_1];
      int noOfNewJobs_2 = add_multipole_jobs_for_2_boxes_branches_recursive(boxIndex_1,
									    branchIndex_1,
									    boxList,
									    boxIndex_2,
									    branchIndex_2,
									    numberOfLevels,
									    currLevel,
									    secondPtr,
									    maxNoOfJobs_multipole - noOfNewJobs_1
									    );
      if(noOfNewJobs_1 < 0 || noOfNewJobs_2 < 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in add_multipole_jobs_for_2_boxes_branches_recursive");
	  return -1;
	}
      *noOfNewJobs_standard = 0;
      *noOfNewJobs_multipole = noOfNewJobs_1 + noOfNewJobs_2;

      return 0;
    }

  // Multipoles could not be used. We must either go to the next level or compute integrals explicitly.
  if(currLevel == numberOfLevels-1)
    {
      // We are at the level of smallest boxes. Add standard job to job list.

      if(boxIndex_1 == boxIndex_2 && branchIndex_1 > branchIndex_2)
	{
	  *noOfNewJobs_standard = 0;
	  *noOfNewJobs_multipole = 0;
	  return 0;
	}
      if(maxNoOfJobs_standard <= 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_joblists_J_for_two_boxes_recursive: (maxNoOfJobs_standard <= 0)");
	  return -1;
	}
      if(jobList_standard != NULL)
	{
	  jobList_standard[0].boxIndex_1 = boxIndex_1;
	  jobList_standard[0].branchIndex_1 = branchIndex_1;
	  jobList_standard[0].boxIndex_2 = boxIndex_2;
	  jobList_standard[0].branchIndex_2 = branchIndex_2;
	}
      *noOfNewJobs_standard = 1;
      *noOfNewJobs_multipole = 0;
      return 0;      
    }

  // Go to next level. Do interaction between all pairs of children of the two boxes.
  int noOfChildren_1 = boxList[boxIndex_1].basicBox.noOfChildBoxes;
  int noOfChildren_2 = boxList[boxIndex_2].basicBox.noOfChildBoxes;

  if(noOfChildren_1 <= 0 || noOfChildren_2 <= 0)
    exit(EXIT_FAILURE);

  int noOfNewJobs_standard_count = 0;
  int noOfNewJobs_multipole_count = 0;
  for(int i = 0; i < noOfChildren_1; i++)
    {
      int start_j = 0;
      if(boxIndex_1 == boxIndex_2)
	start_j = i;
      for(int j = start_j; j < noOfChildren_2; j++)
	{
	  int childIndex_1 = boxList[boxIndex_1].basicBox.firstChildBoxIndex + i;
	  int childIndex_2 = boxList[boxIndex_2].basicBox.firstChildBoxIndex + j;
	  int nJobs_standard = 0;
	  int nJobs_multipole = 0;
	  job_list_multipole_entry_J_struct* jobList_multipole_mod = NULL;
	  if(jobList_multipole != NULL)
	    jobList_multipole_mod = &jobList_multipole[noOfNewJobs_multipole_count];
	  job_list_standard_entry_J_struct* jobList_standard_mod = NULL;
	  if(jobList_standard != NULL)
	    jobList_standard_mod = &jobList_standard[noOfNewJobs_standard_count];
	  if(get_joblists_J_for_two_boxes_recursive(integralInfo,
						    threshold,
						    boxList,
						    numberOfLevels,
						    currLevel + 1,
						    childIndex_1,
						    childIndex_2,
						    branchIndex_1,
						    branchIndex_2,

						    jobList_standard_mod,
						    maxNoOfJobs_standard - noOfNewJobs_standard_count,
						    &nJobs_standard,

						    jobList_multipole_mod,
						    maxNoOfJobs_multipole - noOfNewJobs_multipole_count,
						    &nJobs_multipole
						    ) != 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_joblists_J_for_two_boxes_recursive for child boxes");
	      return -1;
	    }
	  noOfNewJobs_standard_count += nJobs_standard;
	  noOfNewJobs_multipole_count += nJobs_multipole;
	} // END FOR j
    } // END FOR i
  
  *noOfNewJobs_standard = noOfNewJobs_standard_count;
  *noOfNewJobs_multipole = noOfNewJobs_multipole_count;
  return 0;
}




static int
get_list_of_labeled_distrs_maxLimitingFactor_linear(const BasisInfoStruct & basisInfo,
						    const IntegralInfo & integralInfo,
						    ergo_real threshold,
						    const basis_func_index_pair_struct* basisFuncIndexPairList,
						    int basisFuncIndexPairCount,
						    ergo_real* resultMaxLimitingFactor)
{ 
  IntegratorWithMemory integrator(&integralInfo);

  ergo_real maxLimitingFactor = 0;
  for(int kk = 0; kk < basisFuncIndexPairCount; kk++)
    {
      int i = basisFuncIndexPairList[kk].index_1;
      int j = basisFuncIndexPairList[kk].index_2;
	  
      const int maxCountProduct = POLY_PRODUCT_MAX_DISTRS;
      DistributionSpecStruct psi_list[maxCountProduct];
      /* form product of basisfuncs i and j, store product in psi_list */
      int n_psi = get_product_simple_primitives(basisInfo, i,
						basisInfo, j,
						psi_list,
						maxCountProduct,
						0);
      if(n_psi < 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_product_simple_primitives");
	  return -1;
	}
      for(int k = 0; k < n_psi; k++)
	{
	  ergo_real limitingFactor = template_blas_sqrt(integrator.do_2e_integral(&psi_list[k]));
	  if(limitingFactor > maxLimitingFactor)
	    maxLimitingFactor = limitingFactor;
	} // END FOR k
    } // END FOR kk
  *resultMaxLimitingFactor = maxLimitingFactor;

  return 0;
}




static int
get_list_of_labeled_distrs_linear(const BasisInfoStruct & basisInfo,
				  const IntegralInfo & integralInfo,
				  ergo_real threshold,
				  DistributionSpecStructLabeled* resultList,
				  int maxCountDistrs,
				  ergo_real maxLimitingFactor,
				  const basis_func_index_pair_struct* basisFuncIndexPairList,
				  int basisFuncIndexPairCount,
				  const ergo_real* D_list)
{
  ergo_real maxDensityMatrixElement = get_max_abs_vector_element(basisFuncIndexPairCount, D_list);

  IntegratorWithMemory integrator(&integralInfo);

  // create list of product primitives, with labels
  int distrCount = 0;
  for(int kk = 0; kk < basisFuncIndexPairCount; kk++)
    {
      int i = basisFuncIndexPairList[kk].index_1;
      int j = basisFuncIndexPairList[kk].index_2;

      ergo_real dmatElement = D_list[kk];

      const int maxCountProduct = POLY_PRODUCT_MAX_DISTRS;
      DistributionSpecStruct psi_list[maxCountProduct];
      /* form product of basisfuncs i and j, store product in psi_list */
      int n_psi = get_product_simple_primitives(basisInfo, i,
						basisInfo, j,
						psi_list,
						maxCountProduct,
						0);
      if(n_psi < 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_product_simple_primitives");
	  return -1;
	}
      for(int k = 0; k < n_psi; k++)
	{
	  ergo_real limitingFactor = template_blas_sqrt(integrator.do_2e_integral(&psi_list[k]));
	  if(limitingFactor*maxLimitingFactor*maxDensityMatrixElement > threshold)
	    {
	      if(maxCountDistrs > 0 && distrCount >= maxCountDistrs)
		{
		  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_list_of_labeled_distrs: (maxCountDistrs > 0 && distrCount >= maxCountDistrs)");
		  return -1;
		}
	      if(resultList != NULL)
		{
		  resultList[distrCount].distr = psi_list[k];
		  resultList[distrCount].basisFuncIndex_1 = i;
		  resultList[distrCount].basisFuncIndex_2 = j;
		  resultList[distrCount].pairIndex = kk;
		  resultList[distrCount].limitingFactor = limitingFactor;
		  resultList[distrCount].dmatElement = dmatElement;
		}
	      distrCount++;
	    } // END IF above threshold
	} // END FOR k
    } // END FOR kk

  return distrCount;
}





static int
compare_multipole_jobs(const void* p1, const void* p2)
{
  job_list_multipole_entry_J_struct* job_1 = (job_list_multipole_entry_J_struct*)p1;
  job_list_multipole_entry_J_struct* job_2 = (job_list_multipole_entry_J_struct*)p2;
  if(job_1->boxIndex > job_2->boxIndex)
    return 1;
  if(job_1->boxIndex < job_2->boxIndex)
    return -1;
  // now we know that boxIndex is the same for both
  if(job_1->branchIndex > job_2->branchIndex)
    return 1;
  if(job_1->branchIndex < job_2->branchIndex)
    return -1;
  // now we know that boxIndex and branchIndex are the same for both
  if(job_1->multipoleBoxIndex > job_2->multipoleBoxIndex)
    return 1;
  if(job_1->multipoleBoxIndex < job_2->multipoleBoxIndex)
    return -1;
  // now we know that boxIndex and branchIndex and multipoleBoxIndex are the same for both
  // we do not care about the order of different multipoleBranchIndex
  return 0;
}





static void
get_largest_and_smallest_extent_for_list_of_distributions(int n, 
							  const DistributionSpecStructLabeled* distrList, 
							  ergo_real* result_extent_min, 
							  ergo_real* result_extent_max)
{
  ergo_real extent_min = distrList[0].distr.extent;
  ergo_real extent_max = distrList[0].distr.extent;
  for(int i = 0; i < n; i++)
    {
      ergo_real extent = distrList[i].distr.extent;
      if(extent > extent_max)
	extent_max = extent;
      if(extent < extent_min)
	extent_min = extent;
    }
  *result_extent_min = extent_min;
  *result_extent_max = extent_max;
}



static int
get_branch_splitter_info(ergo_real* branchSplitterList,
			 int maxNoOfBranches,
			 const JK::Params& J_K_params,
			 ergo_real toplevelBoxSize,
			 ergo_real extent_max)
{
  int noOfBranches = 0;
  if(J_K_params.fmm_no_of_branches > 0)
    {
      // Use branches as specified in input parameters
      noOfBranches = J_K_params.fmm_no_of_branches;
      if(noOfBranches >= maxNoOfBranches)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_branch_splitter_info: (noOfBranches >= maxNoOfBranches)");
	  return -1;
	}
      for(int i = 0; i < noOfBranches-1; i++)
	{
	  ergo_real splitterValue = 0;
	  switch(i)
	    {
	    case 0: splitterValue = J_K_params.fmm_branch_splitter_extent_1; break;
	    case 1: splitterValue = J_K_params.fmm_branch_splitter_extent_2; break;
	    case 2: splitterValue = J_K_params.fmm_branch_splitter_extent_3; break;
	    case 3: splitterValue = J_K_params.fmm_branch_splitter_extent_4; break;
	    case 4: splitterValue = J_K_params.fmm_branch_splitter_extent_5; break;
	    default:
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_branch_splitter_info: default reached.");
	      return -1;
	    }
	  branchSplitterList[i] = splitterValue;
	}  
    }
  else
    {
      // Use default branch settings based on box size and extent_max
      ergo_real splitterValue = toplevelBoxSize / 2;
      noOfBranches = 2;
      while(splitterValue < extent_max)
	{
	  noOfBranches++;
	  splitterValue *= 2;
	}
      // Now we know how many branches we need. Create splitter list.
      int count = 0;
      splitterValue = 0;
      branchSplitterList[noOfBranches-2-count] = splitterValue;
      count++;
      splitterValue = toplevelBoxSize / 2;
      for(count = 1; count < noOfBranches-1; count++)
	{
	  branchSplitterList[noOfBranches-2-count] = splitterValue;
	  splitterValue *= 2;
	}
    }
  return noOfBranches;
}









static int
create_branches(int noOfBranches,
		const ergo_real* branchSplitterList,
		int distrCount,
		DistributionSpecStructLabeled* distrListOrdered,
		int noOfBoxesTopLevel,
		box_struct* boxListTopLevel
		)
{
  // Start by finding out largest number of distrs per box.
  int maxNoOfDistrsPerBox = 0;
  for(int i = 0; i < noOfBoxesTopLevel; i++) {
    int distrCountCurrBox = boxListTopLevel[i].basicBox.noOfItems;
    if(distrCountCurrBox > maxNoOfDistrsPerBox)
      maxNoOfDistrsPerBox = distrCountCurrBox;
  }
  
  std::vector<int> branchBucketIndexList[MAX_NO_OF_BRANCHES];
  int branchBucketCountList[MAX_NO_OF_BRANCHES];
  for(int i = 0; i < noOfBranches; i++)
    branchBucketIndexList[i].resize(maxNoOfDistrsPerBox);
  
  std::vector<DistributionSpecStructLabeled> distrListTemp(maxNoOfDistrsPerBox);
  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating distrListTemp");

  for(int i = 0; i < noOfBoxesTopLevel; i++) {
    DistributionSpecStructLabeled* distrListCurrBox  = &distrListOrdered[boxListTopLevel[i].basicBox.firstItemIndex];
    int distrCountCurrBox = boxListTopLevel[i].basicBox.noOfItems;
    memcpy(&distrListTemp[0], distrListCurrBox, distrCountCurrBox*sizeof(DistributionSpecStructLabeled));
    DistributionSpecStructLabeled* distrListCurrBox2 = &distrListTemp[0];
    for(int branchIndex = 0; branchIndex < noOfBranches; branchIndex++)
      branchBucketCountList[branchIndex] = 0;
    for(int j = 0; j < distrCountCurrBox; j++) {
      int branchIndex; // declare here because value is used after loop after break.
      for(branchIndex = noOfBranches-1; branchIndex > 0; branchIndex--) {
	ergo_real extent = distrListCurrBox[j].distr.extent;
	ergo_real width = boxListTopLevel[i].basicBox.width;
	// get minWallDist : minimum wall distance
	ergo_real minWallDist = width;
	for(int coordIndex = 0; coordIndex< 3; coordIndex++) {
	  // get wall distance for this coordinate
	  ergo_real dx = distrListCurrBox[j].distr.centerCoords[coordIndex] - boxListTopLevel[i].basicBox.centerCoords[coordIndex];
	  ergo_real wallDist = width - template_blas_fabs(dx);
	  if(wallDist < minWallDist)
	    minWallDist = wallDist;
	} // END FOR coordIndex
	if(minWallDist < 0) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: (minWallDist < 0)");
	  return -1;
	}
	if((extent - minWallDist) < branchSplitterList[branchIndex-1])
	  break;
      }
      branchBucketIndexList[branchIndex][branchBucketCountList[branchIndex]] = j;
      branchBucketCountList[branchIndex]++;
    } // END FOR j
    int newCount = 0;
    for(int branchIndex = 0; branchIndex < noOfBranches; branchIndex++) {
      boxListTopLevel[i].branchIndexListForJ[branchIndex] = boxListTopLevel[i].basicBox.firstItemIndex + newCount;	  
      boxListTopLevel[i].branchCountListForJ[branchIndex] = branchBucketCountList[branchIndex];
      for(int k = 0; k < branchBucketCountList[branchIndex]; k++) {
	distrListCurrBox[newCount] = distrListCurrBox2[branchBucketIndexList[branchIndex][k]];
	newCount++;
      }
    } // END FOR branchIndex
  } // END FOR i divide distrs into branches according to extent.

  return 0;
}









static int
execute_joblist_J_std_serial(int noOfJobs_J_standard,
			     const job_list_standard_entry_J_struct* jobList_J_standard,
			     const IntegralInfo & integralInfo,
			     int maxNoOfMonomials,
			     ergo_real* result_J_list,
			     const box_struct* boxList,
			     ergo_real threshold)
{
  Util::TimeMeter timeMeter;

  JK_contribs_buffer_struct bufferStruct;
  allocate_buffers_needed_by_integral_code(integralInfo, maxNoOfMonomials, 0, &bufferStruct);

  for(int jobIndex = 0; jobIndex < noOfJobs_J_standard; jobIndex++)
    {
      int boxIndex_1 = jobList_J_standard[jobIndex].boxIndex_1;
      int boxIndex_2 = jobList_J_standard[jobIndex].boxIndex_2;
      int branchIndex_1 = jobList_J_standard[jobIndex].branchIndex_1;
      int branchIndex_2 = jobList_J_standard[jobIndex].branchIndex_2;
      int self = 0;
      if(boxIndex_1 == boxIndex_2 && branchIndex_1 == branchIndex_2)
	self = 1;
      if(get_J_contribs_from_2_interacting_boxes(integralInfo,
						 result_J_list,
						 NULL,
						 boxList[boxIndex_1].branchListForJ[branchIndex_1].org,
						 boxList[boxIndex_2].branchListForJ[branchIndex_2].org,
						 self,
						 threshold,
						 &bufferStruct) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_J_contribs_from_2_interacting_boxes");
	  return -1;
	}
    } // END FOR jobIndex

  free_buffers_needed_by_integral_code(&bufferStruct);

  timeMeter.print(LOG_AREA_INTEGRALS, "execute_joblist_J_std_serial");

  return 0;
}


struct J_std_joblist_thread_struct
{
  pthread_t thread;
  const IntegralInfo* integralInfo;
  ergo_real* result_J_list;
  int maxNoOfMonomials;
  ergo_real threshold;
  const box_struct* boxList;
  const job_list_standard_entry_J_struct* jobList_J_standard;
  int noOfJobs_J_standard;
  int thread_ID;
  int noOfThreads;
  int resultCode;
  explicit J_std_joblist_thread_struct() {}
};


static void*
execute_joblist_J_std_thread_func(void* arg)
{
  J_std_joblist_thread_struct* params = (J_std_joblist_thread_struct*)arg;

  JK_contribs_buffer_struct bufferStruct;
  allocate_buffers_needed_by_integral_code(*params->integralInfo, params->maxNoOfMonomials, 0, &bufferStruct);

  const box_struct* boxList = params->boxList;

  for(int jobIndex = 0; jobIndex < params->noOfJobs_J_standard; jobIndex++)
    {
      if(jobIndex % params->noOfThreads != params->thread_ID)
	continue;

      int boxIndex_1 = params->jobList_J_standard[jobIndex].boxIndex_1;
      int boxIndex_2 = params->jobList_J_standard[jobIndex].boxIndex_2;
      int branchIndex_1 = params->jobList_J_standard[jobIndex].branchIndex_1;
      int branchIndex_2 = params->jobList_J_standard[jobIndex].branchIndex_2;
      int self = 0;
      if(boxIndex_1 == boxIndex_2 && branchIndex_1 == branchIndex_2)
	self = 1;
      if(get_J_contribs_from_2_interacting_boxes(*params->integralInfo,
						 params->result_J_list,
						 NULL,
						 boxList[boxIndex_1].branchListForJ[branchIndex_1].org,
						 boxList[boxIndex_2].branchListForJ[branchIndex_2].org,
						 self,
						 params->threshold,
						 &bufferStruct) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_J_contribs_from_2_interacting_boxes");
	  params->resultCode = -1;
	  return NULL;
	}
    } // END FOR jobIndex

  free_buffers_needed_by_integral_code(&bufferStruct);

  params->resultCode = 0;
  return NULL;
}


static int
execute_joblist_J_std_threaded(int noOfThreads,
			       int noOfJobs_J_standard,
			       const job_list_standard_entry_J_struct* jobList_J_standard,
			       const IntegralInfo & integralInfo,
			       int maxNoOfMonomials,
			       ergo_real* result_J_list,
			       int noOfBasisFuncIndexPairs,
			       const box_struct* boxList,
			       ergo_real threshold)
{
  Util::TimeMeter timeMeter;

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "execute_joblist_J_std_threaded, noOfThreads = %2i", noOfThreads);
  
  J_std_joblist_thread_struct* threadParamsList[noOfThreads];
  
  // Set common parameters for all threads
  for(int i = 0; i < noOfThreads; i++)
    {
      threadParamsList[i] = new J_std_joblist_thread_struct();
      threadParamsList[i]->integralInfo = &integralInfo;
      threadParamsList[i]->maxNoOfMonomials = maxNoOfMonomials;
      threadParamsList[i]->boxList = boxList;
      threadParamsList[i]->jobList_J_standard = jobList_J_standard;
      threadParamsList[i]->noOfJobs_J_standard = noOfJobs_J_standard;
      threadParamsList[i]->noOfThreads = noOfThreads;
      threadParamsList[i]->resultCode = -1; // initialize to error code
      threadParamsList[i]->threshold = threshold;
    } // END FOR i
  
  // Set result pointer for thread 0
  // Thread 0 uses the original result_J_list pointer.
  threadParamsList[0]->result_J_list = result_J_list;

  // Set result_J_list pointer for other threads
  for(int i = 1; i < noOfThreads; i++)
    {
      threadParamsList[i]->result_J_list = new ergo_real[noOfBasisFuncIndexPairs];
      memset(threadParamsList[i]->result_J_list, 0, noOfBasisFuncIndexPairs * sizeof(ergo_real));
    }

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating memory for threads.");


  // Set ID number for all threads
  for(int i = 0; i < noOfThreads; i++)
    threadParamsList[i]->thread_ID = i;

  /* start threads */
  for(int i = 0; i < noOfThreads; i++)
    {
      if(pthread_create(&threadParamsList[i]->thread, 
			NULL, 
			execute_joblist_J_std_thread_func, 
			threadParamsList[i]) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in pthread_create for thread %i", i);
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "waiting for already created threads..");
	  for(int j = 0; j < i; j++)
	    {
	      if(pthread_join(threadParamsList[j]->thread, NULL) != 0)
		do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in pthread_join for thread %i", j);
	    } /* END FOR j */
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "all threads finished, returning error code");
	  return -1;
	}
    } /* END FOR i */
  
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "%i threads started OK.", noOfThreads);

  /* wait for threads to finish */
  for(int i = 0; i < noOfThreads; i++)
    {
      if(pthread_join(threadParamsList[i]->thread, NULL) != 0)
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in pthread_join for thread %i", i);
    } /* END FOR i */
  
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "all %i threads have finished.", noOfThreads);
  
  /* now all threads have finished, check for errors */
  for(int i = 0; i < noOfThreads; i++)
    {
      if(threadParamsList[i]->resultCode != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in execute_joblist_J_std_thread_func"
		    " for thread %i", i);
	  return -1;
	}
    } /* END FOR i */
  
  
  // add contributions from other threads
  for(int i = 1; i < noOfThreads; i++)
    {
      for(int j = 0; j < noOfBasisFuncIndexPairs; j++)
	result_J_list[j] += threadParamsList[i]->result_J_list[j];
    }
  
  // Free extra result_J_list buffers used by threads.
  // Note that this loop must start with 1, not 0.
  for(int i = 1; i < noOfThreads; i++)
    delete [] threadParamsList[i]->result_J_list;

  for(int i = 0; i < noOfThreads; i++)
    delete threadParamsList[i];

  timeMeter.print(LOG_AREA_INTEGRALS, "execute_joblist_J_std_threaded");
  
  return 0;
}


static int 
sort_list_of_multipole_jobs_fixed_boxIndex(job_list_multipole_entry_J_struct* jobList, int n)
{
  // Start by bucket-sort by branchIndex.
  const int maxNoOfBranches = 10;  
  job_list_multipole_entry_J_struct* bucketList[maxNoOfBranches];
  // Get number of branches
  int branchIndex_min = maxNoOfBranches;
  int branchIndex_max = 0;
  for(int i = 0; i < n; i++)
    {
      int currBranchIndex = jobList[i].branchIndex;
      if(currBranchIndex > branchIndex_max)
	branchIndex_max = currBranchIndex;
      if(currBranchIndex < branchIndex_min)
	branchIndex_min = currBranchIndex;
    }
  assert(branchIndex_min >= 0);
  assert(branchIndex_max < maxNoOfBranches);
  int noOfBranches = branchIndex_max + 1;
  for(int i = 0; i < noOfBranches; i++)
    bucketList[i] = new job_list_multipole_entry_J_struct[n];

  int counterList[maxNoOfBranches];
  for(int i = 0; i < maxNoOfBranches; i++)
    counterList[i] = 0;
  
  for(int i = 0; i < n; i++)
    {
      int currBranchIndex = jobList[i].branchIndex;
      assert(currBranchIndex < noOfBranches);
      int count = counterList[currBranchIndex];
      bucketList[currBranchIndex][count] = jobList[i];
      counterList[currBranchIndex]++;
    }

  // OK, bucket-sort done. Now sort contents of each bucket.
  int count = 0;
  for(int i = 0; i < maxNoOfBranches; i++)
    {
      int currCount = counterList[i];
      if(currCount == 0) continue;

      qsort(bucketList[i], currCount, sizeof(job_list_multipole_entry_J_struct),
            compare_multipole_jobs);

      // check qsort result
      for(int j  = 0; j < currCount-1; j++)
      {
        job_list_multipole_entry_J_struct* curr = &bucketList[i][j];
        job_list_multipole_entry_J_struct* next = &bucketList[i][j+1];
        if(compare_multipole_jobs(curr, next) > 0)
        {
          do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: qsort result not sorted.");
          return -1;
        }
      }

      // Copy result
      memcpy(&jobList[count], bucketList[i], currCount * sizeof(job_list_multipole_entry_J_struct));
      count += currCount;
    }

  for(int i = 0; i < noOfBranches; i++)
    delete [] bucketList[i];
  
  return 0;
}


static int 
sort_list_of_multipole_jobs(job_list_multipole_entry_J_struct* jobList, int n)
{
  if(n == 0)
    return 0;

  Util::TimeMeter timeMeterInit;

  // Start by bucket-sort by boxIndex.
  // Go through list once to find max boxIndex.
  int boxIndex_max = jobList[0].boxIndex;
  for(int i = 0; i < n; i++)
    {
      int currBoxIndex = jobList[i].boxIndex;
      if(currBoxIndex > boxIndex_max)
	boxIndex_max = currBoxIndex;
    }

  // Go through list once more to find maxNoOfJobsWithSameBoxIndex
  int noOfBoxIndexes = boxIndex_max + 1;
  std::vector<int> counterList1(noOfBoxIndexes);
  for(int i = 0; i < noOfBoxIndexes; i++)
    counterList1[i] = 0;
  for(int i = 0; i < n; i++)
    {
      int currBoxIndex = jobList[i].boxIndex;
      counterList1[currBoxIndex]++;
    }
  int maxNoOfJobsWithSameBoxIndex = 0;
  for(int i = 0; i < noOfBoxIndexes; i++)
    {
      if(counterList1[i] > maxNoOfJobsWithSameBoxIndex)
	maxNoOfJobsWithSameBoxIndex = counterList1[i];
    }

  timeMeterInit.print(LOG_AREA_INTEGRALS, "sort_list_of_multipole_jobs init part");

  //std::vector<job_list_multipole_entry_J_struct> bucketList(noOfBoxIndexes*maxNoOfJobsWithSameBoxIndex);
  std::vector< std::vector<job_list_multipole_entry_J_struct> > bucketList(noOfBoxIndexes);
  for(int i = 0; i < noOfBoxIndexes; i++)
    bucketList[i].resize(counterList1[i]);

  std::vector<int> counterList2(noOfBoxIndexes);
  for(int i = 0; i < noOfBoxIndexes; i++)
    counterList2[i] = 0;
  for(int i = 0; i < n; i++)
    {
      int currBoxIndex = jobList[i].boxIndex;
      if(counterList2[currBoxIndex] >= counterList1[currBoxIndex])
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in sort_list_of_multipole_jobs: (counterList2[currBoxIndex] >= counterList1[currBoxIndex])");
	  return -1;
	}
      bucketList[currBoxIndex][counterList2[currBoxIndex]] = jobList[i];
      counterList2[currBoxIndex]++;
    }

  // OK, bucket-sort done. Now sort contents of each bucket.
  int count = 0;
  for(int i = 0; i < noOfBoxIndexes; i++)
    {
      int currCount = counterList2[i];
      if(sort_list_of_multipole_jobs_fixed_boxIndex(&bucketList[i][0], currCount) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in sort_list_of_multipole_jobs_fixed_boxIndex");
	  return -1;
	}
      // Copy result
      memcpy(&jobList[count], &bucketList[i][0], currCount * sizeof(job_list_multipole_entry_J_struct));
      count += currCount;
    }

  // check that list is sorted
  for(int i = 0; i < n-1; i++)
    {
      if(jobList[i].boxIndex > jobList[i+1].boxIndex)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: list not sorted!");
	  return -1;
	}
      if(jobList[i].boxIndex == jobList[i+1].boxIndex)
	{
	  if(jobList[i].branchIndex > jobList[i+1].branchIndex)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: list not sorted!");
	      return -1;
	    }
	  if(jobList[i].branchIndex == jobList[i+1].branchIndex)
	    {
	      if(jobList[i].multipoleBoxIndex > jobList[i+1].multipoleBoxIndex)
		{
		  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: list not sorted!");
		  return -1;
		}
	    }
	}
    }

  timeMeterInit.print(LOG_AREA_INTEGRALS, "sort_list_of_multipole_jobs complete");

  return 0;
}



/** executes given jobList using FMM.
    @param jobIndexLo the first jobindex for which this thread is responsible.
    @param jobIndexHi the last jobindex for which this thread is responsible is jobIndexHi-1.
    @param integralInfo info needed for evaluation of integrals of Gaussian functions.
    @param J_K_params includes various parameters for J and K matrix construction.
    @param jobList_J_multipole list of multipole-jobs.
    @param boxList list of boxes.
    @param maxnoOfMinimalDistrsPerBoxBranch needed to determine size of work buffer.
    @param result_J_list the list of matrix elements to be updated.
    @param largest_L_used largest L-value used (output).
*/
static int
execute_joblist_J_fmm_shared(int jobIndexLo, int jobIndexHi,
                             const IntegralInfo& integralInfo,
                             const JK::Params& J_K_params,
                             const job_list_multipole_entry_J_struct* jobList_J_multipole,
                             const box_struct *boxList,
                             int maxnoOfMinimalDistrsPerBoxBranch,
                             ergo_real* result_J_list,
                             int* largest_L_used)
{
  // Execute multipole job list for J
  int boxIndexSaved = -1;
  int branchIndexSaved = -1;

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "execute_joblist_J_fmm_shared: Allocating multipoleList_4, maxnoOfMinimalDistrsPerBoxBranch = %9i", 
	    maxnoOfMinimalDistrsPerBoxBranch);
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "execute_joblist_J_fmm_shared: jobIndexLo = %12d, jobIndexHi = %12d", 
	    jobIndexLo, jobIndexHi);

  std::vector<multipole_struct_small> multipoleList_4(maxnoOfMinimalDistrsPerBoxBranch);

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating multipoleList_4");

  MMInteractor interactor(integralInfo.GetMultipolePrep());
  
  int jobIndex = jobIndexLo;
  while(jobIndex<jobIndexHi)
    {
      // check how many of the following jobs that differ only in multipoleBranchIndex
      int jobIndex2 = jobIndex;

      int boxIndex = jobList_J_multipole[jobIndex].boxIndex;
      int branchIndex = jobList_J_multipole[jobIndex].branchIndex;
      int multipoleBoxIndex = jobList_J_multipole[jobIndex].multipoleBoxIndex;

      while(++jobIndex2 < jobIndexHi)
        {
          if(jobList_J_multipole[jobIndex2].boxIndex != boxIndex)
            break;
          if(jobList_J_multipole[jobIndex2].branchIndex != branchIndex)
            break;
          if(jobList_J_multipole[jobIndex2].multipoleBoxIndex != multipoleBoxIndex)
            break;
        }
      int nJobs = jobIndex2 - jobIndex;

      // check if we need to create new list of multipoles
      if(boxIndex != boxIndexSaved || branchIndex != branchIndexSaved) {
	if(create_list_of_multipoles_for_box(integralInfo, boxList[boxIndex].branchListForJ[branchIndex].org, &multipoleList_4[0]) != 0) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in create_list_of_multipoles_for_box().");
	  return -1;
	}
	// save these boxIndex and branchIndex values, so that we do not need to recompute multipoles until we reach the next box/branch
	boxIndexSaved = boxIndex;
	branchIndexSaved = branchIndex;
      } // END IF need to create new list of multipoles

      // OK, now we have nJobs, which is at least 1
      int first_multipoleBranchIndex = jobList_J_multipole[jobIndex].multipoleBranchIndex;
      multipole_struct_large multipoleSum = boxList[multipoleBoxIndex].branchListForJ[first_multipoleBranchIndex].org_mm.data.multipole;
      memset(multipoleSum.momentList, 0, MAX_NO_OF_MOMENTS_PER_MULTIPOLE*sizeof(ergo_real));
      for(int jobNo = 0; jobNo < nJobs; jobNo++)
	{
	  int multipoleBranchIndex = jobList_J_multipole[jobIndex+jobNo].multipoleBranchIndex;
	  const multipole_struct_large* branchMultipole = &boxList[multipoleBoxIndex].branchListForJ[multipoleBranchIndex].org_mm.data.multipole;
	  for(int mm = 0; mm < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; mm++)
	    multipoleSum.momentList[mm] += branchMultipole->momentList[mm];
	} // END FOR jobNo
      setup_multipole_maxAbsMomentList(&multipoleSum);

      if(do_multipole_interaction_between_2_boxes_branches(boxList[boxIndex].branchListForJ[branchIndex],
							   multipoleSum,
							   &multipoleList_4[0],
							   result_J_list,
							   NULL,
							   J_K_params.threshold_J * J_K_params.multipole_threshold_factor,
							   largest_L_used,
							   interactor,
							   integralInfo.GetMMLimitTable()
							   ) != 0) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in do_multipole_interaction_between_2_boxes_branches().");
	return -1;
      }

      jobIndex = jobIndex2;
      
    } // END WHILE (jobIndex < jobIndexHi)

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "execute_joblist_J_fmm_shared: done!");
  
  return 0;
}

static int
execute_joblist_J_fmm_serial(const IntegralInfo& integralInfo,
                             const JK::Params& J_K_params,
                             int noOfJobs_J_multipole,
                             const job_list_multipole_entry_J_struct* jobList_J_multipole,
                             const box_struct *boxList,
                             int maxnoOfMinimalDistrsPerBoxBranch,
                             ergo_real* result_J_list)
{
  Util::TimeMeter timeMeterJmul;
  int largest_L_used = 0;
  int rc = execute_joblist_J_fmm_shared(0, noOfJobs_J_multipole,
                                        integralInfo,
                                        J_K_params,
                                        jobList_J_multipole,
                                        boxList,
                                        maxnoOfMinimalDistrsPerBoxBranch,
                                        result_J_list,
                                        &largest_L_used);
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "multipole job list for J executed, largest L used: %2i", 
	    largest_L_used);
  timeMeterJmul.print(LOG_AREA_INTEGRALS, "Executing multipole job list for J");
  return rc;
}

struct JFMMWorkerData {
  const IntegralInfo* integralInfo;
  const JK::Params* J_K_params;
  const job_list_multipole_entry_J_struct* jobList_J_multipole;
  const box_struct *boxList;
  ergo_real* result_J_list;
  pthread_t threadID;
  int jobIndexLo, jobIndexHi;
  int maxnoOfMinimalDistrsPerBoxBranch;
  int result;
  int largest_L_used;
};

static void*
execute_J_fmm_worker(void *arg)
{
  JFMMWorkerData *data = static_cast<JFMMWorkerData*>(arg);
  
  data->result = 
    execute_joblist_J_fmm_shared(data->jobIndexLo, data->jobIndexHi,
                                 *data->integralInfo,
				 *data->J_K_params,
                                 data->jobList_J_multipole,
                                 data->boxList,
                                 data->maxnoOfMinimalDistrsPerBoxBranch,
                                 data->result_J_list,
                                 &data->largest_L_used);
  return NULL;
}



/** Compute the FMM part of the Coulomb matrix using threads. 0th
    thread reuses result_J_list, all the other threads need to have temporary
    memory allocated.
*/
static int
execute_joblist_J_fmm_thread(int noOfThreads, int noOfBasisFuncIndexPairs,
                             const IntegralInfo& integralInfo,
                             const JK::Params& J_K_params,
                             int noOfJobs_J_multipole,
                             const job_list_multipole_entry_J_struct* jobList_J_multipole,
                             const box_struct *boxList,
                             int maxnoOfMinimalDistrsPerBoxBranch,
                             ergo_real* result_J_list)
{
  std::vector<JFMMWorkerData> threadData(noOfThreads);
  int lastJob = 0;
  Util::TimeMeter timeMeterJmul;

  int th; // declare here because value used after loop after break.
  for(th = 0; th < noOfThreads; th++) {
    threadData[th].integralInfo = &integralInfo;
    threadData[th].J_K_params   = &J_K_params;
    threadData[th].jobList_J_multipole  = jobList_J_multipole;
    threadData[th].boxList = boxList;
    threadData[th].maxnoOfMinimalDistrsPerBoxBranch = maxnoOfMinimalDistrsPerBoxBranch;
    threadData[th].jobIndexLo = lastJob;
    /* Now we want to compute jobIndexHi as
       (((th+1)*noOfJobs_J_multipole)/noOfThreads) but we need to be
       careful to avoid integer overflow if noOfJobs_J_multipole is
       large. */
    size_t noOfJobs_J_multipole_as_size_t = noOfJobs_J_multipole;
    threadData[th].jobIndexHi = lastJob = ((th+1)*noOfJobs_J_multipole_as_size_t)/noOfThreads;
    threadData[th].largest_L_used = 0;
    if(th) {
      threadData[th].result_J_list = new ergo_real[noOfBasisFuncIndexPairs];
      if( threadData[th].result_J_list == NULL) {
        do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error allocating data for thread %i", th);
        do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "waiting for already created threads..");
        th++; /* Correct the wait loop upper index. */
        break;
      }
      memset(threadData[th].result_J_list, 0, noOfBasisFuncIndexPairs*sizeof(ergo_real));
    } else {
      threadData[th].result_J_list = result_J_list;
    }
    if(pthread_create(&threadData[th].threadID, NULL,
                      execute_J_fmm_worker, &threadData[th]) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in pthread_create for thread %i", th);
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "waiting for already created threads..");
      th++; /* Correct the wait loop upper index. */
      break;
    }
  }
  int myResult = 0, largest_L_used = 0;
  for(int i = 0; i < th; i++) {
    if(pthread_join(threadData[i].threadID, NULL) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in pthread_join for thread %i", i);
      myResult++;
    }
    myResult += threadData[i].result;
    if(i == 0) {
      threadData[i].result_J_list = NULL;
    } else {
      for(int idx=0; idx<noOfBasisFuncIndexPairs; idx++)
        result_J_list[idx] += threadData[i].result_J_list[idx];
      delete [] threadData[i].result_J_list;
    }
    if(threadData[i].largest_L_used > largest_L_used)
      largest_L_used = threadData[i].largest_L_used;
  }

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "multipole job list for J executed, largest L used: %2i",
	    largest_L_used);
  timeMeterJmul.print(LOG_AREA_INTEGRALS, "Executing multipole job list for J");
  return myResult;
}


/** Computes the Coulomb interaction.
    @param basisInfo
    @param integralInfo
    @param J_K_params the evaluation parameters, thresholds and all.
    @param basisFuncIndexPairList
    @param basisFuncIndexPairCount the length of basisFuncIndexPairList.
    @param D_list basisFuncIndexPairCount elements, with indices
    matching basisFuncIndexPairList.
    @param result_J_list preallocated list that will contain the results.
    @param noOfBasisFuncIndexPairs the length of result_J_list.
    happens to be always equal to basisFuncIndexPairCount
 */
int
compute_J_by_boxes_linear(const BasisInfoStruct & basisInfo,
			  const IntegralInfo & integralInfo,
			  const JK::Params& J_K_params,
			  const basis_func_index_pair_struct* basisFuncIndexPairList,
			  int basisFuncIndexPairCount,
			  const ergo_real* D_list,
			  ergo_real* result_J_list,
			  int noOfBasisFuncIndexPairs)
{
  Util::TimeMeter timeMeterTot;
  
  int n = basisInfo.noOfBasisFuncs;

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "entering compute_J_by_boxes_linear, no of basis funcs = %5i, threshold_J = %7.3g", 
	    n, (double)J_K_params.threshold_J);
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "use_fmm = %i, fmm_box_size = %6.2f", 
	    J_K_params.use_fmm, (double)J_K_params.fmm_box_size);

  output_current_memory_usage(LOG_AREA_INTEGRALS, "beginning of compute_J_by_boxes_linear");
  
  Util::TimeMeter timeMeterDistrList;

  // get largest limiting factor
  ergo_real maxLimitingFactor = 0;
  if(get_list_of_labeled_distrs_maxLimitingFactor_linear(basisInfo,
							 integralInfo,
							 J_K_params.threshold_J,
							 basisFuncIndexPairList,
							 basisFuncIndexPairCount,
							 &maxLimitingFactor) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_list_of_labeled_distrs_maxLimitingFactor_linear");
      return -1;
    }

  // Get number of distributions
  int distrCount = get_list_of_labeled_distrs_linear(basisInfo,
						     integralInfo,
						     J_K_params.threshold_J,
						     NULL,
						     0,
						     maxLimitingFactor,
						     basisFuncIndexPairList,
						     basisFuncIndexPairCount,
						     D_list);
  if(distrCount <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_J_by_boxes_linear: (distrCount <= 0)");
      return -1;
    }

  std::vector<DistributionSpecStructLabeled> distrList(distrCount);

  // create list of product primitives, with labels
  int distrCountTemp = get_list_of_labeled_distrs_linear(basisInfo,
							 integralInfo,
							 J_K_params.threshold_J,
							 &distrList[0],
							 distrCount,
							 maxLimitingFactor,
							 basisFuncIndexPairList,
							 basisFuncIndexPairCount,
							 D_list);
  if(distrCountTemp != distrCount)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_J_by_boxes_linear: (distrCountTemp != distrCount)");
      return -1;
    }

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating list of primitive distributions");

  ergo_real maxDensityMatrixElement = get_max_abs_vector_element(basisFuncIndexPairCount, D_list);

  // compute extent for all distrs
  Util::TimeMeter timeMeterComputeExtentForAllDistrs;
  compute_extent_for_list_of_distributions(distrCount, 
					   &distrList[0], 
					   J_K_params.threshold_J,
					   maxLimitingFactor,
					   maxDensityMatrixElement);
  timeMeterComputeExtentForAllDistrs.print(LOG_AREA_INTEGRALS, "Compute extent for all distrs");

  // get largest and smallest extent
  ergo_real extent_min, extent_max;
  get_largest_and_smallest_extent_for_list_of_distributions(distrCount, &distrList[0], &extent_min, &extent_max);
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "extent_min = %8.3f, extent_max = %8.3f", (double)extent_min, (double)extent_max);


  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "Creating list of distributions etc done, distrCount = %9i", distrCount);
  timeMeterDistrList.print(LOG_AREA_INTEGRALS, "Creating list of distributions etc");


  //
  // This is where we start to worry about the box system
  //

  Util::TimeMeter timeMeterBoxes;
  BoxSystem boxSystem;
  const ergo_real toplevelBoxSize = J_K_params.fmm_box_size;
  if(create_box_system_and_reorder_distrs(distrCount, 
					  &distrList[0],
					  toplevelBoxSize,
					  boxSystem) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system_and_reorder_distrs");
      return -1;
    }



  // Now we have the box system.
  // Create new list of boxes (more advanced boxes this time)

  std::vector<box_struct> boxList(boxSystem.totNoOfBoxes);

  for(int i = 0; i < boxSystem.totNoOfBoxes; i++)
    boxList[i].basicBox = boxSystem.boxList[i];


  int numberOfLevels = boxSystem.noOfLevels;
  int levelCounterList[numberOfLevels];
  int levelStartIndexList[numberOfLevels];
  for(int i = 0; i < numberOfLevels; i++)
    {
      levelCounterList[i] = boxSystem.levelList[i].noOfBoxes;
      levelStartIndexList[i] = boxSystem.levelList[i].startIndexInBoxList;
    }

  
  // OK, boxes created.

  timeMeterBoxes.print(LOG_AREA_INTEGRALS, "Creating boxes");


  int noOfBoxesTopLevel = levelCounterList[numberOfLevels-1];
  box_struct* boxListTopLevel = &boxList[levelStartIndexList[numberOfLevels-1]];
  

  // within each box, divide distrs into branches according to how far they penetrate outside the box.
  ergo_real branchSplitterList[MAX_NO_OF_BRANCHES];

  int noOfBranches = get_branch_splitter_info(branchSplitterList,
					      MAX_NO_OF_BRANCHES,
					      J_K_params,
					      toplevelBoxSize,
					      extent_max);
  if(noOfBranches <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_branch_splitter_info");
      return -1;
    }

  char s[888];
  s[0] = '\0';
  for(int i = 0; i < noOfBranches-1; i++)
    {
      char ss[888];
      sprintf(ss, " %5.2f", (double)branchSplitterList[i]);
      strcat(s, ss);
    }
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "noOfBranches = %i, splitters: %s", noOfBranches, s);


  if(create_branches(noOfBranches,
		     branchSplitterList,
		     distrCount,
		     &distrList[0],
		     noOfBoxesTopLevel,
		     boxListTopLevel) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_branches");
      return -1;
    }


  Util::TimeMeter timeMeterJorg;
  int groupCount = 0;
  int maxnoOfMinimalDistrsPerBoxBranch = 0;
  int maxNoOfMonomials = 0;
  for(int i = 0; i < noOfBoxesTopLevel; i++) {
    for(int branchIndex = 0; branchIndex < noOfBranches; branchIndex++) {
      DistributionSpecStructLabeled* distrListCurrBox = &distrList[boxListTopLevel[i].branchIndexListForJ[branchIndex]];
      int distrCountCurrBox = boxListTopLevel[i].branchCountListForJ[branchIndex];
      if(organize_distributions(integralInfo,
				distrListCurrBox, 
				distrCountCurrBox,
				&boxListTopLevel[i].branchListForJ[branchIndex].org,
				boxListTopLevel[i].basicBox.centerCoords,
				boxListTopLevel[i].basicBox.width) != 0) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in organize_distributions for box %i branch %i", i, branchIndex);
	return -1;
      }
      groupCount += boxListTopLevel[i].branchListForJ[branchIndex].org.groupList.size();
      int minimalDistrCount = boxListTopLevel[i].branchListForJ[branchIndex].org.minimalDistrList.size();
      if(minimalDistrCount > maxnoOfMinimalDistrsPerBoxBranch)
	maxnoOfMinimalDistrsPerBoxBranch = minimalDistrCount;
      int maxNoOfMonomialsCurr = boxListTopLevel[i].branchListForJ[branchIndex].org.data.maxNoOfMonomials;
      if(maxNoOfMonomialsCurr > maxNoOfMonomials)
	maxNoOfMonomials = maxNoOfMonomialsCurr;
    } // END FOR branchIndex
  } // END FOR i
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "J org done, groupCount = %8i", groupCount);
  timeMeterJorg.print(LOG_AREA_INTEGRALS, "J org");


  // Generate multipole for each group, and find center-of-charge for each branch
  Util::TimeMeter timeMeterGenerateGr;

  ergo_real totChargeWholeSystem = 0;
  for(int i = 0; i < noOfBoxesTopLevel; i++) {
    ergo_real averagePosList[3];
    for(int kk = 0; kk < 3; kk++)
      averagePosList[kk] = 0;
    int avgPosCounter = 0;
    for(int branchIndex = 0; branchIndex < noOfBranches; branchIndex++) {
      if(generate_multipoles_for_groups(integralInfo,
					boxListTopLevel[i].branchListForJ[branchIndex].org,
					boxListTopLevel[i].branchListForJ[branchIndex].org_mm,
					averagePosList,
					avgPosCounter // FIXME ELIAS: remove this if it is known to be same as groupCount
					) != 0) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in generate_multipoles_for_groups.");
	return -1;
      }
      totChargeWholeSystem += boxListTopLevel[i].branchListForJ[branchIndex].org_mm.data.chargeSum;
    } // END FOR branchIndex
    ergo_real multipolePointCurrBox[3];
    if(get_multipole_pt_for_box(boxListTopLevel[i].basicBox.centerCoords, boxListTopLevel[i].basicBox.width,
				averagePosList, avgPosCounter, multipolePointCurrBox) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in determine_multipole_pt_for_box.");
      return -1;
    }
    // Now we have determined multipolePointCurrBox
    // Copy it to each branch
    for(int branchIndex = 0; branchIndex < noOfBranches; branchIndex++) {
      for(int kk = 0; kk < 3; kk++)
	boxListTopLevel[i].branchListForJ[branchIndex].org_mm.data.multipolePoint[kk] = multipolePointCurrBox[kk];
    }
  } // END FOR i

  timeMeterGenerateGr.print(LOG_AREA_INTEGRALS, "Generate group multipoles");
  
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "totChargeWholeSystem = %22.11f", (double)totChargeWholeSystem);

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "Generating multipole for each branch at top level, MAX_MULTIPOLE_DEGREE = %2i", (int)MAX_MULTIPOLE_DEGREE);

  // Generate multipole for each branch at top level (smallest boxes)
  Util::TimeMeter timeMeterTranslate1;
  MMTranslator translator(integralInfo.GetMultipolePrep());
  for(int i = 0; i < noOfBoxesTopLevel; i++) {
    for(int branchIndex = 0; branchIndex < noOfBranches; branchIndex++) {
      if(translate_multipoles_for_box(boxListTopLevel[i].branchListForJ[branchIndex].org_mm,
				      boxListTopLevel[i].branchListForJ[branchIndex].org,
				      translator) != 0) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in translate_multipoles_for_box.");
	return -1;
      }
    } // END FOR branchIndex
  } // END FOR i
  timeMeterTranslate1.print(LOG_AREA_INTEGRALS, "Translate multipoles (step 1)");


  // OK, multipoles created for top level.
  // Now go through the other levels, joining multipoles from child boxes to a single multipole (one per branch) in parent box

  Util::TimeMeter timeMeterTranslate2;
  for(int levelNumber = numberOfLevels-2; levelNumber >= 0; levelNumber--) {
    int noOfBoxesCurrLevel = levelCounterList[levelNumber];
    box_struct* boxListCurrLevel = &boxList[levelStartIndexList[levelNumber]];
    for(int boxIndex = 0; boxIndex < noOfBoxesCurrLevel; boxIndex++) {
      box_struct* currBox = &boxListCurrLevel[boxIndex];
      int noOfChildren = currBox->basicBox.noOfChildBoxes;
      if(noOfChildren == 0) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "ERROR: (noOfChildren == 0)");
	return -1;
      }
      for(int branchIndex = 0; branchIndex < noOfBranches; branchIndex++) {
	distr_list_description_struct & currBoxBranch = currBox->branchListForJ[branchIndex];
	const distr_list_description_struct* childBoxBranches[8];
	for(int i = 0; i < 8; i++)
	  childBoxBranches[i] = NULL;
	for(int childIndex = 0; childIndex < noOfChildren; childIndex++) {
	  int childIndexInBoxList = currBox->basicBox.firstChildBoxIndex + childIndex;
	  box_struct* childBox = &boxList[childIndexInBoxList];
	  childBoxBranches[childIndex] = &childBox->branchListForJ[branchIndex];
	}
	assert(noOfChildren <= 8);
	if(combine_mm_info_for_child_boxes(currBoxBranch, childBoxBranches, noOfChildren, translator) != 0) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in combine_mm_info_for_child_boxes.");
	  return -1;
	}
      } // END FOR branchIndex
    } // END FOR boxIndex
  } // END FOR levelNumber
  timeMeterTranslate2.print(LOG_AREA_INTEGRALS, "Translate multipoles (step 2)");

  // Set J to zero
  memset(result_J_list, 0, basisFuncIndexPairCount*sizeof(ergo_real));

  // Create job lists for J

  Util::TimeMeter timeMeterJjoblist;
  int noOfJobs_J_standard_firstCount = 0;
  int noOfJobs_J_multipole_firstCount = 0;
  for(int branch_i = 0; branch_i < noOfBranches; branch_i++)
    {
      for(int branch_j = 0; branch_j < noOfBranches; branch_j++)
	{
	  int noOfNewJobs_standard = 0;
	  int noOfNewJobs_multipole = 0;
	  if(get_joblists_J_for_two_boxes_recursive(integralInfo,
						    J_K_params.threshold_J,
						    &boxList[0],
						    numberOfLevels,
						    0,
						    0,
						    0,
						    branch_i,
						    branch_j,
						    NULL,
						    HUGE_INTEGER_NUMBER,
						    &noOfNewJobs_standard,
						    NULL,
						    HUGE_INTEGER_NUMBER,
						    &noOfNewJobs_multipole
						    ) != 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_joblists_J_for_two_boxes_recursive");
	      return -1;
	    }
	  noOfJobs_J_standard_firstCount += noOfNewJobs_standard;
	  noOfJobs_J_multipole_firstCount += noOfNewJobs_multipole;
	}
    }
  std::vector<job_list_standard_entry_J_struct> jobList_J_standard(noOfJobs_J_standard_firstCount);
  std::vector<job_list_multipole_entry_J_struct> jobList_J_multipole(noOfJobs_J_multipole_firstCount);

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating jobLists for J");

  int noOfJobs_J_standard = 0;
  int noOfJobs_J_multipole = 0;
  for(int branch_i = 0; branch_i < noOfBranches; branch_i++)
    {
      for(int branch_j = 0; branch_j < noOfBranches; branch_j++)
	{
	  int noOfNewJobs_standard = 0;
	  int noOfNewJobs_multipole = 0;
	  if(get_joblists_J_for_two_boxes_recursive(integralInfo,
						    J_K_params.threshold_J,
						    &boxList[0],
						    numberOfLevels,
						    0,
						    0,
						    0,
						    branch_i,
						    branch_j,
						    &jobList_J_standard[noOfJobs_J_standard],
						    noOfJobs_J_standard_firstCount - noOfJobs_J_standard,
						    &noOfNewJobs_standard,
						    &jobList_J_multipole[noOfJobs_J_multipole],
						    noOfJobs_J_multipole_firstCount - noOfJobs_J_multipole,
						    &noOfNewJobs_multipole
						    ) != 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_joblists_J_for_two_boxes_recursive");
	      return -1;
	    }
	  noOfJobs_J_standard += noOfNewJobs_standard;
	  noOfJobs_J_multipole += noOfNewJobs_multipole;
	}
    }
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "job lists for J created OK, noOfJobs_J_standard = %8i, noOfJobs_J_multipole = %8i",
	    noOfJobs_J_standard, noOfJobs_J_multipole);
  timeMeterJjoblist.print(LOG_AREA_INTEGRALS, "Creating job lists for J");

  // Execute standard job list for J

  int noOfThreads = J_K_params.noOfThreads_J;
  
  if(noOfThreads <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_J_by_boxes_linear: (noOfThreads <= 0)");
      return -1;
    }
  if(noOfThreads == 1)
    {
      if(execute_joblist_J_std_serial(noOfJobs_J_standard,
				      &jobList_J_standard[0],
				      integralInfo,
				      maxNoOfMonomials,
				      result_J_list,
				      &boxList[0],
				      J_K_params.threshold_J) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in execute_joblist_J_std_serial");
	  return -1;
	}
    }
  else
    {
      if(execute_joblist_J_std_threaded(noOfThreads,
					noOfJobs_J_standard,
					&jobList_J_standard[0],
					integralInfo,
					maxNoOfMonomials,
					result_J_list,
					noOfBasisFuncIndexPairs,
					&boxList[0],
					J_K_params.threshold_J) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in execute_joblist_J_std_threaded");
	  return -1;
	}
    }


  // sort multipole job list by boxindex and branchindex.
  Util::TimeMeter timeMeterJmulSort;
  if(sort_list_of_multipole_jobs(&jobList_J_multipole[0], noOfJobs_J_multipole) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in sort_list_of_multipole_jobs");
      return -1;
    }
  timeMeterJmulSort.print(LOG_AREA_INTEGRALS, "sort_list_of_multipole_jobs");


  /* Execute multipole list */
  if(noOfThreads == 1)
    {
      if( execute_joblist_J_fmm_serial(integralInfo, J_K_params,
                                       noOfJobs_J_multipole, &jobList_J_multipole[0],
                                       &boxList[0], maxnoOfMinimalDistrsPerBoxBranch,
                                       result_J_list) != 0)
        {
          do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in execute_joblist_J_fmm_serial");
          return -1;
        }
    }
  else
    {
      if( execute_joblist_J_fmm_thread(noOfThreads, noOfBasisFuncIndexPairs,
                                       integralInfo, J_K_params,
                                       noOfJobs_J_multipole, &jobList_J_multipole[0],
                                       &boxList[0], maxnoOfMinimalDistrsPerBoxBranch,
                                       result_J_list) != 0)
        {
          do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in execute_joblist_J_fmm_thread");
          return -1;
        }
    }
  
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "compute_J_by_boxes_linear ending OK.");
  timeMeterTot.print(LOG_AREA_INTEGRALS, "compute_J_by_boxes_linear");
  
  return 0;
}




int
compute_J_by_boxes(const BasisInfoStruct & basisInfo,
		   const IntegralInfo & integralInfo,
		   const JK::Params& J_K_params,
		   ergo_real* J,
		   const ergo_real* dens)
{
  int n = basisInfo.noOfBasisFuncs;

  ergo_real maxDensityMatrixElement = get_max_abs_vector_element(n*n, dens);

  std::vector<basis_func_index_pair_struct> basisFuncIndexPairList;

  int noOfBasisFuncIndexPairs = get_basis_func_pair_list_2el(basisInfo,
							     integralInfo,
							     J_K_params.threshold_J,
							     maxDensityMatrixElement,
							     basisFuncIndexPairList);
  if(noOfBasisFuncIndexPairs <= 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_basis_func_pair_list");
    return -1;
  }

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "noOfBasisFuncIndexPairs = %i ==> storing %6.2f %% of a full matrix", 
	    noOfBasisFuncIndexPairs, (double)100*noOfBasisFuncIndexPairs/((double)n*n));

  std::vector<ergo_real> D_list(noOfBasisFuncIndexPairs);
  std::vector<ergo_real> J_list(noOfBasisFuncIndexPairs);

  // Setup D_list
  for(int i = 0; i < noOfBasisFuncIndexPairs; i++)
    {
      int a = basisFuncIndexPairList[i].index_1;
      int b = basisFuncIndexPairList[i].index_2;
      D_list[i] = dens[a*n+b];
    }
  
  if(compute_J_by_boxes_linear(basisInfo,
			       integralInfo,
			       J_K_params,
			       &basisFuncIndexPairList[0],
			       noOfBasisFuncIndexPairs,
			       &D_list[0],
			       &J_list[0],
			       noOfBasisFuncIndexPairs) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_J_by_boxes_linear");
      return -1;
    }

  // Now transfer result from J_list to J
  // First set all of J to zero (otherwise there may be something in the non-relevant part of the matrix)
  for(int i = 0; i < n*n; i++)
    J[i] = 0;
  for(int i = 0; i < noOfBasisFuncIndexPairs; i++)
    {
      int a = basisFuncIndexPairList[i].index_1;
      int b = basisFuncIndexPairList[i].index_2;
      J[a*n+b] = J_list[i];
      J[b*n+a] = J_list[i];
    }

  return 0;
}



/*
compute_J_by_boxes_nosymm does the same as compute_J_by_boxes,
but without assuming the density matrix to be symmetric.
*/
int
compute_J_by_boxes_nosymm(const BasisInfoStruct & basisInfo,
			  const IntegralInfo & integralInfo,
			  const JK::Params& J_K_params,
			  ergo_real* J,
			  const ergo_real* dens)
{
  int n = basisInfo.noOfBasisFuncs;

  // Create symmetrized density matrix P_ab = 0.5*(D_ab + D_ba)
  std::vector<ergo_real> P(n*n);
  for(int a = 0; a < n; a++)
    for(int b = 0; b < n; b++)
      P[a*n+b] = 0.5 * ( dens[a*n+b] + dens[b*n+a] );

  int rc = compute_J_by_boxes(basisInfo,
			      integralInfo,
			      J_K_params,
			      J,
			      &P[0]);
  return rc;
}






