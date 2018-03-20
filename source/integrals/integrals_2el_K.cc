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

/** @file integrals_2el_K.cc

    \brief Code for computing the Hartree-Fock exchange matrix K.

    @author: Elias Rudberg <em>responsible</em>.
*/

#include <string.h>
#include <stdio.h>

#include "integrals_2el_K.h"
#include "integrals_2el_utils.h"
#include "integrals_hermite.h"
#include "mm_limit_table.h"
#include "pi.h"
#include "pthread.h"
#include "utilities.h"
#include "matrix_algebra.h"
#include "integrals_2el_util_funcs.h"
#include "integrals_2el_K_kernel.h"
#include "integrals_2el_K_prep_groups.h"


static const int HUGE_INTEGER_NUMBER = 2000000000;


struct job_list_entry_K_struct {
  int boxIndex_1;
  int boxIndex_2;
  int useMultipole;
  ergo_real distance;
};


static int
create_joblist_exchange_for_two_boxes_recursive(const IntegralInfo & integralInfo,
						int maxNoOfMonomials,
						ergo_real threshold,
						const box_struct* boxList,
						int numberOfLevels,
						const csr_matrix_struct* dmatLimitMatrixCSRList,
						const int* basisFuncGroupCounterList,
						int currLevel,
						int boxIndex_1,
						int boxIndex_2,
						job_list_entry_K_struct* jobList_K,
						int maxNoOfJobs
						) {
  // Check if this pair of boxes can be skipped.
  int noOfRelevantBasisFuncGroups_1 = boxList[boxIndex_1].distrListForK.org.basisFuncGroupInfoListForK.size();
  int noOfRelevantBasisFuncGroups_2 = boxList[boxIndex_2].distrListForK.org.basisFuncGroupInfoListForK.size();
  const csr_matrix_struct* dmatLimitMatrixCSR = &dmatLimitMatrixCSRList[currLevel];

  // start by computing the minimum distance between the boxes.
  // We assume that both boxes have the same width.
  ergo_real dxList[3];
  for(int coordIndex = 0; coordIndex< 3; coordIndex++) {
    ergo_real x1 = boxList[boxIndex_1].basicBox.centerCoords[coordIndex];
    ergo_real x2 = boxList[boxIndex_2].basicBox.centerCoords[coordIndex];
    ergo_real dx = template_blas_fabs(x1 - x2);
    ergo_real width = boxList[boxIndex_1].basicBox.width;
    if(dx > width)
      dxList[coordIndex] = dx - width;
    else
      dxList[coordIndex] = 0;
  }
  ergo_real sumOfSquares = 0;
  for(int coordIndex = 0; coordIndex< 3; coordIndex++)
    sumOfSquares += dxList[coordIndex] * dxList[coordIndex];
  ergo_real distance = template_blas_sqrt(sumOfSquares);
  
  ergo_real maxDistanceOutsideBox_1 = boxList[boxIndex_1].distrListForK.org.data.maxDistanceOutsideBox;
  ergo_real maxDistanceOutsideBox_2 = boxList[boxIndex_2].distrListForK.org.data.maxDistanceOutsideBox;
  
  int useMultipole = 0;
  if(boxIndex_1 != boxIndex_2 && distance > maxDistanceOutsideBox_1 + maxDistanceOutsideBox_2)
    useMultipole = 1;
  
  ergo_real maxValue_CauschySchwartz = 0;
  ergo_real maxValue_multipole = 0;
  for(int i = 0; i < noOfRelevantBasisFuncGroups_1; i++)
    for(int j = 0; j < noOfRelevantBasisFuncGroups_2; j++) {
      ergo_real size_1 = boxList[boxIndex_1].distrListForK.org.basisFuncGroupInfoListForK[i].max_CS_factor;
      int index_1 = boxList[boxIndex_1].distrListForK.org.basisFuncGroupInfoListForK[i].basisFuncGroupIndex;
      ergo_real size_2 = boxList[boxIndex_2].distrListForK.org.basisFuncGroupInfoListForK[j].max_CS_factor;
      int index_2 = boxList[boxIndex_2].distrListForK.org.basisFuncGroupInfoListForK[j].basisFuncGroupIndex;
      if(index_1 < 0 || index_2 < 0) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_joblist_exchange_for_two_boxes_recursive: (index_1 < 0 || index_2 < 0)");
	return -1;
      }
      ergo_real maxDensElement = ergo_CSR_get_element(dmatLimitMatrixCSR, index_1, index_2);
      ergo_real currMax = size_1 * size_2 * maxDensElement;
      if(currMax > maxValue_CauschySchwartz)
	maxValue_CauschySchwartz = currMax;
      if(useMultipole == 1) {
	int degreeNeeded_1 = boxList[boxIndex_1].distrListForK.org.basisFuncGroupInfoListForK[i].maxMultipoleDegree;
	int degreeNeeded_2 = boxList[boxIndex_2].distrListForK.org.basisFuncGroupInfoListForK[j].maxMultipoleDegree;
	ergo_real maxAbsContributionFromMultipole = integralInfo.GetMMLimitTable().get_max_abs_mm_contrib(degreeNeeded_1,
													  boxList[boxIndex_1].distrListForK.org.basisFuncGroupInfoListForK[i].maxMomentVectorNormList,
													  degreeNeeded_2,
													  boxList[boxIndex_2].distrListForK.org.basisFuncGroupInfoListForK[j].maxMomentVectorNormList,
													  distance);
	ergo_real currMaxFromMultipole = maxAbsContributionFromMultipole * maxDensElement;
	if(currMaxFromMultipole > maxValue_multipole)
	  maxValue_multipole = currMaxFromMultipole;
      } // END IF useMultipole
    } // END FOR i j
  
  if(useMultipole == 1 && maxValue_multipole < threshold)
    return 0;
  if(maxValue_CauschySchwartz < threshold)
    return 0;
  if(currLevel == numberOfLevels-1) {
    // We are at the level of smallest boxes. Add job to job list.
    if(maxNoOfJobs <= 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_joblist_exchange_for_two_boxes_recursive: (maxNoOfJobs <= 0)");
      return -1;
    }
    if(jobList_K != NULL) {
      jobList_K[0].boxIndex_1 = boxIndex_1;
      jobList_K[0].boxIndex_2 = boxIndex_2;
      jobList_K[0].useMultipole = useMultipole;
      jobList_K[0].distance = distance;
    }
    return 1;
  }
  // Go to next level. Do interaction between all pairs of children of the two boxes.
  int noOfChildren_1 = boxList[boxIndex_1].basicBox.noOfChildBoxes;
  int noOfChildren_2 = boxList[boxIndex_2].basicBox.noOfChildBoxes;
  int jobCount = 0;
  for(int i = 0; i < noOfChildren_1; i++) {
    int start_j = 0;
    if(boxIndex_1 == boxIndex_2)
      start_j = i;
    for(int j = start_j; j < noOfChildren_2; j++) {
      int childIndex_1 = boxList[boxIndex_1].basicBox.firstChildBoxIndex + i;
      int childIndex_2 = boxList[boxIndex_2].basicBox.firstChildBoxIndex + j;
      job_list_entry_K_struct* jobList_K_mod = NULL;
      if(jobList_K != NULL)
	jobList_K_mod = &jobList_K[jobCount];
      int noOfJobs = create_joblist_exchange_for_two_boxes_recursive(integralInfo,
								     maxNoOfMonomials,
								     threshold,
								     boxList,
								     numberOfLevels,
								     dmatLimitMatrixCSRList,
								     basisFuncGroupCounterList,
								     currLevel + 1,
								     childIndex_1,
								     childIndex_2,
								     jobList_K_mod,
								     maxNoOfJobs - jobCount
								     );
      if(noOfJobs < 0) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_joblist_exchange_for_two_boxes_recursive for child boxes");
	return -1;
      }
      jobCount += noOfJobs;
    } // END FOR j
  } // END FOR i
  return jobCount;
}






struct K_joblist_thread_struct {
  pthread_t thread;
  int nBasisFuncs;
  const IntegralInfo* integralInfo;
  const JK::ExchWeights & CAM_params;
  csr_matrix_struct* K_CSR_shared;
  const csr_matrix_struct* densCSR;
  int maxNoOfMonomials;
  int basisFuncListCount_max;
  ergo_real threshold;
  const box_struct* boxList;
  const job_list_entry_K_struct* jobList_K;
  int noOfJobs_K_total;
  int thread_ID;
  int noOfThreads;
  int resultCode;
  int symmetryFlag;
  K_joblist_thread_struct(int nBasisFuncsIn,
			  const JK::ExchWeights & CAM_paramsIn) :
    nBasisFuncs(nBasisFuncsIn), CAM_params(CAM_paramsIn) { }
};


static void*
execute_joblist_K_thread_func(void* arg) {
  K_joblist_thread_struct* params = (K_joblist_thread_struct*)arg;
  try {
    const box_struct* boxList = params->boxList;
    int threadID = params->thread_ID;
    int noOfThreads = params->noOfThreads;

    JK_contribs_buffer_struct bufferStruct;
    allocate_buffers_needed_by_integral_code(*params->integralInfo, params->maxNoOfMonomials, params->basisFuncListCount_max, &bufferStruct);

    for(int jobIndex = 0; jobIndex < params->noOfJobs_K_total; jobIndex++) {
      if(jobIndex % noOfThreads != threadID)
	continue;
      int self = 0;
      int boxIndex_1 = params->jobList_K[jobIndex].boxIndex_1;
      int boxIndex_2 = params->jobList_K[jobIndex].boxIndex_2;
      if(boxIndex_1 == boxIndex_2)
	self = 1;
      if(get_K_contribs_from_2_interacting_boxes(*params->integralInfo,
						 params->CAM_params,
						 params->maxNoOfMonomials,
						 params->K_CSR_shared,
						 NULL,
						 params->densCSR,
						 params->symmetryFlag,
						 boxList[boxIndex_1].distrListForK.org,
						 boxList[boxIndex_2].distrListForK.org,
						 self,
						 params->threshold,
						 &bufferStruct,
						 params->jobList_K[jobIndex].useMultipole,
						 params->jobList_K[jobIndex].distance
						 ) != 0) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_K_contribs_from_2_interacting_boxes");
	params->resultCode = -1;
	return NULL;
      }
    } // END FOR jobIndex

    free_buffers_needed_by_integral_code(&bufferStruct);

    params->resultCode = 0;
  }
  catch(char const* e) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "=============================================================");
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "char* exception caught in execute_joblist_K_thread_func: '%s'", e);    
    do_output_time(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Time of exception: ");
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "=============================================================");
    params->resultCode = -1;
  }
  catch (std::exception & e) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "=============================================================");
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "std::exception caught in execute_joblist_K_thread_func: '%s'", e.what());
    do_output_time(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Time of exception: ");
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "=============================================================");
  }
  return NULL;
}



static int
execute_joblist_K_threaded(int noOfThreads,
			   csr_matrix_struct* densCSR,
			   int noOfBasisFuncs,
			   const IntegralInfo & integralInfo,
			   const JK::ExchWeights & CAM_params,
			   int maxNoOfMonomials,
			   int basisFuncListCount_max,
			   const box_struct* boxList,
			   const job_list_entry_K_struct* jobList_K,
			   int noOfJobs_K,
			   ergo_real threshold,
			   csr_matrix_struct* K_CSR,
			   int symmetryFlag
			   ) {
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "execute_joblist_K_threaded, noOfThreads = %2i, basisFuncListCount_max = %5i", noOfThreads, basisFuncListCount_max);

  K_joblist_thread_struct* threadParamsList[noOfThreads];
  
  // Set common parameters for all threads
  for(int i = 0; i < noOfThreads; i++) {
    threadParamsList[i] = new K_joblist_thread_struct(noOfBasisFuncs, CAM_params);
    threadParamsList[i]->densCSR = densCSR;
    threadParamsList[i]->integralInfo = &integralInfo;
    threadParamsList[i]->maxNoOfMonomials = maxNoOfMonomials;
    threadParamsList[i]->basisFuncListCount_max = basisFuncListCount_max;
    threadParamsList[i]->boxList = boxList;
    threadParamsList[i]->jobList_K = jobList_K;
    threadParamsList[i]->noOfJobs_K_total = noOfJobs_K;
    threadParamsList[i]->noOfThreads = noOfThreads;
    threadParamsList[i]->resultCode = -1; // initialize to error code
    threadParamsList[i]->threshold = threshold;
    threadParamsList[i]->symmetryFlag = symmetryFlag;
    threadParamsList[i]->K_CSR_shared = K_CSR;
  } // END FOR i
  
  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating memory for threads.");


  // Set ID number for all threads
  for(int i = 0; i < noOfThreads; i++)
    threadParamsList[i]->thread_ID = i;

  /* start threads */
  for(int i = 0; i < noOfThreads; i++) {
    if(pthread_create(&threadParamsList[i]->thread, 
		      NULL, 
		      execute_joblist_K_thread_func, 
		      threadParamsList[i]) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in pthread_create for thread %i", i);
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "waiting for already created threads..");
      for(int j = 0; j < i; j++) {
	if(pthread_join(threadParamsList[j]->thread, NULL) != 0)
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in pthread_join for thread %i", j);
      } /* END FOR j */
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "all threads finished, returning error code");
      return -1;
    }
  } /* END FOR i */

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "%i threads started OK.", noOfThreads);
  time_t workStartTime;
  time(&workStartTime);

  /* wait for threads to finish */
  for(int i = 0; i < noOfThreads; i++) {
    if(pthread_join(threadParamsList[i]->thread, NULL) != 0)
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in pthread_join for thread %i", i);
  } /* END FOR i */

  time_t workEndTime;
  time(&workEndTime);
  int secondsTaken = workEndTime - workStartTime;
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "all %i threads have finished, took %8i wall s.", noOfThreads, secondsTaken);
  
  /* now all threads have finished, check for errors */
  for(int i = 0; i < noOfThreads; i++) {
    if(threadParamsList[i]->resultCode != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in execute_joblist_K_thread_func"
		" for thread %i", i);
      return -1;
    }
  } /* END FOR i */
  
  for(int i = 0; i < noOfThreads; i++)
    delete threadParamsList[i];

  return 0;
}


static int
execute_joblist_K_serial(csr_matrix_struct* densCSR,
			 const IntegralInfo & integralInfo,
			 const JK::ExchWeights & CAM_params,
			 int maxNoOfMonomials,
			 int basisFuncListCount_max,
			 const box_struct* boxList,
			 const job_list_entry_K_struct* jobList_K,
			 int noOfJobs_K,
			 ergo_real threshold,
			 csr_matrix_struct* K_CSR,
			 int symmetryFlag) {
  JK_contribs_buffer_struct bufferStruct;
  allocate_buffers_needed_by_integral_code(integralInfo, maxNoOfMonomials, basisFuncListCount_max, &bufferStruct);
  for(int jobIndex = 0; jobIndex < noOfJobs_K; jobIndex++) {
    int self = 0;
    int boxIndex_1 = jobList_K[jobIndex].boxIndex_1;
    int boxIndex_2 = jobList_K[jobIndex].boxIndex_2;
    if(boxIndex_1 == boxIndex_2)
      self = 1;

    if(get_K_contribs_from_2_interacting_boxes(integralInfo,
					       CAM_params,
					       maxNoOfMonomials,
					       K_CSR,
					       NULL,
					       densCSR,
					       symmetryFlag,
					       boxList[boxIndex_1].distrListForK.org,
					       boxList[boxIndex_2].distrListForK.org,						
					       self,
					       threshold,
					       &bufferStruct,
					       jobList_K[jobIndex].useMultipole,
					       jobList_K[jobIndex].distance
					       ) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_K_contribs_from_2_interacting_boxes");
      return -1;
    }
  } // END FOR jobIndex

  free_buffers_needed_by_integral_code(&bufferStruct);

  return 0;
}



struct basisFuncGroupPairStruct {
  int i1;
  int i2;
};


static int
compare_basisFuncGroupPairs(const void* p1, const void* p2) {
  basisFuncGroupPairStruct* pair_1 = (basisFuncGroupPairStruct*)p1;
  basisFuncGroupPairStruct* pair_2 = (basisFuncGroupPairStruct*)p2;
  if(pair_1->i1 > pair_2->i1)
    return 1;
  if(pair_1->i1 < pair_2->i1)
    return -1;
  if(pair_1->i2 > pair_2->i2)
    return 1;
  if(pair_1->i2 < pair_2->i2)
    return -1;
  return 0;
}



static int
get_basisFuncGroupInfoList_maxsize(int distrCountTot,
				   const DistributionSpecStructLabeled* distrList,
				   int numberOfLevels,
				   const int* levelStartIndexList,
				   const int* levelCounterList,
				   const box_struct* boxList,
				   int** basisFuncGroupList) {
  int basisFuncGroupInfoList_count_max = 0;
  std::vector<basisFuncGroupPairStruct> pairList(2*distrCountTot);
  for(int levelNumber = 0; levelNumber < numberOfLevels; levelNumber++) {
    int pairCount = 0;
    for(int i = levelStartIndexList[levelNumber]; i < levelStartIndexList[levelNumber] + levelCounterList[levelNumber]; i++) {
      // go through all distrs of this box, and update basisFuncGroupInfoList accordingly.
      int distrStartIndex = boxList[i].basicBox.firstItemIndex;
      int distrCountCurrBox = boxList[i].basicBox.noOfItems;
      for(int j = distrStartIndex; j < distrStartIndex + distrCountCurrBox; j++) {
	const DistributionSpecStructLabeled & currDistr = distrList[j];
	int basisFuncGroup_1 = basisFuncGroupList[levelNumber][currDistr.basisFuncIndex_1];
	int basisFuncGroup_2 = basisFuncGroupList[levelNumber][currDistr.basisFuncIndex_2];
	pairList[pairCount].i1 = i;
	pairList[pairCount].i2 = basisFuncGroup_1;
	pairCount++;
	pairList[pairCount].i1 = i;
	pairList[pairCount].i2 = basisFuncGroup_2;
	pairCount++;
      } // END FOR j
    } // END FOR i
      // sort pairList
    qsort(&pairList[0], pairCount, sizeof(basisFuncGroupPairStruct), compare_basisFuncGroupPairs);
    int nn = 0;
    int i = 0;
    while(i < pairCount) {
      // now i should point to a new i1
      int i1 = pairList[i].i1;
      int j = i;
      while(j < pairCount && pairList[j].i1 == i1) {
	nn++;
	int i2 = pairList[j].i2;
	// now skip until another i2 is found.
	while(j < pairCount && pairList[j].i1 == i1 && pairList[j].i2 == i2)
	  j++;
      }
      i = j;
    }
    if(nn > basisFuncGroupInfoList_count_max)
      basisFuncGroupInfoList_count_max = nn;
  } // END FOR levelNumber
  return basisFuncGroupInfoList_count_max;
}




struct dmatElementStruct {
  int i1;
  int i2;
  ergo_real x;
};

static int 
compare_dmatElements(const void* p1, const void* p2) {
  dmatElementStruct* e1 = (dmatElementStruct*)p1;
  dmatElementStruct* e2 = (dmatElementStruct*)p2;
  if(e1->i1 > e2->i1)
    return 1;
  if(e1->i1 < e2->i1)
    return -1;
  if(e1->i2 > e2->i2)
    return 1;
  if(e1->i2 < e2->i2)
    return -1;
  return 0;
}




static int
create_reduced_vector(int nvalues,
		      const std::vector<dmatElementStruct> & dmatElementList,
		      std::vector<dmatElementStruct> & resultVector) {
  resultVector.resize(nvalues);
  int i = 0;
  int nvalues2 = 0;
  int curr_i1 = dmatElementList[0].i1;
  int curr_i2 = dmatElementList[0].i2;
  ergo_real curr_maxAbs = template_blas_fabs(dmatElementList[0].x);
  while(i < nvalues) {
    i++;
    int closeCurr = 0;
    if(i == nvalues)
      closeCurr = 1;
    else {
      // now we know it is safe to access element i
      int i1 = dmatElementList[i].i1;
      int i2 = dmatElementList[i].i2;
      if(i1 != curr_i1 || i2 != curr_i2)
	closeCurr = 1;
      else {
	// now we know this i is just a continuation of the current batch
	ergo_real absx = template_blas_fabs(dmatElementList[i].x);
	if(absx > curr_maxAbs)
	  curr_maxAbs = absx;
      }
    }
    if(closeCurr) {
      resultVector[nvalues2].i1 = curr_i1;
      resultVector[nvalues2].i2 = curr_i2;
      resultVector[nvalues2].x  = curr_maxAbs;
      nvalues2++;
      if(i < nvalues) {
	// Now we know it is safe to access element i. Start new batch.
	curr_i1 = dmatElementList[i].i1;
	curr_i2 = dmatElementList[i].i2;
	curr_maxAbs = template_blas_fabs(dmatElementList[i].x);
      }
    }
  }
  resultVector.resize(nvalues2);
  return nvalues2;
} /* End create_reduced_vector */



static int 
getDmatLimitMatrixCSRList(csr_matrix_struct* dmatLimitMatrixCSRList, 
			  int numberOfLevels,
			  const csr_matrix_struct* densCSR,
			  const int* const* basisFuncGroupList,
			  const int* basisFuncGroupCounterList)
{
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "getDmatLimitMatrixCSRList start.");
  memset(dmatLimitMatrixCSRList, 0, numberOfLevels*sizeof(csr_matrix_struct));
  for(int levelNumber = 0; levelNumber < numberOfLevels; levelNumber++) {
    // Populate dmatElementList with info for this level, one row at a time.
    int nn = basisFuncGroupCounterList[levelNumber];
    std::vector< std::vector<dmatElementStruct> > dmatElementListList(densCSR->n);
    for(int dmatrow = 0; dmatrow < densCSR->n; dmatrow++) {
      int nValuesCurrRow = ergo_CSR_get_nvalues_singlerow(densCSR, dmatrow);
      if(nValuesCurrRow == 0)
	continue;
      std::vector<int> colind(nValuesCurrRow);
      std::vector<ergo_real> values(nValuesCurrRow);
      if(ergo_CSR_get_values_singlerow(densCSR,
				       dmatrow,
				       &colind[0], 
				       &values[0],
				       nValuesCurrRow) != 0) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in ergo_CSR_get_values_singlerow.");
	return -1;
      }
      std::vector<dmatElementStruct> dmatElementList(nValuesCurrRow);
      for(int i = 0; i < nValuesCurrRow; i++) {
	int grrow = basisFuncGroupList[levelNumber][dmatrow];
	int grcol = basisFuncGroupList[levelNumber][colind[i]];
	if(grrow < grcol) {
	  dmatElementList[i].i1 = grrow;
	  dmatElementList[i].i2 = grcol;
	}
	else {
	  dmatElementList[i].i1 = grcol;
	  dmatElementList[i].i2 = grrow;
	}
	dmatElementList[i].x  = values[i];
      }
      // sort list to gather equal i1 i2 pairs together
      qsort(&dmatElementList[0], nValuesCurrRow, sizeof(dmatElementStruct), compare_dmatElements);
      // Create reduced vector.
      std::vector<dmatElementStruct> reducedVector;
      int nReduced = create_reduced_vector(nValuesCurrRow, dmatElementList, reducedVector);
      // Store result for this row in dmatElementListList.
      dmatElementListList[dmatrow].resize(nReduced);
      for(int i = 0; i < nReduced; i++)
	dmatElementListList[dmatrow][i] = reducedVector[i];
    }
    // OK, all rows done. Now create a single long list of all rows.
    int nTot = 0;
    for(int row = 0; row < densCSR->n; row++)
      nTot += dmatElementListList[row].size();
    std::vector<dmatElementStruct> dmatElementList(nTot);
    int count = 0;
    for(int row = 0; row < densCSR->n; row++)
      for(int i = 0; i < (int)dmatElementListList[row].size(); i++)
	dmatElementList[count++] = dmatElementListList[row][i];
    if(count != nTot) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in getDmatLimitMatrixCSRList: (count != nTot).");
      return -1;
    }
    // OK, single long list created.
    // sort list to gather equal i1 i2 pairs together
    qsort(&dmatElementList[0], nTot, sizeof(dmatElementStruct), compare_dmatElements);
    // Create reduced vector.
    std::vector<dmatElementStruct> reducedVector;
    int nReduced = create_reduced_vector(nTot, dmatElementList, reducedVector);
    // Create CSR matrix for this level.
    std::vector<int> rowind2(nReduced);
    std::vector<int> colind2(nReduced);
    std::vector<ergo_real> values2(nReduced);
    for(int i = 0; i < nReduced; i++) {
      rowind2[i] = reducedVector[i].i1;
      colind2[i] = reducedVector[i].i2;
      values2[i] = reducedVector[i].x;
    }
    // Create CSR
    csr_matrix_struct* currCSR = &dmatLimitMatrixCSRList[levelNumber];
    if(ergo_CSR_create(currCSR, 
		       1,
		       nn,
		       nReduced,
		       &rowind2[0],
		       &colind2[0]) != 0)
      {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in ergo_CSR_create for dmatLimitMatrixCSRList.");
	return -1;
      }
    for(int i = 0; i < nReduced; i++) {
      ergo_CSR_add_to_element(currCSR, 
			      rowind2[i],
			      colind2[i],
			      values2[i]);
    }
  }
  return 0;
}



/*
NOTE: This function adds its result to K.
This means that if only K is wanted, it must be set to zero before calling this function.
*/
int
compute_K_by_boxes(const BasisInfoStruct & basisInfo,
		   const IntegralInfo & integralInfo,
		   const JK::ExchWeights & CAM_params_in,
		   const JK::Params& J_K_params,
                   csr_matrix_struct* K_CSR,
		   csr_matrix_struct* densCSR,
		   int symmetryFlag)
{
  Util::TimeMeter timeMeterTot;

  int n = basisInfo.noOfBasisFuncs;

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "entering compute_K_by_boxes, no of basis funcs = %5i, threshold_K = %7.3g, exchange_box_size = %6.2f", 
	    n, (double)J_K_params.threshold_K, (double)J_K_params.exchange_box_size);

  output_current_memory_usage(LOG_AREA_INTEGRALS, "beginning of compute_K_by_boxes");

  const JK::ExchWeights CAM_params(CAM_params_in);


  Util::TimeMeter timeMeterDistrList;

  ergo_real maxDensityMatrixElement = ergo_CSR_get_max_abs_element(densCSR);


  // get largest limiting factor
  ergo_real maxLimitingFactor = 0;
  if(get_list_of_labeled_distrs_maxLimitingFactor(basisInfo,
						  integralInfo,
						  J_K_params.threshold_K,
						  &maxLimitingFactor,
						  maxDensityMatrixElement) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_list_of_labeled_distrs_maxLimitingFactor");
    return -1;
  }

  // Get number of distributions
  int distrCountTot = get_list_of_labeled_distrs(basisInfo,
						 integralInfo,
						 J_K_params.threshold_K,
						 NULL,
						 0,
						 maxLimitingFactor,
						 NULL,
						 maxDensityMatrixElement);
  if(distrCountTot == 0) {
    do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "compute_K_by_boxes: (distrCountTot == 0), skipping.");
    return 0;
  }
  if(distrCountTot <= 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_K_by_boxes: (distrCountTot <= 0)");
    return -1;
  }

  std::vector<DistributionSpecStructLabeled> distrList(distrCountTot);

  // create list of product primitives, with labels
  int distrCountTemp = get_list_of_labeled_distrs(basisInfo,
						  integralInfo,
						  J_K_params.threshold_K,
						  &distrList[0],
						  distrCountTot,
						  maxLimitingFactor,
						  NULL,
						  maxDensityMatrixElement);
  if(distrCountTemp != distrCountTot) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_K_by_boxes:(distrCountTemp != distrCountTot)");
    return -1;
  }

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating list of primitive distributions");


  // compute extent for all distrs
  Util::TimeMeter timeMeterComputeExtentForAllDistrs;
  compute_extent_for_list_of_distributions(distrCountTot, 
					   &distrList[0], 
					   J_K_params.threshold_K,
					   maxLimitingFactor,
					   maxDensityMatrixElement);
  timeMeterComputeExtentForAllDistrs.print(LOG_AREA_INTEGRALS, "Compute extent for all distrs");
  
  
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "Creating list of distributions done, distrCountTot = %9i", distrCountTot);
  timeMeterDistrList.print(LOG_AREA_INTEGRALS, "Creating list of distributions");




  //
  // This is where we start to worry about the box system
  //

  Util::TimeMeter timeMeterBoxes;

  BoxSystem boxSystem;
  if(create_box_system_and_reorder_distrs(distrCountTot,
					  &distrList[0],
					  J_K_params.exchange_box_size,
					  boxSystem) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system_and_reorder_distrs");
    return -1;
  }

  // Create new list of boxes (more advanced boxes this time)
  std::vector<box_struct> boxList(boxSystem.totNoOfBoxes);
  // FIXME: TODO: need to clear contents of boxList here?  
  for(int i = 0; i < boxSystem.totNoOfBoxes; i++)
    boxList[i].basicBox = boxSystem.boxList[i];
  
  int numberOfLevels = boxSystem.noOfLevels;
  int levelCounterList[numberOfLevels];
  int levelStartIndexList[numberOfLevels];
  for(int i = 0; i < numberOfLevels; i++) {
    levelCounterList[i] = boxSystem.levelList[i].noOfBoxes;
    levelStartIndexList[i] = boxSystem.levelList[i].startIndexInBoxList;
  }


  // Set up basisFuncGroups for all levels
  // Create another box system, this time with basis functions as items

  std::vector<box_item_struct> itemListBasisFuncs(n);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < 3; j++)
      itemListBasisFuncs[i].centerCoords[j] = basisInfo.basisFuncList[i].centerCoords[j];
    itemListBasisFuncs[i].originalIndex = i;
  } // END FOR i

  
  const ergo_real maxToplevelBoxSizeBasisFuncs = J_K_params.exchange_box_size;
  
  BoxSystem boxSystemBasisFuncs;

  if(boxSystemBasisFuncs.create_box_system(&itemListBasisFuncs[0],
					   n,
					   maxToplevelBoxSizeBasisFuncs * 0.5) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system");
    return -1;
  }

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after creating second box system");

  if(boxSystemBasisFuncs.noOfLevels < boxSystem.noOfLevels) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: (boxSystemBasisFuncs.noOfLevels < boxSystem.noOfLevels)");
    return -1;
  }
  int noOfLevelsBasisFuncs = boxSystemBasisFuncs.noOfLevels;
  int noOfLevelsDiff = noOfLevelsBasisFuncs - numberOfLevels;
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "noOfLevelsBasisFuncs = %i, noOfLevelsDiff = %i", noOfLevelsBasisFuncs, noOfLevelsDiff);

  // Now we set up basisFuncGroupList which for each basis function
  // contains information about which group that basis function
  // belongs to, at each level. So that later, if we have a basis
  // function index and a level, we can check which group that basis
  // function belongs to.
  int* basisFuncGroupList[numberOfLevels];
  int basisFuncGroupCounterList[numberOfLevels]; // Number of groups at each level
  for(int i = 0; i < numberOfLevels; i++)
    basisFuncGroupList[i] = new int[n];

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating basisFuncGroupList");

  int maxNoOfBasisFuncGroupsPerLevel = 0;
  for(int levelNumber = 0; levelNumber < numberOfLevels; levelNumber++) {
    // Set up basisFuncGroup list for this level.
    int noOfBoxesCurrLevel = boxSystemBasisFuncs.levelList[levelNumber+noOfLevelsDiff].noOfBoxes;
    int startIndex         = boxSystemBasisFuncs.levelList[levelNumber+noOfLevelsDiff].startIndexInBoxList;
    for(int i = startIndex; i < startIndex + noOfBoxesCurrLevel; i++) {
      // assign basis funcs of this box to basisFuncGroup i
      int firstItemIndex = boxSystemBasisFuncs.boxList[i].firstItemIndex;
      for(int j = firstItemIndex; j < firstItemIndex + boxSystemBasisFuncs.boxList[i].noOfItems; j++) {
	int basisFuncIndex = itemListBasisFuncs[j].originalIndex;
	basisFuncGroupList[levelNumber][basisFuncIndex] = i - startIndex;
      } // END FOR j
    } // END FOR i
    basisFuncGroupCounterList[levelNumber] = noOfBoxesCurrLevel;
    if(noOfBoxesCurrLevel > maxNoOfBasisFuncGroupsPerLevel)
      maxNoOfBasisFuncGroupsPerLevel = noOfBoxesCurrLevel;
  } // END FOR levelNumber
  
  // OK, basisFuncGroups done.


  // OK, boxes created.


  timeMeterBoxes.print(LOG_AREA_INTEGRALS, "Creating boxes etc");



  Util::TimeMeter timeMeterGetMultipoleNormVectors;

  // Create list of multipole norm vectors, for later use.
  std::vector<ergo_real> multipoleNormVectorList((MAX_MULTIPOLE_DEGREE_BASIC+1)*distrCountTot);
  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating multipoleNormVectorList");

  std::vector<int> multipoleDegreeList(distrCountTot);
  for(int j = 0; j < distrCountTot; j++) {
    DistributionSpecStructLabeled* currDistr = &distrList[j];
    ergo_real* multipoleNormVectorList_curr = &multipoleNormVectorList[j*(MAX_MULTIPOLE_DEGREE_BASIC+1)];
    multipole_struct_small multipole;
    if(compute_multipole_moments(integralInfo, &currDistr->distr, &multipole) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_multipole_moments");
      return -1;
    }
    multipoleDegreeList[j] = multipole.degree;
    for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
      multipoleNormVectorList_curr[l] = 0;
    for(int l = 0; l <= multipole.degree; l++) {
      int startIndex = l*l;
      int endIndex = (l+1)*(l+1);
      ergo_real sum = 0;
      for(int A = startIndex; A < endIndex; A++)
	sum += multipole.momentList[A]*multipole.momentList[A];
      ergo_real subNorm = template_blas_sqrt(sum);
      multipoleNormVectorList_curr[l] = subNorm;
    }
  } // END FOR j
  
  timeMeterGetMultipoleNormVectors.print(LOG_AREA_INTEGRALS, "getting multipoleNormVectorList");


  int noOfBoxesTopLevel = levelCounterList[numberOfLevels-1];
  box_struct* boxListTopLevel = &boxList[levelStartIndexList[numberOfLevels-1]];

  int basisFuncListCount_max = 0;
  int maxNoOfMonomials = 0;

  // Now call organize_distributions for each top-level box
  Util::TimeMeter timeMeterKorg;
  for(int i = 0; i < noOfBoxesTopLevel; i++) {
    DistributionSpecStructLabeled* distrListCurrBox = &distrList[boxListTopLevel[i].basicBox.firstItemIndex];
    int distrCountCurrBox = boxListTopLevel[i].basicBox.noOfItems;
    if(organize_distributions(integralInfo,
			      distrListCurrBox,
			      distrCountCurrBox,
			      &boxListTopLevel[i].distrListForK.org,
			      boxListTopLevel[i].basicBox.centerCoords,
			      boxListTopLevel[i].basicBox.width) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in organize_distributions for box %i", i);
      return -1;
    }
    int basisFuncListCount_curr = boxListTopLevel[i].distrListForK.org.basisFuncList.size();
    if(basisFuncListCount_curr > basisFuncListCount_max)
      basisFuncListCount_max = basisFuncListCount_curr;
    int maxNoOfMonomialsCurr = boxListTopLevel[i].distrListForK.org.data.maxNoOfMonomials;
    if(maxNoOfMonomialsCurr > maxNoOfMonomials)
      maxNoOfMonomials = maxNoOfMonomialsCurr;
  } // END FOR i
  timeMeterKorg.print(LOG_AREA_INTEGRALS, "K organize_distributions for all boxes");


  // Now go through the other levels, getting info for parent boxes

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
      // We want to get maxDistanceOutsideBox for parent box (use largest value found among the children).
      ergo_real maxDistanceOutsideBox = 0;
      for(int childIndex = 0; childIndex < noOfChildren; childIndex++) {
	int childIndexInBoxList = currBox->basicBox.firstChildBoxIndex + childIndex;
	box_struct* childBox = &boxList[childIndexInBoxList];
	if(childBox->distrListForK.org.data.maxDistanceOutsideBox > maxDistanceOutsideBox)
	  maxDistanceOutsideBox = childBox->distrListForK.org.data.maxDistanceOutsideBox;
      } // END FOR childIndex
      currBox->distrListForK.org.data.maxDistanceOutsideBox = maxDistanceOutsideBox;
    } // END FOR boxIndex
  } // END FOR levelNumber



  // For each box at each level, store information about the largest size distr associated with each basisFuncGroup.
  Util::TimeMeter timeMeterGetLimitsAllLevels;

  // Predict size of basisFuncGroupInfoList
  int basisFuncGroupInfoList_maxcount_predicted = get_basisFuncGroupInfoList_maxsize(distrCountTot,
										     &distrList[0],
										     numberOfLevels,
										     levelStartIndexList,
										     levelCounterList,
										     &boxList[0],
										     basisFuncGroupList);
  if(basisFuncGroupInfoList_maxcount_predicted <= 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_basisFuncGroupInfoList_maxsize");
    return -1;
  }

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "starting loop to setup basisFuncGroupInfoList, with added index checks.");
  
  for(int levelNumber = 0; levelNumber < numberOfLevels; levelNumber++) {
    // Set up basisFuncGroupInfoList for each box at level.
    
    for(int boxIndex = levelStartIndexList[levelNumber]; boxIndex < levelStartIndexList[levelNumber] + levelCounterList[levelNumber]; boxIndex++) {
      if(boxIndex < 0 || boxIndex >= boxSystem.totNoOfBoxes) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error doing basisFuncGroupInfoList: (boxIndex < 0 || boxIndex >= boxSystem.totNoOfBoxes)");
	return -1;
      }
      box_struct & currBox = boxList[boxIndex];

      const ergo_real* multipoleNormVectorList_currBox = &multipoleNormVectorList[currBox.basicBox.firstItemIndex*(MAX_MULTIPOLE_DEGREE_BASIC+1)];
      const int* multipoleDegreeList_currBox = &multipoleDegreeList[currBox.basicBox.firstItemIndex];

      int maxCount = basisFuncGroupInfoList_maxcount_predicted;
      int distrCountCurrBox = currBox.basicBox.noOfItems;
      const DistributionSpecStructLabeled* distrListCurrBox = &distrList[currBox.basicBox.firstItemIndex];

      std::vector<int> basisFuncGroupList1(distrCountCurrBox);
      std::vector<int> basisFuncGroupList2(distrCountCurrBox);
      std::vector<ergo_real> limitingFactorList(distrCountCurrBox);
      for(int jjj = 0; jjj < distrCountCurrBox; jjj++) {
	const DistributionSpecStructLabeled & currDistr = distrListCurrBox[jjj];
	int basisFuncGroup_1 = basisFuncGroupList[levelNumber][currDistr.basisFuncIndex_1];
	int basisFuncGroup_2 = basisFuncGroupList[levelNumber][currDistr.basisFuncIndex_2];
	basisFuncGroupList1[jjj] = basisFuncGroup_1;
	basisFuncGroupList2[jjj] = basisFuncGroup_2;
	limitingFactorList[jjj] = distrListCurrBox[jjj].limitingFactor;
      }

      if(prep_info_for_K(maxCount,
			 currBox.distrListForK.org,
			 distrCountCurrBox,
			 multipoleNormVectorList_currBox,
			 multipoleDegreeList_currBox,
			 &limitingFactorList[0],
			 &basisFuncGroupList1[0],
			 &basisFuncGroupList2[0]) != 0) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error: prep_info_for_K() failed.");
	return -1;
      }
      
    } // END FOR boxIndex
  } // END FOR levelNumber

  // OK, basisFuncGroup info done for all boxes.

  timeMeterGetLimitsAllLevels.print(LOG_AREA_INTEGRALS, "GetLimitsAllLevels");
  


  Util::TimeMeter timeMeterGetDensityMatrixLimitMatrixList;

  // Prepare densityMatrixLimit matrix for each level.
  // For a given level, the "densityMatrixLimit matrix" contains the
  // largest absolute dmat element for pairs of basis functions from
  // two basis function groups at that level.

  csr_matrix_struct dmatLimitMatrixCSRList[numberOfLevels];
  output_current_memory_usage(LOG_AREA_INTEGRALS, "before calling getDmatLimitMatrixCSRList");
  if(getDmatLimitMatrixCSRList(dmatLimitMatrixCSRList, 
			       numberOfLevels, 
			       densCSR, 
			       basisFuncGroupList, 
			       basisFuncGroupCounterList) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in getDmatLimitMatrixCSRList.");
    return -1;
  }

  // OK, densityMatrixLimitMatrixList done.

  timeMeterGetDensityMatrixLimitMatrixList.print(LOG_AREA_INTEGRALS, "getting densityMatrixLimitMatrixList");
  output_current_memory_usage(LOG_AREA_INTEGRALS, "after doing densityMatrixLimitMatrixList");


  // Crete job-list for K
  Util::TimeMeter timeMeterKjoblist;

  // compute number of jobs before allocating list.
  int noOfJobs_K_firstCount = create_joblist_exchange_for_two_boxes_recursive(integralInfo,
									      maxNoOfMonomials,
									      J_K_params.threshold_K,
									      &boxList[0],
									      numberOfLevels,
									      dmatLimitMatrixCSRList,
									      basisFuncGroupCounterList,
									      0,
									      0, 
									      0,
									      NULL,
									      HUGE_INTEGER_NUMBER
									      );
  if(noOfJobs_K_firstCount < 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_joblist_exchange_for_two_boxes_recursive");
    return -1;
  }

  std::vector<job_list_entry_K_struct> jobList_K(noOfJobs_K_firstCount);

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating jobList_K");

  int noOfJobs_K = create_joblist_exchange_for_two_boxes_recursive(integralInfo,
								   maxNoOfMonomials,
								   J_K_params.threshold_K,
								   &boxList[0],
								   numberOfLevels,
								   dmatLimitMatrixCSRList,
								   basisFuncGroupCounterList,
								   0,
								   0, 
								   0,
								   &jobList_K[0],
								   noOfJobs_K_firstCount
								   );
  if(noOfJobs_K < 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_joblist_exchange_for_two_boxes_recursive");
    return -1;
  }
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "job list for K created, %8i jobs", noOfJobs_K);
  timeMeterKjoblist.print(LOG_AREA_INTEGRALS, "creating job list for K");


  // Execute job-list for K

  Util::TimeMeter timeMeterK;

  int noOfThreads = J_K_params.noOfThreads_K;
  if(noOfThreads <= 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: (noOfThreads <= 0)");
    return -1;
  }
  if(noOfThreads == 1) {
    // no threading requested
    if(execute_joblist_K_serial(densCSR,
				integralInfo,
				CAM_params,
				maxNoOfMonomials,
				basisFuncListCount_max,
				&boxList[0],
				&jobList_K[0],
				noOfJobs_K,
				J_K_params.threshold_K,
				K_CSR,
				symmetryFlag) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in execute_joblist_K_serial");
      return -1;
    }
  }
  else
    {
      if(execute_joblist_K_threaded(noOfThreads,
				    densCSR,
				    basisInfo.noOfBasisFuncs,
				    integralInfo,
				    CAM_params,
				    maxNoOfMonomials,
				    basisFuncListCount_max,
				    &boxList[0],
				    &jobList_K[0],
				    noOfJobs_K,
				    J_K_params.threshold_K,
				    K_CSR,
				    symmetryFlag) != 0) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in execute_joblist_K_threaded");
	return -1;
      }
    }

  timeMeterK.print(LOG_AREA_INTEGRALS, "Executing job list for K");
  
  for(int i = 0; i < numberOfLevels; i++)
    delete [] basisFuncGroupList[i];

  for(int levelNumber = 0; levelNumber < numberOfLevels; levelNumber++)
    ergo_CSR_destroy(&dmatLimitMatrixCSRList[levelNumber]);

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after freeing stuff at end of compute_K_by_boxes");

  timeMeterTot.print(LOG_AREA_INTEGRALS, "compute_K_by_boxes");

  return 0;
}


int
compute_K_by_boxes_dense(const BasisInfoStruct & basisInfo,
			 const IntegralInfo & integralInfo,
			 const JK::ExchWeights & CAM_params_in,
			 const JK::Params& J_K_params,
			 ergo_real* K_dense,
			 const ergo_real* D_dense,
			 int symmetryFlag) {
  if(symmetryFlag != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in compute_K_by_boxes_dense: (symmetryFlag != 0) case not implemented.");
    return -1;
  }
  csr_matrix_struct K_CSR;
  csr_matrix_struct D_CSR;
  long n = basisInfo.noOfBasisFuncs;
  long nnz = n*n;
  std::vector<int> rowind(nnz);
  std::vector<int> colind(nnz);
  long count = 0;
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) {
      rowind[count] = i;
      colind[count] = j;
      count++;
    }
  assert(count == nnz);
  // Create zero K_CSR matrix where the result will be stored
  if(ergo_CSR_create(&K_CSR, 0, n, nnz, &rowind[0], &colind[0]) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in compute_K_by_boxes_dense, in ergo_CSR_create.");
    return -1;
  }
  // Create D_CSR and set its values to the values from D_dense
  if(ergo_CSR_create(&D_CSR, 0, n, nnz, &rowind[0], &colind[0]) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in compute_K_by_boxes_dense, in ergo_CSR_create.");
    return -1;
  }
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      ergo_CSR_add_to_element(&D_CSR, i, j, D_dense[i*n+j]);
  if(compute_K_by_boxes(basisInfo,
			integralInfo,
			CAM_params_in,
			J_K_params,
			&K_CSR,
			&D_CSR,
			symmetryFlag) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in compute_K_by_boxes_dense, in compute_K_by_boxes.");
    return -1;
  }
  // Transfer results from K_CSR to K_dense
  std::vector<int> colind_single_row(n);
  for(int i = 0; i < n; i++)
    colind_single_row[i] = i;
  std::vector<ergo_real> values_single_row(n);
  for(int i = 0; i < n; i++) {
    if(ergo_CSR_get_values_singlerow(&K_CSR, i, &colind_single_row[0], &values_single_row[0], n) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in compute_K_by_boxes_dense, in ergo_CSR_get_values_singlerow.");
      return -1;
    }
    for(int k = 0; k < n; k++)
      K_dense[i*n+k] = values_single_row[k];
  }
  return 0;
}


