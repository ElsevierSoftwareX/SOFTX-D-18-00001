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

/** @file integrals_2el_K_prep.cc

    \brief Code for preparing for computation of the Hartree-Fock
    exchange matrix K.

    @author: Elias Rudberg <em>responsible</em>.
*/

#include <cstring>
#include <vector>

#include <pthread.h>

#include "integrals_2el_K_prep.h"
#include "utilities.h"
#include "integrals_2el_utils.h"


typedef struct
{
  int i;
  ergo_real max_CS_factor;
} neighbor_basisfunc_struct;



static int
find_int_in_sorted_list(const int* list, int listLength, int i)
{
  if(listLength == 0)
    return 0;
  /* ELIAS NOTE 2014-03-27: chaning datatype from "int" to "unsigned
     int" for the lo, hi, mid variables here gave a considerable
     performance inprovement. Apparently the computation of mid is
     more efficient for unsigned int.  */
  unsigned int lo = 0;
  unsigned int hi = listLength-1;
  while(lo < hi-1)
    {
      unsigned int mid = (lo + hi) / 2;
      if(list[mid] > i)
	hi = mid;
      else
	lo = mid;
    } // END WHILE
  if(list[lo] == i)
    return 1;
  if(list[hi] == i)
    return 1;
  return 0;
}


static int
find_int_in_list(const int* list, int listLength, int i)
{
  for(int k = 0; k < listLength; k++)
    {
      if(list[k] == i)
	return 1;
    }
  return 0;
}


static void
do_sort_int_list(int* list, int n)
{
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n-i-1; j++)
      {
	if(list[j+1] < list[j])
	  {
	    int temp = list[j];
	    list[j] = list[j+1];
	    list[j+1] = temp;
	  }
      } // END FOR i j
}


static int
find_doubles_in_sorted_list(const int* list, int n)
{
  for(int i = 0; i < n-1; i++)
    {
      if(list[i] == list[i+1])
	return 1;
    }
  return 0;
}



static int
merge_lists(const int* list_1, int list_1_len, const int* list_2, int list_2_len,
            int* resultList)
{
  int i1 = 0;
  int i2 = 0;
  int p = 0;
  while(i1 < list_1_len && i2 < list_2_len)
    {
      if(list_1[i1] < list_2[i2])
	{
	  resultList[p] = list_1[i1];
	  i1++;
	}
      else
	{
	  resultList[p] = list_2[i2];
	  i2++;
	}
      p++;
    }
  while(i1 < list_1_len)
    {
      resultList[p] = list_1[i1];
      i1++;
      p++;
    }
  while(i2 < list_2_len)
    {
      resultList[p] = list_2[i2];
      i2++;
      p++;
    }
  return 0;
}




static int
identify_needed_elements_part(ergo_real threshold, 
			      const csr_matrix_struct* dens_CSR, 
			      const int noOfNeighborsList[],
			      const neighbor_basisfunc_struct* neighborList, 
			      const int maxNoOfNeighbors,
			      int ** longList, 
			      int *longListCounterList,
			      int myIndex,
			      int noOfParts)
{
  const int n = dens_CSR->n;
  int shortListCounterList[n];
  memset(shortListCounterList, 0, n*sizeof(int));
  /* ELIAS NOTE 2014-03-27: the "short list length" value chosen here
     has considerable impact on performance. Increasing it from 20 to
     100 made the identify_needed_elements_part call 10% faster. This
     effect could be machine-dependent, of course.  */
  const int sll = 100; // sll = "short list length"
  std::vector<int> shortList(n*sll); // shortList is n lists of sll elements each.

  // now go through all nonzero dmat elements.
  for(int r = 0; r < n; r++) {
    if(r % noOfParts != myIndex)
      continue;
    // Now take care of dmat elements on row r (and some other related elements in symmetric case)
    int nValuesCurrRow = ergo_CSR_get_nvalues_singlerow(dens_CSR, r);
    std::vector<int> rowind(nValuesCurrRow);
    std::vector<int> colind(nValuesCurrRow);
    std::vector<ergo_real> values(nValuesCurrRow);
    if(ergo_CSR_get_values_singlerow(dens_CSR, r, &colind[0], &values[0], nValuesCurrRow) != 0)
      return -1;
    for(int k = 0; k < nValuesCurrRow; k++)
      rowind[k] = r;
    int nValuesToHandle = nValuesCurrRow;
    if(dens_CSR->symmetryFlag) {
      // Symmetric matrix case. In this case we also want to handle elements on the other side of the matrix diagonal, which were not returned by ergo_CSR_get_values_singlerow().
      int nOffDiagElements = 0;
      for(int k = 0; k < nValuesCurrRow; k++) {
	if(colind[k] != r)
	  nOffDiagElements++;
      }
      nValuesToHandle = nValuesCurrRow + nOffDiagElements;
      rowind.resize(nValuesToHandle);
      colind.resize(nValuesToHandle);
      values.resize(nValuesToHandle);
      for(int k = 0; k < nValuesCurrRow; k++)
	rowind[k] = r;
      int currIdx = nValuesCurrRow;
      for(int k = 0; k < nValuesCurrRow; k++) {
	if(colind[k] != r) {
	  rowind[currIdx] = colind[k];
	  colind[currIdx] = r;
	  values[currIdx] = values[k];
	  currIdx++;
	}
      }
    } // end if dens_CSR->symmetryFlag
    // OK, now we have the nValuesToHandle elements to handle in the rowind colind values vectors.
    for(int idx = 0; idx < nValuesToHandle; idx++) {
      int i = rowind[idx];
      int j = colind[idx];
      ergo_real currDmatElement = values[idx];
      ergo_real absDmatElement = template_blas_fabs(currDmatElement);
      if(absDmatElement == 0)
	continue;
      // OK, we have a non-zero density matrix element.
      // index pair (i,j)
      for(int ii = 0; ii < noOfNeighborsList[i]; ii++) {
	const neighbor_basisfunc_struct* currNeighbor_i = &neighborList[i*maxNoOfNeighbors+ii];
	int index_i = currNeighbor_i->i;
	ergo_real max_CS_factor_i = currNeighbor_i->max_CS_factor;
	for(int jj = 0; jj < noOfNeighborsList[j]; jj++) {
	  const neighbor_basisfunc_struct* currNeighbor_j = &neighborList[j*maxNoOfNeighbors+jj];
	  int index_j = currNeighbor_j->i;
	  if(index_j < index_i)
	    continue;
	  ergo_real max_CS_factor_j = currNeighbor_j->max_CS_factor;
	  ergo_real maxContrib = absDmatElement * max_CS_factor_i * max_CS_factor_j;
	  if(maxContrib < threshold)
	    break;
	  // OK, this index pair must be included.
	  // Check if this index pair is already in long list or in short list.
	  int foundInLongList = find_int_in_sorted_list(longList[index_i], longListCounterList[index_i], index_j);
	  if(foundInLongList)
	    continue;
	  int foundInShortList = find_int_in_list(&shortList[index_i*sll], shortListCounterList[index_i], index_j);
	  if(foundInShortList)
	    continue;
	  // not found, this is a new index pair, must be added.
	  // check if shortList is full.
	  if(shortListCounterList[index_i] == sll) {
	    // allocate new long list
	    int* newLongList = new int[longListCounterList[index_i]+sll];
	    // Merge short list with long list.
	    do_sort_int_list(&shortList[sll*index_i], sll);
	    if(find_doubles_in_sorted_list(longList[index_i], longListCounterList[index_i])
	       || find_doubles_in_sorted_list(&shortList[sll*index_i], sll) ) {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error: doubles found before merge! (1) or (2).");
	      delete [] newLongList; /* Added this to satisfy cppcheck. */
	      return -1;
	    }
	    merge_lists(longList[index_i], longListCounterList[index_i],
			&shortList[sll*index_i], sll, newLongList);
	    if(find_doubles_in_sorted_list(newLongList, longListCounterList[index_i] + sll)) {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error: doubles found after merge!");
	      return -1;
	    }
	    if(longList[index_i])
	      delete [] longList[index_i];
	    longList[index_i] = newLongList;
	    longListCounterList[index_i] += sll;
	    shortListCounterList[index_i] = 0;
	  } // END IF short list full
	  // OK, now we know the short list is not full, so we can safely add the new column index.
	  shortList[sll*index_i+shortListCounterList[index_i]] = index_j;
	  shortListCounterList[index_i]++;
	} // END FOR jj
      } // END FOR ii
    } // end for idx
  } // end for r

  // Now merge contents of short list to long list for each row.
  for(int i = 0; i < n; i++)
    {
      // allocate new long list
      int* newLongList = new int[longListCounterList[i]+shortListCounterList[i]];
      // Merge short list with long list.
      do_sort_int_list(&shortList[sll*i], shortListCounterList[i]);
      merge_lists(longList[i], longListCounterList[i], &shortList[sll*i], shortListCounterList[i], newLongList);
      if(longList[i])
	delete [] longList[i];
      longList[i] = newLongList;
      longListCounterList[i] += shortListCounterList[i];
      shortListCounterList[i] = 0;
    }

  return 0;
}


typedef int* intPtr;

struct listsStruct {
  int** longList;
  int* longListCounterList;
  int n_saved;
  listsStruct() : longList(0), longListCounterList(0), n_saved(0) { }
  ~listsStruct() { 
    for(int i = 0; i < n_saved; i++) {
      if(longList[i])
	delete [] longList[i];
    }
    delete [] longList;
    delete [] longListCounterList;
  }
  void init(int n) {
    longList = new intPtr[n];
    longListCounterList = new int[n];
    memset(longList, 0, n*sizeof(int*));
    memset(longListCounterList, 0, n*sizeof(int));
    n_saved= n;
  }
};


struct identify_needed_elements_thread_struct {
  pthread_t thread;
  int thread_ID;
  int nThreads;
  ergo_real threshold;
  const csr_matrix_struct* dens_CSR;
  const int* noOfNeighborsList;
  const neighbor_basisfunc_struct* neighborList;
  int maxNoOfNeighbors;
  int ** longList;
  int *longListCounterList;
  int resultCode;
  identify_needed_elements_thread_struct(const csr_matrix_struct* dens_CSR_,
					 const int* noOfNeighborsList_,
					 const neighbor_basisfunc_struct* neighborList_)
    : dens_CSR(dens_CSR_),
      noOfNeighborsList(noOfNeighborsList_),
      neighborList(neighborList_),
      resultCode(-1)
  { }
};


static void*
identify_needed_elements_thread_func(void* arg) {
  identify_needed_elements_thread_struct* params = (identify_needed_elements_thread_struct*)arg;
  params->resultCode = identify_needed_elements_part(params->threshold, 
						     params->dens_CSR, 
						     params->noOfNeighborsList,
						     params->neighborList, 
						     params->maxNoOfNeighbors,
						     params->longList,
						     params->longListCounterList,
						     params->thread_ID,
						     params->nThreads);
  return NULL;
}


/** Tries to predict which elements of K will be needed.
    Use two different lists, a "long list" and a "short list" for each row of K.
    The "long list" is always sorted so that we can quickly check if a column index is already present.
*/
static int
identify_needed_elements(ergo_real threshold, 
			 const csr_matrix_struct* dens_CSR, 
			 const int noOfNeighborsList[],
                         const neighbor_basisfunc_struct* neighborList, 
			 int maxNoOfNeighbors,
                         int ** longList, 
			 int *longListCounterList,
			 int nThreads)
{
  if(nThreads == 1) {
    // Single-thread case.
    return identify_needed_elements_part(threshold, 
					 dens_CSR, 
					 noOfNeighborsList,
					 neighborList, 
					 maxNoOfNeighbors,
					 longList,
					 longListCounterList,
					 0,
					 1);
  }
  // Multi-thread case.
  const int n = dens_CSR->n;
  std::vector<listsStruct> listList(nThreads);
  for(int i = 0; i < nThreads; i++)
    listList[i].init(n);

  identify_needed_elements_thread_struct* threadParamsList[nThreads];

  // Set common parameters for all threads
  for(int i = 0; i < nThreads; i++)
    {
      threadParamsList[i] = 
	new identify_needed_elements_thread_struct(dens_CSR, noOfNeighborsList, neighborList);
      threadParamsList[i]->threshold = threshold;
      threadParamsList[i]->maxNoOfNeighbors = maxNoOfNeighbors;
      threadParamsList[i]->longList = listList[i].longList;
      threadParamsList[i]->longListCounterList = listList[i].longListCounterList;
      threadParamsList[i]->thread_ID = i;
      threadParamsList[i]->nThreads = nThreads;
    } // END FOR i

  /* start threads */
  for(int i = 0; i < nThreads; i++)
    {
      if(pthread_create(&threadParamsList[i]->thread, 
			NULL, 
			identify_needed_elements_thread_func, 
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

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "%i threads started OK.", nThreads);

  /* wait for threads to finish */
  for(int i = 0; i < nThreads; i++)
    {
      if(pthread_join(threadParamsList[i]->thread, NULL) != 0)
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in pthread_join for thread %i", i);
    } /* END FOR i */
  
  /* now all threads have finished, check for errors */
  for(int i = 0; i < nThreads; i++)
    {
      if(threadParamsList[i]->resultCode != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in identify_needed_elements_thread_func."
		    " for thread %i", i);
	  return -1;
	}
    } /* END FOR i */  
    
  // Now merge contents of lists to long list for each row.
  for(int i = 0; i < n; i++)
    {
      for(int ti = 0; ti < nThreads; ti++) {
	int count = longListCounterList[i] + listList[ti].longListCounterList[i];
	int* newLongList = new int[count];
	merge_lists(longList[i], longListCounterList[i], 
		    listList[ti].longList[i], listList[ti].longListCounterList[i], newLongList);
	// Now newLongList probably contains duplicates. Count unique entries.
	int nUnique = 0;
	if(count > 0)
	  nUnique++;
	for(int k = 0; k < count-1; k++) {
	  if(newLongList[k] > newLongList[k+1]) {
	    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, 
		      "Error in identify_needed_elements: list not sorted.");
	    delete [] newLongList; /* Added this to satisfy cppcheck. */
	    return -1;
	  }
	  if(newLongList[k] < newLongList[k+1])
	    nUnique++;
	}
	if(longList[i])
	  delete [] longList[i];
	longList[i] = new int[nUnique];
	longListCounterList[i] = nUnique;
	if(count > 0)
	  longList[i][0] = newLongList[0];
	int nUnique2 = 0;
	if(count > 0)
	  nUnique2++;
	for(int k = 0; k < count-1; k++) {
	  if(newLongList[k] < newLongList[k+1])
	    longList[i][nUnique2++] = newLongList[k+1];
	}
	delete [] newLongList;
      }
    }

  // Cleanup
  for(int i = 0; i < nThreads; i++)
    delete threadParamsList[i];

  return 0;
}





struct distr_idxs_and_factor_struct {
  int i1;
  int i2;
  ergo_real limitingFactor;
};

static int
compare_distr_idxs_and_factor_structs(const void* p1in, const void* p2in) {
  const distr_idxs_and_factor_struct* p1 = (const distr_idxs_and_factor_struct*)p1in;
  const distr_idxs_and_factor_struct* p2 = (const distr_idxs_and_factor_struct*)p2in;
  if(p1->i1 > p2->i1)
    return 1;
  if(p1->i1 < p2->i1)
    return -1;
  if(p1->i2 > p2->i2)
    return 1;
  if(p1->i2 < p2->i2)
    return -1;
  if(p1->limitingFactor < p2->limitingFactor)
    return 1;
  if(p1->limitingFactor > p2->limitingFactor)
    return -1;
  return 0;
}


int
create_CSR_for_K(const BasisInfoStruct & basisInfo,
		 const IntegralInfo & integralInfo,
		 const JK::Params& J_K_params,
		 csr_matrix_struct* dens_CSR,
		 csr_matrix_struct* K_CSR,
		 int symmetryFlag)
{
  Util::TimeMeter timeMeterTot;
  Util::TimeMeter timeMeterDistrs;

  int n = basisInfo.noOfBasisFuncs;

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "entering create_CSR_for_K, no of basis funcs = %5i, threshold_K = %7.3g", 
	    n, (double)J_K_params.threshold_K);

  output_current_memory_usage(LOG_AREA_INTEGRALS, "beginning of create_CSR_for_K");

  // compute list of distributions, with max_CS_factor for each distr.

  ergo_real maxDensityMatrixElement = ergo_CSR_get_max_abs_element(dens_CSR);


  // get largest limiting factor
  ergo_real maxLimitingFactor = 0;
  if(get_list_of_labeled_distrs_maxLimitingFactor(basisInfo,
						  integralInfo,
						  J_K_params.threshold_K,
						  &maxLimitingFactor,
						  maxDensityMatrixElement) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_list_of_labeled_distrs_maxLimitingFactor");
      return -1;
    }

  // Get number of distributions
  int distrCount = get_list_of_labeled_distrs(basisInfo,
					      integralInfo,
					      J_K_params.threshold_K,
					      NULL,
					      0,
					      maxLimitingFactor,
					      NULL,
					      maxDensityMatrixElement);
  if(distrCount < 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_CSR_for_K: (distrCount < 0)");
      return -1;
    }

  std::vector<DistributionSpecStructLabeled> distrList(distrCount);

  // create list of product primitives, with labels
  int distrCountTemp = get_list_of_labeled_distrs(basisInfo,
						  integralInfo,
						  J_K_params.threshold_K,
						  &distrList[0],
						  distrCount,
						  maxLimitingFactor,
						  NULL,
						  maxDensityMatrixElement);
  if(distrCountTemp != distrCount)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_CSR_for_K:(distrCountTemp != distrCount)");
      return -1;
    }

  output_current_memory_usage(LOG_AREA_INTEGRALS, "after allocating list of primitive distributions");


  // Create simple list of distrs, where each entry in the list only
  // holds the two indexes and the limitingFactor.
  std::vector<distr_idxs_and_factor_struct> distrListSimple(distrCount*2);
  int distrCount2 = 0;
  for(int i = 0; i < distrCount; i++) {
    DistributionSpecStructLabeled* curr = &distrList[i];
    distrListSimple[distrCount2].i1 = curr->basisFuncIndex_1;
    distrListSimple[distrCount2].i2 = curr->basisFuncIndex_2;
    distrListSimple[distrCount2].limitingFactor = curr->limitingFactor;
    distrCount2++;
    if(curr->basisFuncIndex_1 != curr->basisFuncIndex_2) {
      distrListSimple[distrCount2].i1 = curr->basisFuncIndex_2;
      distrListSimple[distrCount2].i2 = curr->basisFuncIndex_1;
      distrListSimple[distrCount2].limitingFactor = curr->limitingFactor;
      distrCount2++;
    }
  }
  // Sort list according to indexes and factor.
  qsort(&distrListSimple[0], distrCount2, sizeof(distr_idxs_and_factor_struct), compare_distr_idxs_and_factor_structs);

  // Check that list of distrs is properly sorted.
  for(int i = 0; i < distrCount2-1; i++) {
    // compare this distr with next one.
    distr_idxs_and_factor_struct* curr = &distrListSimple[i];
    distr_idxs_and_factor_struct* next = &distrListSimple[i+1];
    if(next->i1 < curr->i1) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_CSR_for_K: distr list not sorted.");
      return -1;
    }
    if(next->i1 == curr->i1) {
      if(next->i2 < curr->i2) {
        do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_CSR_for_K: distr list not sorted.");
        return -1;
      }
      if(next->i2 == curr->i2) {
        if(next->limitingFactor > curr->limitingFactor) {
          do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_CSR_for_K: distr list not sorted.");
          return -1;
        }
      }
    }
  }

  timeMeterDistrs.print(LOG_AREA_INTEGRALS, "create_CSR_for_K getting list of distrs");

  // Now use sorted list to get number of neighbors for each basis function.
  int noOfNeighborsList[n];
  for(int i = 0; i < n; i++)
    noOfNeighborsList[i] = 0;
  int iWhile = 0;
  int basisFuncIndex = -1;
  while(iWhile < distrCount2) {
    distr_idxs_and_factor_struct* curr = &distrListSimple[iWhile];
    // Now find last entry for this i1 value.
    int ii = iWhile;
    while(ii < distrCount2) {
      distr_idxs_and_factor_struct* curr2 = &distrListSimple[ii];
      if(curr2->i1 != curr->i1)
        break;
      ii++;
    }
    int noOfDistrsCurrIdx = ii - iWhile;
    int noOfNeighbors = 1;
    for(int k = 0; k < noOfDistrsCurrIdx-1; k++) {
      if(distrListSimple[iWhile+k].i2 != distrListSimple[iWhile+k+1].i2)
        noOfNeighbors++;
    }
    noOfNeighborsList[curr->i1] = noOfNeighbors;
    iWhile += noOfDistrsCurrIdx;
  }
  int maxNoOfNeighbors = 0;
  for(int i = 0; i < n; i++) {
    if(noOfNeighborsList[i] > maxNoOfNeighbors)
      maxNoOfNeighbors = noOfNeighborsList[i];
  }

  Util::TimeMeter timeMeterNeighbors;


  // Create list of neighbor basis funcs, with max_CS_factor for each.
  std::vector<neighbor_basisfunc_struct> neighborList(n*maxNoOfNeighbors);
  for(int i = 0; i < n*maxNoOfNeighbors; i++) {
    neighborList[i].i = -1;
    neighborList[i].max_CS_factor = 0;
  }
  iWhile = 0;
  basisFuncIndex = -1;
  while(iWhile < distrCount2) {
    // now i should point to start of a new basis func
    basisFuncIndex = distrListSimple[iWhile].i1;
    int noOfNeighbors = 0;
    // get all neighbors of this basis func.
    while(iWhile < distrCount2 && distrListSimple[iWhile].i1 == basisFuncIndex) {
      // now i should point to a new neighbor basis func.
      int neighborIndex = distrListSimple[iWhile].i2;
      neighbor_basisfunc_struct* currNeighbor = &neighborList[basisFuncIndex*maxNoOfNeighbors+noOfNeighbors];
      currNeighbor->i = neighborIndex;
      ergo_real max_CS_factor = 0;
      while(iWhile < distrCount2 && distrListSimple[iWhile].i1 == basisFuncIndex && distrListSimple[iWhile].i2 == neighborIndex) {
        ergo_real CS_factor = distrListSimple[iWhile].limitingFactor;
        if(CS_factor > max_CS_factor)
          max_CS_factor = CS_factor;
        iWhile++;
      }
      currNeighbor->max_CS_factor = max_CS_factor;
      noOfNeighbors++;
    }
    if(noOfNeighbors != noOfNeighborsList[basisFuncIndex]) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_CSR_for_K: noOfNeighbors mismatch.");
      return -1;
    }      
  }


  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "list of neighbors complete.");

  timeMeterNeighbors.print(LOG_AREA_INTEGRALS, "create_CSR_for_K getting list of neighbors");




  // Sort each list of neighbors by max_CS_factor, so we can skip out of loops earlier.
  Util::TimeMeter timeMeterSortNeighbors;
  for(int i = 0; i < n; i++)
    {
      neighbor_basisfunc_struct* currList = &neighborList[i*maxNoOfNeighbors]; 
      int count = noOfNeighborsList[i];
      // Bubble sort
      for(int ii = 0; ii < count-1; ii++)
        for(int jj = 0; jj < count-1-ii; jj++) 
          {
            if(currList[jj].max_CS_factor < currList[jj+1].max_CS_factor)
              {
                // switch
                neighbor_basisfunc_struct temp = currList[jj];
                currList[jj] = currList[jj+1];
                currList[jj+1] = temp;
              }
          }
      // check that list is really sorted
      for(int ii = 0; ii < count-1; ii++)
        {
          if(currList[ii].max_CS_factor < currList[ii+1].max_CS_factor)
            {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "ERROR: neighbor list not sorted!");
	      return -1;
            }
        }
    } // END FOR i sort each list of neighbors
  timeMeterSortNeighbors.print(LOG_AREA_INTEGRALS, "create_CSR_for_K sort each list of neighbors");


  Util::TimeMeter timeMeterIdentifyNeededElements;

  std::vector<int*> longList(n);
  std::vector<int> longListCounterList(n);
  memset(&longList[0], 0, n*sizeof(int*));
  memset(&longListCounterList[0], 0, n*sizeof(int));

  if(identify_needed_elements(J_K_params.threshold_K, dens_CSR, noOfNeighborsList,
                              &neighborList[0], maxNoOfNeighbors,
                              &longList[0], &longListCounterList[0], J_K_params.noOfThreads_K) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in identify_needed_elements.");
      return -1;
    }

  timeMeterIdentifyNeededElements.print(LOG_AREA_INTEGRALS, "identify_needed_elements()");

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "computing nnz..");

  // compute nnz;
  size_t nnz = 0;
  for(int i = 0; i < n; i++)
    nnz += longListCounterList[i];

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "create_CSR_for_K: predicted nnz = %12i, <-> %6.2f%% of a full matrix", 
	    nnz, ((double)nnz*100) / ((double)n*n));

  Util::TimeMeter timeMeterLastPart;

  std::vector<int> rowind(nnz);
  std::vector<int> colind(nnz);
  int count = 0;
  for(int i = 0; i < n; i++)
    {
      for(int j = 0; j < longListCounterList[i]; j++)
	{
	  rowind[count] = i;
	  colind[count] = longList[i][j];
	  count++;
	}
    }
  // Now all info we need is in the vectors rowind and colind.
  // Free all other memory.
  distrList.clear();
  neighborList.clear();
  for(int i = 0; i < n; i++)
    {
      if(longList[i])
	{
	  delete [] longList[i];
	  longList[i] = NULL;
	}
    }
  
  // Now use the vectors rowind and colind to create CSR structure.
  if(symmetryFlag == 1) {
    // Symmetric case.
    if(ergo_CSR_create(K_CSR, 
		       1,
		       n,
		       nnz,
		       &rowind[0],
		       &colind[0]) != 0)
      {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in ergo_CSR_create for density matrix.");
	return -1;
      }
  }
  else {
    // Non-symmetric case.
    std::vector<int> rowind2(nnz*2);
    std::vector<int> colind2(nnz*2);
    size_t count = 0;
    for(size_t i = 0; i < nnz; i++) {
      int row = rowind[i];
      int col = colind[i];
      rowind2[count] = row;
      colind2[count] = col;
      count++;
      if(row != col) {
        rowind2[count] = col;
        colind2[count] = row;
        count++;
      }
    }
    size_t nnz2 = count;
    if(ergo_CSR_create(K_CSR, 
                       0,
                       n,
                       nnz2,
                       &rowind2[0],
                       &colind2[0]) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in ergo_CSR_create for density matrix.");
      return -1;
    }
  }
  
  
  timeMeterLastPart.print(LOG_AREA_INTEGRALS, "create_CSR_for_K last part");

  timeMeterTot.print(LOG_AREA_INTEGRALS, "create_CSR_for_K total");

  return 0;  
}


