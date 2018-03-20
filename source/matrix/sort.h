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

/** @file sort.h

    \brief A quicksort implementation.
*/

#ifndef MAT_SORT
#define MAT_SORT
namespace mat {

#if 1

  template<class Treal>
    void quicksort(const Treal* value, int* index, int low, int high) 
    {
    if(high >= low)
      {
	int i = low;
	int j = high;
	int tmp;
	Treal pivot = value[index[(low + high) / 2]];     /* Introduce the pivot */
	do {                                     /* Permute elements so that all */
	  while (value[index[i]] < pivot) i++;   /* elements end up in one of    */
	  while (value[index[j]] > pivot) j--;   /* two groups, one where the    */
	  if (i <= j) {                        /* elements have a value in value */
	    tmp = index[i];                      /* smaller than the pivot and   */
	    index[i] = index[j];                 /* one group with a value larger*/
	    index[j] = tmp;                      /* than the pivot               */
	    i++;
	    j--;
	  }
	} while (i <= j);
	if(low < j)  quicksort(value, index, low, j); /* Sort the two groups     */
	if(i < high) quicksort(value, index, i, high);
      }
  }

#else

  
  template<typename Treal, typename Tfun>
    void quicksort(Tfun const & fun, int* index, int low, int high) 
    throw(std::exception){
    int i = low;
    int j = high;
    int tmp;
    Treal pivot = fun.value(index[(low + high) / 2]); /* Introduce pivot */
    do {                                     /* Permute elements so that all */
      while (fun.value(index[i]) < pivot) i++;/* elements end up in one of*/
      while (fun.value(index[j]) > pivot) j--;/* two groups, one where the*/
      if (i <= j) {                        /* elements have a value in value */
	tmp = index[i];                      /* smaller than the pivot and   */
	index[i] = index[j];                 /* one group with a value larger*/
	index[j] = tmp;                      /* than the pivot               */
	i++;
	j--;
      }
    } while (i <= j);
    /* Sort the two groups     */
    if(low < j)  quicksort<Treal>(fun, index, low, j); 
    if(i < high) quicksort<Treal>(fun, index, i, high);
  }

  template<class Treal>
    void quicksort(const Treal* value, int* index, int low, int high) {
    int i = low;
    int j = high;
    int tmp;
    Treal pivot = value[index[(low + high) / 2]];     /* Introduce the pivot */
    do {                                     /* Permute elements so that all */
      while (value[index[i]] < pivot) i++;   /* elements end up in one of    */
      while (value[index[j]] > pivot) j--;   /* two groups, one where the    */
      if (i <= j) {                        /* elements have a value in value */
	tmp = index[i];                      /* smaller than the pivot and   */
	index[i] = index[j];                 /* one group with a value larger*/
	index[j] = tmp;                      /* than the pivot               */
	i++;
	j--;
      }
    } while (i <= j);
    if(low < j)  quicksort(value, index, low, j); /* Sort the two groups     */
    if(i < high) quicksort(value, index, i, high);
  }
  
#endif
  
} /* end namespace mat */

#endif
