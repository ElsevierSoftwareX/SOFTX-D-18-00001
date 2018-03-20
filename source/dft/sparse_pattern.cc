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

/** @file sparse_pattern.cc

    @brief Class that can be used to store sparse matrix patterns.

    @author: Pawel Salek <em>responsible</em>
*/

#include "output.h"
#include "dft_common.h"
#include "sparse_pattern.h"

BEGIN_NAMESPACE(Dft)

/** Add interval { i: lo <= i < hi } to the list. The list is specific
    to given column. */

void SparsePattern::Column::addInterval(int lo, int hi)
{
  if(list.empty()) {
    list.push_back(SparsePattern::Interval(lo, hi));
    return;
  }

  /* There are four cases: The interval is disjoint to its neighbours:
     insert it; It can overlap with it precessor: extend the
     precessor; it can overlap with successor only: extend the
     successor. It can also overlap with many intervals: extend the
     first one and keep removing following until a disjoint or the end
     of list is found. */
  for(SparsePattern::IntervalList::iterator i = list.begin();
      i != list.end(); ++i) {
    if(hi < i->lo) {
      list.insert(i,SparsePattern::Interval(lo, hi));
      return;
    }
    if(lo <= i->hi) { /* Here we go, first one in the chain localized! */

      if(lo < i->lo)
        i->lo = lo;

      if(hi<= i->hi) /* Nothing left to do. */
        return;

      /* Now, we have to only figure out where it ends... */
      do {
        SparsePattern::IntervalList::iterator j = i;
        ++j;
        if(j == list.end() || /* There is no next one, or... */
           hi < j->lo) {      /* It's too far up. */
          i->hi = hi;
          return;
        }
        /* OK, it needs to be joined... */
        if(hi <= j->hi) {
          i->hi = j->hi;
          /* here, the story ends... */
          list.erase(j);
          return;
        }
        /* hi apparently goes higher... */
        list.erase(j);
      } while(true);
      return;
    }
  }
  /* We reached the end of list and still no hit, time to append
     stuff. */
  list.push_back(SparsePattern::Interval(lo, hi));
}

void
SparsePattern::Column::addIntervals(int nIntervals, int (*intervals)[2])
{
  int currentInterval = 0;
  SparsePattern::IntervalList::iterator i = list.begin();

  /* There are four cases: The interval is disjoint to its neighbours:
     insert it; It can overlap with it precessor: extend the
     precessor; it can overlap with successor only: extend the
     successor. It can also overlap with many intervals: extend the
     first one and keep removing following until a disjoint or the end
     of list is found. */
  // printf("Begin, nIntervals = %d current pattern length %d\n",  nIntervals, list.size());
  
  while(i != list.end() && currentInterval<nIntervals) {  
    int lo = intervals[currentInterval][0];
    int hi = intervals[currentInterval][1];
    //printf("Begin loop current interval %d %d\n", lo, hi);
    if(hi < i->lo) {
      i = list.insert(i,SparsePattern::Interval(lo, hi));
      ++currentInterval;
      continue;
    }
    if(lo <= i->hi) { /* Here we go, first one in the chain localized! */

      if(lo < i->lo)
        i->lo = lo;

      if(hi<= i->hi) { /* Nothing left to do for this one. */
        ++currentInterval;
        continue;
      }

      /* Now, we have to only figure out where it ends... */
      do {
        SparsePattern::IntervalList::iterator j = i;
        ++j;
        if(j == list.end() || /* There is no next one, or... */
           hi < j->lo) {      /* It's too far up. */
          i->hi = hi;
          ++currentInterval;
          break;
        }
        /* OK, it needs to be joined... */
        if(hi <= j->hi) {
          i->hi = j->hi;
          /* here, the story ends... */
          list.erase(j);
          ++currentInterval;
          break;
        }
        /* hi apparently goes higher... */
        list.erase(j);
      } while(true);
      continue;
    }

    ++i;
  }

  /* We reached the end of list and still no hit, time to append
     the remaining stuff. */
  while(currentInterval<nIntervals) {
    list.push_back(SparsePattern::Interval(intervals[currentInterval][0],
                                           intervals[currentInterval][1]));
    ++currentInterval;
  }
}

void
SparsePattern::add(int nRanges, const int (*shellRanges)[2])
{
  /* Works in a number of steps:
     a) translate shell ranges to basis function ranges.
  */
  int nOrbs;
  int (*orbRanges)[2] = new int[nRanges][2];

  ergoShellsToOrbs(&nRanges, shellRanges, &nOrbs, orbRanges, bis);

  for(int colBlock=0; colBlock<nOrbs; colBlock++) 
    for(int col=orbRanges[colBlock][0]; col<orbRanges[colBlock][1]; col++) {
      SparsePattern::Column& column = ranges[col];
      //printf("Adding intervals to col %d\n", col);
#if 0
      for(int rowBlock=0; rowBlock<nOrbs; rowBlock++) {
        column.addInterval(orbRanges[rowBlock][0], orbRanges[rowBlock][1]);
      }
#else
      column.addIntervals(nOrbs, orbRanges);
#endif
    }

  delete []orbRanges;
}

/** Load itself from the specified stream. */
void SparsePattern::load(FILE *f)
{
  int nBasis;
  Interval tmp(0,0);

  if(fread(&nBasis, sizeof(int), 1, f) != 1)
    throw "SparsePattern::load, point 1";
  if(nBasis != bis.noOfBasisFuncs)
    throw "SparsePattern::load, size misalignment";
  for(int i=0; i<nBasis; i++) {
    IntervalList& list = ranges[i].list;
    int intervalCnt;
    if(fread(&intervalCnt, sizeof(int), 1, f) != 1)
      throw "Sparse::Pattern::load, interval cnt read";
    list.clear();
    for(int interval=0; interval<intervalCnt; interval++) {
      if(fread(&tmp, sizeof(Interval), 1, f) != 1)
        throw "Sparse::Pattern::load, interval read";
      list.push_back(tmp);
    }
  }

  int sumSize = 0, maxWidth=0;
  for(int i=0; i<bis.noOfBasisFuncs; i++) {
    int width = getColumnSize(i);
    sumSize += width;
    if(width>maxWidth)
      maxWidth = width;
  }
#if 0
  printf("Read sparse pattern has %d elemts. Width max %d, avg. %4d Size: %f G\n",
         sumSize, maxWidth, int(sumSize/float(bis.noOfBasisFuncs)),
         sumSize*double(sizeof(ergo_real)/(1024.0*1024.0*1024.0)));
#endif
}

/** Save itself to the specified stream. */
void SparsePattern::save(FILE *f) const
{
  int sumSize = 0, maxWidth=0;
  for(int i=0; i<bis.noOfBasisFuncs; i++) {
    int width = getColumnSize(i);
    sumSize += width;
    if(width>maxWidth)
      maxWidth = width;
  }
#if 0
  printf("Saved sparse pattern has %d elemts. Width max %d, avg. %4d Size: %f G\n",
         sumSize, maxWidth, int(sumSize/float(bis.noOfBasisFuncs)),
         sumSize*double(sizeof(ergo_real)/(1024.0*1024.0*1024.0)));
#endif         
  do_output(LOG_CAT_INFO, LOG_AREA_DFT, 
            "Sparse pattern has %d elemts. Width max %d, avg. %4d Size: %f G",
            sumSize, maxWidth, int(sumSize/float(bis.noOfBasisFuncs)),
            (double)(sumSize*double(sizeof(ergo_real)/(1024.0*1024.0*1024.0))));

  if(fwrite(&bis.noOfBasisFuncs, sizeof(int), 1, f) != 1)
    throw "Cannot save sparsity pattern";
  for(int col=0; col<bis.noOfBasisFuncs; col++) {
    IntervalList& l = ranges[col].list;
#if 1
    int cnt = l.size();
    if(fwrite(&cnt, sizeof(int), 1, f) != 1) throw "Save size";
    for(IntervalList::const_iterator i = l.begin();
        i != l.end(); ++i) {
      if(fwrite( &(*i), sizeof(Interval), 1, f) != 1) throw "Save interval";
    }
#else
    int cnt = 1;
    if(fwrite(&cnt, sizeof(int), 1, f) != 1) throw "Save size";
    //Interval i(l.begin()->lo, l.rbegin()->hi);
    Interval i(0, bis.noOfBasisFuncs);
    
    if(fwrite( &i, sizeof(Interval), 1, f) != 1) throw "Save interval";
#endif
  }
}

int
SparsePattern::sizeTotal() const
{
  int sumSize = 0;
  for(int i=0; i<bis.noOfBasisFuncs; i++) {
    int width = getColumnSize(i);
    sumSize += width;
  }
  return sumSize;
}

/** Prepares the AO map given a shell map. */
static void
prepareAOMap(const BasisInfoStruct& bis, const int *shellMap, int *aoMap)
{
  int newIdx = 0;
  for(int i=0; i< bis.noOfShells; i++) {
    const ShellSpecStruct& shell = bis.shellList[shellMap[i]];
    for(int ao=0; ao<shell.noOfBasisFuncs; ao++) {
      aoMap[newIdx+ao] = shell.startIndexInMatrix+ao;
      //printf("%3d %3d\n", newIdx+ao, aoMap[newIdx+ao]);
    }
    newIdx += shell.noOfBasisFuncs;
  }
}

/** prepares a shell map and matching AO map permuted to optimize
    performance of the XC code. We cannot reuse the permutations used
    in the matrix library because they can are AO based and not shell
    based.
    
    @param shellMap - previously allocated vector that will be filled
    with the permutation data.

    @param aoMap - corresponding AO permutation vector, preallocated.

    The code uses a variant of the Cuthill-McKee algorithm to
    determine the shell map. The AO map is trivially generated from
    the shell map (perhaps it could be a separate function?).
    Generation of the optimal performance is in general a complex
    matter but since we use a discrete selection criteria to determine
    the shell radius, Cuthill-McKee will do.
*/
#define USE_CUTHILL_MCKEE 0
#if USE_CUTHILL_MCKEE
void
setupShellMap(const BasisInfoStruct& bis, int *shellMap, int *aoMap)
{
  static const ergo_real THR=1e-2;
  /* BEGIN setupShellMap */
  int nShells = bis.noOfShells;
  std::vector<int> result;
  
  /* prepare lists of neighbours */
  std::vector<NeighbourList> shellNeighbours;

  for(int iShell=0; iShell<bis.noOfShells; iShell++)
    shellNeighbours.push_back(NeighbourList(&bis.shellList[iShell], THR));

  /* WARNING: suspicious quadratic loop. */
  for(int iShell=0; iShell<bis.noOfShells; iShell++)
    shellNeighbours[iShell].setOverlappingWith(shellNeighbours);

  /* Core of Cuthill-McKee  */
  int icore = findExtremeShell(bis.noOfShells, bis.shellList);
  std::set<int> processedShells;
  result.push_back(icore);
  processedShells.insert(icore);
  assert(result.size() >0);

  int perhapsAvailableShellNo=0;
  for(int element=1; element<nShells; element++) {
    int iShell = result.at(element-1);
    NeighbourList& neighbours = shellNeighbours.at(iShell);

    /* Find elements that have not been processed yet. */

    std::list< std::pair<ergo_real,int> > unprocessedNeighbours;
    for(std::list<int>::iterator i= neighbours.begin();
        i != neighbours.end(); i++) {
      if( processedShells.find(*i) == processedShells.end()) {
        ergo_real r = std::sqrt(sqDist(bis.shellList[iShell].centerCoords,
                                       bis.shellList[*i].centerCoords));
        unprocessedNeighbours.push_back( std::pair<ergo_real,int>(r,*i) );
      }
    }
#if 0
    printf("\nShell %d has %d neighbours, %d unprocessed. Total processed: %d\n",
           iShell, neighbours.size(), unprocessedNeighbours.size(),
           processedShells.size());
#endif
    /* Sort wrt the shell radius */
    unprocessedNeighbours.sort();
    
    for(std::list< std::pair<ergo_real,int> >::iterator
          i = unprocessedNeighbours.begin();
        i != unprocessedNeighbours.end(); i++) {
      result.push_back(i->second);
      processedShells.insert(i->second);
    }

    /* Make sure we did not hit a disjoint element here... */
    if(result.size()<=(unsigned)element) {
      for(;perhapsAvailableShellNo<nShells; perhapsAvailableShellNo++) {
        if( processedShells.find(perhapsAvailableShellNo) ==
            processedShells.end() ) {
          result.push_back(perhapsAvailableShellNo);
          processedShells.insert(perhapsAvailableShellNo);
          ++perhapsAvailableShellNo;
          break;
        }
      }
    }
    
  } /* END FOR element<nShells */

  assert(result.size() == (unsigned)bis.noOfShells);
  std::copy(result.begin(), result.end(), shellMap);
  
  prepareAOMap(bis, shellMap, aoMap);
  //printf("F: Cuthill-McKee is done.\n");
}
#else
static void
clusterShells(const ShellSpecStruct *shells, const Box& box,
              const std::vector<int>& inputList,
              std::vector<int>& result, int depth)
{
  static const ergo_real BOX_SIZE = 1.5;

  if(inputList.empty())
    return;

  int dividingDim = box.getMaxDim();
  ergo_real dividingSize = box.hi[dividingDim]-box.lo[dividingDim];
  

  int n = inputList.size();
  if(n<=2 || dividingSize < BOX_SIZE) {
    for(std::vector<int>::const_iterator i= inputList.begin();
        i != inputList.end(); ++i) {
      result.push_back(*i);
    }
    return;
  }

  std::vector<int> lessThanList, greaterList;

  lessThanList.reserve(n);
  greaterList.reserve(n);

  ergo_real dividingValue = 0.5*(box.hi[dividingDim]+box.lo[dividingDim]);
  for(std::vector<int>::const_iterator i= inputList.begin();
      i != inputList.end(); ++i) {
    int shellIdx = *i;
    if(shells[shellIdx].centerCoords[dividingDim] < dividingValue)
      lessThanList.push_back(shellIdx);
    else
      greaterList.push_back(shellIdx);
  }

  ++depth;
  Box bb(box);
  bb.hi[dividingDim] = dividingValue;
  clusterShells(shells, bb, lessThanList, result, depth);
  bb.lo[dividingDim] = dividingValue;
  bb.hi[dividingDim] = box.hi[dividingDim];
  clusterShells(shells, bb, greaterList, result, depth);
}

void
setupShellMap(const BasisInfoStruct& bis, int *shellMap, int *aoMap)
{
  /* this uses a tree based clustering... */
  Box box;
  unsigned i;
  for(i=0; i<3; i++) {
    box.lo[i] = box.hi[i] = bis.shellList[0].centerCoords[i];
  }

  std::vector<int> inputList(bis.noOfShells);
  for(int shell=0; shell< bis.noOfShells; shell++) {
    inputList[shell] = shell;
    /* modify box here */
    for(i=0; i<3; i++) {
      ergo_real c = bis.shellList[shell].centerCoords[i];

      if(c < box.lo[i])
        box.lo[i] = c;
      else if(c > box.hi[i])
        box.hi[i] = c;
    }
  }

  std::vector<int> result;
  result.reserve(bis.noOfShells);
  clusterShells(bis.shellList, box, inputList, result, true);

  assert(result.size() == (unsigned)bis.noOfShells);
  std::copy(result.begin(), result.end(), shellMap);

  prepareAOMap(bis, shellMap, aoMap);
}
#endif

END_NAMESPACE(Dft)
