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

/** @file box_system.cc

    @brief BoxSystem class representing a hierarchical data structure
    of boxes in 3D space (an oct-tree).

    The idea is that you have a list of items at different points in space,
    and you want a hierarchical system of boxes containing those items.

    You give a list of items, and the function create_box_system will
    create a system of boxes for you.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <cmath>
#include <stdlib.h>
#include <vector>
#include "box_system.h"
#include "output.h"
#include "memorymanag.h"
#include "utilities.h"
#include "mat_gblas.h"


BoxSystem::BoxSystem()
{
  boxList = NULL;
}

BoxSystem::~BoxSystem()
{
  if(boxList)
    delete [] boxList;
}


/** Creates the box system.

    @param itemList list of items to create the box structure for.
    @param noOfItems their number.
    @param toplevelBoxSize
*/

int
BoxSystem::create_box_system(box_item_struct* itemList,
			     int noOfItems,
			     ergo_real toplevelBoxSize)
{
  // Allocate resultBoxList with just one item to begin with,
  // It will be expanded later as needed.
  int maxNoOfBoxes = 1;
  boxList = new box_struct_basic[maxNoOfBoxes];

  // create "mother box" containing all distrs.
  int currBoxIndex = 0;
  if(currBoxIndex >= maxNoOfBoxes)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system: (currBoxIndex >= maxNoOfBoxes)");
      return -1;
    }
  box_struct_basic* motherBox = &boxList[currBoxIndex];
  currBoxIndex++;
  
  // first get dimensions of box
  const ergo_real HUGE_NUMBER = 888888888;
  ergo_real xminList[3];
  ergo_real xmaxList[3];
  ergo_real xdiffList[3];
  for(int kk = 0; kk < 3; kk++)
    {
      xminList[kk] =  HUGE_NUMBER;
      xmaxList[kk] = -HUGE_NUMBER;
    }
  for(int i = 0; i < noOfItems; i++)
    {
      for(int kk = 0; kk < 3; kk++)
	{
	  ergo_real x = itemList[i].centerCoords[kk];
	  if(x < xminList[kk])
	    xminList[kk] = x;
	  if(x > xmaxList[kk])
	    xmaxList[kk] = x;
	}
    } // END FOR i
  int bestCoordIndex = 0;
  for(int kk = 0; kk < 3; kk++)
    {
      xdiffList[kk] = xmaxList[kk] - xminList[kk];
      if(xdiffList[kk] > xdiffList[bestCoordIndex])
	bestCoordIndex = kk;
    }
  for(int kk = 0; kk < 3; kk++) {
    if(xdiffList[kk] < 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system: (xdiffList[kk] < 0).");
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system: kk = %d, xdiffList[kk] = %9.4g < 0.\n", kk, xdiffList[kk]);
      return -1;
    }
  }

  ergo_real largestDiff = xdiffList[bestCoordIndex];

  // compute number of levels and size of mother box
  int numberOfLevels = 1;
  ergo_real width = toplevelBoxSize;
  while(width < largestDiff)
    {
      width *= 2;
      numberOfLevels++;
    }

  if(numberOfLevels >= MAX_NO_OF_BOX_LEVELS)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system: (numberOfLevels >= MAX_NO_OF_BOX_LEVELS)");
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system: numberOfLevels = %d, MAX_NO_OF_BOX_LEVELS = %d.", numberOfLevels, (int)MAX_NO_OF_BOX_LEVELS);
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system: toplevelBoxSize = %7.3g, largestDiff = %7.3g.", (double)toplevelBoxSize, (double)largestDiff);
      return -1;
    }

  motherBox->width = width;
  // Set center of mother box.
  for(int kk = 0; kk < 3; kk++)
    motherBox->centerCoords[kk] = xminList[kk] + motherBox->width / 2;
  motherBox->firstItemIndex = 0;
  motherBox->noOfItems = noOfItems;
  
  // OK, mother box done.
  // Now create boxes on other levels, one level at a time.
  
  this->levelList[0].noOfBoxes = 1;
  this->levelList[0].startIndexInBoxList = 0;

  int* itemIndexBucketList[2][2][2];
  for(int ix = 0; ix < 2; ix++)
    for(int iy = 0; iy < 2; iy++)
      for(int iz = 0; iz < 2; iz++)
	itemIndexBucketList[ix][iy][iz] = new int[noOfItems];
  int itemCounterList[2][2][2];
  // Allocate temporary itemList for use when reordering items.
  std::vector<box_item_struct> itemListTemp(noOfItems);

  for(int levelNumber = 1; levelNumber < numberOfLevels; levelNumber++)
    {
      levelList[levelNumber].startIndexInBoxList = currBoxIndex;
      int currLevelNoOfBoxes_1 = 0;
      // go through boxes of previous level, and create new boxes at this level if needed.
      // We go through boxes of previous level twice, the first time to check how many
      // new boxes must be allocated.
      int startIndex = levelList[levelNumber-1].startIndexInBoxList;
      for(int i = startIndex; i < startIndex + levelList[levelNumber-1].noOfBoxes; i++)
	{
	  // now resultBoxList[i] is the box whose children we are creating.
	  boxList[i].firstChildBoxIndex = currBoxIndex;
	  int noOfChildren = 0;
	  if(boxList[i].noOfItems == 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "ERROR: (boxList[i].noOfItems == 0)");
	      return -1;
	    }
	  // create 2*2*2 = 8 new boxes
	  for(int ix = 0; ix < 2; ix++)
	    for(int iy = 0; iy < 2; iy++)
	      for(int iz = 0; iz < 2; iz++)
		{
		  itemCounterList[ix][iy][iz] = 0;
		} // END FOR ix iy iz
	  // now go through the points of the parent box to decide which of the 8 children each point belongs to.
	  for(int j = boxList[i].firstItemIndex; j < boxList[i].firstItemIndex + boxList[i].noOfItems; j++)
	    {
	      int ix = 0;
	      int iy = 0;
	      int iz = 0;
	      if(itemList[j].centerCoords[0] > boxList[i].centerCoords[0])
		ix = 1;
	      if(itemList[j].centerCoords[1] > boxList[i].centerCoords[1])
		iy = 1;
	      if(itemList[j].centerCoords[2] > boxList[i].centerCoords[2])
		iz = 1;
	      itemCounterList[ix][iy][iz]++;
	    } // END FOR j
	  // OK, all items counted, counters updated to indicate how many items 
	  // belong to each of the 8 candidate child boxes.
	  // Now count new boxes.
	  for(int ix = 0; ix < 2; ix++)
	    for(int iy = 0; iy < 2; iy++)
	      for(int iz = 0; iz < 2; iz++)
		{
		  if(itemCounterList[ix][iy][iz] > 0)
		    {
		      currLevelNoOfBoxes_1++;
		      noOfChildren++;
		    } // END IF (itemCounterList[ix][iy][iz] > 0)
		} // END FOR ix iy iz
	} // END FOR i
      
      // OK, now we know how many new boxes are needed.
      int maxNoOfBoxesNew = maxNoOfBoxes + currLevelNoOfBoxes_1;
      box_struct_basic* boxListNew = new box_struct_basic[maxNoOfBoxesNew];
      // copy previous contents of resultBoxList
      memcpy(boxListNew, boxList, maxNoOfBoxes*sizeof(box_struct_basic));
      // free old list, and set pointer to new list instead.
      delete [] boxList;
      boxList = boxListNew;
      maxNoOfBoxes = maxNoOfBoxesNew;

      // Now go through level again, this time creating new boxes and reordering items.
      int currLevelNoOfBoxes_2 = 0;
      for(int i = startIndex; i < startIndex + levelList[levelNumber-1].noOfBoxes; i++)
	{
	  // now resultBoxList[i] is the box whose children we are creating.
	  boxList[i].firstChildBoxIndex = currBoxIndex;
	  int noOfChildren = 0;
	  if(boxList[i].noOfItems == 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "ERROR: (resultBoxList[i].noOfItems == 0)");
	      return -1;
	    }
	  // create 2*2*2 = 8 new boxes
	  box_struct_basic tempBoxList[2][2][2];
	  for(int ix = 0; ix < 2; ix++)
	    for(int iy = 0; iy < 2; iy++)
	      for(int iz = 0; iz < 2; iz++)
		{
		  ergo_real newWidth = boxList[i].width / 2;
		  ergo_real x = boxList[i].centerCoords[0] + (ix*2-1)*0.5*newWidth;
		  ergo_real y = boxList[i].centerCoords[1] + (iy*2-1)*0.5*newWidth;
		  ergo_real z = boxList[i].centerCoords[2] + (iz*2-1)*0.5*newWidth;
		  tempBoxList[ix][iy][iz].centerCoords[0] = x;
		  tempBoxList[ix][iy][iz].centerCoords[1] = y;
		  tempBoxList[ix][iy][iz].centerCoords[2] = z;
		  tempBoxList[ix][iy][iz].width = newWidth;
		  itemCounterList[ix][iy][iz] = 0;
		} // END FOR ix iy iz
	  // now go through the points of the parent box to decide which of the 8 children each point belongs to.
	  for(int j = boxList[i].firstItemIndex; j < boxList[i].firstItemIndex + boxList[i].noOfItems; j++)
	    {
	      int ix = 0;
	      int iy = 0;
	      int iz = 0;
	      if(itemList[j].centerCoords[0] > boxList[i].centerCoords[0])
		ix = 1;
	      if(itemList[j].centerCoords[1] > boxList[i].centerCoords[1])
		iy = 1;
	      if(itemList[j].centerCoords[2] > boxList[i].centerCoords[2])
		iz = 1;
	      itemIndexBucketList[ix][iy][iz][itemCounterList[ix][iy][iz]] = j;
	      itemCounterList[ix][iy][iz]++;
	    } // END FOR j
	  // OK, all items copied to itemBucketList.
	  // Now add new boxes, and order the items accordingly.
	  int currStartIndex = boxList[i].firstItemIndex;
	  for(int ix = 0; ix < 2; ix++)
	    for(int iy = 0; iy < 2; iy++)
	      for(int iz = 0; iz < 2; iz++)
		{
		  if(itemCounterList[ix][iy][iz] > 0)
		    {
		      if(currBoxIndex >= maxNoOfBoxes)
			{
			  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system: (currBoxIndex >= maxNoOfBoxes)");
			  return -1;
			}
		      boxList[currBoxIndex] = tempBoxList[ix][iy][iz];
		      boxList[currBoxIndex].firstItemIndex = currStartIndex;
		      boxList[currBoxIndex].noOfItems = itemCounterList[ix][iy][iz];
		      for(int k = 0; k < itemCounterList[ix][iy][iz]; k++)
			itemListTemp[currStartIndex+k] = itemList[itemIndexBucketList[ix][iy][iz][k]];
		      currStartIndex += itemCounterList[ix][iy][iz];
		      currBoxIndex++;
		      currLevelNoOfBoxes_2++;
		      noOfChildren++;
		    } // END IF (itemCounterList[ix][iy][iz] > 0)
		} // END FOR ix iy iz
	  boxList[i].noOfChildBoxes = noOfChildren;
	} // END FOR i
      if(currLevelNoOfBoxes_2 != currLevelNoOfBoxes_1)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system: (currLevelNoOfBoxes_2 != currLevelNoOfBoxes_1)");
	  return -1;
	}
      levelList[levelNumber].noOfBoxes = currLevelNoOfBoxes_2;
      memcpy(itemList, &itemListTemp[0], noOfItems*sizeof(box_item_struct));
    } // END FOR levelNumber
  
  // OK, boxes created.

  totNoOfBoxes = currBoxIndex;
  noOfLevels = numberOfLevels;

  for(int ix = 0; ix < 2; ix++)
    for(int iy = 0; iy < 2; iy++)
      for(int iz = 0; iz < 2; iz++)
	delete [] itemIndexBucketList[ix][iy][iz];

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "create_box_system end OK, toplevelBoxSize: %5.1f, #levels: %2i, #boxes at top level: %6i", 
	    (double)toplevelBoxSize, numberOfLevels, levelList[noOfLevels-1].noOfBoxes);
  
  return 0;
}



ergo_real
get_min_distance_from_point_to_box(const ergo_real* boxCenterCoords,
				   ergo_real halfwidth,
				   const ergo_real* point)
{
  ergo_real dxList[3];
  for(int k = 0; k < 3; k++)
    {
      ergo_real dx = template_blas_fabs(boxCenterCoords[k] - point[k]);
      if(dx > halfwidth)
	dxList[k] = dx - halfwidth;
      else
	dxList[k] = 0;
    }
  ergo_real sum = 0;
  for(int k = 0; k < 3; k++)
    sum += dxList[k] * dxList[k];
  return template_blas_sqrt(sum);
}

int BoxSystem::get_items_near_point_recursive(const box_item_struct* itemList,
					      const ergo_real* coords, 
					      ergo_real distance, 
					      int* resultOrgIndexList,
					      int level,
					      int boxIndex) const
{
  const box_struct_basic* box = &boxList[boxIndex];
  // Check if this entire box is so far away that it can be skipped.
  ergo_real min_distance_from_point_to_box = get_min_distance_from_point_to_box(box->centerCoords, box->width/2, coords);
  if(min_distance_from_point_to_box > distance)
    return 0;
  // No, we could not skip. Take care of box contents.
  if(level == noOfLevels-1)
    {
      // We are at top level. 
      // Go through all points in box and add the relevant ones to result list.
      int count = 0;
      for(int i = 0; i < box->noOfItems; i++)
	{
	  const box_item_struct* currItem = &itemList[box->firstItemIndex+i];
	  ergo_real sum = 0;
	  for(int coordNo = 0; coordNo < 3; coordNo++)
	    {
	      ergo_real d = coords[coordNo] - currItem->centerCoords[coordNo];
	      sum += d*d;
	    }
	  ergo_real currDist = template_blas_sqrt(sum);
	  if(currDist < distance)
	    {
	      // Add to result.
	      resultOrgIndexList[count] = currItem->originalIndex;
	      count++;
	    }
	} // END FOR i
      return count;
    }
  // Not top level. Go to next level.
  int count = 0;
  for(int i = 0; i < box->noOfChildBoxes; i++)
    {
      int childBoxIndex = box->firstChildBoxIndex + i;
      int nNew = get_items_near_point_recursive(itemList,
						coords, 
						distance, 
						&resultOrgIndexList[count],
						level+1,
						childBoxIndex);
      count += nNew;
    }
  return count;
}



static int
compare_ints(const void* p1, const void* p2)
{
  int i1 = *((int*)p1);
  int i2 = *((int*)p2);
  if(i1 > i2)
    return 1;
  if(i1 < i2)
    return -1;
  return 0;
}



/** Goes through existning box system to find all items within 
    specified distance from given reference point.

    @param itemList the list of items for which the box system was created.
    @param coords list of 3 coordinates for reference point.
    @param distance the distance to find items within.
    @param resultOrgIndexList preallocated list of resulting org indexes.
*/

int BoxSystem::get_items_near_point(const box_item_struct* itemList,
				    const ergo_real* coords, 
				    ergo_real distance, 
				    int* resultOrgIndexList) const
{
  if(!boxList)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in BoxSystem::get_items_near_point: (!boxList). Must create box system first!");
      return -1;
    }
  int noOfItems = get_items_near_point_recursive(itemList, coords, distance, resultOrgIndexList, 0, levelList[0].startIndexInBoxList);
  // sort resultOrgIndexList
  qsort(resultOrgIndexList, noOfItems, sizeof(int), compare_ints);
  return noOfItems;
}


