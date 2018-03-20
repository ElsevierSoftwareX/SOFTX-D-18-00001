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

/** @file box_system.h

    @brief BoxSystem class representing a hierarchical data structure
    of boxes in 3D space (an oct-tree).

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef BOX_SYSTEM_HEADER
#define BOX_SYSTEM_HEADER


#include "realtype.h"


typedef struct
{
  ergo_real centerCoords[3];
  int originalIndex;
} box_item_struct;


typedef struct
{
  ergo_real centerCoords[3];
  ergo_real width;
  int noOfItems;
  int firstItemIndex;
  int noOfChildBoxes;
  int firstChildBoxIndex;
} box_struct_basic;

typedef struct
{
  int noOfBoxes;
  int startIndexInBoxList;
} box_level_struct;

#define MAX_NO_OF_BOX_LEVELS 30

class BoxSystem
{
 public:
  int totNoOfBoxes;
  int noOfLevels;
  box_level_struct levelList[MAX_NO_OF_BOX_LEVELS];
  box_struct_basic* boxList;
  BoxSystem();
  ~BoxSystem();
  int create_box_system(box_item_struct* itemList,
			int noOfItems,
			ergo_real toplevelBoxSize);
  int get_items_near_point(const box_item_struct* itemList,
			   const ergo_real* coords, 
			   ergo_real distance, 
			   int* resultOrgIndexList) const;
 private:
  int get_items_near_point_recursive(const box_item_struct* itemList,
				     const ergo_real* coords, 
				     ergo_real distance, 
				     int* resultOrgIndexList,
				     int level,
				     int boxIndex) const;
};

ergo_real
get_min_distance_from_point_to_box(const ergo_real* boxCenterCoords,
				   ergo_real halfwidth,
				   const ergo_real* point);



#endif
