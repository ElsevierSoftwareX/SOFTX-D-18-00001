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

/** @file atom_labels.cc

    @brief Provides a way to map atom labels to their charges. The
    main procedure provided by this file is
    get_charge_int_from_atom_label().

    @author: Elias Rudberg <em>responsible</em>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <string>

#include "atom_labels.h"

#define kMaxNoOfAtomTypes 200


static void setupLabelList(std::string labelList []) {
  int i;
  for(i = 0; i < kMaxNoOfAtomTypes; i++)
    labelList[i] = "";
  labelList[1]  = "H";
  labelList[2]  = "He";
  labelList[3]  = "Li";
  labelList[4]  = "Be";
  labelList[5]  = "B";
  labelList[6]  = "C";
  labelList[7]  = "N";
  labelList[8]  = "O";
  labelList[9]  = "F";
  labelList[10] = "Ne";
  labelList[11] = "Na";
  labelList[12] = "Mg";
  labelList[13] = "Al";
  labelList[14] = "Si";
  labelList[15] = "P";
  labelList[16] = "S";
  labelList[17] = "Cl";
  labelList[18] = "Ar";
  labelList[19] = "K";
  labelList[20] = "Ca";
  labelList[21] = "Sc";
  labelList[22] = "Ti";
  labelList[23] = "V";
  labelList[24] = "Cr";
  labelList[25] = "Mn";
  labelList[26] = "Fe";
  labelList[27] = "Co";
  labelList[28] = "Ni";
  labelList[29] = "Cu";
  labelList[30] = "Zn";
  labelList[31] = "Ga";
  labelList[32] = "Ge";
  labelList[33] = "As";
  labelList[34] = "Se";
  labelList[35] = "Br";
  labelList[36] = "Kr";
  labelList[37] = "Rb";
  labelList[38] = "Sr";
  labelList[39] = "Y";
  labelList[40] = "Zr";
  labelList[41] = "Nb";
  labelList[42] = "Mo";
  labelList[43] = "Tc";
  labelList[44] = "Ru";
  labelList[45] = "Rh";
  labelList[46] = "Pd";
  labelList[47] = "Ag";
  labelList[48] = "Cd";
  labelList[49] = "In";
  labelList[50] = "Sn";
  labelList[51] = "Sb";
  labelList[52] = "Te";
  labelList[53] = "I";
  labelList[54] = "Xe";  
}


int get_charge_int_from_atom_label(const char* atomLabel)
{
  std::string labelList[kMaxNoOfAtomTypes];
  int i;

  if(atomLabel == NULL)
    return -1;

  if(atomLabel[0] == '\0')
    return -1;
  
  /* handle all-digit atomLabels */
  for(i=0; atomLabel[i] && isdigit(atomLabel[i]); i++)
    ;

  if( atomLabel[i] == '\0') /* all digit atom label */
    return atoi(atomLabel);

  setupLabelList(labelList);

  for(i = 1; i < kMaxNoOfAtomTypes; i++)
    {
      if(strcasecmp(atomLabel, labelList[i].c_str()) == 0)
	return i;
    }

  return -1;
}


int get_atom_label_from_charge_int(int charge, char* atomLabelString, size_t bufferSize)
{
  std::string labelList[kMaxNoOfAtomTypes];
  memset(atomLabelString, 0, bufferSize);
  if(charge <= 0 || charge >= kMaxNoOfAtomTypes)
    return -1;
  setupLabelList(labelList);
  if(labelList[charge].length() >= bufferSize)
    return -1;
  strcpy(atomLabelString, labelList[charge].c_str());
  return 0;
}

