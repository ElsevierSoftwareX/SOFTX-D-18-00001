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

/** @file xyz_file_parser.cc

    @brief Functionality for parsing a file in xyz file format,
    storing the result as a molecule object.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <string.h>
#include <cassert>
#include "xyz_file_parser.h"
#include "output.h"
#include "memorymanag.h"
#include "units.h"
#include "atom_labels.h"
#include "utilities.h"


int 
readMoleculeFileInXyzFormat(Molecule& result, 
			    const char* fileName, 
			    int netCharge,
			    bool expectPlainCharges)
{
  FILE* f = fopen(fileName, "rt");
  if(f == NULL) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error opening file '%s'", fileName);
    return -1;
  }

  long int  fileSize = get_file_size(fileName);
  if(fileSize < 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error getting file size for file '%s'", fileName);
    return -1;    
  }
  size_t bufSize = fileSize + 10000;
  std::vector<char> buf(bufSize);
  memset(&buf[0], 0, bufSize);
  
  int nBytes = (int)fread(&buf[0], sizeof(char), bufSize, f);
  fclose(f);
  if(nBytes <= 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error reading file '%s'", fileName);
    return -1;
  }

  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "read_molecule_file_in_xyz_format: File '%s' read OK, nBytes = %i", fileName, nBytes);

  const char* p = &buf[0];

  // Get number of atoms from first line
  // Skip blanks
  while(*p == ' ' || *p == '\t')
    p++;
  int nAtoms = atoi(p);
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "read_molecule_file_in_xyz_format: nAtoms = %6i", nAtoms);

  // Skip rest of line
  while(*p != '\n' && *p != '\0')
    p++;
  p++;

  // Skip one line
  while(*p != '\n' && *p != '\0')
    p++;
  p++;

  // Now p should point to the line with the first atom
  
  int atomIndex;
  for(atomIndex = 0; atomIndex < nAtoms; atomIndex++)
    {
      // Skip blanks
      while(*p == ' ' || *p == '\t')
	p++;
      // Now p should point to atom label.
      char labelString[88];
      const char* q = p;
      // Move q forwarn until blank
      while(*q != ' ' && *q != '\t' && *q != '\n' && *q != '\0')
	q++;
      int labelLen = q - p;
      if(labelLen > 22)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in read_molecule_file_in_xyz_format: label too long.");
	  return -1;
	}
      memcpy(labelString, p, labelLen);
      labelString[labelLen] = '\0';
      p = q;

      // Now we have the label string null-terminated in labelString.
      // We want to handle the case when the labels have some numbers on each label, stuff like "N17" and "Fe31".
      // Do this by simply replacing any digits in labelString with null-characters.
      for(int i = 0; i < labelLen; i++) {
        if(labelString[i] >= '0' && labelString[i] <= '9')
          labelString[i] = '\0';
      }

      // Now p should point to first blank after label.
      if(*p != ' ' && *p != '\t')
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in read_molecule_file_in_xyz_format: blank not found after label.");
	  return -1;
	}
      // Skip blanks
      while(*p == ' ' || *p == '\t')
	p++;
      // Now p should point to x coordinate.
      ergo_real coord_x = atof(p);
      // Skip coordinate
      while((*p >= '0' && *p <= '9') || *p == '.' || *p == '-')
	p++;
      // Skip blanks
      while(*p == ' ' || *p == '\t')
	p++;
      // Now p should point to y coordinate.
      ergo_real coord_y = atof(p);
      // Skip coordinate
      while((*p >= '0' && *p <= '9') || *p == '.' || *p == '-')
	p++;
      // Skip blanks
      while(*p == ' ' || *p == '\t')
	p++;
      // Now p should point to z coordinate.
      ergo_real coord_z = atof(p);
      // Skip coordinate
      while((*p >= '0' && *p <= '9') || *p == '.' || *p == '-')
	p++;
      // Skip blanks
      while(*p == ' ' || *p == '\t')
	p++;
      // Now p should point to newline character.
      if(*p != '\n')
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in read_molecule_file_in_xyz_format: newline not found after coordinates.");
	  printf("p points to '%s'\n", p);
	  return -1;
	}
      // Skip to next line
      p++;
      
      // OK, now we have labelString and coords.

      // Depending on the expectPlainCharges parameter, we interpret
      // the labelString as a charge or as an atom label.
      ergo_real atomCharge = 0;
      if(expectPlainCharges) {
	atomCharge = atof(labelString);
      }
      else {
	// Get charge corresponding to labelString.
	int chargeInt = get_charge_int_from_atom_label(labelString);
	if(chargeInt <= 0) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, 
		    "error in read_molecule_file_in_xyz_format: label '%s' not recognized as an atom type.", 
		    labelString);
	  return -1;
	}
	atomCharge = chargeInt;
      }

      result.addAtom(atomCharge,
		     coord_x * UNIT_one_Angstrom,
		     coord_y * UNIT_one_Angstrom,
		     coord_z * UNIT_one_Angstrom);
      
    } // END FOR atomIndex

  // accept only blank space and newlines after last atom
  while(*p == ' ' || *p == '\t' || *p == '\n')
    p++;
  // Now p should point to end of buffer
  if(*p != '\0')
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, 
		"error in read_molecule_file_in_xyz_format: garbage found after last atom.");
      return -1;
    }
  
  assert(result.getNoOfAtoms() == nAtoms);
  result.setNetCharge(netCharge);
  
  return 0;
}

