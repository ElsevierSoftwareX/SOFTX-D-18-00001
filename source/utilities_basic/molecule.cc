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

/** @file molecule.cc

    @brief Class representing a molecule as a set of atoms with
    assiciated coordinates and charges of the atomic nuclei.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <string.h>
#include <cassert>
#include "molecule.h"
#include "xyz_file_parser.h"
#include "output.h"
#include "memorymanag.h"
#include "units.h"
#include "utilities.h"


static ergo_real
get_distance_between_atoms(const Atom& atomA, const Atom& atomB)
{
  int i;
  ergo_real sum = 0;
  for(i = 0; i < 3; i++)
    {
      ergo_real dx = atomB.coords[i] - atomA.coords[i];
      sum += dx*dx;
    }
  return template_blas_sqrt(sum);
}

void
Molecule::getExtremeInternuclearDistancesQuadratic(ergo_real & minDist, ergo_real & maxDist) const
{
  minDist = maxDist = 0; // This matters if there is only one atom.
  bool firstTime = true;
  for(int A = 0; A < noOfAtoms; A++)
    for(int B = A+1; B < noOfAtoms; B++) {
      const Atom& atomA = atoms[A];
      const Atom& atomB = atoms[B];
      ergo_real RAB = get_distance_between_atoms(atomA, atomB);
      if(firstTime) {
	minDist = maxDist = RAB;
	firstTime = false;
      }
      else {
	minDist = RAB < minDist ? RAB : minDist;
	maxDist = RAB > maxDist ? RAB : maxDist;
      }
    }
}

ergo_real
Molecule::getNuclearRepulsionEnergyQuadratic() const
{
  int A, B;
  ergo_real sum;
  int N = noOfAtoms;
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, 
	    "get_nuclear_repulsion_energy, N = %i", N);
  sum = 0;
  for(A = 0; A < N; A++)
    for(B = A+1; B < N; B++)
      {
	const Atom& atomA = atoms[A];
	const Atom& atomB = atoms[B];
	ergo_real ZA = atomA.charge;
	ergo_real ZB = atomB.charge;
	ergo_real RAB = get_distance_between_atoms(atomA, atomB);
	sum += ZA * ZB / RAB;
      }
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, 
	    "get_nuclear_repulsion_energy returning energy %33.22f", (double)sum);
  return sum;
}

void
Molecule::getNuclearRepulsionEnergyGradientContribQuadratic(ergo_real* resultGradient) const 
{
  int N = noOfAtoms;
  int A, B;
  for(A = 0; A < N; A++)
    for(B = A+1; B < N; B++) {
      const Atom& atomA = atoms[A];
      const Atom& atomB = atoms[B];
      ergo_real ZA = atomA.charge;
      ergo_real ZB = atomB.charge;
      ergo_real dx = atomB.coords[0] - atomA.coords[0];
      ergo_real dy = atomB.coords[1] - atomA.coords[1];
      ergo_real dz = atomB.coords[2] - atomA.coords[2];
      ergo_real r = template_blas_sqrt(dx*dx + dy*dy + dz*dz);
      ergo_real derivative_x = ZA * ZB * (-0.5) * 2.0 * dx / (r*r*r);
      ergo_real derivative_y = ZA * ZB * (-0.5) * 2.0 * dy / (r*r*r);
      ergo_real derivative_z = ZA * ZB * (-0.5) * 2.0 * dz / (r*r*r);
      resultGradient[A*3+0] += -1 * derivative_x;
      resultGradient[A*3+1] += -1 * derivative_y;
      resultGradient[A*3+2] += -1 * derivative_z;
      resultGradient[B*3+0] +=  1 * derivative_x;
      resultGradient[B*3+1] +=  1 * derivative_y;
      resultGradient[B*3+2] +=  1 * derivative_z;
    }    
}

ergo_real
Molecule::getNuclearElectricFieldEnergy(const Vector3D& electricField) const
{
  int N = noOfAtoms;
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, 
	    "get_nuclear_electric_field_energy, N = %i", N);
  ergo_real sum = 0;
  int A;
  for(A = 0; A < N; A++)
    {
      const Atom* atom = &atoms[A];
      sum += atom->charge * atom->coords[0] * electricField.v[0];
      sum += atom->charge * atom->coords[1] * electricField.v[1];
      sum += atom->charge * atom->coords[2] * electricField.v[2];
    } // END FOR A
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, 
	    "get_nuclear_electric_field_energy returning energy %33.22f", (double)sum);
  return sum;
}

int 
Molecule::getNumberOfElectrons() const
{
  int i;
  ergo_real noOfElectrons;
  ergo_real sum = 0;
  for(i = 0; i < noOfAtoms; i++)
    sum += atoms[i].charge;
  noOfElectrons = sum - netCharge;
  return (int)noOfElectrons;
}


static int 
readMoleculeFileInMolFormat(Molecule* result, const char* fileName, 
                            int netCharge, char **basisfilename)
{
  int nAtoms, nAtomsCurrType;
  int i, j, ii, lineLen, noOfAtomTypes, atomTypeNo, can_read_basset;
  char* p;
  ergo_real unitFactor, atomCharge;
  long int  fileSize = get_file_size(fileName);
  if(fileSize < 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error getting file size for file '%s'", fileName);
    return -1;    
  }
  FILE* f = fopen(fileName, "rt");
  if(f == NULL)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error opening file '%s'", fileName);
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

  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "File '%s' read OK, nBytes = %i", fileName, nBytes);

  nAtoms = 0;

  p = &buf[0];

  // skip 3 or 4 lines
  if(strncmp(p, "BASIS", 5) == 0) { i = 3; can_read_basset = 1; }
  else if(strncmp(p, "ATOMDF", 6) == 0 ||
          strncmp(p, "ATOMBASIS", 9) == 0) { i = 2; can_read_basset = 0; }
  else 
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "Unsupported molecule input file type");
      return -1;
    }
  while((*p != '\n') && (*p != '\0')) p++; /* skip first line */
  p++;
  // if basisfilename is NULL we skip getting basis set string.
  if(basisfilename)
    {
      if(can_read_basset && !*basisfilename) {
	char *eol = strchr(p, '\n');
	if(eol) {
	  while (--eol>p && *eol == ' ')
	    ;
	  if (eol>p) {
	    size_t len = eol-p+1;
	    *basisfilename = (char*)ergo_malloc(len+1);
	    strncpy(*basisfilename, p, len);
	    (*basisfilename)[len] = '\0';
	  }
	}
      }
    }
  while(i--) {
    while((*p != '\n') && (*p != '\0')) p++;
    p++;
  }
  if(*p == '\0') return -2;
  
  // the molecule input can follow in a fixed or a free format.
  unitFactor = 1.0;
  if( sscanf(p, "%d", &noOfAtomTypes) == 1)
    { 
      // check if Angstrom flag is set
      lineLen = 0;
      for(ii = 0; ii <= 19; ii++)
        {
          if((p[ii] == '\n') || (p[ii] == '\0'))
            break;
          lineLen++;
        }
      if(lineLen > 19)
        {
          if(p[19] != ' ')
            unitFactor = UNIT_one_Angstrom;
        }
      // skip till end of line
      while(*p && *p != '\n') p++;
    }
  else
    { /* "free" format */
      do {
        if(strncasecmp(p, "Atomtypes=", 10) == 0)
          sscanf(p+10, "%d", &noOfAtomTypes);
        else if(strncmp(p, "Angstrom", 8) == 0)
          unitFactor = UNIT_one_Angstrom;
        while(*p && *p != '\n' && *p != ' ') p++;
        if(*p != ' ') break;
        p++;
      } while(1);
    }
  if(*p == '\n') p++; 
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN,
	    "noOfAtomTypes = %i dist. unit=%f au", noOfAtomTypes, unitFactor);
  if(noOfAtomTypes <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "Error: (noOfAtomTypes <= 0)");
      return -3;
    }
  
  for(atomTypeNo = 0; atomTypeNo < noOfAtomTypes; atomTypeNo++)
    {
      // now p should point to beginning of new atom type
      float charge;
      if(*p == '\0')
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in readMoleculeFile for atom type %i", 
		    atomTypeNo+1);
	  return -4;
	}

      do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "atomTypeNo = %i", atomTypeNo);

      // skip blanks
      while(*p == ' ') p++;
      if(*p == '\0') return -5;
      // if the first character on the line is a digit, we assume fixed format
      if(*p >= '0' && *p <= '9')
	{
	  charge = (float)atof(p);
	  while(*p != ' ' && *p != '\n' && *p != '\0') p++;
	  // skip blanks
	  while(*p == ' ') p++;
	  // now we expect another digit
	  if(*p < '0' || *p > '9')
	    return -6;
	  nAtomsCurrType = atoi(p);
	  while(*p != ' ' && *p != '\n' && *p != '\0') p++;
	  // skip blanks
	  while(*p == ' ') p++;
	  // now we expect end of line
	  if(*p != '\n')
	    return -7;
	  p++;
	}
      else
        { /* "free" format */
	  charge = -1;
	  nAtomsCurrType = -1;
          do {
            if(strncmp(p, "Charge=", 7) == 0)
              sscanf(p+7, "%f", &charge);
            else if(strncmp(p, "Atoms=", 6) == 0)
              sscanf(p+6, "%d", &nAtomsCurrType);
            while(*p && *p != '\n' && *p != ' ') p++;
            if(*p != ' ') break;
            p++;
          } while(1);
	  if(charge < 0) {
	    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, 
		      "Invalid charge = %5.2f", (double)charge);
            return -8;
	  }
	  if(nAtomsCurrType < 0) {
	    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, 
		      "Invalid number of atoms = %i", nAtomsCurrType);
		       return -8;
	  }
        }

      do_output(LOG_CAT_INFO, LOG_AREA_MAIN, 
		"charge = %5.2f, nAtomsCurrType = %i", (double)charge, nAtomsCurrType);
      
      atomCharge = charge;
      if(atomCharge <= 0)
        return -9;
          
      for(i = 0; i < nAtomsCurrType; i++)
	{
	  // skip atom label
	  while(*p == ' ') p++;
	  while((*p != ' ') && (*p != '\0')) p++;
	  if(*p == '\0') return -10;
	  while(*p == ' ') p++;
	  // read coordinates
	  ergo_real coords[3];
	  for(j = 0; j < 3; j++)
	    {
	      // skip blanks
	      while(*p == ' ') p++;
	      if(*p == '\0') return -11;
	      coords[j] = atof(p) * unitFactor;
	      while((*p != ' ') && (*p != '\n') && (*p != '\0')) p++;
	      if(*p == '\0') return -12;
	    }

	  // skip rest of line
	  while((*p != '\n') && (*p != '\0'))
	    p++;
	  p++;

	  result->addAtom(atomCharge, coords[0], coords[1], coords[2]);
	  nAtoms++;
	} // END FOR i
    } // END FOR atomTypeNo

  // done. now there should be nothing more in the file.
  while(*p != '\0')
    {
      if(*p != ' ' && *p != '\n')
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in readMoleculeFile: "
		    "non-blank character found after last atom.");
	  return -14;
	}
      p++;
    }

#if 0
  for(i = 0; i < nAtoms; i++)
    {
      printf("%5.2f %22.11f %22.11f %22.11f\n", 
	     result->atoms[i].charge,
	     result->atoms[i].coords[0],
	     result->atoms[i].coords[1],
	     result->atoms[i].coords[2]);
    }
#endif

  assert(result->getNoOfAtoms() == nAtoms);
  result->setNetCharge(netCharge);
  
  return 0;
}



int 
Molecule::setFromMoleculeFile(const char* fileName, 
                              int netCharge, char **basisfilename)
{
  // Check filename extension
  int len = strlen(fileName);
  if(len < 5)
    {
      // Too short to have a 3-letter extension after dot.
      // Assume mol format.
      return readMoleculeFileInMolFormat(this, fileName, netCharge,
                                         basisfilename);
    }
  const char* extensionPtr = &fileName[len-4];
  if(strcasecmp(extensionPtr, ".xyz") == 0)
    return readMoleculeFileInXyzFormat(*this, fileName, netCharge, false);

  // Not xyz, assume mol format.
  return readMoleculeFileInMolFormat(this, fileName, netCharge,
                                     basisfilename);
}


