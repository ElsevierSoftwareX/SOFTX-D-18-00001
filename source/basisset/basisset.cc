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

/** @file basisset.cc

    \brief Code for representing basis set information for Gaussian
    basis sets, and for parsing a text file specifying such a
    basisset.

    @author: Elias Rudberg <em>responsible</em>.
*/

/* -*-mode:c; indent-tabs-mode: nil -*- */
/* basisset.c: provides read_basisset_file() which creates a
   basisset_struct from a data contained in a file */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <stdexcept>

#include "basisset.h"
#include "output.h"
#include "memorymanag.h"


basisset_info::basisset_info() {
  atoms.resize(MAX_NO_OF_ATOM_TYPES);
  clear();
}

void basisset_info::clear() {
  memset(&atoms[0], 0x00, MAX_NO_OF_ATOM_TYPES*sizeof(basisset_atom_struct));
}

void basisset_info::write_to_buffer ( char * dataBuffer, size_t const bufferSize ) const {
  if(bufferSize < get_size())
    throw std::runtime_error("Error in basisset_info::write_to_buffer: bufferSize too small.");
  memcpy(dataBuffer, &atoms[0], MAX_NO_OF_ATOM_TYPES*sizeof(basisset_atom_struct));
}

size_t basisset_info::get_size() const {
  return MAX_NO_OF_ATOM_TYPES*sizeof(basisset_atom_struct);
}

void basisset_info::assign_from_buffer ( char const * dataBuffer, size_t const bufferSize) {
  if(bufferSize != get_size())
    throw std::runtime_error("Error in basisset_info::assign_from_buffer: wrong bufferSize.");
  memcpy(&atoms[0], dataBuffer, MAX_NO_OF_ATOM_TYPES*sizeof(basisset_atom_struct));
}


static void
remove_zeros(basisset_atom_struct* currAtom,
             int shellBaseIndex, int noOfShellsCurrBatch) {
  /*  Remove zero coefficients. */
  for(int i = 0; i < noOfShellsCurrBatch; i++) {
    int count = 0;
    ergo_real *coeffList = currAtom->shells[shellBaseIndex+i].coeffList;
    ergo_real *expList   = currAtom->shells[shellBaseIndex+i].exponentList;
    for(int j = 0; j < currAtom->shells[shellBaseIndex+i].contrCount; j++) {
      ergo_real currCoeff    = coeffList[j];
      ergo_real currExponent = expList[j];
      if(currCoeff != 0) {
	coeffList[count] = currCoeff;
	expList[count] = currExponent;
	count++;
      }
    } /*  END FOR j */
    currAtom->shells[shellBaseIndex+i].contrCount = count;
  } /*  END FOR i         */
}

/* read_basisset_file: reads a basis set from fileName. The basis set
   exponents and contraction coefficients are placed in result.  The
   reading procedure is bit convoluted because the basis set file
   follows the Fortran syntax, with wrapping and skipping empty
   elements. We basically need to emulate fortran `format'
   statement. There is one "but" though: AhlrichsDenFit does not
   follow the format syntax so it will/may be misread by eg. dalton. What
   a mess...
   
   The parser is implemented as a state machine. It still cannot read
   ANO-type basis sets...
*/
int 
read_basisset_file(basisset_info & result, 
		   const char* fileName,
                   int dirc, 
		   const char *dirv[],
                   int print_raw)
{
  enum { END_PARSING, ATOM_EXPECTED, SHELL_EXPECTED,
         SHELL_OR_ATOM_EXPECTED, CONTRACTION_BLOCK_EXPECTED } state;
  int uncontracted = 0;
  char line[512];
  basisset_atom_struct* currAtom = NULL;
  int spdf = -1;
  int shellBaseIndex = -1;
  int expNo = -1;
  FILE* f = NULL;

  if(!fileName) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: fileName == NULL.");
    return -1;
  }
  if(fileName[0] == '/')
    f = fopen(fileName, "rt");
  else {
    for(int i = 0; i < dirc; i++) {
      char ffname[256];
      int len = strlen(dirv[i]);
      strncpy(ffname, dirv[i], sizeof(ffname));
      if(len>0 && ffname[len-1] != '/')
        strncat(ffname, "/", sizeof(ffname)-1-len);
      strncat(ffname, fileName, sizeof(ffname)-2-len);
      do_output(LOG_CAT_WARNING, LOG_AREA_INTEGRALS, 
                "Trying basis set file '%s'...", ffname);
      if( (f=fopen(ffname, "rt")) != NULL)
        break;
    }
  }
      
  if(f == NULL)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error opening file '%s'", fileName);
      return -1;
    }
  
  /*  now create basis set by reading buf2 */
  result.clear();
  int noOfAtomTypes = 0;
  state = ATOM_EXPECTED;
  int lineNo = 0;
  int lineConsumed = 1;
  ergo_real currExponent = 0;
  /* start global parsing loop */
  do {
    int dummy;
    if(lineConsumed) {
      for(unsigned int cc = 0; cc < sizeof(line); cc++)
	line[cc] = '\0';
      if(fgets(line, sizeof(line), f) == NULL) {
        state = END_PARSING;
        break;
      }
      // This only works if sizeof(line) is large enough. Here we check that sizeof(line) is really large enough.
      if(strlen(line) >= sizeof(line)-1) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "Error in read_basisset_file: (strlen(line) >= sizeof(line)-1). "
		  "This indicates that the basis set file '%s' included a line that was too long, longer than %d characters.", fileName, (int)(sizeof(line)-1));
	return -1;
      }
      lineConsumed = 0;
      lineNo++;
    }

    for(int cc = strlen(line)-1; cc>=0 && isspace(line[cc]); cc--)
      line[cc] = '\0';
        
    if(line[0] == '$' || line[0] == '!' || line[0] == '#'||
       line[0] == '*' || line[0] == '\0'|| line[0] == '\n') {
      lineConsumed = 1; /* skip the comment and move on */
      continue;
    }
    switch(state) {
    case ATOM_EXPECTED:
      if(line[0] == 'a' || line[0] == 'A') {
        noOfAtomTypes++;
        int atomType = atoi(line+1);
        if(print_raw) 
          do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "Basis set for atom of Z=%d", atomType);
        state = SHELL_EXPECTED;
        if(atomType <= 0)
          {
            do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: (atomType <= 0) "
                      " in line %d %s", lineNo, fileName);
            return -1;
          }
        if(atomType >= MAX_NO_OF_ATOM_TYPES)
          {
            do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: "
                      "(atomType >= MAX_NO_OF_ATOM_TYPES) in line %d %s",
                      lineNo, fileName);
            return -1;
          }
        currAtom = &result.atoms[atomType];
        spdf = 0;
        shellBaseIndex = 0;
      }
      lineConsumed = 1;
      break;

    case SHELL_OR_ATOM_EXPECTED:
      if(line[0] == 'a' || line[0] == 'A') {
        state = ATOM_EXPECTED;
        /* fininalize current atom data */
        if(currAtom == NULL || shellBaseIndex < 0) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: (currAtom == NULL || shellBaseIndex < 0) in line %d %s", lineNo, fileName);
          return -1;
	}
        currAtom->noOfShells = shellBaseIndex;
        break;
      } /* else fall through */

    case SHELL_EXPECTED:
      if(shellBaseIndex < 0) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: (shellBaseIndex < 0) in line %d %s", lineNo, fileName);
        return -1;
      }
      int noOfExponents, noOfShellsCurrBatch;
      if(sscanf(line, "%d %d %d",
                &noOfExponents, &noOfShellsCurrBatch, &dummy ) != 3)
        {
          do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: "
                    "Shell data expected in line %d:\n%s", lineNo, line);
          return -1;
        }
      if(noOfExponents <= 0)
        {
          do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: "
                    "(noOfExponents <= 0) in line %d %s", lineNo, fileName);
          return -1;
        }
      if(noOfExponents >= MAX_NO_OF_CONTR)
        {
          do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: "
                    "(noOfExponents >= MAX_NO_OF_CONTR) in line %d",
                    lineNo);
          return -1;
        }
      if(noOfShellsCurrBatch < 0)
        {
          do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: "
                    "(noOfShellsCurrBatch < 0) in line %d", lineNo);
          return -1;
        }
      if(noOfShellsCurrBatch == 0) {
        /*  special case: uncontracted, expect only one column */
        noOfShellsCurrBatch = noOfExponents;
        uncontracted = 1;
      } else uncontracted = 0;

      if(shellBaseIndex + noOfShellsCurrBatch >= MAX_NO_OF_SHELLS_PER_ATOM)
        {
          do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: too many shells.");
          return -1;
        }
      /* initialize shell data. Set the contraction count to its upper
         limit. remove_zeros() will later check for a better, lower
         value. */
      for(int i = 0; i < noOfShellsCurrBatch; i++) {
	if(currAtom == NULL || spdf < 0) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: (currAtom == NULL || spdf < 0) in line %d %s", lineNo, fileName);
	  return -1;
	}
	currAtom->shells[shellBaseIndex+i].type = spdf;
	currAtom->shells[shellBaseIndex+i].contrCount = noOfExponents;
      }
      expNo = 0;
      state = CONTRACTION_BLOCK_EXPECTED;
      if(print_raw)
        do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
                  "Block for L=%d primitives: %d contracted: %d",
                  spdf, noOfExponents, noOfShellsCurrBatch);

      lineConsumed = 1;
      break;

    case CONTRACTION_BLOCK_EXPECTED:
      currExponent = atof(line);
      if(currExponent <= 0) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error: (currExponent <= 0) in line %d", lineNo);
	return -1;
      }
      if(currAtom == NULL || shellBaseIndex < 0 || expNo < 0) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: (currAtom == NULL || shellBaseIndex < 0 || expNo < 0) in line %d %s", lineNo, fileName);
        return -1;
      }
      if(uncontracted) {
        for(int i = 0; i < noOfShellsCurrBatch; i++) {
          currAtom->shells[shellBaseIndex+i].exponentList[expNo] =
            currExponent;
          currAtom->shells[shellBaseIndex+i].coeffList[expNo] =
            i == expNo ? 1.0 : 0.0;
        }
      } else {
        int idx = 0;
        /* skip exponent */
        while(line[idx] && isspace(line[idx]))  idx++;
        for(int i = 0; i < noOfShellsCurrBatch; i++) {
          currAtom->shells[shellBaseIndex+i].exponentList[expNo] =
            currExponent;
          while(line[idx] && !isspace(line[idx])) idx++;
          while(line[idx] && isspace(line[idx]))  idx++;
          if( !line[idx] ) {
	    /* Second line begins when we are about to read 7th
	       contraction coefficient (i=6), third line for 14th
	       (i=13), fourth for i=20 etc. If this pattern does not
	       match, warn the user. */
	    if( i != 6 && (i+1) % 7 != 0 )
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "WARN: line %d has trailing data: '%s'"
			"non-conformant basis set file for i=%d.",
			lineNo, line + idx, i);
	    if(fgets(line, sizeof(line), f) == NULL) {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "reading error when continuing shell data.");
	      return -1;
	    }
	    lineNo++;
	    idx = 0;
	    while(line[idx] && isspace(line[idx]))  idx++;
	  }
          currAtom->shells[shellBaseIndex+i].coeffList[expNo] =
            atof(line + idx);
        }  /*  END FOR i */
      }
      if(print_raw) {
        char line[256], eee[20];
        line[0] = '\0';
        for(int i = 0; i<noOfShellsCurrBatch; i++) {
          sprintf(eee, "%10.5f",
                  (double)currAtom->shells[shellBaseIndex+i].coeffList[expNo]);
          strcat(line, eee);
        }
        do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
                  "%d %12.6f: %s", expNo, (double)currExponent, line);
      }
      if(++expNo == noOfExponents) {
        remove_zeros(currAtom, shellBaseIndex, noOfShellsCurrBatch);
        shellBaseIndex += noOfShellsCurrBatch;
        spdf++;
        state = SHELL_OR_ATOM_EXPECTED;
      }
      lineConsumed = 1;
      break;
    case END_PARSING:
      /*  do nothing */
      break;
    }
  } while(state != END_PARSING);
  fclose(f);

  /* fininalize current atom data */
  if(currAtom == NULL || shellBaseIndex < 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file: (currAtom == NULL || shellBaseIndex < 0) in line %d %s", lineNo, fileName);
    return -1;
  }
  currAtom->noOfShells = shellBaseIndex;

  /*  Postprocessing... */
  /*  set shell ID for each shell in basis set */
  int currShellID = 0;
  for(int i = 0; i < MAX_NO_OF_ATOM_TYPES; i++) {
    int noOfShells = result.atoms[i].noOfShells;
    for(int j = 0; j < noOfShells; j++) {
      basisset_shell_struct* currShell = &result.atoms[i].shells[j];
      currShellID++;
      currShell->shell_ID = currShellID;
    } /*  END FOR j */
  } /*  END FOR i */
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "total number of shells in basis set: %i", currShellID);
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "Basis set file '%s' processed OK, noOfAtomTypes = %i",
            fileName, noOfAtomTypes);

  return 0;
}
