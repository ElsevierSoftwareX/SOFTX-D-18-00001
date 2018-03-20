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

/** @file density_description_file.h

    @brief An interface file for writing and reading density matrices
    to/from a file, including basis set information.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef DENSITY_DESC_FILE
#define  DENSITY_DESC_FILE 1

#include "basisinfo.h"


int ddf_writeShellListAndDensityMatricesToFile(const BasisInfoStruct & basisInfo,
					       int noOfDensityMatrices,
					       ergo_real** densityMatrixList,
					       const char* fileName);

typedef struct
{
  long nvalues;
  int* rowind;
  int* colind;
  ergo_real* values;
} matrix_description_struct;

/** Writes basisInfo and sparse matrices in a format that can be later
    read by ddf_load_density.
*/

int ddf_writeShellListAndDensityMatricesToFile_sparse(const BasisInfoStruct & basisInfo,
						      int noOfDensityMatrices,
						      matrix_description_struct* densityMatrixList,
						      const char* fileName);


/** Function opens fileName, fills in basisInfo (which has
   to be allocated and nullified), allocates densityMatrixList and
   reads density matrix or at most two matrices and puts it/them in
   densityMatrixList.
*/
int ddf_load_density(const char *densityFileName,
		     int noOfDensityMatrices,
		     const IntegralInfo& integralInfo,
		     BasisInfoStruct **basisInfo,
		     ergo_real **densityMatrixList);

/** Function opens fileName, fills in basisInfo (which has
   to be allocated and nullified), allocates densityMatrixList and
   reads density matrix or at most two matrices and puts it/them in
   densityMatrixList.
*/
int ddf_load_density_sparse(const char *densityFileName,
			    const IntegralInfo& integralInfo,
			    BasisInfoStruct **basisInfo,
			    int *noOfDensitiesRead,
			    int** rowindList,
			    int** colindList,
			    ergo_real** valuesList,
			    long* nvaluesList);


#endif /*  DENSITY_DESC_FILE */
