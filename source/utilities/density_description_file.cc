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

/** @file density_description_file.cc

    @brief An interface file for writing and reading density matrices
    to/from a file, including basis set information.

    @author: Elias Rudberg <em>responsible</em>
*/

#define _LARGEFILE_SOURCE 1
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <errno.h>
#include <string>

#include "density_description_file.h"
#include "densfromf_full.h"
#include "integrals_general.h"
#include "matrix_algebra.h"
#include "memorymanag.h"
#include "output.h"
#include "utilities.h"

// Version number.
// This is written as part of the header.
// Should be changed if the implementation is changed in such a way
// that old density files will no longer work.
const int DENSITY_FILE_VERSION_NUMBER = 10002;


// Matrix storage types
#define MATRIX_STORAGE_TYPE_FULL      1
#define MATRIX_STORAGE_TYPE_TRIANGLE  2
#define MATRIX_STORAGE_TYPE_VECTORS   3

const ergo_real THRESHOLD_FOR_VECTOR_STORAGE = 1e-12;


// File header.
// We want this to contain an even number of integers,
// to avoid compilers inserting padding.
// We use double instead of size_t for sizes to ensure compatibility of density files.
typedef struct
{
  int densityFileVersion;
  int typeOfMatrixStorage;
  int noOfShells;
  int noOfBasisFuncs;
  int noOfDensityMatrices;
  int padding_int;
  double matrixSize_1;
  double matrixSize_2;
  double fileSizeInBytes;
} densityFileHeaderStruct;


/** stores the upper triangle of a square matrix in the specified
    memory block.  
    @param p memory block
    @param n matrix dimension
    @param matrix square matrix
*/
static void
ddf_store_triangular_matrix(char* p, int n, const ergo_real* matrix)
{
  ergo_real* resultPtr = (ergo_real*)p;
  long count = 0;
  for(long i = 0; i < n; i++)
    for(long j = i; j < n; j++)
      {
	resultPtr[count] = matrix[i*n+j];
	count++;
      }
}


/** stores the upper triangle of the matrix given by the
    matrix_description_struct in the specified memory block.
    @param p memory block
    @param n matrix dimension
    @param matrix matrix
*/
static void
ddf_store_triangular_matrix_sparse(char* p, long n, matrix_description_struct* matrix)
{
  ergo_real* resultPtr = (ergo_real*)p;
  long totCount = n*(n+1) / 2;
  memset(resultPtr, 0, totCount*sizeof(ergo_real));
  for(long i = 0; i < matrix->nvalues; i++)
    {
      long row = matrix->rowind[i];
      long col = matrix->colind[i];
      ergo_real value = matrix->values[i];
      // Make sure (row,col) refers to upper triangle.
      if(col < row)
	{
	  long tmp = row;
	  row = col;
	  col = tmp;
	}
      long index = row*n+col-row*(row+1)/2;
      resultPtr[index] = value;
    }
}


static void
ddf_get_triangular_matrix_from_storage(const char* p, long n, ergo_real* resultMatrix)
{
  const ergo_real* sourcePtr = (ergo_real*)p;
  long count = 0;
  for(long i = 0; i < n; i++)
    for(long j = i; j < n; j++)
      {
	ergo_real a = sourcePtr[count];
	resultMatrix[i*n+j] = a;
	resultMatrix[j*n+i] = a;
	count++;
      }  
}

static void
ddf_get_triangular_matrix_from_storage_sparse(const char* p, 
					      int n, 
					      int* rowind,
					      int* colind,
					      ergo_real* values)
{
  const ergo_real* sourcePtr = (ergo_real*)p;
  long count = 0;
  for(long i = 0; i < n; i++)
    for(long j = i; j < n; j++)
      {
	ergo_real a = sourcePtr[count];
	rowind[count] = i;
	colind[count] = j;
	values[count] = a;
	count++;
      }
}

static int
ddf_get_nvalues_symm_matrix(int n, const ergo_real* matrix)
{
  long nvalues = 0;
  for(long i = 0; i < n; i++)
    for(long j = i; j < n; j++)
      {
	if(template_blas_fabs(matrix[i*n+j]) > THRESHOLD_FOR_VECTOR_STORAGE)
	  nvalues++;
      }
  return nvalues;
}

static int
ddf_store_matrix_as_vectors(char* p, int n, const ergo_real* matrix, size_t sizeInBytes)
{
  long nvalues = ddf_get_nvalues_symm_matrix(n, matrix);
  if(sizeInBytes != nvalues*(sizeof(ergo_real)+2*sizeof(int)))
    return -1;
  int* rowind = (int*)p;
  int* colind = &rowind[nvalues];
  // The values follow as a sequence of ergo_real after rowind and colind, but we use a char ptr and memcpy to avoid memory alignment problems e.g. in quad precision case.
  char* valuesPtr = (char*)(&rowind[2*nvalues]);
  long count = 0;
  for(long i = 0; i < n; i++)
    for(long j = i; j < n; j++)
      {
	if(template_blas_fabs(matrix[i*n+j]) > THRESHOLD_FOR_VECTOR_STORAGE)
	  {
	    if(count >= nvalues)
	      return -1;
	    rowind[count] = i;
	    colind[count] = j;
	    memcpy(valuesPtr+count*sizeof(ergo_real), &matrix[i*n+j], sizeof(ergo_real)); // Use memcpy to avoid memory alignment problems e.g. in quad precision case.
	    count++;
	  }
      } // END FOR i j
  return 0;
}


static int
ddf_store_matrix_as_vectors_sparse(char* p, int n, matrix_description_struct* matrix, size_t sizeInBytes)
{
  long nvalues = matrix->nvalues;
  if(sizeInBytes != nvalues*(sizeof(ergo_real)+2*sizeof(int)))
    return -1;
  int* rowind = (int*)p;
  int* colind = &rowind[nvalues];
  ergo_real* values = (ergo_real*)(&rowind[2*nvalues]);
  for(long i = 0; i < nvalues; i++)
    {
      rowind[i] = matrix->rowind[i];
      colind[i] = matrix->colind[i];
      values[i] = matrix->values[i];
    }
  return 0;
}


static int
ddf_get_matrix_from_vector_storage(const char* p, size_t sizeInBytes, int n, ergo_real* resultMatrix)
{
  int nBytesPerElement = sizeof(ergo_real)+2*sizeof(int);
  if(sizeInBytes % nBytesPerElement != 0)
    return -1;
  long nvalues = sizeInBytes / nBytesPerElement;
  int* rowind = (int*)p;
  int* colind = &rowind[nvalues];
  // The values follow as a sequence of ergo_real after rowind and colind, but we use a char ptr and memcpy to avoid memory alignment problems e.g. in quad precision case.
  char* valuesPtr = (char*)(&rowind[2*nvalues]);
  memset(resultMatrix, 0, n*n*sizeof(ergo_real));
  for(long i = 0; i < nvalues; i++)
    {
      long row = rowind[i];
      long col = colind[i];
      ergo_real value = 0;
      memcpy(&value, &valuesPtr[i*sizeof(ergo_real)], sizeof(ergo_real)); // Use memcpy to avoid memory alignment problems e.g. in quad precision case.
      resultMatrix[row*n+col] = value;
      resultMatrix[col*n+row] = value;
    }
  return 0;
}

static int
ddf_get_matrix_from_vector_storage_sparse(const char* p, size_t sizeInBytes, int n, int* rowind2, int* colind2, ergo_real* values2)
{
  int nBytesPerElement = sizeof(ergo_real)+2*sizeof(int);
  if(sizeInBytes % nBytesPerElement != 0)
    return -1;
  long nvalues = sizeInBytes / nBytesPerElement;
  int* rowind = (int*)p;
  int* colind = &rowind[nvalues];
  ergo_real* values = (ergo_real*)(&rowind[2*nvalues]);
  for(long i = 0; i < nvalues; i++)
    {
      rowind2[i] = rowind[i];
      colind2[i] = colind[i];
      values2[i] = values[i];
    }
  return 0;
}

static size_t
ddf_get_matrix_storage_size(int storageType, int n, const ergo_real* matrix)
{
  size_t storageSize = 0;
  switch(storageType)
    {
    case MATRIX_STORAGE_TYPE_FULL:
      storageSize = n*n*sizeof(ergo_real);
      break;
    case MATRIX_STORAGE_TYPE_TRIANGLE:
      storageSize = (n * (n+1) / 2) * sizeof(ergo_real);
      break;
    case MATRIX_STORAGE_TYPE_VECTORS:
      storageSize = ddf_get_nvalues_symm_matrix(n, matrix) * (sizeof(ergo_real) + 2*sizeof(int));
      break;
    default:
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_get_matrix_storage_size: unknown storage type %i", storageType);
      return 0;
    }
  return storageSize;
}


static size_t
ddf_get_matrix_storage_size_sparse(int storageType, int n, matrix_description_struct* matrix)
{
  size_t storageSize = 0;
  switch(storageType)
    {
    case MATRIX_STORAGE_TYPE_FULL:
      storageSize = n*n*sizeof(ergo_real);
      break;
    case MATRIX_STORAGE_TYPE_TRIANGLE:
      storageSize = (n * (n+1) / 2) * sizeof(ergo_real);
      break;
    case MATRIX_STORAGE_TYPE_VECTORS:
      storageSize = matrix->nvalues * (sizeof(ergo_real) + 2*sizeof(int));
      break;
    default:
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_get_matrix_storage_size_sparse: unknown storage type %i", storageType);
      return 0;
    }
  return storageSize;
}


static int
ddf_store_matrix(char* p, size_t sizeInBytes, int n, const ergo_real* matrix, int storageType)
{
  size_t sizeNeeded = 0;
  switch(storageType)
    {
    case MATRIX_STORAGE_TYPE_FULL:
      sizeNeeded = n*n*sizeof(ergo_real);
      if(sizeNeeded != sizeInBytes)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_store_matrix: (sizeNeeded != sizeInBytes)");
	  return -1;
	}
      memcpy(p, matrix, sizeNeeded);
      break;
    case MATRIX_STORAGE_TYPE_TRIANGLE:
      sizeNeeded = (n * (n+1) / 2) * sizeof(ergo_real);
      if(sizeNeeded != sizeInBytes)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_store_matrix: (sizeNeeded != sizeInBytes)");
	  return -1;
	}
      ddf_store_triangular_matrix(p, n, matrix);
      break;
    case MATRIX_STORAGE_TYPE_VECTORS:
      if(ddf_store_matrix_as_vectors(p, n, matrix, sizeInBytes) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_store_matrix_as_vectors.");
	  return -1;
	}
      break;
    default:
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_store_matrix: unknown storage type %i", storageType);
      return -1;
    }
  return 0;
}


static int
ddf_store_matrix_sparse(char* p, 
			size_t sizeInBytes, 
			int n, 
			matrix_description_struct* matrix,
			int storageType)
{
  size_t sizeNeeded = 0;
  switch(storageType)
    {
    case MATRIX_STORAGE_TYPE_FULL:
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_store_matrix_sparse: MATRIX_STORAGE_TYPE_FULL not supported.");
      return -1;
      break;
    case MATRIX_STORAGE_TYPE_TRIANGLE:
      sizeNeeded = (n * (n+1) / 2) * sizeof(ergo_real);
      if(sizeNeeded != sizeInBytes)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_store_matrix: (sizeNeeded != sizeInBytes)");
	  return -1;
	}
      ddf_store_triangular_matrix_sparse(p, n, matrix);
      break;
    case MATRIX_STORAGE_TYPE_VECTORS:
      if(ddf_store_matrix_as_vectors_sparse(p, n, matrix, sizeInBytes) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_store_matrix_as_vectors.");
	  return -1;
	}
      break;
    default:
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_store_matrix: unknown storage type %i", storageType);
      return -1;
    }
  return 0;
}


static int
ddf_get_matrix_nvalues_from_storage(const char* p, size_t sizeInBytes, int n,
				    long* result_nvalues, int storageType)
{
  switch(storageType)
    {
    case MATRIX_STORAGE_TYPE_FULL:
      // In this case the needed nvalues is n*n
      *result_nvalues = n*n;
      break;
    case MATRIX_STORAGE_TYPE_TRIANGLE:
      *result_nvalues = (n * (n+1) / 2);
      break;
    case MATRIX_STORAGE_TYPE_VECTORS:
      {
	int nBytesPerElement = sizeof(ergo_real)+2*sizeof(int);
	if(sizeInBytes % nBytesPerElement != 0)
	  return -1;
	*result_nvalues = sizeInBytes / nBytesPerElement;
      }
      break;
    default:
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_get_matrix_nvalues_from_storage: unknown storage type %i", storageType);
      return -1;
    }
  return 0;
}


static int
ddf_get_matrix_from_storage(const char* p, size_t sizeInBytes, int n,
                            ergo_real* resultMatrix, int storageType)
{
  size_t sizeNeeded = 0;
  switch(storageType)
    {
    case MATRIX_STORAGE_TYPE_FULL:
      sizeNeeded = n*n*sizeof(ergo_real);
      if(sizeNeeded != sizeInBytes)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_get_matrix_from_storage: (sizeNeeded != sizeInBytes)");
	  return -1;
	}
      memcpy(resultMatrix, p, sizeNeeded);
      break;
    case MATRIX_STORAGE_TYPE_TRIANGLE:
      sizeNeeded = (n * (n+1) / 2) * sizeof(ergo_real);
      if(sizeNeeded != sizeInBytes)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_store_matrix: (sizeNeeded != sizeInBytes)");
	  return -1;
	}
      ddf_get_triangular_matrix_from_storage(p, n, resultMatrix);
      break;
    case MATRIX_STORAGE_TYPE_VECTORS:
      if(ddf_get_matrix_from_vector_storage(p, sizeInBytes, n, resultMatrix) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_get_matrix_from_vector_storage.");
	  return -1;
	}	
      break;
    default:
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_get_matrix_from_storage: unknown storage type %i", storageType);
      return -1;
    }
  return 0;
}


static int
ddf_get_matrix_from_storage_sparse(const char* p, 
				   size_t sizeInBytes, 
				   int n,
				   int* rowind,
				   int* colind,
				   ergo_real* values,
				   int storageType)
{
  size_t sizeNeeded = 0;
  switch(storageType)
    {
    case MATRIX_STORAGE_TYPE_FULL:
      {
	sizeNeeded = n*n*sizeof(ergo_real);
	if(sizeNeeded != sizeInBytes)
	  {
	    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_get_matrix_from_storage: (sizeNeeded != sizeInBytes)");
	    return -1;
	  }
	ergo_real* A = (ergo_real*)p;
	long count = 0;
	for(long i = 0; i < n; i++)
	  for(long j = 0; j < n; j++)
	    {
	      rowind[count] = i;
	      colind[count] = j;
	      values[count] = A[i*n+j];
	      count++;
	    }
      }
      break;
    case MATRIX_STORAGE_TYPE_TRIANGLE:
      sizeNeeded = (n * (n+1) / 2) * sizeof(ergo_real);
      if(sizeNeeded != sizeInBytes)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_store_matrix: (sizeNeeded != sizeInBytes)");
	  return -1;
	}
      ddf_get_triangular_matrix_from_storage_sparse(p, n, rowind, colind, values);
      break;
    case MATRIX_STORAGE_TYPE_VECTORS:
      if(ddf_get_matrix_from_vector_storage_sparse(p, sizeInBytes, n, rowind, colind, values) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_get_matrix_from_vector_storage.");
	  return -1;
	}
      break;
    default:
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_get_matrix_from_storage: unknown storage type %i", storageType);
      return -1;
    }
  return 0;
}


int 
ddf_writeShellListAndDensityMatricesToFile(const BasisInfoStruct & basisInfo,
					   int noOfDensityMatrices,
					   ergo_real** densityMatrixList,
					   const char* fileName) {
  int n = basisInfo.noOfBasisFuncs;
  int nShells = basisInfo.noOfShells;
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "ddf_writeShellListAndDensityMatricesToFile, n = %d, nShells = %d", n, nShells);
  if(n <= 0 || nShells <= 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in writeShellListAndDensityMatricesToFile: "
	      "(n <= 0 || nShells <= 0)");
    return -1;
  }
  if(noOfDensityMatrices <= 0 || noOfDensityMatrices > 2) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in writeShellListAndDensityMatricesToFile: "
	      "(noOfDensityMatrices <= 0 || noOfDensityMatrices > 2)");
    return -1;
  }
  FILE* f = fopen(fileName, "wb");
  if(f == NULL) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error opening file '%s' for writing", fileName);
    return -1;
  }

  // decide which storage format to use
  int matrixStorageType;
  if(n < 22000) {
    if(ddf_get_matrix_storage_size(MATRIX_STORAGE_TYPE_TRIANGLE, n, densityMatrixList[0]) < ddf_get_matrix_storage_size(MATRIX_STORAGE_TYPE_VECTORS, n, densityMatrixList[0])) {
      matrixStorageType = MATRIX_STORAGE_TYPE_TRIANGLE;
      do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "Choosing matrixStorageType = MATRIX_STORAGE_TYPE_TRIANGLE.");
    }
    else
      {
	matrixStorageType = MATRIX_STORAGE_TYPE_VECTORS;
	do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "Choosing matrixStorageType = MATRIX_STORAGE_TYPE_VECTORS.");
      }
  }
  else {
    // triangle storage not possible, use vector storage.
    matrixStorageType = MATRIX_STORAGE_TYPE_VECTORS;
    do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "Choosing matrixStorageType = MATRIX_STORAGE_TYPE_VECTORS.");
  } 

  // Get sizes needed for matrix storage.
  size_t matrixStorageSizeList[2];
  matrixStorageSizeList[0] = 0;
  matrixStorageSizeList[1] = 0;
  for(int i = 0; i < noOfDensityMatrices; i++) {
    size_t sizeInBytes = ddf_get_matrix_storage_size(matrixStorageType, n, densityMatrixList[i]);
    if(sizeInBytes == 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_get_matrix_storage_size");
      return -1;
    }
    matrixStorageSizeList[i] = sizeInBytes;
  }  

  size_t matrixStorageSizeTotal = 0;
  for(int i = 0; i < noOfDensityMatrices; i++)
    matrixStorageSizeTotal += matrixStorageSizeList[i];

  // Compute file size
  size_t fileSize = 
    sizeof(densityFileHeaderStruct) +
    nShells * sizeof(ShellSpecStruct) +
    matrixStorageSizeTotal;
  densityFileHeaderStruct fileHeader;

  // Create header
  fileHeader.densityFileVersion = DENSITY_FILE_VERSION_NUMBER;
  fileHeader.typeOfMatrixStorage = matrixStorageType;
  fileHeader.noOfShells = nShells;
  fileHeader.noOfBasisFuncs = n;
  fileHeader.noOfDensityMatrices = noOfDensityMatrices;
  fileHeader.padding_int = 0; /* Set padding_int also to make valgrind happy. */
  fileHeader.matrixSize_1 = matrixStorageSizeList[0];
  fileHeader.matrixSize_2 = matrixStorageSizeList[1];
  fileHeader.fileSizeInBytes = fileSize;

  // prepare data to write
  char* buffer = (char*)ergo_malloc(fileSize);
  char* p = buffer;
  memcpy(p, &fileHeader, sizeof(densityFileHeaderStruct));
  p += sizeof(densityFileHeaderStruct);
  memcpy(p, basisInfo.shellList, nShells * sizeof(ShellSpecStruct));
  p += nShells * sizeof(ShellSpecStruct);
  for(int i = 0; i < noOfDensityMatrices; i++) {
    if(ddf_store_matrix(p, matrixStorageSizeList[i], n, densityMatrixList[i], matrixStorageType) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_store_matrix");
      return -1;
    }
    p += matrixStorageSizeList[i];
  } // END FOR i
  if(size_t(p - buffer) != fileSize) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error: ((p - buffer) != fileSize)");
    return -1;
  }

  // write to file
  size_t noOfBytesWritten = fwrite(buffer, 1, fileSize, f);
  fclose(f);  
  if(noOfBytesWritten != fileSize) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in fwrite, fileSize = %i, noOfBytesWritten = %i",
	      fileSize, noOfBytesWritten);
    return -1;
  }

  ergo_free(buffer);

  // return 0 to indicate success
  return 0;
}



int 
ddf_writeShellListAndDensityMatricesToFile_sparse(const BasisInfoStruct & basisInfo,
						  int noOfDensityMatrices,
						  matrix_description_struct* densityMatrixList,
						  const char* fileName)
{
  int n = basisInfo.noOfBasisFuncs;

  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "ddf_writeShellListAndDensityMatricesToFile_sparse, n = %6i, densityMatrixList[0].nvalues = %11i", n, densityMatrixList[0].nvalues);

  int nShells = basisInfo.noOfShells;
  if(n <= 0 || nShells <= 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in writeShellListAndDensityMatricesToFile: "
	      "(n <= 0 || nShells <= 0)");
    return -1;
  }
  if(noOfDensityMatrices <= 0 || noOfDensityMatrices > 2) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in writeShellListAndDensityMatricesToFile: "
	      "(noOfDensityMatrices <= 0 || noOfDensityMatrices > 2)");
    return -1;
  }

  /* We create a tmp file and rename it on success. This is useful
     when a file of identical name exists: this method guarantees that
     the created density file does not get corrupted when the process
     is interrupted mid-stream. */
  std::string tmpFName(fileName);
  tmpFName += ".XXXXXX";
  int tmpFD = mkstemp(&tmpFName[0]);
  FILE* f;

  if(tmpFD == -1 || (f = fdopen(tmpFD, "wb")) == NULL)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
                "error opening temporary file '%s' for writing",
                tmpFName.c_str());
      return -1;
    }
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN,
	    "opened temporary file '%s' for writing",
	    tmpFName.c_str());
  
  // decide which storage format to use
  int matrixStorageType;
  if(n < 22000)
    {
      if(ddf_get_matrix_storage_size_sparse(MATRIX_STORAGE_TYPE_TRIANGLE, n, &densityMatrixList[0]) < ddf_get_matrix_storage_size_sparse(MATRIX_STORAGE_TYPE_VECTORS, n, &densityMatrixList[0]))
        {
          matrixStorageType = MATRIX_STORAGE_TYPE_TRIANGLE;
          do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "Choosing matrixStorageType = MATRIX_STORAGE_TYPE_TRIANGLE.");
        }
      else
        {
          matrixStorageType = MATRIX_STORAGE_TYPE_VECTORS;
          do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "Choosing matrixStorageType = MATRIX_STORAGE_TYPE_VECTORS.");
        }
    }
  else
    {
      // triangle storage not possible, use vector storage.
      matrixStorageType = MATRIX_STORAGE_TYPE_VECTORS;
      do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "Choosing matrixStorageType = MATRIX_STORAGE_TYPE_VECTORS.");
    } 

  // Get sizes needed for matrix storage.
  size_t matrixStorageSizeList[2];
  matrixStorageSizeList[0] = 0;
  matrixStorageSizeList[1] = 0;
  for(int i = 0; i < noOfDensityMatrices; i++)
    {
      size_t sizeInBytes = ddf_get_matrix_storage_size_sparse(matrixStorageType, n, &densityMatrixList[i]);
      if(sizeInBytes <= 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_get_matrix_storage_size_sparse");
	  return -1;
	}
      matrixStorageSizeList[i] = sizeInBytes;
    }
  
  size_t matrixStorageSizeTotal = 0;
  for(int i = 0; i < noOfDensityMatrices; i++)
    matrixStorageSizeTotal += matrixStorageSizeList[i];

  // Compute file size
  size_t fileSize = 
    sizeof(densityFileHeaderStruct) +
    nShells * sizeof(ShellSpecStruct) +
    matrixStorageSizeTotal;

  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "fileSize = %15.0f bytes = %10.3f MegaBytes", (double)fileSize, (double)fileSize / 1000000);

  densityFileHeaderStruct fileHeader;

  // Create header
  fileHeader.densityFileVersion = DENSITY_FILE_VERSION_NUMBER;
  fileHeader.typeOfMatrixStorage = matrixStorageType;
  fileHeader.noOfShells = nShells;
  fileHeader.noOfBasisFuncs = n;
  fileHeader.noOfDensityMatrices = noOfDensityMatrices;
  fileHeader.padding_int = 0; /* Set padding_int also to make valgrind happy. */
  fileHeader.matrixSize_1 = matrixStorageSizeList[0];
  fileHeader.matrixSize_2 = matrixStorageSizeList[1];
  fileHeader.fileSizeInBytes = fileSize;

  // prepare data to write
  char* buffer = (char*)ergo_malloc(fileSize);
  char* p = buffer;
  memcpy(p, &fileHeader, sizeof(densityFileHeaderStruct));
  p += sizeof(densityFileHeaderStruct);
  memcpy(p, basisInfo.shellList, nShells * sizeof(ShellSpecStruct));
  p += nShells * sizeof(ShellSpecStruct);
  for(int i = 0; i < noOfDensityMatrices; i++) {
    if(ddf_store_matrix_sparse(p, matrixStorageSizeList[i], n, &densityMatrixList[i], matrixStorageType) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_store_matrix");
      return -1;
    }
    p += matrixStorageSizeList[i];
  } // END FOR i
  if(size_t(p - buffer) != fileSize) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error: ((p - buffer) != fileSize)");
    return -1;
  }

  // write to file
  size_t noOfBytesWritten = fwrite(buffer, 1, fileSize, f);
  if(fsync(fileno(f)) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in fsync on %s: %s",
	      fileName, strerror(errno));
    return -1;
  }
  if(fclose(f) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in fclose on %s: %s",
	      fileName, strerror(errno));
    return -1;
  }

  if(noOfBytesWritten != fileSize) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in fwrite, fileSize = %i, noOfBytesWritten = %i",
	      fileSize, noOfBytesWritten);
    return -1;
  }
  if(rename(tmpFName.c_str(), fileName) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error renaming %s to %s: %s\n",
	      tmpFName.c_str(), fileName, strerror(errno));
    return -1;
  }

  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "renamed file '%s' to '%s'",
	    tmpFName.c_str(), fileName);

  /** Data loss was observed with large files on AFS, we do extra verification to detect it early on... */
  struct stat fileStats;
  if(stat(fileName, &fileStats) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "Cannot stat saved file %s: %s\n",
	      fileName, strerror(errno));
    return -1;
  }
  if( size_t(fileStats.st_size) != fileSize) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "File %s size got altered from %lu to %lu\n",
	      fileName, (unsigned long)fileSize, (unsigned long)fileStats.st_size);
    return -1;
  }
    
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "ddf_writeShellListAndDensityMatricesToFile_sparse freeing buffer.");
  ergo_free(buffer);

  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "ddf_writeShellListAndDensityMatricesToFile_sparse returning OK.");
  // return 0 to indicate success
  return 0;
}



static int 
ddf_load_density_getSizes(const char* fileName,
                          int* result_noOfShells,
                          int* result_noOfBasisFuncs,
			  int* result_noOfDensitiesOnFile,
			  long* result_noOfValuesList)
{
  FILE* f = fopen(fileName, "rb");
  if(f == NULL)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error opening file '%s'", fileName);
      return -1;
    }
  
  // check file size
  fseeko(f, 0L, SEEK_END);
  size_t fileSize = size_t(ftello(f));
  if(fileSize < sizeof(densityFileHeaderStruct))
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
		"error in ddf_load_density_getSizes: (fileSize < sizeof(densityFileHeaderStruct)) %lu < %lu",
		(unsigned long)fileSize, (unsigned long)sizeof(densityFileHeaderStruct));
      return -1;
    }
  fseeko(f, 0L, SEEK_SET);
  char* buffer = (char*)ergo_malloc(fileSize);
  if(fread(buffer, 1, fileSize, f) != fileSize)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error reading file '%s'", fileName);
      return -1;
    }
  fclose(f);  
  
  // OK, file contents are now in buffer.
  densityFileHeaderStruct* headerPtr = (densityFileHeaderStruct*)buffer;
  if(headerPtr->densityFileVersion != DENSITY_FILE_VERSION_NUMBER)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_load_density_getSizes: (headerPtr->densityFileVersion != DENSITY_FILE_VERSION_NUMBER)");
      return -1;
    }
  if(headerPtr->fileSizeInBytes != fileSize)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
		"error in ddf_load_density_getSizes: (headerPtr->fileSizeInBytes != fileSize) %lu %lu",
		(unsigned long)headerPtr->fileSizeInBytes, (unsigned long)fileSize);
      return -1;
    }
  int n = headerPtr->noOfBasisFuncs;
  int nShells = headerPtr->noOfShells;

  size_t fileSize_expected = 
    sizeof(densityFileHeaderStruct) + 
    nShells * sizeof(ShellSpecStruct) +
    (size_t)headerPtr->matrixSize_1 + (size_t)headerPtr->matrixSize_2;
  if(fileSize_expected != fileSize)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_load_density_getSizes: (fileSize_expected != fileSize) %lu != %lu",
		(unsigned long)fileSize_expected, (unsigned long)fileSize);
      return -1;
    }
  
  *result_noOfShells = nShells;
  *result_noOfBasisFuncs = n;
  *result_noOfDensitiesOnFile = headerPtr->noOfDensityMatrices;
  if(result_noOfValuesList != NULL)
    {
      // Get result_noOfValuesList.
      char* p = buffer + sizeof(densityFileHeaderStruct);
  
      // skip shell list
      p += nShells*sizeof(ShellSpecStruct);

      // get nvalues for stored matrices
      size_t matrixSizeList[2];
      matrixSizeList[0] = (size_t)headerPtr->matrixSize_1;
      matrixSizeList[1] = (size_t)headerPtr->matrixSize_2;
      for(int i = 0; i < headerPtr->noOfDensityMatrices; i++)
	{
	  size_t currSize = matrixSizeList[i];
	  if(ddf_get_matrix_nvalues_from_storage(p, currSize, n, &result_noOfValuesList[i],
						 headerPtr->typeOfMatrixStorage) != 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_get_matrix_nvalues_from_storage");
	      return -1;
	    }
	  p += currSize;
	} // END FOR i
  
      if(size_t(p - buffer) != fileSize)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error: ((p - buffer) != fileSize)");
	  return -1;
	}
    }
  
  ergo_free(buffer);
  return 0;
}


/** ddf_read_shells_and_density_matrices() reads the basis set
    information and requested number of density matrices from a
    specified file. basisInfo needs to be allocated and zeroed in
    advance. densityMatrixList must be properly allocated as well. */
static int
ddf_read_shells_and_density_matrices(BasisInfoStruct* basisInfo,
                                     int noOfDensityMatrices,
                                     ergo_real** densityMatrixList,
                                     const char* fileName)
{
  FILE* f = fopen(fileName, "rb");
  if(f == NULL)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error opening file '%s'", fileName);
      return -1;
    }

  // check file size
  fseeko(f, 0L, SEEK_END);
  size_t fileSize = size_t(ftello(f));
  if(fileSize < sizeof(densityFileHeaderStruct))
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_read_shells_and_density_matrices: (fileSize < sizeof(densityFileHeaderStruct))");
      return -1;
    }
  fseeko(f, 0L, SEEK_SET);
  char* buffer = (char*)ergo_malloc(fileSize);
  size_t readsz = fread(buffer, 1, fileSize, f);
  fclose(f);
  if(readsz != fileSize)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error reading file '%s'", fileName);
      ergo_free(buffer);
      return -1;
    }
  
  // OK, file contents are now in buffer.
  densityFileHeaderStruct* headerPtr = (densityFileHeaderStruct*)buffer;
  if(headerPtr->densityFileVersion != DENSITY_FILE_VERSION_NUMBER)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_read_shells_and_density_matrices: (headerPtr->densityFileVersion != DENSITY_FILE_VERSION_NUMBER)");
      return -1;
    }
  if(headerPtr->fileSizeInBytes != fileSize)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_read_shells_and_density_matrices: (headerPtr->fileSizeInBytes != fileSize)");
      return -1;
    }
  if(headerPtr->noOfDensityMatrices != noOfDensityMatrices)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_read_shells_and_density_matrices: (headerPtr->noOfDensityMatrices != noOfDensityMatrices)");
      return -1;
    }
  int n = headerPtr->noOfBasisFuncs;
  int nShells = headerPtr->noOfShells;

  size_t fileSize_expected = 
    sizeof(densityFileHeaderStruct) + 
    nShells * sizeof(ShellSpecStruct) +
    (size_t)headerPtr->matrixSize_1 + (size_t)headerPtr->matrixSize_2;
  if(fileSize_expected != fileSize)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_read_shells_and_density_matrices: (fileSize_expected != fileSize)");
      return -1;
    }
  
  char* p = buffer + sizeof(densityFileHeaderStruct);
  
  // get shell list
  basisInfo->noOfShells = nShells;
  memcpy(basisInfo->shellList, p, nShells*sizeof(ShellSpecStruct));
  p += nShells*sizeof(ShellSpecStruct);

  // get density matrices
  size_t matrixSizeList[2];
  matrixSizeList[0] = (size_t)headerPtr->matrixSize_1;
  matrixSizeList[1] = (size_t)headerPtr->matrixSize_2;
  for(int i = 0; i < noOfDensityMatrices; i++)
    {
      long currSize = matrixSizeList[i];
      if(ddf_get_matrix_from_storage(p, currSize, n, densityMatrixList[i],
                                     headerPtr->typeOfMatrixStorage) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_get_matrix_from_storage");
	  return -1;
	}
      p += currSize;
    } // END FOR i

  if(size_t(p - buffer) != fileSize)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error: ((p - buffer) != fileSize)");
      return -1;
    }
  
  ergo_free(buffer);
  return 0;
}



static int
ddf_read_shells_and_density_matrices_sparse(BasisInfoStruct** basisInfo,
					    int noOfDensityMatrices,
					    int** rowindList,
					    int** colindList,
					    ergo_real** valuesList,
					    const char* fileName)
{
  FILE* f = fopen(fileName, "rb");
  if(f == NULL)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error opening file '%s'", fileName);
      return -1;
    }

  // check file size
  fseeko(f, 0L, SEEK_END);
  size_t fileSize = size_t(ftello(f));
  if(fileSize < sizeof(densityFileHeaderStruct))
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_read_shells_and_density_matrices: (fileSize < sizeof(densityFileHeaderStruct))");
      return -1;
    }
  fseeko(f, 0L, SEEK_SET);
  char* buffer = (char*)ergo_malloc(fileSize);
  size_t readsz = fread(buffer, 1, fileSize, f);
  fclose(f);
  if(readsz != fileSize)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error reading file '%s'", fileName);
      ergo_free(buffer);
      return -1;
    }
  
  // OK, file contents are now in buffer.
  densityFileHeaderStruct* headerPtr = (densityFileHeaderStruct*)buffer;
  if(headerPtr->densityFileVersion != DENSITY_FILE_VERSION_NUMBER)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_read_shells_and_density_matrices: (headerPtr->densityFileVersion != DENSITY_FILE_VERSION_NUMBER)");
      return -1;
    }
  if(headerPtr->fileSizeInBytes != fileSize)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_read_shells_and_density_matrices: (headerPtr->fileSizeInBytes != fileSize)");
      return -1;
    }
  if(headerPtr->noOfDensityMatrices != noOfDensityMatrices)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_read_shells_and_density_matrices: (headerPtr->noOfDensityMatrices != noOfDensityMatrices)");
      return -1;
    }
  int n = headerPtr->noOfBasisFuncs;
  int nShells = headerPtr->noOfShells;

  size_t fileSize_expected = 
    sizeof(densityFileHeaderStruct) + 
    nShells * sizeof(ShellSpecStruct) +
    (size_t)headerPtr->matrixSize_1 + (size_t)headerPtr->matrixSize_2;
  if(fileSize_expected != fileSize)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_read_shells_and_density_matrices: (fileSize_expected != fileSize)");
      return -1;
    }
  
  char* p = buffer + sizeof(densityFileHeaderStruct);
  
  // get shell list
  ShellSpecStruct* shellList = new ShellSpecStruct[nShells];
  memcpy(shellList, p, nShells*sizeof(ShellSpecStruct));
  // Now check if use_6_d_funcs was used for the basis set.
  int use_6_d_funcs = 0;
  for(int i = 0; i < nShells; i++) {
    if(shellList[i].noOfBasisFuncs == 6)
      use_6_d_funcs = 1;
  }
  // Now create basisInfo object.
  *basisInfo = new BasisInfoStruct(use_6_d_funcs);
  (*basisInfo)->noOfShells = nShells;
  (*basisInfo)->shellList = shellList;
  // Move pointer past shell list
  p += nShells*sizeof(ShellSpecStruct);

  // get density matrices
  size_t matrixSizeList[2];
  matrixSizeList[0] = (size_t)headerPtr->matrixSize_1;
  matrixSizeList[1] = (size_t)headerPtr->matrixSize_2;
  for(int i = 0; i < noOfDensityMatrices; i++)
    {
      size_t currSize = matrixSizeList[i];
      if(ddf_get_matrix_from_storage_sparse(p, currSize, n, rowindList[i], colindList[i], valuesList[i],
					    headerPtr->typeOfMatrixStorage) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_get_matrix_from_storage_sparse");
	  return -1;
	}
      p += currSize;
    } // END FOR i

  if(size_t(p - buffer) != fileSize)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error: ((p - buffer) != fileSize)");
      return -1;
    }
  
  ergo_free(buffer);
  return 0;
}




int
ddf_load_density(const char *densityFileName,
                 int noOfDensityMatrices,
                 const IntegralInfo& integralInfo,
                 BasisInfoStruct **basisInfo,
                 ergo_real **densityMatrix)
{
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, 
	    "getting starting guess from file, noOfDensityMatrices = %i",
            noOfDensityMatrices);
  
  // get basisInfo and density matrix (or matrices) for starting guess
  // FIXME: here it is assumed that use_6_d_funcs was not used for the saved density.
  int use_6_d_funcs = 0;
  *basisInfo = new BasisInfoStruct(use_6_d_funcs);

  int noOfShells_sg = 0;
  int noOfBasisFuncs_sg = 0;
  int noOfDensitiesOnFile = 0;
  if(ddf_load_density_getSizes(densityFileName,
                               &noOfShells_sg,
                               &noOfBasisFuncs_sg,
			       &noOfDensitiesOnFile,
			       NULL) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_load_density_getSizes");
      return -1;
    }
  if(noOfDensitiesOnFile != noOfDensityMatrices)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN,
		"error: (headerPtr->noOfDensityMatrices != noOfDensityMatrices)");
      return -1;
    }


  (*basisInfo)->shellList = new ShellSpecStruct[noOfShells_sg];
  for(int i = 0; i < noOfDensityMatrices; i++)
    densityMatrix[i] = 
      ergo_new(noOfBasisFuncs_sg*noOfBasisFuncs_sg,ergo_real);

  if(ddf_read_shells_and_density_matrices(*basisInfo,
                                          noOfDensityMatrices,
                                          densityMatrix,
                                          densityFileName) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_load_density");
      return -1;
    }
  if((*basisInfo)->get_basis_funcs() != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in get_basis_funcs");
      return -1;
    }
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "get_basis_funcs for starting guess basis set returned OK,"
	    " number of basis funcs: %i",
            (*basisInfo)->noOfBasisFuncs);
  if((*basisInfo)->getSimplePrimitivesAll(integralInfo) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in getSimplePrimitivesAll");
      return -1;
    }
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, 
	    "getSimplePrimitivesAll for starting guess basis set "
            "returned OK, n = %i",
            (*basisInfo)->noOfSimplePrimitives);
  return 0;
}



int 
ddf_load_density_sparse(const char *densityFileName,
			const IntegralInfo& integralInfo,
			BasisInfoStruct **basisInfo,
			int* noOfDensitiesRead,
			int** rowindList,
			int** colindList,
			ergo_real** valuesList,
			long* nvaluesList)
{
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, 
	    "ddf_load_density_sparse(\"%s\")", densityFileName);
  
  // get basisInfo and density matrix (or matrices) for starting guess
  // FIXME: here it is assumed that use_6_d_funcs was not used for the saved density.
  int use_6_d_funcs = 0;
  *basisInfo = new BasisInfoStruct(use_6_d_funcs);

  int noOfShells_sg = 0;
  int noOfBasisFuncs_sg = 0;
  if(ddf_load_density_getSizes(densityFileName,
                               &noOfShells_sg,
                               &noOfBasisFuncs_sg,
			       noOfDensitiesRead,
			       nvaluesList) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_load_density_getSizes");
      return -1;
    }
  
  for(int i = 0; i < *noOfDensitiesRead; i++)
    {
      long nvalues = nvaluesList[i];
      rowindList[i] = new int[nvalues];
      colindList[i] = new int[nvalues];
      valuesList[i] = new ergo_real[nvalues];
    }  
  
  if(ddf_read_shells_and_density_matrices_sparse(basisInfo,
						 *noOfDensitiesRead,
						 rowindList,
						 colindList,
						 valuesList,
						 densityFileName) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in ddf_read_shells_and_density_matrices_sparse");
      return -1;
    }

  if((*basisInfo)->get_basis_funcs() != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in get_basis_funcs");
      return -1;
    }
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "get_basis_funcs for starting guess basis set returned OK,"
	    " number of basis funcs: %i",
            (*basisInfo)->noOfBasisFuncs);
  if((*basisInfo)->getSimplePrimitivesAll(integralInfo) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in getSimplePrimitivesAll");
      return -1;
    }
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, 
	    "getSimplePrimitivesAll for starting guess basis set "
            "returned OK, n = %i",
            (*basisInfo)->noOfSimplePrimitives);
  return 0;
}
