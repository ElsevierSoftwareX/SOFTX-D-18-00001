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

/** @file matrix_utilities.cc

    @brief Utilities related to the hierarchical matrix library (HML),
    including functions for setting up permutations of basis functions
    to increase data locality in the hierarchical matrix data
    structure.
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include "matrix_utilities.h"
#include "output.h"


mat::SizesAndBlocks prepareMatrixSizesAndBlocks(int n_basis_functions,
						int sparse_block_size,
						int factor1, 
						int factor2, 
						int factor3) {
  int bSizeVecTmp[5];
  bSizeVecTmp[4] = 1; 
  bSizeVecTmp[3] = sparse_block_size;
  bSizeVecTmp[2] = bSizeVecTmp[3] * factor1;
  bSizeVecTmp[1] = bSizeVecTmp[2] * factor2;
  bSizeVecTmp[0] = bSizeVecTmp[1] * factor3;
  std::vector<int> blockSizeVector(bSizeVecTmp, bSizeVecTmp + 5);
  do_output(LOG_CAT_INFO, LOG_AREA_UNDEFINED, "creating matrix SizesAndBlocks using blockSizeVector:");
  for(int i = 0; i < 5; i++)
    do_output(LOG_CAT_INFO, LOG_AREA_UNDEFINED, "blockSizeVector[%i] = %12i", i, blockSizeVector[i]);
  return mat::SizesAndBlocks(blockSizeVector, n_basis_functions);
}


/* ** Permutation help functions 
   Note that the following functions are used to create the inverse 
   permutation compared to the permutation used in the ergo main program.    
*/

template<typename RandomAccessIterator>
struct CompareClass {
  RandomAccessIterator first;
  explicit CompareClass(RandomAccessIterator firstel)
    : first(firstel){}
  bool operator() (int i,int j) { return (*(first + i) < *(first + j));}
};


template<typename Treal, typename TIndexIterator>
void sortCoord(std::vector<Treal> const & xpos, 
	       std::vector<Treal> const & ypos, 
	       std::vector<Treal> const & zpos,
	       TIndexIterator first, 
	       TIndexIterator last) {
  CompareClass<typename std::vector<Treal>::const_iterator> 
    compareX(xpos.begin());
  CompareClass<typename std::vector<Treal>::const_iterator> 
    compareY(ypos.begin());
  CompareClass<typename std::vector<Treal>::const_iterator>  
    compareZ(zpos.begin());
  Treal xmin = xpos[*std::min_element(first, last, compareX)];
  Treal xmax = xpos[*std::max_element(first, last, compareX)];
  Treal ymin = ypos[*std::min_element(first, last, compareY)];
  Treal ymax = ypos[*std::max_element(first, last, compareY)];
  Treal zmin = zpos[*std::min_element(first, last, compareZ)];
  Treal zmax = zpos[*std::max_element(first, last, compareZ)];
  Treal xrange = xmax - xmin;
  Treal yrange = ymax - ymin;
  Treal zrange = zmax - zmin;
  /* Sort in direction with largest range */
  if (xrange>=yrange && xrange>=zrange) 
    std::sort (first, last, compareX);
  else if (yrange>zrange) 
    std::sort (first, last, compareY);
  else 
    std::sort (first, last, compareZ);
}

template<typename Treal>
void permuteAndRecurse(std::vector<Treal> const & xpos, 
		       std::vector<Treal> const & ypos, 
		       std::vector<Treal> const & zpos,
		       std::vector<int> & index,
		       int const first, 
		       int const last,
		       std::vector<int> const & blockSizes,
		       int bSizeIndex) {
  if (last - first > blockSizes[bSizeIndex]) {
    sortCoord(xpos, ypos, zpos, 
	      index.begin() + first, 
	      index.begin() + last);
    int sizeBox1 = 0;
    while (sizeBox1 < (last - first) / 2)
      sizeBox1 += blockSizes[bSizeIndex];
    permuteAndRecurse(xpos, ypos, zpos, 
		      index, first, first + sizeBox1,
		      blockSizes, bSizeIndex);
    permuteAndRecurse(xpos, ypos, zpos, 
		      index, first + sizeBox1, last,
		      blockSizes, bSizeIndex);    
  }
  else {
    ++bSizeIndex;
    if (bSizeIndex < (int)blockSizes.size()) {
      permuteAndRecurse(xpos, ypos, zpos, 
			index, first, last,
			blockSizes, bSizeIndex);
    }
  }
}


/* Added by Anastasia, neeeded for creation of
 large artificial hamiltonian matrix. */
template<typename Treal>
void permuteAndRecurseStart(std::vector<Treal> const & xpos, 
		       std::vector<Treal> const & ypos, 
		       std::vector<Treal> const & zpos,
		       std::vector<int> & index,
		       int const first, 
		       int const last,
		       std::vector<int> const & blockSizes,
		       int bSizeIndex) {

    sortCoord(xpos, ypos, zpos, 
	      index.begin() + first, 
	      index.begin() + last);

		int sizeBox1 = index.size()/7;
		for(int i = 0; i < 7; ++i)
    permuteAndRecurse(xpos, ypos, zpos, 
		      index, first+sizeBox1*i, first + sizeBox1*(i+1),
		      blockSizes, bSizeIndex);   
}


template<typename Treal>
void getPermutation(std::vector<Treal> const & xpos, 
		    std::vector<Treal> const & ypos, 
		    std::vector<Treal> const & zpos,
		    std::vector<int> & permutation,
		    std::vector<int> const & blockSizes) {
  permutation.resize(xpos.size());
  for (unsigned int ind = 0; ind < permutation.size(); ++ind) {
    permutation[ind] = (int)ind;
  }
  permuteAndRecurse(xpos, ypos, zpos, 
		    permutation, 
		    0, permutation.size(), 
		    blockSizes, 0);
}


/* ** End of permutation help functions*/

void getMatrixPermutation(const BasisInfoStruct& basisInfo,
			  int sparse_block_size,
			  int factor1, 
			  int factor2, 
			  int factor3,
			  std::vector<int> & permutation,
			  std::vector<int> & inversePermutation) {
  
  int n = basisInfo.noOfBasisFuncs;
  std::vector<ergo_real> xlong(n);
  std::vector<ergo_real> ylong(n);
  std::vector<ergo_real> zlong(n);
  for(int i = 0; i < n; i++) {
    xlong[i] = basisInfo.basisFuncList[i].centerCoords[0];
    ylong[i] = basisInfo.basisFuncList[i].centerCoords[1];
    zlong[i] = basisInfo.basisFuncList[i].centerCoords[2];
  }

  int bSizeVecTmp[5];
  bSizeVecTmp[4] = 1;
  bSizeVecTmp[3] = sparse_block_size;
  bSizeVecTmp[2] = bSizeVecTmp[3] * factor1;
  bSizeVecTmp[1] = bSizeVecTmp[2] * factor2;
  bSizeVecTmp[0] = bSizeVecTmp[1] * factor3;
  std::vector<int> blockSizeVector(bSizeVecTmp, bSizeVecTmp + 5);
  do_output(LOG_CAT_INFO, LOG_AREA_UNDEFINED, "creating matrix permutation using blockSizeVector:");
  for(int i = 0; i < 5; i++)
    do_output(LOG_CAT_INFO, LOG_AREA_UNDEFINED, "blockSizeVector[%i] = %12i", i, blockSizeVector[i]);
  getPermutation(xlong, ylong, zlong, 
		 inversePermutation,
		 blockSizeVector);    
  permutation.resize(inversePermutation.size());
  for (unsigned int ind = 0; ind < inversePermutation.size(); ++ind)
    permutation[inversePermutation[ind]] = ind;
}

void getMatrixPermutation(const BasisInfoStruct& basisInfo,
			  int sparse_block_size,
			  int factor1, 
			  int factor2, 
			  int factor3,
			  std::vector<int> & permutation) {
  std::vector<int> inversePermutationDummy;
  getMatrixPermutation(basisInfo,
		       sparse_block_size,
		       factor1, 
		       factor2, 
		       factor3,
		       permutation,
		       inversePermutationDummy);
}


void getMatrixPermutationOnlyFactor2(const std::vector<ergo_real> & xcoords,
				     const std::vector<ergo_real> & ycoords,
				     const std::vector<ergo_real> & zcoords,
				     int sparse_block_size_lowest,
				     int first_factor_in, // this factor may be different from 2, all other factors are always 2.
				     std::vector<int> & permutation,
				     std::vector<int> & inversePermutation) {
  int first_factor = first_factor_in;
  // If the given first_factor parameter is 1, then we proceed as if first_factor=2 anyway. This anyway gives us the permutation we want in that case.
  if(first_factor == 1)
    first_factor = 2;
  // Check how many levels are needed.
  int n = xcoords.size();
  int nLevels = 2;
  int nTmp = n / sparse_block_size_lowest;
  bool first = true;
  while(nTmp > 1) {
    int currFactor = 2;
    if(first) {
      currFactor = first_factor;
      first = false;
    }
    nTmp /= currFactor;
    nLevels++;
  }
  std::vector<int> bSizeVecTmp(nLevels);
  bSizeVecTmp[nLevels-1] = 1;
  bSizeVecTmp[nLevels-2] = sparse_block_size_lowest;
  if(nLevels >= 3)
    bSizeVecTmp[nLevels-3] = sparse_block_size_lowest * first_factor;
  for(int i = nLevels-4; i >= 0; i--)
    bSizeVecTmp[i] = bSizeVecTmp[i+1] * 2;
  std::vector<int> blockSizeVector(&bSizeVecTmp[0], &bSizeVecTmp[nLevels]);
  do_output(LOG_CAT_INFO, LOG_AREA_UNDEFINED, "creating matrix permutation using blockSizeVector (%2d levels):", nLevels);
  for(int i = 0; i < nLevels; i++)
    do_output(LOG_CAT_INFO, LOG_AREA_UNDEFINED, "blockSizeVector[%i] = %12i", i, blockSizeVector[i]);
  getPermutation(xcoords, ycoords, zcoords,
		 inversePermutation,
		 blockSizeVector);
  permutation.resize(inversePermutation.size());
  for (unsigned int ind = 0; ind < inversePermutation.size(); ++ind)
    permutation[inversePermutation[ind]] = ind;
}


void getMatrixPermutationOnlyFactor2(const BasisInfoStruct& basisInfo,
				     int sparse_block_size_lowest,
				     int first_factor, // this factor may be different from 2, all other factors are always 2.
				     std::vector<int> & permutation,
				     std::vector<int> & inversePermutation) {
  int n = basisInfo.noOfBasisFuncs;
  std::vector<ergo_real> xcoords(n);
  std::vector<ergo_real> ycoords(n);
  std::vector<ergo_real> zcoords(n);
  for(int i = 0; i < n; i++) {
    xcoords[i] = basisInfo.basisFuncList[i].centerCoords[0];
    ycoords[i] = basisInfo.basisFuncList[i].centerCoords[1];
    zcoords[i] = basisInfo.basisFuncList[i].centerCoords[2];
  }
  getMatrixPermutationOnlyFactor2(xcoords,
				  ycoords,
				  zcoords,
				  sparse_block_size_lowest,
				  first_factor,
				  permutation,
				  inversePermutation);
}


void 
fill_matrix_with_random_numbers(int n, symmMatrix & M)
{
#if 1
  M.random();
#else
  ergo_real* full = new ergo_real[n*n];
  for(int i = 0; i < n; i++)
    for(int j = i; j < n; j++)
      {
	ergo_real a = rand() / (ergo_real)RAND_MAX;
	full[i*n+j] = a;
	full[j*n+i] = a;
      }
  M.assign_from_full(full, n, n);
  delete []full;
#endif
}


static ergo_real rand_minus1_to_1()
{
  ergo_real a = rand() / (ergo_real)RAND_MAX;
  // now a is between 0 and 1
  a *= 2;
  // now a is between 0 and 2
  a -= 1;
  // now a is between -1 and 1
  return a;
}

void 
add_random_diag_perturbation(int n, 
			     symmMatrix & M, 
			     ergo_real eps)
{
  std::vector<ergo_real> randomVector(n);
  std::vector<int> rowIndVector(n);
  std::vector<int> colIndVector(n);
  for(int i = 0; i < n; i++)
    {
      rowIndVector[i] = i;
      colIndVector[i] = i;
      randomVector[i] = eps * rand_minus1_to_1();
    }
  /* No permutation needed since this is a diagonal random element add. */
  M.add_values(rowIndVector, colIndVector, randomVector); 
}


/** This function is supposed to check if a matrix contains any strange numbers such as "inf" or "nan". 
 *  The function returns true is any strange numbers are found, and false if the matrix seems ok. */
bool
check_if_matrix_contains_strange_elements(const symmMatrix & M,
                                          std::vector<int> const & inversePermutationHML)
{
  std::vector<int> rowind;
  std::vector<int> colind;
  std::vector<ergo_real> values;
  M.get_all_values(rowind,
                   colind,
                   values,
                   inversePermutationHML,
                   inversePermutationHML);
  int nvalues = values.size();
  for(int i = 0; i < nvalues; i++) {
    ergo_real x = values[i];
    bool ok1 = false;
    if(x > -template_blas_get_num_limit_max<ergo_real>())
      ok1 = true;
    bool ok2 = false;
    if(x < template_blas_get_num_limit_max<ergo_real>())
      ok2 = true;
    if( ! (ok1 && ok2) )
      return true;
  }
  return false;
}


void 
output_matrix(int n, const ergo_real* matrix, const char* matrixName)
{
  int nn = n;
  printf("output_matrix for matrix '%s', n = %i:\n", matrixName, n);
  if(n > 15) {
    printf("output_matrix showing truncated matrix\n");
    nn = 15;
  }
  for(int i = 0; i < nn; i++)
    {
      for(int j = 0; j < nn; j++)
	printf("%9.4f ", (double)matrix[i*n+j]);
      printf("\n");
    }
}





