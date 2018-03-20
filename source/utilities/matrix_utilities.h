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

/** @file matrix_utilities.h

    @brief Utilities related to the hierarchical matrix library (HML),
    including functions for setting up permutations of basis functions
    to increase data locality in the hierarchical matrix data
    structure.
*/

#ifndef MATRIX_UTILITIES_HEADER
#define MATRIX_UTILITIES_HEADER

#include "matrix_typedefs.h"
#include "basisinfo.h"

#if 0
/** prepare_matrix_permutation creates a perm object as required by
    further matrix manipulation routines. */
Perm* prepare_matrix_permutation(const BasisInfoStruct& basisInfo,
				 int sparse_block_size,
				 int factor1, int factor2, int factor3);
#else

mat::SizesAndBlocks prepareMatrixSizesAndBlocks(int n_basis_functions,
						int sparse_block_size,
						int factor1, 
						int factor2, 
						int factor3);

void getMatrixPermutation(const BasisInfoStruct& basisInfo,
			  int sparse_block_size,
			  int factor1, 
			  int factor2, 
			  int factor3,
			  std::vector<int> & permutation);
void getMatrixPermutation(const BasisInfoStruct& basisInfo,
			  int sparse_block_size,
			  int factor1, 
			  int factor2, 
			  int factor3,
			  std::vector<int> & permutation,
			  std::vector<int> & inversePermutation);
void getMatrixPermutationOnlyFactor2(const std::vector<ergo_real> & xcoords,
				     const std::vector<ergo_real> & ycoords,
				     const std::vector<ergo_real> & zcoords,
				     int sparse_block_size_lowest,
				     int first_factor, // this factor may be different from 2, all other factors are always 2.
				     std::vector<int> & permutation,
				     std::vector<int> & inversePermutation);
void getMatrixPermutationOnlyFactor2(const BasisInfoStruct& basisInfo,
				     int sparse_block_size_lowest,
				     int first_factor, // this factor may be different from 2, all other factors are always 2.
				     std::vector<int> & permutation,
				     std::vector<int> & inversePermutation);

#endif
void fill_matrix_with_random_numbers(int n, symmMatrix & M);

void add_random_diag_perturbation(int n, 
				  symmMatrix & M, 
				  ergo_real eps);

bool check_if_matrix_contains_strange_elements(const symmMatrix & M,
                                               std::vector<int> const & inversePermutationHML);

void output_matrix(int n, const ergo_real* matrix, const char* matrixName);

template<class Tmatrix>
ergo_real compute_maxabs_sparse(const Tmatrix & M)
{
  return M.maxAbsValue();

}

template<typename RandomAccessIterator>
struct matrix_utilities_CompareClass {
  RandomAccessIterator first;
  explicit matrix_utilities_CompareClass(RandomAccessIterator firstel)
    : first(firstel){}
  bool operator() (int i,int j) { return (*(first + i) < *(first + j));}
};

template<typename Tmatrix>
void get_all_nonzeros( Tmatrix const & A, 
		       std::vector<int> const & inversePermutation,
		       std::vector<int> & rowind,
		       std::vector<int> & colind,
		       std::vector<ergo_real> & values) {
  rowind.resize(0);
  colind.resize(0);
  values.resize(0);
  size_t nvalues = 0;
  size_t nvalues_tmp = A.nvalues();
  std::vector<int>       rowind_tmp; rowind_tmp.reserve(nvalues_tmp);
  std::vector<int>       colind_tmp; colind_tmp.reserve(nvalues_tmp);
  std::vector<ergo_real> values_tmp; values_tmp.reserve(nvalues_tmp);
  A.get_all_values(rowind_tmp,
		   colind_tmp,
		   values_tmp,
		   inversePermutation,
		   inversePermutation);
  // Count the number of nonzeros
  for(size_t i = 0; i < nvalues_tmp; i++) {
    nvalues += ( values_tmp[i] != 0 );
  }
  rowind.reserve(nvalues);
  colind.reserve(nvalues);
  values.reserve(nvalues);
  // Extract all nonzeros
  for(size_t i = 0; i < nvalues_tmp; i++) {
    if ( values_tmp[i] != 0 ) {
      rowind.push_back( rowind_tmp[i] );
      colind.push_back( colind_tmp[i] );
      values.push_back( values_tmp[i] );	
    }
  }
} // end get_all_nonzeros(...)


template<typename Tmatrix>
void write_matrix_in_matrix_market_format( Tmatrix const & A, 
					   std::vector<int> const & inversePermutation,
					   std::string filename,
					   std::string identifier,
					   std::string method_and_basis)
{

  // Get all matrix elements
  size_t nvalues = 0;
  std::vector<int>       rowind;
  std::vector<int>       colind;
  std::vector<ergo_real> values;
  get_all_nonzeros( A, inversePermutation, rowind, colind, values);
  nvalues = values.size();
  // Now we have all matrix elements
  // Open file stream
  std::string mtx_filename = filename + ".mtx";
  std::ofstream os(mtx_filename.c_str());

  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  std::string matrix_market_matrix_type = "general";
  bool matrixIsSymmetric = (A.obj_type_id() == "MatrixSymmetric");
  if (matrixIsSymmetric)
    matrix_market_matrix_type = "symmetric";
  os << "%%MatrixMarket matrix coordinate real " << matrix_market_matrix_type << std::endl
     << "%===============================================================================" << std::endl
     << "% Generated by the Ergo quantum chemistry program version " << VERSION << " (www.ergoscf.org)" << std::endl
     << "% Date      : " << asctime (timeinfo) // newline added by asctime
     << "% ID-string : " << identifier << std::endl
     << "% Method    : " << method_and_basis << std::endl
     << "%" << std::endl
     << "% MatrixMarket file format:" << std::endl
     << "% +-----------------" << std::endl
     << "% | % comments" << std::endl
     << "% | nrows ncols nentries" << std::endl
     << "% | i_1        j_1        A(i_1,j_1)" << std::endl
     << "% | i_2        j_2        A(i_2,j_2)" << std::endl
     << "% | ..." << std::endl
     << "% | i_nentries j_nentries A(i_nentries,j_nentries) " << std::endl
     << "% +----------------" << std::endl
     << "% Note that indices are 1-based, i.e. A(1,1) is the first element." << std::endl
     << "%" << std::endl
     << "%===============================================================================" << std::endl;
  os << A.get_nrows() << "  " << A.get_ncols() << "  " << nvalues << std::endl;
  if (matrixIsSymmetric)
    for(size_t i = 0; i < nvalues; i++) {
      // Output lower triangle
      if ( rowind[i] < colind[i] )
	os << colind[i]+1 << "  " << rowind[i]+1 << "  " << std::setprecision(10) << (double)values[i] << std::endl;
      else
	os << rowind[i]+1 << "  " << colind[i]+1 << "  " << std::setprecision(10) << (double)values[i] << std::endl;
    }      
  else
    for(size_t i = 0; i < nvalues; i++) {
      os << rowind[i]+1 << "  " << colind[i]+1 << "  " << std::setprecision(10) << (double)values[i] << std::endl;
    }  
  os.close();
} // end write_matrix_in_matrix_market_format(...)


#endif
