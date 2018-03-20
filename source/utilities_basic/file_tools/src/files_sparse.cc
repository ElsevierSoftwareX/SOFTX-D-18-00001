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

/** @file files_sparse.cc

    @brief File containing definitions of functions for reading/writing sparse matrices from/to mtx (MatrixMarket format) files.

    @author Anastasia Kruchinina <em>responsible</em>
*/


#include "files_sparse.h"

using namespace std;

typedef ergo_real real;




/* READ SPARSE MATRIX FROM THE MATRIX MARKET FILE */
int read_matrix_from_mtx(const char* filename, std::vector<int> &I, vector<int> &J, vector<real> &val, int &N, int &M)
{
  
  int ret_code;
  MM_typecode matcode;
  FILE *f;
  int nz;   
  int i;
  
  if ((f = fopen(filename, "r")) == NULL) 
  {
    printf("Error opening file!\n");
    return -1;
  }
  
  
  if (mm_read_banner(f, &matcode) != 0)
  {
    throw std::runtime_error("Error in read_matrix_from_mtx: could not process Matrix Market banner.\n");
  }
  /*  This is how one can screen matrix types if their application */
  /*  only supports a subset of the Matrix Market data types.      */
  
  if (mm_is_complex(matcode) && mm_is_matrix(matcode))
  {
    std::ostringstream err;
    err << "Error in read_matrix_from_mtx: sorry, this application does not support ";
    err << "Market Market type: [" << mm_typecode_to_str(matcode) << "]  \n";
    throw std::runtime_error(err.str());
  }
  
  if (!mm_is_sparse(matcode))
  {
    throw std::runtime_error("Error in read_matrix_from_mtx: sorry, this application does not support non sparse matrix types");
  }
  
  
  /* find out size of sparse matrix .... */
  
  
  if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
  return -1;
  
  /* reserve memory for matrices */
  
  I.resize(nz);// = (int *) malloc(nz * sizeof(int));
  J.resize(nz);// = (int *) malloc(nz * sizeof(int));
  val.resize(nz);// = (double *) malloc(nz * sizeof(double));
  
  /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
  /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
  /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
  
  // if matrix is symmetric
  if (mm_is_symmetric(matcode))
  {
    // we need the upper triangle
    int Ii, Ji;
    real vali;
    for(i = 0; i < nz; ++i)
    {
      // Use vali_double here to avoid problem with %lg for long double precision.
      double vali_double = 0;
      int nvalues_from_fscanf = fscanf(f, "%d %d %lg\n", &Ii, &Ji, &vali_double);
      vali = vali_double;
      assert(nvalues_from_fscanf == 3);
      // Matrix Market store in the lower triangle, so we transpose it
      I[i] = Ji-1; J[i] = Ii-1; 
      val[i] = vali; 
    }
  }   
  else
  {
    int Ii, Ji;
    real vali;
    for(i = 0; i < nz; ++i)
    {
      // Use vali_double here to avoid problem with %lg for long double precision.
      double vali_double = 0;
      int nvalues_from_fscanf = fscanf(f, "%d %d %lg\n", &Ii, &Ji, &vali_double);
      vali = vali_double;
      assert(nvalues_from_fscanf == 3);
      I[i] = Ii-1; J[i] = Ji-1; 
      val[i] = vali;
    }
  }
  
  
  if (f !=stdin) fclose(f);
  
  // print matrix
  /*
  mm_write_banner(stdout, matcode);
  mm_write_mtx_crd_size(stdout, M, N, nz);
  for (i=0; i<nz; i++)
  fprintf(stdout, "%d %d %20.19g\n", I[i]+1, J[i]+1, val[i]);
  */
  return 1;
  
}


/* WRITE SPARSE MATRIX FROM THE MATRIX MARKET FILE - SYMMETRIC MATRIX*/
int write_matrix_to_mtx(const char* filename, const vector<int> &I, const vector<int> &J, const vector<real> &val, const int &N)
{
  assert(I.size() == J.size());
  assert(I.size() == val.size());
  
  MM_typecode matcode;                         
  size_t NNZ = I.size();
  
  mm_initialize_typecode(&matcode);
  mm_set_matrix(&matcode);
  mm_set_coordinate(&matcode);
  mm_set_real(&matcode);
  mm_set_symmetric(&matcode);
  mm_set_sparse(&matcode);
  
  FILE *f = fopen(filename, "w");
  if (f == NULL)
  {
    throw std::runtime_error("Error in write_matrix_to_mtx: error opening file!\n");
  }
  
  mm_write_banner(f, matcode); 
  mm_write_mtx_crd_size(f, N, N, NNZ);
  
  /* NOTE: matrix market files use 1-based indices, i.e. first element
  of a vector has index 1, not 0.  */
  for(size_t i = 0; i < NNZ; ++i)
  {
    // input is an upper triangle, so transpose the matrix before saving
    // Store just lower triangle
    assert(J[i] >= I[i]);
    fprintf(f, "%d %d %10.16g\n", J[i]+1, I[i]+1, (double)val[i]);
  }
  
  fclose(f);
  
  return 0;
}


/* WRITE SPARSE MATRIX TO THE MATRIX MARKET FILE - UNSYMMETRIC MATRIX*/
int write_matrix_to_mtx_nonsymm(const char* filename, const vector<int> &I, const vector<int> &J, const vector<real> &val, const int &N, const int &M)
{
  assert(I.size() == J.size());
  assert(I.size() == val.size());
  
  MM_typecode matcode;                         
  size_t NNZ = I.size();
  
  mm_initialize_typecode(&matcode);
  mm_set_matrix(&matcode);
  mm_set_coordinate(&matcode);
  mm_set_real(&matcode);
  mm_set_general(&matcode);
  mm_set_sparse(&matcode);
  
  FILE *f = fopen(filename, "w");
  if (f == NULL)
  {
    throw std::runtime_error("Error in write_matrix_to_mtx: error opening file!\n");
  }
  
  mm_write_banner(f, matcode); 
  mm_write_mtx_crd_size(f, N, M, NNZ);
  
  /* NOTE: matrix market files use 1-based indices, i.e. first element
  of a vector has index 1, not 0.  */
  for(size_t i = 0; i < NNZ; ++i)
  fprintf(f, "%d %d %10.16g\n", I[i]+1, J[i]+1, (double)val[i]);
  
  fclose(f);
  
  return 0;
}




















