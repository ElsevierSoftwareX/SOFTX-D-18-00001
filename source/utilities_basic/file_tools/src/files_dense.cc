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

/** @file files_dense.cc

    @brief File containing definition of functions for reading/writing dense matrices and vectors.

    @author Anastasia Kruchinina <em>responsible</em>
*/


#include "files_dense.h"

using namespace std;

typedef ergo_real real;
#define PRECISION 16



/* READ DENSE MATRIX FROM THE FILE (TEXT OR BINARY, ALL MATRIX ELEMENTS ARE STORED) */
int read_matrix(const char *filename, vector<real> &A, int &N, int &M, bool is_binary)
{  
  ifstream f;
  if(!is_binary)
    f.open (filename, ios::in);
  else
    f.open (filename, ios::in | std::ios::binary);

  if (!f.is_open())
    {
      printf("Error: cannot open file\n");
      return -1;
    }

  // Make sure to avoid int overflow when computing N*M
  size_t N_x_M = (size_t)N * (size_t)M;

  if(!is_binary)
  {
    // read size
    f >> N;
    f >> M;

    A.clear();
    A.resize(N_x_M);

  }


  if(!is_binary)
    { 
      for(size_t i = 0; i < N_x_M; ++i ) {
	long double tmp;
	f >> tmp;
	A[i] = tmp;
      }
    }
  else
    {
    // get length of file:
      f.seekg (0, f.end);
      size_t length = f.tellg(); // number of bytes
      f.seekg (0, f.beg);

      N = std::sqrt(length/8);
      M = N;

      A.clear();
      A.resize(N_x_M);

      std::cout << "Reading " << length << " characters   ==   " << N_x_M << " double numbers... ";
      // read data as a block:
      f.read((char *) &A[0], A.size()*sizeof(real)); 
      if (f)
	std::cout << "all characters read successfully.";
      else
	std::cout << "error: only " << f.gcount() << " could be read";
    }

  f.close();

  return 1;
}



/* READ VECTOR FROM THE FILE (TEXT OR BINARY)*/
int read_vector(const char *filename, vector<real> &A, int &N, bool is_binary)
{  
  ifstream f;
  if(!is_binary)
    f.open (filename, ios::in);
  else
    f.open (filename, ios::in | std::ios::binary);

  if (!f.is_open())
    {
      printf("Error: cannot open file\n");
      return -1;
    }

  if(!is_binary)
  {
    // read size
    f >> N;

    A.clear();
    A.resize(N);
  }

  if(!is_binary)
    { 
      for(int i = 0; i < N; ++i ) {
	long double tmp;
	f >> tmp;
	A[i] = tmp;
      }
    }
  else
    {
    // get length of file:
      f.seekg (0, f.end);
      size_t length = f.tellg(); // number of bytes
      f.seekg (0, f.beg);

      N = std::sqrt(length/8);

      A.clear();
      A.resize(N);

      std::cout << "Reading " << length << " characters   ==   " << N << " double numbers... ";
      // read data as a block:
      f.read((char *) &A[0], A.size()*sizeof(real)); 
      if (f)
	std::cout << "all characters read successfully.";
      else
	std::cout << "error: only " << f.gcount() << " could be read";
    }

  f.close();

  return 1;
}







/* WRITE DENSE MATRIX TO THE FILE (NOT BINARY, WRITE ELEMENTWISE) */
int write_matrix(const char *filename, const vector<real> &A, int N, int M)
{
  ofstream f;
  f.open (filename, ios::out);
  if (!f.is_open())
    {
      printf("Error: cannot open file\n");
      return -1;
    }


  f.precision(PRECISION);

  f << N << "  " << M << endl;
  
  for(int i = 0; i < N; ++i )
    {
      for(int j = 0; j < M; ++j)
	f << (double)A[i*M+j] << " ";
      f << endl;
    }
  
  return 1;
}


/* WRITE VECTOR TO THE FILE (NOT BINARY) */
int write_vector(const char *filename, const vector<real> &v)
{
  ofstream f;
  f.open (filename, ios::out);
  if (!f.is_open())
    {
      printf("Error: cannot open file\n");
      return -1;
    }


  f.precision(PRECISION);

  int N = v.size();

  for(int i = 0; i < N; ++i )
    f << (double)v[i] << " ";
  
  return 1;
}

