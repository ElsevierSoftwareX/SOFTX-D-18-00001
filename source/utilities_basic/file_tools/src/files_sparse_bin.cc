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

/** @file files_sparse_bin.cc

    @brief File containing definitions of functions for reading/writing sparse matrices from/to binary files.

    @author Anastasia Kruchinina <em>responsible</em>
*/

#include "files_sparse_bin.h"




// I, J, val CONTAIN ROWS, COLUMNS AND VALUES OF THE UPPER TRIANGLE OF THE MATRIX
void write_matrix_to_bin(const char* filename, const std::vector<int> &I, const std::vector<int> &J, const std::vector<real> &val, const int &N)
{
  assert(I.size() == J.size());
  assert(I.size() == val.size());                   
  size_t NNZ = I.size();

  char * buffer;
  // compute expected file size (NNZ, N, N, I, J, val)
  size_t filesize = sizeof(size_t) + 2*sizeof(int) + sizeof(int)*NNZ*2 + 
                    sizeof(real)*NNZ;
    

  // open binary file for writing
  fstream file (filename, ios::out | ios::binary);

  if(!file)
  {
      throw std::runtime_error("write_matrix_to_bin: error writing matrix to the binary file.");
  }

  buffer = new char[filesize];
  size_t pos = 0;
  memcpy(&buffer[pos], &NNZ, sizeof(size_t));
  pos += sizeof(size_t);
  memcpy(&buffer[pos], &N, sizeof(int));
  pos += sizeof(int);  
  memcpy(&buffer[pos], &N, sizeof(int));
  pos += sizeof(int);  
  memcpy(&buffer[pos], &I[0], sizeof(int)*I.size());
  pos += sizeof(int)*I.size();  
  memcpy(&buffer[pos], &J[0], sizeof(int)*J.size());
  pos += sizeof(int)*J.size();  
  memcpy(&buffer[pos], &val[0], sizeof(real)*val.size());

  // write data to the file
  file.write(buffer, filesize);
  delete[] buffer;
  
  // get length of file:
  file.seekg (0, file.end);
  size_t length = file.tellg();
  file.seekg (0, file.beg);

  if(filesize != length)
   throw std::runtime_error("write_matrix_to_bin: filesize != length");

  // close the file
  file.close();

}

/**
   Read data from the binary file. Matrix is stored in the format:
   row column value row column value row column value ...
 */
void read_matrix_from_bin_Elias_format(const char* filename, std::vector<int> &I, std::vector<int> &J, std::vector<real> &val)
{
  // open binary file for reading                                                                   
  fstream file (filename, ios::in | ios::binary);
  
  if(!file)
    {
      throw std::runtime_error("read_matrix_from_bin: error reading matrix from the binary file.");
    }

  // get length of file:                                                                                    
  file.seekg (0, file.end);
  size_t filesize = file.tellg();
  file.seekg (0, file.beg);

  char *buffer = new char[filesize];

  // read data from the file                                                                                
  file.read(buffer, filesize);

  size_t NNZ = filesize/16; // 16 = 4+4+8
  size_t pos = 0;

  I.resize(NNZ);
  J.resize(NNZ);
  val.resize(NNZ);  

  for(unsigned long int i = 0; i < NNZ; ++i)
  {
    memcpy(&I[i], &buffer[pos], sizeof(int));
    pos += sizeof(int); 
    memcpy(&J[i], &buffer[pos], sizeof(int));
    pos += sizeof(int);
    memcpy(&val[i], &buffer[pos], sizeof(double));
    pos += sizeof(double);
  }

  file.close();
  delete[] buffer;

}




/**
   Read data from the binary file. Matrix is stored in the format:
   NNZ N M [rows] [columns] [values]
 */
void read_matrix_from_bin(const char* filename, std::vector<int> &I, std::vector<int> &J, std::vector<real> &val, int &N, int &M)
{
  // open binary file for reading
  fstream file (filename, ios::in | ios::binary);

  if(!file)
  {
      throw std::runtime_error("read_matrix_from_bin: error reading matrix from the binary file.");
  }

  // get length of file:
  file.seekg (0, file.end);
  size_t filesize = file.tellg();
  file.seekg (0, file.beg);
  
  char *buffer = new char[filesize];
    
  // read data from the file
  file.read(buffer, filesize);
  
  size_t NNZ;
  
  size_t pos = 0;
  memcpy(&NNZ, &buffer[pos], sizeof(size_t));
  pos += sizeof(size_t);
  
  size_t expected_filesize =  sizeof(size_t) + 2*sizeof(int) + sizeof(int)*NNZ*2 + 
                              sizeof(real)*NNZ;
  if(filesize != expected_filesize)
    throw std::runtime_error("write_matrix_to_bin: filesize != expected_filesize");

  memcpy(&N, &buffer[pos], sizeof(int));
  pos += sizeof(int);  
  memcpy(&M, &buffer[pos], sizeof(int));
  pos += sizeof(int);  
  
  I.resize(NNZ);
  J.resize(NNZ);
  val.resize(NNZ);
  
  // copy matrix data
  memcpy(&I[0], &buffer[pos], sizeof(int)*I.size());
  pos += sizeof(int)*I.size();  
  memcpy(&J[0], &buffer[pos], sizeof(int)*J.size());
  pos += sizeof(int)*J.size();  
  memcpy(&val[0], &buffer[pos], sizeof(real)*val.size());
  
  // for (size_t i = 0; i < val.size(); i++) {
  //   std::cout << val[i] << " ";
  // }
  //   std::cout << '\n';
  
  file.close();
  delete[] buffer;
}

