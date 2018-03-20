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

/** @file file_tools/test/test.cc

    @brief File containg tests for reading/writing sparse matrices from/to mtx files and dense matrices from/to binary files.

    @author Anastasia Kruchinina <em>responsible</em>
*/


#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include "files_dense.h"
#include "files_sparse.h"
#include "files_sparse_bin.h"

using namespace std;


int main(int argc, char *argv[])
{
   // Create test symmetric sparse matrix
   int  nnz   = 10;
   int  N     = 10;
   int  I[]   = { 1,4,2,3,1,5,6,6,2,6};
   int  J[]   = { 4,4,6,6,8,8,8,9,10,10};
   real val[] = { 2.3505,-1.0616,0.4882,-0.7648,-1.5998,-1.4023,0.7481,-0.6156,0.8886,-0.1924};

   vector<int>  Iv(I, I + sizeof(I) / sizeof(int));
   vector<int>  Jv(J, J + sizeof(J) / sizeof(int));
   vector<real> valv(val, val + sizeof(val) / sizeof(real));


   {
      std::cout << "Test mtx files" << '\n';

      vector<int>  Iv2;
      vector<int>  Jv2;
      vector<real> valv2;
      int          N2;

      if (write_matrix_to_mtx("test_sparse.mtx", Iv, Jv, valv, N) == -1)
      {
         printf("Error in write_matrix_to_mtx.\n");
         return -1;
      }
      read_matrix_from_mtx("test_sparse.mtx", Iv2, Jv2, valv2, N2, N2);

      assert(N == N2);
      assert((int)Iv2.size() == nnz);

      assert(Iv == Iv2);
      assert(Jv == Jv2);
      assert(valv == valv2);
   }

////////////////////////////

{
   std::cout << "Test binary files" << '\n';
   
   vector<int>  Iv2;
   vector<int>  Jv2;
   vector<real> valv2;
   int          N2;
   
   try
   {
      write_matrix_to_bin("test_sparse.bin", Iv, Jv, valv, N);
      read_matrix_from_bin("test_sparse.bin", Iv2, Jv2, valv2, N2, N2);
   }
   catch (const std::exception& e)
   {
     std::cout << "Error:" << e.what() << '\n';
   }
   
   assert(N == N2);
   assert((int)Iv2.size() == nnz);

   assert(Iv == Iv2);
   assert(Jv == Jv2);
   assert(valv == valv2);

}

   cout << "All tests finished OK!" << endl;

   return 1;
}
