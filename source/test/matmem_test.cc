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

/** @file matmem_test.cc Tests matrix library memory usage. */

#include <stdio.h>
#include <unistd.h>
#include <memory>
#include <limits>
#include <vector>
#include <sstream>

#include "matrix_utilities.h"
#include "integral_matrix_wrappers.h"
#include "output.h"


int main(int argc, char *argv[])
{
  enable_memory_usage_output();

#ifdef _OPENMP
  int defThreads;
  const char *env = getenv("OMP_NUM_THREADS");
  if ( !(env && (defThreads=atoi(env)) > 0) ) {
    defThreads = 1;
  }
  mat::Params::setNProcs(defThreads);
  mat::Params::setMatrixParallelLevel(2);
  std::cout<<"OpenMP is used, number of threads set to "
	   <<mat::Params::getNProcs()<<". Matrix parallel level: "
	   <<mat::Params::getMatrixParallelLevel()<<"."<<std::endl;
#endif


  output_current_memory_usage(LOG_AREA_MAIN, "Beginning");

  int sizeOfVectorObject = sizeof(std::vector<char>);
  int sizeOfPointer = sizeof(char*);
  int sizeOfSizet = sizeof(size_t);
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "sizeOfVectorObject = %d, sizeOfPointer = %d, sizeOfSizet = %d", sizeOfVectorObject, sizeOfPointer, sizeOfSizet);

  for(int k = 0; k < 6; k++) {
    do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "k = %2d", k);
    {
      const int Nvectors = 1200000;
      std::vector< char* > vectorList(Nvectors);
      int totSize = 0;
      for(int i = 0; i < Nvectors; i++) {
	//for(int i = Nvectors-1; i >= 0; i--) {
	int size = (i % 2) * 500 + 12;
	vectorList[i] = new char[size];
	totSize += size;
      }
      do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "totSize = %d = %9.5f GB", totSize, (double)totSize / 1000000000);
      output_current_memory_usage(LOG_AREA_MAIN, "After allocating vectors");
      //for(int i = 0; i < Nvectors; i++)
      for(int i = Nvectors-1; i >= 0; i--)
	delete [] vectorList[i];    
      output_current_memory_usage(LOG_AREA_MAIN, "After deleting vectors");
    }
    output_current_memory_usage(LOG_AREA_MAIN, "After vectors scope ended");
#if 1
    {
      do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "Allocating single vector of 350 MB.");
      std::vector<char> tmpVect(350000000);
      output_current_memory_usage(LOG_AREA_MAIN, "After allocating single vector of 350 MB");
    }
#endif
  }



  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "Doing vector allocation stuff again!");

  {
    const int Nvectors = 3000000;
    std::vector< char* > vectorList(Nvectors);
    int totSize = 0;
    for(int i = 0; i < Nvectors; i++) {
      //for(int i = Nvectors-1; i >= 0; i--) {
      int size = (i % 3) * 250 + 10;
      vectorList[i] = new char[size];
      totSize += size;
    }
    do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "totSize = %d = %9.5f GB", totSize, (double)totSize / 1000000000);
    output_current_memory_usage(LOG_AREA_MAIN, "After allocating vectors");
    //for(int i = 0; i < Nvectors; i++)
    for(int i = Nvectors-1; i >= 0; i--)
      delete [] vectorList[i];    
    output_current_memory_usage(LOG_AREA_MAIN, "After deleting vectors");
  }
  output_current_memory_usage(LOG_AREA_MAIN, "After vectors scope ended");








  
  IntegralInfo integralInfo(true);
  BasisInfoStruct bis;

  Molecule m;

  int nx, ny, nz;
  nx = 10;
  ny = 5;
  nz = 7;


  const ergo_real space = 8.0;
  int atomCount = 0;
  for(int ix = 0; ix < nx; ix++)
    for(int iy = 0; iy < ny; iy++)
      for(int iz = 0; iz < nz; iz++)
	{
	  ergo_real x = ix*space + 0.4*std::cos((ix+iy+iz)*0.2+0.0)*space;
	  ergo_real y = iy*space + 0.4*std::cos((ix+iy+iz)*0.2+0.3)*space;
	  ergo_real z = iz*space + 0.4*std::cos((ix+iy+iz)*0.2+0.6)*space;
	  /* Use a mix of charges: H, C, Zn.
	     It is good to have some Zn there so we check also usage
	     of basis functions of f type. */
	  int charge = 1;
	  if(atomCount%3 == 0)
	    charge = 6;
	  if(atomCount%9 == 0)
	    charge = 30;
	  m.addAtom(charge, x, y, z);
	  atomCount++;
	}
	
  if(bis.addBasisfuncsForMolecule(m, ERGO_SPREFIX "/basis/6-31Gss",
                                   0, NULL, integralInfo, 0, 0, 0) != 0) {
    printf("bis.addBasisfuncsForMolecule failed.\n");
    return 1;
  }
  
  mat::SizesAndBlocks matrix_size_block_info =
    prepareMatrixSizesAndBlocks(bis.noOfBasisFuncs,
				20, 8, 8, 8);
  std::vector<int> permutationHML(bis.noOfBasisFuncs);
  getMatrixPermutation(bis, 20, 8, 8, 8, permutationHML);


  symmMatrix Aempty;
  Aempty.resetSizesAndBlocks(matrix_size_block_info,
			     matrix_size_block_info);

  
  symmMatrix Adiag;
  Adiag.resetSizesAndBlocks(matrix_size_block_info,
			    matrix_size_block_info);
  {
    const int n = bis.noOfBasisFuncs;
    const int nvalues = n;
    std::vector<int> idxrow(nvalues);
    std::vector<int> idxcol(nvalues);
    std::vector<ergo_real> values(nvalues);
    int count = 0;
    for(int i=0; i<n; i++) {
      idxrow[count] = i;
      idxcol[count] = i;
      values[count] = 1.0;
      count++;
    }
    idxrow.resize(count);
    idxcol.resize(count);
    values.resize(count);
    output_current_memory_usage(LOG_AREA_MAIN, "Before Adiag.add_values");
    Adiag.add_values(idxrow, idxcol, values);
    output_current_memory_usage(LOG_AREA_MAIN, "After Adiag.add_values");
    ergo_real memUsage = (double)Adiag.memory_usage();
    ergo_real memUsageGB = memUsage / 1000000000;
    do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "n = %d, Adiag mem usage: %12.0f = %9.5f GB", n, memUsage, memUsageGB);
  }
  output_current_memory_usage(LOG_AREA_MAIN, "After Adiag.add_values scope ended");
  {
    symmMatrix AdiagCopy;
    AdiagCopy = Adiag;
  }
  output_current_memory_usage(LOG_AREA_MAIN, "After AdiagCopy scope ended");



  output_current_memory_usage(LOG_AREA_MAIN, "Before creating A");  

  symmMatrix A;
  A.resetSizesAndBlocks(matrix_size_block_info,
			matrix_size_block_info);
  
  {
    const int n = bis.noOfBasisFuncs;
    const int nvalues = n*n;
    std::vector<int> idxrow(nvalues);
    std::vector<int> idxcol(nvalues);
    std::vector<ergo_real> values(nvalues);
    int count = 0;
    for(int i=0; i<n; i++)
      for(int j=i; j<n; j++) {
	idxrow[count] = i;
	idxcol[count] = j;
	values[count] = 1.0;
	count++;
      }
    idxrow.resize(count);
    idxcol.resize(count);
    values.resize(count);
    output_current_memory_usage(LOG_AREA_MAIN, "Before A.add_values");
#if 0
    A.add_values(idxrow, idxcol, values, 
		 permutationHML,
		 permutationHML);
#else
    A.add_values(idxrow, idxcol, values);
#endif
    output_current_memory_usage(LOG_AREA_MAIN, "After A.add_values");
    ergo_real memUsage = (double)A.memory_usage();
    ergo_real memUsageGB = memUsage / 1000000000;
    do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "n = %d, A mem usage: %12.0f = %9.5f GB", n, memUsage, memUsageGB);
  }
  output_current_memory_usage(LOG_AREA_MAIN, "After A.add_values scope ended");

  {
    output_current_memory_usage(LOG_AREA_MAIN, "Before creating AcopyList");
    symmMatrix AcopyList[88];
    output_current_memory_usage(LOG_AREA_MAIN, "After creating AcopyList");

    {
      std::vector<char> tmpVector(100000000);
      output_current_memory_usage(LOG_AREA_MAIN, "After creating tmpVector of 100000000 bytes");
    }
    output_current_memory_usage(LOG_AREA_MAIN, "After tmpVector scope ended");
    
    for(int i = 0; i < 1; i++)
      {
	std::stringstream ss1;
	ss1 << "Before creating Acopy " << i;
	output_current_memory_usage(LOG_AREA_MAIN, ss1.str().c_str());
	{
	  std::vector<char> tmpVector(100000000);
	  output_current_memory_usage(LOG_AREA_MAIN, "After creating tmpVector of 100000000 bytes");
	  AcopyList[i] = Aempty;
	}
	output_current_memory_usage(LOG_AREA_MAIN, "After AcopyList[i] = Aempty scope ended");
	AcopyList[i] = A;
	std::stringstream ss2;
	ss2 << "After creating Acopy " << i;
	output_current_memory_usage(LOG_AREA_MAIN, ss2.str().c_str());
      }
  }
  output_current_memory_usage(LOG_AREA_MAIN, "After Acopy scope ended");


  printf("Matrix mem usage test OK\n");
  return 0;
}
