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

/** @file xcmat_sparse_test.cc Tests the sparse XC matrix construction. 
*/

#include <stdio.h>
#include <unistd.h>
#include <memory>
#include <limits>
#include <vector>

#include "integrals_1el_potential.h"
#include "integrals_2el.h"
#include "memorymanag.h"
#include "dft_common.h"
#include "grid_reader.h"
#include "xc_matrix_sparse.h"
#include "matrix_utilities.h"
#include "config.h" // Needed to get the PRECISION_SINGLE macro

static const bool PRINT_TIME = false;

static void
calculation_shared(const IntegralInfo& ii, const Molecule& mol,
                   const char *funcName, int blSize, int blFactor,
                   symmMatrix& xcMat, ergo_real *energy,
		   std::vector<int> & permutationHML,
		   bool useHiCu)
{
  time_t tm; time(&tm);
  BasisInfoStruct bis;
  if(bis.addBasisfuncsForMolecule(mol, ERGO_SPREFIX "/basis/4-31G",
                                   0, NULL, ii, 0, 0, 0) != 0) {
    printf("bis.addBasisfuncsForMolecule failed.\n");
    throw "addBasisfuncs failed";
  }
  if(dft_setfunc(funcName) == 0)
    {
      printf("error in dft_setfunc\n");
      throw "dft functional setup failed";
    }
  static Dft::GridParams gridParams(1e-7, 6, 35);
  if(useHiCu)
    gridParams.gridType = Dft::GridParams::TYPE_HICU;
  
  int nElectrons = mol.getNumberOfElectrons();

  mat::SizesAndBlocks matrix_size_block_info =
    prepareMatrixSizesAndBlocks(bis.noOfBasisFuncs, blSize,
				blFactor, blFactor, blFactor);
  getMatrixPermutation(bis, blSize,
		       blFactor, blFactor, blFactor, 
		       permutationHML);
  
  symmMatrix dmat;
  dmat.resetSizesAndBlocks(matrix_size_block_info,
				  matrix_size_block_info);
  xcMat.resetSizesAndBlocks(matrix_size_block_info,
				   matrix_size_block_info);

  {
    std::vector<int> idx(bis.noOfBasisFuncs);
    std::vector<ergo_real> values(bis.noOfBasisFuncs);
    for(int i=0; i<bis.noOfBasisFuncs; i++) {
      idx[i] = i;
      values[i] = 1.0;
    }
    dmat.add_values(idx, idx, values, permutationHML, permutationHML);
  }

  Dft::getXC_mt(bis, ii, mol, gridParams, nElectrons, dmat,
		xcMat, energy, permutationHML);

  if(PRINT_TIME)
    printf("Stop %lu s wall time\n", ((unsigned long)time(NULL))-tm);
}

static bool
small_calculation_core(const IntegralInfo& ii,
		       const char *functionalName,
		       const long double (*xcRef)[2], long double xcERef,
		       bool useHiCu)
{
  bool failed = false;
  
  Molecule m;
  /* The code later will change the order of atoms, this is why the
     reference table may seem strange at the first sight. */
  m.addAtom(2, 0,0,0);
  m.addAtom(1, 0,0,1.5);

  static const int BL_SIZE = 4;
  static const int BL_FACTOR = 2;

  symmMatrix xcMat;
  ergo_real dftEnergy;
  std::vector<int> permutationHML;
  calculation_shared(ii, m, functionalName, BL_SIZE, BL_FACTOR,
		     xcMat, &dftEnergy, permutationHML, useHiCu);

  /* We give some room to accumulation error. For long double calculation,
   * accumulation factor for xc matrix element equal to 200 would suffice.
   * Energy comparison need to be looser, this why we take 2100. */
  /* Since the reference values are computed using long double we
     cannot check more accurately than that. */
  ergo_real EPS_accurate = mat::getMachineEpsilon<ergo_real>();
  if(EPS_accurate < mat::getMachineEpsilon<long double>())
    EPS_accurate = mat::getMachineEpsilon<long double>();
  ergo_real extraFactor = 100;
  if(EPS_accurate <= mat::getMachineEpsilon<long double>())
    extraFactor = 2100;
  EPS_accurate *= extraFactor;
#ifdef PRECISION_SINGLE
  static const ergo_real EPS_sloppy = 3e-5;
#else
  static const ergo_real EPS_sloppy = 2e-5;
#endif
  static ergo_real EPS = EPS_accurate;
  /* Allow larger error for HiCu grid. */
  if(useHiCu)
    EPS = EPS_sloppy;

  std::vector<int> rowind(1);
  std::vector<int> colind(1);
  std::vector<ergo_real> values(1);
  for(int row=0; row<2; row++) {
    for(int col=0; col<2; col++) {
      rowind[0] = row;
      colind[0] = col;
      xcMat.get_values(rowind, colind, values, permutationHML, permutationHML);
      if (template_blas_fabs(values[0] - xcRef[row][col])>EPS) {
        printf(" (%d,%d): ref: %28.25Lf got: %28.25Lf diff: %g eps: %g\n",
               row, col,
               static_cast<long double>(xcRef[row][col]),
               static_cast<long double>(values[0]), 
               double(template_blas_fabs(values[0] - xcRef[row][col])),
	       double(EPS));
        failed = true;
      }
    }
  }
  std::string gridStr = "std ";
  if(useHiCu)
    gridStr = "HiCu";
  if(template_blas_fabs(xcERef - dftEnergy) > EPS)
    {
      printf("Sparse XC %s (grid %s) test failed: could not reproduce the same energy.\n"
	     "Computed: %25.22Lf diff: %g eps: %g\n", functionalName, gridStr.c_str(),
             static_cast<long double>(dftEnergy), double(xcERef-dftEnergy),
             double(EPS));
      return false;
    }
  
  if(!failed) {
    printf("Sparse XC %-10s (grid %s) test OK\n", functionalName, gridStr.c_str());
    unlink("ergoscf.out");
  } else {
      printf("Sparse XC %-10s (grid %s) test FAILED\n", functionalName, gridStr.c_str());
  }
  return !failed;
}

static bool
small_calculation(const IntegralInfo& ii)
{
#if 0
  /* these used to work at some point in time */
  static const long double XCRefBLYP[2][2] = {
    { -0.6469105968311400356582017L,  -0.3406940239203784543807908L  },
    { -0.3406940239203784543807908L,  -0.3377037854748100635876931L  }
  };
  static const long double REF_XC_ENERGY_BLYP = -1.3018204657660243451014L;

  static const long double XCRefSVWN5[2][2] = {
    { -0.6219879322708015512524184L,  -0.3367149426990971624717303L  },
    { -0.3367149426990971624717303L,  -0.3428892519010255774368490L  }
  };
  static const long double REF_XC_ENERGY_SVWN5 = -1.2452936316020187427307L;
#else
  /* gcc version 4.5.1 20100924 (Red Hat 4.5.1-4) (GCC) on x86_64
     yields different long double results! */
  static const long double XCRefBLYP[2][2] = {
    { -0.6469105968311401637566883L,  -0.3406940239203784613196847L  },
    { -0.3406940239203784613196847L,  -0.3377037854748099993487144L  }
  };
  static const long double REF_XC_ENERGY_BLYP = -1.3018204657660245257295L;

  static const long double XCRefSVWN5[2][2] = {
    { -0.6219879322708016678583620L,  -0.3367149426990971693564141L  },
    { -0.3367149426990971693564141L,  -0.3428892519010255078039644L  }
  };
  static const long double REF_XC_ENERGY_SVWN5 = -1.2452936316020189159862L;
#endif
  int errors = 0;

  /* Tests with standard grid.  */

  if (!small_calculation_core(ii, "BLYP", XCRefBLYP, REF_XC_ENERGY_BLYP, false))
    errors++;

  if (!small_calculation_core(ii, "SVWN5", XCRefSVWN5, REF_XC_ENERGY_SVWN5, false))
    errors++;

  /* Remove grid files to make sure new grid is generated. */
  grid_free_files();

  /* Tests with HiCu grid.  */

  if (!small_calculation_core(ii, "BLYP", XCRefBLYP, REF_XC_ENERGY_BLYP, true))
    errors++;

  if (!small_calculation_core(ii, "SVWN5", XCRefSVWN5, REF_XC_ENERGY_SVWN5, true))
    errors++;
    
  return errors == 0;
}

static bool
benchmark_calculation(const IntegralInfo& ii, int sideLength)
{
  static const int CHARGE = 2;
  static const double REF_XC_ENERGY 
    = -0.448990406907508*sideLength*sideLength;
  static const double DISTANCE_BETWEEN_ATOMS = 3.5;
  bool failed = false;
  int nElectrons = 0;
  
  Molecule mol;
  for(int i=0; i<sideLength; i++) {
    for(int j=0; j<sideLength; j++)
      mol.addAtom(CHARGE, 0,
                   i*DISTANCE_BETWEEN_ATOMS,j*DISTANCE_BETWEEN_ATOMS);
    nElectrons += CHARGE*sideLength;
  }

  static const int BL_SIZE = 32;
  static const int BL_FACTOR = 8;

  symmMatrix xcMatrix;

  ergo_real dftEnergy;
  std::vector<int> permutationHML;
  calculation_shared(ii, mol, "SVWN5", BL_SIZE, BL_FACTOR,
		     xcMatrix, &dftEnergy, permutationHML, false);
#if 1
  /* We give some room to accumulation error. */
  static const ergo_real EPS = 1e-5;
  if(template_blas_fabs(REF_XC_ENERGY - dftEnergy) > EPS*sideLength*sideLength)
    {
      printf("DFT XC test failed: could not reproduce the same energy.\n");
      printf("Computed: %25.22Lf Reference: %25.2Lf diff: %g\n",
             static_cast<long double>(dftEnergy),
             static_cast<long double>(REF_XC_ENERGY),
             double(dftEnergy-REF_XC_ENERGY));
      return false;
    }
#endif
  //unlink("ergoscf.out");
  return !failed;
}

static bool
mol_calculation(const IntegralInfo& ii, const char *fname)
{
  bool failed = false;
  
  Molecule mol;

  char *basisSetFile = NULL;
  int res = mol.setFromMoleculeFile(fname, 0, &basisSetFile);

  if(res != 0)
    return false;

  static const int BL_SIZE = 32;
  static const int BL_FACTOR = 8;

  symmMatrix xcMatrix;

  ergo_real dftEnergy;
  std::vector<int> permutationHML;
  calculation_shared(ii, mol, "SVWN5", BL_SIZE, BL_FACTOR,
		     xcMatrix, &dftEnergy, permutationHML, false);

  //unlink("ergoscf.out");
  return !failed;
}

int main(int argc, char *argv[])
{
  bool success = false;
  try {
    static const int PROBLEM_SQUARE_SIDE_LENGTH = 60;
    IntegralInfo ii(true);
    dft_init();
    const char *tmpdir = getenv("TMPDIR");
    tmpdir = tmpdir ? tmpdir : "/tmp";
    grid_set_tmpdir(tmpdir);
    if(getenv("RUN_BENCHMARK")) {
      printf("Running an XC benchmark, tmpdir=%s.\n", tmpdir);
      success = benchmark_calculation(ii, PROBLEM_SQUARE_SIDE_LENGTH);
    } else {
      if(argc>1) {
	int side = strtol(argv[1], NULL, 10);
	if(side)
	  success = benchmark_calculation(ii, side);
	else
	  success = mol_calculation(ii, argv[1]);
      } else 
	success = small_calculation(ii);
    }
  }
  catch(const std::exception& e) {
    std::cerr << "Error: exception caught, e.what() = " << e.what() << std::endl;
    return EXIT_FAILURE;
  } 
  catch(char const* e) {
    std::cerr << "Error: (char*) exception caught: " << e << std::endl;
    return EXIT_FAILURE;
  }
  catch(...) {
    printf("Error: excpetion caught.\n");
    return EXIT_FAILURE;
  }
  printf("Success: %s\n", success ? "YES" : "NO");
  return  success ? 0 : 1;
}
