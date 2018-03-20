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

/** @file xcmat_nan_inf_test.cc Tests that the DFT XC matrix
    construction does not result in "nan" or "inf" values. */

#include <stdio.h>
#include <unistd.h>
#include <memory>
#include <limits>

#include "integrals_1el_potential.h"
#include "integrals_2el.h"
#include "memorymanag.h"
#include "grid_reader.h"
#include "dft_common.h"
#include "xc_matrix.h"

static bool value_seems_like_nan_or_inf(ergo_real x) {
  bool ok1 = false;
  if(x > -template_blas_get_num_limit_max<ergo_real>())
    ok1 = true;
  bool ok2 = false;
  if(x < template_blas_get_num_limit_max<ergo_real>())
    ok2 = true;
  return ( ! (ok1 && ok2) );
}

static bool
compare_matrices(char mat_name,
                 const real *computed, const long double *ref, int sz,
                 ergo_real eps)
{
  bool failed = false;

  for(int row=0; row<sz; row++) {
    for(int col=0; col<sz; col++) {
      ergo_real theDiff = computed[row + col*sz]- ref[row+col*sz];
      bool nanOrInfFlag = value_seems_like_nan_or_inf(theDiff);
      if(nanOrInfFlag) {
        printf("Error! nan/inf found in compare_matrices().\n");
        failed = true;
      }
      if (template_blas_fabs(theDiff)>eps) {
        printf("%c (%d,%d): ref: %28.25Lf got: %28.25Lf diff: %12g eps: %g\n",
               mat_name, row, col,
               (long double)ref[row + col*sz],
               (long double)computed[row + col*sz],
               (double)(computed[row + col*sz]- ref[row+col*sz]),
               (double)eps);
        failed = true;
      }
    }
  }
  return failed;
}

static int
test_small(const IntegralInfo& ii, const char *functional,
           const Dft::GridParams::RadialScheme& gridScheme,
           const char *gridSchemeName,
	   const int *charges, const real (*coords)[3],
           const long double (*XCRef)[2])
{
  // Skip tests for some functionals for single precision because they cannot handle it due to overflow or similar problems.
#ifdef PRECISION_SINGLE
  std::vector<std::string> functionalsToSkip;
  functionalsToSkip.push_back("PW91X");
  functionalsToSkip.push_back("B3LYP");
  functionalsToSkip.push_back("B3LYP-G");
  functionalsToSkip.push_back("B3P86");
  functionalsToSkip.push_back("B3P86-G");
  functionalsToSkip.push_back("B3PW91");
  functionalsToSkip.push_back("BLYP");
  functionalsToSkip.push_back("BP86");
  functionalsToSkip.push_back("BPW91");
  functionalsToSkip.push_back("Camb3lyp");
  functionalsToSkip.push_back("HSE");
  functionalsToSkip.push_back("KT3");
  functionalsToSkip.push_back("OLYP");
  functionalsToSkip.push_back("PBE0");
  functionalsToSkip.push_back("PBE");
  bool skip_for_single_precision = false;
  for(unsigned int i = 0; i < functionalsToSkip.size(); i++) {
    if(strcmp(functional, functionalsToSkip[i].c_str()) == 0)
      skip_for_single_precision = true;
  }
  if(skip_for_single_precision) {
    printf("test_small(): skipping functional '%s' for PRECISION_SINGLE case.\n", functional);
    return 0;
  }
#endif

  BasisInfoStruct* bis = new BasisInfoStruct();
  Molecule m;
  /* The code later will change the order of atoms, this is why the
     reference table may seem strange at the first sight. */
  for(int i=0; i<2; i++) {
    m.addAtom(charges[i], coords[i][0],coords[i][1],coords[i][2]);
  }

  if(bis->addBasisfuncsForMolecule(m, ERGO_SPREFIX "/basis/STO-3G",
                                   0, NULL, ii, 0, 0, 0) != 0) {
    printf("bis->addBasisfuncsForMolecule failed.\n");
    return 1;
  }

  int n = bis->noOfBasisFuncs;

  /* set up density matrix */
  ergo_real *dmat= ergo_new(n*n, ergo_real);
  dmat[0*n+0] = 1.1; dmat[0*n+1] = 0.2;
  dmat[1*n+0] = 0.2; dmat[1*n+1] = 1.3;

  dft_init();
  if(dft_setfunc(functional) == 0)
    {
      printf("error in dft_setfunc\n");
      return 1;
    }
  grid_set_tmpdir("/tmp");
  static const ergo_real GRID_CELL_SIZE = 2.5;
  Dft::GridParams gridParams(1e-5, 6, 7, GRID_CELL_SIZE);
  gridParams.radialGridScheme = gridScheme;

  ergo_real *xcmat= ergo_new(n*n, ergo_real);
  ergo_real *xca = ergo_new(n*n, ergo_real);
  ergo_real *xcb = ergo_new(n*n, ergo_real);
  ergo_real *dmata = ergo_new(n*n, ergo_real);
  for(int i=n*n-1; i>=0; --i) dmata[i] = 0.5*dmat[i];

  int noOfElectrons = 2;
  char mode;
  ergo_real dftEnergy = 0;
  dft_get_xc_mt(noOfElectrons, dmat, bis, &m, gridParams, xcmat, &dftEnergy);
  if(value_seems_like_nan_or_inf(dftEnergy)) {
    printf("Error: dftEnergy from dft_get_xc_mt() is nan/inf for functional '%s'.\n", functional);
    return 1;
  }
  /* We give some room to accumulation error. */
  static const ergo_real EPS = 0.3;
  int nrepeat = 2;
  bool failed = false;
  for(int i = 0; i < nrepeat; i++)
    {
      mode = 'R';
      ergo_real dftEnergyAgain = 0, electronsR, electronsU, dftEnergyU;
      memset(xcmat, 0, n*n*sizeof(ergo_real));
      electronsR = dft_get_xc_mt(noOfElectrons, dmat, bis, &m, gridParams,
                                 xcmat, &dftEnergyAgain);
      if(value_seems_like_nan_or_inf(dftEnergyAgain)) {
	printf("Error: dftEnergyAgain from dft_get_xc_mt() is nan/inf for functional '%s'.\n", functional);
	return 1;
      }
      failed = compare_matrices('R', xcmat, &XCRef[0][0], n, EPS);
      if(template_blas_fabs(dftEnergyAgain - dftEnergy) > EPS)
	{
	  printf("%s/%s energy repeatability test failed.\n",
		 selected_func->is_gga() ? "GGA" : "LDA", functional);
	  printf("i = %5i of %5i: computed: %20.19f diff: %g\n", 
		 i, nrepeat,
                 (double)dftEnergyAgain, (double)(dftEnergy-dftEnergyAgain));
          failed = true;
	}
      if(failed)
	break;

      mode = 'U';
      memset(xca, 0, n*n*sizeof(ergo_real));
      memset(xcb, 0, n*n*sizeof(ergo_real));
      electronsU = dft_get_uxc_mt(noOfElectrons,
                                  dmata, dmata,
                                  bis, &m, gridParams,
                                  xca, xcb, &dftEnergyU);
      if(value_seems_like_nan_or_inf(dftEnergyU)) {
	printf("Error: dftEnergyU from dft_get_uxc_mt() is nan/inf for functional '%s'.\n", functional);
	return 1;
      }
      failed = compare_matrices('A', xca, &XCRef[0][0], n, EPS)
        || compare_matrices('B', xcb, &XCRef[0][0], n, EPS);
      if (template_blas_fabs(electronsU - electronsR) > EPS) {
          printf("%s/%s Electrons restricted %28.25Lg unrestricted %28.25Lg\n",
                 selected_func->is_gga() ? "GGA" : "LDA", functional,
                 (long double)electronsR,
                 (long double)electronsU);
      }   
      if(failed)
        break;      
    }

  ergo_free(dmat);
  ergo_free(dmata);
  ergo_free(xcmat);
  ergo_free(xca);
  ergo_free(xcb);
  grid_free_files();
  delete bis;
  printf("%cXC %s %s/%s simple inf/nan test %s\n", failed ? mode : ' ',
         gridSchemeName,
	 selected_func->is_gga() ? "GGA" : "LDA",
	 functional, failed ? "failed" : "OK"); 
  if(!failed)
    unlink("ergoscf.out");
  return  failed ? 1 : 0;
}

static int test_functional(const IntegralInfo & ii, const char* funcName) {
  static const int sys1Z[2] = { 2, 1 };
  static const ergo_real sys1C[2][3] = { { 0, 0, 0 }, { 0, 0, 1.5 } };

  static const long double XCRefSys1BP86_TURBO[2][2] = {
    { -0.4844723531473195618241717L, -0.2847608922553022361067940L },
    { -0.2847608922553022361067940L, -0.6584790455338922763911698L }
  };
  return test_small(ii, funcName, Dft::GridParams::TURBO, "Turbo", 
                    sys1Z, sys1C, &XCRefSys1BP86_TURBO[0]);
}

static int
test_small_many()
{
  int res = 0;
  IntegralInfo ii(true);

  //  res += test_functional(*ii, "Becke");
  //  res += test_functional(*ii, "KT");
  res += test_functional(ii, "LB94");
  //  res += test_functional(*ii, "LYP");
  //  res += test_functional(*ii, "OPTX");
  //  res += test_functional(*ii, "P86c");
  //  res += test_functional(*ii, "PW86x");

  res += test_functional(ii, "PW91X");
  //  res += test_functional(*ii, "PW91c");
  //  res += test_functional(*ii, "PW92c");
  //  res += test_functional(*ii, "PZ81");
  //  res += test_functional(*ii, "PBEC");
  res += test_functional(ii, "Pbex");
  res += test_functional(ii, "Slater");

  //  res += test_functional(*ii, "SVWNI");
  //  res += test_functional(*ii, "SVWN3I");
  //  res += test_functional(*ii, "SVWN");
  //  res += test_functional(*ii, "XAlpha");

  res += test_functional(ii, "B3LYP");
  res += test_functional(ii, "B3LYP-G");
  res += test_functional(ii, "B3P86");
  res += test_functional(ii, "B3P86-G");
  res += test_functional(ii, "B3PW91");

  //  res += test_functional(*ii, "BHandH");
  //  res += test_functional(*ii, "BHandHLYP");

  res += test_functional(ii, "BLYP");
  res += test_functional(ii, "BP86");
  res += test_functional(ii, "BPW91");
  res += test_functional(ii, "Camb3lyp");
  //  res += test_functional(*ii, "Cam");

  res += test_functional(ii, "HSE");
  res += test_functional(ii, "KT1");
  res += test_functional(ii, "KT2");
  res += test_functional(ii, "KT3");
  res += test_functional(ii, "LDA");

  res += test_functional(ii, "OLYP");
  res += test_functional(ii, "PBE0");
  res += test_functional(ii, "PBE");
  res += test_functional(ii, "SVWN3");
  res += test_functional(ii, "SVWN5");

  return res;
}

int main(int argc, char *argv[])
{
  int result = test_small_many();
  if(result == 0)
    printf("xcmat_nan_inf_test finished OK\n");
  else
    printf("Error in xcmat_nan_inf_test: result = %d\n", result);
  return result;
}
