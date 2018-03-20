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

/** @file xcmat_test.cc Tests the DFT XC matrix construction. 
    This test computes the XC energy many times and checks that the 
    resulting energy is the same every time. If this fails, it is
    probably because of some bug related to synchronization 
    of threads. 
*/

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

static bool
compare_matrices(char mat_name,
                 const real *computed, const long double *ref, int sz,
                 ergo_real eps)
{
  bool failed = false;
  for(int row=0; row<sz; row++) {
    for(int col=0; col<sz; col++) {
      real currdiff = computed[row + col*sz]- ref[row+col*sz];
      if (template_blas_fabs(currdiff)>eps) {
        printf("%c (%d,%d): ref: %28.25Lf got: %28.25Lf diff: %12g eps: %g\n",
               mat_name, row, col,
               (long double)ref[row + col*sz],
               (long double)computed[row + col*sz],
               (double)currdiff,
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
  /* We give some room to accumulation error. */
  /* Since the reference values are computed using long double we
     cannot check more accurately than that. */
  ergo_real EPS = mat::getMachineEpsilon<ergo_real>();
  if(EPS < mat::getMachineEpsilon<long double>())
    EPS = mat::getMachineEpsilon<long double>();
  ergo_real extraFactor = 20;
  if(EPS <= mat::getMachineEpsilon<long double>())
    extraFactor = 230;
  /* ELIAS NOTE 2016-09-06: earlier the factor used here was selected like this:
     (sizeof(ergo_real) == sizeof(ergo_long_real) ? 230 : 20)
     unclear why that would be needed.
     FIXME: find out why long double does not give as high accuracy as
     expected here. */
  EPS *= extraFactor;
  int nrepeat = 2;
  bool failed = false;
  for(int i = 0; i < nrepeat; i++)
    {
      mode = 'R';
      ergo_real dftEnergyAgain = 0, electronsR, electronsU, dftEnergyU;
      memset(xcmat, 0, n*n*sizeof(ergo_real));
      electronsR = dft_get_xc_mt(noOfElectrons, dmat, bis, &m, gridParams,
                                 xcmat, &dftEnergyAgain);
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
  printf("%cXC %s %s/%s test %s\n", failed ? mode : ' ',
         gridSchemeName,
	 selected_func->is_gga() ? "GGA" : "LDA",
	 functional, failed ? "failed" : "OK"); 
  if(!failed)
    unlink("ergoscf.out");
  return  failed ? 1 : 0;
}

static int
test_small_both()
{
  int res = 0;
  IntegralInfo ii(true);

  static const int sys1Z[2] = { 2, 1 };
  static const ergo_real sys1C[2][3] = { { 0, 0, 0 }, { 0, 0, 1.5 } };
  static const int sys2Z[2] = { 2, 1 };
  static const ergo_real sys2C[2][3] = { { 0, 0, 0 }, { 0, 0, 20.0 } };

  /* ELIAS NOTE 2016-09-06: new reference values were computed (using
     gcc 6.1.1) because the test failed after the previous fix adding
     the L suffix for the pi constant in pi.h which is apparently used
     by the XC computations. So the old reference values were wrong
     because they were computed using an inaccurate pi value. Of
     course it's a bad idea to have this kind of test, it just
     verifies that we get the same result as earlier, it does not
     really check if it is correct. But now it is like this, and it is
     still better than no test at all. */
  static const long double XCRefSys1Svwn5_GC2[2][2] = {
    { -0.46846890237096027839L, -0.27741047198177963827L },
    { -0.27741047198177963827L, -0.63052487439961731318L }
  };
  static const long double XCRefSys1Svwn5_Turbo[2][2] = {
    { -0.46837701264321619278L, -0.27730481689558706597L },
    { -0.27730481689558706597L, -0.63065863444345805105L }
  };
  static const long double XCRefSys1Svwn5_LMG[2][2] = {
    { -0.46847089939399400585L, -0.27732916610776639027L },
    { -0.27732916610776639027L, -0.63074970313439271847L }
  };
  static const long double XCRefSys1BP86_LMG[2][2] = {
    { -0.48456321202299733965L, -0.28477354319527885414L },
    { -0.28477354319527885414L, -0.65856888979121373437L }
  };
  static const long double XCRefSys1BP86_TURBO[2][2] = {
    { -0.48447235314731952526L, -0.28476089225530221456L },
    { -0.28476089225530221456L, -0.65847904553389222668L }
  };
  static const long double XCRefSys2Combine[2][2] = {
    { -0.41581586869051081036L,  0.00000000000000000000L },
    {  0.00000000000000000000L, -0.58371746046633447966L }
  };
  static const long double XCRefSys2Blyp[2][2] = {
    { -0.43551596812630800643L,  0.00000000000000000000L },
    {  0.00000000000000000000L, -0.61589810422545549677L }
  };

  res += test_small(ii, "SVWN5", Dft::GridParams::GC2, "GC2  ",
                    sys1Z, sys1C, &XCRefSys1Svwn5_GC2[0]);

  res += test_small(ii, "SVWN5", Dft::GridParams::TURBO, "Turbo",
                    sys1Z, sys1C, &XCRefSys1Svwn5_Turbo[0]);

  res += test_small(ii, "SVWN5", Dft::GridParams::LMG, "LMG  ",
                    sys1Z, sys1C, &XCRefSys1Svwn5_LMG[0]);

  res += test_small(ii, "BP86", Dft::GridParams::LMG, "LMG  ",
                    sys1Z, sys1C, &XCRefSys1BP86_LMG[0]);

  res += test_small(ii, "BP86", Dft::GridParams::TURBO, "Turbo",
                    sys1Z, sys1C, &XCRefSys1BP86_TURBO[0]);

  res += test_small(ii, "Combine Slater=1 PZ81=1",
                    Dft::GridParams::LMG, "LMG  ",
		    sys2Z, sys2C, &XCRefSys2Combine[0]);

  res += test_small(ii, "BLYP", Dft::GridParams::LMG, "LMG  ",
                    sys2Z, sys2C, &XCRefSys2Blyp[0]);
  return res;
}

static int
test_mol(const char *mol_fname, const char *basisSet, const char *xcFunc)
{
  unlink("ergoscf.out");
  dft_init();

  Molecule m;

  if(m.setFromMoleculeFile(mol_fname, 0, NULL) != 0) {
    printf("Molecule::setFromMoleculeFile failed.\n");
    return 1;
  }

  IntegralInfo integralInfo(true);
  BasisInfoStruct *bis = new BasisInfoStruct();

  if(bis->addBasisfuncsForMolecule(m, basisSet, 0, NULL,
                                   integralInfo, 0, 1, 0) != 0) {
    printf("bis->addBasisfuncsForMolecule failed.\n");
    delete bis;
    return 1;
  }

  int n = bis->noOfBasisFuncs;

  /* Set up density matrix. Equivalent to use_simple_starting_guess. */
  ergo_real *dmat= ergo_new(n*n, ergo_real);
  int noOfElectrons = m.getNumberOfElectrons();
  real diag = noOfElectrons/ergo_real(n);
  for(int col=0; col<n; col++) {
    for(int row=0; row<n; row++)
      dmat[row+n*col] = 0.0;
    dmat[col+n*col] = diag;
  }

  if(dft_setfunc(xcFunc) == 0) {
    fprintf(stderr, "Error in dft_setfunc(%s)\n", xcFunc);
    return 1;
  }
  grid_set_tmpdir("/tmp");
  static const int ANGMIN = 6;
  static const int ANGINT = 35;
  static const real RADINT = 1e-7;

  Dft::GridParams gridParams(RADINT, ANGMIN, ANGINT);

  ergo_real *xcmat= ergo_new(n*n, ergo_real);

  ergo_real dftEnergy = 0;
  ergo_real integratedNoOfElectrons = 
    dft_get_xc_mt(noOfElectrons, dmat, bis, &m, gridParams, xcmat, &dftEnergy);

  ergo_free(dmat);
  ergo_free(xcmat);
  grid_free_files();
  delete bis;
  printf("%s/%s benchmark executed. "
         "Expected %d electrons. Integrated %f\n",
         selected_func->is_gga() ? "GGA" : "LDA",
	 xcFunc, noOfElectrons, (double)integratedNoOfElectrons); 
  return  0;
}

int main(int argc, char *argv[])
{
  const char *DEFAULT_XC_FUNC = "BP86";
  if(argc<=1)
    return test_small_both();
  else if(argc==2) 
    return test_mol(argv[1],  ERGO_SPREFIX "/basis/Turbomole-SVP",
		    DEFAULT_XC_FUNC);
  else if(argc==3)
    return test_mol(argv[1],  argv[2], DEFAULT_XC_FUNC);
  else if(argc==4)
    return test_mol(argv[1],  argv[2], argv[3]);
  else {
    fputs("Usage: xcmat_test [MOL_FILE [BASIS_SET]].\n", stderr);
    return 1;
  }
  /* Not reached */
  return 0;
}
