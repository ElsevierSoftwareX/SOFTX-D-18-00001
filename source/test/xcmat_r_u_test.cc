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

/** @file xcmat_r_u_test.cc Tests the DFT XC matrix
    construction. Calls both restricted and unrestricted versions of
    the XC matrix construction code, and verifies that the end results
    are close enough. */

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
                 const real *computed, const real *ref, int sz,
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
	   const int *charges, const real (*coords)[3])
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
  printf("Calling dft_setfunc() for functional '%s'\n", functional);
  if(dft_setfunc(functional) == 0) {
    printf("error in dft_setfunc\n");
    return 1;
  }
  grid_set_tmpdir("/tmp");
  static const ergo_real GRID_CELL_SIZE = 2.2;
  Dft::GridParams gridParams(2e-5, 6, 7, GRID_CELL_SIZE);
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
  if(EPS <= mat::getMachineEpsilon<long double>()) extraFactor = 230; // FIXME: is this really needed?
  EPS *= extraFactor;
  ergo_real tol_nel = EPS*10; // tolerance for electron count comparison
  int nrepeat = 1;
  bool failed = false;
  for(int i = 0; i < nrepeat; i++) {
    // Restricted case
    mode = 'R';
    ergo_real dftEnergyAgain = 0, electronsR, electronsU, dftEnergyU;
    memset(xcmat, 0, n*n*sizeof(ergo_real));
    electronsR = dft_get_xc_mt(noOfElectrons, dmat, bis, &m, gridParams,
			       xcmat, &dftEnergyAgain);
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
    // Unrestricted case
    mode = 'U';
    memset(xca, 0, n*n*sizeof(ergo_real));
    memset(xcb, 0, n*n*sizeof(ergo_real));
    electronsU = dft_get_uxc_mt(noOfElectrons,
				dmata, dmata,
				bis, &m, gridParams,
				xca, xcb, &dftEnergyU);
    failed = compare_matrices('A', xca, xcmat, n, EPS)
      || compare_matrices('B', xcb, xcmat, n, EPS);
    if (template_blas_fabs(electronsU - electronsR) > tol_nel) {
      printf("%s/%s Electrons restricted %28.25Lg unrestricted %28.25Lg ( diff %g  tol %g )\n",
	     selected_func->is_gga() ? "GGA" : "LDA", functional,
	     (long double)electronsR,
	     (long double)electronsU,
	     (double)template_blas_fabs(electronsU - electronsR),
	     (double)tol_nel);
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

  res += test_small(ii, "SVWN5", Dft::GridParams::GC2, "GC2  ", sys1Z, sys1C);
  res += test_small(ii, "SVWN5", Dft::GridParams::TURBO, "Turbo", sys1Z, sys1C);
  res += test_small(ii, "SVWN5", Dft::GridParams::LMG, "LMG  ", sys1Z, sys1C);
  res += test_small(ii, "BP86", Dft::GridParams::LMG, "LMG  ", sys1Z, sys1C);
  res += test_small(ii, "BP86", Dft::GridParams::TURBO, "Turbo", sys1Z, sys1C);
  res += test_small(ii, "Combine Slater=1 PZ81=1", Dft::GridParams::LMG, "LMG  ", sys2Z, sys2C);
  res += test_small(ii, "BLYP", Dft::GridParams::LMG, "LMG  ", sys2Z, sys2C);

  return res;
}

int main(int argc, char *argv[])
{
  int result = test_small_both();
  if(result != 0)
    printf("xcmat_r_u_test failed.\n");
  else
    printf("xcmat_r_u_test finished OK.\n");
  return result;
}
