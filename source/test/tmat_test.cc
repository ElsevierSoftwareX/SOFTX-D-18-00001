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

/** @file tmat_test.cc Tests the kinetic energy matrix
    construction. The purpose of the test in its current form is
    mostly to verify compilation correctness. */

#include <stdio.h>
#include <unistd.h>
#include <memory>
#include <limits>

#include "integrals_1el_kinetic.h"
#include "memorymanag.h"
#include "matInclude.h"

int main(int argc, char *argv[])
{
  IntegralInfo integralInfo(true);
  BasisInfoStruct bis;
  Molecule m;
  static const long double T[2][2] = {
    {   0.760031854530854163621168L, 0.383253665253394198238115L }, 
    {   0.383253665253394198238115L, 0.760031854530854163621168L }
  };
  int verbose = getenv("VERBOSE") != NULL;

  m.addAtom(1, 0,0,0);
  m.addAtom(1, 0,0,1);

  if(bis.addBasisfuncsForMolecule(m, ERGO_SPREFIX "/basis/STO-3G",
                                   0, NULL, integralInfo, 0, 0, 0) != 0) {
    printf("bis.addBasisfuncsForMolecule failed.\n");
    return 1;
  }
  int n = bis.noOfBasisFuncs;
  ergo_real *mat= ergo_new(n*n, ergo_real);
  ergo_real EPS = mat::getMachineEpsilon<ergo_real>()*10;
  /* Ugly fix because the reference values are only accurate ato
     double precision so we cannot compare long double to anything
     more accurate than that. */
  ergo_real EPS_double = mat::getMachineEpsilon<double>()*10;
  if(EPS < EPS_double)
    EPS = EPS_double;
  if (compute_T_matrix_full(bis, EPS, mat)) {
    printf("compute_T_matrix_full failed.\n");
    return 1;
  }
  int failed = 0;
  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) 
      if ( template_blas_fabs(mat[j+i*n]-T[j][i]) > EPS) {
        printf("(%d %d): reference: %27.24Lf got: %27.24Lf diff: %12g\n", j, i,
	       static_cast<long double>(T[j][i]),
               static_cast<long double>(mat[j+i*n]),
               static_cast<double>(T[j][i]-mat[j+i*n]));
	failed++;
      } else {
        if(verbose)
          printf("%d %d OK\n", i, j);
      }
  }
  ergo_free(mat);
  if (!failed) { puts("T test succeeded."); unlink("ergoscf.out"); }
  return failed;
}
