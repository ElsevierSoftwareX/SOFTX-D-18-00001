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

/** @file explicit_integral_test.cc Tests the explicit computation of
    2-electron integrals by moving basis functions by small distances
    and verifying that the computed 2-el integrals vary smoothly. */

#include <stdio.h>
#include <unistd.h>
#include <memory>
#include <limits>
#include <vector>

#include "integrals_2el_explicit.h"
#include "matInclude.h"

int try_integral_diffs(const IntegralInfo & integralInfo, ergo_real h) {

  Molecule m;

  int nx = 2;
  int ny = 1;
  int nz = 1;

  const ergo_real space = 8.8;
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
	  //	  if(atomCount%9 == 0)
	  //	    charge = 30;
	  m.addAtom(charge, x, y, z);
	  atomCount++;
	}

  ergo_real moleculeBaseCoord = m.getAtom(0).coords[0];

  ergo_real machine_epsilon = mat::getMachineEpsilon<ergo_real>();

  // OK, now we have a molecule. Now do a loop where we move one atom
  // by a small distance and save resulting 2-el integrals for each
  // case.
  const int noOfCases = 55;
  std::vector< std::vector< std::vector< std::vector< std::vector<ergo_real> > > > > integralList(noOfCases);
  int n = 0; // Will be set later.
  for(int caseIdx = 0; caseIdx < noOfCases; caseIdx++) {
    // Move first atom by small distance in x direction.
    Atom atom = m.getAtom(0);
    atom.coords[0] = moleculeBaseCoord + caseIdx * h;
    m.replaceAtom(0, atom);
    // New scope here so that the BasisInfoStruct pointer is freed each time.
    {
      BasisInfoStruct bis;
      if(bis.addBasisfuncsForMolecule(m, ERGO_SPREFIX "/basis/3-21G",
				       0, NULL, integralInfo, 0, 0, 0) != 0) {
	printf("bis.addBasisfuncsForMolecule failed.\n");
	return -1;
      }
      // OK, now we have a basis set for the current version of the
      // molecule. Now compute integrals.
      n = bis.noOfBasisFuncs;
      integralList[caseIdx].resize(n);
      for(int i = 0; i < n; i++) {
	integralList[caseIdx][i].resize(n);
	for(int j = i; j < n; j++) {
	  integralList[caseIdx][i][j].resize(n);
	  for(int k = 0; k < n; k++) {
	    integralList[caseIdx][i][j][k].resize(n);
	    for(int l = k; l < n; l++) {
	      integralList[caseIdx][i][j][k][l] = do_2e_integral(i, k, l, j, bis, integralInfo);
	    }
	  }
	}
      }
    } // End of scope for temporary BasisInfoStruct
  } // END FOR caseIdx
  // OK, now we have integrals computed for all cases. Now check if
  // each integral value varies smoothly.
  int count = 0;
  for(int i = 0; i < n; i++)
    for(int j = i; j < n; j++)
      for(int k = 0; k < n; k++)
	for(int l = k; l < n; l++) {
	  // Now we are interested in the integral (ij|kl). First
	  // check if this integral values seems to be approximately
	  // linearly dependent on the modufied atomic coordinate.
	  ergo_real totDiff  = integralList[noOfCases-1][i][j][k][l] - integralList[0][i][j][k][l];
	  ergo_real halfDiff = integralList[noOfCases/2][i][j][k][l] - integralList[0][i][j][k][l];
	  if(halfDiff != 0 && template_blas_fabs(totDiff) > machine_epsilon) {
	    ergo_real kvot = totDiff / halfDiff;
	    // Now if kvot is approximately 2 we think this integral behaves linearly.
	    if(kvot > 1.95 && kvot < 2.05) {
	      count++;
	      for(int caseIdx = 1; caseIdx < noOfCases; caseIdx++) {
		ergo_real diff = integralList[caseIdx][i][j][k][l] - integralList[caseIdx-1][i][j][k][l];
		ergo_real expectedDiff = totDiff / (noOfCases-1);
		if(diff/expectedDiff < 0.9 || diff/expectedDiff > 1.1) {
		  printf("Error for integral i j k l : %d %d %d %d, expectedDiff = %9.5g, diff = %9.5g\n", 
			 i, j, k, l, (double)expectedDiff, (double)diff);
		  return -1;
		} // END IF
	      } // END FOR caseIdx
	    } // END IF
	  } // END IF (halfDiff != 0)
	} // END FOR i j k l

  unlink("ergoscf.out");

  printf("try_integral_diffs finished OK, count = %d.\n", count);

  if(count == 0)
    exit(0);

  return 0;
}


int main(int argc, char *argv[]) {
  IntegralInfo integralInfo(true);
  integralInfo.init();

  ergo_real machine_epsilon = mat::getMachineEpsilon<ergo_real>();
  ergo_real h = 10000*machine_epsilon;
  printf("machine_epsilon = %9.5g, using h = %9.5g\n", (double)machine_epsilon, (double)h);

  if(try_integral_diffs(integralInfo, h) != 0)
    return -1;

  return 0;
}
