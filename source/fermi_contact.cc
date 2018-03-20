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

/** @file fermi_contact.cc implements "Fermi contact" integrals. The
    formulas are give in the hyperfine project.

    @author: Pawel Salek <em>responsible</em>
*/

#include <memory>

static const double GE = 2.0023193044;

#include "memorymanag.h"
#include "basisinfo.h"
#include "aos.h"
#include "density_description_file.h"
#include "Matrix.h"
#include "MatrixSymmetric.h"

struct FCAccumulator {
  const ergo_real *basisFuncValues; /**< vector of basis function values at
                                     * given point in space. */
  explicit FCAccumulator(const ergo_real *bfs) : basisFuncValues(bfs) {}
  ergo_real accumulate(const ergo_real& dij, int i, int j) const
  {
    return dij*basisFuncValues[i]*basisFuncValues[j];
  }
};

template<class Accumulator>
ergo_real
accumulate(int n, const ergo_real *spinMat, const Accumulator& ac)
{
  ergo_real res = 0;
  for(int col=0; col<n; col++)
    for(int row=0; row<n; row++) {
      ergo_real c = ac.accumulate(spinMat[row + col*n], row, col);
      res += c;
    }
  return res;
}

/** computeFermiContact computes the Fermi contact interaction for
    given molecule and specified spin density.

    @param bis basis set specification.
    @param spinDensity spin density defined as D_alpha-D_beta.

    @param R the cartesian coordinates to which we compute interaction
    to.
    @param result will contain the interaction if the function succeeds.

    @return 0 on success, -1 on failure.
*/
int
computeFermiContact(const BasisInfoStruct& bis,
                    const ergo_real* spinDensity,
                    const Vector3D& R, ergo_real& result)
{
  ergo_real *bfs = new ergo_real[bis.noOfBasisFuncs];

  memset(bfs, 0, sizeof(bfs[0])*bis.noOfBasisFuncs);

  int bfBlock[2] = { 0, bis.noOfShells };

  dft_get_orbs(1, bfs, &R.v, 1, &bfBlock, 0, bis);

  result = accumulate(bis.noOfBasisFuncs, spinDensity, FCAccumulator(bfs));

  delete []bfs;
  return 0;
}

#if 0
static void
printmat(const char*label, int n, const ergo_real *m)
{
  for(int row=0; row<n; row++) {
    for(int col=0; col<n; col++)
      printf("%14f ", (double)m[row + col*n]);
    puts("");
  }
}
#endif

int main(int argc, char *argv[])
{
  static char usage[] =
    "fermi_contact SPIN_DENSITY_FILE [MOLECULE]\n"
    "\tComputes Fermi contact term for given unrestricted density file.\n"
    "\tThe first file is assumed to contain the alpha and beta densities.\n"
    "\tThe optional second file contains molecular geometry.\n"
    "\n"
    "Spin couplings to all basis function centers are computed\n"
    "if only one file is specified. If the molecule file is specified\n"
    "as well, the labels from there are used to tag the atoms.\n";

  if(argc <= 1) {
    fputs(usage, stderr);
    return 1;
  }

  IntegralInfo integralInfo(true);
  
  ergo_real       *densMatrix[2] = { NULL, NULL };
  BasisInfoStruct *basisRead  = NULL;
  if (ddf_load_density(argv[1], 2, integralInfo,
                       &basisRead, densMatrix)) {
    fprintf(stderr, "Loading unrestricted densities from '%s' failed. "
            "Calculation aborted.\n", argv[1]);
    return 1;
  }

#if 0
  printmat("ALPHA", basisRead->noOfBasisFuncs, densMatrix[0]);
  printmat("BETA",  basisRead->noOfBasisFuncs, densMatrix[1]);
#endif

  /* Create spin density */
  int nElements = basisRead->noOfBasisFuncs*basisRead->noOfBasisFuncs;
  for(int i=0; i<nElements; i++)
    densMatrix[0][i] -= densMatrix[1][i];

  if (argc>=3) {
    Molecule molecule;
    printf("Loading molecule from %s\n", argv[2]);
    int res = molecule.setFromMoleculeFile
      (argv[2], 0,  /* we are guessing the net charge here */
       NULL);
    if(res) {
      fprintf(stderr,
              "Molecule file '%s' specified but could not be loaded.\n",
              argv[2]);
      return 1;
    }
    printf("%-6s %6s %12s %12s      : %5s\n", "Charge",
           "X", "Y", "Z", "FC Coupling");
    for(int iAtom=0; iAtom<molecule.getNoOfAtoms(); iAtom++) {
      const Atom &atom = molecule.getAtom(iAtom);
      ergo_real fc;
      Vector3D pos(atom.coords[0], atom.coords[1], atom.coords[2]);
      if (computeFermiContact(*basisRead, densMatrix[0], pos, fc) ) {
        fprintf(stderr, "Calculation of FC SS failed for atom %d\n",
                iAtom+1);
        break;
      } else {
        printf("%-6.1f %12.6f %12.6f %12.6f : %17.10g\n", (double)atom.charge,
               (double)pos[0], (double)pos[1], (double)pos[2], 
	       (double)(4*M_PI*GE*fc/3.0));
      }
    }
  } else { /* No molecule file specified - use basis function centers */
    printf("FC Couplings\n"
           "%6s %12s %12s       : %s\n",
           "X", "Y", "Z", "FC Coupling");
    Vector3D lastPos(-12345e6, -12345e6, -12345e6);
    for(int iShell=0; iShell<basisRead->noOfShells; iShell++) {
      const ShellSpecStruct &shell = basisRead->shellList[iShell];
      ergo_real fc;
      Vector3D pos(shell.centerCoords[0], shell.centerCoords[1],
                   shell.centerCoords[2]);
      if(pos.dist(lastPos) >0.1) {
        if (computeFermiContact(*basisRead, densMatrix[0], pos, fc) ) {
          fprintf(stderr, "Calculation of FC SS failed for shell %d\n",
                  iShell+1);
          break;
        } else {
          printf("%12.6f %12.6f %12.6f : %17.10g\n",
                 (double)pos[0], (double)pos[1], (double)pos[2], 
		 (double)(4*M_PI*GE*fc/3.0));
        }
      }
      lastPos = pos;
    }
  }
  ergo_free(densMatrix[0]);
  ergo_free(densMatrix[1]);
  delete basisRead;

  return 0;
}
