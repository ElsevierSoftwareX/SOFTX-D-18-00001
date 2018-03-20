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

/** @file grid_test.cc Tests the DFT grid generation.  This
    test generates the grid, possibly several times to detect problems
    with eg. thread synchronisation.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <memory>

#include "dft_common.h"
#include "grid_reader.h"
#include "grid_stream.h"

#define NATOMS      2
#define N_BF_SHELLS 2
#define N_BFS       4

static const struct {
  ergo_real position[3];
  int charge;
} Atoms[] = {
  { { 0, 0, 0 }, 1 },
  { { 0, 0, 1 }, 2 }
};

static ergo_real ShellRadii[] = { 1, 2 };

class MyMolInfo : public  GridGenMolInfo {
public:
  MyMolInfo() : GridGenMolInfo(NATOMS, N_BFS, N_BF_SHELLS) {}
  virtual void getAtom(int icent, int *cnt, real (*coor)[3],
                       int *charge, int *mult) const
  {
    *cnt = 1;
    coor[0][0] = Atoms[icent].position[0];
    coor[0][1] = Atoms[icent].position[1];
    coor[0][2] = Atoms[icent].position[2];
    *charge    = Atoms[icent].charge;
    *mult = 1;
  }

  virtual void setShellRadii(real *shellRadii) const
  {
    memcpy(shellRadii, ShellRadii, N_BF_SHELLS*sizeof(ergo_real));
  }

  virtual void getBlocks(const real *center, real cellsz, const real *rshell,
                         int *nblcnt, int (*iblcks)[2]) const
  {
    *nblcnt = 1;
    if (center[0]*center[0] + center[1]*center[1] + center[2]*center[2]
        < 3) {
      iblcks[0][0] = 0;
      iblcks[0][1] = N_BF_SHELLS;
    } else {
      iblcks[0][0] = 1;
      iblcks[0][1] = 2;
    }
  }

  virtual void getExps(int *maxl, int **bascnt, real (**aa)[2]) const
  {
    static const int lda = 2;
    maxl[0] = 2;
    *bascnt = (int*)calloc(lda*NATOMS, sizeof(int));
    *aa     = (real(*)[2])calloc(2*lda*NATOMS, sizeof(real));
    (*bascnt)[0 + 0*lda] = 1; /* s functions on first atom */
    (*bascnt)[1 + 0*lda] = 0; /* p functions on the first atom */
    (*bascnt)[0 + 1*lda] = 0; /* s functions on first atom */
    (*bascnt)[1 + 1*lda] = 1; /* p functions on the first atom */
  
    /* Range of exponents on the first atom */
    (*aa)[0][0] = 0.5;
    (*aa)[0][1] = 0.5;
    (*aa)[1][0] = 0.0; /* no p functions */
    (*aa)[1][1] = 0.0; /* no p functions */

    /* Range of exponents on the second atom */
    (*aa)[2][0] = 0.0;
    (*aa)[2][1] = 0.0;
    (*aa)[3][0] = 0.2;
    (*aa)[3][1] = 0.2;
  }
};

static const MyMolInfo MolInfo;

static bool
pattern_to_ps(Dft::SparsePattern& p, const char* fName)
{
//  const char * PS_PREFIX = "";
#if 0
    "%%!PS-Adobe-2.0 EPSF-2.0\n"
    "%%%%Title: Test\n"
    "%%%%Creator: Pawel Salek\n"
    "%%%%Magnification: 1.00\n"
    "%%%%Orientation: Portrait\n"
    "%%%%EndComments\n\n"
    "/inch { 72 mul } def\n"
    "/b { 0 rmoveto currentpoint currentpoint newpath moveto\n"
    "1 0 rlineto 0 1 rlineto -1 0 rlineto closepath fill moveto } def\n"
    "%% now run it!\n"
    "1 inch 1 inch moveto\n";
#endif

//  const char * PS_SUFFIX =
//    "\nshowpage\n"
//    "%%EOF\n";

  FILE *f=fopen(fName, "wt");
  if (!f) {
    perror("fopen");
    return false;
  }
//  fprintf(f, PS_PREFIX);
  fprintf(f, "%f %f scale %% size = %d\n", 72*6.0/p.size(), 72*6.0/p.size(),
	  p.size());
  for(int column=0; column<p.size(); ++column) {
    int lastbox = 0;
    fprintf(f,"currentpoint\n");
    for (Dft::SparsePattern::Column::Iterator i = p[column].begin();
	 i != p[column].end();
	 ++i) {
      fprintf(f,"%d b\n", *i-lastbox);
      //for(int x=lastbox; x<(*i)-1; ++x) putchar(' '); putchar('x');
      lastbox = *i;

    }
    fprintf(f,"moveto 0 1 rmoveto\n");
    //puts("");
  }
//  fprintf(f, PS_SUFFIX);
  fclose(f);
  return true;
}
/** This routine tests the sparsity pattern generation scalability
    properties.  At some point in time, water boxes of increasing
    sizes had irregular sparse pattern. This test helps to debug such
    cases.
 */ 
static bool
grid_test_scaling(const char *fName)
{
  time_t tm; time(&tm);

  Molecule mol;
  int res = mol.setFromMoleculeFile(fName, 0, NULL);

  if(res != 0) {
    fprintf(stderr, "Loading molecule from %s failed.\n", fName);
    return false;
  }

  IntegralInfo ii(true);
  BasisInfoStruct bisOrig;
  if(bisOrig.addBasisfuncsForMolecule(mol, ERGO_SPREFIX "/basis/STO-1G",
				       0, NULL, ii, 0, 0, 0) != 0) {
    printf("bis->addBasisfuncsForMolecule failed.\n");
    throw "addBasisfuncs failed";
  }
  Dft::GridParams gridParams(1e-7, 6, 25);
  gridParams.radialGridScheme = Dft::GridParams::LMG;
  
  int *shellMap = new int[bisOrig.noOfShells];
  int *aoMap = new int[bisOrig.noOfBasisFuncs];

  Dft::setupShellMap(bisOrig, shellMap, aoMap);
  BasisInfoStruct * bisPermuted = bisOrig.permuteShells(shellMap, ii);
  delete []shellMap;
  delete []aoMap;

  ErgoMolInfo molInfo(*bisPermuted, mol);
  /* The important structure of this test */
  Dft::SparsePattern pattern(*bisPermuted); 

  char grid_file_name[] = "scaling_test_XXXXXX";
  close(mkstemp(grid_file_name));
  ErgoGridStream *egStream = grid_stream_new(gridParams, molInfo);
  grid_stream_set_sparse_pattern(egStream, &pattern);
  grid_stream_generate(egStream, grid_file_name, dft_get_num_threads());
  grid_stream_free(egStream);

  printf("Stop %lu s wall time %d nonzero elements\n",
	 ((unsigned long)time(NULL))-tm, pattern.sizeTotal());
  unlink(grid_file_name);

  pattern_to_ps(pattern, "pattern.ps");
  return true;
}

static void
grid_test_synchronisation()
{
  static const int NREPEAT = 500;
  static const int BATCH_LENGTH = 50000;

  // Elias note: earlier, coor and weight were placed directly on the stack like this:
  //    real coor[BATCH_LENGTH][3];
  //    real weight[BATCH_LENGTH];
  // but that gave stack overflow on some systems, so it was changed
  // to the dynamic allocation below.

  real (*coor)[3];
  std::vector<real> coorData(3*BATCH_LENGTH);
  coor = (real (*)[3])(&coorData[0]);
  std::vector<real> weight(BATCH_LENGTH);

  Dft::GridParams gridParams(1e-5, 6, 7);

  Molecule mol;
  for(int i = 0; i < NATOMS; i++)
    mol.addAtom(Atoms[i].charge, Atoms[i].position[0], Atoms[i].position[1], Atoms[i].position[2]);

  IntegralInfo ii(true);
  BasisInfoStruct bis;
  if(bis.addBasisfuncsForMolecule(mol, ERGO_SPREFIX "/basis/STO-1G",
				   0, NULL, ii, 0, 0, 0) != 0) {
    printf("bis.addBasisfuncsForMolecule failed.\n");
    throw "addBasisfuncs failed";
  }
  
  for(int i = 0; i < NREPEAT; i++) {
    DftGridReader* g = grid_open_full(&MolInfo, gridParams, NULL, NULL, bis);
    int np;
    unsigned cnt = 0;
    while ( (np = grid_getchunk_blocked(g, BATCH_LENGTH,
                                        NULL, NULL, 
                                        coor, &weight[0])) >= 0){
      cnt += np;
    }
    grid_close(g);
    grid_free_files();
  }
  printf("Grid generation %i times succeeded.\n", NREPEAT); 
}

int
main(int argc, char *argv[])
{  
  const char *tmpdir = getenv("TMPDIR");
  tmpdir = tmpdir ? tmpdir : "/tmp";
  grid_set_tmpdir(tmpdir);

  switch (argc) {
  case 2:
    grid_test_scaling(argv[1]);
    break;
  default:
    grid_test_synchronisation();
  }

  return 0;
}
