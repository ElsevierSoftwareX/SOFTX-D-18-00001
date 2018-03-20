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

/** @file stream_test.cc Tests the streaming grid generator.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "grid_reader.h"
#include "grid_stream.h"

class MyMolInfo : public  GridGenMolInfo {
public:
  MyMolInfo(int a, int b, int s) : GridGenMolInfo(a, b, s) {}
  void getAtom(int icent, int *cnt, real (*coor)[3],
               int *charge, int *mult) const
  {
    *cnt = 1;
    coor[0][0] = 0;
    coor[0][1] = 0;
    coor[0][2] = icent*100;
    *charge = 1;
    *mult = 1;
  }

  void setShellRadii(real *shellRadii) const
  {
    for(int i=0; i<noOfShells; i++) 
      shellRadii[i] = 5.0;
  }

  void getBlocks(const real *center, real cellsz,
                 const real *rshell,
                 int *nblcnt, int (*iblcks)[2]) const
  {
    *nblcnt = 1;
    iblcks[0][0] = 0;
    iblcks[0][1] = noOfShells;
  }

  void getExps(int *maxl, int **bascnt, real (**aa)[2]) const
  {
    static const int lda = 1;
    *maxl = 1;
    *bascnt = (int*)calloc(lda*noOfAtoms, sizeof(int));
    *aa     = (real(*)[2])calloc(2*lda*noOfShells, sizeof(real));
    for(int i=0; i<noOfAtoms; ++i) {
      (*bascnt)[0 + i*lda] = 1; /* 1 s function */

      /* Range of exponents */
      (*aa)[i][0] = 0.5;
      (*aa)[i][1] = 0.5;
    }
  }
};

int main(int argc, char *argv[])
{
  int nAtoms;

  if(argc<=1 || (nAtoms = atoi(argv[1])) <1)
    nAtoms = 3;
  MyMolInfo mmi(nAtoms,0,3*nAtoms);
  static const Dft::GridParams ggs(1e-11, 9, 35);

  try {
    const char fName[] = "TST.grid";
    ErgoGridStream *egs = grid_stream_new(ggs, mmi);
    const char *str = getenv("OMP_NUM_THREADS");
    grid_stream_generate(egs, fName,
			 str ? strtol(str, NULL, 10) : 1)
      && unlink(fName);
    grid_stream_free(egs);
  } catch(const char *s){
    printf("%s\n", s);
  }
  return 0;
}
