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

/** @file bin2m.cc A program for conversion of binary matrix file to a
    matlab format file. It accepts a file name as the input and prints
    the text form on the standard output. The file can be then fed to
    Matlab and will create matrix called "m" which can be then renamed
    within Matlab.
*/

#include <stdio.h>
#include <memory>

#include "basisinfo.h"
#include "density_description_file.h"
#include "memorymanag.h"

int
main(int argc, char *argv[])
{
  static const int COMPRESSED_DIM = 1200;
  static const char usage[] = 
    "Usage: bin2m [-c] density.bin\n"
    "Program will produce on standard output a matlab script\n"
    "loading given density.\n"
    "Option -c compresses the matrix for visualiation purposes\n"
    "to not exceed 1200x1200.\n";
  if(argc<2) {
    fputs(usage, stderr);
    return 1;
  }

  bool compress = false;
  int fname_index = 1;

  if(strcmp(argv[1], "-c") == 0) {
    fname_index++;
    compress = true;
  }

  if(fname_index >=argc) {
    fputs(usage, stderr);
    return 1;
  }

  IntegralInfo integralInfo(true);
  ergo_real *matrix = NULL;
  BasisInfoStruct *basis = NULL;

  if(ddf_load_density(argv[fname_index], 1, integralInfo,
                      &basis, &matrix)) {
    fprintf(stderr, "Loading a matrix from '%s' failed.\n",
            argv[fname_index]);
    return -1;
  }

  printf("m=[\n");
  if(compress && basis->noOfBasisFuncs > COMPRESSED_DIM) {

    for(int row=0; row<COMPRESSED_DIM; row++) {
      int rowLo = (basis->noOfBasisFuncs*row)/COMPRESSED_DIM;
      int rowHi = (basis->noOfBasisFuncs*(row+1))/COMPRESSED_DIM;
      for(int col=0; col<COMPRESSED_DIM; col++) {
        double sum = 0;
        int colLo = (basis->noOfBasisFuncs*col)/COMPRESSED_DIM;
        int colHi = (basis->noOfBasisFuncs*(col+1))/COMPRESSED_DIM;
        
        for(int i=rowLo; i<rowHi; i++)
          for(int j=colLo; j<colHi; j++)
            sum += matrix[j+i*basis->noOfBasisFuncs];

        printf("%lg ", sum);
      }
      puts( row+1== COMPRESSED_DIM ? "];" : ";");
    }

  } else {

    for(int row=0; row<basis->noOfBasisFuncs; row++) {
      for(int col=0; col<basis->noOfBasisFuncs; col++)
        printf("%lg ", (double)matrix[col + row*basis->noOfBasisFuncs]);
      puts( row+1== basis->noOfBasisFuncs ? "];" : ";");
    }
  }

  ergo_free(matrix);
  delete basis;
  
  return 0;
}

