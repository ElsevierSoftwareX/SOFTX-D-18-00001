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

/** @file general.h

    \brief Some general utilities used by other parts of the
    hierarchical matrix library.
*/

#ifndef MAT_GENERAL
#define MAT_GENERAL
#include <cassert>
namespace mat {

  

    template<class Treal>
    static Treal maxdiff(const Treal* f1,const Treal* f2,int size) {
    Treal diff = 0;
    Treal tmpdiff;
    for(int i = 0; i < size * size; i++) {
      tmpdiff = template_blas_fabs(f1[i] - f2[i]);
      if (tmpdiff > 0)
	diff = (diff > tmpdiff ? diff : tmpdiff);
    }
    return diff;
  }

  template<class Treal>
    static Treal maxdiff_tri(const Treal* f1,const Treal* f2,int size) {
    Treal diff = 0;
    Treal tmpdiff;
    for (int col = 0; col < size; col++)
      for (int row = 0; row < col + 1; row++) {
	tmpdiff = template_blas_fabs(f1[col * size + row] - f2[col * size + row]);
	diff = (diff > tmpdiff ? diff : tmpdiff);
      }
    return diff;
  }


  template<class Treal>
    static Treal frobdiff(const Treal* f1,const Treal* f2,int size) {
    Treal diff = 0;
    Treal tmp;
    for(int i = 0; i < size * size; i++) {
      tmp = f1[i] - f2[i];
      diff += tmp * tmp;
    }
    return template_blas_sqrt(diff);
  }

#if 0
  template<class T>
    static void fileread(T *ptr,int size,FILE*) {
    std::cout<<"error reading file"<<std::endl;
  }
  template<>
    void fileread<double>(double *ptr,int size,FILE* file) {
    fread(ptr,sizeof(double),size*size,file);
  }
  template<>
    void fileread<float>(float *ptr,int size,FILE* file) {
    double* tmpptr=new double [size*size];
    fread(tmpptr,sizeof(double),size*size,file);
    for (int i=0;i<size*size;i++)
      {
	ptr[i]=(float)tmpptr[i];
      }
    delete[] tmpptr;
  }
#else
  template<typename Treal, typename Trealonfile>
    static void fileread(Treal *ptr, int size, FILE* file) {
    if (sizeof(Trealonfile) == sizeof(Treal))
      fread(ptr,sizeof(Treal),size,file);
    else {
      Trealonfile* tmpptr=new Trealonfile[size];
      fread(tmpptr,sizeof(Trealonfile),size,file);
      for (int i = 0; i < size; i++) {
	ptr[i]=(Treal)tmpptr[i];
      }
      delete[] tmpptr;
    }
  }
#endif

  template<typename Treal, typename Tmatrix>
    static void read_matrix(Tmatrix& A, 
			    char const * const matrixPath, 
			    int const size) {
    FILE* matrixfile=fopen(matrixPath,"rb");
    if (!matrixfile) {
      throw Failure("read_matrix: Cannot open inputfile");
    }
    Treal* matrixfull = new Treal [size*size];
    fileread<Treal, double>(matrixfull, size*size, matrixfile);
    /* A must already have built data structure */
    A.assign_from_full(matrixfull, size, size);
    delete[] matrixfull;
    return;
  }

  template<typename Treal, typename Trealonfile, typename Tmatrix>
    static void read_sparse_matrix(Tmatrix& A, 
				   char const * const rowPath,
				   char const * const colPath,
				   char const * const valPath,
				   int const nval) {
    FILE* rowfile=fopen(rowPath,"rb");
    if (!rowfile) {
      throw Failure("read_matrix: Cannot open inputfile rowfile");
    }
    FILE* colfile=fopen(colPath,"rb");
    if (!colfile) {
      throw Failure("read_matrix: Cannot open inputfile colfile");
    }
    FILE* valfile=fopen(valPath,"rb");
    if (!valfile) {
      throw Failure("read_matrix: Cannot open inputfile valfile");
    }
    int* row = new int[nval];
    int* col = new int[nval];
    Treal* val = new Treal[nval];
    fileread<int, int>(row, nval, rowfile);
    fileread<int, int>(col, nval, colfile);
    fileread<Treal, Trealonfile>(val, nval, valfile);

    /* A must already have built data structure */
    A.assign_from_sparse(row, col, val, nval);
#if 0
    Treal* compval = new Treal[nval];
    A.get_values(row, col, compval, nval);
    Treal maxdiff = 0;
    Treal diff;
    for (int i = 0; i < nval; i++) {
      diff = template_blas_fabs(compval[i] - val[i]);
      maxdiff = diff > maxdiff ? diff : maxdiff;
    }
    std::cout<<"Maxdiff: "<<maxdiff<<std::endl;
#endif
    delete[] row;
    delete[] col;
    delete[] val;
    return;
  } 
  
  template<typename Treal>
    static void read_xyz(Treal* x, Treal* y, Treal* z, 
			 char * atomsPath, 
			 int const natoms,
			 int const size) {
    char* atomfile(atomsPath);
    std::ifstream input(atomfile);
    if (!input) {
      throw Failure("read_xyz: Cannot open inputfile");
    }
    input >> std::setprecision(10);
    Treal* xtmp = new Treal[natoms];  
    Treal* ytmp = new Treal[natoms];
    Treal* ztmp = new Treal[natoms];
    int* atomstart = new int[natoms+1];  
    for(int i = 0 ; i < natoms ; i++) {
      input >> x[i];
      input >> y[i];
      input >> z[i];
      input >> atomstart[i];    
    }
    atomstart[natoms] = size;  
    for (int atom = 0; atom < natoms; atom++)
      for (int bf = atomstart[atom]; bf < atomstart[atom + 1]; bf++) {
	x[bf] = x[atom];
	y[bf] = y[atom];
	z[bf] = z[atom];
      }
    delete[] xtmp;
    delete[] ytmp;
    delete[] ztmp;
    delete[] atomstart;
  }
} /* end namespace mat */

#endif
