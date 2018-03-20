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


/** @file random_matrices.cc

    @brief File containing definitions of functions required for testing purposes.
    Functions include generation of the random dense matrices, random sparse symmetric matrices, initialization of the 
    hierarchical matrix structure, work with files and printing matrix to the screen.

    @author Anastasia Kruchinina 
    @sa random_matrices.h
*/


#include "random_matrices.h"



void print_matrix(std::vector<ergo_real> const &A)
{
  int N = sqrt(A.size());
  for(int i = 0; i < N*N; ++i )
    {
      if( i%N == 0 && i != 0 ) printf("\n");
      printf("%lf ", (double)A[i]);
    }
  printf("\n");    
}



void get_random_matrix(int N, MatrixTypeInner &X)
{
  init_matrix<MatrixTypeInner>(X, N);

  assert(X.get_nrows()*X.get_ncols() == N*N);
  // Create vectors I, J, val.
  vector<int> I(N*N);
  vector<int> J(N*N);
  vector<ergo_real> val(N*N);
  int count = 0;
  for(int i = 0; i < N; i++)
    for(int j = 0; j <= i; j++) {
      real rand0to1 = ((real)rand()) / RAND_MAX;
      real rand0to2 = rand0to1 * 2;
      real rand_m1_to_p1 = rand0to2 - 1;
      // Use denominator that becomes large when far away from diagonal, to get some sparsity from truncation
      real denominator = std::pow((i+1-j), 5);
      I[count] = i;
      J[count] = j;
      val[count] = rand_m1_to_p1 / denominator;
      count++;
    }
  I.resize(count);
  J.resize(count);
  val.resize(count);
  X.assign_from_sparse(J, I, val);
}




void get_all_eigenvalues_of_matrix(std::vector<ergo_real> & eigvalList, const MatrixTypeInner & M) {
  assert(M.get_nrows() == M.get_ncols());
  int n = M.get_nrows();
  assert(n > 0);
  eigvalList.resize(n);
  std::vector<ergo_real> M_full;
  M.fullMatrix(M_full);
  int lwork = 3*n*n;
  std::vector<ergo_real> work(lwork);
  std::vector<ergo_real> A(n*n);
  memcpy(&A[0], &M_full[0], n*n*sizeof(ergo_real));
  int info = -1;
  mat::syev("N", "U", &n, &A[0],
	    &n, &eigvalList[0], &work[0], &lwork, 
	    &info);
  assert(info == 0);  
}







void sprandsym(int N, MatrixTypeInner &X, MatrixGeneral &Q, vector<ergo_real> &D, const double MATRIX_SPARSITY)
{
  // determine numit
  int numit = N;

  real val;
  int count;
  int k, l;
  real c, s, theta, a_kk, a_kl, a_ll;

  init_matrix<MatrixTypeInner>(X, N);
  assert(X.get_nrows()*X.get_ncols() == N*N);

  // for temp values
  vector<int> Ig;
  vector<int> Jg;
  vector<real> Vg;

  // save values to add to the matrix
  vector<int> I(N);
  vector<int> J(N);
  vector<real> V;

  for(int i = 0; i < N; ++i)
    {I[i] = i; J[i] = i;}
  X.add_values(I, J, D); // set diagonal

  // matrix with eigenvectors
  init_matrix<MatrixGeneral>(Q, N);
  Q.add_identity(1);

  for(int it = 0; it < numit; ++it)
    {
      MatrixTypeInner Y;
      init_matrix<MatrixTypeInner>(Y, N);
  

      // choose angle between 0 and 2*pi
      theta = 2*PI*(real)rand()/RAND_MAX;
      s = template_blas_sin(theta);
      c = template_blas_cos(theta);
      if( N == 2)
	{ k = 0; l = 1;}
      else 
	{
	  // choose rotation plane
	  k = rand() % (N-1);   // random integer in the range 0 to N-1 
	  l = k;
	  while(l == k) l = rand() % (N-1);   // random integer in the range 0 to N-1 
	}
    
      // make l larger than k
      if(k > l) {real t = l; l = k; k = t;} 

      I.clear();
      J.clear();
      V.clear();
      I.resize(2*N);
      J.resize(2*N);
      V.resize(2*N);
      count = 0;

     Ig.resize(2);  Jg.resize(2);  Vg.resize(2);



      // A(i, k) = c*A(i, k) - s*A(i, l);
      // A(i, l) = s*A(i, k) + c*A(i, l);
    
      // i = 0:k-1, j = k
      // i = 0:k-1, j = l

     vector<real> Vk;
     Ig.resize(k);
     Jg.resize(k);
     Vk.resize(k);
     for(int i = 0; i < k; i++)
       {
	 Ig[i] = i;
	 Jg[i] = k;
       }
     X.get_values(Ig, Jg, Vk);

     vector<real> Vl;
     Vl.resize(k);
     for(int i = 0; i < k; i++)
       {
	 Jg[i] = l;
       }
     X.get_values(Ig, Jg, Vl);


      for(int i = 0; i < k; i++)
      	{
      	  // Ig[0] = i; Ig[1] = i;
      	  // Jg[0] = k; Jg[1] = l;
      	  // X.get_values(Ig, Jg, Vg); 

	  Vg[0] = Vk[i]; Vg[1] = Vl[i]; 	  

      	  I[count] = i;
      	  J[count] = k;
      	  V[count] = c*Vg[0] - s*Vg[1];
	  if(V[count] != 0)
	    count++;
      	  I[count] = i;
      	  J[count] = l;
      	  V[count] = s*Vg[0] + c*Vg[1];
      	  if(V[count] != 0)
	    count++;
      	}


      // i = k, j = l+1:N-1
      // i = l, j = l+1:N-1

      Ig.resize(N-(l+1));
      Jg.resize(N-(l+1));
      Vk.resize(N-(l+1));
      for(int i = 0; i < N-(l+1); i++)
	{
	  Ig[i] = k;
	  Jg[i] = i + (l+1);
	}
      X.get_values(Ig, Jg, Vk);
      
      Vl.resize(N-(l+1));
      for(int i = 0; i < N-(l+1); i++)
	{
	  Ig[i] = l;
	}
      X.get_values(Ig, Jg, Vl);
      

     for(int j = l+1; j < N; j++)
     	{
     	  // Ig[0] = k; Ig[1] = l;
     	  // Jg[0] = j; Jg[1] = j;
     	  // X.get_values(Ig, Jg, Vg); 

	  Vg[0] = Vk[j-(l+1)]; Vg[1] = Vl[j-(l+1)];
	  
     	  I[count] = k;
     	  J[count] = j;
     	  V[count] = c*Vg[0] - s*Vg[1];
     	  if(V[count] != 0)
	    count++;
     	  I[count] = l;
     	  J[count] = j;
     	  V[count] = s*Vg[0] + c*Vg[1];
     	  if(V[count] != 0)
	    count++;
     	}


     // i = k, j = k+1:l-1
     // i = k+1:l-1, j = l

      Ig.resize(l-(k+1));
      Jg.resize(l-(k+1));
      Vk.resize(l-(k+1));
      for(int i = 0; i < l-(k+1); i++)
	{
	  Ig[i] = k;
	  Jg[i] = i + (k+1);
	}
      X.get_values(Ig, Jg, Vk);
      
      Vl.resize(l-(k+1));
      Ig = Jg;
      for(int i = 0; i < l-(k+1); i++)
	{
	  Jg[i] = l;
	}
      X.get_values(Ig, Jg, Vl);
      

     for(int j = k+1; j < l; j++)
     	{
     	  // Ig[0] = k; Ig[1] = j;
     	  // Jg[0] = j; Jg[1] = l;
     	  // X.get_values(Ig, Jg, Vg); 

	  Vg[0] = Vk[j-(k+1)]; Vg[1] = Vl[j-(k+1)];
	  
     	  I[count] = k;
     	  J[count] = j;
     	  V[count] = c*Vg[0] - s*Vg[1];
     	  if(V[count] != 0)
	    count++;
     	  I[count] = j;
     	  J[count] = l;
     	  V[count] = s*Vg[0] + c*Vg[1];
     	  if(V[count] != 0)
	    count++;
     	}


     Ig.resize(3);  Jg.resize(3);  Vg.resize(3);
      Ig[0] = k; Ig[1] = k; Ig[2] = l;
      Jg[0] = l; Jg[1] = k; Jg[2] = l;
      X.get_values(Ig, Jg, Vg);
      a_kl = Vg[0]; a_kk = Vg[1]; a_ll = Vg[2];

      // val = (c^2 - s^2)*a_kl + s*c*(a_kk - a_ll);
      // A(k, l) = val;
      val =  (c*c - s*s)*a_kl + s*c*(a_kk - a_ll);
      I[count] = k;
      J[count] = l;
      V[count] = val;
      if(V[count] != 0)
	    count++;

      // A(k, k) = s^2*a_ll + c^2 * a_kk - 2*s*c*a_kl;
      // A(l, l) = c^2*a_ll + s^2 * a_kk + 2*s*c*a_kl;
      I[count] = k;
      J[count] = k;
      V[count] = s*s*a_ll + c*c * a_kk - 2*s*c*a_kl;
      if(V[count] != 0)
	    count++;
      I[count] = l;
      J[count] = l;
      V[count] = c*c*a_ll + s*s * a_kk + 2*s*c*a_kl;
      if(V[count] != 0)
	    count++;


      I.resize(count);
      J.resize(count);
      V.resize(count);


      vector<real> Vtmp(count);
      X.get_values(I, J, Vtmp);
      Y.assign_from_sparse(I, J, Vtmp);
      Y *= -1;
      Y.add_values(I, J, V);
      X += Y; // x_ij = x_ij + (newx_ij - x_ij)


    // construct Q'
      
      MatrixGeneral Qtmp;
      init_matrix<MatrixGeneral>(Qtmp, N);
      Ig.resize(4);  Jg.resize(4);  Vg.resize(4);
      Ig[0] = k; Ig[1] = l; Ig[2] = k; Ig[3] = l;  
      Jg[0] = k; Jg[1] = l; Jg[2] = l; Jg[3] = k;
      Vg[0] = c-1; Vg[1] = c-1; Vg[2] = -s; Vg[3] = s;

      mat::SizesAndBlocks rows;
      mat::SizesAndBlocks cols;
      Qtmp.getRows(rows);
      Qtmp.getCols(cols);
      Qtmp.assign_from_sparse(Ig, Jg, Vg, rows, cols);
      Qtmp.add_identity(1);

      MatrixGeneral QQ;
      init_matrix<MatrixGeneral>(QQ, N);
      QQ = Qtmp*Q; // QQ cannot be Qtmp or Q
      Q = QQ;

      // check sparsity
      size_t nnz = X.nnz();
      double sparsity = (double)nnz/(N*N) * 100;
      if(it % 10 == 0)
      if(sparsity > MATRIX_SPARSITY) return;
    }


}





// filename - the name of the .mtx file
int get_matrix_from_sparse(char *filename, MatrixTypeInner &X)
{
  vector<int> I, J;
  vector<ergo_real> val;
  int N, M;
  if(read_matrix_from_mtx(filename, I, J, val, N, M) == -1) return -1;
  assert(N==M);
  init_matrix<MatrixTypeInner>(X, N);
  assert(X.get_nrows()*X.get_ncols() == N*N);
  X.assign_from_sparse(I, J, val);
  return 1;
}



// filename - the name of the .mtx file
int get_matrix_from_sparse_vec(char *filename, std::vector<int> &I, std::vector<int> &J, std::vector<real> &val)
{
  I.clear();
  J.clear();
  val.clear();

  int N, M;
  if(read_matrix_from_mtx(filename, I, J, val, N, M) == -1) return -1;

  assert(N==M);

  return 1;
}


// filename - the name of the binary file
int get_matrix_from_binary(char *filename, MatrixTypeInner &X)
{
  int N;
 
  /*********** Read matrix from file */
  
  vector<ergo_real> A;
  if(read_matrix(filename, A, N, N, 1) == -1) return -1;

  init_matrix<MatrixTypeInner>(X, N);
  assert(X.get_nrows()*X.get_ncols() == N*N);

  X.assignFromFull(A);
  
  return 1;
}


// filename - the name of the binary file
int get_matrix_from_binary_vec(char *filename, std::vector<int> &I, std::vector<int> &J, std::vector<real> &val, int &N)
{
  /*********** Read matrix from file */
  vector<ergo_real> A;
  if(read_matrix(filename, A, N, N, 1) == -1) return -1;

    I.resize((N*N + N)/2);
    J.resize((N*N + N)/2);
    val.resize((N*N + N)/2);
    int count = 0;
    for (int r = 0; r < N; r++) 
      for (int c = r; c < N; c++){
	I[count] = r;
	J[count] = c;
	val[count] = A[r*N+c];
	count++;
      }  

  return 1;
}



// filename - the name of the .txt file
int get_matrix_from_full(char * filename, MatrixTypeInner &X)
{
  int N;
 
  /*********** Read matrix from file */
  
  vector<ergo_real> A;
  if(read_matrix(filename, A, N, N, 0) == -1) return -1;

  init_matrix<MatrixTypeInner>(X, N);

  assert(X.get_nrows()*X.get_ncols() == N*N);

  X.assignFromFull(A);

  return 1;
}

