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

/** @file recexp_test.cc

    @brief  Test serial recursive expansion on a random symmetric matrix or
            matrix from a given binary file. Matrix in a binary file should contain only the upper triangle. Note: to get homo-lumo gap all matrix eigenvalues are computed.

    @author Anastasia Kruchinina <em>responsible</em>
*/

#ifndef USE_CHUNKS_AND_TASKS

#include "purification_sp2.h"
#include "purification_sp2acc.h"
#include "matrix_typedefs.h" // definitions of matrix types and interval type (source)
#include "realtype.h"   // definitions of types (utilities_basic)
#include "matrix_utilities.h"
#include "integral_matrix_wrappers.h"
#include "SizesAndBlocks.h"
#include "Matrix.h"
#include "Vector.h"
#include "MatrixSymmetric.h"
#include "MatrixTriangular.h"
#include "MatrixGeneral.h"
#include "VectorGeneral.h"
#include "output.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>

#include "random_matrices.h"
#include "get_eigenvectors.h"

typedef ergo_real real;
typedef symmMatrix MatrixType;
typedef MatrixType::VectorType VectorType;

#define SQRT_EPSILON_REAL    template_blas_sqrt(mat::getMachineEpsilon<real>())
#define MAX_DOUBLE std::numeric_limits<real>::max()
#define MIN_DOUBLE std::numeric_limits<real>::min()

real TOL_ERR_SUBS_DEFAULT = SQRT_EPSILON_REAL;
real TOL_TRACE_ERROR_DEFAULT = SQRT_EPSILON_REAL;

#ifdef PRECISION_SINGLE
real TOL_EIGENSOLVER_ACC_DEFAULT = 1e-6;
#elif PRECISION_DOUBLE
real TOL_EIGENSOLVER_ACC_DEFAULT = 1e-12;
#elif PRECISION_LONG_DOUBLE
real TOL_EIGENSOLVER_ACC_DEFAULT = 1e-16;
#elif PRECISION_QUAD_FLT128
real TOL_EIGENSOLVER_ACC_DEFAULT = 1e-24;
#endif



int main(int argc, char *argv[])
{
  printf("Program performing the recursive expansion on a given matrix.\n");
  printf("Written by Anastasia Kruchinina, Nov 2017\n");
  printf("Possible usage:\n");
  printf("    %s \n", argv[0]);
  printf("              Create a random symmetric matrix F and run recursive expansion. Matrix size is N=100, number of occupied orbitals is N/2 and allowed error in the invariant subspace is given. Homo-lumo gap is computed by explicit eigenvalue computation. \n");
  printf("    %s N_occ err_sub N \n", argv[0]);
  printf("              Create a random symmetric matrix F of size N and run recursive expansion. Number of occupied orbitals is N_occ, allowed error in the invariant subspace is set to err_sub. Homo-lumo gap is computed by explicit eigenvalue computation. \n");
  printf("    %s N_occ err_sub N filenameF \n", argv[0]);
  printf("              Read symmetric matrix F of size N from a file filenameF and run recursive expansion. Homo-lumo gap is computed by explicit eigenvalue computation. \n");
  printf("    %s N_occ err_sub N filenameF homo_lower homo_upper lumo_lower lumo_upper\n", argv[0]);
  printf("              Read symmetric matrix F of size N from a file filenameF and run recursive expansion. Bounds for homo of F are [homo_lower,homo_upper], bounds for lumo of F are [lumo_lower,lumo_upper]. \n");
  
  printf("\n");
  
  #ifdef _OPENMP
  int defThreads;
  const char *env = getenv("OMP_NUM_THREADS");
  if ( !(env && (defThreads=atoi(env)) > 0) ) {
    defThreads = 1;
  }
  
  mat::Params::setNProcs(defThreads);
  mat::Params::setMatrixParallelLevel(2);
  std::cout<<"OpenMP is used, number of threads set to "
  <<mat::Params::getNProcs()<<". Matrix parallel level: "
  <<mat::Params::getMatrixParallelLevel()<<"."<<std::endl;
  #endif
  
  int N = 100;
  int N_occ = N/2;
  real err_sub = TOL_ERR_SUBS_DEFAULT;
  char filename_F[100];
  bool use_random_matrix_F = false;
  real homo_in = MAX_DOUBLE;
  real homo_out = MIN_DOUBLE;
  real lumo_in = MIN_DOUBLE;
  real lumo_out = MAX_DOUBLE;
  bool homo_lumo_bounds_are_given = false;
  
  MatrixType F;
  
  printf("Use spectral norm for the truncation and the stopping criterion.\n");
  mat::normType normPuri = mat::euclNorm;
  mat::normType normPuriStopCrit = mat::euclNorm;
  
  printf("\n\n");
  
  if(argc == 1)
  {
    use_random_matrix_F = true;
  }
  else if(argc == 4)
  {
    use_random_matrix_F = true;
    N_occ = atoi(argv[1]);
    err_sub = atof(argv[2]);
    N = atoi(argv[3]);
  }
  else if(argc == 5)
  {
    N_occ = atoi(argv[1]);
    err_sub = atof(argv[2]);
    N = atoi(argv[3]);
    strcpy(filename_F, argv[4]);
  }
  else if(argc == 9)
  {
    N_occ = atoi(argv[1]);
    err_sub = atof(argv[2]);
    N = atoi(argv[3]);
    strcpy(filename_F, argv[4]);
    homo_out = atof(argv[5]);
    homo_in = atof(argv[6]);
    lumo_in = atof(argv[7]);
    lumo_out = atof(argv[8]);
    homo_lumo_bounds_are_given = true;
  }
  else
  {
    printf("Wrong usage!\n");
    return EXIT_FAILURE;
  }
  
  if(use_random_matrix_F)
  {
    get_random_matrix(N, F);
    printf("Created random symmetric matrix F.\n");
  }
  else
  {
    // READ MATRIX FROM FILE
    printf("Reading matrix %s\n", filename_F);
    try
    {
      std::vector<int> I;
      std::vector<int> J;
      std::vector<real> val;
      read_matrix_from_bin(filename_F, I, J, val, N, N);
      init_matrix<MatrixType>(F, N);
      F.assign_from_sparse(I, J, val);
    }
    catch(const runtime_error& error)
    {
      printf("%s\n", error.what());
      return EXIT_FAILURE;
    }
    printf("Done with reading matrix F.\n");
    assert(N == F.get_nrows());
  }
  
  
  if(!homo_lumo_bounds_are_given)
  {
    // Get all eigenvalues of F. We need this to get bounds for homo and lumo for F.
    std::vector<ergo_real> eigvalList;
    get_all_eigenvalues_of_matrix(eigvalList, F);
    
    ergo_real homo = eigvalList[N_occ-1];
    ergo_real lumo = eigvalList[N_occ  ];
    
    ergo_real epsilon_for_homo_lumo_intervals = 1e-4;
    
    homo_out = homo-epsilon_for_homo_lumo_intervals;
    homo_in  = homo+epsilon_for_homo_lumo_intervals;
    lumo_in  = lumo-epsilon_for_homo_lumo_intervals;
    lumo_out = lumo+epsilon_for_homo_lumo_intervals;
  }
  
  // SET HOMO AND LUMO BOUNDS
  IntervalType homo_bounds(homo_out, homo_in);
  IntervalType lumo_bounds(lumo_in, lumo_out);
  
  
  printf("N = %d\n", N);
  printf("N_occ = %d\n", N_occ);
  printf("err_sub = %e\n", (double)err_sub);
  printf("homo bounds: [%lf, %lf]\n", (double)homo_out, (double)homo_in);
  printf("lumo bounds: [%lf, %lf]\n", (double)lumo_in, (double)lumo_out);
  
  
  // JUST A CHECK...
  if( homo_bounds.empty() )
  {
    printf("Interval homo_bounds is empty.\n");
    return EXIT_FAILURE;
  }
  if ( lumo_bounds.empty() )
  {
    printf("Interval lumo_bounds is empty.\n");
    return EXIT_FAILURE;
  }
  
  
  
  
  int maxit = 100;
  real err_eig = 1e-3; // it is not needed for the new stopping criterion
  
  
  Purification_sp2<MatrixType> Puri;
  //Purification_sp2acc<MatrixType> Puri;
  
  string eigenvectors_method = "square";
  string eigenvectors_iterative_method = "lanczos";
  real eigensolver_accuracy = TOL_EIGENSOLVER_ACC_DEFAULT;
  int eigensolver_maxiter = 200;
  VectorType eigVecHOMO;
  VectorType eigVecLUMO;
  
  IntervalType dummy;
  
  Puri.set_eigenvectors_params(eigenvectors_method,
    eigenvectors_iterative_method,
    eigensolver_accuracy,
    eigensolver_maxiter,
    0,
    0,
    0,
    &eigVecLUMO,
    &eigVecHOMO
  );
  
  //Puri.set_compute_eigenvectors_in_each_iteration();
  
  // enable_printf_output(); // write more debug info about each iteration
  
  Puri.initialize(F,
    lumo_bounds,
    homo_bounds,
    maxit,
    err_sub,
    err_eig,
    1, // 1 = new, 0 = old stopping criterion
    normPuri,
    normPuriStopCrit,
    N_occ);
    
    
    
    // RUN RECURSIVE EXPANSION
    printf("Start recursive expansion...\n");
    Puri.PurificationStart();
    Puri.info.print_collected_info_printf();
    
    // CHECK RESULT OF THE RECURSIVE EXPANSION
    if (Puri.info.converged != 1)
    {
      throw std::runtime_error("SP2 did not converge!");
    }
    MatrixType X(Puri.X);
    ergo_real traceX = X.trace();
    if (template_blas_fabs(traceX - N_occ) > TOL_TRACE_ERROR_DEFAULT)
    {
      throw std::runtime_error("SP2: Wrong value of trace! (abs(traceX - N_occ) > TOL_TRACE_ERROR_DEFAULT)");
    }
    
    printf("\n");
    
    real homo_eigenvalue = eigvec::compute_rayleigh_quotient<real>(F, eigVecHOMO);
    real lumo_eigenvalue = eigvec::compute_rayleigh_quotient<real>(F, eigVecLUMO);
    printf("HOMO eigenvalue is %.8lf, homo bounds: [%lf, %lf]\n", (double)homo_eigenvalue, (double)homo_out, (double)homo_in);
    printf("LUMO eigenvalue is %.8lf, lumo bounds: [%lf, %lf]\n", (double)lumo_eigenvalue, (double)lumo_in, (double)lumo_out);
    
    printf("\n");
    
    // OUTPUT FIGURES
    // ostringstream name;
    // name << "out_norm" << ".m";
    // Puri.gen_matlab_file_norm_diff(name.str().c_str());
    // name.str("");
    //
    // name << "out_tau" << ".m";
    // Puri.gen_matlab_file_threshold(name.str().c_str());
    // name.str("");
    //
    // name << "out_nnz"  << ".m";
    // Puri.gen_matlab_file_nnz(name.str().c_str());
    // name.str("");
    //
    // name << "out_eigs"  << ".m";
    // Puri.gen_matlab_file_eigs(name.str().c_str());
    // name.str("");
    //
    // name << "out_cond"  << ".m";
    // Puri.gen_matlab_file_cond_num(name.str().c_str());
    // name.str("");
    //
    // name << "out_time"  << ".m";
    // Puri.gen_matlab_file_time(name.str().c_str());
    
    
    // SAVE COMPUTED DENSITY MATRIX TO FILE
    
    // vector<int>  I, J;
    // vector<real> V;
    // Puri.X.get_all_values(I, J, V);
    // ostringstream name; name << "Fout.bin";
    // write_matrix_to_bin(name.str().c_str(), I, J, V, F.get_nrows());
    
    
    F.clear();
    Puri.clear();
    
    printf("DONE!\n");
    
    
    return EXIT_SUCCESS;
}

#endif
