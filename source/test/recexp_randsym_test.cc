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



/** @file recexp_randsym_test.cc

    \brief Test serial recursive expansion on a random matrix. Generate a       
           random sparse matrix with a given spectrum.  (Note, there is no decay of elements.)

    @author Anastasia Kruchinina <em>responsible</em>
*/

#include "purification_sp2acc.h"
#include "purification_sp2.h"
#include "matrix_typedefs.h" // definitions of matrix types and interval type (source)
#include "matrix_typedefs_chtml.h" // definitions of matrix types and interval type (source)
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
#include "transform.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>

#include "random_matrices.h"

using namespace std;


typedef ergo_real real;

#define SQRT_EPSILON_REAL template_blas_sqrt(mat::getMachineEpsilon<real>())

real TOL_ERR_SUBS_DEFAULT = SQRT_EPSILON_REAL;
real TOL_TRACE_ERROR_DEFAULT = SQRT_EPSILON_REAL;


typedef symmMatrix MatrixType;
typedef symmMatrixWrap MatrixTypeWrap;
typedef MatrixType::VectorType VectorType;


int main(int argc, char *argv[])
{
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
    
  int use_default = 0; 
  if (argc == 1) use_default = 1;
      
  int N;
  int N_occ;
  int rand_seed;
  double err_sub;

  if(use_default == 0)
    {
      if(argc != 5)
	{
	  printf("Usage: %s N N_occ rand_seed err_sub\n", argv[0]);
	  return EXIT_FAILURE;
	}
      N          = atoi(argv[1]);
      N_occ      = atoi(argv[2]);
      rand_seed  = atoi(argv[3]);
      err_sub      = atof(argv[4]);  
    }
  else // use default
    {
      N          = 242;
      N_occ      = 93;
      rand_seed  = 15217;
      err_sub   = TOL_ERR_SUBS_DEFAULT; // do some truncation
    }

  srand(rand_seed);
     
  // Create random symmetric matrix F
  MatrixType F;
  get_random_matrix(N, F);
      
  // Get all eigenvalues of F. We need this to get bounds for homo and lumo for F.
  std::vector<ergo_real> eigvalList;
  get_all_eigenvalues_of_matrix(eigvalList, F);

  printf("Data for the matrix F:\n");
  printf("N         = %4d\n", N        );
  printf("N_occ     = %4d\n", N_occ    );
  printf("rand_seed = %4d\n", rand_seed);
  printf("err_sub   = %e\n", err_sub);

  ergo_real homo = eigvalList[N_occ-1];
  ergo_real lumo = eigvalList[N_occ  ];
  printf("gap       = %lf\n", (double)(lumo-homo));
      
  ergo_real epsilon_for_homo_lumo_intervals = 1e-3;
  IntervalType homo_bounds(homo-epsilon_for_homo_lumo_intervals, homo+epsilon_for_homo_lumo_intervals);
  IntervalType lumo_bounds(lumo-epsilon_for_homo_lumo_intervals, lumo+epsilon_for_homo_lumo_intervals);

      
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
 

    #ifdef USE_CHUNKS_AND_TASKS
      int leavesSizeMax = 11;
      int blocksize = 4;
      printf("leavesSizeMax       = %d\n", leavesSizeMax);
      printf("blocksize           = %d\n", blocksize);
        #ifdef USE_CHUNKS_AND_TASKS_BSM
              ParamsType params(leavesSizeMax, blocksize, N, N);
        #else
              ParamsType params(leavesSizeMax, N, N);
        #endif
    #else
            ParamsType params;
    #endif

      MatrixTypeWrap Fw;
      transform_matrix_from_to(F, Fw, params);


  mat::normType normPuri = mat::mixedNorm;
  mat::normType normPuriStopCrit = mat::mixedNorm;
 
  int maxit = 100;
  real error_eig = 0; // is not needed for new stopping criterion
  real error_sub = err_sub;


  cout << "*******************" << endl;
  cout << "   Run SP2..." << endl;
  cout << "*******************" << endl;
  try 
    {
      Purification_sp2<MatrixTypeWrap> Puri;
      Puri.initialize(Fw,
		      lumo_bounds,  
		      homo_bounds, 
		      maxit, 
		      error_sub, 
		      error_eig,
		      1, // 1 = new, 0 = old stopping criterion 
		      normPuri, 
		      normPuriStopCrit,
		      N_occ);

      Puri.PurificationStart();
      Puri.info.print_collected_info_printf();
 
      MatrixTypeWrap X(Puri.X);

      // CHECK RESULT

      if(Puri.info.converged != 1)
	throw std::runtime_error("SP2 did not converge!");

      real traceX = X.trace();
      if(template_blas_fabs(traceX - N_occ) > TOL_TRACE_ERROR_DEFAULT)
	throw std::runtime_error("SP2: Wrong value of trace! (abs(traceX - N_occ) > TOL_TRACE_ERROR_DEFAULT)");
    }
  catch(char const* e)
    {
      cerr << e << endl;
      return EXIT_FAILURE;
    }
  catch(const std::exception& e) 
    {
      cerr << e.what() << endl;
      return EXIT_FAILURE;
    } 



  cout << "*******************" << endl;
  cout << "   Run SP2ACC..." << endl;
  cout << "*******************" << endl;
  try
    {
      Purification_sp2acc<MatrixTypeWrap> Puri;
      Puri.initialize(Fw,
		      lumo_bounds,  
		      homo_bounds, 
		      maxit, 
		      error_sub, 
		      error_eig,
		      1, // 1 = new, 0 = old stopping criterion 
		      normPuri, 
		      normPuriStopCrit,
		      N_occ);

      Puri.PurificationStart();
      Puri.info.print_collected_info_printf();
 
      MatrixTypeWrap X(Puri.X);

      // CHECK RESULT
      if(Puri.info.converged != 1)
        throw std::runtime_error("SP2ACC did not converge!");
      
      real traceX = X.trace();
      if(template_blas_fabs(traceX - N_occ) > TOL_TRACE_ERROR_DEFAULT)
        throw std::runtime_error("SP2ACC: Wrong value of trace! (abs(traceX - N_occ) > TOL_TRACE_ERROR_DEFAULT)");
    }
    catch(char const* e)
    {
      cerr << e << endl;
      return EXIT_FAILURE;
    }
    catch(const std::exception& e) 
    {
      cerr << e.what() << endl;
      return EXIT_FAILURE;
    } 


  cout << "Recursive expansion test on random matrix finished OK!" << endl;


  return EXIT_SUCCESS;
}
