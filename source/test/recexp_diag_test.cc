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



/** @file recexp_diag_test.cc

    \brief Test serial recursive expansion on a diagonal matrix.

    @author Anastasia Kruchinina <em>responsible</em>
*/


#include "purification_sp2.h"
#include "purification_sp2acc.h"
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
#include <string.h>

#include "random_matrices.h"

using namespace std;

typedef ergo_real real;
typedef symmMatrix MatrixType;
typedef symmMatrixWrap MatrixTypeWrap;
typedef MatrixType::VectorType VectorType; 


#define SQRT_EPSILON_REAL    template_blas_sqrt(mat::getMachineEpsilon<real>())

real TOL_ERR_SUBS_DEFAULT = 0;  // no truncation
real TOL_TRACE_ERROR_DEFAULT = SQRT_EPSILON_REAL;


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
    double gap;
    double gap_around;
    int eig_option;

    if(use_default == 0)
      {
	if(argc < 6 || argc > 7)
	  {
	    printf("Usage: %s N N_occ rand_seed gap gap_around [option]\n", argv[0]);
	    printf("       where option is: \n");
	    printf("          1 - equidistant eigenvalues outside gap (default) in [0,1] \n");
	    printf("          2 - random spectrum in [0,1]\n");
	    return EXIT_FAILURE;
	  }


	N          = atoi(argv[1]);
	N_occ      = atoi(argv[2]);
	rand_seed  = atoi(argv[3]);
	gap     = atof(argv[4]);
	gap_around = atof(argv[5]);
	if(argc == 7)
	  eig_option = atoi(argv[6]);
	else
	  eig_option = 1;

      }
    else // use default
      {
	N          = 10;
	N_occ      = 5;
	rand_seed  = 9187;
	gap        = 0.01;
	gap_around = 0.5;
	eig_option = 1;
      }
    

    srand(rand_seed);

    // Create random symmetric matrix F with eigenvalues eigvalList  
    std::vector<ergo_real> eigvalList(N);
     
    if(eig_option == 1)
      {
	// [0, gap_around-gap/2]
	for(int i = 1; i < N_occ; ++i)
	  eigvalList[i] = (double)i/(N_occ-1)*(gap_around-gap/2);

	for(int i = 0; i < N-N_occ; ++i)
	  eigvalList[N_occ+i] = (double)i/(N-N_occ-1)*(1-(gap_around+gap/2)) + gap_around+gap/2;
      }

    if(eig_option == 2)
      {
	// [0, gap_around-gap/2]
	for(int i = 1; i < N_occ; ++i)
	  eigvalList[i] = (double)rand()/RAND_MAX*(gap_around-gap/2);
	eigvalList[N_occ-1] = gap_around-gap/2;

	eigvalList[N_occ] = gap_around+gap/2;
	for(int i = 1; i < N-N_occ-1; ++i)
	  eigvalList[N_occ+i] = (double)rand()/RAND_MAX*(1-(gap_around+gap/2)) + gap_around+gap/2;
	eigvalList[N-1] = 1;

	sort(eigvalList.begin(), eigvalList.end());  
      }


    printf("Data for the matrix F:\n");
    printf("N         = %4d\n", N        );
    printf("N_occ     = %4d\n", N_occ    );
    printf("rand_seed = %4d\n", rand_seed);
    printf("gap       = %lf\n", gap);
    printf("gap_around %lf\n", gap_around);

    real eigmin, eigmax;
    std::vector<real>::iterator result;
    result = std::min_element(eigvalList.begin(), eigvalList.end());
    eigmin = *result;
    result = std::max_element(eigvalList.begin(), eigvalList.end());
    eigmax = *result;
    printf("Spectrum is in [%lf,  %lf]\n", (double)eigmin, (double)eigmax);

    if(eig_option == 1)
      printf("Equidistant eigenvalues in [0,1] outside gap\n");
    else printf("Random spectrum in [0,1]\n");
  
    printf("Generating matrix...\n");


  MatrixType Fin;
  init_matrix<MatrixType>(Fin, N);

  std::vector<int> rows(N), cols(N);
  for(int i = 0; i < N; ++i)
    {rows[i] = i; cols[i] = i;}

  Fin.assign_from_sparse(rows, cols, eigvalList);

  printf("Matrix generated!\n");

  ergo_real homo = eigvalList[N_occ-1];
  ergo_real lumo = eigvalList[N_occ  ];
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

  MatrixTypeWrap F;
  transform_matrix_from_to(Fin, F, params);

  mat::normType normPuri = mat::euclNorm;
  mat::normType normPuriStopCrit = mat::euclNorm;
  int maxit = 100;
  ergo_real error_eig = 0; // is not needed for new stopping criterion
  ergo_real error_sub = TOL_ERR_SUBS_DEFAULT; // no truncation

  cout << "*******************" << endl;
  cout << "   Run SP2..." << endl;
  cout << "*******************" << endl;
  try
    {
      Purification_sp2<MatrixTypeWrap> Puri;
      Puri.initialize(F,
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
      //Puri.info.print_collected_info();


      MatrixTypeWrap X(Puri.X);

      // CHECK RESULT

      if(Puri.info.converged != 1)
	throw std::runtime_error("SP2 did not converge!");

      ergo_real traceX = X.trace();
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
      Purification_sp2acc<MatrixTypeWrap> PuriACC;
      PuriACC.initialize(F,
			 lumo_bounds, 
			 homo_bounds, 
			 maxit, 
			 error_sub, 
			 error_eig,
			 1, // 1 = new, 0 = old stopping criterion 
			 normPuri, 
			 normPuriStopCrit,
			 N_occ);

      PuriACC.PurificationStart(); 
      //Puri.info.print_collected_info();


      MatrixTypeWrap Xacc(PuriACC.X);


      // CHECK RESULT

      if(PuriACC.info.converged != 1)
	throw std::runtime_error("SP2ACC did not converge!");


      ergo_real traceXacc = Xacc.trace();

      if(template_blas_fabs(traceXacc - N_occ) > TOL_TRACE_ERROR_DEFAULT)
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



  cout << "Recursive expansion test on diagonal matrix finished OK!" << endl;



  return EXIT_SUCCESS;
}

