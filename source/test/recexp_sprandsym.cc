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

 /** @file recexp_sprandsym.cc
 
    \brief   Test serial recursive expansion and computation of homo and lumo
             eigenvectors of a sparse matrix with a given eigenspectrum. Matrix is generated using Givens rotations starting from a diagonal matrix with a required eigenvalues, so eigenvectors of a result matrix are explicitly known and saved into the matrix Q. 
   
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

#define SQRT_EPSILON_REAL    template_blas_sqrt(mat::getMachineEpsilon<real>())
#define MAX_DOUBLE std::numeric_limits<real>::max()
#define MIN_DOUBLE std::numeric_limits<real>::min()

real TOL_ERR_SUBS_DEFAULT = SQRT_EPSILON_REAL;
real TOL_TRACE_ERROR_DEFAULT = SQRT_EPSILON_REAL;
real SCALAR_TOL = SQRT_EPSILON_REAL;

#ifdef PRECISION_SINGLE
real TOL_EIGENSOLVER_ACC_DEFAULT = 1e-6;
#elif PRECISION_DOUBLE
real TOL_EIGENSOLVER_ACC_DEFAULT = 1e-12;
#elif PRECISION_LONG_DOUBLE
real TOL_EIGENSOLVER_ACC_DEFAULT = 1e-16;
#elif PRECISION_QUAD_FLT128
real TOL_EIGENSOLVER_ACC_DEFAULT = 1e-24;
#endif


static bool abs_compare(real a, real b)
{return (template_blas_fabs(a) < template_blas_fabs(b));}

typedef ergo_real real;
typedef symmMatrix MatrixType;
typedef normalMatrix MatrixTypeGeneral;
typedef MatrixType::VectorType VectorType;


int main(int argc, char *argv[])
{
        printf("Program performing the recursive expansion on a given matrix.\n");
        printf("Written by Anastasia Kruchinina, Nov 2017\n");
        printf("Usage: %s N N_occ rand_seed gap gap_around sparsity err_sub [option filename]\n", argv[0]);
        printf("       where option is: \n");
        printf("          1 - equidistant eigenvalues outside gap (default) in [0,1] \n");
        printf("          2 - random spectrum in [0,1]\n");
        printf("          3 - read vector with eigs in [a,b] from file 'filename'\n");
        printf("Note: matrix sparsity is given in %%.\n");

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


        int N          = 200;
        int N_occ      = 50;
        int rand_seed  = 10000;
        double gap     = 0.01;
        double gap_around = 0.5;
        double MATRIX_SPARSITY = 20;
        double err_sub   = TOL_ERR_SUBS_DEFAULT;
        int eig_option = 1;

        real homo_in = MAX_DOUBLE;
        real homo_out = MIN_DOUBLE;
        real lumo_in = MIN_DOUBLE;
        real lumo_out = MAX_DOUBLE;

        char filename[100];


        if(argc == 9 || argc == 10)
        {
                N          = atoi(argv[1]);
                N_occ      = atoi(argv[2]);
                rand_seed  = atoi(argv[3]);
                gap     = atof(argv[4]);
                gap_around = atof(argv[5]);
                MATRIX_SPARSITY = atof(argv[6]);
                err_sub   = atof(argv[7]);
                if(argc >= 9)
                        eig_option = atoi(argv[8]);
                else
                        eig_option = 1;

                if(argc == 10)
                        strcpy(filename, argv[9]);

                if(eig_option == 3 && argc < 10)
                {
                        printf("Filename is not specified!\n");
                        return EXIT_FAILURE;
                }
        }

        printf("Data for the matrix F:\n");
        printf("N         = %4d\n", N        );
        printf("N_occ     = %4d\n", N_occ    );
        printf("rand_seed = %4d\n", rand_seed);
        printf("gap       = %lf\n", gap);
        printf("gap_around  %lf\n", gap_around);
        printf("err_sub   = %e\n", (double)err_sub);
        printf("sparsity  = %lf\n", MATRIX_SPARSITY);
        if(eig_option == 1)
                printf("Equidistant eigenvalues in [0,1] outside gap\n");
        else if(eig_option == 2)
                printf("Random spectrum in [0,1]\n");
        else if(eig_option == 3)
                printf("Read vector with eigenvalues from file %s\n", filename);
        else
                printf("Random spectrum\n");


        srand(rand_seed);

        // Create random symmetric matrix F with eigenvalues eigvalList
        std::vector<real> eigvalList(N);

        if(eig_option == 1)
        {
                // [0, gap_around-gap/2]
                for(int i = 1; i < N_occ; ++i)
                        eigvalList[i] = (real)i/(N_occ-1)*(gap_around-gap/2);

                for(int i = 0; i < N-N_occ; ++i)
                        eigvalList[N_occ+i] = (real)i/(N-N_occ-1)*(1-(gap_around+gap/2)) + gap_around+gap/2;

        }

        if(eig_option == 2)
        {
                // [0, gap_around-gap/2]
                for(int i = 1; i < N_occ; ++i)
                        eigvalList[i] = (real)rand()/RAND_MAX*(gap_around-gap/2);
                eigvalList[N_occ-1] = gap_around-gap/2;

                eigvalList[N_occ] = gap_around+gap/2;
                for(int i = 1; i < N-N_occ-1; ++i)
                        eigvalList[N_occ+i] = (real)rand()/RAND_MAX*(1-(gap_around+gap/2)) + gap_around+gap/2;
                eigvalList[N-1] = 1;

                sort(eigvalList.begin(), eigvalList.end());

                write_vector("vector_eigs.txt", eigvalList);
        }


        if(eig_option == 3)
        {
                if(read_vector(filename, eigvalList, N, false) == -1)
                        return EXIT_FAILURE;

                gap = eigvalList[N_occ] - eigvalList[N_occ-1];
                gap_around = (eigvalList[N_occ] + eigvalList[N_occ-1])/2;
        }



        // for(int i = 0; i < N; ++i)
        //   printf("%lf ", (double)eigvalList[i]);


        printf("Generating matrix...\n");

        MatrixType F;
        MatrixTypeGeneral Q;
        sprandsym(N, F, Q, eigvalList, MATRIX_SPARSITY);

        printf("Matrices generated:\n");
        printf("   F - Fock matrix, Q - matrix of eigenvectors of F\n");
        printf("NNZ in F: %lu, sparsity is %lf %%\n", F.nnz(), (double)F.nnz()/(N*N) * 100);

        // // write F to the file
        // char filenameF[] = "F.bin";
        // std::vector<real> vals;
        // F.fullMatrix(vals);
        // ofstream f (filenameF, ios::out | ios::binary);
        // f.write ((char*)&vals[0], vals.size() * sizeof(real));
        // f.close();
        // 
        // // write Q to the file
        // char filenameQ[] = "Q.bin";
        // vals.clear();
        // Q.fullMatrix(vals);
        // ofstream fQ (filenameQ, ios::out | ios::binary);
        // fQ.write ((char*)&vals[0], vals.size() * sizeof(real));
        // fQ.close();
        // 
        // printf("Matrices F and Q are saved into files.\n");

        printf("Use spectral norm for the truncation and the stopping criterion.\n");
        mat::normType normPuri = mat::euclNorm;
        mat::normType normPuriStopCrit = mat::euclNorm;

        printf("\n\n");


        ergo_real homo = eigvalList[N_occ-1];
        ergo_real lumo = eigvalList[N_occ  ];

        ergo_real epsilon_for_homo_lumo_intervals = 1e-4;
        homo_out = homo-epsilon_for_homo_lumo_intervals;
        homo_in  = homo+epsilon_for_homo_lumo_intervals;
        lumo_in  = lumo-epsilon_for_homo_lumo_intervals;
        lumo_out = lumo+epsilon_for_homo_lumo_intervals;

        // SET HOMO AND LUMO BOUNDS
        IntervalType homo_bounds(homo_out, homo_in);
        IntervalType lumo_bounds(lumo_in, lumo_out);

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

        real homo_computed = eigvec::compute_rayleigh_quotient<real>(F, eigVecHOMO);
        real lumo_computed = eigvec::compute_rayleigh_quotient<real>(F, eigVecLUMO);
        printf("Computed HOMO eigenvalue is %.8lf, exact HOMO is %.8lf, homo bounds: [%lf, %lf]\n", (double)homo_computed, (double)homo, (double)homo_out, (double)homo_in);
        printf("Computed LUMO eigenvalue is %.8lf, exact LUMO is %.8lf, lumo bounds: [%lf, %lf]\n", (double)lumo_computed, (double)lumo, (double)lumo_in, (double)lumo_out);

        assert(template_blas_fabs(lumo_computed - lumo) < SCALAR_TOL);
        assert(template_blas_fabs(homo_computed - homo) < SCALAR_TOL);

        printf("\n");


        // CHECK EIGENVECTORS

        vector<int> I(N);
        vector<int> JH(N), JL(N);
        vector<real> eigVecHOMOExact(N), eigVecLUMOExact(N);
        // get exact homo vector
        for(int i = 0; i < N; ++i)
        {
                I[i] = i;
                JH[i] = N_occ-1; JL[i] = N_occ;
        }
        Q.get_values(I, JH, eigVecHOMOExact);
        Q.get_values(I, JL, eigVecLUMOExact);

        vector<real> eigVecHOMO_vec, eigVecLUMO_vec;
        eigVecHOMO.fullvector(eigVecHOMO_vec);
        eigVecLUMO.fullvector(eigVecLUMO_vec);

        // let maximum component will be positive
        {
                vector<real>::iterator ind;
                ind = max_element(eigVecLUMOExact.begin(), eigVecLUMOExact.end(), abs_compare);
                real a = eigVecLUMOExact[distance(eigVecLUMOExact.begin(), ind)];
                if( a < 0 )
                        std::transform(eigVecLUMOExact.begin(), eigVecLUMOExact.end(), eigVecLUMOExact.begin(), std::negate<real>());

                ind = max_element(eigVecHOMOExact.begin(), eigVecHOMOExact.end(), abs_compare);
                a = eigVecHOMOExact[distance(eigVecHOMOExact.begin(), ind)];
                if( a < 0 )
                        std::transform(eigVecHOMOExact.begin(), eigVecHOMOExact.end(), eigVecHOMOExact.begin(), std::negate<real>());

                ind = max_element(eigVecLUMO_vec.begin(), eigVecLUMO_vec.end(), abs_compare);
                a = eigVecLUMO_vec[distance(eigVecLUMO_vec.begin(), ind)];
                if( a < 0 )
                        std::transform(eigVecLUMO_vec.begin(), eigVecLUMO_vec.end(), eigVecLUMO_vec.begin(), std::negate<real>());

                ind = max_element(eigVecHOMO_vec.begin(), eigVecHOMO_vec.end(), abs_compare);
                a = eigVecHOMO_vec[distance(eigVecHOMO_vec.begin(), ind)];
                if( a < 0 )
                        std::transform(eigVecHOMO_vec.begin(), eigVecHOMO_vec.end(), eigVecHOMO_vec.begin(), std::negate<real>());


        }

        printf("Checking eigenvectors...\n");

        for(int i = 0; i < N; ++i)
        {
                assert(template_blas_fabs(eigVecHOMO_vec[i] - eigVecHOMOExact[i]) < SCALAR_TOL);
                assert(template_blas_fabs(eigVecLUMO_vec[i] - eigVecLUMOExact[i]) < SCALAR_TOL);
        }

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


        F.clear();
        Puri.clear();

        printf("DONE!\n");


        return EXIT_SUCCESS;
}

#endif
