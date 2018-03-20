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
#ifndef EIGENVECTORS_HEADER
#define EIGENVECTORS_HEADER

/************************************************/
/*           EIGENVECTORS COMPUTATIONS          */
/************************************************/


/** @file get_eigenvectors.h
 *
 *  \brief Defined namespace eigvec containing functions
 *         for computing largest eigenvalues and corresponding eigenvectors
 *         using the power method or the Lanczos method. 
 *         See function computeEigenvectors.
 *
 *  @author Anastasia Kruchinina 
 */

#include "matrix_utilities.h"
#include "integral_matrix_wrappers.h"
#include "SizesAndBlocks.h"
#include "output.h"

#include <iostream>
#include <string.h>

#include "LanczosSeveralLargestEig.h"


namespace eigvec
{
/** Get Rayleigh quotient: A = (y'Ay)/(y'y), y = eigVecPtr. */
template<typename Treal, typename MatrixType, typename VectorType>
Treal compute_rayleigh_quotient(const MatrixType& A, const VectorType& eigVec)
{
   VectorType y, Ay;

   y = eigVec;
   Treal ONE = (Treal)1.0;
   y *= ONE / y.eucl();              // y = y/norm(y)
   Ay = ONE * A * y;                 // Ay = A*y
   Treal lambda = transpose(y) * Ay; // lambda = y'*Ay
   return lambda;
}


/** Use Lanzcos method for computing eigenvectors. See function
 *  computeEigenvectors for the meaning of parameters. */
template<typename Treal, typename MatrixType, typename VectorType>
void lanczos_method(const MatrixType&        A,
                    std::vector<Treal>&      eigVal,
                    std::vector<VectorType>& eigVec,
                    int                      number_of_eigenvalues,
                    const Treal              TOL,
                    std::vector<int>&        num_iter,
                    int                      maxit = 200,
                    bool                     do_deflation = false)
{
   assert(eigVal.size() == eigVec.size());
   assert(eigVal.size() == num_iter.size());
   assert(number_of_eigenvalues >= 1);

   if (eigVec[0].is_empty())
   {
      throw std::runtime_error("Error in eigvec::lanczos_method : eigVec[0].is_empty()");
   }

   const Treal ONE = 1.0;

   if (!do_deflation)
   {
      try
      {
         VectorType y;
         y  = eigVec[0];
         y *= (ONE / y.eucl());                               // normalization


         mat::arn::LanczosSeveralLargestEig<Treal, MatrixType, VectorType> lan
            (A, y, number_of_eigenvalues, maxit);
         lan.setAbsTol(TOL);
         lan.setRelTol(TOL);
         lan.run();
         Treal acc = 0;
         lan.get_ith_eigenpair(1, eigVal[0], eigVec[0], acc);

         VectorType resVec(eigVec[0]);                                 // residual
         resVec *= eigVal[0];
         resVec += -ONE * A * eigVec[0];

         /* if(number_of_eigenvalues == 2) */
         /*   { */
         /*     lan.get_ith_eigenpair(2, eigVal[1], eigVec[1], acc); */
         /*     VectorType resVec2(eigVec[1]); // residual */
         /*     resVec2 *= eigVal[1]; */
         /*     resVec2 += -ONE * A.MATRIX * eigVec[1]; */
         /*   } */
         num_iter[0] = lan.get_num_iter();
      }
      catch (std::exception& e)
      {
         num_iter[0] = maxit;                                 // lanczos did not converge in maxIter iterations
      }
   }
   else                               // do_deflation
   {
      // use the vector stored in eigVec[0]
      if (number_of_eigenvalues > 1)
      {
         VectorType y, v1;
         v1 = eigVec[0];

         // get initial guess
         if (eigVec[1].is_empty())
         {
            throw std::runtime_error("Error in eigvec::lanczos_method : eigVec[1].is_empty()");
         }
         y  = eigVec[1];
         y *= (ONE / y.eucl());                               // normalization

         try
         {
            int num_eig = 1;                                  // just one eigenpair should be computed
            // find bounds of the spectrum
            Treal eigmin, eigmax;
            A.gershgorin(eigmin, eigmax);
            Treal sigma = eigVal[0] - eigmin;                                  // out eigenvalue to the uninteresting end of the spectrum
            mat::arn::LanczosSeveralLargestEig<Treal, MatrixType, VectorType> lan
               (A, y, num_eig, maxit, 100, &v1, sigma);
            lan.setAbsTol(TOL);
            lan.setRelTol(TOL);
            lan.run();
            Treal acc = 0;
            lan.get_ith_eigenpair(1, eigVal[1], eigVec[1], acc);

            VectorType resVec(eigVec[1]);                                  // residual
            resVec *= eigVal[1];
            resVec += -ONE * A * eigVec[1];

            num_iter[1] = lan.get_num_iter();
         }
         catch (std::exception& e)
         {
            num_iter[1] = maxit;                                  // lanczos did not converge in maxIter iterations
         }
      }
      else
      {
         throw std::runtime_error("Error in eigvec::lanczos_method :  number_of_eigenvalues <= 1");
      }
   }
}


/** Use power method for computing eigenvectors. See function
 *  computeEigenvectors for the meaning of parameters. */
template<typename Treal, typename MatrixType, typename VectorType>
void power_method(const MatrixType& A,
                  Treal&            eigVal,
                  VectorType&       eigVec,
                  const Treal       TOL,
                  int&              num_iter,
                  int               maxit = 200)
{
   VectorType  y;
   VectorType  Ay;
   VectorType  residual;
   VectorType  temp;
   Treal       lambda;
   const Treal ONE  = 1.0;
   const Treal MONE = -1.0;

   y  = eigVec;           // initial guess
   y *= (ONE / y.eucl()); // normalization

   // init
   Ay       = y;
   residual = y;
   temp     = y;

   int it = 1;
   Ay = ONE * A * y;                             // Ay = A*y

   while (it == 1 || (residual.eucl() / template_blas_fabs(lambda) > TOL && it <= maxit))
   {
      y      = Ay;
      y     *= ONE / Ay.eucl();                               // y = Ay/norm(Ay)
      Ay     = ONE * A * y;                            // Ay = A*y
      lambda = transpose(y) * Ay;                             // lambda = y'*Ay

      // r = A*y - lambda*y
      residual  = Ay;
      residual += (MONE * lambda) * y;
      //printf("residual.eucl() = %e\n", residual.eucl());

      ++it;
   }

   printf("Power method required %d iterations.\n", it - 1);

   eigVal   = lambda;
   eigVec   = y;
   num_iter = it - 1;
}


/** Function for choosing method for computing eigenvectors. */
template<typename Treal, typename MatrixType, typename VectorType>
int computeEigenvectors(const MatrixType&        A,                                /**< [in] Matrix for which to compute eigenvectors. */
                        Treal                    tol,                              /**< [in] Eigensolver tolerance. */
                        std::vector<Treal>&      eigVal,                           /**< [out] Eigenvalue(s). */
                        std::vector<VectorType>& eigVec,                           /**< [in/out] Eigenvector(s). */
                        int                      number_of_eigenvalues_to_compute, /**< [in] Number of eigenvalues which Lanczos should compute. */
                        std::string              method,                           /**< [in] Chosen eigensolver (power or Lanczos). */
                        std::vector<int>&        num_iter,                         /**< [out] Actual number of iterations (now just num_iter[0] is used). */
                        int                      maxit = 200,                      /**< [in] Maximum number of iterations. */
                        bool                     do_deflation = false              /**< [in] Use deflation with eigVec[0]. */
                        )
{
   assert(number_of_eigenvalues_to_compute >= 1);
   assert(eigVal.size() >= 1);   // note: number_of_eigenvalues may not be equal to eigVal.size()
   assert(eigVec.size() == eigVal.size());
   assert(eigVec.size() == num_iter.size());

   if (method == "power")
   {
      if (eigVal.size() > 1)
      {
         throw std::runtime_error("Error in eigvec::computeEigenvectors: computation of more "
                                  "than 1 eigenpair is not implemented for the power method.");
      }
      if (do_deflation)
      {
         throw std::runtime_error("Error in eigvec::computeEigenvectors: deflation is not implemented for the power method.");
      }
      power_method(A, eigVal[0], eigVec[0], tol, num_iter[0], maxit);
   }
   else if (method == "lanczos")
   {
      lanczos_method(A, eigVal, eigVec, number_of_eigenvalues_to_compute, tol, num_iter, maxit, do_deflation);
   }
   else
   {
      throw std::runtime_error("Error in eigvec::computeEigenvectors: unknown method.");
   }
   return 0;
}
} // namespace

#endif // EIGENVECTORS_HEADER
