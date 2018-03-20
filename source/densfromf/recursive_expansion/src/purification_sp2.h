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

/** @file purification_sp2.h

    @brief SP2 recursive density matrix expansion (or density matrix
    purification).

    @author Anastasia Kruchinina <em>responsible</em>
*/


#ifndef HEADER_PURIFICATION_SP2
#define HEADER_PURIFICATION_SP2

#include "purification_general.h"

//#define DEBUG_OUTPUT

/** Purification_sp2acc is a class which provides an interface for SP2
 * recursive expansion.
 *
 * \tparam MatrixType Type of a matrix (ex. symmMatrix). */
template<typename MatrixType>
class Purification_sp2 : public PurificationGeneral<MatrixType>
{
public:

   typedef typename PurificationGeneral<MatrixType>::real             real;
   typedef typename PurificationGeneral<MatrixType>::IntervalType     IntervalType;
   typedef typename PurificationGeneral<MatrixType>::NormType         NormType;

   typedef typename PurificationGeneral<MatrixType>::VectorTypeInt    VectorTypeInt;
   typedef typename PurificationGeneral<MatrixType>::VectorTypeReal   VectorTypeReal;

   typedef generalVector VectorType;

   Purification_sp2() : PurificationGeneral<MatrixType>() {}

   void set_init_params()
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Chosen SP2 purification method");
      this->info.method = 1;

      this->gammaStopEstim = (3 - template_blas_sqrt((real)5)) / 2.0;
   }

   void get_poly(const int it, int& poly);
   void set_poly(const int it);

   void truncate_matrix(real& threshold, int it);

   void estimate_number_of_iterations(int& numit);
   void purify_X(const int it);
   void purify_bounds(const int it);
   void save_other_iter_info(IterationInfo& iter_info, int it);
   void apply_inverse_poly_vector(const int it, VectorTypeReal& bounds_from_it);

   void return_constant_C(const int it, real& Cval);

   //real apply_inverse_poly(const int it, real x);
   real apply_poly(const int it, real x);
   void apply_poly_to_eigs(const int it, real& homo, real& lumo);
   real compute_derivative(const int it, real x, real& DDf);
};

template<typename MatrixType>
void Purification_sp2<MatrixType>::set_poly(const int it)
{
   assert((int)this->VecPoly.size() > it);

   // if cannot compute polynomial using homo and lumo eigevalues, compute using trace
   if (this->VecPoly[it] == -1)
   {
      real Xtrace   = this->X.trace();
      real Xsqtrace = this->Xsq.trace();

      if ((template_blas_fabs(Xsqtrace - this->nocc) <
           template_blas_fabs(2 * Xtrace - Xsqtrace - this->nocc))
          ||
          (it % 2
           &&
           (template_blas_fabs(Xsqtrace - this->nocc) ==
            template_blas_fabs(2 * Xtrace - Xsqtrace - this->nocc))
          ))
      {
         this->VecPoly[it] = 1;
      }
      else
      {
         this->VecPoly[it] = 0;
      }
   }
}


template<typename MatrixType>
void Purification_sp2<MatrixType>::get_poly(const int it, int& poly)
{
   assert((int)this->VecPoly.size() > it);
   assert(this->VecPoly[it] != -1);
   poly = this->VecPoly[it];
}


template<typename MatrixType>
void Purification_sp2<MatrixType>::purify_X(const int it)
{
#ifdef DEBUG_OUTPUT
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Purify X...");
#endif
   int poly;

   set_poly(it);
   get_poly(it, poly);

   /* It may happen that X2 has many more nonzeros than X, for
    * example 5 times as many.  Therefore it makes sense to try
    * having only one "big" matrix in memory at a time. However,
    * file operations have proved to be quite expensive and should
    * be avoided if possible. Hence we want to achieve having only
    * one big matrix in memory without unnecessary file
    * operations. We are currently hoping that it will be ok to add
    * a "small" matrix to a "big" one, that the memory usage after
    * that operation will be like the memory usage for one big
    * matrix + one small matrix. Therefore we are adding X to X2 (X
    * is truncated, a "small" matrix) instead of the opposite when
    * the 2*X-X*X polynomial is evaluated.
    */

   if (poly == 0)
   {
      this->Xsq *= ((real) - 1.0);
      this->X *= (real)2.0;
      this->Xsq += this->X;    // Xsq = -Xsq + 2X
   }

   this->Xsq.transfer(this->X); // clear Xsq and old X
}


template<typename MatrixType>
void Purification_sp2<MatrixType>::purify_bounds(const int it)
{
#ifdef DEBUG_OUTPUT
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Change homo and lumo bounds according to the chosen polynomial VecPoly = %d", this->VecPoly[it]);
#endif
   real homo_low, homo_upp, lumo_upp, lumo_low;
   int  poly;

   get_poly(it, poly);

   if (poly == 1)
   {
      // update bounds
      homo_low          = 2 * this->homo_bounds.low() - this->homo_bounds.low() * this->homo_bounds.low(); // 2x - x^2
      homo_upp          = 2 * this->homo_bounds.upp() - this->homo_bounds.upp() * this->homo_bounds.upp(); // 2x - x^2
      lumo_low          = this->lumo_bounds.low() * this->lumo_bounds.low();                               // x^2
      lumo_upp          = this->lumo_bounds.upp() * this->lumo_bounds.upp();                               // x^2
      this->homo_bounds = IntervalType(homo_low, homo_upp);
      this->lumo_bounds = IntervalType(lumo_low, lumo_upp);
   }
   else
   {
      // update bounds
      lumo_low          = 2 * this->lumo_bounds.low() - this->lumo_bounds.low() * this->lumo_bounds.low(); // 2x - x^2
      lumo_upp          = 2 * this->lumo_bounds.upp() - this->lumo_bounds.upp() * this->lumo_bounds.upp(); // 2x - x^2
      homo_low          = this->homo_bounds.low() * this->homo_bounds.low();                               // x^2
      homo_upp          = this->homo_bounds.upp() * this->homo_bounds.upp();                               // x^2
      this->homo_bounds = IntervalType(homo_low, homo_upp);
      this->lumo_bounds = IntervalType(lumo_low, lumo_upp);
   }

   IntervalType zero_one(0, 1);
   this->homo_bounds.intersect(zero_one);
   this->lumo_bounds.intersect(zero_one);

#ifdef DEBUG_OUTPUT
   if (this->homo_bounds.empty())
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Interval homo_bounds is empty.");
   }
   if (this->lumo_bounds.empty())
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Interval lumo_bounds is empty.");
   }

   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "1-homo: [ %lf , %lf ],", this->homo_bounds.low(), this->homo_bounds.upp());
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "lumo:   [ %lf , %lf ].", this->lumo_bounds.low(), this->lumo_bounds.upp());
#endif
}


/**************************************************************************************/

template<typename MatrixType>
void Purification_sp2<MatrixType>::return_constant_C(const int it, real& Cval)
{
   Cval = C_SP2;
}


/**************************************************************************************/

template<typename MatrixType>
void Purification_sp2<MatrixType>::truncate_matrix(real& threshold, int it)
{
   real allowed_error = this->error_per_it;

   threshold = 0;
   if (allowed_error == 0)
   {
      return;
   }

   assert((int)this->VecGap.size() > it);
#ifdef DEBUG_OUTPUT
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Truncate matrix X: ");
#endif

   real tau;                 // threshold for the error matrix

   if (this->VecGap[it] > 0) // TODO:  zero or some small value?
   {
#ifdef DEBUG_OUTPUT
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "VecGap[ %d ] = %e , ", it, this->VecGap[it]);
#endif
      tau = (allowed_error * this->VecGap[it]) / (1 + allowed_error); // control error in the occupied subspace
   }
   else // if gap is not known
   {
      tau = allowed_error * 0.01;
   }

   // return the actual introduced error
#ifdef USE_CHUNKS_AND_TASKS
   threshold = (this->X).thresh_frob(tau);
#else
   threshold = (this->X).thresh(tau, this->normPuriTrunc);
#endif

#ifdef DEBUG_OUTPUT
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "tau = %e", tau);
#endif
}


/**************************************************************************************/


template<typename MatrixType>
void Purification_sp2<MatrixType>::estimate_number_of_iterations(int& estim_num_iter)
{
   int  it = 1;
   int  maxit_tmp = this->maxit + this->additional_iterations;
   real x, y;
   real epsilon = this->get_epsilon();

   // we need values on iteration it-2 for the stopping criterion
   this->check_stopping_criterion_iter = 2;

   int max_size = maxit_tmp + 1 + this->additional_iterations + 2; // largest possible size

   this->VecPoly.clear();
   this->VecPoly.resize(max_size, -1);

   this->VecGap.clear();
   this->VecGap.resize(max_size, -1);

   // we are interested in the inner bounds of gap
   x = this->lumo_bounds.upp(); // = lumo
   y = this->homo_bounds.upp(); // = 1 - homo


   // if VecGap is zero
   if (1 - x - y <= 0)
   {
#ifdef DEBUG_OUTPUT
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "VecGap cannot be computed. Set estimated number of iteration to the maxit.");
#endif
      estim_num_iter = this->maxit;
      return;
   }

   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "INIT LUMO: [ %.12lf , %.12lf ]", (double)this->lumo_bounds.low(), (double)this->lumo_bounds.upp());
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "INIT HOMO: [ %.12lf , %.12lf ]", (double)this->homo_bounds.low(), (double)this->homo_bounds.upp());



   this->VecPoly[0] = -1;
   this->VecGap[0]  = 1 - x - y;

   estim_num_iter = -1;

   while (it <= maxit_tmp)
   {
      // note: avoid not-stopping in case of idempotent matrix
      if ((x > y) || (it % 2 && (x == y))) // lumo > 1-homo
      {
         x *= x;
         y  = 2 * y - y * y;
         this->VecPoly[it] = 1;
      }
      else
      {
         x  = 2 * x - x * x;
         y *= y;
         this->VecPoly[it] = 0;
      }

      this->VecGap[it] = 1 - x - y;

      // maybe we wish to perform some more iterations, then stopping criterion suggest
      if ((estim_num_iter == -1) &&
          (x - x * x < epsilon) && (y - y * y < epsilon) &&
          (this->VecPoly[it] != this->VecPoly[it - 1])) // to apply stopping criterion, polynomials must alternate
      {
         estim_num_iter = it;
         maxit_tmp      = it + this->additional_iterations;

         if (this->normPuriStopCrit == mat::frobNorm)
         {
            estim_num_iter += 2;
            maxit_tmp      += 2;
         }
      }

      ++it;
   }  //while



   if ((estim_num_iter == -1) && (it == maxit_tmp + 1))
   {
#ifdef DEBUG_OUTPUT
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "maxit = %d number of iterations is reached in estimate_number_of_iteration()", this->maxit);
#endif
      estim_num_iter = this->maxit;
   }


   this->VecPoly.resize(maxit_tmp + 1); // make it less if needed
   this->VecGap.resize(maxit_tmp + 1);

#ifdef DEBUG_OUTPUT
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Sequence of polynomials VecPoly: ");
   for (int i = 0; i < (int)this->VecPoly.size(); ++i)
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "%d ", this->VecPoly[i]);
   }
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "");
#endif
}


/************ SAVE INFORMATION ABOUT ITERATIONS SPECIFIC FOR SP2 PURIFICATION  *************/

template<typename MatrixType>
void Purification_sp2<MatrixType>::save_other_iter_info(IterationInfo& iter_info, int it)
{
   assert((int)this->VecPoly.size() > it);
   assert((int)this->VecGap.size() > it);

   iter_info.poly = this->VecPoly[it];
   iter_info.gap  = this->VecGap[it];
}


/************ APPLY INVERSE POLYNOMIAL (FOR ESTIMATION OF EIGENVALUES) ************/


template<typename MatrixType>
void Purification_sp2<MatrixType>::apply_inverse_poly_vector(const int it, VectorTypeReal& bounds_from_it)
{
   int  poly;
   for (int i = it; i >= 1; i--)
   {
      get_poly(i, poly);

      if (poly == 1)
      {
         bounds_from_it[0] = template_blas_sqrt(bounds_from_it[0]);
         bounds_from_it[1] = template_blas_sqrt(bounds_from_it[1]);

         bounds_from_it[2] = bounds_from_it[2] / (1 + template_blas_sqrt(1 - bounds_from_it[2]));
         bounds_from_it[3] = bounds_from_it[3] / (1 + template_blas_sqrt(1 - bounds_from_it[3]));
      }
      else
      {
         bounds_from_it[0] = bounds_from_it[0] / (1 + template_blas_sqrt(1 - bounds_from_it[0]));
         bounds_from_it[1] = bounds_from_it[1] / (1 + template_blas_sqrt(1 - bounds_from_it[1]));

         bounds_from_it[2] = template_blas_sqrt(bounds_from_it[2]);
         bounds_from_it[3] = template_blas_sqrt(bounds_from_it[3]);
      }
   }
}


/*
 *
 *
 * template<typename MatrixType>
 * typename Purification_sp2<MatrixType>::real
 * Purification_sp2<MatrixType>::apply_inverse_poly(const int it, real x)
 * {
 * if( it == 0 ) return -1; // no polynomials was applied in 0 iteration
 *
 * int poly;
 *
 * real finvx = x;
 * for(int i = it; i >= 1; i--)
 * {
 * get_poly(i, poly);
 *
 * if(poly == 1)
 * finvx = template_blas_sqrt(finvx);
 * else
 * finvx = finvx/(1+template_blas_sqrt(1-finvx));
 * }
 * return finvx;
 * }
 */

template<typename MatrixType>
typename Purification_sp2<MatrixType>::real
Purification_sp2<MatrixType>::apply_poly(const int it, real x)
{
   assert(it >= 0);
   if (it == 0)
   {
      return x;
   }

   real fx;
   int  poly;
   get_poly(it, poly);

   if (poly == 1)
   {
      fx = x * x;
   }
   else
   {
      fx = 2 * x - x * x;
   }

   return fx;
}


template<typename MatrixType>
void Purification_sp2<MatrixType>::apply_poly_to_eigs(const int it, real& homo, real& lumo)
{
   assert(it >= 0);
   if (it == 0)
   {
      return;
   }

   int poly;
   get_poly(it, poly);

   if (poly == 1)
   {
      homo = 2 * homo - homo * homo;
      lumo = lumo * lumo;
   }
   else
   {
      homo = homo * homo;
      lumo = 2 * lumo - lumo * lumo;
   }
}


template<typename MatrixType>
typename Purification_sp2<MatrixType>::real
Purification_sp2<MatrixType>::compute_derivative(const int it, real x, real& DDf)
{
   assert(it > 0);

   real Df;
   real a;
   int  poly;

   a   = x;
   Df  = 1;
   DDf = 1;

   for (int i = 1; i <= it; i++)
   {
      get_poly(i, poly);

      if (poly == 1)
      {
         DDf = 2 * Df * Df + (2 * a) * DDf;
         Df *= 2 * a;
         a   = a * a;
      }
      else
      {
         DDf = -2 * Df * Df + (2 - 2 * a) * DDf;
         Df *= 2 - 2 * a;
         a   = 2 * a - a * a;
      }
   }

   //do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"Derivative (it = %d) = %lf\n", it, Df);

   return Df;
}


#endif //HEADER_PURIFICATION_SP2
