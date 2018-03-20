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

/** @file purification_sp2acc.h

    @brief SP2ACC (SP2 accelerated) recursive density matrix expansion (or density matrix
    purification).

    @author Anastasia Kruchinina <em>responsible</em>
*/


#ifndef HEADER_PURIFICATION_SP2ACC
#define HEADER_PURIFICATION_SP2ACC

#include "purification_general.h"

//#define DEBUG_OUTPUT

/** Purification_sp2acc is a class which provides an interface for
 * SP2ACC recursive expansion.
 *
 * \tparam MatrixType Type of a matrix (ex. symmMatrix). */
template<typename MatrixType>
class Purification_sp2acc : public PurificationGeneral<MatrixType>
{
public:

   typedef typename PurificationGeneral<MatrixType>::real             real;
   typedef typename PurificationGeneral<MatrixType>::IntervalType     IntervalType;
   typedef typename PurificationGeneral<MatrixType>::NormType         NormType;

   typedef typename PurificationGeneral<MatrixType>::VectorTypeInt    VectorTypeInt;
   typedef typename PurificationGeneral<MatrixType>::VectorTypeReal   VectorTypeReal;

   typedef generalVector   VectorType;

   Purification_sp2acc() : PurificationGeneral<MatrixType>() {}

   virtual void set_init_params()
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Chosen SP2 ACCELERATED purification method");
      this->info.method = 2;

      this->gammaStopEstim = 6 - 4 * template_blas_sqrt((real)2);

      this->check_stopping_criterion_iter = -1; // will be changed during purification
   }

   virtual void get_poly(const int it, int& poly, real& alpha);
   virtual void set_poly(const int it);

   virtual void truncate_matrix(real& threshold, int it);

   virtual void estimate_number_of_iterations(int& numit);
   virtual void purify_X(const int it);
   virtual void purify_bounds(const int it);
   virtual void save_other_iter_info(IterationInfo& iter_info, int it);
   virtual void apply_inverse_poly_vector(const int it, VectorTypeReal& bounds_from_it);

   virtual void return_constant_C(const int it, real& Cval);

   virtual real apply_poly(const int it, real x);
   virtual void apply_poly_to_eigs(const int it, real& homo, real& lumo);
   virtual real compute_derivative(const int it, real x, real& DDf);


   /* PARAMETERS */

   VectorTypeReal VecAlpha;

   // defined the iteration when we turn off acceleration
   static const real deltaTurnOffAcc;
};

template<typename MatrixType>
const typename Purification_sp2acc<MatrixType>::real
Purification_sp2acc<MatrixType>::deltaTurnOffAcc = 0.01;



template<typename MatrixType>
void Purification_sp2acc<MatrixType>::set_poly(const int it)
{
   assert((int)this->VecPoly.size() > it);

   // if cannot compute polynomial using homo and lumo eigevalues, compute using trace
   if (this->VecPoly[it] == -1)
   {
      real Xtrace   = this->X.trace();
      real Xsqtrace = this->Xsq.trace();

      real delta = deltaTurnOffAcc;

      // Should we turn off acceleration or not
      if ((this->check_stopping_criterion_iter == -1) && (this->lumo_bounds.low() < delta) && (this->homo_bounds.low() < delta))
      {
#ifdef DEBUG_OUTPUT
         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Outer bounds of homo and lumo are less then %e: ", delta);
         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "lumo_out = %e, homo_out = %e ", this->lumo_bounds.low(), this->homo_bounds.low());
#endif

         this->lumo_bounds = IntervalType(0, this->lumo_bounds.upp());
         this->homo_bounds = IntervalType(0, this->homo_bounds.upp());

         // start to check stopping criterion
         if (it == 1)
         {
            this->check_stopping_criterion_iter = it + 1; // in the it=0 we had the same eigenvalue bounds
         }
         else
         {
            this->check_stopping_criterion_iter = it + 2;
         }
         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Start to check stopping criterion on iteration %d", this->check_stopping_criterion_iter);
      }

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
         VecAlpha[it]      = 2 / (2 - this->lumo_bounds.low());
      }
      else
      {
         this->VecPoly[it] = 0;
         VecAlpha[it]      = 2 / (2 - this->homo_bounds.low());
      }
#ifdef DEBUG_OUTPUT
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Acceleration parameter: alpha = %lf", VecAlpha[it]);
#endif
   }
}


template<typename MatrixType>
void Purification_sp2acc<MatrixType>::get_poly(const int it, int& poly, real& alpha)
{
   assert((int)this->VecPoly.size() > it);
   assert(this->VecPoly[it] != -1);

   //check also if alpha is computed
   assert(this->VecAlpha[it] != -1);

   poly  = this->VecPoly[it];
   alpha = VecAlpha[it];
}


template<typename MatrixType>
void Purification_sp2acc<MatrixType>::purify_X(const int it)
{
#ifdef DEBUG_OUTPUT
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Purify X...");
#endif
   real alpha_tmp;
   int  poly;

   set_poly(it);

   get_poly(it, poly, alpha_tmp);

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
    * is truncated, a "small" matrix) instead of the opposite.
    */

   if (poly == 1)
   {
      if (alpha_tmp != 1)
      {
         // (1-a+a*x)^2 = (1-a)^2 + 2*(1-a)*a*x + a^2*x^2
        //  this->X.mult_scalar((real)2.0 * (1 - alpha_tmp) * alpha_tmp);
        //  this->X.add_identity((real)(1 - alpha_tmp) * (1 - alpha_tmp));
        //  this->Xsq.mult_scalar((real)alpha_tmp * alpha_tmp);
        //  this->Xsq.add(this->X);  // Xsq = (1-a+a*X)^2

        this->X *= ((real)2.0 * (1 - alpha_tmp) * alpha_tmp);
        this->X.add_identity((real)(1 - alpha_tmp) * (1 - alpha_tmp));
        this->Xsq *= ((real)alpha_tmp * alpha_tmp);
        this->Xsq += this->X;  // Xsq = (1-a+a*X)^2


      }
      else
      {
         // DO NOTHING
      }
   }
   else
   {
      if (alpha_tmp != 1)
      {
         this->X *= ((real) - 2.0 * alpha_tmp);
         this->Xsq *= ((real) - alpha_tmp * alpha_tmp);
         this->Xsq -= this->X;  // Xsq = 2*a*X - (a*X)^2
      }
      else
      {
        this->Xsq *= ((real) - 1.0);
        this->X *= (real)2.0;
        this->Xsq += this->X;    // Xsq = -Xsq + 2X


      }
   }  // if poly == 1

   this->Xsq.transfer(this->X); // clear Xsq and old X
}


template<typename MatrixType>
void Purification_sp2acc<MatrixType>::purify_bounds(const int it)
{
#ifdef DEBUG_OUTPUT
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Change homo and lumo bounds according to the chosen polynomial VecPoly = %d", this->VecPoly[it]);
#endif

   real homo_low, homo_upp, lumo_upp, lumo_low;
   real alpha_tmp;
   int  poly;

   get_poly(it, poly, alpha_tmp);

   if (poly == 1)
   {
      // update bounds
      homo_low  = 2 * alpha_tmp * this->homo_bounds.low() - alpha_tmp * alpha_tmp * this->homo_bounds.low() * this->homo_bounds.low(); // 2*a*x - (a*x)^2
      homo_upp  = 2 * alpha_tmp * this->homo_bounds.upp() - alpha_tmp * alpha_tmp * this->homo_bounds.upp() * this->homo_bounds.upp(); // 2*a*x - (a*x)^2
      lumo_low  = (1 - alpha_tmp + alpha_tmp * this->lumo_bounds.low());                                                               // (1-a+a*x)^2
      lumo_low *= lumo_low;
      lumo_upp  = (1 - alpha_tmp + alpha_tmp * this->lumo_bounds.upp());                                                               // (1-a+a*x)^2
      lumo_upp *= lumo_upp;

      this->homo_bounds = IntervalType(homo_low, homo_upp);
      this->lumo_bounds = IntervalType(lumo_low, lumo_upp);
   }
   else
   {
      // update bounds
      lumo_low  = 2 * alpha_tmp * this->lumo_bounds.low() - alpha_tmp * alpha_tmp * this->lumo_bounds.low() * this->lumo_bounds.low(); // 2*a*x - (a*x)^2
      lumo_upp  = 2 * alpha_tmp * this->lumo_bounds.upp() - alpha_tmp * alpha_tmp * this->lumo_bounds.upp() * this->lumo_bounds.upp(); // 2*a*x - (a*x)^2
      homo_low  = (1 - alpha_tmp + alpha_tmp * this->homo_bounds.low());                                                               // (1-a+a*x)^2
      homo_low *= homo_low;
      homo_upp  = (1 - alpha_tmp + alpha_tmp * this->homo_bounds.upp());                                                               // (1-a+a*x)^2
      homo_upp *= homo_upp;

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


   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "1-homo: [ %g , %g ],", this->homo_bounds.low(), this->homo_bounds.upp());
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "lumo:   [ %g , %g ].", this->lumo_bounds.low(), this->lumo_bounds.upp());
#endif
}


/*****************************************************/

template<typename MatrixType>
void Purification_sp2acc<MatrixType>::return_constant_C(const int it, real& Cval)
{
   assert(it >= 1);

   real alpha1 = VecAlpha[it - 1];
   real alpha2 = VecAlpha[it];

   Cval = -1;

#ifdef DEBUG_OUTPUT
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "alpha1 = %.4e , alpha2 = %.4e", alpha1, alpha2);
#endif

   if (it < 2)
   {
      return; // -1
   }
   // no acceleration
   if (((alpha1 == 1) && (alpha2 == 1)))
   {
#ifdef DEBUG_OUTPUT
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Check SP2 stopping criterion.");
#endif
      Cval = C_SP2;
      return;
   }
#ifdef DEBUG_OUTPUT
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Do not check stopping criterion.");
#endif
   // exit and return C = -1


   // If we want to compute the constant C on every iterations, we can use the following:
#if 0
   // get bounds from the iteration it-2
   real homo_low = this->info.Iterations[it - 2].homo_bound_low;
   real homo_upp = this->info.Iterations[it - 2].homo_bound_upp;
   real lumo_low = this->info.Iterations[it - 2].lumo_bound_low;
   real lumo_upp = this->info.Iterations[it - 2].lumo_bound_upp;
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "lumo [%.16e, %e], homo [%.16e, %e]", lumo_low, lumo_upp, homo_low, homo_upp);


   if ((homo_upp > 0.5) || (lumo_upp > 0.5))
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Inner bounds interval do not contain 0.5. Skip iteration. ");
      Cval = -1;
      return;
   }

   a = std::max(lumo_low, homo_low);

   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Chosen a = %g", a);

   if (a <= 0)
   {
      // just in case
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Cannot compute constant C since a = %g when alpha1 = %g"
                                                  " and  alpha2 = %g", a, alpha1, alpha2);
      Cval = -1;
      return;
   }

   real C1;
   C1 = -7.88 + 11.6 * alpha1 + 0.71 * alpha2;
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Local maximum of g1:  %g", C1);

   Cval = C1;

   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "*********** C = %g  ************", Cval);
#endif
}


/**************************************************************************************/

template<typename MatrixType>
void Purification_sp2acc<MatrixType>::truncate_matrix(real& threshold, int it)
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
   real tau;

   if (this->VecGap[it] > 0) // TODO:  zero or some small value?
   {
#ifdef DEBUG_OUTPUT
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "VecGap[ %d ] = %e , ", it, this->VecGap[it]);
#endif
      tau = (allowed_error * this->VecGap[it]) / (1 + allowed_error);  // control error in the occupied subspace
   }
   else // if gap is not known
   {
      tau = allowed_error * 0.01;
   }


// TODO
   #ifdef USE_CHUNKS_AND_TASKS
      threshold = (this->X).thresh_frob(tau);
   #else
      threshold = (this->X).thresh(tau, this->normPuriTrunc);
   #endif

#ifdef DEBUG_OUTPUT
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "tau = %e", tau);
#endif
}


/****************************************************************************************/



template<typename MatrixType>
void Purification_sp2acc<MatrixType>::estimate_number_of_iterations(int& estim_num_iter)
{
   int  it = 1;
   int  maxit_tmp = this->maxit + this->additional_iterations;
   real x, y, x_out, y_out;
   real alpha_tmp;
   real epsilon = this->get_epsilon();

   int max_size = maxit_tmp + 1 + this->additional_iterations + 2; // largest possible size

   this->VecPoly.clear();
   this->VecPoly.resize(max_size, -1);

   this->VecGap.clear();
   this->VecGap.resize(max_size, -1);

   VecAlpha.clear();
   VecAlpha.resize(max_size, -1);

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


   x_out = this->lumo_bounds.low();
   y_out = this->homo_bounds.low();



   real delta = deltaTurnOffAcc;

   this->VecPoly[0] = -1;
   this->VecGap[0]  = 1 - x - y;

   estim_num_iter = -1;

   while (it <= maxit_tmp)
   {
      // note: avoid not-stopping in case of idempotent matrix
      if ((x > y) || (it % 2 && (x == y))) // lumo > 1-homo
      {
         alpha_tmp = 2 / (2 - x_out);

         x  = (1 - alpha_tmp + alpha_tmp * x);
         x *= x;
         y  = 2 * alpha_tmp * y - alpha_tmp * alpha_tmp * y * y;

         x_out  = (1 - alpha_tmp + alpha_tmp * x_out);
         x_out *= x_out;
         y_out  = 2 * alpha_tmp * y_out - alpha_tmp * alpha_tmp * y_out * y_out;

         this->VecPoly[it] = 1;
      }
      else
      {
         alpha_tmp = 2 / (2 - y_out);

         x  = 2 * alpha_tmp * x - alpha_tmp * alpha_tmp * x * x;
         y  = (1 - alpha_tmp + alpha_tmp * y);
         y *= y;

         x_out  = 2 * alpha_tmp * x_out - alpha_tmp * alpha_tmp * x_out * x_out;
         y_out  = (1 - alpha_tmp + alpha_tmp * y_out);
         y_out *= y_out;

         this->VecPoly[it] = 0;
      }

      VecAlpha[it]     = alpha_tmp;
      this->VecGap[it] = 1 - x - y;

      // find iteration where x_out < delta && y_out < delta
      if ((x_out < delta) && (y_out < delta) && (this->check_stopping_criterion_iter == -1))
      {
         // start to check stopping criterion
         if (it == 1)
         {
            this->check_stopping_criterion_iter = it + 1; // in the it=0 we had the same eigenvalue bounds
         }
         else
         {
            this->check_stopping_criterion_iter = it + 2;
         }
         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Start to check stopping criterion on iteration %d", this->check_stopping_criterion_iter);
         x_out = 0;
         y_out = 0;
      }


      // maybe we wish to perform some more iterations, then stopping criterion suggest
      if ((estim_num_iter == -1) && (it >= this->check_stopping_criterion_iter) &&
          (x - x * x < epsilon) && (y - y * y < epsilon) && // the eucledian norm is less then epsilon
          (this->VecPoly[it] != this->VecPoly[it - 1]))     // to apply stopping criterion, polynomials must alternate
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
   } //while

   if ((estim_num_iter == -1) && (it == maxit_tmp + 1))
   {
#ifdef DEBUG_OUTPUT
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "maxit = %d number of iterations is reached in estimate_number_of_iteration()", this->maxit);
#endif
      estim_num_iter = this->maxit;
   }

   this->VecPoly.resize(maxit_tmp + 1); // make it less if needed
   this->VecGap.resize(maxit_tmp + 1);
   this->VecAlpha.resize(maxit_tmp + 1);
}


/************ SAVE INFORMATION ABOUT ITERATIONS SPECIFIC FOR ACC SP2 PURIFICATION **********x***/

template<typename MatrixType>
void Purification_sp2acc<MatrixType>::save_other_iter_info(IterationInfo& iter_info, int it)
{
   assert((int)this->VecPoly.size() > it);
   assert((int)this->VecGap.size() > it);
   assert((int)this->VecAlpha.size() > it);

   iter_info.poly  = this->VecPoly[it];
   iter_info.gap   = this->VecGap[it];
   iter_info.alpha = VecAlpha[it];
}


/************ APPLY INVERSE POLYNOMIAL (FOR ESTIMATION OF EIGENVALUES) ************/

template<typename MatrixType>
void Purification_sp2acc<MatrixType>::apply_inverse_poly_vector(const int it, VectorTypeReal& bounds_from_it)
{
   real tau;
   real alpha_tmp;
   int  poly;

   for (int i = it; i >= 1; i--)
   {
      tau = 0;//this->info.Iterations[i].threshold_X;

      get_poly(i, poly, alpha_tmp);

      if (poly == 1)
      {
         bounds_from_it[0] = template_blas_sqrt(bounds_from_it[0]);
         bounds_from_it[0] = (bounds_from_it[0] - 1 + alpha_tmp) / alpha_tmp - tau;
         bounds_from_it[1] = template_blas_sqrt(bounds_from_it[1]);
         bounds_from_it[1] = (bounds_from_it[1] - 1 + alpha_tmp) / alpha_tmp + tau;

         bounds_from_it[2] = bounds_from_it[2] / (1 + template_blas_sqrt(1 - bounds_from_it[2]));
         bounds_from_it[2] = bounds_from_it[2] / alpha_tmp - tau;
         bounds_from_it[3] = bounds_from_it[3] / (1 + template_blas_sqrt(1 - bounds_from_it[3]));
         bounds_from_it[3] = bounds_from_it[3] / alpha_tmp + tau;
      }
      else
      {
         bounds_from_it[0] = bounds_from_it[0] / (1 + template_blas_sqrt(1 - bounds_from_it[0]));
         bounds_from_it[0] = bounds_from_it[0] / alpha_tmp - tau;
         bounds_from_it[1] = bounds_from_it[1] / (1 + template_blas_sqrt(1 - bounds_from_it[1]));
         bounds_from_it[1] = bounds_from_it[1] / alpha_tmp + tau;

         bounds_from_it[2] = template_blas_sqrt(bounds_from_it[2]);
         bounds_from_it[2] = (bounds_from_it[2] - 1 + alpha_tmp) / alpha_tmp - tau;
         bounds_from_it[3] = template_blas_sqrt(bounds_from_it[3]);
         bounds_from_it[3] = (bounds_from_it[3] - 1 + alpha_tmp) / alpha_tmp + tau;
      }
   }
}


template<typename MatrixType>
typename Purification_sp2acc<MatrixType>::real
Purification_sp2acc<MatrixType>::apply_poly(const int it, real x)
{
   assert(it >= 0);
   if (it == 0)
   {
      return x;
   }

   real fx;
   int  poly;
   real alpha_tmp;
   get_poly(it, poly, alpha_tmp);

   if (poly == 1)
   {
      fx = (1 - alpha_tmp + alpha_tmp * x) * (1 - alpha_tmp + alpha_tmp * x);
   }
   else
   {
      fx = 2 * alpha_tmp * x - alpha_tmp * alpha_tmp * x * x;
   }

   return fx;
}


template<typename MatrixType>
void Purification_sp2acc<MatrixType>::apply_poly_to_eigs(const int it, real& homo, real& lumo)
{
   assert(it >= 0);
   if (it == 0)
   {
      return;
   }

   int  poly;
   real alpha_tmp;
   get_poly(it, poly, alpha_tmp);

   if (poly == 1)
   {
      homo = 2 * alpha_tmp * homo - alpha_tmp * alpha_tmp * homo * homo;
      lumo = (1 - alpha_tmp + alpha_tmp * lumo) * (1 - alpha_tmp + alpha_tmp * lumo);
   }
   else
   {
      homo = (1 - alpha_tmp + alpha_tmp * homo) * (1 - alpha_tmp + alpha_tmp * homo);
      lumo = 2 * alpha_tmp * lumo - alpha_tmp * alpha_tmp * lumo * lumo;
   }
}


template<typename MatrixType>
typename Purification_sp2acc<MatrixType>::real
Purification_sp2acc<MatrixType>::compute_derivative(const int it, real x, real& DDf)
{
   assert(it > 0);

   real Df;
   real temp, a, b;
   int  poly;
   real alpha;

   a   = x;
   Df  = 1;
   DDf = -1; // TODO

   for (int i = 1; i <= it; i++)
   {
      temp = a;

      get_poly(i, poly, alpha);

      if (poly == 1)
      {
         a = ((1 - alpha) + alpha * temp) * ((1 - alpha) + alpha * temp);
         b = 2 * alpha * ((1 - alpha) + alpha * temp);
      }
      else
      {
         a = 2 * alpha * temp - alpha * alpha * temp * temp;
         b = 2 * alpha - 2 * alpha * alpha * temp;
      }
      Df *= b;
   }

   return Df;
}


#endif //HEADER_PURIFICATION_SP2ACC
