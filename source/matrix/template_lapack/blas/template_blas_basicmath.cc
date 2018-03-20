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
 
 /* This file belongs to the template_lapack part of the Ergo source 
  * code. The source files in the template_lapack directory are modified
  * versions of files originally distributed as CLAPACK, see the
  * Copyright/license notice in the file template_lapack/COPYING.
  */
 

/* We need to include config.h to get macros PRECISION_QUAD_FLT128 and HAVE_SQRTL etc. */
#include "config.h"

#ifdef PRECISION_QUAD_FLT128
#include <quadmath.h>
#endif

#include <cmath>
#include <math.h>
#include <stdio.h>

#include "template_blas_basicmath.h"

/* fabs function */

#ifdef HAVE_FABSF
template<> float template_blas_fabs<float>(float x) { return fabsf(x); }
#else
template<> float template_blas_fabs<float>(float x) { return fabs(x); }
#endif

template<> double template_blas_fabs<double>(double x) { return fabs(x); }

#ifdef HAVE_FABSL
template<> long double template_blas_fabs<long double>(long double x) { return fabsl(x); }
#else
template<> long double template_blas_fabs<long double>(long double x) { return fabs(x); }
#endif

#ifdef PRECISION_QUAD_FLT128
#ifdef HAVE_FABSQ
template<> __float128 template_blas_fabs<__float128>(__float128 x) { return fabsq(x); }
#else
template<> __float128 template_blas_fabs<__float128>(__float128 x) { return fabs(x); }
#endif
#endif



/* sqrt function */

#ifdef HAVE_SQRTF
template<> float template_blas_sqrt<float>(float x) { return sqrtf(x); }
#else
template<> float template_blas_sqrt<float>(float x) { return sqrt(x); }
#endif

template<> double template_blas_sqrt<double>(double x) { return sqrt(x); }

#ifdef HAVE_SQRTL
template<> long double template_blas_sqrt<long double>(long double x) { return sqrtl(x); }
#else
template<> long double template_blas_sqrt<long double>(long double x) { return sqrt(x); }
#endif

#ifdef PRECISION_QUAD_FLT128
#ifdef HAVE_SQRTQ
template<> __float128 template_blas_sqrt<__float128>(__float128 x) { return sqrtq(x); }
#else
template<> __float128 template_blas_sqrt<__float128>(__float128 x) { return sqrt(x); }
#endif
#endif



/* exp function */

#ifdef HAVE_EXPF
template<> float template_blas_exp<float>(float x) { return expf(x); }
#else
template<> float template_blas_exp<float>(float x) { return exp(x); }
#endif

template<> double template_blas_exp<double>(double x) { return exp(x); }

#ifdef HAVE_EXPL
template<> long double template_blas_exp<long double>(long double x) { return expl(x); }
#else
template<> long double template_blas_exp<long double>(long double x) { return exp(x); }
#endif

#ifdef PRECISION_QUAD_FLT128
#ifdef HAVE_EXPQ
template<> __float128 template_blas_exp<__float128>(__float128 x) { return expq(x); }
#else
template<> __float128 template_blas_exp<__float128>(__float128 x) { return exp(x); }
#endif
#endif



/* log function */

#ifdef HAVE_LOGF
template<> float template_blas_log<float>(float x) { return logf(x); }
#else
template<> float template_blas_log<float>(float x) { return log(x); }
#endif

template<> double template_blas_log<double>(double x) { return log(x); }

#ifdef HAVE_LOGL
template<> long double template_blas_log<long double>(long double x) { return logl(x); }
#else
template<> long double template_blas_log<long double>(long double x) { return log(x); }
#endif

#ifdef PRECISION_QUAD_FLT128
#ifdef HAVE_LOGQ
template<> __float128 template_blas_log<__float128>(__float128 x) { return logq(x); }
#else
template<> __float128 template_blas_log<__float128>(__float128 x) { return log(x); }
#endif
#endif



/* log10 function */

#ifdef HAVE_LOG10F
template<> float template_blas_log10<float>(float x) { return log10f(x); }
#else
template<> float template_blas_log10<float>(float x) { return log10(x); }
#endif

template<> double template_blas_log10<double>(double x) { return log10(x); }

#ifdef HAVE_LOG10L
template<> long double template_blas_log10<long double>(long double x) { return log10l(x); }
#else
template<> long double template_blas_log10<long double>(long double x) { return log10(x); }
#endif

#ifdef PRECISION_QUAD_FLT128
#ifdef HAVE_LOG10Q
template<> __float128 template_blas_log10<__float128>(__float128 x) { return log10q(x); }
#else
template<> __float128 template_blas_log10<__float128>(__float128 x) { return log10(x); }
#endif
#endif



/* error function erf */

#ifdef HAVE_ERFF
template<> float template_blas_erf<float>(float x) { return erff(x); }
#else
template<> float template_blas_erf<float>(float x) { return erf(x); }
#endif

template<> double template_blas_erf<double>(double x) { return erf(x); }

#ifdef HAVE_ERFL
template<> long double template_blas_erf<long double>(long double x) { return erfl(x); }
#else
template<> long double template_blas_erf<long double>(long double x) { return erf(x); }
#endif

#ifdef PRECISION_QUAD_FLT128
#ifdef HAVE_ERFQ
template<> __float128 template_blas_erf<__float128>(__float128 x) { return erfq(x); }
#else
template<> __float128 template_blas_erf<__float128>(__float128 x) { return erf(x); }
#endif
#endif



/* complementary error function erfc */

#ifdef HAVE_ERFCF
template<> float template_blas_erfc<float>(float x) { return erfcf(x); }
#else
template<> float template_blas_erfc<float>(float x) { return erfc(x); }
#endif

template<> double template_blas_erfc<double>(double x) { return erfc(x); }

#ifdef HAVE_ERFCL
template<> long double template_blas_erfc<long double>(long double x) { return erfcl(x); }
#else
template<> long double template_blas_erfc<long double>(long double x) { return erfc(x); }
#endif

#ifdef PRECISION_QUAD_FLT128
#ifdef HAVE_ERFCQ
template<> __float128 template_blas_erfc<__float128>(__float128 x) { return erfcq(x); }
#else
template<> __float128 template_blas_erfc<__float128>(__float128 x) { return erfc(x); }
#endif
#endif



/* sine function sin */

#ifdef HAVE_SINF
template<> float template_blas_sin<float>(float x) { return sinf(x); }
#else
template<> float template_blas_sin<float>(float x) { return sin(x); }
#endif

template<> double template_blas_sin<double>(double x) { return sin(x); }

#ifdef HAVE_SINL
template<> long double template_blas_sin<long double>(long double x) { return sinl(x); }
#else
template<> long double template_blas_sin<long double>(long double x) { return sin(x); }
#endif

#ifdef PRECISION_QUAD_FLT128
#ifdef HAVE_SINQ
template<> __float128 template_blas_sin<__float128>(__float128 x) { return sinq(x); }
#else
template<> __float128 template_blas_sin<__float128>(__float128 x) { return sin(x); }
#endif
#endif



/* cosine function cos */

#ifdef HAVE_COSF
template<> float template_blas_cos<float>(float x) { return cosf(x); }
#else
template<> float template_blas_cos<float>(float x) { return cos(x); }
#endif

template<> double template_blas_cos<double>(double x) { return cos(x); }

#ifdef HAVE_COSL
template<> long double template_blas_cos<long double>(long double x) { return cosl(x); }
#else
template<> long double template_blas_cos<long double>(long double x) { return cos(x); }
#endif

#ifdef PRECISION_QUAD_FLT128
#ifdef HAVE_COSQ
template<> __float128 template_blas_cos<__float128>(__float128 x) { return cosq(x); }
#else
template<> __float128 template_blas_cos<__float128>(__float128 x) { return cos(x); }
#endif
#endif



/* power function pow */

#ifdef HAVE_POWF
template<> float template_blas_pow<float>(float x, float y) { return powf(x, y); }
#else
template<> float template_blas_pow<float>(float x, float y) { return pow(x, y); }
#endif

template<> double template_blas_pow<double>(double x, double y) { return pow(x, y); }

#ifdef HAVE_POWL
template<> long double template_blas_pow<long double>(long double x, long double y) { return powl(x, y); }
#else
template<> long double template_blas_pow<long double>(long double x, long double y) { return pow(x, y); }
#endif

#ifdef PRECISION_QUAD_FLT128
#ifdef HAVE_POWQ
template<> __float128 template_blas_pow<__float128>(__float128 x, __float128 y) { return powq(x, y); }
#else
template<> __float128 template_blas_pow<__float128>(__float128 x, __float128 y) { return pow(x, y); }
#endif
#endif




