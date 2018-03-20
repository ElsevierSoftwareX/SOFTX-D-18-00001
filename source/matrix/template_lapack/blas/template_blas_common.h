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
 

//#define sqrt HULAhulaHULAsqrt
//#define log  HULAhulaHULAlog

typedef int integer;
typedef bool logical;
typedef int ftnlen;
typedef char *address;

#define maxMACRO(a,b) (a >= b ? a : b)
#define minMACRO(a,b) (a <= b ? a : b)
#define absMACRO(x) (x >= 0 ? x : (-x))

#if 0
// Pawel suggests this variant:
template<typename A>
A maxMACRO(const A& a, const A&b) { return a >= b ? a : b; }
maxMACRO(f(a), b)
#endif

// include math.h to get sqrt etc
#include <math.h>


#include "template_blas_basicmath.h"


logical template_blas_lsame(const char *ca, const char *cb);
int template_blas_erbla(const char *srname, integer *info);
void template_blas_s_cat(char *lp, char *rpp[], ftnlen rnp[], ftnlen *np, ftnlen ll);

#include "template_blas_axpy.h"
#include "template_blas_scal.h"
#include "template_blas_dot.h"
#include "template_blas_spr.h"
#include "template_blas_spr2.h"
#include "template_blas_gemv.h"
#include "template_blas_gemm.h"
#include "template_blas_trmm.h"
#include "template_blas_trsm.h"
#include "template_blas_syrk.h"
#include "template_blas_syr2.h"
#include "template_blas_syr2k.h"
#include "template_blas_symv.h"
#include "template_blas_symm.h"
#include "template_blas_tpsv.h"
#include "template_blas_tpmv.h"
#include "template_blas_spmv.h"
#include "template_blas_trsv.h"
#include "template_blas_trmv.h"
#include "template_blas_swap.h"
#include "template_blas_nrm2.h"
#include "template_blas_copy.h"
#include "template_blas_ger.h"
#include "template_blas_idamax.h"
#include "template_blas_rot.h"
#include "template_blas_asum.h"
