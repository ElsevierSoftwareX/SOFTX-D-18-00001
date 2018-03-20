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
 * Kohnâˆ’Sham Density Functional Theory Electronic Structure Calculations 
 * with Linearly Scaling Computational Time and Memory Usage,
 * Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek,
 * J. Chem. Theory Comput. 7, 340 (2011),
 * <http://dx.doi.org/10.1021/ct100611z>
 * 
 * For further information about Ergo, see <http://www.ergoscf.org>.
 */

/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/** @file functionals.h Functional library interface.
 The variables, structures and functions related to computation
   of functional and their derivatives.
   (c) Pawel Salek, pawsa@theochem.kth.se. 2001.07.13

   NOTE1: the derivatives are computed with respect to the density,
   and SQUARE of the density gradient. This is a choice. It stems from
   the fact that the factors involved in the derivative vector
   distribution depend on the square of the density gradient.

   NOTE2: C version is included once per file, Fortran version -
   multiple times.
*/
#ifndef _FUNCTIONALS_H_
#define _FUNCTIONALS_H_

/* We need to include config.h to get macro PRECISION_QUAD_FLT128 */
#include "config.h"

#ifdef PRECISION_QUAD_FLT128
#include <quadmath.h>
#endif

#include "realtype.h"

#ifdef __cplusplus
#define EXTERN_C extern "C"
#else
#define EXTERN_C
#endif

typedef ergo_real real;

#if defined(FUNC_PRECISION) && FUNC_PRECISION == 1

/* Use single-precision math functions */
#define ASINH asinhf
#define ATAN atanf
#define ERF erff
#define EXP expf
#define FABS fabsf
#define LOG logf
#define POW powf
#define SQRT sqrtf

#elif defined(FUNC_PRECISION) && FUNC_PRECISION == 2

/* Use double-precision math functions */
#define ASINH asinh
#define ATAN atan
#define ERF erf
#define EXP exp
#define FABS fabs
#define LOG log
#define POW pow
#define SQRT sqrt

#else

#ifdef PRECISION_QUAD_FLT128
/* Use quad-precision math functions */
#define ASINH asinhq
#define ATAN atanq
#define ERF erfq
#define EXP expq
#define FABS fabsq
#define LOG logq
#define POW powq
#define SQRT sqrtq
#else
/* Use "long double"-precision math functions */
#define ASINH asinhl
#define ATAN atanl
#define ERF erfl
#define EXP expl
#define FABS fabsl
#define LOG logl
#define POW powl
#define SQRT sqrtl
#endif

#endif

/* FirstDrv: matrix of first order derivatives with respect to two
 * parameters: density rho and SQUARE of the gradient of density grho.  
 * zeta_i = |\nabla\rho_i|²
 * mu     = |\nabla\rho_\alpha||\nabla\rho_\beta|
 */
typedef struct {
    real df1000;  /* d/drho F     */
    real df0100;
    real df0010;  /* d/zeta F     */
    real df0001;
    real df00001;
} FunFirstFuncDrv;

/* SecondFuncDrv: this structure is used by functional derivative
 * evaluation procedures. Do not include "triplet" transformation.
 */
typedef struct {
    real df1000;  /* d/drho_alpha F               */
    real df0100;  /* d/drho_beta F                */
    real df0010;  /* d/|zeta_alpha| F             */
    real df0001;  /* d/|zeta_beta| F              */
    real df00001;
    real df2000;  /* d/drho_alpha^2 F             */
    real df1100;  /* d/(drho_alpha drho_beta) F   */
    real df1010;  /* d/drho_alpha d/dzeta_alpha F */
    real df1001;  /* d/drho_alpha d/dzeta_beta F  */
    real df10001;
    real df0200;  /* d/drho_beta^2 F              */
    real df0110;  /* d/drho_beta d/dzeta_alpha F  */ 
    real df0101;  /* d/drho_beta d/dzeta_beta F   */
    real df01001;
    real df0020;  /* d/dzeta_alpha^2 F          */
    real df0011;  /* d2/dzeta_a zeta_b F          */
    real df00101;
    real df0002;  /* d/dzeta_beta^2 F             */
    real df00011;
    real df00002;
} FunSecondFuncDrv;


/* ThirdFuncDrv: matrix of third derivatives with respect to five
   parameters: density rho_alpha and SQUARE of the density gradient
   zeta.  and mu.
*/

typedef struct {
    real df1000;   /* d/drho F          */
    real df0100; 
    real df0010;   /* d/|zeta| F        */
    real df0001;
    real df00001;
     
    real df2000;  /* d/drho_alpha^2 F             */
    real df1100;  /* d/(drho_alpha drho_beta) F   */
    real df1010;  /* d/drho_alpha d/dzeta_alpha F */
    real df1001;  /* d/drho_alpha d/dzeta_beta F  */
    real df10001;
    real df0200;  /* d/drho_beta^2 F              */
    real df0110;  /* d/drho_beta d/dzeta_alpha F  */ 
    real df0101;  /* d/drho_beta d/dzeta_beta F   */
    real df01001;
    real df0020;  /* d/dzeta_alpha^2 F          */
    real df0011;  /* d2/dzeta_a zeta_b F          */
    real df00101;
    real df0002;  /* d/dzeta_beta^2 F             */
    real df00011;
    real df00002;

    real df3000;
    real df2100;
    real df2010;
    real df2001;
    real df20001;
    real df1200;
    real df1110;
    real df1101;
    real df11001;
    real df1020;
    real df1011;
    real df10101;
    real df1002;
    real df10011;
    real df10002;
    real df0300;
    real df0210;
    real df0201;
    real df02001;
    real df0120;
    real df0111;
    real df01101;
    real df0102;
    real df01011;
    real df01002;
    real df0030;
    real df0021;
    real df00201;
    real df0012;
    real df00111;
    real df00102;
    real df0003;
    real df00021;
    real df00012;
    real df00003;
} FunThirdFuncDrv;


typedef struct {
  
  /* First order derivatives with respect to all 5 variables */

    real df1000;
    real df0100;
    real df0010;
    real df0001;
    real df00001;
	
  /* Second order mixed derivatives with respect to all 5 variables */

    real df2000;
    real df1100;
    real df1010;
    real df1001;
    real df10001;
    real df0200;
    real df0110;
    real df0101;
    real df01001;
    real df0020;
    real df0011;
    real df00101;
    real df0002;
    real df00011;
    real df00002;
	
  /* Third order mixed derivatives with respect to all 5 variables */

    real df3000;
    real df2100;
    real df2010;
    real df2001;
    real df20001;
    real df1200;
    real df1110;
    real df1101;
    real df11001;
    real df1020;
    real df1011;
    real df10101;
    real df1002;
    real df10011;
    real df10002;
    real df0300;
    real df0210;
    real df0201;
    real df02001;
    real df0120;
    real df0111;
    real df01101;
    real df0102;
    real df01011;
    real df01002;
    real df0030;
    real df0021;
    real df00201;
    real df0012;
    real df00111;
    real df00102;
    real df0003;
    real df00021;
    real df00012;
    real df00003;
    
  /* Fourth order mixed derivatives with respect to all 5 variables */

    real df4000;
    real df3100;
    real df3010;
    real df3001;
    real df30001;
    real df2200;
    real df2110;
    real df2101;
    real df21001;
    real df2020;
    real df2011;
    real df20101;
    real df2002;
    real df20011;
    real df20002;
    real df1300;
    real df1210;
    real df1201;
    real df12001;
    real df1120;
    real df1111;
    real df11101;
    real df1102;
    real df11011;
    real df11002;
    real df1030;
    real df1021;
    real df10201;
    real df1012;
    real df10111;
    real df10102;
    real df1003;
    real df10021;
    real df10012;
    real df10003;
    real df0400;
    real df0310;
    real df0301;
    real df03001;
    real df0220;
    real df0211;
    real df02101;
    real df0202;
    real df02011;
    real df02002;
    real df0130;
    real df0121;
    real df01201;
    real df0112;
    real df01111;
    real df01102;
    real df0103;
    real df01021;
    real df01012;
    real df01003;
    real df0040;
    real df0031;
    real df00301;
    real df0022;
    real df00211;
    real df00202;
    real df0013;
    real df00121;
    real df00112;
    real df00103;
    real df0004;
    real df00031;
    real df00022;
    real df00013;
    real df00004;    
} FunFourthFuncDrv;


typedef struct Functional_ Functional;

enum FunError { FUN_OK, FUN_UNKNOWN, FUN_CONF_ERROR };
EXTERN_C enum FunError fun_select_by_name(const char *conf_string);
extern Functional *selected_func;
extern int (*fun_printf)(const char *fmt, ...);
extern void (*fun_set_hf_weight)(real w);
extern real (*fun_get_hf_weight)(void);
extern void (*fun_set_cam_param)(real w, real b);
int dft_get_cam_param(real *alpha, real *beta, real *mu);

/* FunDensProp structure contains properties of the density that are
   needed for functional evaluation and possibly other purposes.
*/
typedef struct FunDensProp_ {
    real rhoa,  rhob;
    real grada, gradb; /* norms of the density gradient, not squares */
    real gradab;       /* scalar product of grada and gradb */
    /* real current[3] or something may come in the future :-) */
} FunDensProp;

/* EnergyFunc: the function returning the energy for given densities
   and gradients. Note that some functionals(like LYP) depend explicitely
   on separately alpha and beta densities
*/
typedef int (*IsGGAFunc)(void);
typedef int (*ReadInputFunc)(const char* conf_string);
typedef void (*ReportFunc)(void);
typedef real (*EnergyFunc)(const FunDensProp* dens_prop);
typedef void (*FirstOrderFun)(FunFirstFuncDrv *ds, real factor,
                              const FunDensProp* dns_prp);

typedef void (*SecondOrderFun)(FunSecondFuncDrv *ds, real factor,
                               const FunDensProp* dens_prop);

typedef void (*ThirdOrderFun)(FunThirdFuncDrv *ds, real factor,
                              const FunDensProp* dens_prop);
typedef void (*FourthOrderFun)(FunFourthFuncDrv *ds, real factor,
                               const FunDensProp *dens_prop);

struct Functional_ {
    const char* name; /* descriptive functional name (usually 5 characters) */
    IsGGAFunc       is_gga;
    ReadInputFunc   read;
    ReportFunc      report;
   /* Only unrestricted implementations are needed. A benchmark for
     * a CO molecule with 28 basis function reveals a 4% time difference.
     * This difference will only decrease for larger systems. */
    EnergyFunc      func;
    FirstOrderFun   first;
    SecondOrderFun  second;
    ThirdOrderFun   third;
    FourthOrderFun  fourth;
};

EXTERN_C void drv1_clear(FunFirstFuncDrv* gga);  /* set all components to 0 */
EXTERN_C void drv2_clear(FunSecondFuncDrv* gga); /* set all components to 0 */
EXTERN_C void drv3_clear(FunThirdFuncDrv* gga);  /* set all components to 0 */
EXTERN_C void drv4_clear(FunFourthFuncDrv* gga); /* set all components to 0 */

/* The list of functionals */
/* sorted list of generic functionals */
extern Functional BeckeFunctional;
extern Functional KTFunctional;
extern Functional LB94Functional;
extern Functional LYPFunctional;
extern Functional OPTXFunctional;
extern Functional P86cFunctional;
extern Functional PW86xFunctional;
extern Functional Pw91xFunctional;
extern Functional Pw91cFunctional;
extern Functional PW92cFunctional;
extern Functional PZ81Functional;
extern Functional PbecFunctional;
extern Functional PbexFunctional;
extern Functional SlaterFunctional;
extern Functional VWN3Functional;
extern Functional VWN5Functional;
extern Functional VWNIFunctional;
extern Functional VWN3IFunctional;
extern Functional VWNFunctional;
extern Functional XAlphaFunctional;

/* sorted list of mixed functionals */
extern Functional B3LYPFunctional;
extern Functional B3LYPGaussFunctional;
extern Functional B3P86Functional;
extern Functional B3P86GFunctional;
extern Functional B3PW91Functional;
extern Functional BHandHFunctional;
extern Functional BHandHLYPFunctional;
extern Functional BLYPFunctional;
extern Functional BP86Functional;
extern Functional BPW91Functional;
extern Functional Camb3lypFunctional;
extern Functional CamFunctional;
extern Functional CombineFunctional;
extern Functional GGAKeyFunctional;
extern Functional HseFunctional;
extern Functional KT1Functional;
extern Functional KT2Functional;
extern Functional KT3Functional;
extern Functional LDAFunctional;
extern Functional OLYPFunctional;
extern Functional PBE0Functional;
extern Functional PBEFunctional;
extern Functional SVWN3Functional;
extern Functional SVWN5Functional;

/* the list of the functionals */
extern Functional* available_functionals[];

extern int fun_true(void);
extern int fun_false(void);
/* fortran (and not only) functional stub routines */
EXTERN_C void dftlistfuncs_(void);
EXTERN_C void dftlistfuncs_using_printf_(void);
int dft_isgga_(void);
int dft_isgga__(void);

EXTERN_C void dftreport_(void);
EXTERN_C int fun_get_cam_param(real *alpha, real *beta, real *mu);

#endif /* _FUNCTIONALS_H_ */
