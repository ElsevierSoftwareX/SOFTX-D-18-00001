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

/** @file fun-lb94.c LB94 implementation.
   Implementation of Exchange-correlation potential with correct 
   asymptotic behavior by R. van Leeuwen and E. J. Baerends:

   [ van Leeuwen and EJ Baerends, Phys Rev A 49, 2421 (1994)]
   See also comments in Gisbergen et al, JCP 105(8) 3142.

   (c) P. Salek, oct 2003 - the working implementation.
*/

#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int  lb94_isgga(void) { return 1; }
static int  lb94_read(const char* conf_line);
static real lb94_energy(const FunDensProp* dens_prop);
static void lb94_first(FunFirstFuncDrv *ds, real factor, 
                       const FunDensProp* dens_prop);
static void lb94_second(FunSecondFuncDrv *ds, real factor,
                        const FunDensProp* dens_prop);

static void lb94_third(FunThirdFuncDrv *ds, real factor,
                       const FunDensProp* dens_prop);

#ifdef FOURTH_ORDER_DERIVATIVES
static void lb94_fourth(FourthFuncDrv *ds, real factor,
                        const FunDensProp* dens_prop);
#endif

Functional LB94Functional = {"LB94",      /* name */
                             lb94_isgga,  /* gga-corrected */
                             lb94_read,   /* set common blocks */
                             NULL,         /* reporter */
                             lb94_energy, 
                             lb94_first,
                             lb94_second,
                             lb94_third
#ifdef FOURTH_ORDER_DERIVATIVES
                             ,lb94_fourth,
#endif
};

/* IMPLEMENTATION PART */

static int
lb94_read(const char* conf_line)
{
    fun_set_hf_weight(0.0);
    return 1;
}

static const real BETA = 0.05;

static real
lb94_energy(const FunDensProp* dp)
{
    return SlaterFunctional.func(dp)+VWNFunctional.func(dp);
}

static void
lb94_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real rho    = dp->rhoa + dp->rhob;
    real rho13 = POW(rho, 1.0/3.0);
    real grad = dp->grada + dp->gradb;
    real rho43=rho*rho13;
    real scaled_grad, sg2, vx;
    scaled_grad = grad/(rho43>1e-13 ? rho43 : 1e-13);
    sg2   = scaled_grad*scaled_grad;

    vx = -BETA*rho13*sg2/
        (1+3*BETA*scaled_grad*ASINH(scaled_grad));

    ds->df1000 += vx*factor;
    ds->df0100 += vx*factor;

    SlaterFunctional.first(ds, factor, dp);
    VWNFunctional.first(ds,   factor, dp);
}

static void
lb94_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
/* according to the authors, it is equivalent to ALDA for LR and higher.
 * See comments in Gisbergen et al, JCP 105(8) 3142. */
#if 0
    real
        x1 = POW(dp->gradb,2.0)+2.0*dp->gradab+POW(dp->grada,2.0),
        x2 = dp->rhob+dp->rhoa,
        x3 = 1/POW(x2,2.333333333333334),
        x4 = sqrt(x1),
        x5 = 1/POW(x2,1.333333333333333),
        x6 = asinh(x4*x5),
        x7 = 0.15*x4*x5*x6+1.0,
        x8 = 1/x7,
        x9 = -0.05*x1*x3*x8,
        x10 = 1/POW(x2,2.666666666666667),
        x11 = 1/sqrt(x1*x10+1.0),
        x12 = 1/POW(x7,2.0),
        x13 = 0.11666666666667*x1*x8/POW(x2,3.333333333333334)
        +0.05*x1*x12*x3*(-0.2*x4*x3*x6-0.2*x1*x11/POW(x2,3.666666666666667)),
        x14 = 1/x4,
        x15 = 0.05*x1*x3*(0.15*dp->grada*x14*x5*x6+0.15*dp->grada*x10*x11)*x12
        -0.1*dp->grada*x3*x8,
        x16 = 0.05*x1*x3*(0.15*dp->gradb*x14*x5*x6+0.15*dp->gradb*x10*x11)*x12
        -0.1*dp->gradb*x3*x8;
    
    ds->df1000 += x9*factor;
    ds->df0100 += x9*factor;

    ds->df2000 += x13*factor;
    ds->df0200 += x13*factor;
    ds->df1100 += x13*factor;

    ds->df1010 += x15*factor;
    ds->df1001 += x16*factor;
    ds->df10001 += (0.05*x1*x3*(0.15*x14*x5*x6+0.15*x10*x11)*x12-0.1*x3*x8)
        *factor;
    ds->df0110 += x15*factor;
    ds->df0101 += x16*factor;
#endif    
    SlaterFunctional.second(ds, factor, dp);
    VWNFunctional.second(ds,   factor, dp);
}

 
static void
lb94_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    SlaterFunctional.third(ds, factor, dp);
    VWNFunctional.third(ds,   factor, dp);
}

#ifdef FOURTH_ORDER_DERIVATIVES  
static void
lb94_fourth(FourthFuncDrv *ds, real factor, const FunDensProp* dp)
{
    SlaterFunctional.fourth(ds, factor, dp);
    VWNFunctional.fourth(ds,   factor, dp);
}
#endif
