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

/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/** @file fun-slater.c
   Implementation of Slater functional and its derivatives .
   (c), Pawel Salek, pawsa@theochem.kth.se, aug 2001
   Z. Rinkevicius adapted for open shell systems: energy, first derivatives.
   NOTE:
   this file may seem unnecessarily complex but the structure really pays off
   when implementing multiple functionals depending on different parameters.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int slater_isgga(void) { return 0; }
static int slater_read(const char* conf_line);
static real slater_energy(const FunDensProp* dp);
static void slater_first(FunFirstFuncDrv *ds,   real fac, const FunDensProp*);
static void slater_second(FunSecondFuncDrv *ds, real fac, const FunDensProp*);
static void slater_third(FunThirdFuncDrv *ds,   real fac, const FunDensProp*);
static void slater_fourth(FunFourthFuncDrv *ds, real fac, const FunDensProp*);

Functional SlaterFunctional = {
  "Slater",       /* name */
  slater_isgga,   /* gga-corrected */
  slater_read, 
  NULL,
  slater_energy, 
  slater_first,
  slater_second,
  slater_third,
  slater_fourth
};

/* IMPLEMENTATION PART */
static int
slater_read(const char* conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}

/* SLATER_THRESHOLD Only to avoid numerical problems due to raising 0
 * to a fractional power. */
static const real SLATER_THRESHOLD = 1e-20;
static real
slater_energy(const FunDensProp* dp)
{
  real ea = 0.0, eb = 0.0; 
  const real PREF= -3.0/4.0*POW(6/M_PI, 1.0/3.0);
  if (dp->rhoa >SLATER_THRESHOLD) {
    real powValue = POW(dp->rhoa,4.0/3.0);
    /* ELIAS NOTE 2011-05-03: When using long double precision it
       seems like the function POW (powl) sometimes gives extremely
       wrong results. For example:
       dp->rhoa = 2.058e-20
       POW(dp->rhoa,4.0/3.0) =    0.3969
       That example occurred when running the test_dft_hicu test case
       on a Intel Core i5 CPU, with gcc version 4.5.1.
       The if statement below was added to expose such powl() problems.
     */
    /* ELIAS NOTE 2011-06-02: This error is probably caused by a bug
       in the glibc powl() implementation, see bug report at
       http://gcc.gnu.org/bugzilla/show_bug.cgi?id=49031 and
       http://sourceware.org/bugzilla/show_bug.cgi?id=12775
     */
    if(dp->rhoa < 1e-5 && powValue > 1e-5) {
      printf("Error in slater_energy: POW gives crazy result. "
	     "This error is probably caused by a bug in the glibc "
	     "powl() implementation, see bug report at "
	     "http://gcc.gnu.org/bugzilla/show_bug.cgi?id=49031 "
	     "and http://sourceware.org/bugzilla/show_bug.cgi?id=12775\n");
      printf("dp->rhoa = %9.4g\n", (double)dp->rhoa);
      printf("POW(dp->rhoa,4.0/3.0) = %9.4g\n", (double)powValue);
      exit(-1);
    }
    ea= PREF*powValue;
  }
  if (dp->rhob >SLATER_THRESHOLD)
      eb= PREF*POW(dp->rhob,4.0/3.0);   
  return ea+eb;  
}

static void
slater_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
  if (dp->rhoa>SLATER_THRESHOLD)
     ds->df1000 += -POW(6.0/M_PI*dp->rhoa, 1.0/3.0)*factor;
  if (dp->rhob>SLATER_THRESHOLD)
     ds->df0100 += -POW(6.0/M_PI*dp->rhob, 1.0/3.0)*factor;
}
static void
slater_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
  const real PREF = POW(6.0/M_PI, 1.0/3.0);
  if (dp->rhoa>SLATER_THRESHOLD) {
    ds->df1000 += -PREF*POW(dp->rhoa,  1.0/3.0)*factor;
    ds->df2000 += -PREF*POW(dp->rhoa, -2.0/3.0)/3*factor;
  }
  if (dp->rhob>SLATER_THRESHOLD) {
    ds->df0100 += -PREF*POW(dp->rhob,  1.0/3.0)*factor;
    ds->df0200 += -PREF*POW(dp->rhob, -2.0/3.0)/3*factor;
  }
}

/* slater_third:
   Slater functional derivatives.
*/
static void
slater_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
  const real PREF = POW(6.0/M_PI, 1.0/3.0);
  if (dp->rhoa>SLATER_THRESHOLD) {
    ds->df1000 += -PREF*POW(dp->rhoa,  1.0/3.0)*factor;
    ds->df2000 += -PREF*POW(dp->rhoa, -2.0/3.0)/3*factor;
    ds->df3000 +=  PREF*POW(dp->rhoa, -5.0/3.0)*2.0/9.0*factor;
  }
  if (dp->rhob>SLATER_THRESHOLD) {
    ds->df0100 += -PREF*POW(dp->rhob,  1.0/3.0)*factor;
    ds->df0200 += -PREF*POW(dp->rhob, -2.0/3.0)/3*factor;
    ds->df0300 +=  PREF*POW(dp->rhob, -5.0/3.0)*2.0/9.0*factor;
  }
}

/* slater_fourth:
   Dirac functional fourth derivatives.
   by B. Jansik
*/
static void
slater_fourth(FunFourthFuncDrv *ds, real factor, const FunDensProp *dp)
{
    const real PREF = (-3.0/4.0)*POW(6.0/M_PI, 1.0/3.0);/* Dirac G prefactor */
    const real DPREF = 40.0/81.0;	  /* Prefactor from 4th derivative */
    const real JPREF = DPREF*PREF*factor; /* Joined prefactor */
    const real ROEXP = -8.0/3.0;	  /* Exponent on density (from 4th deriv.) */
    FunThirdFuncDrv ds_third;
   
   
    /* set up lower order derivatives */
    /* dirac_third contain third and also lower order derivatives	 */

   drv3_clear(&ds_third);
   slater_third(&ds_third, factor, dp);
   
   ds->df1000 += ds_third.df1000;
   ds->df2000 += ds_third.df2000;
   ds->df3000 += ds_third.df3000;

   ds->df0100 += ds_third.df0100;
   ds->df0200 += ds_third.df0200;
   ds->df0300 += ds_third.df0300;

   if (dp->rhoa > SLATER_THRESHOLD)
     ds->df4000 += JPREF*POW(dp->rhoa, ROEXP);
   
   if (dp->rhob > SLATER_THRESHOLD)
     ds->df0400 += JPREF*POW(dp->rhob, ROEXP);
}
