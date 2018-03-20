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

/** @file fun-pw86x.c
   The PW86 exchange functional and its derivative.
   Contributed by Olav Fossgaard, olav@chem.uit.no, May 2002
 
   Reference: Phys. Rev. B 33. 8800 (1986)
*/

#include <math.h>
#include <stdio.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int pw86x_isgga(void) { return 1; }
static int pw86x_read(const char* conf_line);
static real pw86x_energy(const FunDensProp* dp);
static void pw86x_first (FunFirstFuncDrv *ds,  real factor, const FunDensProp*dp);

Functional PW86xFunctional = {
  "PW86x",       /* name */
  pw86x_isgga,   /* gga-corrected */
  pw86x_read, 
  NULL,
  pw86x_energy, 
  pw86x_first,
  NULL,
  NULL
};

/* IMPLEMENTATION PART */
static int
pw86x_read(const char* conf_line)
{
  fun_set_hf_weight(0);
  return 1;
}

static real
pw86x_energy(const FunDensProp* dp)
{
/* Use density functional form. In case of spin polarization,
   this function will have to be called twice with arguments
   rho=2rhoa and rho=2rhob, respectively. The total energy is then
   half the sum of the returned values.
*/
  const real a = 1.0;
  const real b = 1.296;
  const real c = 14.0;
  const real d = 0.20;
/* Closed shell (See eq. (25) in reference) */
  real rho = dp->rhoa+dp->rhob, grad = dp->grada+dp->gradb;

  const real Ax = -POW(3.0/M_PI,1.0/3.0)*3.0/4.0;
  const real kf = POW(3.0*POW(M_PI,2.0)*rho,1.0/3.0);
  real s = grad/(2.0*kf*rho);
  real F = POW(a+b*POW(s,2.0)+c*POW(s,4.0)+d*POW(s,6.0),1.0/15.0); 
  return Ax*POW(rho,4.0/3.0)*F;
}

static void
pw86x_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
/* The energy expression is the integral int(Ax*rho**(4/3)*F). We first 
   calculate d(F)/d(rho) and d(F)/d(grad_rho) and differentiate the
   product in the last step.
*/
  const real a = 1.0;
  const real b = 1.296;
  const real c = 14.0;
  const real d = 0.20;
/* Closed shell (See eq. (25) in reference) */
  real rho = dp->rhoa+dp->rhob, grad = dp->grada+dp->gradb;

  const real Ax= -POW(3.0/M_PI,1.0/3.0)*3.0/4.0;
  const real kf= POW(3.0*M_PI*M_PI*rho,1.0/3.0);
  real  s = grad/(2.0*kf*rho);
  real  F = POW(a+b*POW(s,2.0)+c*POW(s,4.0)+d*POW(s,6.0),1.0/15.0);

  real F1 = 1.0/15.0*POW(a+b*POW(s,2.0)+c*POW(s,4.0)+d*POW(s,6.0),-14.0/15.0)
            *(2.0*b*s+4.0*c*POW(s,3.0)+6.0*d*POW(s,5.0)); /* dF/ds */

  real s1 = -4.0*s/(3.0*rho); /* ds/d(rho) */
  real s2 = 1.0/(2.0*kf*rho); /* d(s)/d(grad) */
  real G1 = F1*s1; /* dF/d(rho) */ 
  real G2 = F1*s2; /* dF/d(grad) */
   
  ds->df1000 += Ax*((4.0/3.0)*POW(rho,1.0/3.0)*F + POW(rho,4.0/3.0)*G1 )*factor;
  ds->df0010 += Ax*POW(rho,4.0/3.0)*G2*factor;
}
