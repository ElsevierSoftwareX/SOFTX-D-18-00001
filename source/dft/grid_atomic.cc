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

/** @file grid_atomic.cc Implements radial grid generators. */

#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include "grid_atomic.h"

/** vector of atoms' Bragg radii. It is indexed by atomic number. */
const real BraggRadii[] = {
  /* dummy         */
  0.75,
  /* H      He*   */
  0.35,  0.35,  
  /* Li     Be     B      C      N      O      F      Ne*  */
  1.45,  1.05,  0.85,  0.70,  0.65,  0.60,  0.50,  0.45,  
  /*Na     Mg     Al     Si     P      S      Cl     Ar*  */
  1.80,  1.50,  1.25,  1.10,  1.00,  1.00,  1.00,  1.00,  
  /* K      Ca     Sc     Ti     V      Cr     Mn     Fe     Co   */
  2.20,  1.80,  1.60,  1.40,  1.35,  1.40,  1.40,  1.40,  1.35,  
  /* Ni     Cu     Zn     Ga     Ge     As     Se     Br     Kr*  */
  1.35,  1.35,  1.35,  1.30,  1.25,  1.15,  1.15,  1.15,  1.10,  
  /* Rb     Sr     Y      Zr     Nb     Mo     Tc     Ru     Rh  */
  2.35,  2.00,  1.80,  1.55,  1.45,  1.45,  1.35,  1.30,  1.35,  
  /* Pd     Ag     Cd     In     Sn     Sb     Te     I      Xe*  */
  1.40,  1.60,  1.55,  1.55,  1.45,  1.45,  1.40,  1.40,  1.40,  
  /* Cs     Ba     La      */
  2.60,  2.15,  1.95,  
  /* Ce     Pr     Nd     Pm     Sm     Eu     Gd  */
  1.85,  1.85,  1.85,  1.85,  1.85,  1.85,  1.80,  
  /* Tb     Dy     Ho     Er     Tm     Yb     Lu  */
  1.75,  1.75,  1.75,  1.75,  1.75,  1.75,  1.75,  
  /* Hf     Ta     W      Re     Os     Ir     Pt     Au     Hg  */
  1.55,  1.45,  1.35,  1.30,  1.30,  1.35,  1.35,  1.35,  1.50,  
  /* Tl     Pb*    Bi     Po     At*    Rn*  */
  1.90,  1.75,  1.60,  1.90,  1.50,  1.50,  
  /* Fr*    Ra     Ac       */
  2.15,  2.15,  1.95,  
  /* rad(U): 1.75 --> 1.37D0  */
  /*Th     Pa     U      Np     Pu     Am     Cm*       */
  1.80,  1.80,  1.37,  1.75,  1.75,  1.75,  1.75,  
  /* Bk*    Cf*    Es*    Fm*    Md*    No*    Lw*  */
  1.75,  1.75,  1.75,  1.75,  1.75,  1.75,  1.75
};
/** Number of defined elements in BraggRadii array */
const unsigned BraggSize = sizeof(BraggRadii)/sizeof(BraggRadii[0]);

/* ===================================================================
 *             RADIAL QUADRATURES
 * the quadratore has to fill in grid->pnt with number of points
 * and set grid->rad.
 * =================================================================== */


/** Initializes RadialSchemeGC2 grid generator.  Determinates number
 * of radial points to be used for Gauss-Chebyshev quadrature of
 * second kind needed to integrate atom of specified Z number to
 * specified threshold thrl.
 */
void
RadialSchemeGC2::init(int myNumber, int Z, real thrl)
{
    static const int MIN_RAD_PT = 20;
    int ta, ri;

    if(Z<=2) ta=0;
    else if(Z<=10) ta=1;
    else if(Z<=18) ta=2;
    else if(Z<=36) ta=3;
    else if(Z<=54) ta=4;
    else if(Z<=86) ta=5;
    else ta=6;

    thrl = 1e-1*template_blas_sqrt(thrl); /* Fudge factor. */
    ri = int( -5.0*(3*template_blas_log10(thrl)-ta+8) );
    gridSize = ri>MIN_RAD_PT ? ri : MIN_RAD_PT;
}

/** Generates grid point positions and weights using Gauss-Chebyshev
   quadrature of second kind. The rad and wght arrays are filled
   in. */
void 
RadialSchemeGC2::generate(real *rad, real *wght)
{
    /* constants */
    static const real pi_2 = 2.0/M_PI;  
    static const real sfac = 2.0/3.0;
    const real rfac = 1.0/template_blas_log(static_cast<real>(2.0));
    real n_one, n_pi, wfac;
    /* variables */
    real x = 0.0, angl = 0.0, w = 0.0;
    int i;

    n_one = gridSize+1.0;
    n_pi  = M_PI/n_one;
    wfac = 16.0/(3*n_one);
    /* radial points */ 
    for (i=0; i<gridSize; i++) {
        real sinangl, sinangl2, r;
        x = (gridSize-1-2*i)/n_one;
        angl = n_pi*(i+1);
        sinangl = template_blas_sin(angl); 
        sinangl2 = sinangl*sinangl;
        x += pi_2*(1.0+sfac*sinangl2)*template_blas_cos(angl)*sinangl;
        r = rfac*template_blas_log( static_cast<real>(2.0/(1.0-x)) );
        w = wfac*sinangl2*sinangl2;
        wght[gridSize-i-1] = w*rfac/(1.0-x)*r*r;
        rad[gridSize-i-1] = r;
        /* transformation factor accumulated in weight */
    }
}

/** This quadrature follows [JCP 102, 346 (1995)].
    That is T2 quadrature with M4 mapping of r.
 */
void
RadialSchemeTurbo::init(int myNumber, int Z, real thrl)
{
  static const real zetas[] = 
    {/* H */ 0.8, /* He */ 0.9,
     /* Li */ 1.8, /* Be */ 1.4, /* B */ 1.3,  /* C */ 1.1,
     /* N */ 0.9,  /* O */  0.9, /* F */ 0.9,  /* Ne */ 0.9,
     /* Na */ 1.4, /* Mg */ 1.3, /* Al */ 1.3, /* Si */ 1.2,
     /* P */ 1.1,  /* S */  1.0, /* Cl */ 1.0, /* Ar */ 1.0,
     /* K */ 1.5,  /* Ca */ 1.4, /* Sc */ 1.3, /* Ti */ 1.2, /* V */ 1.2, 
     /* Cr */ 1.2, /* Mn */ 1.2, /* Fe */ 1.2, /* Co */ 1.2, /* Ni */ 1.1,
     /* Cu */ 1.1, /* Zn */ 1.1, /* Ga */ 1.1, /* Ge */ 1.0, /* As */ 0.9,
     /* Se */ 0.9, /* Br */ 0.9, /* Kr */ 0.9
    };

  int ta, accuracy_correction, z_correction;

  if(Z<=2) ta=0;
  else if(Z<=10) ta=1;
  else if(Z<=18) ta=2;
  else if(Z<=36) ta=3;
  else if(Z<=54) ta=4;
  else if(Z<=86) ta=5;
  else ta=6;

  /* thrl = 1e-5 maps to 0, 1e-13 -> 25, following Table III */
  accuracy_correction = int( (-template_blas_log10(thrl)-5.0)*3.0 );
  if(accuracy_correction<0) accuracy_correction = 0;
  z_correction = ta*5;

  static const int MIN_RAD_PT = 20;

  gridSize = MIN_RAD_PT + accuracy_correction + z_correction;
  zeta = Z >=1 && Z <= int(sizeof(zetas)/sizeof(zetas[0]))
    ?  zetas[Z-1] : 0.9;
}

/** Actual generation of the radial quadrature.
*/
void 
RadialSchemeTurbo::generate(real *rad, real *wght)
{
  const real piOverN  = M_PI/gridSize;
  static const real a = 1.0;

  const real rfac = zeta/M_LN2;

  /* radial points */ 
  for (int i=0; i<gridSize; i++) {
    real angle = (i+0.5)*piOverN;
    real x = template_blas_cos(angle);
    real s = template_blas_sin(angle);
    real w = piOverN * s;
    real aPlusX06 = template_blas_pow(a+x, (ergo_real)0.6);
    real logAPlus1Over1MinusX = template_blas_log( (a+1.0)/(1.0-x) );
    real r = rfac*aPlusX06*logAPlus1Over1MinusX;
    real rdiff = rfac*(aPlusX06/(1.0-x) +
                       0.6*logAPlus1Over1MinusX/template_blas_pow(a+x, (ergo_real)0.4));
    wght[gridSize-i-1] = w*rdiff*r*r;
    rad[gridSize-i-1] = r;
  }
}

/* gen_lmg_quad:
 *  As proposed by Roland Lindh, Per-Aake Malmqvist and Laura
 *  Gagliardi. */

RadialSchemeLMG::RadialSchemeLMG(const GridGenMolInfo& ggmi_)
  : RadialScheme("LMG scheme"), ggmi(ggmi_)
{
  ggmi.getExps(&maxL, &nucorb, &aa);
}

/**                                                                  
 *  diserr() provides grid spacing h for given angular momentum L and
 *  discretization error RD Based on eqs. (17) and (18) of R. Lindh,
 *  P.-Aa. Malmqvist and L. Gagliardi * "Molecular integrals by
 *  numerical quadrature", * Theor. Chem. Acc. 106 (2001) 178-187
 *
 * The array CF(4,L) contains coefficients of a 3rd order polynomial
 * fit to provide start values for the determination of H by a
 * Newton-Raphson search.
 *
 * Based on Fortran-77 code by T. Saue July 2002
 * This code DOES need double precision or precision dependent ACC factor.
 */
static real
diserr(int l, real rd)
{
    static const double ACC=1.0e-7;
    static const double cf[][4] = {
        { 0.91570e0,0.78806e-1,0.28056e-2,3.4197e-05 }, 
        { 0.74912e0,0.61502e-1,0.21558e-2,2.6100e-05 },
        { 0.65449e0,0.52322e-1,0.18217e-2,2.2004e-05 },
        { 0.59321e0,0.46769e-1,0.16261e-2,1.9649e-05 },
        { 0.55125e0,0.43269e-1,0.15084e-2,1.8270e-05 } }; 
    long_real fac, rdlog, res, x, htlog;
    int ifac, i, it, lm;

    fac  = 4*M_SQRT2;
    ifac = 1;
    for(i = 1; i<=l; i++) {
        fac   *= 2;
        ifac  = ifac*(2*i+1);
    }

    fac = fac/ifac;
    lm = l>4 ? 4 : l;
    rdlog = template_blas_log(rd);
    //res = polval(3,CF(1,LM),RDLOG);
    res = cf[lm][0]; x = rdlog;
    for(i=1; i<4; i++) {
        res += cf[lm][i]*x;
        x *= rdlog;
    }
    htlog = template_blas_log(res);
    // Newton-Raphson search
    for(it = 0; it<20; it++) {
        long_real u0, u1, f0, f1, dx;
        long_real pih  = M_PI/res;
        long_real pihl = pih;
        long_real piex = M_PI*pih*0.5;
        for(i = 0; i<l; i++)
            pihl = pihl*pih;
        u0   = fac*pihl*template_blas_exp(-piex);
        u1   = u0*((piex/res)-(l+1)/pih);
        f0   = template_blas_log(u0)-rdlog;
        f1   = res*u1/u0;
        dx = f0/f1;
        htlog = htlog - dx;
        res = template_blas_exp(htlog);
        if(template_blas_fabs(dx)<ACC) return res;
    }
    puts("diserr never reached"); return 0.1;
}

/** outerr() provides outer grid point for given angular momentum L
 * outer exponent AL and discretization error RD
 
 * Based on eq. (19) of R. Lindh, P.-Aa. Malmqvist and L. Gagliardi
 * "Molecular integrals by numerical quadrature",
 * Theor. Chem. Acc. 106 (2001) 178-187
                                                                    
 * The variable U = AL*R*R is found by a Newton-Raphson search.
 * Based on a F77 code by T. Saue July 2002
 */
static real
outerr(real al, int l, real rd)
{
    real tolen, fac, a;
    long_real aln, expl, u, rln;
    int  it;

    tolen = 2;
    fac = 1;
    for(it=1; it<=l; it++) {
        tolen *= 2;
        fac   = fac*(2*it+1);
    }
    expl = 0.5*(2*l+1);
    a = 2*fac/(tolen*M_2_SQRTPI);
    aln = template_blas_log(a);
    rln = template_blas_log(rd);
    u = 35.0;
    //Newton-Raphson search
    for(it = 0; it<8; it++) {
        long_real f0hln = aln+expl*template_blas_log(u)-u-rln;
        long_real f1hln = expl/u-1.0;
        long_real dx = f0hln/f1hln;
        u = u - dx;
        if(template_blas_fabs(dx)<1e-8) return template_blas_sqrt(u/al);
    }
    puts("outerr never reached");
    return 10.0;
}

/** Initializes the LMG radial grid generator for given atom charge
    and acceptable error threshold. */
void 
RadialSchemeLMG::init(int myNumber, int charge, real thrl)
{
    int  *lNucorb;
    real (*lAA)[2];
    real ah, rh;
    int l;


    /*
     *     Grid spacing to H and inner grid point to AH
     */
    lNucorb = nucorb+myNumber*maxL;
    lAA     = aa    +myNumber*maxL;
    h  = 1e30;
    ah = 0.0;
    for(l=0; l<maxL; l++) {
        if(lNucorb[l] > 0) {
            real htmp = diserr(l, thrl);
            if(htmp<h) h = htmp;
            //printf("Spacing for l=%d: %f\n", l, htmp);
        }
        if(lAA[l][1]>ah) ah = lAA[l][1];
        //printf("ah=%g %g\n", ah, aa[l][1]);
    }
    eph = template_blas_exp(h);
    ah *= 2;
    rl = (1.9+template_blas_log(thrl))/3.0 - 0.5*template_blas_log(ah);
    rl = template_blas_exp(rl);
    //printf("* Inner grid point: %10.5g ah=%f\n", rl, ah);
    //...  Outer point
    rh = 0.0;
    for(l=0; l<maxL; l++) {
        if(lAA[l][0]>0) {
            real rhtmp = outerr(2*lAA[l][0], l, thrl);
            if(rh<rhtmp) rh = rhtmp;
        }
    }
    grdc = rl/(eph-1.0);
    gridSize = (int) (template_blas_log(1.0 + rh/grdc)/h);
}

/** Generates grid point positions and associated weights using LMG
    method. */
void 
RadialSchemeLMG::generate(real *radposn, real *radwght)
{
  radposn[0] = rl;
  radwght[0]  = (rl+grdc)*rl*rl*h;
  for(int ip=1; ip<gridSize; ip++) {
    radposn[ip] = (radposn[ip-1]+grdc)*eph-grdc;
    radwght[ip] = (radposn[ip]  +grdc)*radposn[ip]*radposn[ip]*h;
  }
}


RadialSchemeLMG::~RadialSchemeLMG() 
{
    free(nucorb);
    free(aa);
}
