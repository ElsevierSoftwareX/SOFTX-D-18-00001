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

/** @file integrals_hermite.cc

    @brief Code for computation of Coulomb integrals of Hermite
    Gaussians, using the the McMurchie-Davidson scheme as described in
    the book "Molecular electronic-structure theory" by Trygve
    Helgaker, Poul Jorgensen, and Jeppe Olsen.

    @author: Elias Rudberg <em>responsible</em>
*/

#include "integrals_hermite.h"
#include "boysfunction.h"
#include "mat_gblas.h"
#include <cmath>
#include <stdio.h>

/* ELIAS NOTE 2014-07-12: The variables R1000 R1001 etc in this file
   were renamed to R_val_1000 etc to avoid R3000/R4000 identifiers,
   since those identifiers apparently cannot be used on Mips/Mips64
   architectures. The variables were renamed using the standalone
   utility program standalone/rename_Rxxxx_in_code.cc.  */

int
get_related_integrals_hermite(const IntegralInfo & integralInfo,
			      const JK::ExchWeights & paramsCAM,
			      int n1max, int noOfMonomials_1,
			      int n2max, int noOfMonomials_2,
			      ergo_real dx0, 
			      ergo_real dx1, 
			      ergo_real dx2, 
			      ergo_real alpha0,
			      ergo_real resultPreFactor,
			      ergo_real* primitiveIntegralList)
{
  
  int Nmax = n1max + n2max;
  ergo_real R_12_squared = dx0*dx0 + dx1*dx1 + dx2*dx2;
  
  ergo_real BoysList[Nmax+1];

  
  /* Get Boys function values and store them in BoysList.
     NOTE: If CAM params are used, the values in BoysList are
     not simply Boys function values but are modified 
     according to the CAM params.  */

  if(paramsCAM.computeRangeSeparatedExchange)
    {
      ergo_real BoysFunctionList_std[Nmax+1];
      ergo_real BoysFunctionList_mod[Nmax+1];
      ergo_real mu = paramsCAM.mu;
      ergo_real v1_squared = mu * mu / (mu * mu + alpha0);
      ergo_real v1 = template_blas_sqrt(v1_squared);
      /* Prepare BoysFunctionList_std */
      /* Use downward recursion to get Boys function values */
      ergo_real arg1 =  alpha0 * R_12_squared;
      ergo_real expMinusArg1 = template_blas_exp(-arg1);
      BoysFunctionList_std[Nmax] = integralInfo.BoysFunction(Nmax, arg1);
      for(int i = Nmax-1; i >= 0; i--)
	BoysFunctionList_std[i] = (2*arg1*BoysFunctionList_std[i+1] + expMinusArg1) / (2*i+1);
      /* Prepare BoysFunctionList_mod */
      /* Use downward recursion to get Boys function values */
      ergo_real arg2 =  alpha0 * R_12_squared * v1_squared;
      ergo_real expMinusArg2 = template_blas_exp(-arg2);
      BoysFunctionList_mod[Nmax] = integralInfo.BoysFunction(Nmax, arg2);
      for(int i = Nmax-1; i >= 0; i--)
	BoysFunctionList_mod[i] = (2*arg2*BoysFunctionList_mod[i+1] + expMinusArg2) / (2*i+1);
      // rescale
      for(int i = 0; i <= Nmax; i++)
	BoysFunctionList_mod[i] *= v1 * template_blas_pow(v1_squared, (ergo_real)i);  /* TODO: avoid using pow() here! */
      // add BoysFunctionList_std and BoysFunctionList_mod using weights given by cam_param_alpha and cam_param_beta
      for(int i = 0; i <= Nmax; i++)
	BoysList[i] = paramsCAM.alpha * BoysFunctionList_std[i] + paramsCAM.beta * BoysFunctionList_mod[i];
    }
  else
    {
      /* Compute all Boys function values needed */
      /* Use downward recursion to get Boys function values */
      ergo_real arg = alpha0 * R_12_squared;
      BoysList[Nmax] = integralInfo.BoysFunction(Nmax, arg);
      if(Nmax > 0)
	{
	  ergo_real expMinusArg = template_blas_exp(-arg); 
	  for(int i = Nmax-1; i >= 0; i--)
	    BoysList[i] = (2*arg*BoysList[i+1] + expMinusArg) / (2*i+1);
	}
    }



#if 1
  if(n1max == 0)
    {
      if(n2max == 0)
	{
	  primitiveIntegralList[0] = resultPreFactor * BoysList[0];
	  return 0;
	} // END IF (n2max == 0)
      if(n2max == 1)
	{
	  ergo_real R_val_0000 = BoysList[0] * resultPreFactor;
	  ergo_real R_val_1000 = -2*alpha0*BoysList[1] * resultPreFactor;
	  
	  ergo_real R_val_0100 = dx0 * R_val_1000;
	  ergo_real R_val_0010 = dx1 * R_val_1000;
	  ergo_real R_val_0001 = dx2 * R_val_1000;

	  primitiveIntegralList[0] = R_val_0000;

	  primitiveIntegralList[1] = R_val_0001;
	  primitiveIntegralList[2] = R_val_0010;
	  primitiveIntegralList[3] = R_val_0100;

	  return 0;
	} // END IF (n2max == 1)
      if(n2max == 2)
	{
	  ergo_real R_val_0000 = BoysList[0] * resultPreFactor;
	  ergo_real R_val_1000 = -2*alpha0*BoysList[1] * resultPreFactor;
	  ergo_real R_val_2000 =  4*alpha0*alpha0*BoysList[2] * resultPreFactor;

	  ergo_real R_val_1100 = dx0 * R_val_2000;
	  ergo_real R_val_1010 = dx1 * R_val_2000;
	  ergo_real R_val_1001 = dx2 * R_val_2000;

	  ergo_real R_val_0100 = dx0 * R_val_1000;
	  ergo_real R_val_0010 = dx1 * R_val_1000;
	  ergo_real R_val_0001 = dx2 * R_val_1000;

	  ergo_real R_val_0200 = R_val_1000 + dx0 * R_val_1100;
	  ergo_real R_val_0020 = R_val_1000 + dx1 * R_val_1010;
	  ergo_real R_val_0002 = R_val_1000 + dx2 * R_val_1001;

	  ergo_real R_val_0110 = dx0 * R_val_1010;
	  ergo_real R_val_0101 = dx0 * R_val_1001;
	  ergo_real R_val_0011 = dx1 * R_val_1001;

	  primitiveIntegralList[0] = R_val_0000;

	  primitiveIntegralList[1] = R_val_0001;
	  primitiveIntegralList[2] = R_val_0010;
	  primitiveIntegralList[3] = R_val_0100;

	  primitiveIntegralList[4] = R_val_0002;
	  primitiveIntegralList[5] = R_val_0011;
	  primitiveIntegralList[6] = R_val_0020;
	  primitiveIntegralList[7] = R_val_0101;
	  primitiveIntegralList[8] = R_val_0110;
	  primitiveIntegralList[9] = R_val_0200;

	  return 0;
	} // END IF (n2max == 1)
    }
  if(n1max == 1)
    {
      if(n2max == 0)
	{
	  ergo_real R_val_0000 = BoysList[0] * resultPreFactor;
	  ergo_real R_val_1000 = -2*alpha0*BoysList[1] * resultPreFactor;

	  ergo_real R_val_0100 = dx0 * R_val_1000;
	  ergo_real R_val_0010 = dx1 * R_val_1000;
	  ergo_real R_val_0001 = dx2 * R_val_1000;

	  primitiveIntegralList[0] = R_val_0000;

	  primitiveIntegralList[1] = -R_val_0001;
	  primitiveIntegralList[2] = -R_val_0010;
	  primitiveIntegralList[3] = -R_val_0100;

	  return 0;
	} // END IF (n2max == 0)
      if(n2max == 1)
	{
	  ergo_real R_val_0000 = BoysList[0] * resultPreFactor;
	  ergo_real R_val_1000 = -2*alpha0*BoysList[1] * resultPreFactor;
	  ergo_real R_val_2000 =  4*alpha0*alpha0*BoysList[2] * resultPreFactor;

	  ergo_real R_val_1100 = dx0 * R_val_2000;
	  ergo_real R_val_1010 = dx1 * R_val_2000;
	  ergo_real R_val_1001 = dx2 * R_val_2000;

	  ergo_real R_val_0100 = dx0 * R_val_1000;
	  ergo_real R_val_0010 = dx1 * R_val_1000;
	  ergo_real R_val_0001 = dx2 * R_val_1000;

	  ergo_real R_val_0200 = R_val_1000 + dx0 * R_val_1100;
	  ergo_real R_val_0020 = R_val_1000 + dx1 * R_val_1010;
	  ergo_real R_val_0002 = R_val_1000 + dx2 * R_val_1001;

	  ergo_real R_val_0110 = dx0 * R_val_1010;
	  ergo_real R_val_0101 = dx0 * R_val_1001;
	  ergo_real R_val_0011 = dx1 * R_val_1001;

	  // i1 = 0 (000)
	  primitiveIntegralList[ 0] =  R_val_0000; // i2 000
	  primitiveIntegralList[ 1] =  R_val_0001; // i2 001
	  primitiveIntegralList[ 2] =  R_val_0010; // i2 010
	  primitiveIntegralList[ 3] =  R_val_0100; // i2 100
	  // i1 = 1 (001)
	  primitiveIntegralList[ 4] = -R_val_0001; // i2 000
	  primitiveIntegralList[ 5] = -R_val_0002; // i2 001
	  primitiveIntegralList[ 6] = -R_val_0011; // i2 010
	  primitiveIntegralList[ 7] = -R_val_0101; // i2 100
	  // i1 = 2 (010)
	  primitiveIntegralList[ 8] = -R_val_0010; // i2 000
	  primitiveIntegralList[ 9] = -R_val_0011; // i2 001
	  primitiveIntegralList[10] = -R_val_0020; // i2 010
	  primitiveIntegralList[11] = -R_val_0110; // i2 100
	  // i1 = 3 (100)
	  primitiveIntegralList[12] = -R_val_0100; // i2 000
	  primitiveIntegralList[13] = -R_val_0101; // i2 001
	  primitiveIntegralList[14] = -R_val_0110; // i2 010
	  primitiveIntegralList[15] = -R_val_0200; // i2 100

	  return 0;
	} // END IF (n2max == 1)
      if(n2max == 2)
	{
	  ergo_real R_val_0000 = BoysList[0] * resultPreFactor;
	  ergo_real R_val_1000 = -2*alpha0*BoysList[1] * resultPreFactor; 
	  ergo_real R_val_2000 =  4*alpha0*alpha0*BoysList[2] * resultPreFactor;
	  ergo_real R_val_3000 = -8*alpha0*alpha0*alpha0*BoysList[3] * resultPreFactor;

	  ergo_real R_val_2100 = dx0 * R_val_3000;
	  ergo_real R_val_2010 = dx1 * R_val_3000;
	  ergo_real R_val_2001 = dx2 * R_val_3000;

	  ergo_real R_val_1100 = dx0 * R_val_2000;
	  ergo_real R_val_1010 = dx1 * R_val_2000;
	  ergo_real R_val_1001 = dx2 * R_val_2000;

	  ergo_real R_val_1200 = R_val_2000 + dx0 * R_val_2100;
	  ergo_real R_val_1020 = R_val_2000 + dx1 * R_val_2010;
	  ergo_real R_val_1002 = R_val_2000 + dx2 * R_val_2001;

	  ergo_real R_val_1110 = dx0 * R_val_2010;
	  ergo_real R_val_1101 = dx0 * R_val_2001;
	  ergo_real R_val_1011 = dx1 * R_val_2001;

	  ergo_real R_val_0100 = dx0 * R_val_1000;
	  ergo_real R_val_0010 = dx1 * R_val_1000;
	  ergo_real R_val_0001 = dx2 * R_val_1000;

	  ergo_real R_val_0200 = R_val_1000 + dx0 * R_val_1100;
	  ergo_real R_val_0020 = R_val_1000 + dx1 * R_val_1010;
	  ergo_real R_val_0002 = R_val_1000 + dx2 * R_val_1001;

	  ergo_real R_val_0110 = dx0 * R_val_1010;
	  ergo_real R_val_0101 = dx0 * R_val_1001;
	  ergo_real R_val_0011 = dx1 * R_val_1001;

	  ergo_real R_val_0111 =             dx0 * R_val_1011;
	  ergo_real R_val_0300 = 2 * R_val_1100 + dx0 * R_val_1200;
	  ergo_real R_val_0030 = 2 * R_val_1010 + dx1 * R_val_1020;
	  ergo_real R_val_0003 = 2 * R_val_1001 + dx2 * R_val_1002;
	  ergo_real R_val_0210 =     R_val_1010 + dx0 * R_val_1110;
	  ergo_real R_val_0201 =     R_val_1001 + dx0 * R_val_1101;
	  ergo_real R_val_0120 =     R_val_1100 + dx1 * R_val_1110;
	  ergo_real R_val_0021 =     R_val_1001 + dx1 * R_val_1011;
	  ergo_real R_val_0102 =     R_val_1100 + dx2 * R_val_1101;
	  ergo_real R_val_0012 =     R_val_1010 + dx2 * R_val_1011;

	  // i1 = 0 (000)
	  primitiveIntegralList[ 0] =  R_val_0000; // i2 000
	  primitiveIntegralList[ 1] =  R_val_0001; // i2 001
	  primitiveIntegralList[ 2] =  R_val_0010; // i2 010
	  primitiveIntegralList[ 3] =  R_val_0100; // i2 100
	  primitiveIntegralList[ 4] =  R_val_0002; // i2 002
	  primitiveIntegralList[ 5] =  R_val_0011; // i2 011
	  primitiveIntegralList[ 6] =  R_val_0020; // i2 020
	  primitiveIntegralList[ 7] =  R_val_0101; // i2 101
	  primitiveIntegralList[ 8] =  R_val_0110; // i2 110
	  primitiveIntegralList[ 9] =  R_val_0200; // i2 200
	  // i1 = 1 (001)
	  primitiveIntegralList[10] = -R_val_0001; // i2 000
	  primitiveIntegralList[11] = -R_val_0002; // i2 001
	  primitiveIntegralList[12] = -R_val_0011; // i2 010
	  primitiveIntegralList[13] = -R_val_0101; // i2 100
	  primitiveIntegralList[14] = -R_val_0003; // i2 002
	  primitiveIntegralList[15] = -R_val_0012; // i2 011
	  primitiveIntegralList[16] = -R_val_0021; // i2 020
	  primitiveIntegralList[17] = -R_val_0102; // i2 101
	  primitiveIntegralList[18] = -R_val_0111; // i2 110
	  primitiveIntegralList[19] = -R_val_0201; // i2 200
	  // i1 = 2 (010)
	  primitiveIntegralList[20] = -R_val_0010; // i2 000
	  primitiveIntegralList[21] = -R_val_0011; // i2 001
	  primitiveIntegralList[22] = -R_val_0020; // i2 010
	  primitiveIntegralList[23] = -R_val_0110; // i2 100
	  primitiveIntegralList[24] = -R_val_0012; // i2 002
	  primitiveIntegralList[25] = -R_val_0021; // i2 011
	  primitiveIntegralList[26] = -R_val_0030; // i2 020
	  primitiveIntegralList[27] = -R_val_0111; // i2 101
	  primitiveIntegralList[28] = -R_val_0120; // i2 110
	  primitiveIntegralList[29] = -R_val_0210; // i2 200
	  // i1 = 3 (100)
	  primitiveIntegralList[30] = -R_val_0100; // i2 000
	  primitiveIntegralList[31] = -R_val_0101; // i2 001
	  primitiveIntegralList[32] = -R_val_0110; // i2 010
	  primitiveIntegralList[33] = -R_val_0200; // i2 100
	  primitiveIntegralList[34] = -R_val_0102; // i2 002
	  primitiveIntegralList[35] = -R_val_0111; // i2 011
	  primitiveIntegralList[36] = -R_val_0120; // i2 020
	  primitiveIntegralList[37] = -R_val_0201; // i2 101
	  primitiveIntegralList[38] = -R_val_0210; // i2 110
	  primitiveIntegralList[39] = -R_val_0300; // i2 200

	  return 0;
	} // END IF (n2max == 2)
    }
  if(n1max == 2)
    {
      if(n2max == 0)
	{
	  ergo_real R_val_0000 = BoysList[0] * resultPreFactor;
	  ergo_real R_val_1000 = -2*alpha0*BoysList[1] * resultPreFactor;
	  ergo_real R_val_2000 =  4*alpha0*alpha0*BoysList[2] * resultPreFactor;

	  ergo_real R_val_1100 = dx0 * R_val_2000;
	  ergo_real R_val_1010 = dx1 * R_val_2000;
	  ergo_real R_val_1001 = dx2 * R_val_2000;

	  ergo_real R_val_0100 = dx0 * R_val_1000;
	  ergo_real R_val_0010 = dx1 * R_val_1000;
	  ergo_real R_val_0001 = dx2 * R_val_1000;

	  ergo_real R_val_0200 = R_val_1000 + dx0 * R_val_1100;
	  ergo_real R_val_0020 = R_val_1000 + dx1 * R_val_1010;
	  ergo_real R_val_0002 = R_val_1000 + dx2 * R_val_1001;

	  ergo_real R_val_0110 = dx0 * R_val_1010;
	  ergo_real R_val_0101 = dx0 * R_val_1001;
	  ergo_real R_val_0011 = dx1 * R_val_1001;

	  // i1 = 0 (000)
	  primitiveIntegralList[0] =  R_val_0000;
	  // i1 = 1 (001)
	  primitiveIntegralList[1] = -R_val_0001;
	  // i1 = 2 (010)
	  primitiveIntegralList[2] = -R_val_0010;
	  // i1 = 3 (100)
	  primitiveIntegralList[3] = -R_val_0100;
	  // i1 = 4 (002)
	  primitiveIntegralList[4] =  R_val_0002;
	  // i1 = 5 (011)
	  primitiveIntegralList[5] =  R_val_0011;
	  // i1 = 6 (020)
	  primitiveIntegralList[6] =  R_val_0020;
	  // i1 = 7 (101)
	  primitiveIntegralList[7] =  R_val_0101;
	  // i1 = 8 (110)
	  primitiveIntegralList[8] =  R_val_0110;
	  // i1 = 9 (200)
	  primitiveIntegralList[9] =  R_val_0200;

	  return 0;
	} // END IF (n2max == 0)
      if(n2max == 1)
	{
	  ergo_real R_val_0000 = BoysList[0] * resultPreFactor;
	  ergo_real R_val_1000 = -2*alpha0*BoysList[1] * resultPreFactor;
	  ergo_real R_val_2000 =  4*alpha0*alpha0*BoysList[2] * resultPreFactor;
	  ergo_real R_val_3000 = -8*alpha0*alpha0*alpha0*BoysList[3] * resultPreFactor;

	  ergo_real R_val_2100 = dx0 * R_val_3000;
	  ergo_real R_val_2010 = dx1 * R_val_3000;
	  ergo_real R_val_2001 = dx2 * R_val_3000;

	  ergo_real R_val_1100 = dx0 * R_val_2000;
	  ergo_real R_val_1010 = dx1 * R_val_2000;
	  ergo_real R_val_1001 = dx2 * R_val_2000;

	  ergo_real R_val_1200 = R_val_2000 + dx0 * R_val_2100;
	  ergo_real R_val_1020 = R_val_2000 + dx1 * R_val_2010;
	  ergo_real R_val_1002 = R_val_2000 + dx2 * R_val_2001;

	  ergo_real R_val_1110 = dx0 * R_val_2010;
	  ergo_real R_val_1101 = dx0 * R_val_2001;
	  ergo_real R_val_1011 = dx1 * R_val_2001;

	  ergo_real R_val_0100 = dx0 * R_val_1000;
	  ergo_real R_val_0010 = dx1 * R_val_1000;
	  ergo_real R_val_0001 = dx2 * R_val_1000;

	  ergo_real R_val_0200 = R_val_1000 + dx0 * R_val_1100;
	  ergo_real R_val_0020 = R_val_1000 + dx1 * R_val_1010;
	  ergo_real R_val_0002 = R_val_1000 + dx2 * R_val_1001;

	  ergo_real R_val_0110 = dx0 * R_val_1010;
	  ergo_real R_val_0101 = dx0 * R_val_1001;
	  ergo_real R_val_0011 = dx1 * R_val_1001;

	  ergo_real R_val_0111 =             dx0 * R_val_1011;
	  ergo_real R_val_0300 = 2 * R_val_1100 + dx0 * R_val_1200;
	  ergo_real R_val_0030 = 2 * R_val_1010 + dx1 * R_val_1020;
	  ergo_real R_val_0003 = 2 * R_val_1001 + dx2 * R_val_1002;
	  ergo_real R_val_0210 =     R_val_1010 + dx0 * R_val_1110;
	  ergo_real R_val_0201 =     R_val_1001 + dx0 * R_val_1101;
	  ergo_real R_val_0120 =     R_val_1100 + dx1 * R_val_1110;
	  ergo_real R_val_0021 =     R_val_1001 + dx1 * R_val_1011;
	  ergo_real R_val_0102 =     R_val_1100 + dx2 * R_val_1101;
	  ergo_real R_val_0012 =     R_val_1010 + dx2 * R_val_1011;

	  // i1 = 0 (000)
	  primitiveIntegralList[ 0] =  R_val_0000; // i2 000
	  primitiveIntegralList[ 1] =  R_val_0001; // i2 001
	  primitiveIntegralList[ 2] =  R_val_0010; // i2 010
	  primitiveIntegralList[ 3] =  R_val_0100; // i2 100
	  // i1 = 1 (001)
	  primitiveIntegralList[ 4] = -R_val_0001; // i2 000
	  primitiveIntegralList[ 5] = -R_val_0002; // i2 001
	  primitiveIntegralList[ 6] = -R_val_0011; // i2 010
	  primitiveIntegralList[ 7] = -R_val_0101; // i2 100
	  // i1 = 2 (010)
	  primitiveIntegralList[ 8] = -R_val_0010; // i2 000
	  primitiveIntegralList[ 9] = -R_val_0011; // i2 001
	  primitiveIntegralList[10] = -R_val_0020; // i2 010
	  primitiveIntegralList[11] = -R_val_0110; // i2 100
	  // i1 = 3 (100)
	  primitiveIntegralList[12] = -R_val_0100; // i2 000
	  primitiveIntegralList[13] = -R_val_0101; // i2 001
	  primitiveIntegralList[14] = -R_val_0110; // i2 010
	  primitiveIntegralList[15] = -R_val_0200; // i2 100
	  // i1 = 4 (002)
	  primitiveIntegralList[16] =  R_val_0002; // i2 000
	  primitiveIntegralList[17] =  R_val_0003; // i2 001
	  primitiveIntegralList[18] =  R_val_0012; // i2 010
	  primitiveIntegralList[19] =  R_val_0102; // i2 100
	  // i1 = 5 (011)
	  primitiveIntegralList[20] =  R_val_0011; // i2 000
	  primitiveIntegralList[21] =  R_val_0012; // i2 001
	  primitiveIntegralList[22] =  R_val_0021; // i2 010
	  primitiveIntegralList[23] =  R_val_0111; // i2 100
	  // i1 = 6 (020)
	  primitiveIntegralList[24] =  R_val_0020; // i2 000
	  primitiveIntegralList[25] =  R_val_0021; // i2 001
	  primitiveIntegralList[26] =  R_val_0030; // i2 010
	  primitiveIntegralList[27] =  R_val_0120; // i2 100
	  // i1 = 7 (101)
	  primitiveIntegralList[28] =  R_val_0101; // i2 000
	  primitiveIntegralList[29] =  R_val_0102; // i2 001
	  primitiveIntegralList[30] =  R_val_0111; // i2 010
	  primitiveIntegralList[31] =  R_val_0201; // i2 100
	  // i1 = 8 (110)
	  primitiveIntegralList[32] =  R_val_0110; // i2 000
	  primitiveIntegralList[33] =  R_val_0111; // i2 001
	  primitiveIntegralList[34] =  R_val_0120; // i2 010
	  primitiveIntegralList[35] =  R_val_0210; // i2 100
	  // i1 = 9 (200)
	  primitiveIntegralList[36] =  R_val_0200; // i2 000
	  primitiveIntegralList[37] =  R_val_0201; // i2 001
	  primitiveIntegralList[38] =  R_val_0210; // i2 010
	  primitiveIntegralList[39] =  R_val_0300; // i2 100

	  return 0;
	} // END IF (n2max == 0)
      if(n2max == 2)
	{
	  ergo_real R_val_0000 = BoysList[0] * resultPreFactor;
	  ergo_real R_val_1000 = -2*alpha0*BoysList[1] * resultPreFactor;
	  ergo_real R_val_2000 =  4*alpha0*alpha0*BoysList[2] * resultPreFactor;
	  ergo_real R_val_3000 = -8*alpha0*alpha0*alpha0*BoysList[3] * resultPreFactor;
	  ergo_real R_val_4000 = 16*alpha0*alpha0*alpha0*alpha0*BoysList[4] * resultPreFactor;

	  ergo_real R_val_3100 = dx0 * R_val_4000;
	  ergo_real R_val_3010 = dx1 * R_val_4000;
	  ergo_real R_val_3001 = dx2 * R_val_4000;

	  ergo_real R_val_2100 = dx0 * R_val_3000;
	  ergo_real R_val_2010 = dx1 * R_val_3000;
	  ergo_real R_val_2001 = dx2 * R_val_3000;

	  ergo_real R_val_2200 = R_val_3000 + dx0 * R_val_3100;
	  ergo_real R_val_2020 = R_val_3000 + dx1 * R_val_3010;
	  ergo_real R_val_2002 = R_val_3000 + dx2 * R_val_3001;

	  ergo_real R_val_2110 = dx0 * R_val_3010;
	  ergo_real R_val_2101 = dx0 * R_val_3001;
	  ergo_real R_val_2011 = dx1 * R_val_3001;

	  ergo_real R_val_1100 = dx0 * R_val_2000;
	  ergo_real R_val_1010 = dx1 * R_val_2000;
	  ergo_real R_val_1001 = dx2 * R_val_2000;

	  ergo_real R_val_1200 = R_val_2000 + dx0 * R_val_2100;
	  ergo_real R_val_1020 = R_val_2000 + dx1 * R_val_2010;
	  ergo_real R_val_1002 = R_val_2000 + dx2 * R_val_2001;

	  ergo_real R_val_1110 = dx0 * R_val_2010;
	  ergo_real R_val_1101 = dx0 * R_val_2001;
	  ergo_real R_val_1011 = dx1 * R_val_2001;

	  ergo_real R_val_1111 =             dx0 * R_val_2011;
	  ergo_real R_val_1300 = 2 * R_val_2100 + dx0 * R_val_2200;
	  ergo_real R_val_1030 = 2 * R_val_2010 + dx1 * R_val_2020;
	  ergo_real R_val_1003 = 2 * R_val_2001 + dx2 * R_val_2002;
	  ergo_real R_val_1210 =     R_val_2010 + dx0 * R_val_2110;
	  ergo_real R_val_1201 =     R_val_2001 + dx0 * R_val_2101;
	  ergo_real R_val_1120 =     R_val_2100 + dx1 * R_val_2110;
	  ergo_real R_val_1021 =     R_val_2001 + dx1 * R_val_2011;
	  ergo_real R_val_1102 =     R_val_2100 + dx2 * R_val_2101;
	  ergo_real R_val_1012 =     R_val_2010 + dx2 * R_val_2011;

	  ergo_real R_val_0100 = dx0 * R_val_1000;
	  ergo_real R_val_0010 = dx1 * R_val_1000;
	  ergo_real R_val_0001 = dx2 * R_val_1000;

	  ergo_real R_val_0200 = R_val_1000 + dx0 * R_val_1100;
	  ergo_real R_val_0020 = R_val_1000 + dx1 * R_val_1010;
	  ergo_real R_val_0002 = R_val_1000 + dx2 * R_val_1001;

	  ergo_real R_val_0110 = dx0 * R_val_1010;
	  ergo_real R_val_0101 = dx0 * R_val_1001;
	  ergo_real R_val_0011 = dx1 * R_val_1001;

	  ergo_real R_val_0111 =             dx0 * R_val_1011;
	  ergo_real R_val_0300 = 2 * R_val_1100 + dx0 * R_val_1200;
	  ergo_real R_val_0030 = 2 * R_val_1010 + dx1 * R_val_1020;
	  ergo_real R_val_0003 = 2 * R_val_1001 + dx2 * R_val_1002;
	  ergo_real R_val_0210 =     R_val_1010 + dx0 * R_val_1110;
	  ergo_real R_val_0201 =     R_val_1001 + dx0 * R_val_1101;
	  ergo_real R_val_0120 =     R_val_1100 + dx1 * R_val_1110;
	  ergo_real R_val_0021 =     R_val_1001 + dx1 * R_val_1011;
	  ergo_real R_val_0102 =     R_val_1100 + dx2 * R_val_1101;
	  ergo_real R_val_0012 =     R_val_1010 + dx2 * R_val_1011;

	  ergo_real R_val_0400 = 3 * R_val_1200 + dx0 * R_val_1300;
	  ergo_real R_val_0040 = 3 * R_val_1020 + dx1 * R_val_1030;
	  ergo_real R_val_0004 = 3 * R_val_1002 + dx2 * R_val_1003;
	  ergo_real R_val_0310 = 2 * R_val_1110 + dx0 * R_val_1210;
	  ergo_real R_val_0301 = 2 * R_val_1101 + dx0 * R_val_1201;
	  ergo_real R_val_0130 = 2 * R_val_1110 + dx1 * R_val_1120;
	  ergo_real R_val_0031 = 2 * R_val_1011 + dx1 * R_val_1021;
	  ergo_real R_val_0103 = 2 * R_val_1101 + dx2 * R_val_1102;
	  ergo_real R_val_0013 = 2 * R_val_1011 + dx2 * R_val_1012;
	  ergo_real R_val_0220 =     R_val_1020 + dx0 * R_val_1120;
	  ergo_real R_val_0202 =     R_val_1002 + dx0 * R_val_1102;
	  ergo_real R_val_0022 =     R_val_1002 + dx1 * R_val_1012;
	  ergo_real R_val_0211 =     R_val_1011 + dx0 * R_val_1111;
	  ergo_real R_val_0121 =     R_val_1101 + dx1 * R_val_1111;
	  ergo_real R_val_0112 =     R_val_1110 + dx2 * R_val_1111;

	  // i1 = 0 (000)
	  primitiveIntegralList[ 0] =  R_val_0000; // i2 000
	  primitiveIntegralList[ 1] =  R_val_0001; // i2 001
	  primitiveIntegralList[ 2] =  R_val_0010; // i2 010
	  primitiveIntegralList[ 3] =  R_val_0100; // i2 100
	  primitiveIntegralList[ 4] =  R_val_0002; // i2 002
	  primitiveIntegralList[ 5] =  R_val_0011; // i2 011
	  primitiveIntegralList[ 6] =  R_val_0020; // i2 020
	  primitiveIntegralList[ 7] =  R_val_0101; // i2 101
	  primitiveIntegralList[ 8] =  R_val_0110; // i2 110
	  primitiveIntegralList[ 9] =  R_val_0200; // i2 200
	  // i1 = 1 (001)
	  primitiveIntegralList[10] = -R_val_0001; // i2 000
	  primitiveIntegralList[11] = -R_val_0002; // i2 001
	  primitiveIntegralList[12] = -R_val_0011; // i2 010
	  primitiveIntegralList[13] = -R_val_0101; // i2 100
	  primitiveIntegralList[14] = -R_val_0003; // i2 002
	  primitiveIntegralList[15] = -R_val_0012; // i2 011
	  primitiveIntegralList[16] = -R_val_0021; // i2 020
	  primitiveIntegralList[17] = -R_val_0102; // i2 101
	  primitiveIntegralList[18] = -R_val_0111; // i2 110
	  primitiveIntegralList[19] = -R_val_0201; // i2 200
	  // i1 = 2 (010)
	  primitiveIntegralList[20] = -R_val_0010; // i2 000
	  primitiveIntegralList[21] = -R_val_0011; // i2 001
	  primitiveIntegralList[22] = -R_val_0020; // i2 010
	  primitiveIntegralList[23] = -R_val_0110; // i2 100
	  primitiveIntegralList[24] = -R_val_0012; // i2 002
	  primitiveIntegralList[25] = -R_val_0021; // i2 011
	  primitiveIntegralList[26] = -R_val_0030; // i2 020
	  primitiveIntegralList[27] = -R_val_0111; // i2 101
	  primitiveIntegralList[28] = -R_val_0120; // i2 110
	  primitiveIntegralList[29] = -R_val_0210; // i2 200
	  // i1 = 3 (100)
	  primitiveIntegralList[30] = -R_val_0100; // i2 000
	  primitiveIntegralList[31] = -R_val_0101; // i2 001
	  primitiveIntegralList[32] = -R_val_0110; // i2 010
	  primitiveIntegralList[33] = -R_val_0200; // i2 100
	  primitiveIntegralList[34] = -R_val_0102; // i2 002
	  primitiveIntegralList[35] = -R_val_0111; // i2 011
	  primitiveIntegralList[36] = -R_val_0120; // i2 020
	  primitiveIntegralList[37] = -R_val_0201; // i2 101
	  primitiveIntegralList[38] = -R_val_0210; // i2 110
	  primitiveIntegralList[39] = -R_val_0300; // i2 200
	  // i1 = 4 (002)
	  primitiveIntegralList[40] =  R_val_0002; // i2 000
	  primitiveIntegralList[41] =  R_val_0003; // i2 001
	  primitiveIntegralList[42] =  R_val_0012; // i2 010
	  primitiveIntegralList[43] =  R_val_0102; // i2 100
	  primitiveIntegralList[44] =  R_val_0004; // i2 002
	  primitiveIntegralList[45] =  R_val_0013; // i2 011
	  primitiveIntegralList[46] =  R_val_0022; // i2 020
	  primitiveIntegralList[47] =  R_val_0103; // i2 101
	  primitiveIntegralList[48] =  R_val_0112; // i2 110
	  primitiveIntegralList[49] =  R_val_0202; // i2 200
	  // i1 = 5 (011)
	  primitiveIntegralList[50] =  R_val_0011; // i2 000
	  primitiveIntegralList[51] =  R_val_0012; // i2 001
	  primitiveIntegralList[52] =  R_val_0021; // i2 010
	  primitiveIntegralList[53] =  R_val_0111; // i2 100
	  primitiveIntegralList[54] =  R_val_0013; // i2 002
	  primitiveIntegralList[55] =  R_val_0022; // i2 011
	  primitiveIntegralList[56] =  R_val_0031; // i2 020
	  primitiveIntegralList[57] =  R_val_0112; // i2 101
	  primitiveIntegralList[58] =  R_val_0121; // i2 110
	  primitiveIntegralList[59] =  R_val_0211; // i2 200
	  // i1 = 6 (020)
	  primitiveIntegralList[60] =  R_val_0020; // i2 000
	  primitiveIntegralList[61] =  R_val_0021; // i2 001
	  primitiveIntegralList[62] =  R_val_0030; // i2 010
	  primitiveIntegralList[63] =  R_val_0120; // i2 100
	  primitiveIntegralList[64] =  R_val_0022; // i2 002
	  primitiveIntegralList[65] =  R_val_0031; // i2 011
	  primitiveIntegralList[66] =  R_val_0040; // i2 020
	  primitiveIntegralList[67] =  R_val_0121; // i2 101
	  primitiveIntegralList[68] =  R_val_0130; // i2 110
	  primitiveIntegralList[69] =  R_val_0220; // i2 200
	  // i1 = 7 (101)
	  primitiveIntegralList[70] =  R_val_0101; // i2 000
	  primitiveIntegralList[71] =  R_val_0102; // i2 001
	  primitiveIntegralList[72] =  R_val_0111; // i2 010
	  primitiveIntegralList[73] =  R_val_0201; // i2 100
	  primitiveIntegralList[74] =  R_val_0103; // i2 002
	  primitiveIntegralList[75] =  R_val_0112; // i2 011
	  primitiveIntegralList[76] =  R_val_0121; // i2 020
	  primitiveIntegralList[77] =  R_val_0202; // i2 101
	  primitiveIntegralList[78] =  R_val_0211; // i2 110
	  primitiveIntegralList[79] =  R_val_0301; // i2 200
	  // i1 = 8 (110)
	  primitiveIntegralList[80] =  R_val_0110; // i2 000
	  primitiveIntegralList[81] =  R_val_0111; // i2 001
	  primitiveIntegralList[82] =  R_val_0120; // i2 010
	  primitiveIntegralList[83] =  R_val_0210; // i2 100
	  primitiveIntegralList[84] =  R_val_0112; // i2 002
	  primitiveIntegralList[85] =  R_val_0121; // i2 011
	  primitiveIntegralList[86] =  R_val_0130; // i2 020
	  primitiveIntegralList[87] =  R_val_0211; // i2 101
	  primitiveIntegralList[88] =  R_val_0220; // i2 110
	  primitiveIntegralList[89] =  R_val_0310; // i2 200
	  // i1 = 9 (200)
	  primitiveIntegralList[90] =  R_val_0200; // i2 000
	  primitiveIntegralList[91] =  R_val_0201; // i2 001
	  primitiveIntegralList[92] =  R_val_0210; // i2 010
	  primitiveIntegralList[93] =  R_val_0300; // i2 100
	  primitiveIntegralList[94] =  R_val_0202; // i2 002
	  primitiveIntegralList[95] =  R_val_0211; // i2 011
	  primitiveIntegralList[96] =  R_val_0220; // i2 020
	  primitiveIntegralList[97] =  R_val_0301; // i2 101
	  primitiveIntegralList[98] =  R_val_0310; // i2 110
	  primitiveIntegralList[99] =  R_val_0400; // i2 200

	  return 0;
	}
    }
#endif

  ergo_real R[Nmax+1][Nmax+1][Nmax+1][Nmax+1];

  int n;
  ergo_real factor = 1;
  for(n = 0; n <= Nmax; n++)
    {
      R[n][0][0][0] = factor * BoysList[n];
      factor *= -2*alpha0;
    }

  int minus1topowList[Nmax+1];
  int ifactor = 1;
  for(n = 0; n <= Nmax; n++)
    {
      minus1topowList[n] = ifactor;
      ifactor *= -1;
    }
  
  // Use recurrences to get remaining R values
  for(n = Nmax - 1; n >= 0; n--)
    {
      int nn = Nmax - n;
      int noOfMonomials = integralInfo.monomial_info.no_of_monomials_list[nn];
      int i;
      for(i = 1; i < noOfMonomials; i++)
	{
	  int ix = integralInfo.monomial_info.monomial_list[i].ix;
	  int iy = integralInfo.monomial_info.monomial_list[i].iy;
	  int iz = integralInfo.monomial_info.monomial_list[i].iz;
	  if(ix > 0)
	    {
	      ergo_real Rval = dx0 * R[n+1][ix-1][iy][iz];
	      if(ix > 1)
		Rval += (ix - 1) * R[n+1][ix-2][iy][iz];
	      R[n][ix][iy][iz] = Rval;
	    }
	  else if(iy > 0)
	    {
	      ergo_real Rval = dx1 * R[n+1][ix][iy-1][iz];
	      if(iy > 1)
		Rval += (iy - 1) * R[n+1][ix][iy-2][iz];
	      R[n][ix][iy][iz] = Rval;
	    }
	  else if(iz > 0)
	    {
	      ergo_real Rval = dx2 * R[n+1][ix][iy][iz-1];
	      if(iz > 1)
		Rval += (iz - 1) * R[n+1][ix][iy][iz-2];
	      R[n][ix][iy][iz] = Rval;
	    }
	} // END FOR i
    } // END FOR n

  int i1, i2;
  for(i1 = 0; i1 < noOfMonomials_1; i1++)
    {
      int ix1 = integralInfo.monomial_info.monomial_list[i1].ix;
      int iy1 = integralInfo.monomial_info.monomial_list[i1].iy;
      int iz1 = integralInfo.monomial_info.monomial_list[i1].iz;
      int n1 = ix1+iy1+iz1;
      ergo_real prefactor = minus1topowList[n1] * resultPreFactor;
      for(i2 = 0; i2 < noOfMonomials_2; i2++)
	{
	  int ix2 = integralInfo.monomial_info.monomial_list[i2].ix;
	  int iy2 = integralInfo.monomial_info.monomial_list[i2].iy;
	  int iz2 = integralInfo.monomial_info.monomial_list[i2].iz;
	  primitiveIntegralList[i1*noOfMonomials_2+i2] = prefactor * R[0][ix1+ix2][iy1+iy2][iz1+iz2];
	} // END FOR i2  
    }
  
  return 0;
}


