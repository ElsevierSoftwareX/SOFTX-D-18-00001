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

/** @file cubature_rules.cc

    @brief Code for various cubature rules that can be used for
    numerical integration.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <stdexcept>

#include "cubature_rules.h"


int
use_cubature_rule(int maxlen,
		  real (*coor)[3],
		  real *weight,
		  BoxStruct* box,
		  int ruleNumber)
{
  real volume, diff0, diff1, diff2;
  real c0, c1, c2, a, b;
  real currCoords[3];
  int Ngrid, currIndex;
  int i, j, k, ii;
  real a0, a1, a2;

  volume = 1;
  for(i = 0; i < NO_OF_DIMENSIONS; i++)
    volume *= (box->max[i] - box->min[i]);
  
  switch(ruleNumber)
    {
    case 1: /* single point in center of box */
      Ngrid = 1;
      if(Ngrid >= maxlen)
	throw std::runtime_error("error in use_cubature_rule: (Ngrid >= maxlen)");
      for(i = 0; i < NO_OF_DIMENSIONS; i++)
	  coor[0][i] = (box->max[i] + box->min[i]) / 2;
      weight[0] = volume;
      break;
    case 2: /* eight points towards corners of box */
      Ngrid = 8;
      if(Ngrid >= maxlen)
	throw std::runtime_error("error in use_cubature_rule: (Ngrid >= maxlen)");
      for(i = 0; i < Ngrid; i++)
	weight[i] = volume / 8;
      diff0 = box->max[0] - box->min[0];
      diff1 = box->max[1] - box->min[1];
      diff2 = box->max[2] - box->min[2];
      currIndex = 0;
      for(i = 0; i < 2; i++)
	{
	  currCoords[0] = box->min[0] + 0.25*diff0 + 0.5*diff0*i;
	  for(j = 0; j < 2; j++)
	    {
	      currCoords[1] = box->min[1] + 0.25*diff1 + 0.5*diff1*i;
	      for(k = 0; k < 2; k++)
		{
		  currCoords[2] = box->min[2] + 0.25*diff2 + 0.5*diff2*i;
		  for(ii = 0; ii < 3; ii++)
		    {
		      coor[currIndex][ii] = currCoords[ii];
		    } /* END FOR ii */
		  currIndex++;
		} /* END FOR k */
	    } /* END FOR j */
	} /* END FOR i */
      break;
    case 3: /* 14 point, degree 5 rule (Stroud 1971) */
      Ngrid = 14;
      if(Ngrid >= maxlen)
	throw std::runtime_error("error in use_cubature_rule: (Ngrid >= maxlen)");
      for(i = 0; i < 6; i++)
	weight[i] = volume *  0.88642659279778393 / 8;
      for(i = 6; i < 14; i++)
	weight[i] = volume * 0.33518005540166204 / 8;
      diff0 = box->max[0] - box->min[0];
      diff1 = box->max[1] - box->min[1];
      diff2 = box->max[2] - box->min[2];
      c0 = (box->max[0] + box->min[0]) / 2;
      c1 = (box->max[1] + box->min[1]) / 2;
      c2 = (box->max[2] + box->min[2]) / 2;
      a = 0.79582242575422146 * 0.5;
      b = 0.75878691063932814 * 0.5;

#define MACRO_3VECT(v,x,y,z) v[0]=x; v[1]=y; v[2]=z;

      MACRO_3VECT(coor[0], c0-a*diff0, c1,         c2        );
      MACRO_3VECT(coor[1], c0+a*diff0, c1,         c2        );
      MACRO_3VECT(coor[2], c0        , c1-a*diff1, c2        );
      MACRO_3VECT(coor[3], c0        , c1+a*diff1, c2        );
      MACRO_3VECT(coor[4], c0        , c1        , c2-a*diff2);
      MACRO_3VECT(coor[5], c0        , c1        , c2+a*diff2);

      MACRO_3VECT(coor[ 6], c0-b*diff0, c1-b*diff1, c2-b*diff2);
      MACRO_3VECT(coor[ 7], c0-b*diff0, c1-b*diff1, c2+b*diff2);
      MACRO_3VECT(coor[ 8], c0-b*diff0, c1+b*diff1, c2-b*diff2);
      MACRO_3VECT(coor[ 9], c0-b*diff0, c1+b*diff1, c2+b*diff2);
      MACRO_3VECT(coor[10], c0+b*diff0, c1-b*diff1, c2-b*diff2);
      MACRO_3VECT(coor[11], c0+b*diff0, c1-b*diff1, c2+b*diff2);
      MACRO_3VECT(coor[12], c0+b*diff0, c1+b*diff1, c2-b*diff2);
      MACRO_3VECT(coor[13], c0+b*diff0, c1+b*diff1, c2+b*diff2);

      break;

    case 4: /* 25 point, degree 5 rule (Stroud 1971) */
      Ngrid = 25;
      if(Ngrid >= maxlen)
	throw std::runtime_error("error in use_cubature_rule: (Ngrid >= maxlen)");
      weight[0] = volume * 1.6842105263157894 / 8;
      for(i = 1; i < 25; i++)
	weight[i] = volume *  0.26315789473684210 / 8;

      diff0 = box->max[0] - box->min[0];
      diff1 = box->max[1] - box->min[1];
      diff2 = box->max[2] - box->min[2];
      c0 = (box->max[0] + box->min[0]) / 2;
      c1 = (box->max[1] + box->min[1]) / 2;
      c2 = (box->max[2] + box->min[2]) / 2;
      a = 0.47800981191507060 * 0.5;
      b = 0.89982215247931316 * 0.5;

      MACRO_3VECT(coor[0], c0, c1, c2);

      ii = 1;
      a0 = a;
      a1 = a;
      a2 = b;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;
      a0 = a;
      a1 = b;
      a2 = a;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;
      a0 = b;
      a1 = a;
      a2 = a;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;

      break;

    case 5: /* 27 point, degree 7 rule (Stroud 1971) */
      Ngrid = 27;
      if(Ngrid >= maxlen)
	throw std::runtime_error("error in use_cubature_rule: (Ngrid >= maxlen)");
      weight[0] = volume * 0.78807348274421057 / 8;
      for(i = 1; i < 7; i++)
	weight[i] = volume *  0.49936900230772032 / 8;
      for(i = 7; i < 19; i++)
	weight[i] = volume *  0.032303742334037395 / 8;
      for(i = 19; i < 27; i++)
	weight[i] = volume *  0.47850844942512734 / 8;

      diff0 = box->max[0] - box->min[0];
      diff1 = box->max[1] - box->min[1];
      diff2 = box->max[2] - box->min[2];
      c0 = (box->max[0] + box->min[0]) / 2;
      c1 = (box->max[1] + box->min[1]) / 2;
      c2 = (box->max[2] + box->min[2]) / 2;

      MACRO_3VECT(coor[0], c0, c1, c2);
      a = 0.84841801147225245 * 0.5;
      ii = 1;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1        , c2        ); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1        , c2        ); ii++;
      MACRO_3VECT(coor[ii], c0        , c1-a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0        , c1+a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0        , c1        , c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1        , c2+a*diff2); ii++;
      a = 1.1064128986267175 * 0.5;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1-a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1+a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1-a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1+a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1        , c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1        , c2+a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1        , c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1        , c2+a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1-a*diff1, c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1-a*diff1, c2+a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1+a*diff1, c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1+a*diff1, c2+a*diff2); ii++;
      a0 = a1 = a2 = 0.65281647210169120 * 0.5;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;

      break;

    case 6: /* 32 point, degree 7 rule (Beckers 1992) */
      Ngrid = 32;
      if(Ngrid >= maxlen)
	throw std::runtime_error("error in use_cubature_rule: (Ngrid >= maxlen)");
      for(i = 0; i < 6; i++)
	weight[i] = volume * 0.14098806933910446 / 8;
      for(i = 6; i < 12; i++)
	weight[i] = volume * 0.53332245896607639 / 8;
      for(i = 12; i < 24; i++)
	weight[i] = volume *  0.049451452995044458 / 8;
      for(i = 24; i < 32; i++)
	weight[i] = volume *  0.42008992427854766 / 8;

      diff0 = box->max[0] - box->min[0];
      diff1 = box->max[1] - box->min[1];
      diff2 = box->max[2] - box->min[2];
      c0 = (box->max[0] + box->min[0]) / 2;
      c1 = (box->max[1] + box->min[1]) / 2;
      c2 = (box->max[2] + box->min[2]) / 2;

      ii = 0;

      a = 1.0 * 0.5;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1        , c2        ); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1        , c2        ); ii++;
      MACRO_3VECT(coor[ii], c0        , c1-a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0        , c1+a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0        , c1        , c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1        , c2+a*diff2); ii++;
      a = 0.66289786904352112 * 0.5;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1        , c2        ); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1        , c2        ); ii++;
      MACRO_3VECT(coor[ii], c0        , c1-a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0        , c1+a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0        , c1        , c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1        , c2+a*diff2); ii++;
      a = 1.0306143700994171 * 0.5;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1-a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1+a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1-a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1+a*diff1, c2        ); ii++;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1        , c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a*diff0, c1        , c2+a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1        , c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a*diff0, c1        , c2+a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1-a*diff1, c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1-a*diff1, c2+a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1+a*diff1, c2-a*diff2); ii++;
      MACRO_3VECT(coor[ii], c0        , c1+a*diff1, c2+a*diff2); ii++;
      a0 = a1 = a2 = 0.66713797405746656 * 0.5;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0-a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1-a1*diff1, c2+a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2-a2*diff2); ii++;
      MACRO_3VECT(coor[ii], c0+a0*diff0, c1+a1*diff1, c2+a2*diff2); ii++;

      break;

    default:
      throw std::runtime_error("error: unknown cubature rule.");
    } /* END SWITCH */

  /*  
  real testSum = 0;
  for(i = 0; i < Ngrid; i++)
    testSum += weight[i];
  printf("testSum = %22.11f\n", testSum);
  printf("volume  = %22.11f\n\n", volume);
  */

  return Ngrid;
} /* END use_cubature_rule */


