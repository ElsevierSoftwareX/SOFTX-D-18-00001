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
 


#include "template_lapack_common.h"


int dummyroutine()
{
  return 0;
}









/* assign strings:  a = b */

#ifdef KR_headers
VOID s_copy(a, b, la, lb) register char *a, *b; ftnlen la, lb;
#else
void s_copy(register char *a, const register char *b, ftnlen la, ftnlen lb)
#endif
{
const register char *aend, *bend;

aend = a + la;

if(la <= lb)
        while(a < aend)
                *a++ = *b++;

else
        {
        bend = b + lb;
        while(b < bend)
                *a++ = *b++;
        while(a < aend)
                *a++ = ' ';
        }
}



/* compare two strings */

integer s_cmp(const register char *a, const register char *b, ftnlen la, ftnlen lb)
{
const register char *aend, *bend;
aend = a + la;
bend = b + lb;

if(la <= lb)
        {
        while(a < aend)
                if(*a != *b)
                        return( *a - *b );
                else
                        { ++a; ++b; }

        while(b < bend)
                if(*b != ' ')
                        return( ' ' - *b );
                else    ++b;
        }

else
        {
        while(b < bend)
                if(*a == *b)
                        { ++a; ++b; }
                else
                        return( *a - *b );
        while(a < aend)
                if(*a != ' ')
                        return(*a - ' ');
                else    ++a;
        }
return(0);
}





typedef double real;






integer ieeeck_(integer *ispec, real *zero, real *one)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1998   


    Purpose   
    =======   

    IEEECK is called from the ILAENV to verify that Infinity and   
    possibly NaN arithmetic is safe (i.e. will not trap).   

    Arguments   
    =========   

    ISPEC   (input) INTEGER   
            Specifies whether to test just for inifinity arithmetic   
            or whether to test for infinity and NaN arithmetic.   
            = 0: Verify infinity arithmetic only.   
            = 1: Verify infinity and NaN arithmetic.   

    ZERO    (input) REAL   
            Must contain the value 0.0   
            This is passed to prevent the compiler from optimizing   
            away this code.   

    ONE     (input) REAL   
            Must contain the value 1.0   
            This is passed to prevent the compiler from optimizing   
            away this code.   

    RETURN VALUE:  INTEGER   
            = 0:  Arithmetic failed to produce the correct answers   
            = 1:  Arithmetic produced the correct answers */
    /* System generated locals */
    integer ret_val;
    /* Local variables */
     real neginf, posinf, negzro, newzro, nan1, nan2, nan3, nan4, nan5, 
	    nan6;


    ret_val = 1;

    posinf = *one / *zero;
    if (posinf <= *one) {
	ret_val = 0;
	return ret_val;
    }

    neginf = -(*one) / *zero;
    if (neginf >= *zero) {
	ret_val = 0;
	return ret_val;
    }

    negzro = *one / (neginf + *one);
    if (negzro != *zero) {
	ret_val = 0;
	return ret_val;
    }

    neginf = *one / negzro;
    if (neginf >= *zero) {
	ret_val = 0;
	return ret_val;
    }

    newzro = negzro + *zero;
    if (newzro != *zero) {
	ret_val = 0;
	return ret_val;
    }

    posinf = *one / newzro;
    if (posinf <= *one) {
	ret_val = 0;
	return ret_val;
    }

    neginf *= posinf;
    if (neginf >= *zero) {
	ret_val = 0;
	return ret_val;
    }

    posinf *= posinf;
    if (posinf <= *one) {
	ret_val = 0;
	return ret_val;
    }




/*     Return if we were only asked to check infinity arithmetic */

    if (*ispec == 0) {
	return ret_val;
    }

    nan1 = posinf + neginf;

    nan2 = posinf / neginf;

    nan3 = posinf / posinf;

    nan4 = posinf * *zero;

    nan5 = neginf * negzro;

    nan6 = nan5 * 0.f;

    if (nan1 == nan1) {
	ret_val = 0;
	return ret_val;
    }

    if (nan2 == nan2) {
	ret_val = 0;
	return ret_val;
    }

    if (nan3 == nan3) {
	ret_val = 0;
	return ret_val;
    }

    if (nan4 == nan4) {
	ret_val = 0;
	return ret_val;
    }

    if (nan5 == nan5) {
	ret_val = 0;
	return ret_val;
    }

    if (nan6 == nan6) {
	ret_val = 0;
	return ret_val;
    }

    return ret_val;
} /* ieeeck_ */





integer template_lapack_ilaenv(const integer *ispec, const char *name__, const char *opts, const integer *n1, 
	const integer *n2, const integer *n3, const integer *n4, ftnlen name_len, ftnlen 
	opts_len)
{

  //printf("entering template_lapack_ilaenv\n");

/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    ILAENV is called from the LAPACK routines to choose problem-dependent   
    parameters for the local environment.  See ISPEC for a description of   
    the parameters.   

    This version provides a set of parameters which should give good,   
    but not optimal, performance on many of the currently available   
    computers.  Users are encouraged to modify this subroutine to set   
    the tuning parameters for their particular machine using the option   
    and problem size information in the arguments.   

    This routine will not function correctly if it is converted to all   
    lower case.  Converting it to all upper case is allowed.   

    Arguments   
    =========   

    ISPEC   (input) INTEGER   
            Specifies the parameter to be returned as the value of   
            ILAENV.   
            = 1: the optimal blocksize; if this value is 1, an unblocked   
                 algorithm will give the best performance.   
            = 2: the minimum block size for which the block routine   
                 should be used; if the usable block size is less than   
                 this value, an unblocked routine should be used.   
            = 3: the crossover point (in a block routine, for N less   
                 than this value, an unblocked routine should be used)   
            = 4: the number of shifts, used in the nonsymmetric   
                 eigenvalue routines   
            = 5: the minimum column dimension for blocking to be used;   
                 rectangular blocks must have dimension at least k by m,   
                 where k is given by ILAENV(2,...) and m by ILAENV(5,...)   
            = 6: the crossover point for the SVD (when reducing an m by n   
                 matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds   
                 this value, a QR factorization is used first to reduce   
                 the matrix to a triangular form.)   
            = 7: the number of processors   
            = 8: the crossover point for the multishift QR and QZ methods   
                 for nonsymmetric eigenvalue problems.   
            = 9: maximum size of the subproblems at the bottom of the   
                 computation tree in the divide-and-conquer algorithm   
                 (used by xGELSD and xGESDD)   
            =10: ieee NaN arithmetic can be trusted not to trap   
            =11: infinity arithmetic can be trusted not to trap   

    NAME    (input) CHARACTER*(*)   
            The name of the calling subroutine, in either upper case or   
            lower case.   

    OPTS    (input) CHARACTER*(*)   
            The character options to the subroutine NAME, concatenated   
            into a single character string.  For example, UPLO = 'U',   
            TRANS = 'T', and DIAG = 'N' for a triangular routine would   
            be specified as OPTS = 'UTN'.   

    N1      (input) INTEGER   
    N2      (input) INTEGER   
    N3      (input) INTEGER   
    N4      (input) INTEGER   
            Problem dimensions for the subroutine NAME; these may not all   
            be required.   

   (ILAENV) (output) INTEGER   
            >= 0: the value of the parameter specified by ISPEC   
            < 0:  if ILAENV = -k, the k-th argument had an illegal value.   

    Further Details   
    ===============   

    The following conventions have been used when calling ILAENV from the   
    LAPACK routines:   
    1)  OPTS is a concatenation of all of the character options to   
        subroutine NAME, in the same order that they appear in the   
        argument list for NAME, even if they are not used in determining   
        the value of the parameter specified by ISPEC.   
    2)  The problem dimensions N1, N2, N3, N4 are specified in the order   
        that they appear in the argument list for NAME.  N1 is used   
        first, N2 second, and so on, and unused problem dimensions are   
        passed a value of -1.   
    3)  The parameter value returned by ILAENV is checked for validity in   
        the calling subroutine.  For example, ILAENV is used to retrieve   
        the optimal blocksize for STRTRI as follows:   

        NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )   
        IF( NB.LE.1 ) NB = MAX( 1, N )   

    ===================================================================== */
    /* Table of constant values */
     integer c__0 = 0;
     real c_b162 = 0.f;
     real c_b163 = 1.f;
     integer c__1 = 1;
    
    /* System generated locals */
    integer ret_val;
    /* Local variables */
     integer i__;
     logical cname, sname;
     integer nbmin;
     char c1[1], c2[2], c3[3], c4[2];
     integer ic, nb;
     integer iz, nx;
     char subnam[6];




    switch (*ispec) {
	case 1:  goto L100;
	case 2:  goto L100;
	case 3:  goto L100;
	case 4:  goto L400;
	case 5:  goto L500;
	case 6:  goto L600;
	case 7:  goto L700;
	case 8:  goto L800;
	case 9:  goto L900;
	case 10:  goto L1000;
	case 11:  goto L1100;
    }

/*     Invalid value for ISPEC */

    ret_val = -1;
    return ret_val;

L100:

/*     Convert NAME to upper case if the first character is lower case. */

    ret_val = 1;
    s_copy(subnam, name__, (ftnlen)6, name_len);
    ic = *(unsigned char *)subnam;
    iz = 'Z';
    if (iz == 90 || iz == 122) {

/*        ASCII character set */

	if (ic >= 97 && ic <= 122) {
	    *(unsigned char *)subnam = (char) (ic - 32);
	    for (i__ = 2; i__ <= 6; ++i__) {
		ic = *(unsigned char *)&subnam[i__ - 1];
		if (ic >= 97 && ic <= 122) {
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
		}
/* L10: */
	    }
	}

    } else if (iz == 233 || iz == 169) {

/*        EBCDIC character set */

      if ( (ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 162 && 
								    ic <= 169) ) {
	    *(unsigned char *)subnam = (char) (ic + 64);
	    for (i__ = 2; i__ <= 6; ++i__) {
		ic = *(unsigned char *)&subnam[i__ - 1];
		if ( (ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 
									      162 && ic <= 169) ) {
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
		}
/* L20: */
	    }
	}

    } else if (iz == 218 || iz == 250) {

/*        Prime machines:  ASCII+128 */

	if (ic >= 225 && ic <= 250) {
	    *(unsigned char *)subnam = (char) (ic - 32);
	    for (i__ = 2; i__ <= 6; ++i__) {
		ic = *(unsigned char *)&subnam[i__ - 1];
		if (ic >= 225 && ic <= 250) {
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
		}
/* L30: */
	    }
	}
    }

    *(unsigned char *)c1 = *(unsigned char *)subnam;
    sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
    cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
    if (! (cname || sname)) {
	return ret_val;
    }
    s_copy(c2, subnam + 1, (ftnlen)2, (ftnlen)2);
    s_copy(c3, subnam + 3, (ftnlen)3, (ftnlen)3);
    s_copy(c4, c3 + 1, (ftnlen)2, (ftnlen)2);

    switch (*ispec) {
	case 1:  goto L110;
	case 2:  goto L200;
	case 3:  goto L300;
    }

L110:

/*     ISPEC = 1:  block size   

       In these examples, separate code is provided for setting NB for   
       real and complex.  We assume that NB will take the same value in   
       single or double precision. */

    nb = 1;

    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	} else if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, 
		"RQF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)
		3, (ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) 
		== 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "PO", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	} else if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
	    nb = 32;
	} else if (sname && s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0) {
	    nb = 64;
	}
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    nb = 64;
	} else if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
	    nb = 32;
	} else if (s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0) {
	    nb = 64;
	}
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nb = 32;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nb = 32;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nb = 32;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nb = 32;
	    }
	}
    } else if (s_cmp(c2, "GB", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		if (*n4 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    } else {
		if (*n4 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    }
	}
    } else if (s_cmp(c2, "PB", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		if (*n2 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    } else {
		if (*n2 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    }
	}
    } else if (s_cmp(c2, "TR", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "LA", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "UUM", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (sname && s_cmp(c2, "ST", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "EBZ", (ftnlen)3, (ftnlen)3) == 0) {
	    nb = 1;
	}
    }
    ret_val = nb;
    return ret_val;

L200:

/*     ISPEC = 2:  minimum block size */

    nbmin = 2;
    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "RQF", (
		ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)3, (
		ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
		 {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	}
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nbmin = 8;
	    } else {
		nbmin = 8;
	    }
	} else if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
	    nbmin = 2;
	}
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
	    nbmin = 2;
	}
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nbmin = 2;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nbmin = 2;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nbmin = 2;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nbmin = 2;
	    }
	}
    }
    ret_val = nbmin;
    return ret_val;

L300:

/*     ISPEC = 3:  crossover point */

    nx = 0;
    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "RQF", (
		ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)3, (
		ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
		 {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	}
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
	if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
	    nx = 32;
	}
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
	    nx = 32;
	}
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nx = 128;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nx = 128;
	    }
	}
    }
    ret_val = nx;
    return ret_val;

L400:

/*     ISPEC = 4:  number of shifts (used by xHSEQR) */

    ret_val = 6;
    return ret_val;

L500:

/*     ISPEC = 5:  minimum column dimension (not used) */

    ret_val = 2;
    return ret_val;

L600:

/*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */

    ret_val = (integer) ((real) minMACRO(*n1,*n2) * 1.6f);
    return ret_val;

L700:

/*     ISPEC = 7:  number of processors (not used) */

    ret_val = 1;
    return ret_val;

L800:

/*     ISPEC = 8:  crossover point for multishift (used by xHSEQR) */

    ret_val = 50;
    return ret_val;

L900:

/*     ISPEC = 9:  maximum size of the subproblems at the bottom of the   
                   computation tree in the divide-and-conquer algorithm   
                   (used by xGELSD and xGESDD) */

    ret_val = 25;
    return ret_val;

L1000:

/*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap   

       ILAENV = 0 */
    ret_val = 1;
    if (ret_val == 1) {
	ret_val = ieeeck_(&c__0, &c_b162, &c_b163);
    }
    return ret_val;

L1100:

/*     ISPEC = 11: infinity arithmetic can be trusted not to trap   

       ILAENV = 0 */
    ret_val = 1;
    if (ret_val == 1) {
	ret_val = ieeeck_(&c__1, &c_b162, &c_b163);
    }
    return ret_val;

/*     End of ILAENV */

} /* ilaenv_ */

