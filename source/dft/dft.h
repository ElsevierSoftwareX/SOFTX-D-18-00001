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
/** @file dft.h Definitions exported by the DFT module.
    Specific to full matrices, containing traces of Fortran influence
    and really deprecated...

   (c) Pawel Salek, pawsa@theochem.kth.se, feb 2002
*/
#ifndef _GENERAL_H_
#define _GENERAL_H_

#include <stdlib.h>

#ifdef __cplusplus
#define EXTERN_C extern "C"
#else
#define EXTERN_C
#endif

#if defined(__cplusplus)
#define RESTRICT
#else
/* We really do want to take advantage of the restrict keyword... */
#define RESTRICT restrict
#endif

#if !defined(__CVERSION)
#define __CVERSION__
#endif

#include "functionals.h"

#include "basisinfo.h"
#include "molecule.h"
#include "grid_reader.h"

/* Match Fortran name mangling. If the Fortran compiler does not
 * mangle names, define NO_UNDERSCORE in CFLAGS.  g77 and compaq fort
 * (cryptically referred to with HAVE_GCPP below) for linux-alpha both
 * insert a second underscore if routine name contains at least one
 * underscore /hjaaj Oct04 */
#ifdef NO_UNDERSCORE
#define FSYM(a) a
#define FSYM2(a) a
#else
#define FSYM(a) a ## _
#if defined(VAR_G77) || defined(HAVE_GCPP)
#define FSYM2(a) a ## __
#else
#define FSYM2(a) a ## _
#endif
#endif

#if defined(VAR_PGF77)
#define __FUNCTION__ "PGI_does_not_define__FUNCTION__"
#endif
#if defined(SYS_SUN)
#define __FUNCTION__ "SUNs CC compiler_does_not_define__FUNCTION__"
#endif
#if defined(SYS_IRIX)
#define __FUNCTION__ "SGIs CC compiler_does_not_define__FUNCTION__"
#endif
#if defined(SYS_DEC)
#define __FUNCTION__ "DEC CC compiler does not define __FUNCTION__"
#endif

#define ELEMENTS(arr) (sizeof(arr)/sizeof(arr[0]))


EXTERN_C void dftpot0_(FirstDrv *ds, const real* weight, const FunDensProp* dp);
EXTERN_C void dftpot1_(SecondDrv *ds, const real* w, const FunDensProp* dp,
		       const int* triplet);

EXTERN_C int dft_setfunc(const char *line);
EXTERN_C void grid_set_tmpdir(const char *tmpdir);


EXTERN_C real dft_get_xc(int nElectrons, const real* dmat,
			 const BasisInfoStruct *bis, const Molecule *mol,
			 const Dft::GridParams& gss,
                         real* ksm, real* edfty,
                         int nThreads);
EXTERN_C real dft_get_uxc(int nElectrons,
                          const real* dmata, const real *dmatb,
			  const BasisInfoStruct *bis, const Molecule *mol,
			  const Dft::GridParams& gss,
			  real* xca,   real *xcb, real* edfty,
                          int nThreads);

/* Property evaluators */
typedef void (*DFTPropEvalMaster)(void);
typedef void (*DFTPropEvalSlave)(real* work, int* lwork, const int* iprint);

extern int (*fort_print)(const char* format, ...);


#if !defined __inline__
/* inline some stuff whenever possible */
#define __inline__
#endif

#endif
