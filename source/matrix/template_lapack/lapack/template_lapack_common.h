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
 

#ifndef TEMPLATE_LAPACK_COMMON_HEADER
#define TEMPLATE_LAPACK_COMMON_HEADER

#include "template_blas_common.h"

#define TRUE_ 1
#define FALSE_ 0

integer template_lapack_ilaenv(const integer *ispec, const char *name__, const char *opts, const integer *n1, 
			       const integer *n2, const integer *n3, const integer *n4, ftnlen name_len, ftnlen 
			       opts_len);



#include "template_lapack_lamch.h"


#include "template_lapack_lae2.h"
#include "template_lapack_lascl.h"
#include "template_lapack_lansy.h"
#include "template_lapack_pptrf.h"
#include "template_lapack_spgst.h"
#include "template_lapack_tptri.h"
#include "template_lapack_potrf.h"
#include "template_lapack_potf2.h"
#include "template_lapack_sygst.h"
#include "template_lapack_sygs2.h"
#include "template_lapack_larfg.h"
#include "template_lapack_latrd.h"
#include "template_lapack_sytd2.h"
#include "template_lapack_sytrd.h"
#include "template_lapack_lanst.h"
#include "template_lapack_lapy2.h"
#include "template_lapack_lasrt.h"
#include "template_lapack_laset.h"
#include "template_lapack_sterf.h"
#include "template_lapack_larft.h"
#include "template_lapack_larfb.h"
#include "template_lapack_larf.h"
#include "template_lapack_org2r.h"
#include "template_lapack_orgqr.h"
#include "template_lapack_org2l.h"
#include "template_lapack_orgql.h"
#include "template_lapack_orgtr.h"
#include "template_lapack_laev2.h"
#include "template_lapack_lasr.h"
#include "template_lapack_lartg.h"
#include "template_lapack_steqr.h"
#include "template_lapack_syev.h"
#include "template_lapack_sygv.h"
#include "template_lapack_trti2.h"
#include "template_lapack_trtri.h"
#include "template_lapack_ladiv.h"
#include "template_lapack_laln2.h"
#include "template_lapack_lapy3.h"
#include "template_lapack_lasv2.h"
#include "template_lapack_lag2.h"
#include "template_lapack_lanhs.h"
#include "template_lapack_geqr2.h"
#include "template_lapack_ggbak.h"
#include "template_lapack_tgevc.h"
#include "template_lapack_hgeqz.h"
#include "template_lapack_gghrd.h"
#include "template_lapack_lacpy.h"
#include "template_lapack_orm2r.h"
#include "template_lapack_ormqr.h"
#include "template_lapack_geqrf.h"
#include "template_lapack_ggbal.h"
#include "template_lapack_labad.h"
#include "template_lapack_lange.h"
#include "template_lapack_ggev.h"
#include "template_lapack_rscl.h"
#include "template_lapack_latrs.h"
#include "template_lapack_lacon.h"
#include "template_lapack_pocon.h"
#include "template_lapack_laruv.h"
#include "template_lapack_laebz.h"
#include "template_lapack_lagts.h"
#include "template_lapack_lagtf.h"
#include "template_lapack_larnv.h"
#include "template_lapack_stein.h"
#include "template_lapack_stebz.h"
#include "template_lapack_stevx.h"
#include "template_lapack_larra.h"
#include "template_lapack_larrb.h"
#include "template_lapack_larrc.h"
#include "template_lapack_larrd.h"
#include "template_lapack_larre.h"
#include "template_lapack_larrf.h"
#include "template_lapack_larrj.h"
#include "template_lapack_larrk.h"
#include "template_lapack_larrr.h"
#include "template_lapack_larrv.h"
#include "template_lapack_lar1v.h"
#include "template_lapack_laneg.h"
#include "template_lapack_isnan.h"
#include "template_lapack_laisnan.h"
#include "template_lapack_lasq2.h"
#include "template_lapack_lasq3.h"
#include "template_lapack_lasq4.h"
#include "template_lapack_lasq5.h"
#include "template_lapack_lasq6.h"
#include "template_lapack_stemr.h"
#include "template_lapack_stevr.h"
#include "template_lapack_laswp.h"
#include "template_lapack_getf2.h"
#include "template_lapack_getrf.h"
#include "template_lapack_getrs.h"
#include "template_lapack_gesv.h"

#endif
