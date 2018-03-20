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
/** @file functionals.c General functional selection and input processing.
   (c) Pawel Salek, pawsa@theochem.kth.se, 2001-08-02

 */

#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define __CVERSION__

#include "functionals.h"

Functional* available_functionals[] = {
    /* generic functionals */
    &BeckeFunctional,
    &KTFunctional,
    &LB94Functional,
    &LYPFunctional,
    &OPTXFunctional,
    &P86cFunctional,
    &PW86xFunctional,
    &Pw91xFunctional,
    &Pw91cFunctional,
    &PW92cFunctional,
    &PZ81Functional,
    &PbecFunctional,
    &PbexFunctional,
    &SlaterFunctional,
    &VWN3Functional,
    &VWN5Functional,
    &VWNIFunctional,
    &VWN3IFunctional,
    &VWNFunctional,
    &XAlphaFunctional,
    /* mixed functionals */
    &B3LYPFunctional,
    &B3LYPGaussFunctional,
    &B3P86Functional,
    &B3P86GFunctional,
    &B3PW91Functional,
    &BHandHFunctional,
    &BHandHLYPFunctional,
    &BLYPFunctional,
    &BP86Functional,
    &BPW91Functional,
    &Camb3lypFunctional,
    &CamFunctional,
    &CombineFunctional,
    &GGAKeyFunctional,
    &HseFunctional,
    &KT1Functional,
    &KT2Functional,
    &KT3Functional,
    &LDAFunctional,
    &OLYPFunctional,
    &PBE0Functional,
    &PBEFunctional,
    &SVWN3Functional,
    &SVWN5Functional,
    NULL
};

static int my_printf(const char *fmt, ...)
{
    int i;va_list ap; va_start(ap, fmt); i= vprintf(fmt, ap); va_end(ap);
    puts("");
    return i;
}
 
static void set_hf_weight(real w)         {}
static real get_hf_weight(void)           {return 0;}
static void set_cam_param(real w, real b) {}

Functional* selected_func = &LDAFunctional;
int (*fun_printf)(const char *fmt, ...) = my_printf;
void (*fun_set_hf_weight)(real w)         = set_hf_weight;
real (*fun_get_hf_weight)(void)           = get_hf_weight;
void (*fun_set_cam_param)(real w, real b) = set_cam_param;

/* =================================================================== */
enum FunError
fun_select_by_name(const char *conf_string)
{
    int ok, i;
    char func_name[40];

    sscanf(conf_string,"%39s", func_name);
    for(i=0; available_functionals[i]; i++)
        if(strcasecmp(available_functionals[i]->name, func_name)==0) {
            selected_func = available_functionals[i];
            ok = selected_func->read ?
                selected_func->read(conf_string+strlen(func_name)) : 1;
            return ok ? FUN_OK : FUN_CONF_ERROR;
        }
    return FUN_UNKNOWN;
}

void
drv1_clear(FunFirstFuncDrv* gga)
{
    memset(gga, 0, sizeof(*gga));
}

void
drv2_clear(FunSecondFuncDrv* gga)
{
    memset(gga, 0, sizeof(*gga));
}

void
drv3_clear(FunThirdFuncDrv* gga)
{
    memset(gga, 0, sizeof(*gga));
}
void
drv4_clear(FunFourthFuncDrv* gga)
{
    memset(gga, 0, sizeof(*gga));
}

int fun_true(void)  { return 1; }
int fun_false(void) { return 0; }

/* =================================================================== */
/*           fortran (and not only) functional stub routines           */
/* =================================================================== */

/** dftreport:
   report the selected functional and its configuration.
*/
void
dftreport_(void)
{
    fun_printf("     This is a DFT calculation of type: %s",
               selected_func->name);
    if(selected_func->report)
        selected_func->report();
}

void
dftlistfuncs_(void)
{
    int i;
    fun_printf("\nAvailable functionals:");
    for(i=0; available_functionals[i]; i++)
        fun_printf(available_functionals[i]->name);
}

void
dftlistfuncs_using_printf_(void)
{
    int i;
    printf("\nAvailable functionals:\n");
    for(i=0; available_functionals[i]; i++)
      printf("%s\n", available_functionals[i]->name);
}

/* declare both known fortran name-mangled variants */
int
dft_isgga_(void)
{ return selected_func->is_gga(); }

int
dft_isgga__(void)
{ return selected_func->is_gga(); }
