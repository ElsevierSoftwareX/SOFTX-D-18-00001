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
/** @file fun-gga.c
   Implementation of a functional being a linear combination of 
   other functionals.
   (c) Pawel Salek, pawsa@theochem.kth.se, sep 2001
*/

#include <ctype.h>
//#include <math.h>
#include <stdio.h>
#include <string.h>
#include <strings.h> /* Needed for strncasecmp() */
#include <stdlib.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int  lda_read(const char* conf_line);
static real lda_energy(const FunDensProp* dp);
static void lda_first(FunFirstFuncDrv *ds,   real fac, const FunDensProp* dp);
static void lda_second(FunSecondFuncDrv *ds, real fac, const FunDensProp* dp);
static void lda_third(FunThirdFuncDrv *ds,   real fac, const FunDensProp* dp);
static void lda_fourth(FunFourthFuncDrv *ds, real fac, const FunDensProp* dp);
static int  ldagauss_read(const char* conf_line);
static int  blyp_read(const char* conf_line);
static int  b3lyp_read(const char* conf_line);
static int  b3lypgauss_read(const char* conf_line);
static int  b3pw91_read(const char* conf_line);
static int  bhandh_read(const char* conf_line);
static int  bhandhlyp_read(const char* conf_line);
static int  bp86_read(const char* conf_line);
static int  bpw91_read(const char* conf_line);
static int  b3p86_read(const char* conf_line);
static int  b3p86g_read(const char* conf_line);
static int  kt1_read(const char* conf_line);
static int  kt2_read(const char* conf_line);
static int  kt3_read(const char* conf_line);
static int  olyp_read(const char* conf_line);
static int  pbe_read(const char* conf_line);
static int  pbe0_read(const char* conf_line);

static int  gga_isgga(void);
static int  xalpha_read(const char* conf_line);
static int  combine_read(const char* conf_line);
static void gga_report(void);
static real gga_energy(const FunDensProp* dp);
static void gga_first(FunFirstFuncDrv *ds,   real fac, const FunDensProp* dp);
static void gga_second(FunSecondFuncDrv *ds, real fac, const FunDensProp* dp);
static void gga_third(FunThirdFuncDrv *ds,   real fac, const FunDensProp* dp);
static void gga_fourth(FunFourthFuncDrv *ds, real fac, const FunDensProp* dp);

#define LDA_FUNCTIONAL(name,read) { (name), \
    fun_false, (read), NULL, lda_energy, lda_first, lda_second, \
    lda_third, lda_fourth }

#define GGA_FUNCTIONAL(name,read) { (name), \
    gga_isgga, (read), gga_report, gga_energy, gga_first, gga_second, \
    gga_third, gga_fourth }

Functional XAlphaFunctional = GGA_FUNCTIONAL("XAlpha", xalpha_read);
Functional LDAFunctional =    LDA_FUNCTIONAL("LDA",     lda_read);
/* SVWN5 aliases LDA */
Functional SVWN5Functional =  LDA_FUNCTIONAL("SVWN5",   lda_read);
Functional SVWN3Functional =  GGA_FUNCTIONAL("SVWN3",   ldagauss_read);
Functional B3LYPFunctional =  GGA_FUNCTIONAL("B3LYP",   b3lyp_read);
Functional B3LYPGaussFunctional = GGA_FUNCTIONAL("B3LYP-G", b3lypgauss_read);
Functional B3P86Functional =  GGA_FUNCTIONAL("B3P86",   b3p86_read);
Functional B3P86GFunctional = GGA_FUNCTIONAL("B3P86-G", b3p86g_read);
Functional B3PW91Functional =  GGA_FUNCTIONAL("B3PW91",  b3pw91_read);
Functional BHandHFunctional = GGA_FUNCTIONAL("BHandH",  bhandh_read);
Functional BHandHLYPFunctional = GGA_FUNCTIONAL("BHandHLYP", bhandhlyp_read);
Functional BLYPFunctional =   GGA_FUNCTIONAL("BLYP",    blyp_read);
Functional BP86Functional =   GGA_FUNCTIONAL("BP86",    bp86_read);
Functional BPW91Functional =  GGA_FUNCTIONAL("BPW91",   bpw91_read);
Functional CombineFunctional =GGA_FUNCTIONAL("Combine", combine_read);
Functional GGAKeyFunctional = GGA_FUNCTIONAL("GGAKey",  combine_read);
Functional KT1Functional =    GGA_FUNCTIONAL("KT1",     kt1_read);
Functional KT2Functional =    GGA_FUNCTIONAL("KT2",     kt2_read);
Functional KT3Functional =    GGA_FUNCTIONAL("KT3",     kt3_read);
Functional OLYPFunctional =   GGA_FUNCTIONAL("OLYP" ,   olyp_read);
Functional PBE0Functional =   GGA_FUNCTIONAL("PBE0",    pbe0_read);
Functional PBEFunctional =    GGA_FUNCTIONAL("PBE",     pbe_read); 

/* MIXED FUNCTIONALS */
typedef struct FuncList_ FuncList;

struct FuncList_ {
    Functional* func;
    real        weight;
    FuncList*   next;
};

FuncList* gga_fun_list = NULL;

static int
gga_isgga(void)
{ 
    int res = 0;

    FuncList* lst;
    for(lst=gga_fun_list; lst && !res; lst=lst->next) 
        res |= lst->func->is_gga();
    return res;
}

static FuncList*
add_functional(FuncList* lst, Functional* f, float weight)
{
    FuncList* n = (FuncList*)malloc(sizeof(FuncList));
    n->func = f; n->weight = weight; n->next = lst;
    return n;
}

static void
free_functionals(FuncList **lst)
{
    
    FuncList *item = *lst;
    while (item) {
        FuncList *n = item->next;
        free(item);
        item = n;
    }
    *lst = NULL;
}
       

static int
xalpha_read(const char* conf_line)
{
    float weight;
    int res = (sscanf(conf_line, "%g", &weight)==1);
    free_functionals(&gga_fun_list);
    if(res) 
        gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 
                                      1.5*weight);
    fun_set_hf_weight(0);
    return res;
}

static int
lda_read(const char* conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}

static real
lda_energy(const FunDensProp* dp)
{
    return SlaterFunctional.func(dp) + VWNFunctional.func(dp);
}

static void
lda_first(FunFirstFuncDrv *ds, real factor,  const FunDensProp* dp)
{
    SlaterFunctional.first(ds, factor, dp);
    VWNFunctional  .first(ds, factor, dp);
}

static void
lda_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    SlaterFunctional.second(ds, factor, dp);
    VWN5Functional  .second(ds, factor, dp);
}

static void
lda_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    SlaterFunctional.third(ds, factor, dp);
    VWN5Functional  .third(ds, factor, dp);
}

static void
lda_fourth(FunFourthFuncDrv *ds, real factor, const FunDensProp* dp)
{
    SlaterFunctional.fourth(ds, factor, dp);
    VWN5Functional  .fourth(ds, factor, dp);
}

static int
ldagauss_read(const char* conf_line)
{
    free_functionals(&gga_fun_list);
    gga_fun_list = add_functional(NULL, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &VWN3Functional,  1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
blyp_read(const char* conf_line)
{
    free_functionals(&gga_fun_list);
    gga_fun_list = add_functional(NULL, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,   1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
b3lyp_read(const char* conf_line)
{
    static const real lypw = 0.81, dirw = 0.8;
    free_functionals(&gga_fun_list);
    gga_fun_list = add_functional(NULL, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  0.72);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    lypw);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,    1-lypw);
    fun_set_hf_weight(1-dirw);
    return 1;
}

static int
b3lypgauss_read(const char* conf_line)
{
    static const real lypw = 0.81, dirw = 0.8;
    free_functionals(&gga_fun_list);
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  0.72);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    lypw);
    gga_fun_list = add_functional(gga_fun_list, &VWN3IFunctional,   1-lypw);
    fun_set_hf_weight(1-dirw);
    return 1;
}

static int
b3pw91_read(const char* conf_line)
{
#if 0
    return combine_read("hf=0.2 slater=0.8 becke=0.72 vwn=0.19 pw91c=0.81");
#else
    static const real corw = 0.81, dirw = 0.8;
    free_functionals(&gga_fun_list);
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,   dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,    0.72);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,      1-corw);
    gga_fun_list = add_functional(gga_fun_list, &Pw91cFunctional,    corw);
    fun_set_hf_weight(1-dirw);
    return 1;
#endif
}

static int
bhandh_read(const char* conf_line)
{
    free_functionals(&gga_fun_list);
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 0.5);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    1.0);
    fun_set_hf_weight(0.5);
    return 1;
}

static int
bhandhlyp_read(const char* conf_line)
{
    free_functionals(&gga_fun_list);
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 0.5);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  0.5);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    1.0);
    fun_set_hf_weight(0.5);
    return 1;
}

static int
bp86_read(const char* conf_line)
{
    free_functionals(&gga_fun_list);
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &PZ81Functional,  1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
b3p86_read(const char* conf_line)
{
    static const real lypw = 0.81, dirw = 0.8;
    free_functionals(&gga_fun_list);
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  0.72);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,   lypw);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,    1);
    fun_set_hf_weight(1-dirw);
    return 1;
}

static int
b3p86g_read(const char* conf_line)
{
    static const real lypw = 0.81, dirw = 0.8;
    free_functionals(&gga_fun_list);
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional,  dirw);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  0.72);
    gga_fun_list = add_functional(gga_fun_list, &P86cFunctional,   lypw);
    gga_fun_list = add_functional(gga_fun_list, &VWN3Functional,   1);
    fun_set_hf_weight(1-dirw);
    return 1;
}

static int
bpw91_read(const char* conf_line)
{
    free_functionals(&gga_fun_list);
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1);
    gga_fun_list = add_functional(gga_fun_list, &BeckeFunctional,  1);
    gga_fun_list = add_functional(gga_fun_list, &Pw91cFunctional,  1);
    fun_set_hf_weight(0);
    return 1;
}

static int
kt1_read(const char* conf_line)
{
    static const real ktgam = -0.006;
    free_functionals(&gga_fun_list);
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &KTFunctional,    ktgam);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,   1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
kt2_read(const char* conf_line)
{   
    static const real dirw = 1.07173, vwnw = 0.576727; 
    static const real ktgam = -0.006;
    free_functionals(&gga_fun_list);
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw);
    gga_fun_list = add_functional(gga_fun_list, &KTFunctional,    ktgam);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,   vwnw);
    fun_set_hf_weight(0);
    return 1;
}

static int
kt3_read(const char* conf_line)
{
    static const real dirw = 1.092, lypw = 0.864409, optw = -0.925452;
    static const real ktgam = -0.004;
    free_functionals(&gga_fun_list);
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw);
    gga_fun_list = add_functional(gga_fun_list, &KTFunctional,    ktgam);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,   lypw);
    gga_fun_list = add_functional(gga_fun_list, &OPTXFunctional, optw);
    fun_set_hf_weight(0);
    return 1;
}

static int
olyp_read(const char* conf_line)
{
    static const real optkw = -1.43169, dirw = 1.05151;
    free_functionals(&gga_fun_list);
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw);
    gga_fun_list = add_functional(gga_fun_list, &OPTXFunctional, optkw);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,   1.0);
    fun_set_hf_weight(0);
    return 1;
}

static int
pbe_read(const char* conf_line)
{
    free_functionals(&gga_fun_list);
    gga_fun_list = add_functional(gga_fun_list, &PbecFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &PbexFunctional, 1.0);
    fun_set_hf_weight(0);
    return 1;
}
 
static int
pbe0_read(const char* conf_line)
{
    free_functionals(&gga_fun_list);
    gga_fun_list = add_functional(gga_fun_list, &PbecFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &PbexFunctional, 0.75);
    fun_set_hf_weight(0.25);
    return 1;
}


/* combine_read:
   read general GGA input: 
*/
static int
combine_read(const char* conf_line)
{
    int res = 1, i;
    float f;
    const char* str = conf_line;

    fun_set_hf_weight(0);
    free_functionals(&gga_fun_list);

    while(*str) {
        while(*str && isspace((int)*str)) str++; /* skip whitespace */
        if(*str =='\0') break; /* line ended by whitespace */
        if(strncasecmp("HF=", str, 3)==0) {
            if(sscanf(str+3,"%g", &f) != 1) {
                fun_printf("Combine: HF not followed by the weight: ",
                           conf_line);
                res = 0;
            } else fun_set_hf_weight(f);
        } else {
            for(i=0; available_functionals[i]; i++) {
                int len = strlen(available_functionals[i]->name);
                if(strncasecmp(available_functionals[i]->name, str, len)==0 &&
                   str[len] == '=') {
                    if(sscanf(str+len+1,"%g", &f) != 1) {
                        fun_printf("GGAKey: keyword '%s' not followed by "
                                   "weight: %s",
                                   available_functionals[i]->name, 
                                   conf_line);
                        res = 0;
                    } else {
                        gga_fun_list = 
                            add_functional(gga_fun_list, 
                                           available_functionals[i], f);
                        break; /* weight properly read, break the 'for' loop */
                    }
                }
            }  
            if(available_functionals[i] == NULL) {
                fun_printf("Combine: functional '%s' not recognised: ", str);
                res = 0;
            }
        }
        while(*str && !isspace((int)*str)) str++; /* skip nonws */
    } return res;
}

static void
gga_report(void)
{
    FuncList* lst;
    fun_printf("Weighted mixed functional:");
    if(fun_get_hf_weight()>0)
      fun_printf("%25s: %10.5f", "HF exchange", (double)fun_get_hf_weight());
    for(lst=gga_fun_list; lst; lst=lst->next) 
      fun_printf("%25s: %10.5f", lst->func->name, (double)lst->weight);
}

static real
gga_energy(const FunDensProp* dp)
{
    real res = 0;
    FuncList* lst;
    for(lst=gga_fun_list; lst; lst=lst->next) {
        real contr = lst->weight*lst->func->func(dp);
/*        fun_printf("[%g,%g,w=%g] %s contributes with %g", dp->rhoa,
                   dp->grada, lst->weight, lst->func->name, contr); */
        res += contr;
    }
    return res;
}

static void
gga_first(FunFirstFuncDrv *ds, real factor,  const FunDensProp* dp)
{
    FuncList* lst;
    for(lst=gga_fun_list; lst; lst=lst->next) {
        lst->func->first(ds, factor*lst->weight, dp);
/*      fun_printf("[%g,%g,w=%g] %s f: %g deriv (%g,%g)", dp->rhoa,
                   dp->grada, lst->weight, lst->func->name, factor,
		   ds->df1000-df10, ds->df0010-df01); */
    }
}

static void
gga_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    FuncList* lst;
    for(lst=gga_fun_list; lst; lst=lst->next) 
        lst->func->second(ds, factor*lst->weight, dp);
}

static void
gga_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    FuncList* lst;
    for(lst=gga_fun_list; lst; lst=lst->next) 
        lst->func->third(ds, factor*lst->weight, dp);
}


static void
gga_fourth(FunFourthFuncDrv *ds, real factor, const FunDensProp* dp)
{
    FuncList* lst;
    for(lst=gga_fun_list; lst; lst=lst->next) 
        lst->func->fourth(ds, factor*lst->weight, dp);
}


