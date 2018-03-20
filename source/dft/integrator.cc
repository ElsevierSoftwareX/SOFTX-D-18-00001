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
/** @file integrator.cc The DFT integrator.
   (c) Pawel Salek, pawsa@theochem.kth.se.
   2001.07.13

   The WRKMEM memory block is not used since it should be deprecated.
   It might be therefore useful to enable memory overcommiting. On linux-2.4.x
   it can be done with echo 1 > /proc/sys/vm/overcommit_memory or a 
   sysctl call. We use it only to pass it to other Fortran routines we call.

   OPTIMIZATIONS: ordinary calculation uses approximately only 4%
   total CPU time in this code. Most likely, the optimizations should
   be sought somewhere else. The simple optimization path is though to
   use block structure of kappa matrices to reduce time by 2 for
   larger matrices.  

   integrator.cc provides dft_integrator() routine. It is passed some
   standard parameters and a table of callbacks and associated
   closures (callback data). The callback gets the grid data for
   current point as well as its own closure.

   The grid file is assumed to be available on call to
   dft_integrator(). Otherwise, it is a black-box implementation.

*/

#include <assert.h>
#include <cmath>
#include <pthread.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#define __CVERSION__
#include "aos.h"
#include "integrator.h"
#include "functionals.h"

#include "output.h"
#include "rho-mat.h"
#include "sparse_matrix.h"
#include "grid_reader.h"
#include "dft_common.h"

/* blocksz_t is a variable that matches the one used by the 
 * fortran runtime library to store the block size. This type
 * is compilator dependent but is usually int or long int.
 */
#if defined(__gnu_linux__)
/* gnu compilers use this */
typedef long blocksz_t;
#else
/* safe default for all the other compilers */
typedef int blocksz_t;
#endif

#define max(a,b) ((a)>(b)? (a):(b))

/** the DFT grid buffer length. grid_getchunk_blocked() will never try
    to read buffers longer than this. */
#define GRID_BUFF_SZ 100000

/* =================================================================== */
/*                     BLOCKED INTEGRATORS                             */
/* =================================================================== */
/* Blocked integrator(s) have altered the block structure to enhance
 * data locality and increase length of internal loops. This should
 * increase performance even for small molecules and reach linear
 * scaling for large ones by enabling vector-like optimization and
 * enhancing data locality.
 */
/* dft_grid_blocked_new:
   initialize grid data.
   ndmat - number of density matrices handled at the same time 
           needed for temporary array.
   bllen - grid point batch length.
*/
DftIntegratorBl*
dft_integrator_bl_new(Functional* f, int ndmat,
                      int bllen, int needlondon, const BasisInfoStruct& bis)
{
    int kmax, nbast;
    DftIntegratorBl* grid = new DftIntegratorBl;

    grid->coor   = (real (*)[3])calloc(3*GRID_BUFF_SZ, sizeof(real));
    grid->weight = (real*)calloc(GRID_BUFF_SZ, sizeof(real));
    grid->dogga  = f->is_gga();
    grid->dfthri = 1e-13;
    grid->needlap= 0;
    grid->needgb = needlondon;
    grid->nsym   = 1;
    grid->ndmat  = ndmat;
    grid->nbast  = nbast = bis.noOfBasisFuncs;

    kmax = bis.noOfShells;
#ifdef DALTON
    geodrv = grid->dogga ? 1 : 0;
    setupsos_(&geodrv, &grid->needgb, &grid->ntypso, &grid->london_off);
    grid->london_off--; /* convert from fortran offset type */
#else
    grid->ntypso = grid->dogga ? 4 : 1;
#endif
    grid->atv   = (real*)calloc(bllen*nbast*grid->ntypso, sizeof(real));
    grid->shlblocks = (int (*)[2])dal_malloc(2*kmax*sizeof(int));
    grid->basblocks = (int (*)[2])dal_malloc(2*kmax*8*sizeof(int));

    /* Allocate memory for rho, taking advantage of the union. */
    grid->r.rho    = dal_new(ndmat*bllen,real);
    grid->g.grad   = (real (*)[3])dal_malloc(ndmat*3*bllen*sizeof(real));
    
    /* and set some aliases in case somebody needed them for open-shell. *
     * Observe that rho aliases with rhoa. */
    if(ndmat == 2) {
        grid->r.ho.b    = grid->r.ho.a  + bllen;
        grid->g.rad.b   = grid->g.rad.a + bllen;
    }
    return grid;
}

void
dft_integrator_bl_free(DftIntegratorBl *res)
{
  free(res->coor);
  free(res->weight);
  free(res->atv);
  free(res->r.rho);
  free(res->g.grad);
  free(res->shlblocks);
  free(res->basblocks);
  delete res;
}

/* grid_blocked_getval:
   evaluates a block of bllen orbitals
*/
void blgetsos_(int *nvclen, real GSO[], real COOR[],
               int *NBLCNT, int IBLCKS[], real WORK[], int *LWORK,
               int *NBAST, int *DOLND, int *DOGGA, real *DFTHRI,
               const int*IPRINT);


static void
output_memory_usage(bool& a)
{
    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;

    pthread_mutex_lock(&m);
    if(!a){
        output_current_memory_usage(LOG_AREA_DFT, "XC integration loop");
        a = true;
    }
    pthread_mutex_unlock(&m);
}

template<typename DensityType>
static real
dft_integrate(int ndmat,
              const DensityType* const * dmat,
              const BasisInfoStruct& bis,
              const Molecule& mol, const Dft::GridParams& gss,
              int nThreads, DftBlockCallback cb, void *cb_data)
{
    int npoints, blocksz;
    int nbast = bis.noOfBasisFuncs;
    real electrons; /* alpha electrons only most of the time */
    DftIntegratorBl* grid;
    real *dmagao;
    DftGridReader* rawgrid;
    int nder;
    static bool firstThreadFlag;

    firstThreadFlag = false;
    ErgoMolInfo mol_info(bis, mol);

    dmagao = dal_new(nbast*DFT_MAX_BLLEN,real);
    grid = dft_integrator_bl_new(selected_func, ndmat,
                                 DFT_MAX_BLLEN, false, bis);

    /* start integration */
    electrons  = 0.0;
    nder = grid->dogga ? 1 : 0;

    Dft::Matrix *mat = createGridMatrix(*dmat[0]);
    rawgrid = grid_open_full(&mol_info, gss, NULL, mat, bis);
    delete mat;

    if(sync_threads(false, nThreads) != 0) throw "Error syncing threads";

    output_memory_usage(firstThreadFlag);
    npoints = 0;
    while( (blocksz=grid_getchunk_blocked(rawgrid, GRID_BUFF_SZ,
                                          &grid->shl_bl_cnt, 
                                          &grid->shlblocks[0][0],
                                          &grid->coor[0],
                                          grid->weight)) >=0) {
        int ipnt;
        ergoShellsToOrbs(&grid->shl_bl_cnt, grid->shlblocks, 
                         grid->bas_bl_cnt, grid->basblocks,
                         bis);

        for(ipnt=0; ipnt<blocksz; ipnt+=DFT_MAX_BLLEN) {
            int i, j, lo, hi;
            int len = ipnt+DFT_MAX_BLLEN<blocksz ? DFT_MAX_BLLEN : blocksz-ipnt;
            grid->curr_point  = ipnt;

            dft_get_orbs(len, grid->atv, (real(*)[3]) &grid->coor[ipnt][0],
                         grid->shl_bl_cnt, (int(*)[2]) &grid->shlblocks[0][0],
                         nder, bis);

            for(i=0; i<ndmat; i++) {
                int roff = i*DFT_MAX_BLLEN;
                if(grid->dogga)
                  getrho_blocked_gga(nbast, *dmat[i], grid->atv,
                                     grid->bas_bl_cnt,
                                     grid->basblocks, grid->shl_bl_cnt,
                                     dmagao, len, grid->r.rho+roff,
                                     grid->g.rad.a+roff);
                else
                  getrho_blocked_lda(nbast, *dmat[i], grid->atv,
                                     grid->bas_bl_cnt,
                                     grid->basblocks, grid->shl_bl_cnt,
                                     dmagao, len, grid->r.rho+roff);
                for(j=0; j<len; j++)
                  electrons += grid->weight[ipnt+j]*grid->r.rho[j+roff];
            }
            lo = 0; hi = len;
            npoints += len;
            /* Consider skipping low-density points at the beginning
               and the end of the batch by modifying lo and hi. */
            if(lo<hi)
                cb(grid, dmagao, len, lo, hi, cb_data);
        }
    }
    grid_close(rawgrid);
#ifdef VAR_MPI
    FSYM(dftintcollect)(&electrons);
#endif
    free(dmagao);
    dft_integrator_bl_free(grid);
    return electrons;
}

BEGIN_NAMESPACE(Dft)
/** reads the grid and calls the callback function for each group of
    grid points. As a courtesy, evaluates first for each batch the
    density.
    
    @param ndmat number of density matrices to evaluate
    @param dmat square density matrices.
    @param bis  a structure describing the used basis set.
    @param mol  a structure describing the molecule.

    @param gss  a structure describing the grid settings.

    @param nThreads - how many threads will execute this function
    simultaneously. Needed for synchronisation purposes.

    @param cb function to be evaluated for each batch of grid points
    @param cb_data its closure.
*/
real
integrate(int ndmat, const FullMatrix * const *dmat,
               const BasisInfoStruct& bis,
               const Molecule& mol, const Dft::GridParams& gss,
               int nThreads, DftBlockCallback cb, void *cb_data)
{
  return dft_integrate(ndmat, dmat, bis, mol, gss, nThreads, cb, cb_data);
}

real
integrate(int nDmat, const SparseMatrix * const *dmat,
          const BasisInfoStruct& bis,
          const Molecule& mol, const Dft::GridParams& gss,
          int nThreads, DftBlockCallback cb, void *cb_data)
{
  return dft_integrate(nDmat, dmat, bis, mol, gss, nThreads, cb, cb_data);
}

END_NAMESPACE(Dft)
