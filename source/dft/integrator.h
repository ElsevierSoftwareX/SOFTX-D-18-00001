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

/** @file integrator.h The DFT integrator interface.
    Pawel Salek.
*/

#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include "basisinfo.h"
#include "matrix_typedefs.h"
#include "grid_stream.h"
#include "functionals.h"

typedef ergo_real      real;
typedef ergo_long_real long_real;

/* =================================================================== */
/*                     BLOCKED INTEGRATORS                             */
/* =================================================================== */

typedef struct DftIntegratorBl_ {
    /* private to integrator */
    real (*coor)[3];
    real* weight;
    real* atv; /* the orbital values and their derivatives at given 
                * grid point. The vector is indexed by dftinf_.kso1, etc
                * the dimensioning is (C syntax) [ntypso][nbast][bllen].
                */
    real dfthri; /* threshold on orbital values */
    int nsym, shl_bl_cnt, bas_bl_cnt[8];
    int (*shlblocks)[2]; /* shell blocks   */
    int (*basblocks)[2]; /* basis function blocks */
#define BASBLOCK(grid,isym) ((grid)->basblocks + (isym)*(grid)->shl_bl_cnt)

    int ntypso; /* how many different vectors are computed for each
                 * (point,orbital) pair. i.e whether only orbital values
                 * are computed (1), orbital values and first derivatives 
                 * (4), etc. */

    int london_off; /* offset of the "london" orbital derivatives */
    /* 1 - only values; 4 - values + (x,y,z) derivatives, etc */

    int ndmat; /* 1 for closed shell, 2 for open shell */
    int nbast; /* number of basis functions */
    /* for closed shell, only rho is set. For open shell, only rhoa and rhob
     * is set. */
    union {
        real *rho; /* total density vector; used in closed shell code. */
        struct {  /* used in open-shell code.    */
            real *a, *b;
        }ho;
    }r;
    union {
        real (*grad)[3]; /*total density gradient; used in closed shell code.*/
        struct {
            real (*a)[3], (*b)[3];
        }rad;
    }g;
    /* public, read only */
    real tgrad[3];/* alpha, also used in closed-shell code */
    int  curr_point;  /* index of the current point */
    real curr_weight; /* the weight at current grid point */
    int dogga, needlap, needgb;
} DftIntegratorBl;

/* dft_integrate_ao_bl:
   numerical integration in atomic orbitals, blocked scheme.
*/
typedef void (*DftBlockCallback)(DftIntegratorBl* grid, real *tmp, 
                                 int bllen, int blstart, int blend,
                                 void* cb_data);

DftIntegratorBl*
dft_integrator_bl_new(Functional* f, int ndmat,
                      int bllen, int needlondon, const BasisInfoStruct& bis);

void
dft_integrator_bl_free(DftIntegratorBl *res);

class Molecule;
namespace Dft {
class FullMatrix;
class SparseMatrix;

real integrate(int ndmat, const FullMatrix * const*dmat,
               const BasisInfoStruct& bis,
               const Molecule& mol, const Dft::GridParams& gss,
               int nThreads, DftBlockCallback cb, void *cb_data);

real integrate(int nDmat, const SparseMatrix * const *dmat,
               const BasisInfoStruct& bis,
               const Molecule& mol, const Dft::GridParams& gss,
               int nThreads, DftBlockCallback cb, void *cb_data);

}

#endif /* _INTEGRATOR_H_ */
