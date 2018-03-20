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

/** @file lin_trans.cc

    @brief Blocked DFT Linear Response contribution evaluator.

    @author: Pawel Salek <em>responsible</em>
*/

/* -*-mode:c; c-basic-offset:4; -*- */

#include "config.h"

#include <stdio.h>
#include <cmath>
#include <string.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>

#include <pthread.h>
static pthread_mutex_t dft_prop_mutex = PTHREAD_MUTEX_INITIALIZER;

#include "aos.h"
#include "dft_common.h"
#include "functionals.h"
#include "integrator.h"
#include "output.h"
#include "grid_matrix.h"
#include "rho-mat.h"
#include "utilities.h"

/* restrict hints should not be necessary... */
#if !defined(restrict)
#define restrict
#endif

inline int
min(int a, int b) {
  return a<b ? a : b;
}

typedef struct {
  const real *kappa;
  real *res;
  real* vt; /* dimensioned [bllen] for LDA, [bllen][4] for GGA */
  int   trplet, nbast, vecs_in_batch;
} LinRespBlData;

static void
lin_resp_cb_b_lda(DftIntegratorBl* grid, real * restrict tmp,
                  int bllen, int blstart, int blend,
                  LinRespBlData* data)
{
  const real * restrict aos = grid->atv;
  real * restrict excmat = data->res;
  real (* restrict vt) = data->vt; /* [bllen][4] */
  int ibl, i, jbl, j, k, isym, ivec;
  int nbast = data->nbast;
  int n2basx = nbast*nbast;
  FunDensProp dp = { 0 };

  for(ivec=0; ivec<data->vecs_in_batch; ivec++) {
    /* compute vector of transformed densities vt */
    getexp_blocked_lda(nbast, data->kappa + ivec*n2basx, grid->atv,
                       grid->bas_bl_cnt,
                       grid->basblocks, grid->shl_bl_cnt,
                       tmp, bllen, vt);

    for(i=blstart; i<blend; i++) {
      SecondDrv vxc;
      real weight = grid->weight[grid->curr_point+i];
      dp.rhoa = dp.rhob = 0.5*grid->r.rho[i];
      dftpot1_(&vxc, &weight, &dp, &data->trplet);
      vt[i] = vxc.fRR*vt[i]*2;
    }
    for(isym=0; isym<grid->nsym; isym++) {
      const real *vPot = vt;
      int (*restrict iblocks)[2] = BASBLOCK(grid,isym);
      int ibl_cnt = grid->bas_bl_cnt[isym];
            
      for(ibl=0; ibl<ibl_cnt; ibl++)
        for(i=iblocks[ibl][0]; i<iblocks[ibl][1]; i++) { 
          int ioff = i*bllen;
          for(k=blstart; k<blend; k++)
            tmp[k+ioff] = aos[k+ioff]*vPot[k];
        }

      for(ibl=0; ibl<ibl_cnt; ibl++) {
        for(i=iblocks[ibl][0]; i<iblocks[ibl][1]; i++) { 
          int ioff = i*nbast + ivec*n2basx;
          int jsym = 0; /* inforb_.muld2h[data->ksymop-1][isym]-1; */
          int (*restrict jblocks)[2] = BASBLOCK(grid,jsym);
          int jbl_cnt = grid->bas_bl_cnt[jsym];
          real *restrict tmpi = &tmp[i*bllen];
          if (isym<jsym) continue;
          for(jbl=0; jbl<jbl_cnt; jbl++) {
            int jtop = min(jblocks[jbl][1],i);
            for(j=jblocks[jbl][0]; j<jtop; j++) { 
              ergo_real s = 0;
              for(k=blstart; k<blend; k++)
                s += aos[k+j*bllen]*tmpi[k];
              excmat[j+ioff] += s;
            }
          }
          for(k=blstart; k<blend; k++)
            excmat[i+ioff] += aos[k+i*bllen]*tmpi[k]*0.5;
        }
      }
    }
  }
}

static void
lin_resp_cb_b_gga(DftIntegratorBl* grid, real * restrict tmp,
                  int bllen, int blstart, int blend,
                  LinRespBlData* data)
{
  int ibl, i, jbl, j, k, isym, ivec;
  int nbast = data->nbast;
  int n2basx = nbast*nbast;
  real (* restrict vt3)[4] = (real(*)[4])data->vt; /* [bllen][4] */
  real * restrict aos = grid->atv;
  real * restrict aox = grid->atv+bllen*nbast;
  real * restrict aoy = grid->atv+bllen*nbast*2;
  real * restrict aoz = grid->atv+bllen*nbast*3;
  real * restrict excmat = data->res;
  FunDensProp dp = { 0 };

  for(ivec=0; ivec<data->vecs_in_batch; ivec++) {
    /* compute vector of transformed densities and dens. gradients vt3 */
    getexp_blocked_gga(nbast, data->kappa + ivec*n2basx, grid->atv,
                       grid->bas_bl_cnt,
                       grid->basblocks, grid->shl_bl_cnt,
                       tmp, bllen, vt3);
    for(i=blstart; i<blend; i++) {
      SecondDrv vxc;
      real facr, facg;
      real weight = grid->weight[grid->curr_point+i];
      real ngrad  = template_blas_sqrt(grid->g.grad[i][0]*grid->g.grad[i][0]+
				       grid->g.grad[i][1]*grid->g.grad[i][1]+
				       grid->g.grad[i][2]*grid->g.grad[i][2]);
      real brg, brz, b0 = vt3[i][0];
      if(ngrad<1e-15|| grid->r.rho[i]<1e-15) {
        vt3[i][0] = vt3[i][1] = vt3[i][2] = vt3[i][3] = 0;
        continue;
      }
      brg = (vt3[i][1]*grid->g.grad[i][0] +
             vt3[i][2]*grid->g.grad[i][1] +
             vt3[i][3]*grid->g.grad[i][2]);
      brz = brg/ngrad;
      dp. rhoa = dp. rhob = 0.5*grid->r.rho[i];
      dp.grada = dp.gradb = 0.5*ngrad;
      dp.gradab = dp.grada*dp.gradb;
      dftpot1_(&vxc, &weight, &dp, &data->trplet);
      facr = vxc.fRZ*b0 + (vxc.fZZ-vxc.fZ/ngrad)*brz + vxc.fZG*brg;
      facr = facr/ngrad + (vxc.fRG*b0+vxc.fZG*brz +vxc.fGG*brg);
      facg = vxc.fZ/ngrad + vxc.fG;
      vt3[i][0] = vxc.fRR*b0 + vxc.fRZ*brz+ vxc.fRG*brg;
      vt3[i][1] = (grid->g.grad[i][0]*facr + facg*vt3[i][1])*2;
      vt3[i][2] = (grid->g.grad[i][1]*facr + facg*vt3[i][2])*2;
      vt3[i][3] = (grid->g.grad[i][2]*facr + facg*vt3[i][3])*2;
    }

    for(isym=0; isym<grid->nsym; isym++) {
      int (*restrict iblocks)[2] = BASBLOCK(grid,isym);
      int ibl_cnt = grid->bas_bl_cnt[isym];

      for(ibl=0; ibl<ibl_cnt; ibl++) {
        for(i=iblocks[ibl][0]; i<iblocks[ibl][1]; i++) { 
          real *restrict g0i = &aos[i*bllen];
          real *restrict gxi = &aox[i*bllen];
          real *restrict gyi = &aoy[i*bllen];
          real *restrict gzi = &aoz[i*bllen];
          int ioff = i*nbast + ivec*n2basx;
          int jsym = 0; /* inforb_.muld2h[data->ksymop-1][isym]-1; */
          int (*restrict jblocks)[2] = BASBLOCK(grid,jsym);
          int jbl_cnt = grid->bas_bl_cnt[jsym];
          for(k=blstart; k<blend; k++)
            tmp[k] = (vt3[k][0]*g0i[k] + 
                      vt3[k][1]*gxi[k] +
                      vt3[k][2]*gyi[k] +
                      vt3[k][3]*gzi[k]);
          for(jbl=0; jbl<jbl_cnt; jbl++) {
            int jtop = jblocks[jbl][1];
            for(j=jblocks[jbl][0]; j<jtop; j++) { 
              real *restrict g0j = &aos[j*bllen];
              real s = 0;
              for(k=blstart; k<blend; k++)
                s += g0j[k]*tmp[k];
              excmat[j+ioff] += s;
            }
          }
        }
      }
    }
  }
}

#if 0
static void
printmat(int n, const ergo_real *m, const char *name)
{
  int i, j;
  printf("Printing matrix %s\n", name);
  for(i=0; i<n; i++) {
    for(j=0; j<n; j++)
      printf("%10.5f", m[i + j*n]);
    puts("");
  }
}
#endif

/** dft_lin_respao performs the transformation of given transition
    density @param vec and the result is stored in @param trans_vec -
    both of which are square matrix A ground state density @param dens
    is required.
    @param bis is the basis set description structure.
    @param mol contains the molecule data (is this strictly needed?)
    @param gss a structure describing the grid settings.
    @param nThreads tells how many threads execute this section
    (needed for grid).
*/
EXTERN_C real
dft_lin_respao(const BasisInfoStruct& bis, const Molecule& mol,
	       const Dft::GridParams& gss,
               const real *dens, const real *vec, real* trans_vec,
               int nThreads)
{
    real electrons = 0;
    LinRespBlData lr_data; /* linear response data */
    int i, j, jvec;
    int nbast = bis.noOfBasisFuncs;
    int n2basx = nbast*nbast;
    Dft::FullMatrix density(dens, nbast);
    const Dft::FullMatrix *densPtr = &density;
    Util::TimeMeter tm;

    lr_data.vt    = dal_new(DFT_MAX_BLLEN*4, real);
    lr_data.kappa = vec;
    lr_data.res   = dal_new(nbast*nbast, real);
    lr_data.trplet= 0;  /* FIXME: should be 1 when finding excitations */
    lr_data.vecs_in_batch = 1;    /*nosim; */
    lr_data.nbast = nbast;
    memset(lr_data.res, 0, nbast*nbast*sizeof(real));
    electrons = Dft::integrate(1, &densPtr, bis, mol, gss, nThreads,
                               (DftBlockCallback)
                               (selected_func->is_gga() ?
                                lin_resp_cb_b_gga : lin_resp_cb_b_lda),
                               &lr_data);
    pthread_mutex_lock(&dft_prop_mutex);
    for(jvec=0; jvec<lr_data.vecs_in_batch; jvec++){
      for(i=0; i<nbast; i++) {
        int ioff = i*nbast + jvec*n2basx;
        for(j=0; j<i; j++) {
          int joff = j*nbast + jvec*n2basx;
          real averag = 0.5*(lr_data.res[i+joff] + lr_data.res[j+ioff]);
          trans_vec[i+joff] += averag;
          trans_vec[j+ioff] += averag;
        }
        trans_vec[i+ioff] += lr_data.res[i+ioff];
      }
    }
    pthread_mutex_unlock(&dft_prop_mutex);
    free(lr_data.res);
    free(lr_data.vt);

    if(nThreads<=1) {
        int nElectrons = mol.getNumberOfElectrons();
        do_output(LOG_CAT_INFO, LOG_AREA_LR, 
                  "LR-DFT*%d finished. Electrons: %f(%9.3g) (serial)", 
                  lr_data.vecs_in_batch, (double)electrons,
                  (double)((electrons-nElectrons)/nElectrons));
	tm.print(LOG_AREA_DFT, __func__);
    }
    return electrons;
}

struct LinData {
    const BasisInfoStruct *bis;
    const Molecule *mol;
    const Dft::GridParams *gss;
    const real     *density;
    const real     *inputVec;
    real           *transformedVec;
    real electrons;
    int nThreads;
};

static void*
dft_lin_resp_worker(void *data)
{
  static const int LINRESP_ERROR = 1;
    LinData *ld = (LinData*)data;
    try {
      ld->electrons =
        dft_lin_respao(*ld->bis, *ld->mol, *ld->gss,
                       ld->density, ld->inputVec, ld->transformedVec,
                       ld->nThreads);
    } catch(const char *s) {
      do_output(LOG_CAT_ERROR, LOG_AREA_DFT, 
		"dft_lin_resp_worker thread caught an exception '%s'", s);
      return (void*)&LINRESP_ERROR;
    } catch(const std::bad_alloc & e) {
      do_output(LOG_CAT_ERROR, LOG_AREA_DFT, 
		"dft_lin_resp_worker thread caught an exception '%s'", e.what());
      return (void*)&LINRESP_ERROR;
    }  catch(...) {
      do_output(LOG_CAT_ERROR, LOG_AREA_DFT, 
		"dft_lin_resp_worker thread caught unexpected exception.");
      return (void*)&LINRESP_ERROR;
    }
    return NULL;
}

EXTERN_C real
dft_lin_resp_mt(const BasisInfoStruct& bis, const Molecule& mol,
		const Dft::GridParams& gss,
                const real *dens, const real *vec, real* trans_vec)
{
    int i, threads;
    real electrons = 0;

    Util::TimeMeter tm;

    threads = dft_get_num_threads();
    std::vector<LinData> data(threads);
    std::vector<pthread_t> pids(threads);

    for(i=0; i<threads; i++) {
        data[i].bis = &bis;
        data[i].mol = &mol;
        data[i].gss = &gss;
        data[i].density = dens;
        data[i].inputVec = vec;
        data[i].transformedVec = trans_vec;
        data[i].nThreads = threads;
        if (pthread_create(&pids[i], NULL, dft_lin_resp_worker, &data[i])) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_DFT, 
		    "Creation of thread # %d failed\n", i);
	  if (i==0)
	    throw "No worker threads could be started";
	  else 
	    break;
	}
    }
    while (--i >= 0) {
        pthread_join(pids[i], NULL);
        electrons += data[i].electrons;
    }
    if(threads>1) {
        int nElectrons = mol.getNumberOfElectrons();
        do_output(LOG_CAT_INFO, LOG_AREA_LR, 
                  "LR-DFT*%d finished. Electrons: %f(%9.3g) (mt)", 
                  1, (double)electrons,
                  (double)((electrons-nElectrons)/nElectrons));
	tm.print(LOG_AREA_DFT, __func__);
    }
    return electrons;
}
