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
/** @file xc_evaluators.h defines LDA and GGA evaluators that work
    both for dense and sparse matrices. */

/** structure describing the data needed by distributors.  This is the
    restricted-orbital formalism, cross alpha-beta terms are not
    needed. */
template<typename Matrix>
struct KsData {
    Matrix *excmat;
    real* dR;
    real* dZ;
    real energy;
  KsData(Matrix *m_, size_t s)
    : excmat(m_), dR(new real[s]), dZ(new real[s]), energy(0)
  {}
  ~KsData() {
    delete []dR;
    delete []dZ;
  }
};


/** distributes a LDA-type xc potential over the XC-matrix elements,
    with optimization for a closed shell case. Matrix is a type that
    provides an operation add(row, col, val).
 */
template<typename Matrix>
struct XCDistributorLda {
  static void distribute(DftIntegratorBl *grid,
                         int bllen, int blstart, int blend,
                         real * restrict tmp, real *restrict dR,
                         Matrix& excmat)
  {
    static const int isym = 0;
    int jbl, j, ibl, i, k;
    real * restrict aos = grid->atv;

    int (*restrict blocks)[2] = BASBLOCK(grid,isym);
    int bl_cnt = grid->bas_bl_cnt[isym];

    for(jbl=0; jbl<bl_cnt; jbl++) {
        for(j=blocks[jbl][0]; j<blocks[jbl][1]; j++) { 
            int joff = j*bllen;
            for(k=blstart; k<blend; k++)
                tmp[k] = aos[k+joff]*dR[k]*0.5;

            for(ibl=0; ibl<bl_cnt; ibl++) {
                int top = blocks[ibl][1] < j
                    ? blocks[ibl][1] : j; 
                for(i=blocks[ibl][0]; i<top; i++) { 
                    real *restrict aosi = aos + i*bllen; /* ith orbital */
                    long_real s = 0;
                    for(k=blstart; k<blend; k++)
                        s += aosi[k]*tmp[k];
                    excmat.add(j,i, s);
                }
            }
            long_real s = 0;
            for(k=blstart; k<blend; k++)
                s += aos[k+j*bllen]*tmp[k];
            excmat.add(j,j, s);
        }
    }
  }
};

/** modifies data->excmat by adding LDA-type contributions from a
    given set of bllen grid points as saved in grid.
 */
template<typename Matrix, typename LDADistributor>
void
xcCallbackLdaR(DftIntegratorBl *grid, real * restrict tmp,
               int bllen, int blstart, int blend,
               KsData<Matrix>* data)
{
    FirstDrv drvs; 
    int k;
    real * restrict dR     = data->dR;
    FunDensProp dp = { 0 };
    
    assert(grid->ntypso >0);
    for(k=blstart; k<blend; k++) {
        real weight = grid->weight[grid->curr_point+k];
        dp.rhoa = dp. rhob = 0.5*grid->r.rho[k];
        data->energy += selected_func->func(&dp)*weight;
        dftpot0_(&drvs, &weight, &dp);
        dR[k] = 2*drvs.fR;
    }
    LDADistributor::distribute(grid, bllen, blstart, blend, tmp, dR,
                               *data->excmat);
}

/** distributes a GGA-type xc potential over the XC-matrix
    elements. Matrix is a type provides an operation add(row, col,
    val).
 */
template<typename Matrix>
struct XCDistributorGga {
  static void distribute(DftIntegratorBl *grid,
                         int bllen, int blstart, int blend,
                         real * restrict tmp,
                         const real *dR, const real *dZ,
                         Matrix& excmat)
  {
    static const int isym = 0;
    int jbl, j, ibl, i, k;
    const real * aox = grid->atv+bllen*grid->nbast;
    const real * aoy = grid->atv+bllen*grid->nbast*2;
    const real * aoz = grid->atv+bllen*grid->nbast*3;
    const real * aos = grid->atv;

    const int (*blocks)[2] = BASBLOCK(grid,isym);
    int nblocks = grid->bas_bl_cnt[isym];
    for(jbl=0; jbl<nblocks; jbl++)
        for(j=blocks[jbl][0]; j<blocks[jbl][1]; j++) { 
            int joff = j*bllen;
            for(k=blstart; k<blend; k++)
                tmp[k+joff] = 
                    (dR[k]* aos[k+j*bllen] +
                     dZ[k]*(aox[k+j*bllen]*grid->g.rad.a[k][0]+
                            aoy[k+j*bllen]*grid->g.rad.a[k][1]+
                            aoz[k+j*bllen]*grid->g.rad.a[k][2]));
        }
        
    for(jbl=0; jbl<nblocks; jbl++) {
        for(j=blocks[jbl][0]; j<blocks[jbl][1]; j++) { 
            const real * tmpj = tmp + j*bllen; /* jth orbital */
            for(ibl=0; ibl<nblocks; ibl++) {
                int top = blocks[ibl][1] < j
                    ? blocks[ibl][1] : j; 
                for(i=blocks[ibl][0]; i<top; i++) { 
                    long_real s = 0;
                    const real * tmpi = tmp + i*bllen; /* ith orbital */
                    for(k=blstart; k<blend; k++)
                        s += aos[k+i*bllen]*tmpj[k] +
                            aos[k+j*bllen]*tmpi[k];
                    excmat.add(i, j, s);
                }
            }
            long_real s = 0;
            for(k=blstart; k<blend; k++)
                s += 2*aos[k+j*bllen]*tmpj[k];
            excmat.add(j,j, s);
        }
    }
  }
};

/** modifies data->excmat by adding GGA-type contributions from a
    given set of bllen grid points as saved in grid.
 */
template<typename Matrix, typename GGADistributor>
void
xcCallbackGgaR(DftIntegratorBl* grid, real * restrict tmp, 
               int bllen, int blstart, int blend,
               KsData<Matrix>* data)
{
    FirstDrv drvs;
    int k;
    real * restrict dR = data->dR;
    real * restrict dZ = data->dZ;
    FunDensProp dp = { 0 };

    assert(grid->ntypso >0);
    for(k=blstart; k<blend; k++) {
        real weight = grid->weight[grid->curr_point+k];
        dp.grada = 0.5*template_blas_sqrt(grid->g.grad[k][0]*grid->g.grad[k][0]+
                                 grid->g.grad[k][1]*grid->g.grad[k][1]+
                                 grid->g.grad[k][2]*grid->g.grad[k][2]);
        dp. rhoa = dp.rhob = 0.5*grid->r.rho[k];
        dp.gradb  = dp.grada;
        dp.gradab = dp.grada*dp.gradb;
        if(dp.rhoa>1e-14) {
            if(dp.grada<1e-35) dp.grada = 1e-35;
            data->energy += selected_func->func(&dp)*weight;
            dftpot0_(&drvs, &weight, &dp);
            dR[k] = 0.5*drvs.fR;
            dZ[k] = 0.5*drvs.fZ/dp.grada;
        }else {
            dR[k] = dZ[k] = 0;
        }
    }

    GGADistributor::distribute(grid, bllen, blstart, blend, tmp, dR, dZ,
                               *data->excmat);
}

/* =================================================================== */
/* Unrestricted evaluators. */

template<typename Matrix>
struct UksData {
  Matrix *exca, *excb;
  real* dRa, *dRb;
  real* dZa, *dZb, *dZab;
  real energy;
  UksData(Matrix * a_, Matrix * b_, size_t s)
    : exca(a_), excb(b_),
      dRa(new real[s]), dRb(new real[s]),
      dZa(new real[s]), dZb(new real[s]), dZab(new real[s]),
      energy(0.0)
  {}
  ~UksData() {
    delete []dRa;    delete []dRb;
    delete []dZa;    delete []dZb;   delete []dZab;
  }
};

template<typename Matrix>
struct XCDistributorGgaU {
  static void
  distribute(DftIntegratorBl *grid, int bllen, int blstart, int blend,
             real * restrict tmp, const real * dR,
             const real *dZ1, const real (*grad1)[3],
             const real *dZ2, const real (*grad2)[3],
             Matrix& excmat);
};

template<typename Matrix>
void
XCDistributorGgaU<Matrix>::distribute(DftIntegratorBl *grid,
                                      int bllen, int blstart, int blend,
                                      real * restrict tmp, const real * dR,
                                      const real *dZ1,
                                      const real (*grad1)[3],
                                      const real *dZ2,
                                      const real (*grad2)[3],
                                      Matrix& excmat)
{
    int isym, jbl, j, ibl, i, k;
    real * restrict aox = grid->atv+bllen*grid->nbast;
    real * restrict aoy = grid->atv+bllen*grid->nbast*2;
    real * restrict aoz = grid->atv+bllen*grid->nbast*3;
    real * restrict aos = grid->atv;

    for(isym=0; isym<grid->nsym; isym++) {
        int (*restrict blocks)[2] = BASBLOCK(grid,isym);
        int nblocks = grid->bas_bl_cnt[isym];
        for(jbl=0; jbl<nblocks; jbl++)
            for(j=blocks[jbl][0]; j<blocks[jbl][1]; j++) { 
                int joff = j*bllen;
                for(k=blstart; k<blend; k++)
                    tmp[k+joff] = 
                        dR[k]* aos[k+j*bllen] +
                        dZ1[k]*(aox[k+j*bllen]*grad1[k][0]+
                                aoy[k+j*bllen]*grad1[k][1]+
                                aoz[k+j*bllen]*grad1[k][2])+
                        dZ2[k]*(aox[k+j*bllen]*grad2[k][0]+
                                aoy[k+j*bllen]*grad2[k][1]+
                                aoz[k+j*bllen]*grad2[k][2]);
        }
        
        for(jbl=0; jbl<nblocks; jbl++) {
            for(j=blocks[jbl][0]; j<blocks[jbl][1]; j++) { 
              real *restrict tmpj = tmp + j*bllen; /* jth orbital */
              for(ibl=0; ibl<nblocks; ibl++) {
//#define FULL_ONLY 1
#if defined(FULL_ONLY)
                for(i=blocks[ibl][0]; i<blocks[ibl][1]; i++) { 
                  long_real s = 0;
                  for(k=blstart; k<blend; k++)
                    s += aos[k+i*bllen]*tmpj[k];
                  excmat.add(i, j, s);
                }
#else /* FULL_ONLY */
                int top = blocks[ibl][1] < j
                  ? blocks[ibl][1] : j; 
                for(i=blocks[ibl][0]; i<top; i++) { 
                  long_real s = 0;
                  const real * tmpi = tmp + i*bllen; /* ith orbital */
                  for(k=blstart; k<blend; k++)
                    s += aos[k+i*bllen]*tmpj[k] +
                      aos[k+j*bllen]*tmpi[k];
                  excmat.add(i, j, s);                      
                }                    
#endif /* FULL_ONLY */
              }
#if !defined(FULL_ONLY)
              long_real s = 0;
              for(k=blstart; k<blend; k++)
                s += 2*aos[k+j*bllen]*tmpj[k];
              excmat.add(j,j, s);
#endif /* FULL_ONLY */
            }
        }
    }
}

/** modifies data->excmat by adding LDA-type xc contributions, for a
    general unrestricted case. Matrix is a type that provides an
    operation add(row, col, val). */
template<typename Matrix, typename LDADistributor>
void
xcCallbackLdaU(DftIntegratorBl *grid, real * restrict tmp,
               int bllen, int blstart, int blend,
               UksData<Matrix>* d)
{
    FunFirstFuncDrv drvs; 
    FunDensProp dp = { 0 };
    int k;

    assert(grid->ntypso >0);
    for(k=blstart; k<blend; k++) {
        real weight = grid->weight[grid->curr_point+k];
        dp.rhoa = grid->r.ho.a[k];
        dp.rhob = grid->r.ho.b[k];
        d->energy += selected_func->func(&dp)*weight;
        drv1_clear(&drvs);
        selected_func->first(&drvs, weight, &dp);
        d->dRa[k] = 2*drvs.df1000;
        d->dRb[k] = 2*drvs.df0100;
    }

    LDADistributor::distribute(grid, bllen, blstart, blend,
                               tmp, d->dRa, *d->exca);
    LDADistributor::distribute(grid, bllen, blstart, blend,
                               tmp, d->dRb, *d->excb);
}

/** modifes data->excmat by adding GGA-type xc contributions, for a
    general unrestricted case. Matrix is a type that provides an
    operation add(row, col, val). */
template<typename Matrix, typename GGADistributor>
static void
xcCallbackGgaU(DftIntegratorBl* grid, real * restrict tmp, 
               int bllen, int blstart, int blend,
               UksData<Matrix>* d)
{
    FunFirstFuncDrv drvs; 
    int k;
    FunDensProp dp = { 0 };

    assert(grid->ntypso >0);
    for(k=blstart; k<blend; k++) {
        const real THR = 1e-40;
        real weight = grid->weight[grid->curr_point+k];
        dp.rhoa = grid->r.ho.a[k];
        dp.rhob = grid->r.ho.b[k];
        dp.grada = template_blas_sqrt(grid->g.rad.a[k][0]*grid->g.rad.a[k][0]+
                             grid->g.rad.a[k][1]*grid->g.rad.a[k][1]+
                             grid->g.rad.a[k][2]*grid->g.rad.a[k][2]);
        dp.gradb = template_blas_sqrt(grid->g.rad.b[k][0]*grid->g.rad.b[k][0]+
                             grid->g.rad.b[k][1]*grid->g.rad.b[k][1]+
                             grid->g.rad.b[k][2]*grid->g.rad.b[k][2]);
        dp.gradab = grid->g.rad.a[k][0]*grid->g.rad.b[k][0]+
            grid->g.rad.a[k][1]*grid->g.rad.b[k][1]+
            grid->g.rad.a[k][2]*grid->g.rad.b[k][2];
        if(dp.rhoa+dp.rhob>1e-17) {
            if(dp.grada<THR) dp.grada = THR;
            if(dp.rhob<THR) dp.rhob = THR;
            if(dp.gradb<THR) dp.gradb = THR;
            if(dp.gradab<THR) dp.gradab = THR;
            d->energy += selected_func->func(&dp)*weight;
            drv1_clear(&drvs);
            selected_func->first(&drvs, weight, &dp);
            d->dRa[k]  = drvs.df1000*0.5;
            d->dRb[k]  = drvs.df0100*0.5;
            d->dZa[k]  = drvs.df0010/dp.grada;
            d->dZb[k]  = drvs.df0001/dp.gradb;
            d->dZab[k] = drvs.df00001;
        } else {
            d->dRa[k] = d->dRb[k] = d->dZa[k] = d->dZb[k] = d->dZab[k] = 0;
        }
    }

    GGADistributor::distribute(grid, bllen, blstart, blend, tmp, d->dRa,
                               d->dZa, grid->g.rad.a, d->dZab,
                               grid->g.rad.b, *d->exca);
    GGADistributor::distribute(grid, bllen, blstart, blend, tmp, d->dRb,
                               d->dZb, grid->g.rad.b, d->dZab,
                               grid->g.rad.a, *d->excb);
}
