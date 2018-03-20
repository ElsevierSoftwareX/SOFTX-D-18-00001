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

/** @file slr.cc contains a Simple Linear Response implementation
    based on the MO orbital picture. Copyright(c) Pawel Salek 2006.

    @author: Pawel Salek <em>responsible</em>

    The solver is generic: it can handle both linear sets of equations
    and eigenvalue problems. The variable parts are implemented by
    subclasses. A subclass must implement following functions:

    - getResidual(trial) - returns a residual vector being a
      difference between true solution and the solution in the
      subspace.

      In case of the set of equations, residual is computed in 3 steps:
      Solution in the subspace:          Xsub = (Asub-freq*Ssub)\\Ysub;
      Solution expanded out fo subspace: X=v*Xsub;
      Residual vector: residualv= (Av-freq*Sv)*Xsub - Y;
      
      In case of the eigenvalue problem, more steps are needed:
      Solution in the subspace: [ Xsub, lambda ] = eig(Asub, Ssub);
      Pick first positive eigenvalue l1 = lambda(step+1);
      Pick corresponding eigenvector: Xsub = Xsub(:,step+1);
      Residual Vector:  residualv = (Av-l1*Sv)*Xsub;

    - getPreconditonerShift() - get the preconditioner shift. In case
      of the set of equations, it is shifted by the constant
      frequency. In case of set of eigenvalues, it is shifted by the
      best approximation to the required solution, obtained in the
      first step from the difference of KS matrix eigenvalues or taken
      as the Ritz value when they are available.

    - all other stuff, like vector transformation and subspace
      extension are handled by the generic solver.

  Example usage:
  <pre>
  MyE2Evaluator e2;
  EigenSolver solver(nbast, nocc, fock_matrix, overlap_matrix);
  solver.solve(e2);
  printf("Lowest eigenvalue: %f\n", solver.getFreq());
  </pre>
  The important stuff will get printed but also a solution should
  probably be returned in some convenient way.
*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <vector>

#include "mat_gblas.h"
#include "slr.h"
#include "solve_lin_eq_syst.h"
#include "utilities.h"
#include "output.h"
#include "units.h"

namespace LR {
/** returns dot_product(a, b) */
inline ergo_real
dot(int len, const ergo_real* a, const ergo_real *b)
{
  ergo_real r=0;
  for(int i=0; i<len; i++) r += a[i]*b[i];
  return r;
}


inline ergo_real operator*(const VarVector& a, const VarVector& b){
  assert(a.nvar == b.nvar);
  assert(a.inMemory); /* or load()? */
  assert(b.inMemory); /* or load()? */
  return dot(a.nvar*2, a.v, b.v);
}


inline ergo_real
operator*(const VarVector &a, const VarVectorProxyOp<false,false>& b){
  assert(a.inMemory);     /* or load()? */
  assert(b.vec.inMemory); /* or load()? */

  return b.scalar*dot(a.nvar*2, a.v, b.vec.v);
}

inline ergo_real
operator*(const VarVector &a, const VarVectorProxyOp<true,false>& b){
  assert(a.inMemory);
  assert(b.vec.inMemory);

  return 2*b.scalar*(dot(a.nvar,a.v, b.vec.v)-
                     dot(a.nvar,a.v+a.nvar, b.vec.v+a.nvar));
}

inline ergo_real
operator*(const VarVector &a, const VarVectorProxyOp<false,true>& b){
  ergo_real r=0;
  assert(a.inMemory);
  assert(b.vec.inMemory);

  for(int i=0; i<a.nvar; i++) {
    r += b.scalar * a[i]        * b.vec[i+a.nvar];
    r += b.scalar * a[i+a.nvar] * b.vec[i];
  }
  return r;
}


inline VarVectorProxyOp<false, false>
operator*(ergo_real s, const VarVector& v) {
  return VarVectorProxyOp<false,false>(v, s);
}

template<bool MultByS, bool SwapXY>
inline VarVectorProxyOp<MultByS, SwapXY>
operator*(ergo_real s, const VarVectorProxyOp<MultByS, SwapXY>& v) {
  return VarVectorProxyOp<MultByS,SwapXY>(v.vec, s*v.scalar);
}

inline VarVector&
operator+=(VarVector& a, const VarVectorProxyOp<false,false>& proxy){
  assert(a.inMemory);
  assert(proxy.vec.inMemory);
  for(int i=0; i<2*a.nvar; i++)
    a[i] += proxy.scalar*proxy.vec[i];
  return a;
}

inline VarVector&
operator+=(VarVector& a, const VarVectorProxyOp<true,false>& proxy){
  assert(a.nvar == proxy.vec.nvar);
  assert(a.inMemory);
  assert(proxy.vec.inMemory);

  for(int i=0; i<a.nvar; i++) {
    a[i]        +=  2*proxy.scalar*proxy.vec[i];
    a[i+a.nvar] += -2*proxy.scalar*proxy.vec[i+a.nvar];
  }
  return a;
}

inline VarVector&
operator+=(VarVector& a, const VarVectorProxyOp<false,true>& proxy){
  assert(a.nvar == proxy.vec.nvar);
  assert(a.inMemory);
  assert(proxy.vec.inMemory);

  for(int i=0; i<a.nvar; i++) {
    a[i]        += proxy.scalar*proxy.vec[i+a.nvar];
    a[i+a.nvar] += proxy.scalar*proxy.vec[i];
  }
  return a;
}

inline VarVector&
operator+=(VarVector& a, const VarVectorProxyOp<true,true>& proxy){
  assert(a.nvar == proxy.vec.nvar);
  assert(a.inMemory);
  assert(proxy.vec.inMemory);
  for(int i=0; i<a.nvar; i++) {
    a[i]        +=  2*proxy.scalar*proxy.vec[i+a.nvar];
    a[i+a.nvar] += -2*proxy.scalar*proxy.vec[i];
  }
  return a;
}

/** returns a proxy object corresponding to a swapped vector. */
const VarVectorProxyOp<false,true>
swapXY(const VarVector& arg)
{ return VarVectorProxyOp<false, true >(arg, 1.0) ; }


/** returns a proxy object corresponding to a vector multiplied by
    S[2], i.e. v -> S[2]*v. */
const VarVectorProxyOp<true,false>
sTimes(const VarVector& arg)
{ return VarVectorProxyOp<true,false>(arg, 1.0) ; }

template<bool SwapXY>
VarVectorProxyOp<true, SwapXY>
sTimes(const VarVectorProxyOp<false, SwapXY>& arg)
{ return VarVectorProxyOp<true, SwapXY >(arg.vec, arg.scalar) ; }



/* VarVector implementation. */
void
VarVector::load(const char* tmpdir)
{
  if(inMemory)
    return;
  if (nvar == 0) {
    inMemory = 1;
    return;
  }

  do_output(LOG_CAT_INFO, LOG_AREA_LR, "load::'%s' with nvar=%d", fName, nvar);
  if (!fName) throw "loading not saved VarVector";

  int fd = open(fName, O_RDONLY, 0);
  if (fd == -1)
    throw "VarVector disappeared from disk";
  if (!v)
    v = new ergo_real[2*nvar];

  ssize_t readAlready = 0, toRead = nvar*2*sizeof(ergo_real);
  do {
    ssize_t ret = read(fd, v + readAlready, toRead);
    if (ret != -1) {
      readAlready += ret;
      toRead      -= ret;
    }
  } while (toRead>0);
  close(fd);
  inMemory = 1;
}

void
VarVector::save(const char* tmpdir)
{
  if(onDisk || nvar == 0)
    return;

  if (!fName) {
    /* NOTE: earlier fName was allocated as new char[strlen(tmpdir) +
       4 + 8 + 1] which was too small on some systems (e.g. luc2)
       where pointer strings are 16 chars long. This gave random
       segfaults and other mysterious program crashes. Elias changed
       the allocation to new char[888] 2008-12-02. */
    fName = new char[888];
    /* FIXME: Consider different creation of temporary file names. */
    sprintf(fName, "%s/LR_%-8p", tmpdir, (void*)this);
  }
  static mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
  int fd = open(fName, O_WRONLY | O_CREAT | O_TRUNC, mode);
  if (fd == -1)
    throw "Cannot save VarVector";

  ssize_t written = 0, toWrite = nvar*2*sizeof(ergo_real);
  do_output(LOG_CAT_INFO, LOG_AREA_LR, "save:'%s' writes %u bytes nvar=%d",
            fName, toWrite, nvar);

  do {
    ssize_t ret = write(fd, v + written, toWrite);
    if (ret != -1) { /* Handle interrrupted system calls... */
      written += ret;
      toWrite -= ret;
    }
  } while (toWrite>0);
  close(fd);
  onDisk = 1;
}

void
VarVector::release(const char* tmpdir)
{
  if(nvar == 0)
    return;
  if(!onDisk)
    save(tmpdir);
  delete []v;
  v = NULL;
  inMemory = 0;
}



/* VarVectorCollection implementation. */
const char *VarVectorCollection::tmpdir = "/tmp";

VarVectorCollection::~VarVectorCollection()
{
  if(vecs) delete []vecs;
  if(ages) delete []ages;
}

void
VarVectorCollection::setSize(int sz)
{
  if(sz > nAllocated) {
    delete []vecs; nAllocated = sz; vecs = new VarVector[sz];
    delete []ages; ages = new unsigned[sz];
  }
  nVecs = sz;
}

VarVector&
VarVectorCollection::operator[](int i)
{
  if (!vecs[i].inMemory) {
    release();
    vecs[i].load(tmpdir);
  }
  ages[i] = currentAge++;
  return vecs[i];
}

void
VarVectorCollection::release()
{
  if(!diskMode)
    return; /* Nothing needs to be done. */

  /** Must allow at least two vectors at the same time in memory or
   *  evil things will happen. */
  int oldestIdx = nVecs, used=0;

  for(int i=1; i<nVecs; i++) {
    if (vecs[i].nvar && vecs[i].inMemory ) {
      used++;
      if (oldestIdx == nVecs || ages[oldestIdx]>ages[i])
        oldestIdx = i;
    }
  }

  if( oldestIdx < nVecs && used>2) {
    do_output(LOG_CAT_INFO, LOG_AREA_LR,
              "ONE: releases vector no. %d %s",
              oldestIdx, vecs[oldestIdx].fName);
    vecs[oldestIdx].release(tmpdir);
  }
}

void
VarVectorCollection::releaseAll()
{
  if(!diskMode)
    return; /* Nothing needs to be done. */
  if (diskMode) {
    do_output(LOG_CAT_INFO, LOG_AREA_LR, "releaseAll called.");
  }
  for(int i=0; i<nVecs; i++) {
    if (vecs[i].nvar > 0 && vecs[i].inMemory ) {
      do_output(LOG_CAT_INFO, LOG_AREA_LR, "ALL: releases vector %d with %d", i, i, vecs[i].nvar);
      
      vecs[i].release(tmpdir);
    }
  }
}

/** increase the dimension of the matrix without losing the data. */
void
SmallMatrix::expand(int newSize)
{
  if(newSize > nsize) {
    ergo_real *m = new ergo_real[newSize*newSize];
    for(int i=0; i<nsize; i++)
      for(int j=0; j<nsize; j++)
        m[j + i*newSize] = mat[j + i*nsize];
    delete []mat;
    mat = m;
    nsize = newSize;
  }
}

/** Initialize the solver by computing the diagonal of the E2 operator
    as needed for preconditioning. */
LRSolver::LRSolver(int aNbast, int aNocc,
                   const ergo_real *fock_matrix,
                   const ergo_real *s)
  : convThreshold(1e-3), maxSubspaceSize(MVEC),
    eSub(MVEC), sSub(MVEC),
    xSub(new ergo_real[MVEC]),
    nbast(aNbast), nocc(aNocc),
    vects(MVEC), Avects(MVEC),
    fdiag(NULL), cmo(NULL)
{
  if(fock_matrix) setE2diag(nbast, nocc, fock_matrix, s);
}


void
VarVector::setFromFull(int nbast, int nocc, const ergo_real *full_mat)
{
  int nvirt = nbast - nocc;
  nvar = nocc*nvirt;
  inMemory = 1;
  onDisk = 0;

  if(nvar == 0) throw "No variables - no excitations";
  if(v) delete []v;
  v = new ergo_real[2*nvar];

  if(!full_mat)
    return;

  int idx =0;
  for(int col=0; col<nocc; col++) {
    for(int row=0; row<nvirt; row++)
      v[idx + row] = full_mat[nocc+row + col*nbast];
    idx += nvirt;
  }

  for(int row=0; row<nocc; row++) {
    for(int col=0; col<nvirt; col++)
      v[idx + col] = full_mat[row + (nocc + col)*nbast];
    idx += nvirt;
  }
}

void
VarVector::setFull(int nbast, int nocc, ergo_real *full_mat) const
{
  int i, j, idx;
  int nvirt = nbast-nocc;

  idx = 0;
  /* first occ columns */
  for(j=0; j<nocc; j++) {
    for(i=0; i<nocc; i++)  full_mat[i + j*nbast] = 0;
    for(i=0; i<nvirt; i++) full_mat[nocc+i + j*nbast] = v[i + idx];
    idx += nvirt;
  }
  /* The "other" block */
  for(i=0; i<nocc; i++) {
    for(j=0; j<nvirt; j++)
      full_mat[i + (nocc+j)*nbast] = v[j + idx];
    idx += nvirt;
  }
  /* Zero the remaining part */
  for(j=nocc; j<nbast; j++)
    for(i=nocc; i<nbast; i++) full_mat[i + j*nbast] = 0;
}

/** Uses symmetric orthogonalization to orthonormalize itself (x y)
    with the swapped vector (y x). It is achieved by performing a
    following transformation: av = [ x y; y x]; s = av'*av [v, e] =
    eig(s); e=diag(e)' a = av* v*diag(e.^-0.5)*v'.

    It may happen that X=Y (ovl=0) - one such case is running Hartree
    approximation. In that case we set Y=0.
*/

/**  */
void VarVector::symorth(void)
{
  ergo_real nrm= dot(nvar*2, v, v);
  ergo_real ovl= 2*dot(nvar, v, v+nvar )/nrm;
  ergo_real x1 = 1+ovl;
  ergo_real x2 = 1-ovl;
  if(template_blas_fabs(x1)>1e-10 && template_blas_fabs(x2)>1e-10) {
    x1 = 0.5 / template_blas_sqrt(x1*nrm);
    x2 = 0.5 / template_blas_sqrt(x2*nrm);
    ergo_real c1 = (x1 + x2);
    ergo_real c2 = (x1 - x2);
    for(int i=0; i<nvar; i++) {
      ergo_real xi = c1*v[i] + c2*v[i+nvar];
      ergo_real yi = c2*v[i] + c1*v[i+nvar];
      v[i]      = xi;
      v[i+nvar] = yi;
    }
  } else {
    do_output(LOG_CAT_INFO, LOG_AREA_LR, "Removing Y part, normfactor = %g", (double)nrm);
    ergo_real f = template_blas_sqrt(2.0/nrm);
    for(int i=0; i<nvar; i++) {
      v[i]      = v[i]*f;
      v[i+nvar] = 0;
    }
  }

  onDisk = 0;
}

/** pre-condition a vector given an approximation of the E[2] operator
    diagonal and a shift of the S[2] operator. */
static void
precondition(VarVector&v, const VarVector& e2diag, ergo_real shift)
{
  ergo_real denom;
  /* Factors 2 here is shaky - compare against exact E2 diagonal. */
  for(int i=0; i<v.nvar; i++) {
    denom = e2diag.x()[i] - shift;
    if(template_blas_fabs(denom)>0) v.x()[i] /= denom;
    denom = e2diag.y()[i] + shift;
    if(template_blas_fabs(denom)>0) v.y()[i] /= denom;
  }
  v.onDisk = 0;
}

/** mat := [mat, D_MO] */
static void
commuteWithDMO(int nbast, int nocc, ergo_real *mat)
{
  int col, row;

  for(col=0; col<nocc; col++) {
    for(row=0;    row<nocc; row++)  mat[row + col*nbast]  = 0;
    for(row=nocc; row<nbast; row++) mat[row + col*nbast] *= 2;
  }
  for(col=nocc; col<nbast; col++) {
    for(row=0;    row<nocc; row++)  mat[row + col*nbast] *= -2;
    for(row=nocc; row<nbast; row++) mat[row + col*nbast]  = 0;
  }
}

#if 0
static void
printmat(int n, const ergo_real *m, const char *name)
{
  printf("Printing matrix %s\n", name);

  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++)
      printf("%10.5f", m[i + j*n]);
    puts("");
  }
}

/** Compute c := beta*c + alpha*op(a)*op(b) */
static void
dgemm(int n, const char *at, const ergo_real *a,
      const char *bt, const ergo_real *b, 
      ergo_real alpha, ergo_real beta, ergo_real *c)
{
  int i, j, k;
  for(j=0; j<n; j++)
    for(i=0; i<n; i++) {
      /* c may be uninitialized here if beta == 0 */
      if(beta == 0.0) c[i + j*n] = 0;
      else c[i + j*n] *= beta;
      switch(*at) {
      case 'N':
      case 'n':
        switch(*bt) {
        case 'N':
        case 'n':
          for(k=0; k<n; k++) c[i + j*n] += alpha*a[i + k*n]*b[k + j*n];
          break;
        case 'T':
        case 't':
          for(k=0; k<n; k++) c[i + j*n] += alpha*a[i + k*n]*b[j + k*n];
            break;
        }
        break;
      case 'T':
      case 't':
        switch(*bt) {
        case 'N':
        case 'n':
          for(k=0; k<n; k++) c[i + j*n] += alpha*a[k + i*n]*b[k + j*n];
            break;
        case 'T':
        case 't':
          for(k=0; k<n; k++) c[i + j*n] += alpha*a[k + i*n]*b[j + k*n];
            break;
        }
      }
    }
}
#endif

static inline void
gemm(int n, const char *at, const ergo_real *a,
      const char *bt, const ergo_real *b, 
      ergo_real alpha, ergo_real beta, ergo_real *c)
{
  mat::gemm(at, bt, &n, &n, &n,
            &alpha, a, &n,
            b, &n, &beta, c, &n);
}

/** res := f*res + [a, b] */
static void
commute(int nbast, const ergo_real *a, const ergo_real *b,
        ergo_real f, ergo_real *res)
{
  gemm(nbast, "N", a, "N", b,  1, f, res);
  gemm(nbast, "N", b, "N", a, -1, 1, res);
}

/* computes ao := cmo*mo*cmo' */
void
LRSolver::mo2ao(int nbast, const ergo_real *mo, ergo_real *ao) const
{
  assert(cmo);
  ergo_real *tmp = new ergo_real[nbast*nbast];
  gemm(nbast, "N", cmo, "N", mo,  1, 0, tmp);
  gemm(nbast, "N", tmp, "T", cmo, 1, 0, ao);
  delete []tmp;
}

/** computes mo := cmo'*ao*cmo */
void
LRSolver::ao2mo(int nbast, const ergo_real *ao, ergo_real *mo) const
{
  assert(cmo);
  ergo_real *tmp = new ergo_real[nbast*nbast];
  gemm(nbast, "T", cmo, "N", ao,  1, 0, tmp);
  gemm(nbast, "N", tmp, "N", cmo, 1, 0, mo);
  delete []tmp;
}

void
LRSolver::load_F_MO(ergo_real *fmat) const
{
  int row, col;

  for(col=0; col<nbast; col++) {
    for(row=0; row<col; row++) fmat[row + col*nbast] = 0;
    fmat[col + col*nbast] = fdiag[col];
    for(row=col+1; row<nbast; row++) fmat[row + col*nbast] = 0;
  }
}

/** extends the subspace with v and its transformed vector Av. The
    eSub and sSub projected matrices are modified as well. We do not
    store explicitely the swapped vectors (y x), only the (x y)
    combination. */
void
LRSolver::addToSpace(VarVectorCollection& v, E2Evaluator& e2)
{
  VarVector Av;
  if(subspaceSize + 2*v.size()>maxSubspaceSize)
    throw "addToSpace: subspace size exceeded";

  int i, j;
  int orig_size = subspaceSize;

  for(int ivec=0; ivec<v.size(); ivec++) {
    precondition(v[ivec], e2diag, getPreconditionerShift(ivec));
    /* orthogonalize */
    for(i=0; i<subspaceSize; i+=2) {
      ergo_real prod = v[ivec] * vects[i];
      if(template_blas_fabs(prod)>1e-10) v[ivec] += -prod*vects[i];
      prod = v[ivec] * swapXY(vects[i]);
      if(template_blas_fabs(prod)>1e-10) v[ivec] += -prod*swapXY(vects[i]);
    }
    v[ivec].symorth();
    v.releaseAll();
    vects.releaseAll();
    if(!lintrans(e2, v[ivec], Av)) throw "E[2]*x failed";
  
    vects[subspaceSize]    = v[ivec];
    Avects[subspaceSize]   = Av;
    //vects[subspaceSize+1]  = swapXY(v[ivec]);
    //Avects[subspaceSize+1] = swapXY(Av);
    subspaceSize += 2;
  }

  /* compute missing elements of the space */
  for(i=orig_size; i<subspaceSize; i+=2)
    for(j=0; j<subspaceSize; j+=2) {
#if 1
      ergo_real eiXjX = dot(Avects[i].nvar, Avects[i].x(), vects[j].x());
      ergo_real eiXjY = dot(Avects[i].nvar, Avects[i].x(), vects[j].y());
      ergo_real eiYjX = dot(Avects[i].nvar, Avects[i].y(), vects[j].x());
      ergo_real eiYjY = dot(Avects[i].nvar, Avects[i].y(), vects[j].y());
      
      ergo_real siXjX = 2*dot(vects[i].nvar, vects[i].x(), vects[j].x());
      ergo_real siXjY = 2*dot(vects[i].nvar, vects[i].x(), vects[j].y());
      ergo_real siYjX = 2*dot(vects[i].nvar, vects[i].y(), vects[j].x());
      ergo_real siYjY = 2*dot(vects[i].nvar, vects[i].y(), vects[j].y());
      
      eSub[i][j]     = eSub[j][i]     = eiXjX + eiYjY;
      eSub[i][j+1]   = eSub[j+1][i]   = eiXjY + eiYjX;
      eSub[i+1][j]   = eSub[j][i+1]   = eiYjX + eiXjY;
      eSub[i+1][j+1] = eSub[j+1][i+1] = eiYjY + eiXjX;

      sSub[i][j]     = sSub[j][i]     = siXjX - siYjY;
      sSub[i][j+1]   = sSub[j+1][i]   = siXjY - siYjX;
      sSub[i+1][j]   = sSub[j][i+1]   = siYjX - siXjY;
      sSub[i+1][j+1] = sSub[j+1][i+1] = siYjY - siXjX;
#else
      ergo_real eij = Avects[i] * vects[j];
      eSub[i][j] = eSub[j][i] = eij;
      ergo_real sij =  vects[i] * sTimes(vects[j]);
      sSub[i][j] = sSub[j][i] = sij;
#endif
    }
}


/** performs the linear transformation of the vector with E[2]
    operator. */
bool
LRSolver::lintrans(E2Evaluator& e2, const VarVector& v, VarVector& Av) const
{
  bool res;
  ergo_real *transDens = new ergo_real[nbast*nbast];

  v.setFull(nbast, nocc, transDens);
  commuteWithDMO(nbast,nocc,transDens); /* transDens := [transDens, D_MO] */
  
  ergo_real *dmat = new ergo_real[nbast*nbast];
  ergo_real *fmat = new ergo_real[nbast*nbast];
  mo2ao(nbast, transDens, dmat);      /* dmat := cmo*tranDens*cmo' */

  res = e2.transform(dmat, fmat);
  ergo_real *f_MO = dmat; dmat = NULL; /* steal another pointer */

  ao2mo(nbast, fmat, f_MO);
  commuteWithDMO(nbast,nocc,f_MO); /* one contribution is ready... */

  load_F_MO(fmat);
  commute(nbast, fmat, transDens, 1, f_MO); /* f_MO += [fmat, transDens] */
  Av.setFromFull(nbast, nocc, f_MO);
  delete []fmat;
  delete []dmat;
  delete []transDens;
  return res;
}

void
LRSolver::increaseSubspaceLimit(int newSize)
{
  if(newSize > maxSubspaceSize) {
    eSub.expand(newSize);
    sSub.expand(newSize);
    ergo_real *v = new ergo_real[newSize];
    for(int i=0; i<maxSubspaceSize; i++) v[i] = xSub[i];
    delete []xSub; xSub = v;
    vects.setSize(newSize);
    Avects.setSize(newSize);
    maxSubspaceSize = newSize;
  }
}

/** solve the problem as defined by the subclass. This involves
    generation of the initial guess, symmetric orthogonalization,
    subspace extension routines, etc. */
bool
LRSolver::solve(E2Evaluator& e2, bool diskMode)
{
  VarVector Av;
  VarVectorCollection v;
  Util::TimeMeter tm;

  vects.setDiskMode(diskMode);
  v.setDiskMode(diskMode);

  getInitialGuess(v);
  do_output(LOG_CAT_INFO, LOG_AREA_LR,
            "LRSolver::solve entered with Nvariables: %d DiskMode: %s",
            v[0].nvar*2, diskMode ? "YES" : "NO");
  //computeExactE2Diag(e2); return false;
  subspaceSize = 0;
  for(int step = 0; step < maxSubspaceSize/2; step++) {
    /* Subspace extension phase */
    addToSpace(v, e2); /* Add the new vectors and precompute the
                        * projected operators as well. */
    /* residue estimation phase */
    if(!getResidual(v))
      break;
    ergo_real nr = template_blas_sqrt(v[0] * v[0]);

    printf("Step %2d: ssize: %3d |r|=%13.6g shift=%10.6f\n",
           step+1, subspaceSize, (double)nr, (double)getPreconditionerShift(0));
    do_output(LOG_CAT_INFO, LOG_AREA_LR,
              "Step %2d: ssize: %3d |r|=%13.6g shift=%10.6f",
	      step+1, subspaceSize, (double)nr, (double)getPreconditionerShift(0));
  }
  tm.print(LOG_AREA_LR, "LRSolver::solve");
  return true;
}

void
LRSolver::computeExactE2Diag(E2Evaluator& e2)
{
  VarVector v, Av;
  int nvar = (nbast-nocc)*nocc;
  v.setSize(nvar); Av.setSize(nvar);

  for(int i=0; i<nvar*2; i++) {
    v = 0; v[i] = 1;
    if(!lintrans(e2, v, Av))
      throw "Linear Transformation with E[2] failed";
    printf("%4i: %20.6g %20.6g\n", i, (double)e2diag[i], (double)Av[i]);

  }
}

ergo_real
LRSolver::setE2diag(int nbast, int nocc,
                    const ergo_real *fock_matrix,
                    const ergo_real *s)
{
  static const int ITYPE = 1;
  assert(fdiag == NULL); /* prevent double-initialization */

  int i, j;
  int lwork = 18*nbast;
  std::vector<ergo_real> work(lwork);
  std::vector<ergo_real> ovl(nbast*nbast);
  ergo_real ret;
  fdiag = new ergo_real[nbast];
  cmo   = new ergo_real[nbast*nbast];

  if(nocc<=0) throw "setE2diag::At least one orbital must be occupied";
  for(i=nbast*nbast-1; i>=0; i--) {
    cmo[i] = fock_matrix[i];
    ovl[i] = s[i];
  }
  //printmat(nbast, cmo, "fock");
  //printmat(nbast, ovl, "ovl");
  mat::sygv(&ITYPE, "V", "L", &nbast, cmo, &nbast, &ovl[0],
            &nbast, fdiag, &work[0], &lwork, &i);
  if(i == 0) {
    std::vector<ergo_real> ediff(nbast*nbast);
    for(i=0; i<nbast; i++) {
      for(j=0; j<i; j++)       ediff[j + nbast*i] = -(fdiag[j]-fdiag[i])*2.0;
      ediff[i + nbast*i] = 0; /* not used */
      for(j=i+1; j<nbast; j++) ediff[j + nbast*i] = +(fdiag[j]-fdiag[i])*2.0;
    }
    e2diag.setFromFull(nbast, nocc, &ediff[0]);
    ret = fdiag[nocc]-fdiag[nocc-1];
    
    do_output(LOG_CAT_INFO, LOG_AREA_LR, "HOMO-LUMO gap: %20.15f eV",
              (double)ret / UNIT_one_eV);

    return ret;
  } else {
    do_output(LOG_CAT_ERROR, LOG_AREA_LR, "setE2diag: mat::sygv failed with info=%d", i);
    throw "setE2diag::dsygv failed";
  }
}

/** get_av_minus_freq_sv scans through transformed vectors creating
    their linear combination and returning r := Av - f*Sv
*/
void
LRSolver::getAvMinusFreqSv(ergo_real f, ergo_real *weights,
                           VarVector& r)
{
  int i;

  r  = weights[0]*Avects[0];
  r += weights[1]*swapXY(Avects[0]);
  for(i=2; i<subspaceSize; i+=2) {
    r += weights[i]  *Avects[i];
    r += weights[i+1]*swapXY(Avects[i]);
  }

  if(f != 0) {
    for(i=0; i<subspaceSize; i+=2) {
      ergo_real w = -f*weights[i];
      r += w*sTimes(vects[i]);
      w = -f*weights[i+1];
      r += w*sTimes(swapXY(vects[i]));
    }
  }
}

/** Projects a full vector onto the reduced space. */
void
LRSolver::projectOnSubspace(const VarVector& full, ergo_real *w)
{
  for(int i=0; i<subspaceSize; i+=2) {
    ergo_real vXrX = dot(full.nvar, vects[i].x(), full.x());
    ergo_real vXrY = dot(full.nvar, vects[i].x(), full.y());
    ergo_real vYrX = dot(full.nvar, vects[i].y(), full.x());
    ergo_real vYrY = dot(full.nvar, vects[i].y(), full.y());
    w[i]   = vXrX + vYrY;
    w[i+1] = vYrX + vXrY;
  }
}

void
LRSolver::buildVector(const ergo_real *w, VarVector& full) /* const */
{
  full.setSize(vects[0].nvar);
  full = 0;
  for(int i=0; i<subspaceSize; i+=2) {
    for(int j=0; j<full.nvar; j++) {
      full[j]           += w[i]*vects[i][j] + w[i+1]*vects[i][j + full.nvar];
      full[j+full.nvar] += w[i]*vects[i][j+full.nvar] + w[i+1]*vects[i][j];
    }
  }
}

void
LRSolver::operToVec(OneElOperator& oper, VarVector& res) const
{
  ergo_real *buf = new ergo_real[nbast*nbast];
  ergo_real *tmp = new ergo_real[nbast*nbast];
  oper.getOper(buf);
  ao2mo(nbast, buf, tmp);
  delete []buf;
  commuteWithDMO(nbast, nocc, tmp);

  res.setFromFull(nbast, nocc, tmp);
  delete []tmp;
}

/* ===================================================================
 *              Solver of Linear Set of Equations
 * =================================================================== */

/** returns the initial guess for the linear set of equations. The
    explicit value is obtained from the diagonal assumption for the
    E[2] operator and is:
    (E[2]- freq*S[2])*g = Y ->
    g = Y./(E[2]-freq*S[2])
*/
int
SetOfEqSolver::getInitialGuess(VarVectorCollection& guess)
{
  int nvar = e2diag.nvar;
  if(rhs.nvar == 0)
    throw "SetOfEqSolver::getInitialGuess() called without RHS set";
  guess.setSize(1);
  VarVector& v = guess[0];
  v.setSize(nvar);
  /* In principle, the best approximation to the solution we can get
     is RHS divided by the approximation of the diagonal but the last
     operation is done when preconditioning which is the first op that
     we do. That's why we set the guess to the right side with
     switched sign for the lower part only. */
  for(int i=0; i<nvar; i++) {
    v[i     ] =  0.5*rhs[i];
    v[i+nvar] = -0.5*rhs[i+nvar];
  }
  return 0;
}

/** multiplies current solution by some vector. If such contractions
    are to be done several times, perhaps a single vector solution
    should be created and only then contracted with rhs vector. */
ergo_real
SetOfEqSolver::multiplyXtimesVec(const VarVector& rhs) /* const */
{
  ergo_real *proj = new ergo_real[subspaceSize];
  projectOnSubspace(rhs, proj);
  ergo_real res = 0;
  for(int i=0; i<subspaceSize; i++)
    res += xSub[i] * proj[i];
  delete []proj;
  return res;
}


/** get the residual of the set of linear equations. This is done in
    two steps: Solution in the subspace: Xsub = (eSub-freq*Ssub)\\Ysub;
    Residual vector is: residualv= (Av-freq*Sv)*Xsub - Y;
*/
bool
SetOfEqSolver::getResidual(VarVectorCollection& residualv)
{
  if(rhs.nvar == 0)
    throw "SetOfEqSolver::getResidual() called without RHS set";

  ergo_real *A = new ergo_real[subspaceSize*subspaceSize];

  for(int j=0; j<subspaceSize; j++) {
    for(int i=0; i<subspaceSize; i++)
      A[i+j*subspaceSize] =
        eSub[j][i] - frequency*sSub[j][i];
  }
  // printmat(subspaceSize, A, "A");
  solve_linear_equation_system(subspaceSize, A, rhsSub, &xSub[0]);
  delete []A;

  residualv.setSize(1);
  getAvMinusFreqSv(frequency, xSub, residualv[0]);
  for(int i=0; i<2*rhs.nvar; i++)
    residualv[0][i] -= rhs[i];
  ergo_real n = template_blas_sqrt(residualv[0]*residualv[0]);
  xTimesRHS = multiplyXtimesVec(rhs);
  do_output(LOG_CAT_INFO, LOG_AREA_LR, 
	    "SetOfEqSolver::getResidual: norm of the residual: %g X*RHS: %g",
            (double)n, (double)xTimesRHS);
  return n >= convThreshold;
}

/** expands above the default limit */
void SetOfEqSolver::increaseSubspaceLimit(int newSize)
{
  if(newSize > maxSubspaceSize) {
    ergo_real *v = new ergo_real[newSize];
    for(int i=0; i<maxSubspaceSize; i++) v[i] = rhsSub[i];
    delete []rhsSub; rhsSub = v;
    LRSolver::increaseSubspaceLimit(newSize);
  }    
}


void
SetOfEqSolver::setRHS(OneElOperator& oper)
{
  ergo_real *res = new ergo_real[nbast*nbast];
  ergo_real *tmp = new ergo_real[nbast*nbast];
  oper.getOper(res);

  ao2mo(nbast, res, tmp);
  commuteWithDMO(nbast, nocc, tmp);
  rhs.setFromFull(nbast, nocc, tmp);
  delete []res;
  delete []tmp;
}

void
SetOfEqSolver::addToSpace(VarVectorCollection& vecs, E2Evaluator& e2)
{
  int orig_size = subspaceSize;
  LRSolver::addToSpace(vecs, e2);
  for(int i=orig_size; i<subspaceSize; i+=2) {
    ergo_real rXX = dot(rhs.nvar, rhs.x(), vects[i].x());
    ergo_real rXY = dot(rhs.nvar, rhs.x(), vects[i].y());
    ergo_real rYX = dot(rhs.nvar, rhs.y(), vects[i].x());
    ergo_real rYY = dot(rhs.nvar, rhs.y(), vects[i].y());
    rhsSub[i]   = rXX + rYY;
    rhsSub[i+1] = rXY + rYX;
  }
}

/** computes polarizability by contracting the response vector with
    specified operator */
ergo_real
SetOfEqSolver::getPolarisability(OneElOperator& oper) /* const */
{
  VarVector x;
  operToVec(oper, x);
  return multiplyXtimesVec(x);
}


/* ===================================================================
 *              Solver of Eigenvalue Problem
 * =================================================================== */
                              
/** get residual of the eigenvalue problem. This is done in following steps:
      Solution in the subspace: [ Xsub, lambda ] = eig(eSub, Ssub);
      Pick first positive eigenvalue l1 = lambda(step+1);
      Pick corresponding eigenvector: Xsub = Xsub(:,step+1);
      Residual Vector:  residualv = (Av-l1*Sv)*Xsub;
*/
bool
EigenSolver::getResidual(VarVectorCollection& residualv)
{
  ergo_real *e = new ergo_real[subspaceSize*subspaceSize];
  ergo_real *s = new ergo_real[subspaceSize*subspaceSize];
  ergo_real *ev = new ergo_real[subspaceSize*subspaceSize];
  ergo_real *alphar = new ergo_real[subspaceSize];
  ergo_real *alphai = new ergo_real[subspaceSize];
  ergo_real *beta   = new ergo_real[subspaceSize];
  int lwork = 16*subspaceSize;
  ergo_real *work = new ergo_real[lwork], *dummy = NULL;
  int i, j;

  for(i=0; i<subspaceSize; i++)
    for(j=0; j<subspaceSize; j++) {
      e[j+i*subspaceSize] = eSub[i][j];
      s[j+i*subspaceSize] = sSub[i][j];
    }
  //printmat(subspaceSize, e, "Reduced E[2]");
  //printmat(subspaceSize, s, "Reduced S[2]");
  mat::ggev("N","V", &subspaceSize, e, &subspaceSize, s, &subspaceSize,
            alphar, alphai, beta,                     /* eigenvalues */
            dummy, &subspaceSize, ev, &subspaceSize, /* eigenvectors */
            work, &lwork, &i);
  delete []work;
  delete []s;
  delete []e;

  if(i == 0) {
    std::vector<int> idx(subspaceSize);
    for(i=0; i<subspaceSize; i++) {
      ritzVals[i] = alphar[i]/beta[i];
      //printf("BEFORE %d: (%g + i %g) / %g -> %10.5f\n", i,
      //         alphar[i], alphai[i], beta[i], ritzVals[i]);
      idx[i] = i;
    }
    for(int ival=0; ival<nStates; ival++) {
      int cidx = ival;
      for(i=ival+1; i< subspaceSize; i++) {
        if( ritzVals[i]>0 &&
            (ritzVals[cidx]<=0 || ritzVals[i] < ritzVals[cidx]) )
          cidx = i;
      }
      ergo_real tr = ritzVals[ival];
      ritzVals[ival] = ritzVals[cidx]; ritzVals[cidx] = tr;
      int ti = idx[ival]; idx[ival] = idx[cidx]; idx[cidx] = ti;
    }

    /* Copy the reduced-space solutions to the private matrix, making
       sure they are normalized properly - very important if
       transition moments are to be computed. ggev returns the vectors
       with somewhat peculiar normalization... */
    if(last_ev) delete []last_ev;
    last_ev = new ergo_real[subspaceSize*nStates]; 
    ergo_real *tmp = new ergo_real[subspaceSize];
    for(i=0; i<nStates; i++) {
      for(int j=0; j<subspaceSize; j++)
        last_ev[j + i*subspaceSize] = ev[j + idx[i]*subspaceSize];
      /* mat:gemv(sSub,  .. */
      for(int j=0; j<subspaceSize; j++) {
        ergo_real s = 0;
        for(int k=0; k<subspaceSize; k++)
          s += sSub[j][k]*last_ev[k+i*subspaceSize];
        tmp[j] = s;
      }
      ergo_real n = 
        1.0/template_blas_sqrt(dot(subspaceSize, tmp, last_ev + i*subspaceSize));
      for(int j=0; j<subspaceSize; j++)
        last_ev[j + i*subspaceSize] *=n;
    }
    delete []tmp;

    /* actually compute the residuals... */
    int nadded=0;
    for(i=nConverged; i<nStates; i++) {
      if(ritzVals[i]<0) {
        printf("CRITICAL!\n");
        throw "Sorting error";
        /* FIXME: no good starting guess!? We should have been
           seeding with more vectors! */
      }
      getAvMinusFreqSv(ritzVals[i], last_ev + i*subspaceSize,
                       residualv[nadded]);
      //ergo_real n = template_blas_sqrt(residualv[0]*residualv[0]);
      ergo_real n = template_blas_sqrt(residualv[nadded]*residualv[nadded]);
      if(n<convThreshold) {
        nConverged++;
        do_output(LOG_CAT_INFO, LOG_AREA_LR, "Converged %d root %g with n=%g",
		  nConverged, (double)ritzVals[i], (double)n);
      } else {
        do_output(LOG_CAT_INFO, LOG_AREA_LR,
                  "Trial vector for %2i th root needed orig %3i, "
                  "n=%g ritz=%f",
		  i+1, idx[i], (double)n, (double)ritzVals[i]);
        nadded++;
      }
    }
    residualv.setSize(nadded);

    delete []alphar;
    delete []alphai;
    delete []beta;
    delete []ev;
    return nConverged<nStates;
  } else {
    delete []alphar;
    delete []alphai;
    delete []beta;
    delete []ev;
    printf("i=%i\n", i); 
    do_output(LOG_CAT_ERROR, LOG_AREA_LR, "dsygv failed with info=%i\n", i);
    throw "Solving projected problem failed. This is not normal.";
  }
}

/** expands above the default limit */
void EigenSolver::increaseSubspaceLimit(int newSize)
{
  if(newSize > maxSubspaceSize) {
    ergo_real *v = new ergo_real[newSize];
    for(int i=0; i<maxSubspaceSize; i++) v[i] = ritzVals[i];
    delete []ritzVals; ritzVals = v;
    v = new ergo_real[newSize];
    for(int i=0; i<maxSubspaceSize; i++) v[i] = transMoms2[i];
    delete []transMoms2; transMoms2 = v;
    LRSolver::increaseSubspaceLimit(newSize);
  }
}


/** generate the starting guess for the HOMO-LUMO excitation by
    placing one in th the right position. Do it the lazy way: creating
    VarVector directly risks inconsistiency when the internal
    representation of VarVector changes. */
int
EigenSolver::getInitialGuess(VarVectorCollection& v)
{
#if 0
  ergo_real *m = new ergo_real[nbast*nbast];
  memset(m, 0, nbast*nbast*sizeof(ergo_real));
  m[nocc + (nocc-1)*nbast] = 1;
  v.setFromFull(nbast, nocc, m);
  delete []m;
#else
  v.setSize(nStates);

  int *idx = new int[e2diag.nvar];
  ergo_real *val = new ergo_real[e2diag.nvar];
  for(int i=0; i<e2diag.nvar; i++)
    val[i] = e2diag[i];

  for(int j=0; j<e2diag.nvar; j++) idx[j] = j;
  for(int i=0; i<nStates; i++) {
    int cidx = i;
    for(int j=i+1; j<e2diag.nvar; j++)
      if(val[j]<val[cidx]) cidx = j;
    ergo_real tr = val[i]; val[i] = val[cidx]; val[cidx] = tr;
    int ti = idx[cidx];  idx[cidx] = idx[i]; idx[i] = ti;
    v[i].setSize(e2diag.nvar);
    v[i] = 0; /* zero the vector */
    v[i][ti] = 1;
  }
  delete []idx;
  delete []val;
#endif
  return 1;
}

void
EigenSolver::computeMoments(OneElOperator& dipx, 
                            OneElOperator& dipy,
                            OneElOperator& dipz)
{
  int state, x;
  ergo_real *proj = new ergo_real[subspaceSize];
  OneElOperator *ops[] = { &dipx, &dipy, &dipz };

  for(state=0; state<nStates; state++) transMoms2[state] = 0;
  for(x=0; x<3; x++) {
    VarVector op;
#if 1
    operToVec(*ops[x], op);
#else
    ergo_real *bufao = new ergo_real[nbast*nbast];
    ergo_real *bufmo = new ergo_real[nbast*nbast];
    ops[x]->getOper(bufao);
    ao2mo(nbast, bufao, bufmo);
    printmat(nbast, bufmo, "OPERATOR IN MO");
    op.setFromFull(nbast, nocc, bufmo);
    delete []bufao;
    delete []bufmo;
#endif
    projectOnSubspace(op, proj);
    for(state=0; state<nStates; state++) {
      ergo_real c = dot(subspaceSize, last_ev + state*subspaceSize, proj);
      transMoms2[state] += c*c;
      printf("STATE %d TRANSITION MOMENT: %20g ENERGY %20f eV\n",
             state+1, (double)c, (double)(ritzVals[state]/UNIT_one_eV));
    }
  }
  delete []proj;
}

} /* end of namespace LR */

