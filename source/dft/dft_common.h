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

/** @file dft_common.h

    @brief Common DFT routines. Mostly functional mixing.

    @author: Pawel Salek <em>responsible</em>
*/

#ifndef _DFT_COMMON_H_
#define _DFT_COMMON_H_

#include <stdlib.h>
#include <vector>

#ifdef __cplusplus
#define EXTERN_C extern "C"
#else
#define EXTERN_C
#endif

#include "realtype.h"
#include "basisinfo.h"
#include "matrix_typedefs.h"
#include "functionals.h"
#include "grid_atomic.h"

/** A vector of first order derivatives with respect to two
 * parameters: density rho and SQUARE of the gradient of density grho.
 * zeta_i = |nabla rho_i|^2
 */

typedef struct {
    real fR;  /* d/drho F     */
    real fZ;  /* d/zeta F     */
} FirstDrv;

/* SecondDrv:  matrix  of  second  order functional  derivatives  with
 * respect  to two  parameters: density  rho and  SQUARE  of the
 * density gradient zeta.  The derivatives are computed for alpha-alpha
 * or beta-beta spin-orbital block (i.e. include triplet flag).
 */
typedef struct {
    real fR; /* d/drho  F */
    real fZ; /* d/dzeta F */
    real fRR; /* d/drho^2 F */
    real fRZ; /* d/(drho dzeta) F */
    real fZZ; /* d/dzeta^2 F */
    /* additional derivatives required by  */
    /* general linear response scheme     */
    real fRG; /* d/(drho dgamma) F */
    real fZG; /* d/(dzeta dgamma) F */
    real fGG; /* d/dgamma^2 F */
    real fG;  /* d/dgamma F */
} SecondDrv;


EXTERN_C void dftpot0_(FirstDrv *ds, const real* weight, const FunDensProp* dp);
EXTERN_C void dftpot1_(SecondDrv *ds, const real* w, const FunDensProp* dp,
		       const int* triplet);

EXTERN_C void dft_init(void);
EXTERN_C int dft_setfunc(const char *line);

class ShellTree;

/** Ergo specific implementation of molecule-grid interface. */
class ErgoMolInfo : public GridGenMolInfo {
  const BasisInfoStruct& bis;
  const Molecule&        molecule;
 public:
  ErgoMolInfo(const BasisInfoStruct& bis_,  const Molecule& mol);
  virtual ~ErgoMolInfo();

  virtual void getAtom(int icent, int *cnt, real (*coor)[3],
                       int *charge, int *mult) const;
  virtual void setShellRadii(real *shellRadii) const;
  virtual void getBlocks(const real *center, real cellsz,
                         const real *rshell,
                         int *nblcnt, int (*iblcks)[2]) const;
  void getBlocks1(const real *center, real cellsz,
                  const real *rshell,
                  int *nblcnt, int (*iblcks)[2]) const;
  virtual void getExps(int *maxl, int **nucbas, real (**aa)[2]) const;
  ShellTree *shellTree;
};

EXTERN_C void ergoShellsToOrbs(const int *nshlbl, const int (*shlblock)[2],
                               int *norbbl, int (*orbblock)[2],
                               const BasisInfoStruct& bis);

EXTERN_C int dft_get_num_threads();
EXTERN_C void dft_set_num_threads(int nThreads);


EXTERN_C void dft_init(void);

#define dal_new(sz,tp) (tp*)dal_malloc_((sz)*sizeof(tp),__FUNCTION__, __LINE__)
void* dal_malloc_(size_t sz, const char *func, unsigned line);
/* dal_malloc: usage discouraged */
#define dal_malloc(sz) dal_malloc_((sz),__FUNCTION__, __LINE__)

/* useful  constants for BLAS interfacing */
extern int  ZEROI, ONEI, THREEI, FOURI;
extern real ZEROR, ONER, TWOR, FOURR;

/** Class Box provides an ability to determine box containing all
    Objects.  The class Object must provide field center[] and method
    radius(). */
class Box {
public:
  real getDistanceTo(const real* v) const;
  int getMaxDim() const;
  real size(int dim) const { return hi[dim]-lo[dim]; }

  bool overlapsWith(const real *center, real radius) const {
    real d = getDistanceTo(center);
    return d < radius;
  }

  /** Determines whether given point is inside the box. In order to
      avoid double counting, the points that are overlap with the
      lower limits are included but those that overlap with the higher
      limit are excluded. */
  bool contains(const real *p) const {
#if 0
    printf("B:(%8.2f %8.2f %8.2f)-(%8.2f %8.2f %8.2f): %8.2f %8.2f %8.2f ",
           lo[0], lo[1], lo[2], hi[0], hi[1], hi[2],
           p[0], p[1], p[2]);
#endif
    for(int i=0; i<3; i++)
      if(p[i]<lo[i] || p[i] >= hi[i]) {
        //puts("F");
        return false;
      }
    //puts("T");
    return true;
  }

  real lo[3];
  real hi[3];
};

template<typename Iterator>
  void getBoundingBox(Box& box, Iterator start, Iterator end)
{
  static const ergo_real OFF = 0.1;
  if(start == end)
    throw "BoundingBox called for empty set";

  real r = start->radius() + OFF;
  for(int i=0; i<3; i++) {
    box.lo[i] = start->center[i]-r;
    box.hi[i] = start->center[i]+r;
  }

  for(++start; start != end; ++start) {
    real r = start->radius() + OFF;
    for(int i=0; i<3; i++) {
      real l = start->center[i]-r; if (l<box.lo[i]) box.lo[i] = l;
      real h = start->center[i]+r; if (h>box.hi[i]) box.hi[i] = h;
    }
  }
}


int sync_threads(bool release, int nThreads);

#endif /* _DFT_COMMON_H_ */
