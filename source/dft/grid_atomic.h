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

/** @file grid_atomic.h Implements shared parts of the grid generation code.
    Examples include radial grids and code for determining active shells
    in given space.
*/
#if !defined(_GRID_ATOMIC_H_)
#define _GRID_ATOMIC_H_ 1

#include "realtype.h"
#include "matrix_typedefs.h"
#include "basisinfo.h"
#include "grid_interface.h"

typedef ergo_real real;
typedef ergo_long_real long_real;

extern const real BraggRadii[];
extern const unsigned BraggSize;


/** RadialScheme describes the radial grid. */
struct RadialScheme {
  const char *name;
  int gridSize;
  explicit RadialScheme(const char *n) : name(n), gridSize(0) {}
  inline int size() const { return gridSize; }
  virtual void init(int myNumber, int charge, real threshold) = 0;
  virtual void generate(real *r, real *w) = 0;
  virtual ~RadialScheme() {}
};

struct RadialSchemeGC2 : public RadialScheme {
  void *quadData;
  RadialSchemeGC2(): RadialScheme("Gauss-Chebychev scheme of second kind")
  {}
  virtual void init(int myNumber, int charge, real threshold);
  virtual void generate(real *r, real *w);
};

struct RadialSchemeTurbo : public RadialScheme {
  real zeta;
  RadialSchemeTurbo(): RadialScheme("Chebychev T2 scheme/M4 mapping (Turbo)")
  {}
  virtual void init(int myNumber, int charge, real threshold);
  virtual void generate(real *r, real *w);
};

struct RadialSchemeLMG : public RadialScheme {
  explicit RadialSchemeLMG(const GridGenMolInfo& ggmi_);

  virtual void init(int myNumber, int charge, real threshold);
  virtual void generate(real *r, real *w);
  virtual ~RadialSchemeLMG();
  private:
  const GridGenMolInfo& ggmi;
  int *nucorb;
  real (*aa)[2];
  int maxL;
  /* grid params */
  real rl, grdc, h, eph;
};



#endif /* _GRID_ATOMIC_H_ */
