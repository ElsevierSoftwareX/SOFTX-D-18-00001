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

/* -*-mode:c; c-style:bsd; c-basic-offset:4;indent-tabs-mode:nil; -*- */
/** @file grid_interface.h Grid Generator interface. */

#if !defined(_GRID_INTERFACE_H_)
#define _GRID_INTERFACE_H_ 1

#include "realtype.h"
typedef ergo_real real;

/** GridGenMolInfo is an abstract class providing information about
 * the molecule so that the grid generator can fetch atom positions
 * and charges, and shell extents. We prefer to provide virtual
 * functions than just store data in order to reduce storage and need
 * no destructor. This abstract interface also allows to share the
 * code between different programs. */
class GridGenMolInfo {
 public:
    int    noOfAtoms;
    int    noOfBasisFuncs;
    int    noOfShells;

    GridGenMolInfo(int a, int b, int s)
      : noOfAtoms(a), noOfBasisFuncs(b), noOfShells(s) {}

    virtual void getAtom(int icent, int *cnt, real (*coor)[3],
                         int *charge, int *mult) const = 0;
    virtual void setShellRadii(real *shellRadii) const = 0;
    virtual void getBlocks(const real *center, real cellsz,
                           const real *rshell,
                           int *nblcnt, int (*iblcks)[2]) const = 0;
    virtual void getExps(int *maxl, int **nucbas, real (**aa)[2]) const = 0;
    virtual ~GridGenMolInfo() {};
};

#endif /* _GRID_INTERFACE_H_ */
