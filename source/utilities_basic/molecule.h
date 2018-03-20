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

/** @file molecule.h

    @brief Class representing a molecule as a set of atoms with
    assiciated coordinates and charges of the atomic nuclei.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef MOLECULE_HEADER
#define MOLECULE_HEADER

#include <cmath>
#include <vector>
#include <cassert>

#include "realtype.h"
#include "mat_gblas.h"

/** Simple atom representation by its charge and cartesian coordinates.
 *
 */
struct Atom {
  ergo_real charge;
  ergo_real coords[3];
};

/**
 * A representation of Vector or point in cartesian space. It provides
 * means to compute distance between two pointsi space.
 */
struct Vector3D {
  ergo_real v[3];
  Vector3D() {}
  Vector3D(ergo_real x, ergo_real y, ergo_real z) {
    v[0] = x; v[1] = y; v[2] = z;
  }
  ergo_real& operator[](unsigned i)       { return v[i]; }
  ergo_real  operator[](unsigned i) const { return v[i]; }
  /** compute square of distance between two points. */
  ergo_real dist2(const ergo_real b[]) const {
    ergo_real d, r;
    d = v[0]-b[0]; r  = d*d;
    d = v[1]-b[1]; r += d*d;
    d = v[2]-b[2]; r += d*d;
    return r;
  }
  /** compute distance between two points. */
  ergo_real dist(const Vector3D& b) const 
  { return template_blas_sqrt(dist2(b.v)); }
  ergo_real dist(const ergo_real b[]) const 
  { return template_blas_sqrt(dist2(b)); }
};

/**
 * Representation of a molecule as a set of nuclei and total
 * charge. It provides I/O methods and basic manipulation routines.
 */
class Molecule {
 private:
  std::vector<Atom> atoms;
  ergo_real netCharge;
  int noOfAtoms;
  
 public:

 Molecule() : atoms(10), netCharge(0), noOfAtoms(0) {}

  void addAtom(ergo_real c, ergo_real x, ergo_real y, ergo_real z) {
    int currListSize = atoms.size();
    if(noOfAtoms >= currListSize)
      atoms.resize(currListSize*2);
    atoms[noOfAtoms].charge = c;
    atoms[noOfAtoms].coords[0] = x;
    atoms[noOfAtoms].coords[1] = y;
    atoms[noOfAtoms].coords[2] = z;
    noOfAtoms++;
  }

  void clear() { noOfAtoms = 0; netCharge = 0; }
  void setNetCharge(ergo_real netCharge_) { netCharge = netCharge_; }
  void replaceAtom(int i, const Atom & atom) { assert(i >= 0 && i < noOfAtoms); atoms[i] = atom; }
  void setAtomList(const std::vector<Atom> atomList) { atoms = atomList; noOfAtoms = atomList.size(); }
  const Atom* getAtomListPtr() const { return &atoms[0]; }
  const Atom & getAtom(int i) const { return atoms[i]; }
  int getNoOfAtoms() const { return noOfAtoms; }
  ergo_real getNetCharge() const { return netCharge; }

  /** Compute smallest and largest internuclear distances. */
  void getExtremeInternuclearDistancesQuadratic(ergo_real & minDist, ergo_real & maxDist) const;
  /** Compute nuclear repulsion energy. */
  ergo_real getNuclearRepulsionEnergyQuadratic() const;
  /** Compute nuclear energy in given electric field. */
  ergo_real getNuclearElectricFieldEnergy(const Vector3D& electricField) const;
  /** Compute total number of electrons. The result is sum of atomic
      charges plus netCharge. */
  int getNumberOfElectrons() const;
  /** Compute gradient of nuclear repulsion energy w.r.t. changes in nuclear coordinates. Result is added to resultGradient vector.  */
  void getNuclearRepulsionEnergyGradientContribQuadratic(ergo_real* resultGradient) const;

  /** Loads molecule from a given file name, assuming given net
      charge. basissetFile will be set if the file contains basis set
      and basissetFile is NULL. */
  int setFromMoleculeFile(const char* fileName,
                          int netCharge,
                          char **basissetFile);

};


#endif /* MOLECULE_HEADER */
