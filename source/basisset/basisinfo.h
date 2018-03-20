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

/** @file basisinfo.h

    \brief Code for setting up basis functions starting from shells.

    @author: Elias Rudberg <em>responsible</em>.
*/

#ifndef BASISINFO_HEADER
#define BASISINFO_HEADER

/* for NULL */
#include <stdlib.h>

#include "realtype.h"
#include "integral_info.h"
/* for Molecule */
#include "molecule.h"

#include "basisset.h"

struct DistributionSpecStruct {
  ergo_real coeff;           /**< Coefficient A */
  ergo_real exponent;        /**< exponent alfa */
  ergo_real extent;          /* FIXME ELIAS: remove this here since the "extent" you want to use depends on what you are doing? */
  ergo_real centerCoords[3]; /**< x0, y0, z0    */
  char monomialInts[4];  /**< nx, ny, nz    */
};

struct DistributionSpecStructLabeled {
  int basisFuncIndex_1;
  int basisFuncIndex_2;
  int pairIndex;
  int groupID;
  ergo_real limitingFactor; // squareroot of repulsion integral of this distr with itself.
  ergo_real dmatElement;
  DistributionSpecStruct distr;
};


#define MAX_NO_OF_CONTR_GAUSSIANS 20

struct ShellSpecStruct {
  ergo_real coeffList[MAX_NO_OF_CONTR_GAUSSIANS];
  ergo_real exponentList[MAX_NO_OF_CONTR_GAUSSIANS];
  ergo_real sizeList[MAX_NO_OF_CONTR_GAUSSIANS];
  ergo_real padding; /* We keep this for compatibility with old density files... */
  ergo_real centerCoords[3]; /* x0, y0, z0 */
  int noOfContr;
  int shellType; 
  int shell_ID;
  int noOfBasisFuncs;
  int startIndexInMatrix; /* start index in density matrix  */
  int dummy; /* padding to make sure the size of this structure is a multiple of 8 bytes */
};

struct BasisFuncStruct {
  int noOfContr;
  ergo_real coeffList[MAX_NO_OF_CONTR_GAUSSIANS];
  ergo_real exponentList[MAX_NO_OF_CONTR_GAUSSIANS];
  ergo_real extent;
  Vector3D centerCoords; /* x0, y0, z0 */
  int shellType; /* 0 <-> 's', 1 <-> 'p', 2 <-> 'd' etc */
  int functionNumber; /* -1,0,1 for 'p', -2,-1,0,1,2 for 'd', etc */
  int noOfSimplePrimitives;
  int simplePrimitiveIndex;
  int noOfTermsInPolynomial;
  basis_func_term_struct poly[MAX_NO_OF_TERMS_IN_BASIS_FUNC_POLY];
};


struct basis_set_range_struct {
  int startAtomIndex;
  int count;
  basisset_info* basisset;
};

struct BasissetNameRange {
  int startAtomIndex;
  int count;
  char* basisSetFileName;
};

struct BasisInfoStruct {
  int use_6_d_funcs; /**< Whether to use 6 d-type basis functions
			instead of the usual 5 functions. This option
			exists to make it possible to get results
			compatible with other codes that have d-type
			functions defined in that way.  */
  int noOfShells;
  ShellSpecStruct* shellList;
  int noOfBasisFuncs;
  BasisFuncStruct* basisFuncList;
  int noOfSimplePrimitives;
  DistributionSpecStruct* simplePrimitiveList;

  /** Initializes all the fields to sane values. */
  BasisInfoStruct(int use_6_d_funcs_ = 0);

  /** Copies values from another BasisInfoStruct. */
  BasisInfoStruct(const BasisInfoStruct & b);

  ~BasisInfoStruct();

  void addBasisfuncsForPoint(ergo_real x,
			     ergo_real y,
			     ergo_real z,
			     int shellType,
			     ergo_real exponent,
			     const IntegralInfo & integralInfo,
			     int print_raw,
			     int do_normalization,
			     int skip_sort_shells);

  void addBasisfuncsForAtomList(const Atom* atomList,
				int noOfAtoms,
				const basisset_info & basissetDefault,
				int noOfRanges,
				const basis_set_range_struct* rangeList,
				const IntegralInfo & integralInfo,
				int print_raw,
				int do_normalization,
				int skip_sort_shells);

  int addBasisfuncsForMolecule(const Molecule& molecule,
				  const char* basisset_filename_default,
				  int noOfRanges,
				  const BasissetNameRange* rangeList,
				  const IntegralInfo& integralInfo,
				  int print_raw,
				  int do_normalization,
				  int skip_sort_shells);

  static int getNoOfBasisFuncsForAtomType(const basisset_info & basisset,
					  const IntegralInfo & integralInfo,
					  int atomCharge,
					  int use_6_d_funcs);

  BasisInfoStruct *permuteShells(const int *shellMap,
                                 const IntegralInfo& ii) const;

  int normalizeShells(const IntegralInfo& integralInfo);

  int get_basis_funcs();

  int getSimplePrimitivesAll(const IntegralInfo& integralInfo);

  // Stuff needed for Chunks&Tasks usage
  void write_to_buffer ( char * dataBuffer, size_t const bufferSize ) const;
  size_t get_size() const;
  void assign_from_buffer ( char const * dataBuffer, size_t const bufferSize);
};



/** Provides temporary storage for
    compute_integral_of_square_of_basis_func.  Stack used to be the
    storage but many operating systems do not like to allocate so
    much space for stack, particularly when many threads are
    present. */
struct SquareFuncIntegrator {
  const int MAX_NO_OF_PRIMS;
  DistributionSpecStruct *list;
  DistributionSpecStruct *productlist;
SquareFuncIntegrator() :  MAX_NO_OF_PRIMS(44444)
  {
    list = new DistributionSpecStruct[MAX_NO_OF_PRIMS];
    productlist = new DistributionSpecStruct[MAX_NO_OF_PRIMS];
  }
  ~SquareFuncIntegrator() 
  {
    delete []list;
    delete []productlist;
  }
  ergo_real computeIntegralOfSquareOfBasisFunc
  (const IntegralInfo& integralInfo, BasisFuncStruct* basisFunc, int use_6_d_funcs);

  ergo_real getShellFactor(const IntegralInfo& integralInfo,
			   ergo_real exponent, int shellType, int use_6_d_funcs);
};




#ifdef ERGO_ENABLE_DEPRECATED

int basisinfo_construct_multi_basis(BasisInfoStruct* result_basisInfo, 
				    const Molecule* molecule,
				    const char* basisset_filename_default,
				    const Molecule* ghostMolecule,
				    const char* ghost_molecule_basisset_filename,
				    int noOfRanges,
				    const BasissetNameRange* rangeList,
				    IntegralInfo* integralInfo,
				    int print_raw,
				    int do_normalization,
				    int skip_sort_shells,
				    int skip_standard_basis);
struct AtomInfoStruct_{
  int charge;
  ergo_real coords[3];
};
typedef struct AtomInfoStruct_ AtomInfoStruct;
#endif

int get_basis_funcs(BasisInfoStruct* basisInfo,
		    const IntegralInfo* integralInfo,
		    int do_normalization);

int get_simple_primitives_all(BasisInfoStruct* basisInfo,
			      const IntegralInfo* integralInfo);

int output_basisinfo(const BasisInfoStruct & basisInfo);

ergo_real getSafeMaxDistance(const BasisInfoStruct & basisInfo);


#endif
