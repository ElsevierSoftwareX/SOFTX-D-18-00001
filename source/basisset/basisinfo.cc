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

/** @file basisinfo.cc

    \brief Code for setting up basis functions starting from shells.

    @author: Elias Rudberg <em>responsible</em>.
*/

/* -*-mode:c; c-style:k&r; indent-tabs-mode: nil -*- */
/* Written by Elias Rudberg, KTH, Stockholm */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <memory.h>
#include <time.h>
#include <stdarg.h>
#include <string.h>

#include "basisinfo.h"
#include "basisset.h"
#include "memorymanag.h"
#include "pi.h"
#include "output.h"
#include "utilities.h"
#include "boysfunction.h"
#include "integral_info.h"
#include "integrals_general.h"
#include "machine_epsilon.h"


int
output_basisinfo(const BasisInfoStruct & basisInfo)
{
  static char shell_names[] = "SPDFGHIJKLMNOP";
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "=============== start of output_basisinfo ===========================");
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "output_basisinfo, basisInfo->noOfShells = %i",
	    basisInfo.noOfShells);
  char line[180];
  for(int i = 0; i < basisInfo.noOfShells; i++) {
    if(basisInfo.shellList[i].shellType >=0 &&
       basisInfo.shellList[i].shellType<(signed)sizeof(shell_names))
      do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "%c-shell at (x y z) = ( %8.4f %8.4f %8.4f )"
		" %d primitive(s)",
		shell_names[basisInfo.shellList[i].shellType],
		(double)basisInfo.shellList[i].centerCoords[0],
		(double)basisInfo.shellList[i].centerCoords[1],
		(double)basisInfo.shellList[i].centerCoords[2],
		basisInfo.shellList[i].noOfContr);
    else
      do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "shell with L=%d at (x y z) = ( %8.4f %8.4f %8.4f )"
		" %d primitive(s)",
		basisInfo.shellList[i].shellType,
		(double)basisInfo.shellList[i].centerCoords[0],
		(double)basisInfo.shellList[i].centerCoords[1],
		(double)basisInfo.shellList[i].centerCoords[2],
		basisInfo.shellList[i].noOfContr);
    int pos = 0;
    for(int j = 0; j < basisInfo.shellList[i].noOfContr; j++) {
      sprintf(line+pos, "%10.6f e^%9.6f,",
	      (double)basisInfo.shellList[i].coeffList[j],
	      (double)basisInfo.shellList[i].exponentList[j]);
      pos = (int)strlen(line);
      if(pos>60) { do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, line); pos = 0; }
    }
    if(pos>0) do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, line);
  }
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "================ end of output_basisinfo ===========================");
  return 0;
}


static void
define_basis_func_poly(BasisFuncStruct* basisFunc, int polyIndex, const IntegralInfo& b)
{
  if(polyIndex >= b.no_of_basis_func_polys)
    throw "Error in define_basis_func_poly: (polyIndex >= b.no_of_basis_func_polys).";
  const basis_func_poly_struct* poly = &b.basis_func_poly_list[polyIndex];
  basisFunc->noOfTermsInPolynomial = poly->noOfTerms;
  for(int i = 0; i < poly->noOfTerms; i++)
    memcpy(&basisFunc->poly[i], &poly->termList[i], sizeof(basis_func_term_struct));
}


static void
define_basis_func_poly_special_6dfuncs(BasisFuncStruct* basisFunc, const IntegralInfo & b)
{
  // SPECIAL CASE: USE 6 d-type functions instead of 5.
  // Use the following 6 functions: x^2 y^2 z^2 xy xz yz
  basisFunc->noOfTermsInPolynomial = 1;
  int i0 = 0, i1 = 0, i2 = 0;
  ergo_real coeff = 0;
  const ergo_real coeff_a = 1;
  const ergo_real coeff_b = template_blas_sqrt((ergo_real)3);
  switch(basisFunc->functionNumber) {
  case 0: i0 = 2; i1 = 0; i2 = 0; coeff = coeff_a; break; // x^2
  case 1: i0 = 0; i1 = 2; i2 = 0; coeff = coeff_a; break; // y^2
  case 2: i0 = 0; i1 = 0; i2 = 2; coeff = coeff_a; break; // y^2
  case 3: i0 = 1; i1 = 1; i2 = 0; coeff = coeff_b; break; // xy
  case 4: i0 = 1; i1 = 0; i2 = 1; coeff = coeff_b; break; // xz
  case 5: i0 = 0; i1 = 1; i2 = 1; coeff = coeff_b; break; // yz
  default: throw "Error: default reached when defining d-type basis function.";
  }
  basisFunc->poly[0].coeff = coeff;
  basisFunc->poly[0].monomialInts[0] = i0;
  basisFunc->poly[0].monomialInts[1] = i1;
  basisFunc->poly[0].monomialInts[2] = i2;
  basisFunc->poly[0].monomialID = b.monomial_info.monomial_index_list[i0][i1][i2];
}


static int 
get_simple_primitives(BasisFuncStruct* currBasisFunc,
		      DistributionSpecStruct* list,
		      int nInput,
		      int nListMax,
		      const IntegralInfo& b,
                      int use_6_d_funcs)
{
  /* first setup polynomial */
  int spd = currBasisFunc->shellType;
  // Note: special case for use_6_d_funcs comes in here.
  if(spd == 2 && use_6_d_funcs == 1) {
    // SPECIAL CASE: USE 6 d-type functions instead of 5.
    // Use the following 6 functions: x^2 y^2 z^2 xy xz yz
    define_basis_func_poly_special_6dfuncs(currBasisFunc, b);
  }
  else {
    int baseIndex = spd*spd;
    define_basis_func_poly(currBasisFunc, baseIndex + currBasisFunc->functionNumber, b);
  }
  int n = nInput;
  int contr = currBasisFunc->noOfContr;
  for(int kk = 0; kk < contr; kk++) {
    for(int ii = 0; ii < currBasisFunc->noOfTermsInPolynomial; ii++) {
      list[n].coeff = currBasisFunc->coeffList[kk] * currBasisFunc->poly[ii].coeff;
      for(int coordNo = 0; coordNo < 3; coordNo++)
	list[n].monomialInts[coordNo] = currBasisFunc->poly[ii].monomialInts[coordNo];
      list[n].exponent = currBasisFunc->exponentList[kk];
      for(int coordNo = 0; coordNo < 3; coordNo++)
	list[n].centerCoords[coordNo] = currBasisFunc->centerCoords[coordNo];
      n++;
      if(n >= nListMax) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_simple_primitives: "
		  "(n >= nListMax)");
	return -1;
      }
    } /* END FOR ii */
  } /* END FOR kk */
  return n - nInput;
} /* END get_simple_primitives */


static int
sort_shells(ShellSpecStruct* list, ShellSpecStruct* listTemp, int n)
{
  int n1, n2, i1, i2;
  if(n == 1)
    return 0;
  n1 = n / 2;
  n2 = n - n1;

  /* sort parts */
  if(sort_shells(&list[0], listTemp, n1) != 0)
    return -1;
  if(sort_shells(&list[n1], listTemp, n2) != 0)
    return -1;

  /* merge to temp list */
  i1 = 0;
  i2 = 0;
  while(i1 < n1 && i2 < n2)
    {
      if(list[i1].shell_ID < list[n1+i2].shell_ID)
	{
	  memcpy(&listTemp[i1+i2], &list[i1], sizeof(ShellSpecStruct));
	  i1++;
	}
      else
	{
	  memcpy(&listTemp[i1+i2], &list[n1+i2], sizeof(ShellSpecStruct));
	  i2++;
	}
    } /* END WHILE */
  while(i1 < n1)
    {
      memcpy(&listTemp[i1+i2], &list[i1], sizeof(ShellSpecStruct));
      i1++;
    }
  while(i2 < n2)
    {
      memcpy(&listTemp[i1+i2], &list[n1+i2], sizeof(ShellSpecStruct));
      i2++;
    }

  /* copy back to result list */
  memcpy(list, listTemp, n * sizeof(ShellSpecStruct));
  return 0;
}


/* FIXME: why was there an anonymous name space here, enclosing the
   SquareFuncIntegrator stuff? */
namespace {

} /* end of anonymous name space */


ergo_real 
SquareFuncIntegrator::computeIntegralOfSquareOfBasisFunc
(const IntegralInfo& integralInfo, BasisFuncStruct* basisFunc, int use_6_d_funcs)
{
  int noOfPrimitives = get_simple_primitives(basisFunc,
					     list,
					     0,
					     MAX_NO_OF_PRIMS,
					     integralInfo,
                                             use_6_d_funcs);
  if(noOfPrimitives == 0)
    throw "error in get_simple_primitives in computeIntegralOfSquareOfBasisFunc (noOfPrimitives == 0)";
  if(noOfPrimitives < 0)
    throw "error in get_simple_primitives in computeIntegralOfSquareOfBasisFunc (noOfPrimitives < 0)";  

  // Compute square of basis function
  int productCount = 0;
  for(int ii = 0; ii < noOfPrimitives; ii++) {
    const DistributionSpecStruct& primA = list[ii];
    for(int jj = 0; jj < noOfPrimitives; jj++) {
      const DistributionSpecStruct& primB = list[jj];
      int nNewPrims = get_product_simple_prims(primA, 
					       primB, 
					       &productlist[productCount],
					       MAX_NO_OF_PRIMS - productCount,
					       0);
      if(nNewPrims < 0)
	throw "Error in computeIntegralOfSquareOfBasisFunc, in get_product_simple_prims.";
      productCount += nNewPrims;
    } // END FOR jj
  } // END FOR ii
  ergo_real sum = 0;
  for(int ii = 0; ii < productCount; ii++)
    sum += compute_integral_of_simple_prim(productlist[ii]);
  if(sum < 0)
    throw "Error in computeIntegralOfSquareOfBasisFunc, norm factor sum < 0.";
  return sum;
}

ergo_real
SquareFuncIntegrator::getShellFactor(const IntegralInfo& integralInfo,
                                     ergo_real exponent,
                                     int shellType,
                                     int use_6_d_funcs)
{
  BasisFuncStruct basisFunc;
  basisFunc.noOfContr = 1;
  basisFunc.coeffList[0] = 1;
  basisFunc.exponentList[0] = exponent;
  for(int kk = 0; kk < 3; kk++)
    basisFunc.centerCoords[kk] = 0;
  basisFunc.shellType = shellType;
  basisFunc.functionNumber = 0;
  // Compute integral of this basis function squared.
  ergo_real integralValue = computeIntegralOfSquareOfBasisFunc(integralInfo, &basisFunc, use_6_d_funcs);
  ergo_real shellFactor = (ergo_real)1.0 / template_blas_sqrt(integralValue);
  return shellFactor;
}



static int
find_range_index(int atomIndex, int noOfRanges, const basis_set_range_struct* rangeList) {
  for(int i = 0; i < noOfRanges; i++) {
    if(atomIndex >= rangeList[i].startAtomIndex && atomIndex < rangeList[i].startAtomIndex + rangeList[i].count)
      return i;
  }
  // Return -1 to indicate range not found.
  return -1;
}


static const basisset_info* 
select_basis_set(int atomIndex, 
		 int noOfRanges,
		 const basis_set_range_struct* rangeList,
		 const basisset_info* basissetDefault)
{
  int rangeIndex = find_range_index(atomIndex, noOfRanges, rangeList);
  if(rangeIndex < 0)
    return basissetDefault;
  else
    return rangeList[rangeIndex].basisset;
}

/** Returns number of shells needed to describe the electronic density
    for given molecule and basis set.

    @param atomList list of atoms

    @param noOfAtoms the length of atomList

    @param basissetDefault the basis set to be used for all atoms but
    those specified by rangeList.

    @param noOfRanges the length of rangeList.

    @param rangeList A list of atoms that should get some other,
    specified basis set.

    @return the number of basis set shells.
*/
static int 
setup_shells_multi_basis_getcount(const Atom* atomList,
				  int noOfAtoms,
				  const basisset_info* basissetDefault,
				  int noOfRanges,
				  const basis_set_range_struct* rangeList)
{
  int noOfShells = 0;
  for(int i = 0; i < noOfAtoms; i++) {
    int z = (int)atomList[i].charge;
    const basisset_info* basissetCurrAtom = select_basis_set(i, noOfRanges, rangeList, basissetDefault);
    int noOfShellsCurrAtom = basissetCurrAtom->atoms[z].noOfShells;
    if(noOfShellsCurrAtom <= 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, 
		"error in setup_shells_multi_basis_getcount: element %i is not supported by selected basis set?", z);
      return -1;
    }
    noOfShells += noOfShellsCurrAtom;
  }
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "setup_shells_multi_basis_getcount, noOfShells = %i", noOfShells);
  return noOfShells;
}


static int
setup_shells_multi_basis(const IntegralInfo& integralInfo,
			 const Atom* atomList,
			 int noOfAtoms,
			 const basisset_info* basissetDefault, 
			 ShellSpecStruct* shell_list,
			 int noOfShells,
			 int noOfRanges,
			 const basis_set_range_struct* rangeList,
                         int use_6_d_funcs)
{
  memset(shell_list, 0, noOfShells*sizeof(ShellSpecStruct));
  
  int count = 0;
  SquareFuncIntegrator sfi;
  for(int i = 0; i < noOfAtoms; i++) {
    int z = (int)atomList[i].charge;
    const basisset_info* basissetCurrAtom = select_basis_set(i, noOfRanges, rangeList, basissetDefault);
    int noOfShellsCurrAtom = basissetCurrAtom->atoms[z].noOfShells;
    if(noOfShellsCurrAtom <= 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in setup_shells_multi_basis: element %i is not supported by selected basis set?", z);
      return -1;
    }
    for(int j = 0; j < noOfShellsCurrAtom; j++) {
      if(count >= noOfShells) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in setup_shells_multi_basis: (count >= noOfShells).");
	return -1;
      }
      const basisset_shell_struct* basissetShell = 
	&basissetCurrAtom->atoms[z].shells[j];
      shell_list[count].shellType = basissetShell->type;

      // Note: special case for use_6_d_funcs comes in here.
      if(basissetShell->type == 2 && use_6_d_funcs == 1)
        shell_list[count].noOfBasisFuncs = 6;
      else
        shell_list[count].noOfBasisFuncs = 1 + basissetShell->type * 2;

      if ( basissetShell->contrCount > MAX_NO_OF_CONTR_GAUSSIANS )
	throw std::runtime_error("basissetShell->contrCount > MAX_NO_OF_CONTR_GAUSSIANS in setup_shells_multi_basis(...)");
      shell_list[count].noOfContr = basissetShell->contrCount;
      shell_list[count].shell_ID  = basissetShell->shell_ID;
      for(int k = 0; k < 3; k++)
	shell_list[count].centerCoords[k] = 
	  atomList[i].coords[k];
      for(int k = 0; k < basissetShell->contrCount; k++) {
	ergo_real exponent = basissetShell->exponentList[k];
	shell_list[count].coeffList[k] = 
	  basissetShell->coeffList[k] * sfi.getShellFactor(integralInfo, 
                                                           exponent,
                                                           basissetShell->type, use_6_d_funcs);
	shell_list[count].exponentList[k] = exponent;
      }
      shell_list[count].startIndexInMatrix = -1; // startIndexInMatrix will be set later.
      count++;
    }
  }
  
  if(count != noOfShells)
    return -1;
  
  return 0;
}



void BasisInfoStruct::addBasisfuncsForPoint(ergo_real x,
					    ergo_real y,
					    ergo_real z,
					    int shellType,
					    ergo_real exponent,
					    const IntegralInfo & integralInfo,
					    int print_raw,
					    int do_normalization,
					    int skip_sort_shells) {
  // Create "dummy atom" at given pos.
  Atom dummyAtom;
  dummyAtom.charge = 1;
  dummyAtom.coords[0] = x;
  dummyAtom.coords[1] = y;
  dummyAtom.coords[2] = z;
  // Create "dummy basis set" defining shell types for atoms with charge 1.
  basisset_info basisset;
  basisset.atoms.resize(2);
  for(int kk = 0; kk < 2; kk++) {
    // Use one shell of each type, up to given shellType.
    int nShells = shellType + 1;
    if(nShells >= MAX_NO_OF_SHELLS_PER_ATOM)
      throw std::runtime_error("Error in BasisInfoStruct::addBasisfuncsForPoint: (nShells >= MAX_NO_OF_SHELLS_PER_ATOM).");
    basisset.atoms[kk].noOfShells = nShells;
    for(int i = 0; i <= shellType; i++) {
      basisset.atoms[kk].shells[i].type = i;
      basisset.atoms[kk].shells[i].contrCount = 1;
      basisset.atoms[kk].shells[i].shell_ID = i;
      basisset.atoms[kk].shells[i].exponentList[0] = exponent;
      basisset.atoms[kk].shells[i].coeffList[0] = 1;
    }
  }
  // OK, now the "dummy basis set" is done.
  addBasisfuncsForAtomList(&dummyAtom, // atomList
			   1, // int noOfAtoms,
			   basisset,
			   0, // noOfRanges,
			   NULL, // rangeList,
			   integralInfo,
			   print_raw,
			   do_normalization,
			   skip_sort_shells);
}


void BasisInfoStruct::addBasisfuncsForAtomList(const Atom* atomList,
					       int noOfAtoms,
					       const basisset_info & basissetDefault,
					       int noOfRanges,
					       const basis_set_range_struct* rangeList,
					       const IntegralInfo& integralInfo,
					       int print_raw,
					       int do_normalization,
					       int skip_sort_shells) {

  int noOfShellsToAdd = setup_shells_multi_basis_getcount(atomList,
							  noOfAtoms,
                                                          &basissetDefault,
                                                          noOfRanges,
                                                          rangeList);
  if(noOfShellsToAdd <= 0)
    throw std::runtime_error("error in setup_shells_multi_basis_getcount");

  int noOfShellsNew = noOfShells + noOfShellsToAdd;
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "BasisInfoStruct::addBasisfuncsForAtomList, noOfShellsNew = %i", noOfShellsNew);
  
  ShellSpecStruct* shell_list_new = new ShellSpecStruct[noOfShellsNew];
  if(this->noOfShells > 0)
    memcpy(shell_list_new, this->shellList, this->noOfShells*sizeof(ShellSpecStruct));
  
  // Setup new shells
  if(setup_shells_multi_basis(integralInfo,
			      atomList,
			      noOfAtoms,
                              &basissetDefault,
                              &shell_list_new[noOfShells],
                              noOfShellsToAdd,
                              noOfRanges,
                              rangeList,
                              use_6_d_funcs) != 0)
    throw std::runtime_error("error in setup_shells_multi_basis");

  for(int i = 0; i < noOfShellsNew; i++) {
    for(int kk = 0; kk < shell_list_new[i].noOfContr; kk++) {
      // calculate size
      shell_list_new[i].sizeList[kk] = 
	template_blas_fabs(template_blas_pow((ergo_real)pi/shell_list_new[i].exponentList[kk], (ergo_real)1.5) * shell_list_new[i].coeffList[kk]);
    }
  } /* END FOR i */
  
  if(skip_sort_shells == 0) {
    /* sort shells by shell ID */
    std::vector<ShellSpecStruct> shellListTemp(noOfShellsNew);
    if(sort_shells(shell_list_new, &shellListTemp[0], noOfShellsNew) != 0)
      throw std::runtime_error("Error in sort_shells.");
  }

  if(this->shellList)
    delete [] this->shellList;
  this->shellList = shell_list_new;
  this->noOfShells = noOfShellsNew;
  
  if(do_normalization) {
    if(this->normalizeShells(integralInfo) != 0)
      throw std::runtime_error("error in normalizeShells");
  }
  if(this->get_basis_funcs() != 0)
    throw std::runtime_error("error in get_basis_funcs");
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "get_basis_funcs returned OK, number of basis funcs: %i",
	    this->noOfBasisFuncs);
  if(this->getSimplePrimitivesAll(integralInfo) != 0)
    throw std::runtime_error("error in getSimplePrimitivesAll");
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "getSimplePrimitivesAll returned OK, n = %i",
	    this->noOfSimplePrimitives);
}





/** Fills in BasisInfoStruct for given molecule and
    basisset_filename. It can be called several times to add basis functions 
    for ghost molecules.
    
    @param molecule contains the description of the molecule geometry.

    @param basisset_filename_default contains the name of the basis
    set that will be used for atoms that have no basis set specified
    in rangeList. A number of directories will be searched for the
    given basis.

    @param noOfRanges the length of rangeList.

    @param rangeList is a list of basis sets associated with ranges of
    atoms that should get non-default basis set.

    @param integralInfo - the core structure for integral
    evaluation, needed for basis set normalization.

    @param print_raw - whether the basis set as read should be printed.

    @param do_normalization - whether the contraction coefficients in
    front of exponentials are to be normalized.
    
    @param skip_sort_shells disable the standard sorting of shells in
    the basis set with respect to atom type and exponent.

    @return 0 on success, -1 on failure.
*/
int BasisInfoStruct::addBasisfuncsForMolecule(const Molecule& molecule,
                                              const char* basisset_filename_default,
                                              int noOfRanges,
                                              const BasissetNameRange* rangeList,
                                              const IntegralInfo& integralInfo,
                                              int print_raw,
                                              int do_normalization,
                                              int skip_sort_shells)
{
  static const char *dirv[] = {
    ".", "basis", "../basis",
    ERGO_DATA_PREFIX "/basis",
    ERGO_DATA_PREFIX,
    ERGO_SPREFIX "/basis",
    ERGO_SPREFIX
  };
  basisset_info* basissetDefault = new basisset_info();
  
  std::vector<basis_set_range_struct> rangeListTemp(noOfRanges);
  memset(&rangeListTemp[0], 0, noOfRanges*sizeof(basis_set_range_struct));
  
  if(noOfRanges > 0 && rangeList == NULL)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in BasisInfoStruct::addBasisfuncsForMolecule: (noOfRanges > 0 && rangeList == NULL).");
      delete basissetDefault;
      return -1;
    }

  if(read_basisset_file(*basissetDefault, basisset_filename_default, 6, dirv,
                        print_raw) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file for standard basis set.");
      delete basissetDefault;
      return -1;
    }

  for(int rangeIndex = 0; rangeIndex < noOfRanges; rangeIndex++) {
    rangeListTemp[rangeIndex].startAtomIndex = rangeList[rangeIndex].startAtomIndex;
    rangeListTemp[rangeIndex].count = rangeList[rangeIndex].count;
    if(rangeList[rangeIndex].count <= 0)
      rangeListTemp[rangeIndex].basisset = NULL;
    else {
      rangeListTemp[rangeIndex].basisset = new basisset_info();
      if(read_basisset_file(*rangeListTemp[rangeIndex].basisset, rangeList[rangeIndex].basisSetFileName, 6, dirv,
			    print_raw) != 0) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in read_basisset_file for rangeIndex = %i", rangeIndex);
	delete basissetDefault;
	return -1;
      }
    }
  }

  addBasisfuncsForAtomList(molecule.getAtomListPtr(),
			   molecule.getNoOfAtoms(),
			   *basissetDefault,
			   noOfRanges,
			   &rangeListTemp[0],
			   integralInfo,
			   print_raw,
			   do_normalization,
			   skip_sort_shells);

  delete basissetDefault;
  
  return 0;
}


int
BasisInfoStruct::getNoOfBasisFuncsForAtomType(const basisset_info & basisset,
					      const IntegralInfo & integralInfo,
					      int atomCharge,
					      int use_6_d_funcs) {
  BasisInfoStruct b(use_6_d_funcs);
  const int do_basis_normalization = 1;
  const int skip_sort_shells = 0;
  // Create dummy atom with the given charge
  Atom atom;
  atom.charge = atomCharge;
  for(int i = 0; i < 3; i++)
    atom.coords[i] = 0;
  b.addBasisfuncsForAtomList(&atom, 1, basisset, 0, NULL, integralInfo, 0, do_basis_normalization, skip_sort_shells);
  return b.noOfBasisFuncs;
}


/** a factory method generating new BasisInfo struct with permuted
    shells and basis functions.
    
    @param shellMap vector defining the permutation of shells.

    newShell(i) = this.shell(shellMap(i));

    @param ii IntegralInfo structure needed to reconstruct the
    primitive gaussian data.
*/
BasisInfoStruct*
BasisInfoStruct::permuteShells(const int *shellMap,
                               const IntegralInfo& ii) const
{
  BasisInfoStruct *res = new BasisInfoStruct(use_6_d_funcs);

  res->noOfShells = noOfShells;
  res->shellList  = new ShellSpecStruct[noOfShells];

  for(int i = 0; i<noOfShells; i++)
    res->shellList[i] = shellList[shellMap[i]];

  if(res->get_basis_funcs() != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_basis_funcs");
    delete res;
    return NULL;
  }
  if(res->getSimplePrimitivesAll(ii) != 0)  {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS,
              "error in getSimplePrimitivesAll");
    delete res;
    return NULL;
  }
  return res;
}



/** Normalizes shells so that the overlap of each basis function with itself will be 1.
    This is done by explicitly generating each basis function in each shell and
    computing the overlap. It is verified that all functions within the same shell 
    have the same normalization factor.
*/
int BasisInfoStruct::normalizeShells(const IntegralInfo& integralInfo)
{
  ergo_real normFactorTot_min = 1e22;
  ergo_real normFactorTot_max = 0;

  // Adapt tolerance to machine accuracy to be able to run with different precision.
  ergo_real max_allowed_difference = template_blas_sqrt(get_machine_epsilon());
  SquareFuncIntegrator sfi;
  for(int i = 0; i < this->noOfShells; i++) {
    ShellSpecStruct* currShell = &this->shellList[i];

    ergo_real normFactorShell_min = 1e22;
    ergo_real normFactorShell_max = 0;

    int nFunctions = currShell->noOfBasisFuncs;
    if(nFunctions <= 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in BasisInfoStruct::normalizeShells: (nFunctions <= 0).");
      return -1;
    }

    for(int j = 0; j < nFunctions; j++) {
      BasisFuncStruct basisFunc;
      basisFunc.noOfContr = currShell->noOfContr;
      for(int kk = 0; kk < currShell->noOfContr; kk++) {
	basisFunc.coeffList[kk] = currShell->coeffList[kk];
	basisFunc.exponentList[kk] = currShell->exponentList[kk];
      } /* END FOR kk */
      for(int kk = 0; kk < 3; kk++)
	basisFunc.centerCoords[kk] = currShell->centerCoords[kk];
      basisFunc.shellType = currShell->shellType;
      basisFunc.functionNumber = j;

      // Compute integral of this basis function squared, for normalization.
      ergo_real integralValue =
        sfi.computeIntegralOfSquareOfBasisFunc(integralInfo, &basisFunc, use_6_d_funcs);
      ergo_real normalizationFactor = (ergo_real)1.0 / template_blas_sqrt(integralValue);

      if(normalizationFactor < normFactorShell_min)
	normFactorShell_min = normalizationFactor;
      if(normalizationFactor > normFactorShell_max)
	normFactorShell_max = normalizationFactor;          
    } /* END FOR j */
      
    ergo_real absdiff = template_blas_fabs(normFactorShell_max - normFactorShell_min);
    if(absdiff > max_allowed_difference) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in BasisInfoStruct::normalizeShells: different norm factors within shell, absdiff = %g.", (double)absdiff);
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "normFactorShell_min = %22.15g  normFactorShell_max = %22.15g", (double)normFactorShell_min, (double)normFactorShell_max);
      return -1;
    }
      
    // Use average. This should not matter, they should be the same anyway.
    ergo_real normalizationFactor = 0.5 * (normFactorShell_max + normFactorShell_min);
      
    if(normalizationFactor < normFactorTot_min)
      normFactorTot_min = normalizationFactor;
    if(normalizationFactor > normFactorTot_max)
      normFactorTot_max = normalizationFactor;
      
    for(int kk = 0; kk < currShell->noOfContr; kk++)
      currShell->coeffList[kk] *= normalizationFactor;

  } /* END FOR i each shell */
  
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "BasisInfoStruct::normalizeShells finished, normalizationFactor min max : %15.11f %15.11f",
	    (double)normFactorTot_min, (double)normFactorTot_max);
  
  return 0;
}






/** creates list of 'basis functions', and set startIndexInMatrix for
    each shell. */
int BasisInfoStruct::get_basis_funcs()
{
  int nShells = this->noOfShells;
  int count = 0;
  for(int i = 0; i < nShells; i++) {
    ShellSpecStruct* currShell = &this->shellList[i];
    currShell->startIndexInMatrix = count;
    count += currShell->noOfBasisFuncs;
  }
  this->noOfBasisFuncs = count;
  this->basisFuncList = new BasisFuncStruct[count];
  count = 0;
  for(int i = 0; i < nShells; i++) {
    ShellSpecStruct* currShell = &this->shellList[i];
    int nFunctions = currShell->noOfBasisFuncs;
    for(int j = 0; j < nFunctions; j++) {
      this->basisFuncList[count].noOfContr = currShell->noOfContr;
      for(int kk = 0; kk < currShell->noOfContr; kk++) {
	this->basisFuncList[count].coeffList[kk] = currShell->coeffList[kk];
	this->basisFuncList[count].exponentList[kk] = 
	  currShell->exponentList[kk];
      } /* END FOR kk */
      for(int kk = 0; kk < 3; kk++)
	this->basisFuncList[count].centerCoords[kk] = 
	  currShell->centerCoords[kk];
      this->basisFuncList[count].shellType = currShell->shellType;
      this->basisFuncList[count].functionNumber = j;
      count++;
    } /* END FOR j */
  } /* END FOR i each shell */
  if(count != this->noOfBasisFuncs) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_basis_funcs: "
	      "(count != this->noOfBasisFuncs)");
    return -1;
  }
  return 0;
}


int BasisInfoStruct::getSimplePrimitivesAll(const IntegralInfo& integralInfo)
{
  BasisFuncStruct* currBasisFunc;
  int nbast = this->noOfBasisFuncs;
  const int MAX_NO_OF_PRIMITIVES_PER_BASIS_FUNC = 44; // Guess here to be able to preallocate large enough buffer
  int maxNoOfSimplePrimsTot = nbast * MAX_NO_OF_PRIMITIVES_PER_BASIS_FUNC;
  DistributionSpecStruct* list = new DistributionSpecStruct[maxNoOfSimplePrimsTot];
  
  /* create list of 'simple primitives' */
  int n = 0;
  for(int i = 0; i < nbast; i++) {
    currBasisFunc = &basisFuncList[i];
    int noOfPrimitives = get_simple_primitives(currBasisFunc,
					       list,
					       n,
					       maxNoOfSimplePrimsTot,
					       integralInfo,
                                               use_6_d_funcs);
    if(noOfPrimitives <= 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_simple_primitives");
      return -1;
    }
    currBasisFunc->noOfSimplePrimitives = noOfPrimitives;
    currBasisFunc->simplePrimitiveIndex = n;
    n += noOfPrimitives;
  } /* END FOR i */
  if(simplePrimitiveList) {
    printf("Releasing old simple primitive list\n");
    delete []simplePrimitiveList;
  }
  this->simplePrimitiveList = new DistributionSpecStruct[n];
  memcpy(this->simplePrimitiveList, list, n * sizeof(DistributionSpecStruct));
  delete [] list;
  this->noOfSimplePrimitives = n;
  return 0;
}

/** Initializes all the fields to sane values. */
BasisInfoStruct::BasisInfoStruct(int use_6_d_funcs_) : 
  use_6_d_funcs(use_6_d_funcs_),
  noOfShells(0),
  shellList(NULL),
  noOfBasisFuncs(0),
  basisFuncList(NULL),
  noOfSimplePrimitives(0),
  simplePrimitiveList(NULL)
{
}

/** Copies values from another BasisInfoStruct. */
BasisInfoStruct::BasisInfoStruct(const BasisInfoStruct & b) : 
  use_6_d_funcs(b.use_6_d_funcs),
  noOfShells(b.noOfShells),
  noOfBasisFuncs(b.noOfBasisFuncs),
  noOfSimplePrimitives(b.noOfSimplePrimitives)
{
  shellList = new ShellSpecStruct[noOfShells];
  memcpy(shellList, b.shellList, noOfShells*sizeof(ShellSpecStruct));
  basisFuncList = new BasisFuncStruct[noOfBasisFuncs];
  memcpy(basisFuncList, b.basisFuncList, noOfBasisFuncs*sizeof(BasisFuncStruct));
  simplePrimitiveList = new DistributionSpecStruct[noOfSimplePrimitives];
  memcpy(simplePrimitiveList, b.simplePrimitiveList, noOfSimplePrimitives*sizeof(DistributionSpecStruct));
}


/** Function needed for Chunks&Tasks usage. */
void BasisInfoStruct::write_to_buffer ( char * dataBuffer, size_t const bufferSize ) const {
  // First store the 4 int numbers.
  char* p = dataBuffer;
  if(bufferSize < get_size())
    throw std::runtime_error("Error in BasisInfoStruct::write_to_buffer: bufferSize too small.");
  memcpy(p, &noOfShells, sizeof(int));
  p += sizeof(int);
  memcpy(p, &noOfBasisFuncs, sizeof(int));
  p += sizeof(int);
  memcpy(p, &noOfSimplePrimitives, sizeof(int));
  p += sizeof(int);
  memcpy(p, &use_6_d_funcs, sizeof(int));
  p += sizeof(int);  
  // There are three lists that need to be stored. Take care of them one by one.
  // shellList
  memcpy(p, shellList, noOfShells * sizeof(ShellSpecStruct));
  p += noOfShells * sizeof(ShellSpecStruct);
  // basisFuncList
  memcpy(p, basisFuncList, noOfBasisFuncs * sizeof(BasisFuncStruct));
  p += noOfBasisFuncs * sizeof(BasisFuncStruct);
  // simplePrimitiveList
  memcpy(p, simplePrimitiveList, noOfSimplePrimitives * sizeof(DistributionSpecStruct));
  p += noOfSimplePrimitives * sizeof(DistributionSpecStruct);
  // DONE!
}

/** Function needed for Chunks&Tasks usage. */
size_t BasisInfoStruct::get_size() const {
  return 4 * sizeof(int) 
    + noOfShells * sizeof(ShellSpecStruct)
    + noOfBasisFuncs * sizeof(BasisFuncStruct)
    + noOfSimplePrimitives * sizeof(DistributionSpecStruct);
}

/** Function needed for Chunks&Tasks usage. */
void BasisInfoStruct::assign_from_buffer ( char const * dataBuffer, size_t const bufferSize) {
  // First get the 4 int numbers.
  const char* p = dataBuffer;
  if(bufferSize < 4 * sizeof(int))
    throw std::runtime_error("Error in BasisInfoStruct::assign_from_buffer: bufferSize too small.");
  memcpy(&noOfShells, p, sizeof(int));
  p += sizeof(int);
  memcpy(&noOfBasisFuncs, p, sizeof(int));
  p += sizeof(int);
  memcpy(&noOfSimplePrimitives, p, sizeof(int));
  p += sizeof(int);
  memcpy(&use_6_d_funcs, p, sizeof(int));
  p += sizeof(int);  
  // There are three lists that need to be set up. Take care of them one by one.
  // shellList
  shellList = new ShellSpecStruct[noOfShells];
  memcpy(shellList, p, noOfShells * sizeof(ShellSpecStruct));
  p += noOfShells * sizeof(ShellSpecStruct);
  // basisFuncList
  basisFuncList = new BasisFuncStruct[noOfBasisFuncs];
  memcpy(basisFuncList, p, noOfBasisFuncs * sizeof(BasisFuncStruct));
  p += noOfBasisFuncs * sizeof(BasisFuncStruct);
  // simplePrimitiveList
  simplePrimitiveList = new DistributionSpecStruct[noOfSimplePrimitives];
  memcpy(simplePrimitiveList, p, noOfSimplePrimitives * sizeof(DistributionSpecStruct));
  p += noOfSimplePrimitives * sizeof(DistributionSpecStruct);
  // DONE!
  if(static_cast<size_t>(p-dataBuffer) > bufferSize)
    throw std::runtime_error("Error: (p > bufferSize).");
}

BasisInfoStruct::~BasisInfoStruct()
{
  if(shellList)                 delete [] shellList;
  if(basisFuncList)             delete [] basisFuncList;
  if(simplePrimitiveList)       delete [] simplePrimitiveList;
}


/** Compute safe upper limit for largest possible distance between any
    two basis functions in given basis set.
*/
ergo_real getSafeMaxDistance(const BasisInfoStruct & basisInfo)
{
  ergo_real minCoords[3];
  ergo_real maxCoords[3];
  for(int coordNo = 0; coordNo < 3; coordNo++)
    {
      minCoords[coordNo] = basisInfo.basisFuncList[0].centerCoords[coordNo];
      maxCoords[coordNo] = basisInfo.basisFuncList[0].centerCoords[coordNo];
    }
  for(int i = 0; i < basisInfo.noOfBasisFuncs; i++)
    {
      for(int coordNo = 0; coordNo < 3; coordNo++)
	{
	  ergo_real curr = basisInfo.basisFuncList[i].centerCoords[coordNo];
	  if(curr < minCoords[coordNo])
	    minCoords[coordNo] = curr;
	  if(curr > maxCoords[coordNo])
	    maxCoords[coordNo] = curr;
	}
    }
  ergo_real sum = 0;
  for(int coordNo = 0; coordNo < 3; coordNo++)
    {
      ergo_real dx = maxCoords[coordNo] - minCoords[coordNo];
      sum += dx*dx;
    }
  return template_blas_sqrt(sum);
}
