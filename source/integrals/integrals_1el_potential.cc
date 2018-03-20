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

/** @file integrals_1el_potential.cc

    @brief Code for 1-electron integrals, computation of
    electron-nuclear potential energy matrix V.

    @author: Elias Rudberg <em>responsible</em>
*/

/* Written by Elias Rudberg. First version was written in 2006 or even
   earlier, at KTH. Then the code has been updated many times after
   that.  */

/* 
   ELIAS NOTE 2014-05-31: gradient computation has been added now. The
   code for computing the gradient uses that the nuclear coordinates
   enter the coputation of the V matrix in two ways: 

   (1) Directly in integrals of the type (ij|A) where ij is a basis
   function pair and A is an atom index. In this case we can get the
   derivatives directly using get_related_integrals_hermite().

   (2) Via the multipole tree. For this kind of contribution we can
   compute the derivatives indirectly, by first computing derivatives
   of th emultipole moments w.r.t. nuclear coordinates, and
   derivatives of the energy w.r.t. the multipole moments. Then, we
   multiply those values to get the final derivatives of the energy
   w.r.t. nuclear coordinates. I.e., we compute the derivative df/dx
   using the chain rule, df/dx = (df/dy)*(dy/dx).
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <memory.h>
#include <time.h>
#include <stdarg.h>
#include <assert.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "integrals_1el.h"
#include "integrals_1el_potential.h"
#include "integrals_1el_potential_prep.h"
#include "memorymanag.h"
#include "pi.h"
#include "output.h"
#include "utilities.h"
#include "boysfunction.h"
#include "integral_info.h"
#include "integrals_general.h"
#include "box_system.h"
#include "multipole.h"
#include "integrals_2el_single.h"
#include "integrals_1el_single.h"
#include "integrals_hermite.h"
#include "matrix_norm.h"
#include "mm_limit_table.h"



typedef struct
{
  box_struct_basic basicBox;
  ergo_real centerOfChargeCoords[3];
  multipole_struct_large multipole;
  ergo_real* multipole_moment_derivatives; // If needed, list of derivatives of multipole moments w.r.t. atom coordinates for all atoms in this box.
  ergo_real* derivatives_wrt_multipole_moments;  // If needed, list of derivatives of the energy w.r.t. the multipole moments for this box.
} atom_box_struct;


static ergo_real
get_distance_3d(const ergo_real* x, const ergo_real* y)
{
  ergo_real sum = 0;
  for(int k = 0; k < 3; k++)
    {
      ergo_real dx = x[k] - y[k];
      sum += dx * dx;
    }
  return template_blas_sqrt(sum);
}


static void get_multipole_contribs_for_atom(multipole_struct_large & boxMultipole, ergo_real* multipolePointCoords, const Atom & currAtom, const MMTranslator & translator) {
  multipole_struct_small multipoleCurrAtom;
  memset(&multipoleCurrAtom, 0, sizeof(multipole_struct_small));
  multipoleCurrAtom.degree = 0;
  multipoleCurrAtom.noOfMoments = 1;
  multipoleCurrAtom.momentList[0] = currAtom.charge;
  ergo_real dx = currAtom.coords[0] - multipolePointCoords[0];
  ergo_real dy = currAtom.coords[1] - multipolePointCoords[1];
  ergo_real dz = currAtom.coords[2] - multipolePointCoords[2];

  ergo_real W[MAX_NO_OF_MOMENTS_PER_MULTIPOLE*MAX_NO_OF_MOMENTS_PER_MULTIPOLE];
  translator.getTranslationMatrix(dx, dy, dz, MAX_MULTIPOLE_DEGREE, 0, W);
		  
  multipole_struct_large translatedMultipole;
  for(int A = 0; A < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; A++)
    {
      ergo_real sum = 0;
      for(int B = 0; B < 1; B++)
	sum += W[A*1+B] * multipoleCurrAtom.momentList[B];
      translatedMultipole.momentList[A] = sum;
    } // END FOR A
  for(int kk = 0; kk < 3; kk++)
    translatedMultipole.centerCoords[kk] = multipolePointCoords[kk];
  translatedMultipole.degree = MAX_MULTIPOLE_DEGREE;
  translatedMultipole.noOfMoments = MAX_NO_OF_MOMENTS_PER_MULTIPOLE;
	  
  // add translated multipole to box multipole
  for(int A = 0; A < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; A++)
    boxMultipole.momentList[A] += translatedMultipole.momentList[A];  
}


static void init_multipole_struct_large(multipole_struct_large & boxMultipole, const ergo_real* multipolePointCoords) {
  memset(&boxMultipole, 0, sizeof(multipole_struct_large));
  for(int kk = 0; kk < 3; kk++)
    boxMultipole.centerCoords[kk] = multipolePointCoords[kk];
  boxMultipole.degree = MAX_MULTIPOLE_DEGREE;
  boxMultipole.noOfMoments = MAX_NO_OF_MOMENTS_PER_MULTIPOLE;
}


static int
create_nuclei_mm_tree(const IntegralInfo& integralInfo,
		      int nAtoms,
		      const Atom* atomList,
		      ergo_real boxSize,
		      BoxSystem & boxSystem,
		      atom_box_struct** return_boxList, // allocated by this routine, must be freed by caller.
		      int *return_numberOfLevels,
		      Atom **return_atomListReordered, // allocated by this routine, must be freed by caller.
		      int* return_atomPermutation, // list of int, length=nAtoms, saying how atoms have been reordered in return_atomListReordered.
		      bool compute_gradient_also
		      )
{
  // Create box system based on atoms
  box_item_struct* itemList = new box_item_struct[nAtoms];
  for(int i = 0; i < nAtoms; i++)
    {
      for(int j = 0; j < 3; j++)
	itemList[i].centerCoords[j] = atomList[i].coords[j];
      itemList[i].originalIndex = i;
    } // END FOR i
  
  const ergo_real maxToplevelBoxSize = boxSize;
  
  if(boxSystem.create_box_system(itemList,
				 nAtoms,
				 maxToplevelBoxSize) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system");
      return -1;
    }

  // Now the itemList has been reordered. Store the ordering info in the return_atomPermutation list.
  for(int i = 0; i < nAtoms; i++)
    return_atomPermutation[i] = itemList[i].originalIndex;

  int numberOfLevels = boxSystem.noOfLevels;
  
  atom_box_struct* boxList = new atom_box_struct[boxSystem.totNoOfBoxes];
  for(int i = 0; i < boxSystem.totNoOfBoxes; i++)
    boxList[i].basicBox = boxSystem.boxList[i];

  // create new list of atoms, where they are ordered box by box at the level of smallest boxes.
  Atom* atomList2 = new Atom[nAtoms];
  for(int i = 0; i < nAtoms; i++)
    atomList2[i] = atomList[itemList[i].originalIndex];

  delete [] itemList;
  itemList = NULL;

  // Find center-of-charge for each box at top level
  atom_box_struct* boxListTopLevel = &boxList[boxSystem.levelList[numberOfLevels-1].startIndexInBoxList];
  int noOfBoxesTopLevel = boxSystem.levelList[numberOfLevels-1].noOfBoxes;
  for(int i = 0; i < noOfBoxesTopLevel; i++)
    {
      atom_box_struct* currBox = &boxListTopLevel[i];
      // Find center-of-charge for this box.
      ergo_real cocSumList[3];
      for(int k = 0; k < 3; k++)
	cocSumList[k] = 0;
      ergo_real chargeSum = 0;
      Atom* atomListCurrBox = &atomList2[currBox->basicBox.firstItemIndex];
      int noOfAtomsCurrBox = currBox->basicBox.noOfItems;
      int nPositiveCharges = 0;
      int nNegativeCharges = 0;
      for(int j = 0; j < noOfAtomsCurrBox; j++)
	{
	  Atom* currAtom = &atomListCurrBox[j];
	  chargeSum += currAtom->charge;
	  if(currAtom->charge > 0)
	    nPositiveCharges++;
	  if(currAtom->charge < 0)
	    nNegativeCharges++;
	  for(int k = 0; k < 3; k++)
	    cocSumList[k] += currAtom->coords[k] * currAtom->charge;
	} // END FOR j Find center-of-charge for this box.
      // We only want to use averaging if all charges are of the same
      // sign, otherwise it may happen that they cancel out and then
      // the "average" has no meaning.
      if(nPositiveCharges == noOfAtomsCurrBox || nNegativeCharges == noOfAtomsCurrBox) {
	// OK, all charges are of the same sign.
	for(int k = 0; k < 3; k++)
	  currBox->centerOfChargeCoords[k] = cocSumList[k] / chargeSum;
      }
      else {
	// All charges are not of the same sign. In this case we just use the center of the box as centerOfChargeCoords.
	for(int k = 0; k < 3; k++)
	  currBox->centerOfChargeCoords[k] = currBox->basicBox.centerCoords[k];	
      }

      if(compute_gradient_also) {
	// Always use center point in this case, to simplify computation of derivatives.
	// FIXME: CHECK IF THIS IS REALLY NEEDED, MAYBE IT WAS REALLY BUGS IN OTHER PLACES THAT CAUSED THE PROBLEMS?
	for(int k = 0; k < 3; k++)
	  currBox->centerOfChargeCoords[k] = currBox->basicBox.centerCoords[k];
      }

    } // END FOR i Find center-of-charge for each box at top level

  // Create multipole for each box at top level (smallest boxes)
  MMTranslator translator(integralInfo.GetMultipolePrep());
  for(int i = 0; i < noOfBoxesTopLevel; i++)
    {
      atom_box_struct* currBox = &boxListTopLevel[i];
      ergo_real* multipolePointCoords = currBox->centerOfChargeCoords;
      int noOfAtomsCurrBox = currBox->basicBox.noOfItems;
      if(compute_gradient_also) {
	currBox->multipole_moment_derivatives = new ergo_real[noOfAtomsCurrBox*MAX_NO_OF_MOMENTS_PER_MULTIPOLE*3];
	currBox->derivatives_wrt_multipole_moments = new ergo_real[MAX_NO_OF_MOMENTS_PER_MULTIPOLE];
	memset(currBox->derivatives_wrt_multipole_moments, 0, MAX_NO_OF_MOMENTS_PER_MULTIPOLE*sizeof(ergo_real));
      }
      multipole_struct_large boxMultipole;
      init_multipole_struct_large(boxMultipole, multipolePointCoords);
      Atom* atomListCurrBox = &atomList2[currBox->basicBox.firstItemIndex];
      // Go through all atoms in this box
      for(int j = 0; j < noOfAtomsCurrBox; j++)
	{
	  const Atom & currAtom = atomListCurrBox[j];
      	  // take multipole for this atom, and translate it to center-of-charge point
	  // the "multipole" for this atom is of course only a monopole.
	  get_multipole_contribs_for_atom(boxMultipole, multipolePointCoords, currAtom, translator);
	  if(compute_gradient_also) {
	    for(int coordIdx = 0; coordIdx < 3; coordIdx++) {
	      multipole_struct_large boxMultipoleTmp1;
	      init_multipole_struct_large(boxMultipoleTmp1, multipolePointCoords);
	      multipole_struct_large boxMultipoleTmp2;
	      init_multipole_struct_large(boxMultipoleTmp2, multipolePointCoords);
	      const ergo_real eps = 1e-5;
	      Atom atomTmp1 = currAtom;
	      atomTmp1.coords[coordIdx] += eps;
	      Atom atomTmp2 = currAtom;
	      atomTmp2.coords[coordIdx] -= eps;
	      get_multipole_contribs_for_atom(boxMultipoleTmp1, multipolePointCoords, atomTmp1, translator);
	      get_multipole_contribs_for_atom(boxMultipoleTmp2, multipolePointCoords, atomTmp2, translator);
	      for(int ii = 0; ii < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; ii++) {
		ergo_real m1 = boxMultipoleTmp1.momentList[ii];
		ergo_real m2 = boxMultipoleTmp2.momentList[ii];
		ergo_real value = (m1 - m2) / (2*eps);
		currBox->multipole_moment_derivatives[j*MAX_NO_OF_MOMENTS_PER_MULTIPOLE*3 + MAX_NO_OF_MOMENTS_PER_MULTIPOLE*coordIdx + ii] = value;
	      }
	    }
	  }
	} // END FOR j Go through all atoms in this box

      currBox->multipole = boxMultipole;
    } // END FOR i Create multipole for each box at top level (smallest boxes)



  // OK, multipoles created for top level.
  // Now go through the other levels, joining multipoles from child boxes to a single multipole in parent box

  for(int levelNumber = numberOfLevels-2; levelNumber >= 0; levelNumber--)
    {
      int noOfBoxesCurrLevel = boxSystem.levelList[levelNumber].noOfBoxes;
      atom_box_struct* boxListCurrLevel = &boxList[boxSystem.levelList[levelNumber].startIndexInBoxList];
      
      for(int boxIndex = 0; boxIndex < noOfBoxesCurrLevel; boxIndex++)
	{
	  atom_box_struct* currBox = &boxListCurrLevel[boxIndex];
	  int noOfChildren = currBox->basicBox.noOfChildBoxes;

	  if(noOfChildren == 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "ERROR: (noOfChildren == 0)");
	      return -1;
	    }

	  multipole_struct_large* newMultipole = &currBox->multipole;
	  memset(newMultipole, 0, sizeof(multipole_struct_large));

	  // get average position of child multipoles
	  ergo_real avgPosList[3];
	  int kk;
	  for(kk = 0; kk < 3; kk++)
	    avgPosList[kk] = 0;
	  int childIndex;
	  for(childIndex = 0; childIndex < noOfChildren; childIndex++)
	    {
	      int childIndexInBoxList = currBox->basicBox.firstChildBoxIndex + childIndex;
	      atom_box_struct* childBox = &boxList[childIndexInBoxList];
	      for(kk = 0; kk < 3; kk++)
		avgPosList[kk] += childBox->multipole.centerCoords[kk];
	    } // END FOR childIndex	  
	  for(kk = 0; kk < 3; kk++)
	    newMultipole->centerCoords[kk] = avgPosList[kk] / noOfChildren;

	  if(compute_gradient_also) {
	    // Always use center point in this case, to simplify computation of derivatives.
	    // FIXME: CHECK IF THIS IS REALLY NEEDED, MAYBE IT WAS REALLY BUGS IN OTHER PLACES THAT CAUSED THE PROBLEMS?
	    for(kk = 0; kk < 3; kk++)
	      newMultipole->centerCoords[kk] = currBox->basicBox.centerCoords[kk];
	  }
	  
	  newMultipole->degree = MAX_MULTIPOLE_DEGREE;
	  newMultipole->noOfMoments = MAX_NO_OF_MOMENTS_PER_MULTIPOLE;

	  // Now translate child multipoles and add to parent multipole
	  for(childIndex = 0; childIndex < noOfChildren; childIndex++)
	    {
	      int childIndexInBoxList = currBox->basicBox.firstChildBoxIndex + childIndex;
	      atom_box_struct* childBox = &boxList[childIndexInBoxList];
	      multipole_struct_large* childMultipole = &childBox->multipole;

	      ergo_real dx = childMultipole->centerCoords[0] - newMultipole->centerCoords[0];
	      ergo_real dy = childMultipole->centerCoords[1] - newMultipole->centerCoords[1];
	      ergo_real dz = childMultipole->centerCoords[2] - newMultipole->centerCoords[2];

	      ergo_real W[MAX_NO_OF_MOMENTS_PER_MULTIPOLE*MAX_NO_OF_MOMENTS_PER_MULTIPOLE];
	      translator.getTranslationMatrix(dx, dy, dz, MAX_MULTIPOLE_DEGREE, MAX_MULTIPOLE_DEGREE, W);

	      multipole_struct_large translatedMultipole;
	      int A, B;
	      for(A = 0; A < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; A++)
		{
		  ergo_real sum = 0;
		  for(B = 0; B < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; B++)
		    sum += W[A*MAX_NO_OF_MOMENTS_PER_MULTIPOLE+B] * childMultipole->momentList[B];
		  translatedMultipole.momentList[A] = sum;
		} // END FOR A
	      for(kk = 0; kk < 3; kk++)
		translatedMultipole.centerCoords[kk] = newMultipole->centerCoords[kk];
	      translatedMultipole.degree = MAX_MULTIPOLE_DEGREE;
	      translatedMultipole.noOfMoments = MAX_NO_OF_MOMENTS_PER_MULTIPOLE;

	      // add translated multipole to parent multipole
	      for(A = 0; A < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; A++)
		newMultipole->momentList[A] += translatedMultipole.momentList[A];
	    } // END FOR childIndex

	  if(compute_gradient_also) {
	    // Get multipole for large box by going through all atoms instead of translating child multipoles
	    ergo_real* multipolePointCoords = newMultipole->centerCoords;
	    multipole_struct_large boxMultipole;
	    init_multipole_struct_large(boxMultipole, multipolePointCoords);
	    Atom* atomListCurrBox = &atomList2[currBox->basicBox.firstItemIndex];
	    int noOfAtomsCurrBox = currBox->basicBox.noOfItems;
	    currBox->multipole_moment_derivatives = new ergo_real[noOfAtomsCurrBox*MAX_NO_OF_MOMENTS_PER_MULTIPOLE*3];
	    currBox->derivatives_wrt_multipole_moments = new ergo_real[MAX_NO_OF_MOMENTS_PER_MULTIPOLE];
	    memset(currBox->derivatives_wrt_multipole_moments, 0, MAX_NO_OF_MOMENTS_PER_MULTIPOLE*sizeof(ergo_real));
	    // Go through all atoms in this box
	    for(int j = 0; j < noOfAtomsCurrBox; j++) {
	      const Atom & currAtom = atomListCurrBox[j];
	      // take multipole for this atom, and translate it to center-of-charge point
	      // the "multipole" for this atom is of course only a monopole.
	      get_multipole_contribs_for_atom(boxMultipole, multipolePointCoords, currAtom, translator);
	      for(int coordIdx = 0; coordIdx < 3; coordIdx++) {
		multipole_struct_large boxMultipoleTmp1;
		init_multipole_struct_large(boxMultipoleTmp1, multipolePointCoords);
		multipole_struct_large boxMultipoleTmp2;
		init_multipole_struct_large(boxMultipoleTmp2, multipolePointCoords);
		const ergo_real eps = 1e-5;
		Atom atomTmp1 = currAtom;
		atomTmp1.coords[coordIdx] += eps;
		Atom atomTmp2 = currAtom;
		atomTmp2.coords[coordIdx] -= eps;
		get_multipole_contribs_for_atom(boxMultipoleTmp1, multipolePointCoords, atomTmp1, translator);
		get_multipole_contribs_for_atom(boxMultipoleTmp2, multipolePointCoords, atomTmp2, translator);
		for(int ii = 0; ii < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; ii++) {
		  ergo_real m1 = boxMultipoleTmp1.momentList[ii];
		  ergo_real m2 = boxMultipoleTmp2.momentList[ii];
		  ergo_real value = (m1 - m2) / (2*eps);
		  currBox->multipole_moment_derivatives[j*MAX_NO_OF_MOMENTS_PER_MULTIPOLE*3 + MAX_NO_OF_MOMENTS_PER_MULTIPOLE*coordIdx + ii] = value;
		}
	      }
	    } // END FOR j Go through all atoms in this box
	  }

	} // END FOR boxIndex
    } // END FOR levelNumber  


  // Prepare info needed later to determine needed multipole degree
  for(int levelNumber = 0; levelNumber < numberOfLevels; levelNumber++)
    {
      int noOfBoxesCurrLevel = boxSystem.levelList[levelNumber].noOfBoxes;
      atom_box_struct* boxListCurrLevel = &boxList[boxSystem.levelList[levelNumber].startIndexInBoxList];      
      for(int boxIndex = 0; boxIndex < noOfBoxesCurrLevel; boxIndex++)
	{
	  atom_box_struct* currBox = &boxListCurrLevel[boxIndex];
	  setup_multipole_maxAbsMomentList(&currBox->multipole);
	}
    }


  *return_boxList = boxList;
  *return_numberOfLevels = numberOfLevels;
  *return_atomListReordered = atomList2;
  return 0;
}



typedef struct {
  DistributionSpecStruct distr;
  int pairIdx;
  int basisFuncIdx1;
  int basisFuncIdx2;
} DistributionSpecStructWithIndexes;


/**
   Take care of interaction between list of distrs and box.
*/
static int
do_interaction_recursive(const IntegralInfo & integralInfo,
			 ergo_real* V_list,
			 int noOfBasisFuncIndexPairs,
			 const basis_func_index_pair_struct_1el* basisFuncIndexPairList,
			 const DistributionSpecStructWithIndexes* list,
			 int nDistrs,
			 const multipole_struct_small* multipoleList,
			 const ergo_real* maxMomentVectorNormForDistrsList,
			 int maxNoOfMomentsForDistrs,
			 int maxDegreeForDistrs,
			 ergo_real distrExtent,
			 const Atom *atomListReordered,
			 const int* atomPermutation, // list of int, length=nAtoms, saying how atoms have been reordered in return_atomListReordered.
			 ergo_real threshold,
			 const atom_box_struct* boxList,
			 MMInteractor & interactor,
			 int boxIndex,
			 int currLevel,
			 int numberOfLevels,
			 bool compute_gradient_also,
			 const ergo_real* D_list,         // used in compute_gradient_also case, NULL otherwise
			 ergo_real* result_gradient_list  // used in compute_gradient_also case, NULL otherwise
			 )
{
  const atom_box_struct* currBox = &boxList[boxIndex];
  const multipole_struct_large* boxMultipole = &currBox->multipole;

  // check if current box is far enough away so that we can use multipole description.
  ergo_real distance = get_distance_3d(list[0].distr.centerCoords, currBox->basicBox.centerCoords);
  ergo_real boxRadius = currBox->basicBox.width * 0.5 * template_blas_sqrt((ergo_real)3);
  ergo_real requiredDistance = boxRadius + distrExtent;

  // Note that the distance to the box multipole is different, since
  // it is not necessarily placed at the box center.
  ergo_real multipoleDistance = get_distance_3d(list[0].distr.centerCoords, boxMultipole->centerCoords);

  int degreeNeeded = integralInfo.GetMMLimitTable().get_minimum_multipole_degree_needed(multipoleDistance,
											boxMultipole,
											maxDegreeForDistrs,
											maxMomentVectorNormForDistrsList,
											threshold);
  if(degreeNeeded < 0)
    return -1;
  degreeNeeded+=2; // We need a couple of extra degrees to handle gradient computation. FIXME: IS THIS REALL NEEDED? WHY? MAYBE IT WAS BUGS IN OTHER PLACES THAT CAUSED THE PROBLEMS?

  bool multipoleDegreeIsSafe = false;
  // Demand at least two degrees margin compared to
  // MAX_MULTIPOLE_DEGREE, otherwise this may fail due to alternating
  // odd/even degrees where only one of them are significant. This has
  // happened in some test cases.
  if(degreeNeeded <= MAX_MULTIPOLE_DEGREE-2)
    multipoleDegreeIsSafe = true;
  
  if(distance > requiredDistance && multipoleDegreeIsSafe)
    {
      // OK, use multipole description of atom charges.
      int boxNeededNoOfMoments = (degreeNeeded+1)*(degreeNeeded+1);
      // create interaction matrix
      ergo_real T[boxNeededNoOfMoments * maxNoOfMomentsForDistrs];
      ergo_real dx = boxMultipole->centerCoords[0] - list[0].distr.centerCoords[0];
      ergo_real dy = boxMultipole->centerCoords[1] - list[0].distr.centerCoords[1];
      ergo_real dz = boxMultipole->centerCoords[2] - list[0].distr.centerCoords[2];

      interactor.getInteractionMatrix(dx, dy, dz, maxDegreeForDistrs, degreeNeeded, T);

      ergo_real tempVector[MAX_NO_OF_MOMENTS_PER_MULTIPOLE];
      for(int A = 0; A < maxNoOfMomentsForDistrs; A++)
	{
	  ergo_real sum = 0;
	  for(int B = 0; B < boxNeededNoOfMoments; B++)
	    {
	      ergo_real momB = boxMultipole->momentList[B];
	      ergo_real Telement = T[A*boxNeededNoOfMoments+B];
	      sum += momB * Telement;
	    }
	  tempVector[A] = sum;
	}

      for(int i = 0; i < nDistrs; i++)
	{
	  ergo_real sum = 0;
	  for(int A = 0; A < multipoleList[i].noOfMoments; A++)
	    sum += tempVector[A] * multipoleList[i].momentList[A];
	  V_list[list[i].pairIdx] += -1 * sum;
	}

      if(compute_gradient_also) {
	// we need to compute the derivatives of the energy with respect to the multipole moments.
	// To reduce number of critical sections when OpenMP threading is used, use temporary derivatives_wrt_multipole_moments list first, and then update the real one.
	ergo_real derivatives_wrt_multipole_moments_tmp[MAX_NO_OF_MOMENTS_PER_MULTIPOLE];
	memset(derivatives_wrt_multipole_moments_tmp, 0, MAX_NO_OF_MOMENTS_PER_MULTIPOLE*sizeof(ergo_real));
	for(int i = 0; i < nDistrs; i++) {
	  ergo_real tempVector2[boxNeededNoOfMoments];
	  for(int B = 0; B < boxNeededNoOfMoments; B++) {
	    ergo_real sum = 0;
	    for(int A = 0; A < multipoleList[i].noOfMoments; A++) {
	      ergo_real Telement = T[A*boxNeededNoOfMoments+B];
	      ergo_real mom = multipoleList[i].momentList[A];
	      sum += mom * Telement;
	    }
	    tempVector2[B] = sum;
	  } // end for B
	  int pairIdx = list[i].pairIdx;
	  ergo_real extraFactor = 1.0;
	  if(basisFuncIndexPairList[pairIdx].index_1 != basisFuncIndexPairList[pairIdx].index_2)
	    extraFactor = 2.0;
	  ergo_real densityMatrixElement = D_list[list[i].pairIdx];
	  for(int B = 0; B < boxNeededNoOfMoments; B++)
	    derivatives_wrt_multipole_moments_tmp[B] += -1 * extraFactor * densityMatrixElement * tempVector2[B];
	} // end for i (loop over distrs)

#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  for(int B = 0; B < boxNeededNoOfMoments; B++)
	    currBox->derivatives_wrt_multipole_moments[B] += derivatives_wrt_multipole_moments_tmp[B];
	}

      } // end if compute_gradient_also

      return 0;
    }

  // No, multipole description could not be used in this case.
  if(currLevel == numberOfLevels-1)
    {
      // We are at top level, must compute explicit interactions

      int Nmax = 0;
      for(int i = 0; i < nDistrs; i++)
	{
	  const DistributionSpecStruct* distr = &list[i].distr;
	  int N = distr->monomialInts[0] + distr->monomialInts[1] + distr->monomialInts[2];
	  if(N > Nmax)
	    Nmax = N;
	}

      const JK::ExchWeights CAM_params_not_used;
      
      const Atom* currAtomList = &atomListReordered[currBox->basicBox.firstItemIndex];
      int noOfAtomsCurrBox = currBox->basicBox.noOfItems;
      for(int ia = 0; ia < noOfAtomsCurrBox; ia++)
	{
	  const Atom* currAtom = &currAtomList[ia];

	  const DistributionSpecStruct* distr1 = &list[0].distr;
	  int noOfMonomials_1 = integralInfo.monomial_info.no_of_monomials_list[Nmax];
	  int noOfMonomials_2 = integralInfo.monomial_info.no_of_monomials_list[0];
	  ergo_real primitiveIntegralList_h[noOfMonomials_1*noOfMonomials_2];
	  ergo_real primitiveIntegralList_2[noOfMonomials_1*noOfMonomials_2];
	  // Let the distr be 1 and the pointcharge 2
	  ergo_real alpha1 = distr1->exponent;
	  ergo_real alpha0 = alpha1;  
	  int n1 = Nmax;
	  int n2 = 0;
	  ergo_real dx0 = currAtom->coords[0] - distr1->centerCoords[0];
	  ergo_real dx1 = currAtom->coords[1] - distr1->centerCoords[1];
	  ergo_real dx2 = currAtom->coords[2] - distr1->centerCoords[2];
	  ergo_real resultPreFactor = 2 * pi / alpha1;
	  get_related_integrals_hermite(integralInfo,
					CAM_params_not_used,
					n1, noOfMonomials_1,
					n2, noOfMonomials_2,
					dx0, 
					dx1, 
					dx2, 
					alpha0,
					resultPreFactor,
					primitiveIntegralList_h);
	  integralInfo.multiply_by_hermite_conversion_matrix_from_right(n1,
									n2,
									1.0/alpha1,
									primitiveIntegralList_h,
									primitiveIntegralList_2);	  
	  for(int id = 0; id < nDistrs; id++)
	    {
	      const DistributionSpecStruct* distr = &list[id].distr;
	      int n1x = distr->monomialInts[0];
	      int n1y = distr->monomialInts[1];
	      int n1z = distr->monomialInts[2];
	      int monomialIndex = integralInfo.monomial_info.monomial_index_list[n1x][n1y][n1z];
	      ergo_real integralValue = currAtom->charge * distr->coeff * primitiveIntegralList_2[monomialIndex];
	      V_list[list[id].pairIdx] += -1 * integralValue;
	    }
	  
	} // END FOR ia

      if(compute_gradient_also) {
	for(int ia = 0; ia < noOfAtomsCurrBox; ia++) {
	  const Atom* currAtom = &currAtomList[ia];
	  const DistributionSpecStruct* distr1 = &list[0].distr;
	  int n1 = Nmax;
	  int n2 = 1;
	  int noOfMonomials_1 = integralInfo.monomial_info.no_of_monomials_list[n1];
	  int noOfMonomials_2 = integralInfo.monomial_info.no_of_monomials_list[n2];
	  ergo_real primitiveIntegralList_h[noOfMonomials_1*noOfMonomials_2];
	  // Let the distr be 1 and the pointcharge 2
	  ergo_real alpha1 = distr1->exponent;
	  ergo_real alpha0 = alpha1;  
	  ergo_real dx0 = currAtom->coords[0] - distr1->centerCoords[0];
	  ergo_real dx1 = currAtom->coords[1] - distr1->centerCoords[1];
	  ergo_real dx2 = currAtom->coords[2] - distr1->centerCoords[2];
	  ergo_real resultPreFactor = 2 * pi / alpha1;
	  get_related_integrals_hermite(integralInfo,
					CAM_params_not_used,
					n1, noOfMonomials_1,
					n2, noOfMonomials_2,
					dx0, 
					dx1, 
					dx2, 
					alpha0,
					resultPreFactor,
					primitiveIntegralList_h);
	  for(int id = 0; id < nDistrs; id++) {
	    const DistributionSpecStruct* distr = &list[id].distr;
	    int n1b = n1;
	    int n2b = 0;
	    ergo_real primitiveIntegralList_h_components[3][noOfMonomials_1];
	    int monomialIndex_x = integralInfo.monomial_info.monomial_index_list[1][0][0];
	    int monomialIndex_y = integralInfo.monomial_info.monomial_index_list[0][1][0];
	    int monomialIndex_z = integralInfo.monomial_info.monomial_index_list[0][0][1];
	    for(int i = 0; i < noOfMonomials_1; i++) {
	      primitiveIntegralList_h_components[0][i] = primitiveIntegralList_h[i*noOfMonomials_2+monomialIndex_x];
	      primitiveIntegralList_h_components[1][i] = primitiveIntegralList_h[i*noOfMonomials_2+monomialIndex_y];
	      primitiveIntegralList_h_components[2][i] = primitiveIntegralList_h[i*noOfMonomials_2+monomialIndex_z];
	    }
	    ergo_real primitiveIntegralList_2_components[3][noOfMonomials_1];
	    for(int i = 0; i < 3; i++)
	      integralInfo.multiply_by_hermite_conversion_matrix_from_right(n1b, n2b, 1.0/alpha1, primitiveIntegralList_h_components[i], primitiveIntegralList_2_components[i]);
	    int n1x = distr->monomialInts[0];
	    int n1y = distr->monomialInts[1];
	    int n1z = distr->monomialInts[2];
	    int monomialIndex = integralInfo.monomial_info.monomial_index_list[n1x][n1y][n1z];
	    ergo_real value_x = currAtom->charge * distr->coeff * primitiveIntegralList_2_components[0][monomialIndex];
	    ergo_real value_y = currAtom->charge * distr->coeff * primitiveIntegralList_2_components[1][monomialIndex];
	    ergo_real value_z = currAtom->charge * distr->coeff * primitiveIntegralList_2_components[2][monomialIndex];
	    int pairIdx = list[id].pairIdx;
	    ergo_real extraFactor = 1.0;
	    if(basisFuncIndexPairList[pairIdx].index_1 != basisFuncIndexPairList[pairIdx].index_2)
	      extraFactor = 2.0;
	    ergo_real densityMatrixElement = D_list[list[id].pairIdx];
	    int atomIndex = atomPermutation[currBox->basicBox.firstItemIndex + ia];
#ifdef _OPENMP
#pragma omp critical
#endif
	    {
	      result_gradient_list[atomIndex*3+0] += -1 * value_x * densityMatrixElement * extraFactor;
	      result_gradient_list[atomIndex*3+1] += -1 * value_y * densityMatrixElement * extraFactor;
	      result_gradient_list[atomIndex*3+2] += -1 * value_z * densityMatrixElement * extraFactor;
	    }
	  }

	}

      }

      return 0;
    }

  // Go to next level
  int noOfChildren = currBox->basicBox.noOfChildBoxes;
  for(int i = 0; i < noOfChildren; i++)
    {
      int childBoxIndex = currBox->basicBox.firstChildBoxIndex + i;
      if(do_interaction_recursive(integralInfo,
				  V_list,
				  noOfBasisFuncIndexPairs,
				  basisFuncIndexPairList,
				  list,
				  nDistrs,
				  multipoleList,
				  maxMomentVectorNormForDistrsList,
				  maxNoOfMomentsForDistrs,
				  maxDegreeForDistrs,
				  distrExtent,
				  atomListReordered,
				  atomPermutation,
				  threshold,
				  boxList,
				  interactor,
				  childBoxIndex,
				  currLevel + 1,
				  numberOfLevels,
				  compute_gradient_also,
				  D_list,
				  result_gradient_list) != 0)
	return -1;
    } // END FOR i 
  return 0;
}



/**
   Take care of interaction between list of distrs and box.
*/
static int
do_interaction_recursive_2(const IntegralInfo & integralInfo,
			   csr_matrix_struct* V_CSR,
			   int noOfBasisFuncIndexPairs,
			   const basis_func_index_pair_struct_1el* basisFuncIndexPairList,
			   const DistributionSpecStructWithIndexes2* list,
			   int nDistrs,
			   const multipole_struct_small* multipoleList,
			   const ergo_real* maxMomentVectorNormForDistrsList,
			   int maxNoOfMomentsForDistrs,
			   int maxDegreeForDistrs,
			   ergo_real distrExtent,
			   const Atom *atomListReordered,
			   const int* atomPermutation, // list of int, length=nAtoms, saying how atoms have been reordered in return_atomListReordered.
			   ergo_real threshold,
			   const atom_box_struct* boxList,
			   MMInteractor & interactor,
			   int boxIndex,
			   int currLevel,
			   int numberOfLevels
			   )
{
  const atom_box_struct* currBox = &boxList[boxIndex];
  const multipole_struct_large* boxMultipole = &currBox->multipole;

  // check if current box is far enough away so that we can use multipole description.
  ergo_real distance = get_distance_3d(list[0].distr.centerCoords, currBox->basicBox.centerCoords);
  ergo_real boxRadius = currBox->basicBox.width * 0.5 * template_blas_sqrt((ergo_real)3);
  ergo_real requiredDistance = boxRadius + distrExtent;

  // Note that the distance to the box multipole is different, since
  // it is not necessarily placed at the box center.
  ergo_real multipoleDistance = get_distance_3d(list[0].distr.centerCoords, boxMultipole->centerCoords);
  int degreeNeeded = integralInfo.GetMMLimitTable().get_minimum_multipole_degree_needed(multipoleDistance,
											boxMultipole,
											maxDegreeForDistrs,
											maxMomentVectorNormForDistrsList,
											threshold);
  if(degreeNeeded < 0)
    return -1;
  degreeNeeded+=2; // We need a couple of extra degrees to handle gradient computation. FIXME: IS THIS REALL NEEDED? WHY? MAYBE IT WAS BUGS IN OTHER PLACES THAT CAUSED THE PROBLEMS?

  bool multipoleDegreeIsSafe = false;
  // Demand at least two degrees margin compared to
  // MAX_MULTIPOLE_DEGREE, otherwise this may fail due to alternating
  // odd/even degrees where only one of them are significant. This has
  // happened in some test cases.
  if(degreeNeeded <= MAX_MULTIPOLE_DEGREE-2)
    multipoleDegreeIsSafe = true;
  
  if(distance > requiredDistance && multipoleDegreeIsSafe)
    {
      // OK, use multipole description of atom charges.
      int boxNeededNoOfMoments = (degreeNeeded+1)*(degreeNeeded+1);
      // create interaction matrix
      ergo_real T[boxNeededNoOfMoments * maxNoOfMomentsForDistrs];
      ergo_real dx = boxMultipole->centerCoords[0] - list[0].distr.centerCoords[0];
      ergo_real dy = boxMultipole->centerCoords[1] - list[0].distr.centerCoords[1];
      ergo_real dz = boxMultipole->centerCoords[2] - list[0].distr.centerCoords[2];

      interactor.getInteractionMatrix(dx, dy, dz, maxDegreeForDistrs, degreeNeeded, T);

      ergo_real tempVector[MAX_NO_OF_MOMENTS_PER_MULTIPOLE];
      for(int A = 0; A < maxNoOfMomentsForDistrs; A++)
	{
	  ergo_real sum = 0;
	  for(int B = 0; B < boxNeededNoOfMoments; B++)
	    {
	      ergo_real momB = boxMultipole->momentList[B];
	      ergo_real Telement = T[A*boxNeededNoOfMoments+B];
	      sum += momB * Telement;
	    }
	  tempVector[A] = sum;
	}

      for(int i = 0; i < nDistrs; i++)
	{
	  ergo_real sum = 0;
	  for(int A = 0; A < multipoleList[i].noOfMoments; A++)
	    sum += tempVector[A] * multipoleList[i].momentList[A];
	  int idx1 = list[i].basisFuncIdx1;
	  int idx2 = list[i].basisFuncIdx2;
	  if(ergo_CSR_add_to_element(V_CSR, idx1, idx2,  -1 * sum) != 0)
	    return -1;
	}

      return 0;
    }

  // No, multipole description could not be used in this case.
  if(currLevel == numberOfLevels-1)
    {
      // We are at top level, must compute explicit interactions

      int Nmax = 0;
      for(int i = 0; i < nDistrs; i++)
	{
	  const DistributionSpecStruct* distr = &list[i].distr;
	  int N = distr->monomialInts[0] + distr->monomialInts[1] + distr->monomialInts[2];
	  if(N > Nmax)
	    Nmax = N;
	}

      const JK::ExchWeights CAM_params_not_used;
      
      const Atom* currAtomList = &atomListReordered[currBox->basicBox.firstItemIndex];
      int noOfAtomsCurrBox = currBox->basicBox.noOfItems;
      for(int ia = 0; ia < noOfAtomsCurrBox; ia++)
	{
	  const Atom* currAtom = &currAtomList[ia];

	  const DistributionSpecStruct* distr1 = &list[0].distr;
	  int noOfMonomials_1 = integralInfo.monomial_info.no_of_monomials_list[Nmax];
	  int noOfMonomials_2 = integralInfo.monomial_info.no_of_monomials_list[0];
	  ergo_real primitiveIntegralList_h[noOfMonomials_1*noOfMonomials_2];
	  ergo_real primitiveIntegralList_2[noOfMonomials_1*noOfMonomials_2];
	  // Let the distr be 1 and the pointcharge 2
	  ergo_real alpha1 = distr1->exponent;
	  ergo_real alpha0 = alpha1;  
	  int n1 = Nmax;
	  int n2 = 0;
	  ergo_real dx0 = currAtom->coords[0] - distr1->centerCoords[0];
	  ergo_real dx1 = currAtom->coords[1] - distr1->centerCoords[1];
	  ergo_real dx2 = currAtom->coords[2] - distr1->centerCoords[2];
	  ergo_real resultPreFactor = 2 * pi / alpha1;
	  get_related_integrals_hermite(integralInfo,
					CAM_params_not_used,
					n1, noOfMonomials_1,
					n2, noOfMonomials_2,
					dx0, 
					dx1, 
					dx2, 
					alpha0,
					resultPreFactor,
					primitiveIntegralList_h);
	  integralInfo.multiply_by_hermite_conversion_matrix_from_right(n1,
									n2,
									1.0/alpha1,
									primitiveIntegralList_h,
									primitiveIntegralList_2);	  
	  for(int id = 0; id < nDistrs; id++)
	    {
	      const DistributionSpecStruct* distr = &list[id].distr;
	      int n1x = distr->monomialInts[0];
	      int n1y = distr->monomialInts[1];
	      int n1z = distr->monomialInts[2];
	      int monomialIndex = integralInfo.monomial_info.monomial_index_list[n1x][n1y][n1z];
	      ergo_real integralValue = currAtom->charge * distr->coeff * primitiveIntegralList_2[monomialIndex];
	      int idx1 = list[id].basisFuncIdx1;
	      int idx2 = list[id].basisFuncIdx2;
	      if(ergo_CSR_add_to_element(V_CSR, idx1, idx2, -1 * integralValue) != 0)
		return -1;
	    }
	  
	} // END FOR ia

      return 0;
    }

  // Go to next level
  int noOfChildren = currBox->basicBox.noOfChildBoxes;
  for(int i = 0; i < noOfChildren; i++)
    {
      int childBoxIndex = currBox->basicBox.firstChildBoxIndex + i;
      if(do_interaction_recursive_2(integralInfo,
				    V_CSR,
				    noOfBasisFuncIndexPairs,
				    basisFuncIndexPairList,
				    list,
				    nDistrs,
				    multipoleList,
				    maxMomentVectorNormForDistrsList,
				    maxNoOfMomentsForDistrs,
				    maxDegreeForDistrs,
				    distrExtent,
				    atomListReordered,
				    atomPermutation,
				    threshold,
				    boxList,
				    interactor,
				    childBoxIndex,
				    currLevel + 1,
				    numberOfLevels) != 0)
	return -1;
    } // END FOR i
  return 0;
}




static int 
get_list_of_distrs_for_V(const BasisInfoStruct& basisInfo,
			 const basis_func_index_pair_struct_1el* basisFuncIndexPairList,
			 int noOfBasisFuncIndexPairs,
			 ergo_real threshold,
			 ergo_real maxCharge,
			 DistributionSpecStructWithIndexes* resultList,
			 int maxCountResult)
{
  int distrCount = 0;
  for(int kk = 0; kk < noOfBasisFuncIndexPairs; kk++)
    {
      int i = basisFuncIndexPairList[kk].index_1;
      int j = basisFuncIndexPairList[kk].index_2;      
      const int maxCountProduct = POLY_PRODUCT_MAX_DISTRS;
      DistributionSpecStruct psi_list[maxCountProduct];
      /* form product of basisfuncs i and j, store product in psi_list */
      int n_psi = get_product_simple_primitives(basisInfo, i,
						basisInfo, j,
						psi_list, maxCountProduct, 0);
      if(n_psi < 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_product_simple_primitives");
	  return -1;
	}
      for(int k = 0; k < n_psi; k++)
	{
	  // now take care of psi_list[k]
	  // Here, we estimate the largest possible contribution to V from this distr as (maxCharge * 2 * pi * coeff / exponent).
	  DistributionSpecStruct* prim = &psi_list[k];
	  ergo_real maxContrib = template_blas_fabs(maxCharge * 2 * pi * prim->coeff / prim->exponent);
	  if(maxContrib > threshold)
	    {
	      if(maxCountResult > 0 && distrCount >= maxCountResult)
		{
		  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_list_of_distrs_for_V: (maxCountResult > 0 && distrCount >= maxCountResult).");
		  return -1;
		}
	      if(resultList != NULL)
		{
		  resultList[distrCount].distr = *prim;
		  resultList[distrCount].pairIdx = kk;
		  resultList[distrCount].basisFuncIdx1 = i;
		  resultList[distrCount].basisFuncIdx2 = j;
		}
	      distrCount++;
	    } // END IF above threshold
	} // END FOR k
    } // END FOR kk
  return distrCount;
}



static ergo_real
get_nucl_repulsion_energy_using_multipoles(const Atom *atomListReordered,
					   const int* atomPermutation, // list of int, length=nAtoms, saying how atoms have been reordered in return_atomListReordered.
					   ergo_real threshold,
					   const atom_box_struct* boxList,
					   MMInteractor & interactor,
					   int boxIndex1,
					   int boxIndex2,
					   int currLevel,
					   int numberOfLevels) {
  const atom_box_struct* currBox1 = &boxList[boxIndex1];
  const atom_box_struct* currBox2 = &boxList[boxIndex2];
  const multipole_struct_large* boxMultipole1 = &currBox1->multipole;
  const multipole_struct_large* boxMultipole2 = &currBox2->multipole;
  // Check if we are at level of smallest boxes
  if(currLevel == numberOfLevels-1) {
    // Do interaction explicitly
    const Atom* currAtomList1 = &atomListReordered[currBox1->basicBox.firstItemIndex];
    int noOfAtomsCurrBox1 = currBox1->basicBox.noOfItems;
    const Atom* currAtomList2 = &atomListReordered[currBox2->basicBox.firstItemIndex];
    int noOfAtomsCurrBox2 = currBox2->basicBox.noOfItems;
    ergo_real sum = 0;
    for(int i1 = 0; i1 < noOfAtomsCurrBox1; i1++)
      for(int i2 = 0; i2 < noOfAtomsCurrBox2; i2++) {
	if(boxIndex1 == boxIndex2 && i1 == i2)
	  continue; // skip interaction of an atom with itself
	if(boxIndex1 == boxIndex2 && i1 > i2)
	  continue; // do not double-count interactions
	const Atom* currAtom1 = &currAtomList1[i1];
	const Atom* currAtom2 = &currAtomList2[i2];
	ergo_real dx0 = currAtom1->coords[0] - currAtom2->coords[0];
	ergo_real dx1 = currAtom1->coords[1] - currAtom2->coords[1];
	ergo_real dx2 = currAtom1->coords[2] - currAtom2->coords[2];
	ergo_real distance = template_blas_sqrt(dx0*dx0+dx1*dx1+dx2*dx2);
	sum += currAtom1->charge * currAtom2->charge / distance;
      }
    return sum;
  }
  // Check if multipole representation can be used
  // check if boxes are far enough apart so that we can consider using multipole description.
  ergo_real distance_between_box_centers = get_distance_3d(currBox1->basicBox.centerCoords, currBox2->basicBox.centerCoords);
  ergo_real boxRadius = currBox1->basicBox.width * 0.5 * template_blas_sqrt((ergo_real)3);
  ergo_real requiredDistance = 2.2*boxRadius; // 2*boxRadius should be enough, use 2.2*boxRadius to have some margin
  if(distance_between_box_centers > requiredDistance) {
    // Try using multipole description of atom charges.
    // Try with different multipole degree and check the difference, if nearly same result is obtained we trust it, otherwise not.
    const int nDegreesToTry = 3;
    ergo_real resultMin = 0;
    ergo_real resultMax = 0;
    const int highest_degree_to_use = MAX_MULTIPOLE_DEGREE / 2; // different choices possible here, pick something that gives good performance
    ergo_real dx = boxMultipole2->centerCoords[0] - boxMultipole1->centerCoords[0];
    ergo_real dy = boxMultipole2->centerCoords[1] - boxMultipole1->centerCoords[1];
    ergo_real dz = boxMultipole2->centerCoords[2] - boxMultipole1->centerCoords[2];
    int boxNeededNoOfMomentsMax = (highest_degree_to_use+1)*(highest_degree_to_use+1);
    ergo_real T[boxNeededNoOfMomentsMax * boxNeededNoOfMomentsMax];
    interactor.getInteractionMatrix(dx, dy, dz, highest_degree_to_use, highest_degree_to_use, T);
    for(int i = 0; i < nDegreesToTry; i++) {
      int degree = highest_degree_to_use - i;
      if(degree < 0)
	throw std::runtime_error("Error in get_nucl_repulsion_energy_using_multipoles: (degree < 0).");
      int boxNeededNoOfMoments = (degree+1)*(degree+1);
      // create interaction matrix
      ergo_real sum = 0;
      for(int A = 0; A < boxNeededNoOfMoments; A++)
	for(int B = 0; B < boxNeededNoOfMoments; B++)
	  sum += T[A*boxNeededNoOfMomentsMax+B] * boxMultipole1->momentList[A] * boxMultipole2->momentList[B];
      if(i == 0) {
	resultMin = sum;
	resultMax = sum;
      }
      if(sum < resultMin)
	resultMin = sum;
      if(sum > resultMax)
	resultMax = sum;
    } // end for i
    ergo_real diff = resultMax - resultMin;
    if(diff < threshold)
      return 0.5*(resultMin+resultMax);
  } // end if distance large enough to consider multipoles
  // Go to next level, smaller boxes
  ergo_real sum = 0;
  int noOfChildren1 = currBox1->basicBox.noOfChildBoxes;
  int noOfChildren2 = currBox2->basicBox.noOfChildBoxes;
  for(int i1 = 0; i1 < noOfChildren1; i1++)
    for(int i2 = 0; i2 < noOfChildren2; i2++) {
      if(boxIndex1 == boxIndex2 && i1 > i2)
	continue; // do not double-count interactions
      int childBoxIndex1 = currBox1->basicBox.firstChildBoxIndex + i1;
      int childBoxIndex2 = currBox2->basicBox.firstChildBoxIndex + i2;
      sum += get_nucl_repulsion_energy_using_multipoles(atomListReordered,
							atomPermutation, // list of int, length=nAtoms, saying how atoms have been reordered in return_atomListReordered.
							threshold,
							boxList,
							interactor,
							childBoxIndex1,
							childBoxIndex2,
							currLevel + 1,
							numberOfLevels);
    } // END FOR i1 i2
  return sum;
}


int compute_V_and_gradient_linear(const BasisInfoStruct& basisInfo,
				  const IntegralInfo& integralInfo,
				  const Molecule& molecule,
				  ergo_real threshold,
				  ergo_real boxSize,
				  const basis_func_index_pair_struct_1el* basisFuncIndexPairList,
				  ergo_real* V_list,
				  int noOfBasisFuncIndexPairs,
				  bool compute_gradient_also,
				  const ergo_real* D_list, // List of corresponding density matrix elemets; used for compute_gradient_also case, NULL otherwise
				  ergo_real* gradient_list, // list of result gradient values; used for compute_gradient_also case, NULL otherwise
				  ergo_real & result_nuclearRepulsionEnergy
				  )
{
  result_nuclearRepulsionEnergy = 0; // to be set later
  int errorCount = 0;
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "compute_V_and_gradient_linear start, compute_gradient_also = %d.", compute_gradient_also);

  Util::TimeMeter timeMeterTotal;
  Util::TimeMeter timeMeterInitPart;

  // Get maxCharge
  ergo_real maxCharge = 0;
  for(int i = 0; i < molecule.getNoOfAtoms(); i++) {
    ergo_real currCharge = molecule.getAtom(i).charge;
    if(currCharge > maxCharge)
      maxCharge = currCharge;
  }

  // Create list of distributions
  int nDistrs = get_list_of_distrs_for_V(basisInfo,
					 basisFuncIndexPairList,
					 noOfBasisFuncIndexPairs,
					 threshold,
					 maxCharge,
					 NULL,
					 0);
  if(nDistrs <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_V_and_gradient_linear: (nDistrs <= 0).");
      return -1;
    }
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "compute_V_and_gradient_linear nDistrs    = %9i", nDistrs);
  std::vector<DistributionSpecStructWithIndexes> list(nDistrs);
  if(get_list_of_distrs_for_V(basisInfo,
			      basisFuncIndexPairList,
			      noOfBasisFuncIndexPairs,
			      threshold,
			      maxCharge,
			      &list[0],
			      nDistrs) != nDistrs)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_V_and_gradient_linear, in get_list_of_distrs_for_V.");
      return -1;
    }

  // Sort list of distrs by x, y, z, exponent.
  // The point of this is to group together distrs that have same center and same exponent.
  sort_distr_list(&list[0], nDistrs);

  // identify groups of distrs that have same center and same exponent.
  // Allocate according to worst case, each distr being a separate group.
  std::vector<group_struct> groupList(nDistrs);
  int ind = 0;
  int currGroupInd = 0;
  int groupCount = 0;
  int maxNDistrsPerGroup = 0;
  while(ind < nDistrs)
    {
      ind++;
      if(ind < nDistrs)
	{
	  if(compare_distrs<DistributionSpecStructWithIndexes>(&list[ind], &list[currGroupInd]) == 0)
	    continue;
	}
      // define new group
      groupList[groupCount].startIndex = currGroupInd;
      groupList[groupCount].count = ind - currGroupInd;
      if (groupList[groupCount].count > maxNDistrsPerGroup)
	maxNDistrsPerGroup = groupList[groupCount].count;
      groupCount++;
      // start next group
      currGroupInd = ind;
    }
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "compute_V_and_gradient_linear groupCount = %9i", groupCount);

  // Note that boxList and atomListReordered are allocated by create_nuclei_mm_tree,
  // we must remember to free them in the end.
  BoxSystem boxSystem;
  std::vector<int> atomPermutation(molecule.getNoOfAtoms());
  atom_box_struct* boxList = NULL;
  Atom *atomListReordered = NULL;
  int numberOfLevels = -1;
  if(create_nuclei_mm_tree(integralInfo,
			   molecule.getNoOfAtoms(), molecule.getAtomListPtr(), boxSize,
			   boxSystem,
			   &boxList, &numberOfLevels,
			   &atomListReordered,
			   &atomPermutation[0],
			   compute_gradient_also) != 0) 
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "create_nuclei_mm_tree failed");
      return -1;
    }

  timeMeterInitPart.print(LOG_AREA_INTEGRALS, "compute_V_and_gradient_linear init part (including create_nuclei_mm_tree)");

  // Compute nuclear repulsion energy using multipole tree.
  Util::TimeMeter timeMeterNuclRep;
  {
    MMInteractor interactor(integralInfo.GetMultipolePrep());
    int boxIndex1 = 0;
    int boxIndex2 = 0;
    int currLevel = 0;
    ergo_real nuclRepEnergy = get_nucl_repulsion_energy_using_multipoles(atomListReordered,
									 &atomPermutation[0],
									 threshold,
									 boxList,
									 interactor,
									 boxIndex1,
									 boxIndex2,
									 currLevel,
									 numberOfLevels);
    do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "get_nucl_repulsion_energy_using_multipoles() gave nuclRepEnergy = %22.11f", nuclRepEnergy);
    result_nuclearRepulsionEnergy = nuclRepEnergy;
  }
  timeMeterNuclRep.print(LOG_AREA_INTEGRALS, "get_nucl_repulsion_energy_using_multipoles");

  Util::TimeMeter timeMeterMainPart;

  const JK::ExchWeights CAM_params_not_used;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
  MMInteractor interactor(integralInfo.GetMultipolePrep());
  ergo_real *private_V_list = V_list;
  ergo_real* private_gradient_list = gradient_list;
  multipole_struct_small* multipoleList =
    new multipole_struct_small[maxNDistrsPerGroup];
#ifdef _OPENMP
  if (omp_get_thread_num() != 0) {
    private_V_list = new ergo_real[noOfBasisFuncIndexPairs];
    if(compute_gradient_also)
      private_gradient_list = new ergo_real[3*molecule.getNoOfAtoms()];
  }
#endif

  memset(private_V_list, 0, noOfBasisFuncIndexPairs*sizeof(ergo_real));
  if(compute_gradient_also)
    memset(private_gradient_list, 0, 3*molecule.getNoOfAtoms()*sizeof(ergo_real));
#ifdef _OPENMP
#pragma omp for
#endif
  for(int groupIndex = 0; groupIndex < groupCount; groupIndex++)
    {
      DistributionSpecStructWithIndexes* currList = &list[groupList[groupIndex].startIndex];
      int nDistrsCurrGroup = groupList[groupIndex].count;

      // Create multipoles for distrs in this group.
      memset(multipoleList, 0, nDistrsCurrGroup*sizeof(multipole_struct_small));
      for(int i = 0; i < nDistrsCurrGroup; i++)
	{
	  compute_multipole_moments(integralInfo,
				    &currList[i].distr,
				    &multipoleList[i]);
	}

      int maxNoOfMoments = 0;
      int maxDegree = 0;
      ergo_real maxMomentVectorNormForDistrsList[MAX_MULTIPOLE_DEGREE_BASIC+1];
      for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
	maxMomentVectorNormForDistrsList[l] = 0;
      for(int i = 0; i < nDistrsCurrGroup; i++)
	{
	  if(multipoleList[i].noOfMoments > maxNoOfMoments)
	    maxNoOfMoments = multipoleList[i].noOfMoments;
	  if(multipoleList[i].degree > maxDegree)
	    maxDegree = multipoleList[i].degree;
	  const multipole_struct_small* distrMultipole = &multipoleList[i];
	  for(int l = 0; l <= distrMultipole->degree; l++)
	    {
	      int startIndex = l*l;
	      int endIndex = (l+1)*(l+1);
	      ergo_real sum = 0;
	      for(int A = startIndex; A < endIndex; A++)
		sum += distrMultipole->momentList[A]*distrMultipole->momentList[A];
	      ergo_real subNorm = template_blas_sqrt(sum);
	      if(subNorm > maxMomentVectorNormForDistrsList[l])
		maxMomentVectorNormForDistrsList[l] = subNorm;
	    }
	}

      // Here we use an extent such that beyond the extent the abs
      // value of any distr is smaller than threshold/maxCharge.
      ergo_real maxabscoeff = 0;
      for(int i = 0; i < nDistrsCurrGroup; i++)
        {
          ergo_real abscoeff = template_blas_fabs(currList[i].distr.coeff);
          if(abscoeff > maxabscoeff)
            maxabscoeff = abscoeff;
        }
      ergo_real exponent = currList[0].distr.exponent; // all exponents in group are equal anyway.
      ergo_real R2 = -1 * (1/exponent) * template_blas_log(threshold/(maxabscoeff*maxCharge));
      ergo_real extent = 0;
      if(R2 > 0) // R2 can become negative, e.g. if maxabscoeff is very small, in such cases we let extent be zero.
	extent = template_blas_sqrt(R2);
      // Take care of interaction of this group with MM tree
      if(do_interaction_recursive(integralInfo,
				  private_V_list, 
				  noOfBasisFuncIndexPairs,
				  basisFuncIndexPairList,
				  currList,
				  nDistrsCurrGroup,
				  multipoleList,
				  maxMomentVectorNormForDistrsList,
				  maxNoOfMoments,
				  maxDegree,
				  extent,
				  atomListReordered,
				  &atomPermutation[0],
				  threshold,
				  boxList,
				  interactor,
				  0,
				  0,
				  numberOfLevels,
				  compute_gradient_also,
				  D_list,
				  private_gradient_list) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in do_interaction_recursive");
#ifdef _OPENMP
#pragma omp atomic
#endif
	  errorCount++;
	  continue;
	}
    }
#ifdef _OPENMP
  if(omp_get_thread_num() != 0) {
    /* collect data from different threads. */
#pragma omp critical
    {
      for(int i=0; i<noOfBasisFuncIndexPairs; i++)
        V_list[i] += private_V_list[i];
      if(compute_gradient_also)
	for(int i = 0; i < 3*molecule.getNoOfAtoms(); i++)
	  gradient_list[i] += private_gradient_list[i];
    }
    delete [] private_V_list;
    delete [] private_gradient_list;
  }
#endif
  delete []multipoleList;
  }

  timeMeterMainPart.print(LOG_AREA_INTEGRALS, "compute_V_and_gradient_linear main part");

  if(compute_gradient_also) {
    Util::TimeMeter timeMeterGradientMultipolePart;
    // Update gradient with multipole contributions
    for(int boxIdx = 0; boxIdx < boxSystem.totNoOfBoxes; boxIdx++) {
      atom_box_struct* currBox = &boxList[boxIdx];
      int noOfAtomsCurrBox = currBox->basicBox.noOfItems;
      // Go through all atoms in this box
      for(int j = 0; j < noOfAtomsCurrBox; j++) {
	int atomIndex = atomPermutation[currBox->basicBox.firstItemIndex + j];
	for(int coordIdx = 0; coordIdx < 3; coordIdx++) {
	  for(int k = 0; k < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; k++) {
	    ergo_real multipole_moment_derivative = currBox->multipole_moment_derivatives[j*MAX_NO_OF_MOMENTS_PER_MULTIPOLE*3 + MAX_NO_OF_MOMENTS_PER_MULTIPOLE*coordIdx + k];
	    ergo_real derivative_wrt_multipole_moment = currBox->derivatives_wrt_multipole_moments[k];
	    ergo_real gradient_contrib = multipole_moment_derivative * derivative_wrt_multipole_moment;
	    gradient_list[atomIndex*3+coordIdx] += gradient_contrib;
	  } // end for k
	} // end for coordIdx
      } // end for j
    } // end for boxIdx
    timeMeterGradientMultipolePart.print(LOG_AREA_INTEGRALS, "compute_V_and_gradient_linear gradient multipole contributions part");

    // Update gradient with nuclear-nuclear energy contribution
    Util::TimeMeter timeMeterNuclearRepulsionEnergyGradientContrib;
    molecule.getNuclearRepulsionEnergyGradientContribQuadratic(gradient_list);
    timeMeterNuclearRepulsionEnergyGradientContrib.print(LOG_AREA_INTEGRALS, "compute_V_and_gradient_linear getNuclearRepulsionEnergyGradientContribQuadratic");

    // Cleanup: delete buffers that were allocated during creation of multipole tree
    for(int boxIdx = 0; boxIdx < boxSystem.totNoOfBoxes; boxIdx++) {
      atom_box_struct* currBox = &boxList[boxIdx];
      delete [] currBox->multipole_moment_derivatives;
      currBox->multipole_moment_derivatives = NULL;
      delete [] currBox->derivatives_wrt_multipole_moments;
      currBox->derivatives_wrt_multipole_moments = NULL;
    }

  } // end if compute_gradient_also    

  delete [] boxList;
  delete [] atomListReordered;
  
  timeMeterTotal.print(LOG_AREA_INTEGRALS, "compute_V_and_gradient_linear total");

  return -errorCount;
}



int compute_V_hierarchical(const BasisInfoStruct& basisInfo,
			   const IntegralInfo& integralInfo,
			   const Molecule& molecule,
			   ergo_real threshold,
			   ergo_real boxSize,
			   const basis_func_index_pair_struct_1el* basisFuncIndexPairList,
			   int noOfBasisFuncIndexPairs,
			   csr_matrix_struct* V_CSR,
			   ergo_real & result_nuclearRepulsionEnergy
			   ) {
  result_nuclearRepulsionEnergy = 0; // computed later
  bool compute_gradient_also = false;
  int errorCount = 0;
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "compute_V_hierarchical start.");

  Util::TimeMeter timeMeterTotal;
  Util::TimeMeter timeMeterInitPart;

  // Get maxCharge
  ergo_real maxCharge = 0;
  for(int i = 0; i < molecule.getNoOfAtoms(); i++) {
    ergo_real currCharge = molecule.getAtom(i).charge;
    if(currCharge > maxCharge)
      maxCharge = currCharge;
  }

  // Create list of distributions
  int nDistrs = get_list_of_distrs_for_V(basisInfo, basisFuncIndexPairList, noOfBasisFuncIndexPairs,
					 threshold, maxCharge, NULL, 0);
  if(nDistrs <= 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_V_hierarchical: (nDistrs <= 0).");
    return -1;
  }
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "compute_V_hierarchical nDistrs    = %9i", nDistrs);
  std::vector<DistributionSpecStructWithIndexes> listTmp(nDistrs);
  if(get_list_of_distrs_for_V(basisInfo, basisFuncIndexPairList, noOfBasisFuncIndexPairs,
			      threshold, maxCharge, &listTmp[0], nDistrs) != nDistrs) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_V_hierarchical, in get_list_of_distrs_for_V.");
    return -1;
  }

  std::vector<DistributionSpecStructWithIndexes2> listTmp2(nDistrs);
  for(int i = 0; i < nDistrs; i++) {
    listTmp2[i].distr = listTmp[i].distr;
    listTmp2[i].basisFuncIdx1 = listTmp[i].basisFuncIdx1;
    listTmp2[i].basisFuncIdx2 = listTmp[i].basisFuncIdx2;
  }
  listTmp.clear();

  SetOfDistrsForV setOfDistrsForV;
  organize_distrs_for_V(integralInfo, setOfDistrsForV, listTmp2, threshold, maxCharge);
  listTmp2.clear();

  // Note that boxList and atomListReordered are allocated by create_nuclei_mm_tree,
  // we must remember to free them in the end.
  BoxSystem boxSystem;
  std::vector<int> atomPermutation(molecule.getNoOfAtoms());
  atom_box_struct* boxList = NULL;
  Atom *atomListReordered = NULL;
  int numberOfLevels = -1;
  if(create_nuclei_mm_tree(integralInfo,
			   molecule.getNoOfAtoms(), molecule.getAtomListPtr(), boxSize,
			   boxSystem, &boxList, &numberOfLevels, &atomListReordered,
			   &atomPermutation[0], compute_gradient_also) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "create_nuclei_mm_tree failed");
    return -1;
  }

  timeMeterInitPart.print(LOG_AREA_INTEGRALS, "compute_V_hierarchical init part (including create_nuclei_mm_tree)");

  // Compute nuclear repulsion energy using multipole tree.
  Util::TimeMeter timeMeterNuclRep;
  {
    MMInteractor interactor(integralInfo.GetMultipolePrep());
    int boxIndex1 = 0;
    int boxIndex2 = 0;
    int currLevel = 0;
    ergo_real nuclRepEnergy = get_nucl_repulsion_energy_using_multipoles(atomListReordered,
									 &atomPermutation[0],
									 threshold,
									 boxList,
									 interactor,
									 boxIndex1,
									 boxIndex2,
									 currLevel,
									 numberOfLevels);
    do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "get_nucl_repulsion_energy_using_multipoles() gave nuclRepEnergy = %22.11f", nuclRepEnergy);
    result_nuclearRepulsionEnergy = nuclRepEnergy;
  }
  timeMeterNuclRep.print(LOG_AREA_INTEGRALS, "get_nucl_repulsion_energy_using_multipoles");

  Util::TimeMeter timeMeterMainPart;

  const JK::ExchWeights CAM_params_not_used;
  MMInteractor interactor(integralInfo.GetMultipolePrep());
  int groupCount = setOfDistrsForV.groupList.size();
  for(int groupIndex = 0; groupIndex < groupCount; groupIndex++) {
    int groupStartIdx = setOfDistrsForV.groupList[groupIndex].startIndex;
    DistributionSpecStructWithIndexes2* currList = &setOfDistrsForV.distrList[groupStartIdx];
    int nDistrsCurrGroup = setOfDistrsForV.groupList[groupIndex].count;
    ergo_real extent = setOfDistrsForV.groupList[groupIndex].maxExtent;
    const multipole_struct_small* multipoleList = &setOfDistrsForV.multipoleList[groupStartIdx];
    const ergo_real* maxMomentVectorNormForDistrsList = setOfDistrsForV.maxMomentVectorNormList[groupIndex].maxMomentVectorNormList;
    int maxNoOfMoments = setOfDistrsForV.groupList[groupIndex].maxNoOfMoments;
    int maxDegree = setOfDistrsForV.groupList[groupIndex].maxDegree;
    // Take care of interaction of this group with MM tree
    if(do_interaction_recursive_2(integralInfo, V_CSR, noOfBasisFuncIndexPairs,
				  basisFuncIndexPairList, currList, nDistrsCurrGroup,
				  multipoleList, maxMomentVectorNormForDistrsList, maxNoOfMoments,
				  maxDegree, extent, atomListReordered,
				  &atomPermutation[0], threshold, boxList,
				  interactor, 0, 0, numberOfLevels) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in do_interaction_recursive");
      errorCount++;
      continue;
    }
  }

  timeMeterMainPart.print(LOG_AREA_INTEGRALS, "compute_V_hierarchical main part");

  delete [] boxList;
  delete [] atomListReordered;
  
  timeMeterTotal.print(LOG_AREA_INTEGRALS, "compute_V_hierarchical total");

  return -errorCount;
}


static ergo_real
simplePrimVintegralSingle(const DistributionSpecStruct & prim,
                          const Atom & atom,
                          const IntegralInfo & integralInfo) {
  return do_1e_repulsion_integral_using_symb_info(prim,
						  atom.charge,
						  atom.coords,
						  integralInfo);
}


ergo_real
simplePrimVintegral_list(const DistributionSpecStruct* list,
                         int nPrims,
			 const Atom & atom,
                         const IntegralInfo & integralInfo) {
  ergo_real sum = 0;
  for(int k = 0; k < nPrims; k++) {
    const DistributionSpecStruct & currDistr = list[k];
    sum += simplePrimVintegralSingle(currDistr, atom, integralInfo);
  }
  return sum;
}


int 
compute_V_matrix_full(const BasisInfoStruct& basisInfo,
		      const IntegralInfo& integralInfo,
		      int nAtoms,
		      const Atom* atomList,
		      ergo_real threshold,
		      ergo_real* result)
{
  int mu, nu, A, j, k, nbast;
  nbast = basisInfo.noOfBasisFuncs;

  for(mu = 0; mu < nbast; mu++)
    {
      BasisFuncStruct* basisFunc_mu = &basisInfo.basisFuncList[mu];
      int n_mu = basisFunc_mu->noOfSimplePrimitives;
      int start_prim_mu = basisFunc_mu->simplePrimitiveIndex;
      DistributionSpecStruct* list_mu = 
        &basisInfo.simplePrimitiveList[start_prim_mu];
      for(nu = 0; nu <= mu; nu++)
        {
          BasisFuncStruct* basisFunc_nu = &basisInfo.basisFuncList[nu];
          int n_nu = basisFunc_nu->noOfSimplePrimitives;
          int start_prim_nu = basisFunc_nu->simplePrimitiveIndex;
          DistributionSpecStruct* list_nu = 
            &basisInfo.simplePrimitiveList[start_prim_nu];
          /* compute matrix element [mu,nu] */
          ergo_real sum = 0;
          for(j = 0; j < n_mu; j++)
            {
              DistributionSpecStruct & prim_mu_j = list_mu[j];
              for(k = 0; k < n_nu; k++)
                {
                  DistributionSpecStruct & prim_nu_k = list_nu[k];
                  const int maxDistrsInTempList = 888;
                  DistributionSpecStruct tempList[maxDistrsInTempList];
                  int nNewPrims = get_product_simple_prims(prim_mu_j, 
                                                           prim_nu_k, 
                                                           tempList,
                                                           maxDistrsInTempList,
                                                           threshold);
                  if(nNewPrims < 0)
                    {
		      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_product_simple_prims");
                      exit(EXIT_FAILURE);
                    }
                  /* now loop over atoms */
                  for(A = 0; A < nAtoms; A++)
                    {
                      sum += simplePrimVintegral_list(tempList,
                                                      nNewPrims,
                                                      atomList[A],
                                                      integralInfo);
                    } /* END FOR A */
                } /* END FOR k */
            } /* END FOR j */
          result[mu*nbast+nu] = -1 * sum;
        } /* END FOR nu */
    } /* END FOR mu */

  // copy values to the other triangle
  for(mu = 0; mu < nbast; mu++)
    for(nu = mu+1; nu < nbast; nu++)
      result[mu*nbast+nu] = result[nu*nbast+mu];

  return 0;
}

