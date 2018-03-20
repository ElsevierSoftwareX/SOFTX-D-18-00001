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

/** @file grid_stream.cc @brief This is a streaming version of the
    linearly-scaling grid generator.

    The grid is generated on the fly, without the in-memory sort
    phase. The disk format is shared with the ordinary grid. The
    separation of the ordinary grid data (coordinates, weights) from
    the auxiliary data (lists of active orbitals) is taken into
    account so that grid can be used for integration in the auxiliary
    basis.
*/

#include <cmath>
#include <stdio.h>
#include <assert.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <vector>

#include "dft_common.h"
#include "grid_stream.h"
#include "grid_atomic.h"
#include "lebedev_laikov.h"
#include "molecule.h"
#include "output.h"
#include "realtype.h"
#include "utilities.h"

/* FIXME: remove this dependency */
#include "sparse_matrix.h"

//#define BEGIN_NAMESPACE(x) namespace x {
// #define END_NAMESPACE(x)   } /* x */

typedef ergo_real real;
/* NAMESPACE NAMESPACE NAMESPACE NAMESPACE NAMESPACE NAMESPACE NAMESPACE */
BEGIN_NAMESPACE(Grid)

/** Ignore all grid points that partitioning scales down by more than
    WEIGHT_THRESHOLD */
static const real WEIGHT_THRESHOLD = 1e-15;

/** A grid describing a radial grid for an atom with a specific
    charge. */
class RadialGrid {
public:
  real *rad;      /**< Array of radial grid points */
  real *weights;  /**< Array of the weights associated with the grid points */
  int  *nAngular; /**< array of sizes of corresponding angular grids. */
  int noOfRadPoints;
  RadialGrid(int charge_, RadialScheme* rs, int angMin, int angMax);

  int getCharge() const  { return charge; }
  real getRadius() const { return rad[noOfRadPoints-1]; }
  unsigned getPointCount() const {
    unsigned s = 0;
    for(int j=0; j<noOfRadPoints; ++j)
      s += nAngular[j];
    return s;
  }
  ~RadialGrid() { 
    delete []rad; delete []weights; delete []nAngular;
  }
protected:
  void setAngularFixed(int minAng, int maxAng);
private:
  int charge;
};


void
RadialGrid::setAngularFixed(int minAng, int maxAng)
{
  static const real BOHR = 1.0/0.529177249;

  real rBragg = BraggRadii[charge]*BOHR;
  /* maxAng specified for carbon-type atoms, corrections for ligher
     and heavier atoms. */
  if(charge<=2)
    maxAng -= 10;
  else if(charge>10) {
    if(charge<=18)
      maxAng += 5;
    else if(charge<=36)
      maxAng += 10;
    else maxAng += 15;
  }

  int currentAng = maxAng;
  int maxAngular = ll_npoint(maxAng);
  for (int i=0; i<noOfRadPoints; i++) {
    if(rad[i]<rBragg) { /* Prune */
      int iAngPoints = int(maxAngular*rad[i]/rBragg);
      currentAng = ll_order(iAngPoints);
      if(currentAng<minAng) currentAng = minAng;
    } /* else current_ang = maxang; */
    nAngular[i] = ll_npoint(currentAng);
  }
}

RadialGrid::RadialGrid(int charge_, RadialScheme* rs,
                       int angMin, int angMax) : charge(charge_)
{
  noOfRadPoints = rs->gridSize;
  rad = new real[noOfRadPoints];
  weights = new real[noOfRadPoints];
  nAngular = new int[noOfRadPoints];
  rs->generate(rad, weights);
  setAngularFixed(angMin, angMax);
}
  
class AtomicGrid {
  const RadialGrid* rGrid;
public:  
  Vector3D center;
  const RadialGrid& getRadialGrid() const { return *rGrid; }
  /** Returns "radius" of an atomic grid. It stretches really bit
      longer than the last grid point. */
  real radius() const {
    return rGrid->getRadius();
  }
  int charge() const { return rGrid->getCharge(); }
  unsigned getPointCount() const { return rGrid->getPointCount(); }
  
  AtomicGrid(const AtomicGrid& a) :
    rGrid(a.rGrid), center(a.center) {}

  AtomicGrid(const real c[3], const RadialGrid* rg) :
    rGrid(rg), center(c[0], c[1], c[2]) { }

  AtomicGrid(const Atom& atom, const RadialGrid* rg) :
    rGrid(rg), center(atom.coords[0],atom.coords[1],atom.coords[2]) { }
};


/** "Block" partitioning is the only one implemented now... We rename
    it here to Box partitioner to avoid name space conflicts.
*/
class BoxPartitioner {
  const real (*coor)[3];
  real *rj;
  long_real *p_kg;
  long_real *vec;
  real *invAtomDistances; /**< a triangular array */
  real *atomFactors;   /**< a triangular array */
  static const int HARDNESS = 11;
  real xpasc[HARDNESS], apasc;

  unsigned maxPointCnt;
  unsigned maxAtomPointCnt;
  unsigned maxRelevantAtoms;
  int LDA; /* leading dimension of rj and p_kg */
  void prepare(const std::vector<AtomicGrid>& atoms,
               unsigned noOfRelevantAtoms, const unsigned *relevantAtoms,
               unsigned pointCnt, const real (*gridPoints)[3]);
  /** return distance between given pair of relevant atoms. Arguments
      i and j specify the number of the atoms on the list of relevant
      atoms. */

  real getInvDistanceBetweenAtoms(int i, int j) const {
    return i>j
      ? invAtomDistances[(i*(i-1))/2 + j] : invAtomDistances[(j*(j-1))/2 +i]; 
  }
 
  real getFactor(int i, int j) const {
    return i>j
      ? atomFactors[(i*(i-1))/2 + j] : atomFactors[(j*(j-1))/2 +i];
  }
    

public:
  BoxPartitioner();
  ~BoxPartitioner();
  unsigned process(unsigned atomNumber,
		   const std::vector<AtomicGrid>& atomGrids,
		   int noOfRelevantAtoms,
		   const unsigned *relevantAtoms,
		   unsigned batchLength,
		   real (*coor)[3], real *w);

};

/** Initializez the BoxPartitioner object. */
BoxPartitioner::BoxPartitioner()
  : coor(NULL), rj(NULL), p_kg(NULL), vec(NULL),
    invAtomDistances(NULL), atomFactors(NULL),
    maxPointCnt(0), maxAtomPointCnt(0), maxRelevantAtoms(0)
{
  int h;
  real facult[HARDNESS], isign = -1;

  facult[0] = 1;
  for(h=1; h<HARDNESS; h++)
    facult[h] = facult[h-1]*h;
  
  for(h=HARDNESS-1; h>=0; h--) {
    isign = -isign;
    xpasc[h] = isign*facult[HARDNESS-1]/
      ((2*h+1)*facult[h]*facult[HARDNESS-1-h]);
  }
  xpasc[0] = 1;
  apasc = 0;
  for(h=0; h<HARDNESS; h++) apasc += xpasc[h];
  apasc = 0.5/apasc;
}


BoxPartitioner::~BoxPartitioner()
{
  if(rj)   delete []rj;
  if(p_kg) delete []p_kg;
  if(vec)  delete []vec;
  if(invAtomDistances) delete []invAtomDistances;
  if(atomFactors)   delete []atomFactors;
}


void
BoxPartitioner::prepare(const std::vector<AtomicGrid>& atoms,
                        unsigned noOfRelevantAtoms,
                        const unsigned *relevantAtoms,
                        unsigned pointCnt, const real (*gridPoints)[3])
{
  unsigned ptno;

  coor = gridPoints;
  LDA  = pointCnt;

  /* FIXME: optimize reallocation strategies to avoid increases by 1? */
  if (noOfRelevantAtoms*pointCnt > maxAtomPointCnt) {
    maxAtomPointCnt = noOfRelevantAtoms*pointCnt;
    if (rj)   delete []rj;
    if (p_kg) delete []p_kg;
    rj   = new real[noOfRelevantAtoms*pointCnt];
    p_kg = new long_real[noOfRelevantAtoms*pointCnt];
  }
  if (pointCnt > maxPointCnt) {
    maxPointCnt = pointCnt;
    if (vec) delete []vec;
    vec  = new long_real[pointCnt];
  }

  if (noOfRelevantAtoms > maxRelevantAtoms) {
    maxRelevantAtoms = noOfRelevantAtoms;
    if(invAtomDistances) {
      delete []invAtomDistances;
      delete []atomFactors;
    }
    invAtomDistances = new real[(noOfRelevantAtoms*(noOfRelevantAtoms-1))/2];
    atomFactors   = new real[(noOfRelevantAtoms*(noOfRelevantAtoms-1))/2];
  }
  memset(vec, 0, pointCnt*sizeof(vec[0]));

  for(unsigned i=0; i<noOfRelevantAtoms; i++) {
    int atno = relevantAtoms[i];
    real iradius = 
      BraggRadii[(atoms[atno].charge()>=(signed)BraggSize
                  ? BraggSize-1 : atoms[atno].charge())];
    
    for(unsigned j=0; j<i; j++) {
      int atno2 = relevantAtoms[j];
      real d = atoms[atno].center.dist(atoms[atno2].center);
      real jradius = 
        BraggRadii[(atoms[atno2].charge()>=(signed)BraggSize
                    ? BraggSize-1 : atoms[atno2].charge())];
      real chi = template_blas_sqrt(iradius/jradius);
      real temp = (chi-1)/(chi+1);
      temp = temp/(temp*temp-1);
      if(temp>0.5) temp = 0.5;      
      else if(temp<-0.5) temp = -0.5;
      atomFactors[(i*(i-1))/2 + j] = temp;
      invAtomDistances[(i*(i-1))/2 + j] = 1.0/d;
    }

    real *rji   = rj   + i*LDA;
    long_real *p_kgi = p_kg + i*LDA;
    for(ptno=0; ptno<pointCnt; ptno++) {
      real dx = atoms[atno].center[0] - gridPoints[ptno][0];
      real dy = atoms[atno].center[1] - gridPoints[ptno][1];
      real dz = atoms[atno].center[2] - gridPoints[ptno][2];

      rji  [ptno] = template_blas_sqrt(dx*dx + dy*dy + dz*dz);
      p_kgi[ptno] = 1.0;
    }
  }
}


/** Applies the partitioning factors to the the grid points associated
   with given atom.

   @param atomNumber index of the atom that the grid
   points are associated with, in the relevantAtoms array.

   @param atomGrids the list of all atom grids.

   @param noOfRelevantAtoms the number of atoms relevant for the box
   being processed. They need to be taken into account when computing
   the partitioning weights.

   @param relevantAtoms the indices of the relevant atoms in the
   atomGrids vector.

   @param batchLength number of the grid points the partitioning
   weights have to be computed for.

   @param coor their cartesian coordinates. Will be modified if the
   grid point compression occurs.

   @param w Their weights - they will be modified.
   @returns new batch length after compression.
*/
unsigned
BoxPartitioner::process(unsigned atomNumber,
                        const std::vector<AtomicGrid>& atomGrids,
                        int noOfRelevantAtoms,
                        const unsigned *relevantAtoms,
                        unsigned batchLength,
                        real (*coor)[3], real *w)
{
  int i, j;
  unsigned ptno;

  prepare(atomGrids, noOfRelevantAtoms, relevantAtoms, batchLength, coor);

  for(i=1; i<noOfRelevantAtoms; i++) {
    const real *distAtomGridI = rj + LDA*i;
    long_real *p_kgI = p_kg + LDA*i;

    for(j=0; j<i; j++) {
      const real *distAtomGridJ = rj + LDA*j;
      long_real *p_kgJ = p_kg + LDA*j;
      long_real mu, mu2, g_mu;

      real invDist = getInvDistanceBetweenAtoms(i, j);
      real bFac = getFactor(i, j);
      for(ptno=0; ptno<batchLength; ptno++) {
        mu =(distAtomGridI[ptno]-distAtomGridJ[ptno])*invDist;
#if 0
        if( std::fabs(mu>1))
          printf("Too large mu: %g by %g ptno: %i\n", mu, 1-std::fabs(mu), ptno);
#endif
        if(mu>1.0) mu = 1.0;
        else if(mu<-1.0) mu = -1.0;

        mu += bFac*(1-mu*mu);
        mu2 = mu*mu;
        g_mu = 0;
        for(int h=0; h<HARDNESS; h++) {
          g_mu += xpasc[h]*mu;
          mu *= mu2;
        }
        long_real fac = apasc*g_mu;
        if(fac>0.5) fac = 0.5;
        else if(fac<-0.5) fac = -0.5;
        p_kgI[ptno] *= 0.5-fac;
        p_kgJ[ptno] *= 0.5+fac;
      }
    }
  }

  /* compute weight normalization factors */
  for(i=0; i<noOfRelevantAtoms; i++) {
    const long_real *p_kgI = p_kg + LDA*i;
    for(ptno=0; ptno<batchLength; ptno++)
      vec[ptno] += p_kgI[ptno];
  }

  /*
   * Apply the computed weights.
   */
  const long_real *myP_kg = p_kg + LDA*atomNumber;
  unsigned dest = 0;
  for(unsigned src=0; src<batchLength; src++) {
    real fac = myP_kg[src]/vec[src];
    if(fac>=WEIGHT_THRESHOLD) {
      w[dest] = w[src]*fac;
      coor[dest][0] = coor[src][0];
      coor[dest][1] = coor[src][1];
      coor[dest][2] = coor[src][2];
      dest++;
    }
  }

  return dest;
}

/** A class that is able to quickly determine the active shells that
    overlap with given box in space. */
class ActiveBfShells {
  const GridGenMolInfo& ggmi;
  real *rShell2;
public:
  explicit ActiveBfShells(const GridGenMolInfo& ggmi_)
    : ggmi(ggmi_), rShell2(new real[ggmi_.noOfShells])
  {
    ggmi.setShellRadii(rShell2);
  }
  int getMaxShells() const {
    return ggmi.noOfShells;
  }
  ~ActiveBfShells() { 
    delete []rShell2;
  }
  void setForBox(const Box& b, int *nBlocks, int (*shlBlocks)[2]) const;
  
   /**< the start and stop+1 indexes. */
  static int getNoOfShells(int nBlocks,
			   int (*shlBlocks)[2]) {
    int sum = 0;
    for(int block=0; block<nBlocks; block++)
      sum += shlBlocks[block][1]-shlBlocks[block][0];
    return sum;
  }
};

void
ActiveBfShells::setForBox(const Box& box,
			  int *nBlocks, int (*shlBlocks)[2]) const
{
  real center[3];
  for(int i=0; i<3; i++)
    center[i] = 0.5*(box.lo[i]+box.hi[i]);
  real cellSize = box.size(box.getMaxDim());

  ggmi.getBlocks(center, cellSize, rShell2, nBlocks, shlBlocks);
}

/** Saves the grid saving context. Contains the information that is
    changed as the grid tree is traversed while it is being saved. */
struct StreamSaveContext {
  FILE *stream;
  int (*shlBlocks)[2];
  BoxPartitioner& partitioner;
  unsigned savedPoints;
  unsigned boxCount;
  unsigned myRank;
  StreamSaveContext(FILE *f, BoxPartitioner& p, unsigned maxShells,
		    unsigned rank)
    : stream(f), shlBlocks(new int[maxShells][2]),
      partitioner(p), savedPoints(0), boxCount(0), myRank(rank) {}
  ~StreamSaveContext()
  { delete []shlBlocks; }
};

/** Streamlined, abstract grid generation class. This class does not
    depend explicitly on the representation of the basis set and
    molecule. All such dependence is abstracted away in case the grid
    generator is to be used with another program. */
class Stream {
public:
  real boxSize;
  bool saveToFile(const char *fname, int noOfThreads);
  Stream(ActiveBfShells& abs, RadialScheme *rs,
         real radint, int angmin, int angmax,
	 real boxSize, bool forceCubic);
  ~Stream();
  unsigned getPointCount() const { return savedPoints; }
  Dft::SparsePattern *sparsePattern;
protected:
  bool forceCubic;
  bool saveAtomsRecursively(StreamSaveContext& ssc, const Box& box,
			    unsigned cnt, const unsigned atoms[],
			    int depth) const;
  unsigned saveAtomGridInBox(unsigned iAtom, const Box& box,
			     BoxPartitioner& partitioner,
			     unsigned cnt, const unsigned atoms[],
			     int (*shlBlocks)[2],
			     FILE *stream) const;
  unsigned noOfAtoms() const;
  const AtomicGrid& getAtomicGrid(unsigned i) const;
  void addAtom(const real coor[3], int charge, int atomNo);
  void addAtom(const Atom& at, int atomNo) {
    addAtom(at.coords, int(at.charge), atomNo);
  }

  unsigned saveBatch(unsigned batchLength, real (*coor)[3],
		     real *weight,
		     unsigned nBlocks, int (*shlBlocks)[2],
		     FILE *f) const;
  static void* saveThread(void *data);
private:
  static pthread_mutex_t fileSaveMutex;
  /** We store just one entry per atom type - there is no reason to
      have thousands identical copies. */
  std::vector<RadialGrid*> atomTypes;
  std::vector<AtomicGrid>  atoms;
  ActiveBfShells&          activeShells;
  RadialScheme*            radialScheme;
  real radialThreshold;
  int angularMin, angularMax;
  /* unsigned long contributionsToCompute; nPoints*nShells*nShells */
  unsigned savedPoints;
  unsigned expectedPoints;
  int noOfThreads;
};

pthread_mutex_t Stream::fileSaveMutex = PTHREAD_MUTEX_INITIALIZER;


const AtomicGrid&
Stream::getAtomicGrid(unsigned i) const 
{
  return atoms[i];
}

inline unsigned
Stream::noOfAtoms() const
{
  return atoms.size();
}

void
Stream::addAtom(const real coor[3], int charge, int atomNo)
{
  std::vector<RadialGrid*>::iterator i;
  for(i= atomTypes.begin(); i != atomTypes.end(); ++i)
    if( (*i)->getCharge() == charge)
      break;
 
  if(i == atomTypes.end()) {
    radialScheme->init(atomNo, charge, radialThreshold);
    RadialGrid *rg = new RadialGrid( charge, radialScheme,
                                     angularMin, angularMax);
    atomTypes.push_back(rg);
    atoms.push_back(AtomicGrid(coor, rg));
    do_output(LOG_CAT_INFO, LOG_AREA_DFT,
              "Grid for charge %d: %d radial points, %d total, r=%f",
	      charge, rg->noOfRadPoints,
	      rg->getPointCount(), (double)rg->getRadius());
#if 0
    printf("Added atom type with charge %f and radius %g, %u points\n",
           charge, rg->getRadius(), rg->getPointCount());
#endif
  } else {
    atoms.push_back(AtomicGrid(coor, *i));
  }
  AtomicGrid &ag = atoms.back();
  unsigned nPoints = ag.getPointCount();
  expectedPoints += nPoints;
}

/** Saves a batch of grid points to given file. */
unsigned
Stream::saveBatch(unsigned batchLength, real (*coor)[3],
                  real *weight,
		  unsigned nBlocks, int (*shlBlocks)[2],
		  FILE *f) const
{  
  pthread_mutex_lock(&fileSaveMutex);
 
  if(fwrite(&batchLength, sizeof(batchLength), 1, f)!=1            ||
     fwrite(&nBlocks,  sizeof(nBlocks), 1, f) != 1                 ||
     fwrite(shlBlocks, sizeof(int), nBlocks*2, f) != nBlocks*2     ||
     fwrite(coor, sizeof(real), 3*batchLength, f) != 3*batchLength ||
     fwrite(weight, sizeof(real), batchLength, f) != batchLength)
    throw "Cannot save grid point batch on file";
  /*unsigned nShells = activeShells.getNoOfShells();
    contributionsToCompute += batchLength*nShells*nShells; */

  if(sparsePattern)
    sparsePattern->add(nBlocks, shlBlocks);
  pthread_mutex_unlock(&fileSaveMutex);
  return batchLength;
}


/** Method saves all grid points associated with specified atom,
    located in the specified box. It will also save the associated
    auxiliary information (usually a list of active orbitals) - this
    is why we pass an atom list. FIXME: this probably needs to be
    thought through: what factor decides the atom radius, really? Is
    it max(auxiliaryRadius, gridRadius)?
*/
unsigned
Stream::saveAtomGridInBox(unsigned iAtom, const Box& box,
                          BoxPartitioner& partitioner,
                          unsigned cnt, const unsigned relevantAtoms[],
			  int (*shlBlocks)[2],
                          FILE *stream) const
{
  /* MAX_BATCH_LENGTH Must be larger than the densest angular grid. */
  /* Elias note 2016-07-08: The MAX_BATCH_LENGTH value here was
     decreased from 8192 to 4096 because the 8192 was apparently too
     large when running on a Mac, it seemed like the threads ran out
     of stack space then. */
  /* Elias note 2017-12-21: The MAX_BATCH_LENGTH value here was
     decreased again, this time from 4096 to 2048, because the 4096
     value was apparently too large when running on a Mac with long
     double precision, there was a "bus error" problem then that
     seemed to be fixed by reducing MAX_BATCH_LENGTH. */
  /* FIXME: there does not seem to be proper checking if the
     MAX_BATCH_LENGTH value is too small, it just gives wrong results
     without detecting the error. */
  static const unsigned MAX_BATCH_LENGTH = 2048;
  const AtomicGrid& atom = atoms[relevantAtoms[iAtom]];
  real dist = box.getDistanceTo(atom.center.v);
  const RadialGrid& rGrid = atom.getRadialGrid();
  assert(rGrid.noOfRadPoints>0);
  assert(rGrid.rad[rGrid.noOfRadPoints-1]>rGrid.rad[0]);

  real coor[MAX_BATCH_LENGTH][3];
  real weight[MAX_BATCH_LENGTH];
  real angX[MAX_BATCH_LENGTH], angY[MAX_BATCH_LENGTH], angZ[MAX_BATCH_LENGTH];
  real angW[MAX_BATCH_LENGTH];
  unsigned usedPoints = 0;
  int nBlocks; 
  unsigned savedPoints = 0;
#if 0
  /* amazingly enough, these memsets quieten valgrind warnings with
     long double calculations. */
  memset(coor, 0, sizeof(real)*3*MAX_BATCH_LENGTH);
  memset(weight, 0, sizeof(real)*MAX_BATCH_LENGTH);
#endif
  activeShells.setForBox(box, &nBlocks, shlBlocks);
  if(nBlocks == 0) {
#if 0
    printf("Got zero AO blocks for box (%f,%f,%f)-(%f,%f,%f).\n",
           box.lo[0], box.lo[1], box.lo[2],
           box.hi[0], box.hi[1], box.hi[2]);
#endif
    return 0;
  }
  for (int i=0; i<rGrid.noOfRadPoints; i++) {
    /* Skip early entire grid point shells that do not reach the box
       in question. Point by point screening is done after grid shell
       generation. */
    if (rGrid.rad[i] < dist)
      continue;
    /* Create the part of the grid associated with the atom that lies
       within the box. */
    int nAngPoints = rGrid.nAngular[i];
    if(nAngPoints > (signed)MAX_BATCH_LENGTH)
      throw "MAX_BATCH_LENGTH too small!";

    /* Add points to the batch here... */
    real radius = rGrid.rad[i];
    real rWeight = rGrid.weights[i]*4.0*M_PI;
    ll_sphere(nAngPoints, angX, angY, angZ, angW);
#if 0
    printf("B %2i: r %8.2f W: %8.2g Ang: %d| %g %g %g\n", i, radius, rWeight,
           nAngPoints, angX[0], angY[0], angZ[0]);
#endif
    for(int pnt=0; pnt<nAngPoints; pnt++) {
      coor[usedPoints][0] = atom.center[0] + angX[pnt]*radius;
      coor[usedPoints][1] = atom.center[1] + angY[pnt]*radius;
      coor[usedPoints][2] = atom.center[2] + angZ[pnt]*radius;
      if(box.contains(coor[usedPoints]) ) {
        weight[usedPoints] = angW[pnt]*rWeight;
        usedPoints++;
        
        if(usedPoints == MAX_BATCH_LENGTH) {
          usedPoints = partitioner.process(iAtom, atoms, cnt, relevantAtoms,
					   usedPoints, coor, weight);
	  if(usedPoints) {
	    savedPoints +=
	      saveBatch(usedPoints, coor, weight, nBlocks, shlBlocks, stream);
	    usedPoints = 0;
	  }
        }
      }
    }
  }

  /* Flush the last batch for this atom... Consider merging points
     from different atoms to a single batch after partitioning. */
  if(usedPoints) {
    usedPoints = partitioner.process(iAtom, atoms, cnt, relevantAtoms,
				     usedPoints, coor, weight);
    if(usedPoints) {
      savedPoints +=
	saveBatch(usedPoints, coor, weight, nBlocks, shlBlocks, stream);
    }
  }
  return savedPoints;
}

/** This is a recursive procedure that generates the grid points that lie
    in the specified bounding box. As an optimization, a list of atoms
    that may overlap with given grid is passed so that linear scaling
    can be achieved. This goal is achieved by recursive division of
    the bounding box until there are no atoms that can overlap with
    it, or the minimal size is achieved. In the latter case, all atoms
    are iterated over and the grid points associated with them that
    lie in the bounding box are saved. An associated auxiliary
    information is saved as well.

    An atom is considered relevant for given box, if its Voronoi
    polyhedra (+safety margin) overlaps with the box.

    @param ssc the saving context containing information about
    selected partitioner and other grid generation specifics.

    @param box    save only the points within the box.

    @param atomCount the number of potentially relevant atoms that have 
    grids overlapping with the box in question..

    @param atomIndices ... and their indices in the global array.

    @param depth the recursion depth for logging purposes.
 */
bool
Stream::saveAtomsRecursively(StreamSaveContext& ssc, const Box& box,
                             unsigned atomCount,
			     const unsigned atomIndices[],
                             int depth) const
{
    /** randomly chosen parameter. We need in general a better way to
        determing whether the voronoi polyhedra associated with a
        given atom overlaps with the box in question... */  
  static const ergo_real POLYHEDRA_SAFETY_MARGIN = 2.0;
  if (atomCount == 0)
    return true;
  unsigned maxDim = box.getMaxDim();
  real maxSize = box.size(maxDim);

#ifdef DEBUG
  for(int d=0; d<depth; d++) printf("  ");
  printf("Recursive: (%f,%f,%f) - (%f,%f,%f) atoms: %i maxSize: %f\n",
         box.lo[0], box.lo[1], box.lo[2],
         box.hi[0], box.hi[1], box.hi[2], atomCount, maxSize);
  for(int d=0; d<depth; d++) printf("  ");
  for(int d=0; d<atomCount; ++d) printf("%d ", atomIndices[d]);
  puts("");
#endif

  /* Now, we divide the boxes really until their smallest size because
     the basis function range may be larger than the atomic grid
     range. */
  if ( /* atomCount > 1 &&*/ maxSize > boxSize) {
    bool res = true;

    std::vector<unsigned> localIndices;
    localIndices.reserve(atomCount);

    real splitPos = box.lo[maxDim] + 0.5*maxSize;
    Box newBox = box;
    newBox.hi[maxDim] = splitPos;
    
    for(unsigned i=0; i<atomCount; i++) {
      const AtomicGrid& g = getAtomicGrid(atomIndices[i]);
      if (newBox.overlapsWith( g.center.v, g.radius()+POLYHEDRA_SAFETY_MARGIN))
        localIndices.push_back(atomIndices[i]);
    }

    if (!localIndices.empty())
      res = saveAtomsRecursively(ssc, newBox,
                                 localIndices.size(),
                                 &(*localIndices.begin()), depth+1);

    newBox.hi[maxDim] = box.hi[maxDim];
    newBox.lo[maxDim] = splitPos;
    localIndices.clear();

    for(unsigned i=0; i<atomCount; i++) {
      const AtomicGrid& g = getAtomicGrid(atomIndices[i]);
      if (newBox.overlapsWith( g.center.v, g.radius()+POLYHEDRA_SAFETY_MARGIN))
        localIndices.push_back(atomIndices[i]);
    }
    if (!localIndices.empty())
      res = saveAtomsRecursively(ssc, newBox,
                                 localIndices.size(),
                                 &(*localIndices.begin()), depth+1);

    return res;
    
  } else { /* small box - or just one atom left... */
    ssc.boxCount++;
    if( ssc.boxCount % noOfThreads != ssc.myRank) return true;
    /* get the list of relevant basis shells here. */
    //    unsigned c = savedPoints;
    for(unsigned i=0; i<atomCount; i++) {
      const AtomicGrid& ga = getAtomicGrid(atomIndices[i]);
      if (!box.overlapsWith(ga.center.v, ga.radius()+POLYHEDRA_SAFETY_MARGIN)) {
        printf("Atom %i -> %i absolute does not overlap\n", i, atomIndices[i]);
        continue;
      }
      /* Save all the grid points of this atoms in the given box, does
         the partitioning as well. */
      ssc.savedPoints +=
	saveAtomGridInBox(i, box, ssc.partitioner,
			  atomCount, atomIndices,
			  ssc.shlBlocks, ssc.stream);
    }
#if 0
    printf("BOX (%f,%f,%f)-(%f,%f,%f) contains %u grid points\n",
           box.lo[0], box.lo[1], box.lo[2],
           box.hi[0], box.hi[1], box.hi[2],
           savedPoints - c);
#endif
  }
  return true;
}

struct ThreadInfo {
  pthread_t tid;
  const Stream* stream;
  FILE *f;
  const Box* box;
  unsigned atomCnt;
  const unsigned *atoms;
  unsigned myRank;
  unsigned savedPoints;
  bool res;
};
void*
Stream::saveThread(void *data) {
  static const int SAVEGRID_ERROR = 0;
  ThreadInfo *ti = (ThreadInfo*)data;
  try {
    BoxPartitioner partitioner;
    StreamSaveContext ssc(ti->f,partitioner,
			  ti->stream->activeShells.getMaxShells(), ti->myRank);
    ti->res =
      ti->stream->saveAtomsRecursively(ssc, *ti->box,
				       ti->atomCnt, ti->atoms, 0);
    ti->savedPoints = ssc.savedPoints;
  } catch(const char *s) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, 
	      "Stream::saveThread: xcWorker thread caught an exception '%s'", s);
    return (void*)&SAVEGRID_ERROR;
  } catch(const std::bad_alloc & e) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, 
	      "Stream::saveThread: xcWorker thread caught an exception '%s'", e.what());
    return (void*)&SAVEGRID_ERROR;
  } catch(const std::runtime_error & e) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, 
	      "Stream::saveThread: xcWorker thread caught an exception '%s'", e.what());
    return (void*)&SAVEGRID_ERROR;
  }  catch(...) {
    do_output(LOG_CAT_ERROR, LOG_AREA_DFT, 
	      "Stream::saveThread: xcWorker thread caught unexpected exception.");
    return (void*)&SAVEGRID_ERROR;
  }
  
  return NULL;
}
/** Generates the grid and saves it to the file with given name. 
    @return Number of saved grid points.
 */
bool
Stream::saveToFile(const char *fname, int noOfThreads)
{
  bool res;
  FILE *f = fopen(fname, "wb");
  if (!f) return false;
  unsigned n = noOfAtoms();
  unsigned *atomIndices = new unsigned[n];

  for(unsigned i=0; i<n; i++) atomIndices[i] = i;

  Box box;
  getBoundingBox(box, atoms.begin(), atoms.end());
  /* FIXME: Correct the box size so that the smallest level boxes are
     about cubic! */
  if (forceCubic) {
    for(int i=0; i<3; ++i) { 
      real center = 0.5* (box.lo[i] + box.hi[i]);
      real dim = 0.5* (box.hi[i] - box.lo[i]);
      int x, dimi = int(std::ceil((long double)(dim/boxSize)));
      for(x=1; x<dimi; x *=2);
      dim = x*boxSize;
      box.lo[i] = center-dim;
      box.hi[i] = center+dim;
    }
  }
  /* End of fix*/

  this->noOfThreads = noOfThreads;
  if(noOfThreads<=1) {
    BoxPartitioner p;
    StreamSaveContext ssc(f, p, activeShells.getMaxShells(), 0);
    res = saveAtomsRecursively(ssc, box, n, atomIndices, 0);
    savedPoints = ssc.savedPoints;
  } else {
    std::vector<ThreadInfo> threads(noOfThreads);
    for(int i=0; i<noOfThreads; i++) {
      threads[i].stream = this;
      threads[i].f = f;
      threads[i].box = &box;
      threads[i].atomCnt = n;
      threads[i].atoms = atomIndices;
      threads[i].myRank = i;
      if(pthread_create(&threads[i].tid, NULL, saveThread, &threads[i]) != 0)
	res = false;
    }
    savedPoints = 0;
    res = true;
    for(int i=0; i<noOfThreads; i++) {
      void *ptr;
      pthread_join(threads[i].tid, &ptr);
      if(ptr != NULL)
	res = false;
      savedPoints += threads[i].savedPoints;
      res = res && threads[i].res;
    }
  }    
  delete []atomIndices;

  if(fclose(f) != 0)
    throw "Cannot close the grid file";

  return res;
}

/** The Stream constructor. Takes over the radial scheme object, which
    has to be allocated dynamically. */
Stream::Stream(ActiveBfShells& abs, RadialScheme *rs,
               real radint, int angmin, int angmax,
	       real boxSize_, bool forceCubic_)
  : boxSize(boxSize_), sparsePattern(NULL), forceCubic(forceCubic_),
    activeShells(abs), radialScheme(rs), radialThreshold(radint),
    angularMin(angmin), angularMax(angmax),
    savedPoints(0), expectedPoints(0), noOfThreads(1)
{
}

Stream::~Stream()
{
  delete radialScheme;
  for(std::vector<RadialGrid*>::iterator i= atomTypes.begin();
      i != atomTypes.end(); ++i) {
    delete *i;
  }

  /* It can happen with GC2 radial grid that some of the grid points
     are ignored because no AO function reaches there. */
  if (savedPoints != expectedPoints) {
    do_output(LOG_CAT_INFO, LOG_AREA_DFT,
              "Grid pruned from %d to %d grid points.",
              expectedPoints, savedPoints);
    /* throw "Grid generator lost some points"; */
  }
}

END_NAMESPACE(Grid)

/* NAMESPACE NAMESPACE NAMESPACE NAMESPACE NAMESPACE NAMESPACE NAMESPACE */

/** Ergo-specific GridStream implementation. */
class ErgoGridStream : public Grid::Stream {
  Grid::ActiveBfShells ergoShells;
public:
  const Dft::GridParams& gsSettings;
  ErgoGridStream(const Dft::GridParams& gss, const GridGenMolInfo& molInfo,
		 RadialScheme *rs);
};

ErgoGridStream::ErgoGridStream(const Dft::GridParams& gss,
			       const GridGenMolInfo& molInfo,
                               RadialScheme *rs) 
  : Stream(ergoShells, rs, gss.radint, gss.angmin, gss.angmax,
	   gss.boxSize, gss.cubicBoxes),
    ergoShells(molInfo), gsSettings(gss)
{
  for(int i=0; i<molInfo.noOfAtoms; i++) {
    int dummy, charge, mult;
    real c[3];
    molInfo.getAtom(i, &dummy, &c, &charge, &mult);
    addAtom(c, charge, i);
  }
}

/** Creates the grid object. The Settings object must have longer
    lifetime than the grid itself - its content is not copied. */
ErgoGridStream*
grid_stream_new(const struct Dft::GridParams& gss,
		const GridGenMolInfo& molInfo)
{
  RadialScheme *radScheme;
  switch(gss.radialGridScheme) {
  default:
  case Dft::GridParams::LMG:
    radScheme = new RadialSchemeLMG(molInfo);
    break;
  case Dft::GridParams::GC2:
    radScheme = new RadialSchemeGC2();
    break;
  case Dft::GridParams::TURBO:
    radScheme = new RadialSchemeTurbo();
    break;
  }
  return new ErgoGridStream(gss, molInfo, radScheme);
}


void
grid_stream_set_sparse_pattern(ErgoGridStream *stream,
                               Dft::SparsePattern *pattern)
{
  stream->sparsePattern = pattern;
}

/** Generate grid for given molecule.
    @param stream The grid object.
    @param fname The file name the grid is to be saved to.
    @param noOfThreads the number of threads that are supposed to be
    created and used for the grid generation.
*/
unsigned
grid_stream_generate(ErgoGridStream *stream, const char *fname,
		     int noOfThreads)
{
  unsigned res;
  Util::TimeMeter tm;
  if(noOfThreads<1)
    noOfThreads = 1;
  try {
    const Dft::GridParams& gss = stream->gsSettings;
    const char *gridType;
    switch(gss.radialGridScheme) {
    case Dft::GridParams::GC2: gridType = "GC2"; break;
    case Dft::GridParams::LMG: gridType = "LMG"; break;
    case Dft::GridParams::TURBO: gridType = "Turbo"; break;
    default: gridType = "Invalid"; break;
    }
    do_output(LOG_CAT_INFO, LOG_AREA_DFT,
              "Generating %s grid Radint = %g AngInt: [ %d %d ]",
	      gridType, (double)gss.radint, gss.angmin, gss.angmax);
    stream->saveToFile(fname, noOfThreads);
    res = stream->getPointCount();
  } catch(const char *s) {
    printf("Error generating the grid: %s\n", s);
    fprintf(stderr, "Error generating the grid: %s\n", s);
    abort();
  }
  tm.print(LOG_AREA_DFT, __func__);
  return res;
}

void
grid_stream_free(ErgoGridStream *stream)
{
  delete stream;
}
