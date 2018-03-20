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

/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/

/** @file dft_common.cc

    @brief Common DFT routines. Mostly functional mixing.

    @author: Pawel Salek <em>responsible</em>

   (c) Pawel Salek, pawsa@theochem.kth.se, 2001-08-02
   NOTES: Adding new functionals:
   a. use fun-slater.c as template.
   b. add 'extern Functional MyFunctional;' to functionals.h
   c. add '&MyFunctional' to available_functionals below.
   d. have a beer. Or some crackers, if you prefer.
*/

#include <ctype.h>
#include <cmath>
#include <pthread.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <list>
#include <map>

#define __CVERSION__

#include "dft_common.h"
#include "functionals.h"
#include "grid_atomic.h"
#include "output.h"

/* C-wide constants */
int  ZEROI = 0,   ONEI = 1, THREEI = 3, FOURI = 4;
real ZEROR = 0.0, ONER = 1.0, TWOR = 2.0, FOURR = 4.0;

void*
dal_malloc_(size_t sz, const char *place, unsigned line)
{
    void* res = malloc(sz);
    if(!res) {
        fprintf(stderr, "dal_malloc(sz=%lu bytes) at %s (line %u) failed.\n",
                (unsigned long)sz, place, line);
        exit(1);
    }
    return res;
}

#if 0
static const ergo_real GET_BLOCKS_FUDGE_FACTOR = 0.5*template_blas_sqrt(3.0);
static const ergo_real SET_SHELL_RADII_ORBITAL_THR = 2e-12;
#else
static const ergo_real GET_BLOCKS_FUDGE_FACTOR = 1.0;
static const ergo_real SET_SHELL_RADII_ORBITAL_THR = 1e-11;
#endif

/** Returns the shortest distance of the border of the box to the
    specified point in space. */
real
Box::getDistanceTo(const real* v) const
{
  /* Identify the contact point and compute the distance to
     it... Observe how nicely the code handles the case of the point
     inside the box! */
  Vector3D contactPoint;

  for(int i=0; i<3; i++) {
    if(v[i] <= lo[i]) {
      contactPoint[i] = lo[i];
    } else if(v[i] >= hi[i]) {
      contactPoint[i] = hi[i];
    } else contactPoint[i] = v[i];
  }
  real dist= contactPoint.dist(v);
  return dist;
}

/** Return the index of the largest Cartesian dimension: 0 for x, 1
    for y and 2 for z. */
int
Box::getMaxDim() const
{
  int maxDim = 0;
  real d = hi[0]-lo[0];
  for(int j=1; j<3; j++) {
    real x = hi[j] - lo[j];
    if(x>d) {
      d = x;
      maxDim = j;
    }
  }
  return maxDim;
}


#include "barrier.h"

/** creates or destroys a barrier for nThreads. 
    
@param release tells whether we are to destroy the barrier (true) or
just sync (false).

@param nThreads informs the code how many threads are supposed to
block on the barrier.
 */
int
sync_threads(bool release, int nThreads)
{
    static pthread_mutex_t barrierInitMutex = PTHREAD_MUTEX_INITIALIZER;
    static ergo_barrier_t barrier;
    static bool barrierInitialized = false;
    if (pthread_mutex_lock(&barrierInitMutex) != 0) return -1;
    if(release) {
        if(barrierInitialized) {
            ergo_barrier_destroy(&barrier);
            barrierInitialized = false;
        }
        if (pthread_mutex_unlock(&barrierInitMutex) != 0) return -1;
    } else {
        int rc = 0;
        if (!barrierInitialized) {
            rc = ergo_barrier_init (&barrier, NULL, nThreads);
            barrierInitialized = true;
        }
        if (pthread_mutex_unlock(&barrierInitMutex) != 0) return -1;
        if(rc) return -1;
        int ret = ergo_barrier_wait(&barrier);
        if ((ret != 0) && (ret != PTHREAD_BARRIER_SERIAL_THREAD)) {
            fprintf (stderr, "Error when waiting for barrier\n");
            return -1;
        }
    }
    return 0;
}

/* =================================================================== */
/* dftinput:

   read DFT functional from given line. The calling convention assumes
   Sun-style passing parameters. ATLAS linear algebra package
   http://www.netlib.org/atlas/ or http://math-atlas.sourceforge.net/
   contains an elaborate discuttion of character type variable passing
   conventions, in particular little bit below
   http://math-atlas.sourceforge.net/errata.html#RH7.0
*/
//static char* DftConfString = NULL;
static real dft_hf_weight = 0.0;
static void dft_set_hf_weight(real w) { dft_hf_weight = w; }
static real dft_get_hf_weight(void)   { return dft_hf_weight; }
int (*fort_print)(const char* format, ...) = printf;

static int
ergo_fort_print(const char *format, ...)
{
  int res;
  va_list v;

  va_start(v, format);
  res = do_voutput(LOG_CAT_INFO, LOG_AREA_DFT, format, v);
  va_end(v);
  return res;
}

EXTERN_C void
dft_init(void)
{
  fort_print = ergo_fort_print;
  fun_printf = ergo_fort_print;
}

static int dft_thread_count = -1;
EXTERN_C int
dft_get_num_threads(void)
{ 
    if(dft_thread_count<0) { /* initialization needed */
        char *env = getenv("OMP_NUM_THREADS");
        if(env) {
	    dft_thread_count=atoi(env);
        } else { /* try opening the /proc/cpuinfo file */
            FILE *f = fopen("/proc/cpuinfo", "rt");
            if(f) {
                char line[256];
                dft_thread_count = 0;
                while(fgets(line, sizeof(line), f))
                    if(strncmp(line, "processor", 9) == 0)
                        dft_thread_count++;
                fclose(f);
            }
        }
        if(dft_thread_count<=0) dft_thread_count = 1;
        fort_print("Default number of threads: %d\n", dft_thread_count);
    }
    return dft_thread_count;
}
EXTERN_C void
dft_set_num_threads(int nThreads)
{ 
  dft_thread_count = nThreads;
}

#if defined(_OLD_FORTRAN_COMPATIBLE_CODE_)
int
FSYM(dftsetfunc)(const char* line, int * inperr, int len)
{
    int i, off;

#if 0
    fun_printf        = fort_print;
    fun_set_hf_weight = dal_set_hf_weight;
    fun_get_hf_weight = dal_get_hf_weight;
    fun_set_cam_param = dal_set_cam_param;
#else
    fun_set_hf_weight = dft_set_hf_weight;
    fun_get_hf_weight = dft_get_hf_weight;
#endif

    for(i=len-1; i>=0 && isspace((int)line[i]); i--)
        ;
    if(DftConfString) free(DftConfString);
    i++;
    for(off=0; line[off] && isspace((int)line[off]); off++)
        ;
    DftConfString = (char*)malloc(i+1-off);
    strncpy(DftConfString, line+off, i-off); 
    DftConfString[i-off] = '\0';
    
    switch(fun_select_by_name(DftConfString)) {
    case FUN_OK: return 1; /* SUCCESS! */
    case FUN_UNKNOWN:
        fun_printf("Unknown functional '%s'. Aborting.\n", DftConfString);
        dftlistfuncs_();
        break;
    case FUN_CONF_ERROR:
        fun_printf("Functional configuration '%s' is not understood. "
                   "Aborting.\n", DftConfString);
        break;
    }
    (*inperr)++;
    return 0; /* failed to find */
}
#endif /* _OLD_FORTRAN_COMPATIBLE_CODE_ */

/* dft_setfunc:
   returns 1 on success and 0 on failure (unknown functional, etc.)
*/
int
dft_setfunc(const char *line)
{
#if 0
    fun_printf        = fort_print;
    fun_set_hf_weight = dal_set_hf_weight;
    fun_get_hf_weight = dal_get_hf_weight;
    fun_set_cam_param = dal_set_cam_param;
#else
    fun_set_hf_weight = dft_set_hf_weight;
    fun_get_hf_weight = dft_get_hf_weight;
#endif

    switch(fun_select_by_name(line)) {
    case FUN_OK: return 1; /* SUCCESS! */
    case FUN_UNKNOWN:
        fun_printf("Unknown functional '%s'. Aborting.\n", line);
        dftlistfuncs_();
        break;
    case FUN_CONF_ERROR:
        fun_printf("Functional configuration '%s' is not understood. "
                   "Aborting.\n", line);
        break;
    }
    return 0; /* failed to find */
}

/** Return atom data. */
void
ErgoMolInfo::getAtom(int icent, int *cnt, real (*coor)[3],
                     int *charge, int *mult) const
{
    *cnt = 1; /* one atom in the batch */
    coor[0][0] = molecule.getAtom(icent).coords[0];
    coor[0][1] = molecule.getAtom(icent).coords[1];
    coor[0][2] = molecule.getAtom(icent).coords[2];
    *charge    = int( molecule.getAtom(icent).charge );
    *mult = 1; /* no atom is multiplied by a symmetry factor */
}

void
ErgoMolInfo::setShellRadii(real *rshell) const
{
    static const real facl[] = {
        1.0, 1.33, 1.6, 1.83, 2.03, 2.22, 2.39, 2.55, 2.70, 2.84 };
    real thlog = template_blas_log(SET_SHELL_RADII_ORBITAL_THR);
    int ishela, iprim;

    for(ishela=0; ishela <bis.noOfShells; ishela++) {
        const ShellSpecStruct& currShell = bis.shellList[ishela];
        rshell[ishela] =  0.0;
        for(iprim=0; iprim<currShell.noOfContr; iprim++) {
            ergo_real absCoeff = template_blas_fabs(currShell.coeffList[iprim]);
            if(absCoeff>0) {
              real r2 = (template_blas_log(absCoeff)-thlog)/
                    currShell.exponentList[iprim];
                if(rshell[ishela]<r2) rshell[ishela] = r2;
            }
        }
    }
    for(ishela=0; ishela <bis.noOfShells; ishela++) {
        rshell[ishela] = template_blas_sqrt(rshell[ishela])
            *facl[bis.shellList[ishela].shellType];
    }
}

/**
  get blocks of active SHELLS in cube of CELLSZ size centered at CENTER.
  
  RSHELL - precomputed shell extents.
  NBLCNT (output) - number of active blocks
  IBLCKS (output) - pairs of (startindex, stopindex)

  This algorithm scales quadratically.
*/
void
ErgoMolInfo::getBlocks1(const real *center, real cellsz, const real *rshell,
                       int *nblcnt, int (*iblcks)[2]) const
{
#if 0
    *nblcnt = 1;
    iblcks[0][0] = 0;
    iblcks[0][1] = bis.noOfShells;
#else
    int iprev, ishela;
    real celldg, px, py, pz, dst;

    *nblcnt = -1;
    iprev = -1111;
    celldg = cellsz*0.5*template_blas_sqrt(3.0);
    for(ishela=0; ishela <bis.noOfShells; ishela++) {
        const ShellSpecStruct& currShell = bis.shellList[ishela];        
        px = center[0]-currShell.centerCoords[0];
        py = center[1]-currShell.centerCoords[1];
        pz = center[2]-currShell.centerCoords[2];
        dst = template_blas_sqrt(px*px + py*py + pz*pz);
        if(dst<rshell[ishela]+celldg) {
            /* accepted... */
            if(ishela != iprev+1)
                iblcks[++(*nblcnt)][0] = ishela;
            iprev = ishela;
            iblcks[*nblcnt][1] = ishela+1;
        }
    }
    (*nblcnt)++;
#endif
}

/** Class that allows to find in NLogN time all shells that overlap
    with a given box. */
class ShellTree {
    const ShellSpecStruct *shells;
    ShellTree *smaller;
    ShellTree *larger;
    ergo_real dividingValue;
    int       dividingDimension;
    std::map<int,ergo_real> ownShells; /**<set only for leaves,
                                        *  i.e. such ShellTree objects
                                        *  that have smaller and
                                        *  larger fields set to
                                        *  NULL. Tree node contains a
                                        *  shell*/
    ergo_real maxRadius; /**< upper limit of ownShell radius. */
 public:
    ShellTree(const BasisInfoStruct& bis_,
              const real *rShells_);
    ShellTree(const BasisInfoStruct& bis_,
              const real *rshell,
              const std::list<int>& activeShells,
              const Box& bb);
    ~ShellTree() {
        if(smaller)
            delete smaller;
        if(larger)
            delete larger;
    }
    void init(const BasisInfoStruct& bis_,
              const real *rShells_,
              const std::list<int>& activeShells,
              const Box& bb);
    void getOverlappingWith(const real *center, real cellsz,
                            std::map<int,int>& res) const;
};

struct Ball {
    ergo_real rad;
    ergo_real center[3]; 
    ergo_real radius() const { return rad; }
    Ball() {}
    Ball(const ergo_real (&c)[3], ergo_real r)
        : rad(r) {
        center[0] = c[0]; center[1] = c[1]; center[2] = c[2];
    }
};

/** root node constructor. It does some initalization work and passes
    on the ball to the leave constructors. */
ShellTree::ShellTree(const BasisInfoStruct& bis_,
                     const real *rShells)
    : shells(bis_.shellList), smaller(NULL), larger(NULL), maxRadius(0.0)
{
    Box box;

    std::vector<Ball> shellBalls(bis_.noOfShells);
    for(int i=0; i<bis_.noOfShells; i++)
        shellBalls[i] = Ball(bis_.shellList[i].centerCoords, 0.1);
    getBoundingBox(box, shellBalls.begin(), shellBalls.end());
    
    std::list<int> shells;
    for(int i=bis_.noOfShells-1; i>=0; i--)
        shells.push_front(i);

    init(bis_, rShells, shells, box);
}

/** Constructs the ShellTree. */
ShellTree::ShellTree(const BasisInfoStruct& bis_,
                     const real *rShells,
                     const std::list<int>& activeShells,
                     const Box& bb)
    : shells(bis_.shellList), smaller(NULL), larger(NULL), maxRadius(0.0)
{
    init(bis_, rShells, activeShells, bb);
}

void
ShellTree::init(const BasisInfoStruct& bis_,
                const real *rShells,
                const std::list<int>& activeShells,
                const Box& bb)
{
    /* When you go down with box size to this value, stop recursion,
       there is too little to win. */
    static const ergo_real MIN_BOX_SIZE = 1.0;
#if defined(DEBUG)
    printf("Puts %d shells in box (%f,%f,%f)-(%f,%f,%f)\n",
           activeShells.size(), bb.lo[0], bb.lo[1], bb.lo[2],
           bb.hi[0], bb.hi[1], bb.hi[2]);
    for(std::list<int>::const_iterator i= activeShells.begin();
        i != activeShells.end(); ++i)
        printf("%d ", *i);
    puts("");
#endif
    dividingDimension = bb.getMaxDim();
    dividingValue = 0.5*(bb.hi[dividingDimension]+bb.lo[dividingDimension]);

    std::list<int> lessThanList, greaterList;

    assert(!activeShells.empty());
    if(activeShells.size() == 1 ||
       bb.hi[dividingDimension]-bb.lo[dividingDimension] < MIN_BOX_SIZE) {
        for(std::list<int>::const_iterator i = activeShells.begin();
            i != activeShells.end(); ++i) {
            ownShells[*i] = rShells[*i];
            if(rShells[*i]>maxRadius)
                maxRadius = rShells[*i];
        }
        return;
    }

    /* Proceed with the recursion... */
    for(std::list<int>::const_iterator i = activeShells.begin();
        i != activeShells.end(); ++i) {
        int shell = *i;
        /* Do the job here! */
        if(shells[shell].centerCoords[dividingDimension]>dividingValue)
            greaterList.push_back(shell);
        else
            lessThanList.push_back(shell);
    }
    maxRadius = 0;
    if(!lessThanList.empty()) {
        Box box = bb; box.hi[dividingDimension] = dividingValue;
        smaller = new ShellTree(bis_, rShells, lessThanList, box);
        if(smaller->maxRadius > maxRadius)
            maxRadius = smaller->maxRadius;
    }
    if(!greaterList.empty()) {
        Box box = bb; box.lo[dividingDimension] = dividingValue;
        larger = new ShellTree(bis_, rShells, greaterList, box);
        if(larger->maxRadius > maxRadius)
            maxRadius = larger->maxRadius;
    }
}


/** Coomputes distance between two points, they do not need to be of
    the Vector3D type. In a way, we lose here some error
    checking. Perhaps it should be avoided. */
static inline ergo_real
distance(const ergo_real *a, const ergo_real *b)
{
    ergo_real r2 = 0;
    for(int i=0; i<3; i++) {
        ergo_real d = a[i]-b[i];
        r2 += d*d;
    }
    return template_blas_sqrt(r2);
}

void
ShellTree::getOverlappingWith(const real *center, real cellsz,
                              std::map<int,int>& res) const
{
#if defined(DEBUG)
    static const char XYZ[] = "XYZ";
    printf("Checking whether the box at %f %f %f (%f) is close to %c %f\n",
           center[0], center[1], center[2], cellsz,
           XYZ[dividingDimension], dividingValue);
#endif
    if(smaller &&
       center[dividingDimension]-cellsz < dividingValue +smaller->maxRadius)
        smaller->getOverlappingWith(center, cellsz, res);
    
    if(larger &&
       center[dividingDimension]+cellsz > dividingValue - larger->maxRadius)
        larger->getOverlappingWith(center, cellsz, res);

    /* "This" node contains at least one shell overlapping with given
       node. We can however tighten up the selection using specified
       shell radii. */
    for(std::map<int,ergo_real>::const_iterator i= ownShells.begin();
        i != ownShells.end(); ++i) {
        int shell = i->first;
        ergo_real r = i->second;
        ergo_real dist = distance(shells[shell].centerCoords, center);
#if defined(DEBUG)
        printf("Shell %d dist %f radius: %f + %f = %f\n",
               shell, dist, r, cellsz, r + cellsz);
#endif
        if(dist <= r + cellsz)
            res[shell] = shell;
    }
}


/** Ther standard constructor. */
ErgoMolInfo::ErgoMolInfo(const BasisInfoStruct& bis_, const Molecule& mol)
  : GridGenMolInfo(mol.getNoOfAtoms(),bis_.noOfBasisFuncs, bis_.noOfShells),
      bis(bis_), molecule(mol), shellTree(NULL)
{ 
    ergo_real *rShell = new ergo_real[bis.noOfShells];
    setShellRadii(rShell);
    shellTree = new ShellTree(bis_, rShell);
    delete []rShell;
}

ErgoMolInfo::~ErgoMolInfo()
{
    if(shellTree)
        delete shellTree;
}

/** same as ergo_get_shlblocks, except it should scale NlogN.
    rshell is not used - we store this information in the tree.
 */
void
ErgoMolInfo::getBlocks(const real *center, real cellsz, const real *rshell,
                        int *nblcnt, int (*iblcks)[2]) const
{
    std::map<int,int> activeShellList;
    /* The factor 0.5*template_blas_sqrt(3.0) does not fit here but is present
       for compatibility with the old code. */
    shellTree->getOverlappingWith(center, cellsz *GET_BLOCKS_FUDGE_FACTOR,
                                  activeShellList);
    
    /* sort list and create the packed structure nblcnt, iblcks */
    /* Not needed for maps... activeShellList.sort(); */
    *nblcnt = -1;
    int iprev = -2;
#if defined(DEBUG)
    printf("Following shell blocks were found to be active around\n"
           "%f %f %f of radius %f:\n",
           center[0], center[1], center[2], cellsz);
#endif
    for(std::map<int,int>::const_iterator i = activeShellList.begin();
        i != activeShellList.end(); ++i) {
#if defined(DEBUG)
        printf("%d (%d)\n", i->first, i->second);
#endif
        if (i->first != iprev+1) /* if noncontinous */
            iblcks[++(*nblcnt)][0] = i->first;
        iprev = i->first;
        iblcks[*nblcnt][1] = i->first+1;
    }
    (*nblcnt)++;

#if defined(DEBUG)
    int (*refBlocks)[2] = new int[noOfShells][2];
    int cnt;
    getBlocks1(center, cellsz, rshell, &cnt, refBlocks);
    if(*nblcnt != cnt) {
	printf("Compressed:\n");
	for(int i=0; i<*nblcnt; i++)
	    printf("%d:%d\n", iblcks[i][0], iblcks[i][1]);
	printf("Correct ones, compressed:\n");
	for(int i=0; i<cnt; i++)
	    printf("%d:%d\n", refBlocks[i][0], refBlocks[i][1]);
    }
    delete []refBlocks;
#endif
}

/** ergo_get_exps() generates a list of exponents for every center as
   in mol_info table: number of gaussians at given center (nucorb),
   and their smallest and largest exponent in aaa. 
   @param maxl max L quantum number - the leading dimension.
   @param bascnt [noOfAtoms][maxl] - no of funcs at center and
   given angular momentum.
   @param aa  [noOfAtoms][maxl][2] min/max exponent.
*/
void
ErgoMolInfo::getExps(int *maxl, int **bascnt, real (**aa)[2]) const
{
  int icent, ishell, il, lda, noOfAtoms = molecule.getNoOfAtoms();
    real (*la)[2];
    const Molecule& ms = molecule;
    
    lda = 0;
    for(ishell=0; ishell<bis.noOfShells; ishell++) {
        int l = bis.shellList[ishell].shellType;
        if(l > lda) lda = l;
    }
    ++lda;
    *maxl = lda;
    *bascnt = (int*)calloc(  lda*noOfAtoms, sizeof(int));
    la = *aa     = (real (*)[2])calloc(2*lda*noOfAtoms, sizeof(real));

    /* init min exponents to largest possible exponent */
    for(icent=0; icent<noOfAtoms; icent++)
        for(il=0; il<lda; il++)
            la[lda*icent + il][0] = 1e30; 


    for(ishell=0; ishell<bis.noOfShells; ishell++) {
        const ShellSpecStruct& shell = bis.shellList[ishell];
        int iexp, lval = shell.shellType;

	/* Go through all atoms to check if this shell is sitting 
	   on an atom center */
        for(icent=0; icent<noOfAtoms; icent++) {
	  real dx = shell.centerCoords[0]-ms.getAtom(icent).coords[0];
	  real dy = shell.centerCoords[1]-ms.getAtom(icent).coords[1];
	  real dz = shell.centerCoords[2]-ms.getAtom(icent).coords[2];
	  // FIXME hard-coded constant 1e-10 here, is that bad?
	  if(dx*dx + dy*dy+dz*dz < 1e-10)
	    break;
        }
	/* If (icent == noOfAtoms) that means this shell is not on any 
	   atom (it is a "ghost" shell).
	   In that case we ant to continue with the next shell. */
	if(icent == noOfAtoms)
	    continue;
	/* There was an assert(icent<noOfAtoms) here but that causes
	   problems when ghost atoms are used.  Elias removed the
	   assertion. Pawel thinks the change is OK. */
	
        (*bascnt)[lval + lda*icent]++;
        for(iexp = 0; iexp<shell.noOfContr; iexp++) {
            int idx = lda*icent+ lval;
            real an = shell.exponentList[iexp];
            if(an < la[idx][0]) la[idx][0] = an;
            if(an > la[idx][1]) la[idx][1] = an;
        }
    }
}

/** transform shell block indices to orbital block indices.  IORIDX
  contains preprocessed information about where given shell begins and
  ends in given symmetry.
*/
void
ergoShellsToOrbs(const int *nshlbl, const int (*shlblock)[2],
                 int *norbbl, int (*orbblock)[2],
                 const BasisInfoStruct& bis)
/*
  INTEGER SHLBLOCK(2,NSHLBL), ORBBLOCK(2,NSHLBL,NSYM)
  INTEGER IORIDX(2,KMAX,NSYM), NORBBL(NSYM)
*/
{
    int ishbl;

    norbbl[0] = 0;
    for(ishbl=0; ishbl< *nshlbl; ishbl++) {
        const ShellSpecStruct& shellLo = bis.shellList[shlblock[ishbl][0]];
        const ShellSpecStruct& shellHi = bis.shellList[shlblock[ishbl][1]-1];
        orbblock[norbbl[0]][0] = shellLo.startIndexInMatrix;
        orbblock[norbbl[0]][1] = shellHi.startIndexInMatrix+
            shellHi.noOfBasisFuncs;
#if 0
        printf("shell [%d %d] to orb [%d %d]\n",
               shlblock[ishbl][0], shlblock[ishbl][1]-1,
               orbblock[norbbl[0]][0],
               orbblock[norbbl[0]][1]);
#endif
        norbbl[0]++;               
    }
}


/* =================================================================== */
/* transformations of functional derivatives from from unrestricted
 * ones to closed shell case. */
real
dftene_(const real *rho, const real *grad)
{
    FunDensProp dp;
    dp.rhoa  = dp.rhob  = *rho  *0.5;
    dp.grada = dp.gradb = *grad *0.5;
    dp.gradab = dp.grada*dp.gradb;
    return selected_func->func(&dp);
}

void
dftptf0_(real *rho, real *grad, real *wght, real *vx)
{
    FunFirstFuncDrv drvs;
    FunDensProp dp;
    dp.rhoa  = dp.rhob  = *rho *0.5;
    dp.grada = dp.gradb = *grad*0.5;
    if(dp.rhoa<1e-13) dp.rhoa = dp.rhob = 1e-13;
    if(dp.grada<1e-13) dp.grada = dp.gradb = 1e-13;
    dp.gradab = dp.grada*dp.gradb;
    drv1_clear(&drvs);
    selected_func->first(&drvs, *wght, &dp);
    vx[0] = drvs.df1000;
    vx[1] = drvs.df0010 + 0.5*drvs.df00001* (*grad);
}

void
dftpot0_(FirstDrv *ds, const real* weight, const FunDensProp* dp)
{
    FunFirstFuncDrv drvs;
    drv1_clear(&drvs);
    selected_func->first(&drvs, *weight, dp);
    ds->fR = drvs.df1000;
    ds->fZ = drvs.df0010 + 0.5*drvs.df00001* (dp->grada+dp->gradb);
}

void
dftpot1_(SecondDrv *ds, const real* w, const FunDensProp* dp, 
         const int* triplet)
{
    
    FunSecondFuncDrv drvs;

    drv2_clear(&drvs);
    if(dp->rhoa + dp->rhob>1e-14)
        selected_func->second(&drvs, *w, dp);
    if (*triplet) { /* triplet */  
        ds->fZ  = drvs.df0010;
        ds->fG  = -0.5*drvs.df00001;
        ds->fRR = 0.5*(drvs.df2000 - drvs.df1100);
        ds->fRZ = 0.5*(drvs.df1010 - drvs.df1001); 
        ds->fRG = 0.0;  
        ds->fZZ = 0.5*(drvs.df0020 - drvs.df0011); 
        ds->fZG = 0.0; 
        ds->fGG = 0.0; 
    } else { /* singlet */
        ds->fR  = 0.5*(drvs.df1000 + drvs.df0100);
        ds->fZ  = drvs.df0010;
        ds->fRR = 0.5*(drvs.df2000 + drvs.df1100);
        ds->fRZ = 0.5*(drvs.df1010 + drvs.df1001);
        ds->fZZ = 0.5*(drvs.df0020 + drvs.df0011); 
        ds->fRG = 0.5*drvs.df10001;   
        ds->fZG = 0.5*drvs.df00101;   
        ds->fGG = 0.25*drvs.df00002; 
        ds->fG  = 0.5*drvs.df00001;  
    }
}


#if defined(_OLD_FORTRAN_COMPATIBLE_CODE_)

/* dftpot2_:
   computes third order derivatives of selected functional with respect
   to rho and zeta=|\nabla\rho|
*/
void
dftpot2_(ThirdDrv *ds, real factor, const FunDensProp* dp, int isgga,
         int triplet)
{
    FunThirdFuncDrv drvs;

    drv3_clear(&drvs);
    selected_func->third(&drvs, factor, dp);
    /* transform to density from density_alpha derivatives below */
    /* This could be a separate function. */
    /* we treat singlet here */
    ds->fR   = (drvs.df1000);
    ds->fRR[0] = (drvs.df2000 + drvs.df1100);
    ds->fRR[1] = (drvs.df2000 - drvs.df1100);
    ds->fRRR[0] = (drvs.df3000 + 3*drvs.df2100);
    ds->fRRR[1] = (drvs.df3000 - drvs.df2100);

    if(isgga) { /* FORMULAE seem ok but were never really tested */
        real grada = dp->grada;
        real grada2= grada*grada;
        real grada3= grada2*grada;
	/* transform to closed shell, second time. Oh, I love this mess. */
        ds->fZ  = drvs.df0010/(2*grada);
        ds->fG = drvs.df00001;
        ds->fZZ[0]  = (drvs.df0020 + drvs.df0011)/(4*grada2) 
	    -drvs.df0010/(4*grada3);
        ds->fZZ[1]  = (drvs.df0020 - drvs.df0011)/(4*grada2) 
	    -drvs.df0010/(4*grada3);
        ds->fRZ[0]  = (drvs.df1010 + drvs.df1001)/(2*grada);
        ds->fRZ[1]  = (drvs.df1010 - drvs.df1001)/(2*grada);
        ds->fRG[0]  = 2*drvs.df10001;
        ds->fRG[1]  = 0.0;  
        ds->fRRZ[0][0] = (drvs.df2010+drvs.df2001+2*drvs.df1110)/(2*grada);
        ds->fRRZ[0][1] = (drvs.df2010+drvs.df2001-2*drvs.df1110)/(2*grada);
        ds->fRRZ[1][0] = (drvs.df2010-drvs.df2001)/(2*grada);
        ds->fRRZ[1][1] = ds->fRRZ[1][0];
        ds->fRRG[0] = drvs.df20001+drvs.df11001;
        ds->fRRG[1] = drvs.df20001-drvs.df11001;
        ds->fRRGX[0][0] = 2*(drvs.df20001+drvs.df11001);
        ds->fRRGX[1][1] = 2*(drvs.df20001-drvs.df11001);
        ds->fRRGX[1][0] = ds->fRRGX[0][1] = 0;
         
        ds->fRZZ[0][0] = (drvs.df1020+drvs.df0120+2*drvs.df1011)/(4*grada2) 
                      - (drvs.df1010+drvs.df1001)/(4*grada3);
        ds->fRZZ[0][1] = (drvs.df1020+drvs.df0120-2*drvs.df1011)/(4*grada2) 
                      - (drvs.df1010+drvs.df1001)/(4*grada3);
        ds->fRZZ[1][0] = (drvs.df1020-drvs.df0120)/(4*grada2) 
                      - (drvs.df1010-drvs.df1001)/(4*grada3);
        ds->fRZZ[1][1] = ds->fRZZ[1][0];
        ds->fZZZ[0] = ((drvs.df0030 + 3*drvs.df0021)/grada3 
                   -3*(drvs.df0020 + drvs.df0011)/(grada2*grada2)
                   +3*drvs.df0010/(grada3*grada2))/8.0; 
        ds->fZZZ[1] = ((drvs.df0030 - drvs.df0021)/grada3 
                   -(3*drvs.df0020 - drvs.df0011)/(grada2*grada2)
                   +3*drvs.df0010/(grada3*grada2))/8.0; 
    } else {
        ds->fZ = ds->fZZ[0] = ds->fZZ[1] = ds->fRZ[0] = ds->fRZ[1] = 0;
        ds->fRRZ[0][0] = ds->fRRZ[0][1]  = ds->fRRZ[1][0] = 0;
        ds->fRRZ[1][1] = ds->fRZZ[0][0] = ds->fRZZ[0][1] = 0;
        ds->fRZZ[1][0] = ds->fRZZ[1][1] = ds->fZZZ[0] = ds->fZZZ[1] = 0; 
        ds->fG = ds->fRG[0] = ds->fRG[1]= ds->fRRG[0] = ds->fRRG[1] = 0;
        ds->fRRGX[0][0] = ds->fRRGX[1][0] = 
	    ds->fRRGX[0][1] = ds->fRRGX[1][1] = 0;  
    }
}

void
dftpot3ab_(FourthDrv *ds, real *factor, const FunDensProp* dp, int *isgga)
{
     FunFourthFuncDrv drvs;
     
     /* Initialize drvs */ 
     drv4_clear(&drvs);
     selected_func->fourth(&drvs, *factor, dp);

     /* Transform derivatives (drvs -> ds) */

     ds->fR = 0.5*(drvs.df1000 + drvs.df0100);
     ds->fRR = 0.25*(drvs.df2000 + 2.0*drvs.df1100 + drvs.df0200);
     ds->fRRR = 0.125*(drvs.df3000 + 3.0*drvs.df2100 + 
                       3.0*drvs.df1200 + drvs.df0300);
     ds->fRRRR = 0.0625*(drvs.df4000 + 4.0*drvs.df3100 + 
                         6.0*drvs.df2200 + 4.0*drvs.df1300 + drvs.df0400);
     
     if (*isgga) {
         real igroa, igrob;
         real igroa_p2, igrob_p2;
         real igroa_p3, igrob_p3;
         real igroa_p4, igrob_p4;
         real igroa_p5, igrob_p5;
         real igroa_p6, igrob_p6;
         real igroa_p7, igrob_p7;
         
         /* 1/groa, 1/grob and its powers */
         igroa = 1.0/(dp->grada);
         igrob = 1.0/(dp->gradb);
         
         igroa_p2=igroa*igroa;
         igroa_p3=igroa_p2*igroa;
         igroa_p4=igroa_p3*igroa;
         igroa_p5=igroa_p4*igroa;
         igroa_p6=igroa_p5*igroa;
         igroa_p7=igroa_p6*igroa;
         
         igrob_p2=igrob*igrob;
         igrob_p3=igrob_p2*igrob;
         igrob_p4=igrob_p3*igrob;
         igrob_p5=igrob_p4*igrob;
         igrob_p6=igrob_p5*igrob;
         igrob_p7=igrob_p6*igrob;
         
         
         /* 1st order */
         
         ds->fZ = 0.125*(igroa*drvs.df0010+2.0*drvs.df00001+igrob*drvs.df0001);
         
         /* 2nd order */
         
         ds->fRZ = 0.0625*(2*drvs.df01001 + 2*drvs.df10001 + drvs.df0110*igroa
                           +drvs.df1010*igroa + drvs.df0101*igrob
                           + drvs.df1001*igrob);
     
         ds->fZZ = 0.015625*(4*drvs.df00002 + 4*drvs.df00101*igroa + 
                     drvs.df0020*igroa_p2 - drvs.df0010*igroa_p3 + 
                     4*drvs.df00011*igrob + 2*drvs.df0011*igroa*igrob + 
                     drvs.df0002*igrob_p2 - drvs.df0001*igrob_p3);

      /* 3rd order */
     
     ds->fRRZ = 0.03125*(2*drvs.df02001 + 4*drvs.df11001 + 2*drvs.df20001 +
		     drvs.df0210*igroa + 2*drvs.df1110*igroa +
                     drvs.df2010*igroa + drvs.df0201*igrob +
                     2*drvs.df1101*igrob + drvs.df2001*igrob);

     ds->fRZZ = 0.0078125*(4*drvs.df01002 + 4*drvs.df10002 + 4*drvs.df01101*igroa + 
	             4*drvs.df10101*igroa + drvs.df0120*igroa_p2 +       
	             drvs.df1020*igroa_p2 - drvs.df0110*igroa_p3 -       
	             drvs.df1010*igroa_p3 + 4*drvs.df01011*igrob +       
	             4*drvs.df10011*igrob + 2*drvs.df0111*igroa*igrob +       
	             2*drvs.df1011*igroa*igrob + drvs.df0102*igrob_p2 +       
	             drvs.df1002*igrob_p2 - drvs.df0101*igrob_p3 -       
	             drvs.df1001*igrob_p3);

     ds->fZZZ = 0.001953125*(8*drvs.df00003 + 12*drvs.df00102*igroa + 
		     6*drvs.df00201*igroa_p2 - 6*drvs.df00101*igroa_p3 +       
		     drvs.df0030*igroa_p3 - 3*drvs.df0020*igroa_p4 +       
		     3*drvs.df0010*igroa_p5 + 12*drvs.df00012*igrob + 
		     12*drvs.df00111*igroa*igrob +3*drvs.df0021*igroa_p2*igrob-
		     3*drvs.df0011*igroa_p3*igrob + 
		     6*drvs.df00021*igrob_p2 + 3*drvs.df0012*igroa*igrob_p2 - 
		     6*drvs.df00011*igrob_p3 + drvs.df0003*igrob_p3 - 
		     3*drvs.df0011*igroa*igrob_p3 - 3*drvs.df0002*igrob_p4 + 
	  	     3*drvs.df0001*igrob_p5);

     /* 4th order */

     ds->fRRRZ = 0.015625*(2*drvs.df03001 + 6*drvs.df12001 + 6*drvs.df21001 + 
		     2*drvs.df30001 + drvs.df0310*igroa + 3*drvs.df1210*igroa +
		     3*drvs.df2110*igroa + drvs.df3010*igroa +
		     drvs.df0301*igrob + 3*drvs.df1201*igrob + 
		     3*drvs.df2101*igrob + drvs.df3001*igrob);
     
     ds->fRRZZ = 0.00390625*(4*drvs.df02002 + 8*drvs.df11002 + 4*drvs.df20002 +
	             4*drvs.df02101*igroa + 8*drvs.df11101*igroa + 
	             4*drvs.df20101*igroa + drvs.df0220*igroa_p2 + 
 	             2*drvs.df1120*igroa_p2 + drvs.df2020*igroa_p2 - 
	             drvs.df0210*igroa_p3 - 2*drvs.df1110*igroa_p3 - 
	             drvs.df2010*igroa_p3 + 4*drvs.df02011*igrob + 
	             8*drvs.df11011*igrob + 4*drvs.df20011*igrob + 
	             2*drvs.df0211*igroa*igrob + 4*drvs.df1111*igroa*igrob + 
	             2*drvs.df2011*igroa*igrob + drvs.df0202*igrob_p2 + 
	             2*drvs.df1102*igrob_p2 + drvs.df2002*igrob_p2 - 
	             drvs.df0201*igrob_p3 - 2*drvs.df1101*igrob_p3 - 
	             drvs.df2001*igrob_p3);

     ds->fRZZZ = 0.0009765625*(8*drvs.df01003 + 8*drvs.df10003 + 12*drvs.df01102*igroa + 
		     12*drvs.df10102*igroa + 6*drvs.df01201*igroa_p2 + 
		     6*drvs.df10201*igroa_p2 - 6*drvs.df01101*igroa_p3 + 
		     drvs.df0130*igroa_p3 - 6*drvs.df10101*igroa_p3 + 
		     drvs.df1030*igroa_p3 - 3*drvs.df0120*igroa_p4 - 
		     3*drvs.df1020*igroa_p4 + 3*drvs.df0110*igroa_p5 + 
		     3*drvs.df1010*igroa_p5 + 12*drvs.df01012*igrob + 
		     12*drvs.df10012*igrob + 12*drvs.df01111*igroa*igrob + 
		     12*drvs.df10111*igroa*igrob +3*drvs.df0121*igroa_p2*igrob+
		     3*drvs.df1021*igroa_p2*igrob - 
		     3*drvs.df0111*igroa_p3*igrob-3*drvs.df1011*igroa_p3*igrob+
		     6*drvs.df01021*igrob_p2 + 6*drvs.df10021*igrob_p2 + 
		     3*drvs.df0112*igroa*igrob_p2 + 3*drvs.df1012*igroa*igrob_p2 -
		     6*drvs.df01011*igrob_p3 + drvs.df0103*igrob_p3 - 
		     6*drvs.df10011*igrob_p3 + drvs.df1003*igrob_p3 - 
		     3*drvs.df0111*igroa*igrob_p3-3*drvs.df1011*igroa*igrob_p3-
		     3*drvs.df0102*igrob_p4 - 3*drvs.df1002*igrob_p4 + 
		     3*drvs.df0101*igrob_p5 + 3*drvs.df1001*igrob_p5);

     ds->fZZZZ = (1.0/4096.0)*(16*drvs.df00004 + 32*drvs.df00103*igroa + 
		     24*drvs.df00202*igroa_p2 - 24*drvs.df00102*igroa_p3 + 
		     8*drvs.df00301*igroa_p3 - 24*drvs.df00201*igroa_p4 + 
		     drvs.df0040*igroa_p4 + 24*drvs.df00101*igroa_p5 - 
		     6*drvs.df0030*igroa_p5 + 15*drvs.df0020*igroa_p6 - 
		     15*drvs.df0010*igroa_p7 + 32*drvs.df00013*igrob + 
		     48*drvs.df00112*igroa*igrob + 24*drvs.df00211*igroa_p2*igrob - 
		     24*drvs.df00111*igroa_p3*igrob +
		     4*drvs.df0031*igroa_p3*igrob - 12*drvs.df0021*igroa_p4*igrob + 
		     12*drvs.df0011*igroa_p5*igrob +
		     24*drvs.df00022*igrob_p2 + 24*drvs.df00121*igroa*igrob_p2+
		     6*drvs.df0022*igroa_p2*igrob_p2 - 6*drvs.df0012*igroa_p3*igrob_p2 -
		     24*drvs.df00012*igrob_p3 + 8*drvs.df00031*igrob_p3 - 
		     24*drvs.df00111*igroa*igrob_p3 + 
		     4*drvs.df0013*igroa*igrob_p3 - 6*drvs.df0021*igroa_p2*igrob_p3 +
		     6*drvs.df0011*igroa_p3*igrob_p3 - 
		     24*drvs.df00021*igrob_p4 + drvs.df0004*igrob_p4 - 
		     12*drvs.df0012*igroa*igrob_p4 + 24*drvs.df00011*igrob_p5 -
		     6*drvs.df0003*igrob_p5 + 12*drvs.df0011*igroa*igrob_p5 + 
		     15*drvs.df0002*igrob_p6 - 15*drvs.df0001*igrob_p7);
     }
}

#if 0
void
daxpy_(const int* cnt, const real* alpha, const real* v1, 
       const int* stride1, real* v2, const int* stride2)
{
    int i;
    real a = *alpha;
    for(i= *cnt - 1; i>=0; i--)
        v2[i * *stride2] += a*v1[i* *stride1];
}
#endif

void
dzero_(real* arr, const int* len)
{
    int i;
    for(i= *len -1; i>=0; i--)
        arr[i] = 0;
}

#endif /* _OLD_FORTRAN_COMPATIBLE_CODE_ */
