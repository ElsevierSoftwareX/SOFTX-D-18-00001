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

/** @file slr.h

    @brief Contains a Simple Linear Response implementation based on
    the MO orbital picture.

    @author: Pawel Salek <em>responsible</em>
*/

#if !defined(SLR_HEADER)
#define SLR_HEADER
/* Copyright(c) Pawel Salek 2006. */

#include <unistd.h>

#include "realtype.h"

namespace LR {
class VarVector;
/** template based proxy object that uses bool-valued policies to
    perform the assignments.  */
template<bool MultByS,bool SwapXY>
  class VarVectorProxyOp {
 public:
  const VarVector& vec;
  ergo_real scalar;
  explicit VarVectorProxyOp(const VarVector& a, ergo_real s=1.0) : vec(a), scalar(s) {}
};


/** Vector of variables parametrising the solution to the linear
    response equations.  It provides some handy operations for
    compact notation. */
class VarVector {
  ergo_real *v;  /**< vector coefficients */
 public:
  char *fName;    /**< Present in secondary storage (i.e. disk) under
                   * given file name. */
  int  nvar;      /**< nvar := nocc*nvirt. Vector length is 2*nvar */
  unsigned onDisk:1;   /**< valid representation on disk */
  unsigned inMemory:1; /**< valid representation in memory */

  VarVector(const VarVector& a) : fName(NULL), nvar(a.nvar),
    onDisk(0), inMemory(1) {
    if(nvar) {
      v = new ergo_real[nvar*2];
      memcpy(v, a.v, 2*nvar*sizeof(ergo_real));
    } else v = NULL;
  }

  VarVector() : v(NULL), fName(NULL), nvar(0), onDisk(0), inMemory(1) {}

  /** Creates a vector from a full matrix. */
  VarVector(int nbast, int nocc, const ergo_real *full_mat)
    : v(NULL), fName(NULL), nvar(0), onDisk(0), inMemory(1)  {
    setFromFull(nbast, nocc, full_mat);
  }

  ~VarVector() {
    if(v) delete v; 
    if(fName) {
      unlink(fName);
      delete fName;
    }
}
  ergo_real* x() const { return v; } /**< return the X part of the vector. */
  ergo_real* y() const { return v + nvar; } /**< returns the Y part. */
  void symorth(void);
  void setSize(int n) { /**< sets the size. Reallocates if needed. */
    if(nvar != n) {
      if(v)
	delete v;
      v = new ergo_real[2*n];
      nvar = n;
      onDisk = 0;
    }
  }
  int size() const { return nvar; }
  void print(const char *comment = NULL) {
    if(comment) puts(comment);
    for(int i=0; i<nvar*2; i++) printf("%15.10f\n", (double)v[i]);
  }
  void setFromFull(int nbast, int nocc, const ergo_real *full_mat);
  void setFull(int nbast, int nocc, ergo_real *full_mat) const;
  const ergo_real& operator[](int i) const { return v[i]; }
  ergo_real& operator[](int i)             { return v[i]; }
  VarVector& operator=(ergo_real scalar) {
    for(int i=0; i<2*nvar; i++) v[i] = scalar;
    onDisk = 0;
    return *this;
  };
  VarVector& operator=(const VarVector& b) {
    if(this == &b) return *this;
    if(nvar != b.nvar) { 
      nvar = b.nvar; 
      if(v) delete v;
      v = new ergo_real[2*nvar];
    }
    memcpy(v, b.v, 2*nvar*sizeof(v[0]));
    onDisk = 0;
    return *this;
  }

    VarVector&
    operator=(const VarVectorProxyOp<false, false >& proxy) {
      if (&proxy.vec == this)
        throw "VarVector self-assignment not-implemented";
      if(nvar != proxy.vec.nvar) {
        if(v) delete v; 
        nvar = proxy.vec.nvar;
        v = new ergo_real[2*nvar];
      }

      for(int i=0; i<2*nvar; i++)
        v[i]      = proxy.scalar*proxy.vec[i];
      onDisk = 0;
      return *this;
    }
    VarVector&
    operator=(const VarVectorProxyOp<false, true >& proxy) {
      if (&proxy.vec == this)
        throw "VarVector self-assignment not-implemented";
      if(nvar != proxy.vec.nvar) {
        if(v) delete v; 
        nvar = proxy.vec.nvar;
        v = new ergo_real[2*nvar];
      }
      for(int i=0; i<nvar; i++) {
        v[i]      = proxy.scalar*proxy.vec[i+nvar];
        v[i+nvar] = proxy.scalar*proxy.vec[i];
      }
      onDisk = 0;
      return *this;
    }
    
    /** Load the object to memory. */
    void load(const char* tmpdir);
    /** Save the object. Probably does not need be done explicitly. */
    void save(const char* tmpdir);
    /** Releases the memory, saving if necessary. */
    void release(const char* tmpdir);

    friend ergo_real operator*(const VarVector& a, const VarVector& b);
    friend ergo_real operator*(const VarVector &a,
                               const VarVectorProxyOp<false,false>& b);
    friend ergo_real operator*(const VarVector &a,
                               const VarVectorProxyOp<true,false>& b);
};

/** E2Evaluator interface provides a way to perform a linear
    transformation of supplied transition density matrix @param
    dmat. The result is returned in @param fmat. */
class E2Evaluator {
 public:
  virtual bool transform(const ergo_real *dmat, ergo_real *fmat) = 0;
  virtual ~E2Evaluator() {}
};

/* ################################################################### */
/** a collection of vectors, usually handled at once. */
class VarVectorCollection {
  VarVector *vecs;
  unsigned *ages;
  unsigned currentAge;
  int nVecs;
  int nAllocated;
  bool diskMode;
 public:
  explicit VarVectorCollection(int nSize=0) : vecs(NULL), ages(NULL), currentAge(0),
    nVecs(0), nAllocated(0), diskMode(false) {
    if(nSize) setSize(nSize);
  }
  ~VarVectorCollection();
  void setSize(int sz);
  VarVector& operator[](int i);
  int size() const { return nVecs; }
  bool getDiskMode() const { return diskMode; }
  void setDiskMode(bool x) { diskMode = x; }
  /** Make sure there is space for one vector. */
  void release(); 
  /** Release all vectors from the memory, saving if necessary. */
  void releaseAll();
  static const char *tmpdir;
};
 
/* ################################################################### */
/** Abstract interface to a one electron operator. */
class OneElOperator {
 public:
  virtual void getOper(ergo_real *result) = 0;
  virtual ~OneElOperator() {}
};

/** a class implementing dynamic resized two dimensional arrays. */
class SmallMatrix {
  ergo_real *mat;
  int nsize;
 protected:
  struct RowProxy {
    ergo_real *toprow;
    explicit RowProxy(ergo_real *r) : toprow(r) {}
    ergo_real& operator[](int j) const {
      //printf("  returning row %i -> %p\n", j, toprow + j);
      return toprow[j]; }
  };
 public:
  explicit SmallMatrix(int sz) : mat(new ergo_real[sz*sz]), nsize(sz) {}
  ~SmallMatrix() { delete [] mat; }
  const RowProxy operator[](int i) {
    //printf("Returning column %i -> %p\n", i, mat + i*nsize);
    return RowProxy(mat + i*nsize); }
  void expand(int newSize);
};


/* ################################################################### */
/** Linear Response iterative solver using a variant of the Davidson
    method. */
class LRSolver {
 public:

  LRSolver(int nbast, int nocc,
           const ergo_real *fock_matrix,
           const ergo_real *s);
  virtual ~LRSolver() {/* FIXME: delete esub etc */
    if(cmo)   delete cmo; 
    if(fdiag) delete fdiag;
    delete [] xSub;
  }

  /** Computes the residual vector. The residual vector is created by
      solving the problem in the subspace and then using the solution
      coefficients to form the approximate solution vector. This trial
      vector is then substituted to the equation and the residual is
      defined as the difference between the transformed trial vector
      and the expected solution.  */
  virtual bool getResidual(VarVectorCollection& residualv) = 0;

  /** Computes the initial vector the subspace is to be seeded
      with. Allocates @param vecs and returns the number of
      vectors. */
  virtual int getInitialGuess(VarVectorCollection& vecs) = 0;

  /** returns the preconditioning shift. Proper preconditioning is
      vital for the quick convergence. */
  virtual ergo_real getPreconditionerShift(int i) const = 0;

  /** expands above the default limit */
  virtual void increaseSubspaceLimit(int newSize);

  /** Solves the problem defined by the subclass. */
  bool solve(E2Evaluator& e, bool diskMode = false);
  void computeExactE2Diag(E2Evaluator& e2);
  ergo_real convThreshold;   /**< iterative method convergence threshold */
  int       maxSubspaceSize; /**< current subspace size limit. */
 protected:
  static const int MVEC = 200; /**< default limit for subspace size */
  VarVector e2diag; /**< approximation to the diagonal of E2 operator */
  int subspaceSize; /**< current subspace size */
  SmallMatrix eSub; /**< E[2] matrix projected onto subspace */
  SmallMatrix sSub; /**< S[2] matrix projected onto subspace */
  ergo_real *xSub;  /**< solution vector projected onto subspace */

  /** Computes a vector built of base vectors with specified vectors.
      r := Av - f*Sv */
  void getAvMinusFreqSv(ergo_real f, ergo_real *weights,
                        VarVector& r);

  /** Projects vector @param full on the reduced subspace, returns the
      result in @param w which is a preallocated vector of projection
      coefficients (weights) with size equal at least to the subspace
      size. */
  void projectOnSubspace(const VarVector& full, ergo_real *w)/* const*/;
  /** Build full fector from the reduced form */
  void buildVector(const ergo_real *w, VarVector& full) /* const */;

  /** Transform square operator to the vector form */
  void operToVec(OneElOperator& oper, VarVector& res) const;

  /** setE2diag is called by the constructor to fill in the
   approximation of the E[2] operator diagonal.  It returns
   E_LUMO-E_HOMO which is useful for other things. */
  ergo_real setE2diag(int nbast, int nocc,
                      const ergo_real *fock_matrix,
                      const ergo_real *s);
  
  int nbast; /**< number of basis functions */
  int nocc;  /**< number of occupied orbitals */
  VarVectorCollection vects;  /**< base vectors */
  virtual void addToSpace(VarVectorCollection& vecs, E2Evaluator& e2);
  void mo2ao(int nbast, const ergo_real *mo, ergo_real *ao) const;
  void ao2mo(int nbast, const ergo_real *ao, ergo_real *mo) const;
 private:
  /** vects and Avects members store the trial vectors and their
      transformed versions. Only every second vector is stored, the
      paired vectors are recovered with help of swapXY() function. */
  VarVectorCollection Avects; /**< transformed base vectors */
  ergo_real *fdiag;       /**< the eigenvalues of the Fock
                             matrix. Used by load_F_MO. */
  ergo_real *cmo;         /**< the MO coefficients. */

  void load_F_MO(ergo_real *fmat) const;
  bool lintrans(E2Evaluator& e2, const VarVector& v, VarVector& Av) const;
};

/* ################################################################### */
/** Iterative Set Of Linear Equations solver, extending the generic
    LRSolver. */
class SetOfEqSolver : public LRSolver {
  ergo_real frequency; /**< frequency for which the SOE is to be solved. */
  VarVector rhs; /**< RHS of the SOE */
 public:
  /** Creates the set-of-equations solver. The KS and overlap matrix
      may be deleted immediately after the object creation. */
  SetOfEqSolver(int nbast, int nocc,
                const ergo_real *fock_matrix,
                const ergo_real *s,
		ergo_real freq)
    : LRSolver(nbast, nocc, fock_matrix, s), frequency(freq),
    rhsSub(new ergo_real[MVEC]) {};
  void setRHS(OneElOperator& op);    /**< initializes the rhs field */
  virtual ~SetOfEqSolver() { delete [] rhsSub; }
  virtual ergo_real getPreconditionerShift(int) const { return frequency; }
  virtual int getInitialGuess(VarVectorCollection& vecs);
  virtual bool getResidual(VarVectorCollection& residualv);
  /** expands above the default limit */
  virtual void increaseSubspaceLimit(int newSize);
  ergo_real getPolarisability(OneElOperator& oper) /* const */;
 protected:
  ergo_real *rhsSub; /**< RHS vector projected onto subspace */
  virtual void addToSpace(VarVectorCollection& vecs, E2Evaluator& e2);
  ergo_real multiplyXtimesVec(const VarVector& rhs) /* const */;
  ergo_real xTimesRHS;
};


/* ################################################################### */
/** Iterative Eigenvalue solver, extending the generic LRSolver. */
class EigenSolver : public LRSolver {
  ergo_real *ritzVals; /**< recent ritz values in the subspace. */
  ergo_real *transMoms2; /**< most recent SQUARED transition moments. */
  int nStates;    /**< number of excited states to compute */
  int nConverged; /**< number of already converged eigenvalues */
  ergo_real *last_ev; /**< most recent eigenvectors in the reduced space */
 public:
  EigenSolver(int nbast, int nocc,
              const ergo_real *fock_matrix,
              const ergo_real *overlap, int n)
    : LRSolver(nbast, nocc, NULL, NULL), 
    ritzVals(new ergo_real[MVEC]), transMoms2(new ergo_real[MVEC]),
    nStates(n), nConverged(0), last_ev(NULL) {
    ritzVals[0] = setE2diag(nbast, nocc, fock_matrix, overlap);
    for(int i=1; i<n; i++) ritzVals[i] = ritzVals[0];
  }
  virtual ~EigenSolver() {
    if(last_ev) delete last_ev;
    delete [] ritzVals;
    delete [] transMoms2;
  }
  virtual ergo_real getPreconditionerShift(int i) const { 
    return ritzVals[nConverged+i];
  }
  virtual int getInitialGuess(VarVectorCollection& vecs);
  virtual bool getResidual(VarVectorCollection& residualv);
  /** expands above the default limit */
  virtual void increaseSubspaceLimit(int newSize);

  ergo_real getFreq(int i) const { return ritzVals[i]; }
  void computeMoments(OneElOperator& dipx,
                      OneElOperator& dipy,
                      OneElOperator& dipz);
  ergo_real getTransitionMoment2(int i) const { return transMoms2[i]; }
};

} /* End of namespace LR */
#endif /* !defined(SLR_HEADER) */
