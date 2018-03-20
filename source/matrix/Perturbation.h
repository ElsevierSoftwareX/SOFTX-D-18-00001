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

/** @file Perturbation.h Perturbation theory class
 *
 * Copyright(c) Emanuel Rubensson 2008
 *
 * @author Emanuel Rubensson
 * @date June 2008
 *
 */
#ifndef MAT_PERTURBATION
#define MAT_PERTURBATION
namespace per {
  template<typename Treal, typename Tmatrix, typename Tvector>
    class Perturbation {
  public: 
    Perturbation(std::vector<Tmatrix *> const & F, 
		 /**< Vector with matrices (input). */
		 std::vector<Tmatrix *> & D, 
		 /**< Vector with matrices (output). */
		 mat::Interval<Treal> const & gap, /**< Band gap. */
		 mat::Interval<Treal> const & allEigs, 
		 /**< Interval containing all eigenvalues of 
		  *   X0 + delta*X1 + delta^2*X2 + ... 
		  *   for all delta in [0, deltaMax]
		  *   for initial X.
		  */
		 Treal const deltaMax, /**< Largest allowed delta. */
		 Treal const errorTol, /**< Error tolerance. */
		 mat::normType const norm, /**< Norm for truncation etc. */
		 Tvector & vect /**< Vector. */
		 );
    void perturb() {
      dryRun();
      run();
    }

    void checkIdempotencies(std::vector<Treal> & idemErrors);
    template<typename TmatNoSymm>
      void checkCommutators(std::vector<Treal> & commErrors, 
			    TmatNoSymm const & dummyMat);
    void checkMaxSubspaceError(Treal & subsError);
    
  protected:
    /* This is input from the beginning */
    std::vector<Tmatrix *> const & F;
    std::vector<Tmatrix *> & X;
    mat::Interval<Treal> gap;
    mat::Interval<Treal> const & allEigs;
    Treal deltaMax;
    Treal errorTol;
    mat::normType const norm;
    Tvector & vect;

    /* These variables are set in the dry run. */
    int nIter;
    std::vector<Treal> threshVal;
    std::vector<Treal> sigma;
    
    /** Dry run to obtain some needed numbers.
     *
     *  After call to this function we know:
     *   - number of iterations (nIter),
     *   - threshold values (threshVal), and 
     *   - polyunomials to choose (sigma = -1 | = 1)
     *
     *  If requested accuracy is too high or gap too small, an
     *  exception is thrown.
     */
    void dryRun();
    void run();
  private:
    
  };
  
  template<typename Treal, typename Tmatrix, typename Tvector>
    Perturbation<Treal, Tmatrix, Tvector>::
    Perturbation(std::vector<Tmatrix *> const & F_in, 
		 std::vector<Tmatrix *> & X_in, 
		 mat::Interval<Treal> const & gap_in,
		 mat::Interval<Treal> const & allEigs_in,
		 Treal const deltaMax_in,
		 Treal const errorTol_in,
		 mat::normType const norm_in,
		 Tvector & vect_in) 
    : F(F_in), X(X_in), gap(gap_in), allEigs(allEigs_in),
    deltaMax(deltaMax_in), errorTol(errorTol_in), norm(norm_in),
    vect(vect_in) {
      if (!X.empty())
	throw "Perturbation constructor: D vector is expected to be empty (size==0)";
      for (unsigned int iMat = 0; iMat < F.size(); ++iMat)
	X.push_back(new Tmatrix(*F[iMat]));
      
      Treal lmin = allEigs.low();
      Treal lmax = allEigs.upp();
      
      /***** Initial linear transformation of matrix sequence. */
      typename std::vector<Tmatrix *>::iterator matIt = X.begin();
      /* Scale to [0, 1] interval and negate */
      (*matIt)->add_identity(-lmax); 
      *(*matIt) *= ((Treal)1.0 / (lmin - lmax));
      matIt++;
      /* ...and derivatives: */
      for ( ; matIt != X.end(); matIt++ ) 
	*(*matIt) *= ((Treal)-1.0 / (lmin - lmax));
      /* Compute transformed gap */
      gap = (gap - lmax) / (lmin - lmax);
    }
  
  template<typename Treal, typename Tmatrix, typename Tvector>
    void Perturbation<Treal, Tmatrix, Tvector>::dryRun() {
    Treal errorTolPerIter;
    int nIterGuess = 0;
    nIter = 1;
    Treal lumo;
    Treal homo;
    Treal m;
    Treal g;
    while (nIterGuess < nIter) {
      nIterGuess++;
      errorTolPerIter = 0.5 * errorTol /nIterGuess;
      nIter = 0;
      mat::Interval<Treal> gapTmp(gap);
      sigma.resize(0);
      threshVal.resize(0);
      while (gapTmp.low() > 0.5 * errorTol || gapTmp.upp() < 0.5 * errorTol) {
	lumo = gapTmp.low();
	homo = gapTmp.upp();
	m = gapTmp.midPoint();
	g = gapTmp.length();
	if (m > 0.5) {
	  lumo = lumo*lumo;
	  homo = homo*homo;
	  sigma.push_back(-1);
	}
	else {
	  lumo = 2*lumo - lumo*lumo;
	  homo = 2*homo - homo*homo;
	  sigma.push_back(1);
	}
	/* Compute threshold value necessary to converge. */
	Treal forceConvThresh = template_blas_fabs(1-2*m) * g / 10;
	/* We divide by 10 > 2 so that this loop converges at some point. */
	/* Compute threshold value necessary to maintain accuracy in subspace.*/
	Treal subspaceThresh = errorTolPerIter * (homo-lumo) / (1+errorTolPerIter);
	/* Choose the most restrictive threshold of the two. */
	threshVal.push_back(forceConvThresh < subspaceThresh ?
			    forceConvThresh : subspaceThresh);
	homo -= threshVal.back();
	lumo += threshVal.back();
	gapTmp = mat::Interval<Treal>(lumo, homo);
	if (gapTmp.empty())
	  throw "Perturbation<Treal, Tmatrix, Tvector>::dryRun() : Perturbation iterations will fail to converge; Gap is too small or desired accuracy too high.";
	nIter++;
      } /* end 2nd while */
    } /* end 1st while */
    /* Now, we have nIter, threshVal, and sigma. */ 
  }
  
  template<typename Treal, typename Tmatrix, typename Tvector>
    void Perturbation<Treal, Tmatrix, Tvector>::run() {
    Treal const ONE = 1.0;
    mat::SizesAndBlocks rowsCopy;
    X.front()->getRows(rowsCopy);
    mat::SizesAndBlocks colsCopy;
    X.front()->getCols(colsCopy);
    Tmatrix tmpMat;
    //    tmpMat.resetSizesAndBlocks(rowsCopy, colsCopy);
    int nMatrices;
    Treal threshValPerOrder;
    Treal chosenThresh;
    for (int iter = 0; iter < nIter; iter++) {
      std::cout<<"\n\nInside outer loop iter = "<<iter
	       <<"  nIter = "<<nIter
	       <<"  sigma = "<<sigma[iter]<<std::endl;
      /* Number of matrices increases by 1 in each iteration: */
      X.push_back(new Tmatrix);
      nMatrices = X.size();
      X[nMatrices-1]->resetSizesAndBlocks(rowsCopy, colsCopy);
      /* Compute threshold value for each order. */
      threshValPerOrder = threshVal[iter] / nMatrices;
      /* Loop through all nonzero orders. */
      std::cout<<"Entering inner loop nMatrices = "<<nMatrices<<std::endl;
      for (int j = nMatrices-1 ; j >= 0 ; --j) {
      std::cout<<"Inside inner loop j = "<<j<<std::endl;
      std::cout<<"X[j]->eucl() (before compute) = "<<X[j]->eucl(vect,1e-7)<<std::endl;
      std::cout<<"X[j]->frob() (before compute) = "<<X[j]->frob()<<std::endl;
      tmpMat = Treal(Treal(1.0)+sigma[iter]) * (*X[j]);
	std::cout<<"tmpMat.eucl() (before for) = "<<tmpMat.eucl(vect,1e-7)<<std::endl;
      std::cout<<"tmpMat.frob() (before for) = "<<tmpMat.frob()<<std::endl;
	for (int k = 0; k <= j; k++) {
	  /* X[j] = X[j] - sigma * X[k] * X[j-k]      */
	  Tmatrix::ssmmUpperTriangleOnly(-sigma[iter], *X[k], *X[j-k],
					 ONE, tmpMat);
	} /* End 3rd for */
	std::cout<<"tmpMat.eucl() (after for) = "<<tmpMat.eucl(vect,1e-7)<<std::endl;
	*X[j] = tmpMat;
	
	/* Truncate tmpMat, remove if it becomes zero. */
	chosenThresh = threshValPerOrder / pow(deltaMax, Treal(j));
	std::cout<<"X[j]->eucl() (before thresh) = "<<X[j]->eucl(vect,1e-7)<<std::endl;
	std::cout<<"Chosen thresh: "<<chosenThresh<<std::endl;
	Treal actualThresh = X[j]->thresh(chosenThresh, vect, norm);
	std::cout<<"X[j]->eucl() (after thresh) = "<<X[j]->eucl(vect,1e-7)<<std::endl;
#if 1
	/*  If the current matrix is zero AND 
	 *  it is the last matrix
	 */
	if (*X[j] == 0 && int(X.size()) == j+1) {
	  std::cout<<"DELETION: j = "<<j<<"  X.size() = "<<X.size()
		   <<"  X[j] = "<<X[j]<< "  X[j]->frob() = "<<X[j]->frob()
		   <<std::endl;
	  delete X[j];
	  X.pop_back();
	}
	else 
	  std::cout<<"NO DELETION: j = "<<j<<"  X.size() = "<<X.size()
		   <<"  X[j] = "<<X[j]<< "  X[j]->frob() = "<<X[j]->frob()
		   <<std::endl;
#endif
	
      } /* End 2nd for (Loop through orders)     */
    }   /* End 1st for (Loop through iterations) */
  }  /* End run() */
  

  template<typename Treal, typename Tmatrix, typename Tvector>
    void Perturbation<Treal, Tmatrix, Tvector>::
    checkIdempotencies(std::vector<Treal> & idemErrors) {
    Tmatrix tmpMat;
    Treal const ONE = 1.0;
    unsigned int j;
    for (unsigned int m = 0; m < X.size(); ++m) {
      tmpMat = (-ONE) * (*X[m]);
      for (unsigned int i = 0; i <= m; ++i) {
	j = m - i;
	/* TMP = TMP + X[i] * X[j]      */
	Tmatrix::ssmmUpperTriangleOnly(ONE, *X[i], *X[j], ONE, tmpMat);
      }
      /* TMP is symmetric! */
      idemErrors.push_back(tmpMat.eucl(vect,1e-10));
    }
  }

  template<typename Treal, typename Tmatrix, typename Tvector>
    template<typename TmatNoSymm>
    void Perturbation<Treal, Tmatrix, Tvector>::
    checkCommutators(std::vector<Treal> & commErrors, 
		     TmatNoSymm const & dummyMat) {
    mat::SizesAndBlocks rowsCopy;
    X.front()->getRows(rowsCopy);
    mat::SizesAndBlocks colsCopy;
    X.front()->getCols(colsCopy);
    TmatNoSymm tmpMat;
    tmpMat.resetSizesAndBlocks(rowsCopy, colsCopy);
    Treal const ONE = 1.0;
    unsigned int j;
    for (unsigned int m = 0; m < X.size(); ++m) {
      tmpMat = 0;
      std::cout<<"New loop\n";
      for (unsigned int i = 0; i <= m && i < F.size(); ++i) {
	j = m - i;
	std::cout<<i<<", "<<j<<std::endl;
	/* TMP = TMP + F[i] * X[j] - X[j] * F[i]    */
	tmpMat += ONE * (*F[i]) * (*X[j]);
	tmpMat += -ONE * (*X[j]) * (*F[i]);
      }
      /* TMP is not symmetric! */
      commErrors.push_back(tmpMat.frob());
    }   
  }
  
  
  template<typename Treal, typename Tmatrix, typename Tvector>
    void Perturbation<Treal, Tmatrix, Tvector>::
    checkMaxSubspaceError(Treal & subsError) {
    Treal const ONE = 1.0;
    Tmatrix XdeltaMax(*F.front());
    for (unsigned int ind = 1; ind < F.size(); ++ind)
      XdeltaMax += pow(deltaMax, Treal(ind)) * (*F[ind]);
    /***** Initial linear transformation of matrix. */
    Treal lmin = allEigs.low();
    Treal lmax = allEigs.upp();
    /* Scale to [0, 1] interval and negate */
    XdeltaMax.add_identity(-lmax); 
    XdeltaMax *= ((Treal)1.0 / (lmin - lmax));
    
    Tmatrix X2;
    for (int iter = 0; iter < nIter; iter++) {
      X2 = ONE * XdeltaMax * XdeltaMax;
      if (sigma[iter] == Treal(1.0)) {
	XdeltaMax *= 2.0; 
	XdeltaMax -= X2;
      }
      else {
	XdeltaMax = X2; 
      }
    }   /* End of for (Loop through iterations) */
    
    Tmatrix DdeltaMax(*X.front());
    for (unsigned int ind = 1; ind < X.size(); ++ind)
      DdeltaMax += pow(deltaMax, Treal(ind)) * (*X[ind]);
    subsError = Tmatrix::eucl_diff(XdeltaMax,DdeltaMax,
				   vect, errorTol *1e-2); 
  }
  


} /* end namespace mat */
#endif

