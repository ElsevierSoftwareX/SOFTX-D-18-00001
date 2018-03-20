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


/** @file GetDensFromFock.cc
 *
 *  @brief Routines for getting density matrix from a given Fock
 *         matrix.
 *
 *  @author Anastasia Kruchinina <em>responsible</em>
 *  @author Elias Rudberg
 */



#include "output.h"
#include <memory.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include "utilities.h"
#include "matrix_utilities.h"
#include "TC2.h"
#include "units.h"
#include "machine_epsilon.h"
#include "AllocatorManager.h"

#include "densfromf_full.h"

#include "purification_general.h"
#include "purification_sp2.h"
#include "purification_sp2acc.h"
#include "GetDensFromFock.h"

typedef generalVector VectorType;

const int       GetDensFromFock::UNDEF_VALUE        = -1;
const int       GetDensFromFock::UNDEF_VALUE_UINT   = -1;
const ergo_real GetDensFromFock::UNDEF_VALUE_REAL   = -1;
const string    GetDensFromFock::UNDEF_VALUE_STRING = "";
const int       GetDensFromFock::SET   = 1;
const int       GetDensFromFock::UNSET = 0;


/** Choose which method to use for computing the density matrix from Fock matrix.
 *
 * Possible alternatives:
 *  - use recursive expansion
 *  - use diagonalization
 */
int GetDensFromFock::get_dens_from_fock(symmMatrix&   Finput,
                                        symmMatrix&   resultDens,
                                        symmMatrix&   F_ort_prev,
                                        generalVector *eigVecLUMO,
                                        generalVector *eigVecHOMO)
{
   Util::TimeMeter timeMeterTot;

   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "get_dens_from_fock_general, n = %i, use_diagonalization = %i, use_diag_on_error = %i",
             n, use_diagonalization, use_diag_on_error);
   resultEntropyTerm = 0; // In nonzero temperature case, this will be set to nonzero value later.

   std::string allocStatsStr1 = mat::AllocatorManager<ergo_real>::instance().getStatistics();
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Before writeAndReadAll(): %s", allocStatsStr1.c_str());

   Util::TimeMeter timeMeterWriteAndReadAll;
   std::string     sizesStr = mat::FileWritable::writeAndReadAll();
   timeMeterWriteAndReadAll.print(LOG_AREA_DENSFROMF, "FileWritable::writeAndReadAll");
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, ((std::string)"writeAndReadAll sizesStr: '" + sizesStr).c_str());

   std::string allocStatsStr2 = mat::AllocatorManager<ergo_real>::instance().getStatistics();
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "After writeAndReadAll(): %s", allocStatsStr2.c_str());

   if (noOfOccupiedOrbs == 0)
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "GetDensFromFock::get_dens_from_fock: (noOfOccupiedOrbs == 0), skipping.");
      resultDens.readFromFile();
      resultDens.clear();
      resultDens.writeToFile();
      return 0;
   }



   int use_diag = 0;
   int purification_has_failed = 0;

   if (use_diagonalization == UNDEF_VALUE)
   {
      throw "Error in get_dens_from_fock (GetDensFromFock class) : use_diagonalization flag is not specified";
   }

   if (use_diagonalization == SET)
   {
      use_diag = 1;
   }
   else
   {
      // Try purification

      if (electronicTemperature != 0)
      {
         throw "Error: (electronicTemperature != 0) not implemented for sparse case.";
      }
      resultDens.readFromFile();
      resultDens.clear();

      int puri_res;

      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,
                "calling get_dens_from_fock_sparse, n = %6i, subspaceErrorLimit = %g",
                n, (double)subspaceErrorLimit);
      puri_res = get_dens_from_fock_sparse(Finput,
                                           resultDens,
                                           F_ort_prev,
                                           eigVecLUMO,
                                           eigVecHOMO);


      if (puri_res != 0)
      {
         // Something was wrong...
         if (use_diag_on_error == UNDEF_VALUE)
         {
            throw "Error in get_dens_from_fock (GetDensFromFock class) : use_diag_on_error flag is not specified";
         }

         if (use_diag_on_error == SET)
         {
            do_output(LOG_CAT_ERROR, LOG_AREA_DENSFROMF, "get_dens_from_fock  (GetDensFromFock class) : error while getting density matrix; trying with diagonalization instead.");
            use_diag = 1;
            purification_has_failed = 1;
         }
         else
         {
            do_output(LOG_CAT_ERROR, LOG_AREA_DENSFROMF, "get_dens_from_fock  (GetDensFromFock class) : error while getting density matrix; aborting.");
            return -1;
         }
      }
      else
      {
         // Purification success!
         do_output(LOG_CAT_ERROR, LOG_AREA_DENSFROMF, "get_dens_from_fock  (GetDensFromFock class) : purification finished OK.");
      }
      resultDens.writeToFile();
   }


   if (use_diag == 1)
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "calling get_dens_from_fock_full, n = %i", n);

      std::vector<ergo_real> F_full(n *n);
      std::vector<ergo_real> S_full(n *n);

      {
         // Create full matrix versions of F and S
         normalMatrix *tmpMat;
         Finput.readFromFile();
         tmpMat = new normalMatrix(Finput);
         Finput.writeToFile();
         tmpMat->fullMatrix(F_full);
         delete tmpMat;
         overlapMatrix.readFromFile();
         tmpMat = new normalMatrix(overlapMatrix);
         overlapMatrix.writeToFile();
         tmpMat->fullMatrix(S_full);
         delete tmpMat;
      }

      std::vector<ergo_real> densityMatrixFull(n *n);
      std::vector<ergo_real> eigVecLUMO_tmp(n);
      std::vector<ergo_real> eigVecHOMO_tmp(n);


      if (store_all_eigenvalues_to_file == UNDEF_VALUE)
      {
         throw "Error in get_dens_from_fock (GetDensFromFock class) : store_all_eigenvalues_to_file flag is not specified";
      }

      ergo_real gap = -1;
      if (get_dens_from_fock_full(n,
                                  noOfOccupiedOrbs,
                                  &densityMatrixFull[0],
                                  &F_full[0],
                                  &S_full[0],
                                  factor,
                                  electronicTemperature,
                                  resultEntropyTerm,
                                  gap,
                                  store_all_eigenvalues_to_file,
                                  &eigVecLUMO_tmp[0],
                                  &eigVecHOMO_tmp[0]) != 0)
      {
         throw "error in get_dens_from_fock_full";
      }

      if (purification_has_failed)
      {
         // Accept purification failure only in case gap is very small.
         ergo_real gapLimit = 1e-4;
         if (gap > gapLimit)
         {
            throw "Error in GetDensFromFock::get_dens_from_fock: purification failed and (gap > gapLimit). Purification should not fail in such cases; something is wrong.";
         }
      }

      resultDens.readFromFile();
      resultDens.assignFromFull(densityMatrixFull);
      resultDens.writeToFile();
      if (eigVecLUMO)
      {
         eigVecLUMO->assign_from_full(eigVecLUMO_tmp, matrixSizesAndBlocks);
      }
      if (eigVecHOMO)
      {
         eigVecHOMO->assign_from_full(eigVecHOMO_tmp, matrixSizesAndBlocks);
      }
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "get_dens_from_fock_full finished");
   }

   timeMeterTot.print(LOG_AREA_DENSFROMF, "get_dens_from_fock");

   return 0;
}


#ifndef USE_CHUNKS_AND_TASKS
static ergo_real get_eucl_diff_with_adapted_accuracy(int n,
						     const symmMatrixWrap & F_w,
						     const symmMatrixWrap & F_ort_prev_w,
						     ergo_real acc) {
  // The symmMatrixWrap::eucl_diff() call may be slow, we use a maxIter param to detect if it is a difficult case, and in such cases use a larger acc value.
  int maxIterForEuclDiff = std::max(n / 10, 500);
  ergo_real maxEigValMovement_eucl = -1; // Value will be set in try/catch code below.
  try {
    Util::TimeMeter timeMeterEuclDiff;
    maxEigValMovement_eucl = symmMatrixWrap::eucl_diff(F_w, F_ort_prev_w, acc, maxIterForEuclDiff) + acc;
    timeMeterEuclDiff.print(LOG_AREA_DENSFROMF, "symmMatrixWrap::eucl_diff for maxEigValMovement_eucl ");
  }
  catch(...) {
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "symmMatrixWrap::eucl_diff() for maxEigValMovement_eucl failed for maxIterForEuclDiff=%d. Calling eucl_diff() again with lower accuracy requirement sqrt(acc).",
	      maxIterForEuclDiff);
    ergo_real acc2 = template_blas_sqrt(acc);
    Util::TimeMeter timeMeterEuclDiff;
    maxEigValMovement_eucl = symmMatrixWrap::eucl_diff(F_w, F_ort_prev_w, acc2) + acc2;
    timeMeterEuclDiff.print(LOG_AREA_DENSFROMF, "symmMatrixWrap::eucl_diff for maxEigValMovement_eucl ");
  }
  return maxEigValMovement_eucl;
}
#endif


/** Use recursive expansion for computing the density matrix from Fock matrix.
 *
 * Construct approximation of the step function by recursive
 * application of low order polynomials. Sparsity is preserved using
 * truncation (see J. Chem. Phys. 128, 074106, 2008), which can be
 * done using spectral, Frobenius or mixed norms (see
 * J. Comput. Chem. 30.6 (2009): 974-977.).
 *
 * Possible alternatives (use_acceleration parameter):
 * - SP2 recursive expansion
 * - SP2 accelerated recursive expansion
 */
int GetDensFromFock::get_dens_from_fock_sparse(symmMatrix&   F,
                                               symmMatrix&   resultDens,
                                               symmMatrix&   F_ort_prev,
                                               generalVector *eigVecLUMO,
                                               generalVector *eigVecHOMO)
{
#ifdef USE_CHUNKS_AND_TASKS
   if (output_homo_and_lumo_eigenvectors == SET)
   {
      throw "Error in get_dens_from_fock_sparse (GetDensFromFock class) : computation of eigenvectors is not implemented with Chunks and Tasks.";
   }
   if ((leavesSizeMax == UNDEF_VALUE_UINT) || (blocksize == UNDEF_VALUE_UINT))
   {
      throw "Error in get_dens_from_fock (GetDensFromFock class) : leavesSizeMax and/or blocksize flags are not specified";
   }
     size_t NRows = n; // defined in GetDensFromFock.h
     size_t NCols = n;
     #ifdef USE_CHUNKS_AND_TASKS_BSM
           ParamsType params(leavesSizeMax, blocksize, NRows, NCols);
     #else
           ParamsType params(leavesSizeMax, NRows, NCols);
     #endif
#else
         ParamsType params;
#endif

   if (use_new_stopping_criterion == UNDEF_VALUE)
   {
      throw "Error in get_dens_from_fock_sparse (GetDensFromFock class) : use_new_stopping_criterion flag is not specified";
   }

   Util::TimeMeter timeMeterTot;
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "get_dens_from_fock_sparse() start!");

   int stopping_criterion;
   if (use_new_stopping_criterion == SET)
   {
      stopping_criterion = 1;
      // If we use the new stopping criterion, then we should stop
      // when the idempotency error cannot be decreased
      // significantly anymore. The error in the subspace
      // accumulates during the iterations, thus we can expect that
      // it would be larger than the error in eigenvalues.
      eigvalueErrorLimit = subspaceErrorLimit;
   }
   else
   {
      stopping_criterion = 0;
   }

#ifdef USE_CHUNKS_AND_TASKS
   if ((truncationNormPurification == mat::euclNorm) || (truncationNormPurification == mat::mixedNorm))
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "get_dens_from_fock_sparse():  set truncation "
                                                  "norm to mat::frobNorm since spectral and mixed norms are not implemented");
      set_truncationNormPurification(mat::frobNorm);
   }
   if ((stopCriterionNormPurification == mat::euclNorm) || (stopCriterionNormPurification == mat::mixedNorm))
   {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "get_dens_from_fock_sparse():  set stopping criterion "
                                                  "norm to mat::frobNorm since spectral and mixed norms are is not implemented");
      set_stopCriterionNormPurification(mat::frobNorm);
   }
#endif

   // Select tolerated errors in the occupied subspace for the three truncations
   // and for purification.
   ergo_real subspaceThr_1    = 0.1 * subspaceErrorLimit;
   ergo_real subspaceThr_Puri = 0.7 * subspaceErrorLimit;
   ergo_real subspaceThr_2    = 0.1 * subspaceErrorLimit;
   ergo_real subspaceThr_3    = 0.1 * subspaceErrorLimit;

   // Select tolerated errors in eigenvalues
   ergo_real eigvalueThr_Puri = 0.7 * eigvalueErrorLimit;
   ergo_real eigvalueThr_2    = 0.15 * eigvalueErrorLimit;
   ergo_real eigvalueThr_3    = 0.15 * eigvalueErrorLimit;

   symmMatrix F_tmp(F);
   F_tmp.readFromFile();  // F stay written in the file
   std::string allocStatsStr3 = mat::AllocatorManager<ergo_real>::instance().getStatistics();
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "After F.readFromFile(): %s", allocStatsStr3.c_str());


   symmMatrixWrap F_w;
   transform_matrix_from_to(F_tmp, F_w, params);



   std::string allocStatsStr4 = mat::AllocatorManager<ergo_real>::instance().getStatistics();
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "After creating wrapper MatrixType F_w(F): %s", allocStatsStr4.c_str());
   F_tmp.clear();  //  we do not need this matrix anymore, it is saved in the wrapper

   // transform to the orthogonal basis
   {
      triangMatrix invCholFactor_tmp(invCholFactor);
      invCholFactor_tmp.readFromFile();
      output_current_memory_usage(LOG_AREA_DENSFROMF, "In get_dens_from_fock_sparse, before F_w to orthogonal basis");

      Util::TimeMeter timeMeterFortTransf;

      triangMatrixWrap invCholFactor_tmp_w;
      transform_matrix_from_to(invCholFactor_tmp, invCholFactor_tmp_w, params);
      invCholFactor_tmp.clear();
    
      F_w = transpose(invCholFactor_tmp_w) * F_w * invCholFactor_tmp_w;

      timeMeterFortTransf.print(LOG_AREA_DENSFROMF, " F_w to orthogonal basis");
      output_current_memory_usage(LOG_AREA_DENSFROMF, "In get_dens_from_fock_sparse,  after  "
                                                      "F_ort = tr(Z) * F_S * Z");
   }

   // Now F_w contains F_ort.

   //Compare F to F_ort_prev to check how far eigenvalues move.
   F_ort_prev.readFromFile();
   output_current_memory_usage(LOG_AREA_DENSFROMF,
                               "In get_dens_from_fock_sparse,  after F_ort_prev.readFromFile()");
                              

  symmMatrixWrap F_ort_prev_w;
  transform_matrix_from_to(F_ort_prev, F_ort_prev_w, params);

   output_current_memory_usage(LOG_AREA_DENSFROMF,
                               "In get_dens_from_fock_sparse,  after symmMatrixWrap F_ort_prev_w(F_ort_prev)");
   F_ort_prev.clear();

   // intervals will be expanded to make sure that they contain the HOMO and LUMO eigenvalues of F_ort
   intervalType lumoInterval_F_ort_prev_expanded;
   intervalType homoInterval_F_ort_prev_expanded;

   {
      ergo_real maxEigValMovement_frob = symmMatrixWrap::frob_diff(F_w, F_ort_prev_w);
      output_current_memory_usage(LOG_AREA_DENSFROMF,
                                  "In get_dens_from_fock_sparse,  after getting maxEigValMovement_frob ");

#ifndef USE_CHUNKS_AND_TASKS
      ergo_real       acc = template_blas_sqrt(get_machine_epsilon());
      Util::TimeMeter timeMeterMixedDiff;
      ergo_real       maxEigValMovement_mixed = symmMatrixWrap::mixed_diff(F_w, F_ort_prev_w, acc) + acc;
      timeMeterMixedDiff.print(LOG_AREA_DENSFROMF, "MatrixType::mixed_diff for maxEigValMovement_mixed");
      output_current_memory_usage(LOG_AREA_DENSFROMF,
                                  "In get_dens_from_fock_sparse,  after getting maxEigValMovement_mixed");
      ergo_real maxEigValMovement_eucl = get_eucl_diff_with_adapted_accuracy(n, F_w, F_ort_prev_w, acc);
      output_current_memory_usage(LOG_AREA_DENSFROMF,
                                  "In get_dens_from_fock_sparse,  after getting maxEigValMovement_eucl ");
#endif

      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "maxEigValMovement_frob  = %22.11f", (double)maxEigValMovement_frob);
#ifndef USE_CHUNKS_AND_TASKS
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "maxEigValMovement_mixed = %22.11f", (double)maxEigValMovement_mixed);
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "maxEigValMovement_eucl  = %22.11f", (double)maxEigValMovement_eucl);
#endif
      // Increase HOMO/LUMO intervals so that they for sure contain the HOMO and LUMO eigenvalues of F_ort

      // Anastasia comment: it may happen for very small cases that
      // bounds for homo and lumo are very good and intervals are
      // actually empty sets, then the increase() function will throw
      // exception
      if (homoInterval_F_ort_prev.low() >= homoInterval_F_ort_prev.upp())
      {
         homoInterval_F_ort_prev = intervalType(homoInterval_F_ort_prev.low(), homoInterval_F_ort_prev.low() + 1e-10);
      }
      if (lumoInterval_F_ort_prev.low() >= lumoInterval_F_ort_prev.upp())
      {
         lumoInterval_F_ort_prev = intervalType(lumoInterval_F_ort_prev.upp() - 1e-10, lumoInterval_F_ort_prev.upp());
      }

#ifdef USE_CHUNKS_AND_TASKS
      ergo_real maxEigValMovement = maxEigValMovement_frob;
#else
      ergo_real maxEigValMovement = maxEigValMovement_eucl;
#endif

      lumoInterval_F_ort_prev_expanded  = lumoInterval_F_ort_prev;
      homoInterval_F_ort_prev_expanded  = homoInterval_F_ort_prev;

      homoInterval_F_ort_prev_expanded.increase(maxEigValMovement);
      lumoInterval_F_ort_prev_expanded.increase(maxEigValMovement);


      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "lumo before truncation: [ %.12lf , %.12lf ]", (double)lumoInterval_F_ort_prev_expanded.low(), (double)lumoInterval_F_ort_prev_expanded.upp());
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "homo before truncation: [ %.12lf , %.12lf ]", (double)homoInterval_F_ort_prev_expanded.low(), (double)homoInterval_F_ort_prev_expanded.upp());
   }



   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Truncate matrix F and update intervals for homo and lumo.");

   /* EMANUEL COMMENT:
    * - truncationNorm can be set to any of
    * mat::frobNorm, mat::mixedNorm, or mat::euclNorm
    * The best choice depends on a trade-off between spending
    * time in truncation and in matrix-matrix multiplication.
    */
   mat::normType truncationNorm = truncationNormPurification;


   // Now, we will truncate F (and update eigenvalue intervals):
   ergo_real truncError_1;
   {
      ergo_real gapMin = lumoInterval_F_ort_prev_expanded.low() - homoInterval_F_ort_prev_expanded.upp();
      ergo_real gapMax = lumoInterval_F_ort_prev_expanded.upp() - homoInterval_F_ort_prev_expanded.low();
      ergo_real threshold_1;
      // We consider the gap to be accurately known if the uncertainty is at most 10 %
      if ((gapMin > 0) && ((gapMax - gapMin) / gapMin < 0.1))
      {
         // Gap is accurately known: we use gapMin
         threshold_1 = subspaceThr_1 * gapMin / (1 + subspaceThr_1);
      }
      else
      {
         // Gap is not accurately known. To avoid choosing a very tight
         // threshold value due to a small lower bound for the gap, we
         // use the largest of 'gap_expected_lower_bound' and calculated
         // 'gapMin':
         threshold_1 = gapMin > gap_expected_lower_bound ?
                       subspaceThr_1 * gapMin / (1 + subspaceThr_1) :
                       subspaceThr_1 * gap_expected_lower_bound / (1 + subspaceThr_1);
      }

      double nnzF_before_trunc_pc = (double)F_w.nnz() * 100 / ((double)n * n);
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Truncating F_ort ( %s ), selected threshold = %10.6g",
                mat::getNormTypeString(truncationNorm).c_str(), (double)threshold_1);
      Util::TimeMeter timeMeterFThresh;
#ifdef USE_CHUNKS_AND_TASKS
      truncError_1 = F_w.thresh_frob(threshold_1);
#else
      truncError_1 = F_w.thresh(threshold_1, truncationNorm);
#endif
      double nnzF_after_trunc_pc = (double)F_w.nnz() * 100 / ((double)n * n);
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,
                "Truncated F_ort ( %s ), selected threshold = %10.6g, returned error = %10.6g, nnz before = %3.4f %%, nnz after = %3.4f %%",
                mat::getNormTypeString(truncationNorm).c_str(), (double)threshold_1, (double)truncError_1, nnzF_before_trunc_pc, nnzF_after_trunc_pc);
      timeMeterFThresh.print(LOG_AREA_DENSFROMF, "Truncation of F_ort");
      puri_stats[stats_prefix + "nnz_percentage_F_ort"] = nnzF_after_trunc_pc;

      // Increase HOMO and LUMO intervals so that they contain the eigenvalues of the truncated matrix:
      homoInterval_F_ort_prev_expanded.increase(truncError_1);
      lumoInterval_F_ort_prev_expanded.increase(truncError_1);
   }


   F_ort_prev_w.clear();

   transform_matrix_from_to(F_w, F_ort_prev, params);
   
   F_ort_prev.writeToFile();



   // The HOMO and LUMO intervals now contain the HOMO and LUMO
   // eigenvalues of F_ort_prev but improved values will hopefully be
   // calculated in purification.

   intervalType homoIntervalSaved = homoInterval_F_ort_prev_expanded;
   intervalType lumoIntervalSaved = lumoInterval_F_ort_prev_expanded;


   if (use_acceleration == UNDEF_VALUE)
   {
      throw "Error in get_dens_from_fock_sparse (GetDensFromFock class) : use_acceleration flag is not specified";
   }

#ifdef USE_CHUNKS_AND_TASKS
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Chosen CHT Wrapper");
#else
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Chosen Ergo Wrapper");
#endif


   PurificationGeneral<symmMatrixWrap> *Puri;  // abstract class

   if (use_acceleration == SET)
   {
      Puri = new Purification_sp2acc<symmMatrixWrap>();
   }
   else
   {
      Puri = new Purification_sp2<symmMatrixWrap>();
   }

   mat::Gblas::timekeeping = true;
   mat::Gblas::time        = 0;

   Puri->initialize(F_w,
                    lumoIntervalSaved,
                    homoIntervalSaved,
                    maxMul,
                    subspaceThr_Puri,
                    eigvalueThr_Puri,
                    stopping_criterion,  // 1 = new, 0 = old
                    truncationNorm,
                    stopCriterionNormPurification,
                    noOfOccupiedOrbs
                    );


   if (output_homo_and_lumo_eigenvectors == UNDEF_VALUE)
   {
      throw "Error in get_dens_from_fock_sparse (GetDensFromFock class) : output_homo_and_lumo_eigenvectors flag is not specified";
   }


   /*
    * eigVecHOMO and eigVecLUMO are pointers to some empty vectors.
    * If we do not want to compute eigenvectors now, set them to NULL.
    * In other case, just send this pointers to needed function.
    */

   if ((output_homo_and_lumo_eigenvectors == UNSET) || (eigenvectors_method == ""))
   {
      eigVecLUMO  = NULL;
      eigVecHOMO  = NULL;
   }
   else
   {
      // check if parameters are specified
      // all these parameters are input parameters to Ergo
      // they set in SCF_general.cc
      if (try_eigv_on_next_iteration_if_fail == UNDEF_VALUE)
      {
         throw "Error in get_dens_from_fock_sparse (GetDensFromFock class) : try_eigv_on_next_iteration_if_fail flag is not specified";
      }
      if (use_prev_vector_as_initial_guess == UNDEF_VALUE)
      {
         throw "Error in get_dens_from_fock_sparse (GetDensFromFock class) : use_prev_vector_as_initial_guess flag is not specified";
      }

      if (puri_compute_eigv_in_each_iteration == UNDEF_VALUE)
      {
         throw "Error in get_dens_from_fock_sparse (GetDensFromFock class) : puri_compute_eigv_in_each_iteration flag is not specified";
      }
      if (run_shift_and_square_method_on_F == UNDEF_VALUE)
      {
         throw "Error in get_dens_from_fock_sparse (GetDensFromFock class) : run_shift_and_square_method_on_F flag is not specified";
      }

      if (eigensolver_maxiter == UNDEF_VALUE_UINT)
      {
         throw "Error in get_dens_from_fock_sparse (GetDensFromFock class) : eigensolver_maxiter value is not specified";
      }
      if (eigensolver_accuracy == UNDEF_VALUE_REAL)
      {
         throw "Error in get_dens_from_fock_sparse (GetDensFromFock class) : eigensolver_accuracy value is not specified";
      }

      if (save_permuted_F_matrix_in_bin == UNDEF_VALUE)
      {
         throw "Error in get_dens_from_fock_sparse (GetDensFromFock class) : save_permuted_F_matrix_in_bin  flag is not specified";
      }

      Puri->set_eigenvectors_params(eigenvectors_method,
                                    eigenvectors_iterative_method,
                                    eigensolver_accuracy,
                                    eigensolver_maxiter,
                                    SCF_step,
                                    use_prev_vector_as_initial_guess,
                                    try_eigv_on_next_iteration_if_fail,
                                    eigVecLUMO,
                                    eigVecHOMO
                                    );

      if (puri_compute_eigv_in_each_iteration == SET)
      {
         Puri->set_compute_eigenvectors_in_each_iteration();
      }


      if (run_shift_and_square_method_on_F == SET)
      {
         // COMPUTE EIGENVALUES FOR F (FOR COMPARISON)
         int eigsolver_maxiter = 5000;
         Puri->compute_eigenvectors_without_diagonalization_on_F(F_w, eigsolver_maxiter);
      }
   }


   // save Hamiltonian in this permutation if needed to compute eigenvectors
   // for plotting error in matlab
   if ((save_permuted_F_matrix_in_bin == SET) && (SCF_step >= 1))
   {
      symmMatrix Y;
      transform_matrix_from_to(F_w, Y, params);
      vector<int>  Itmp, I, Jtmp, J;
      vector<real> Vtmp, V;
      Y.get_all_values(Itmp, Jtmp, Vtmp);
   
      size_t nnz = 0;
      // Count nonzeros
      for (size_t i = 0; i < Itmp.size(); i++)
      {
         nnz += (Vtmp[i] != 0);
      }
   
      I.reserve(nnz);
      J.reserve(nnz);
      V.reserve(nnz);
      // Extract nonzeros
      for (size_t i = 0; i < Itmp.size(); i++)
      {
         if (Vtmp[i] != 0)
         {
            I.push_back(Itmp[i]);
            J.push_back(Jtmp[i]);
            V.push_back(Vtmp[i]);
         }
      }
      ostringstream name;
      name << "F.bin";
      write_matrix_to_bin(name.str().c_str(), I, J, V, Y.get_nrows());
      name.str("");
   }


   F_w.clear();

   output_current_memory_usage(LOG_AREA_DENSFROMF, "Before  Puri->PurificationStart()");
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,
             "calling Puri->PurificationStart(), number of threads = %i, trunc norm '%s'",
             mat::Params::getNProcs(), mat::getNormTypeString(truncationNorm).c_str());
   mat::FileWritable::resetStats();
   time_t puriStartWallTime;
   time(&puriStartWallTime);

   try
   {
      Util::TimeMeter timeMeterPurification;
      Puri->PurificationStart();
      timeMeterPurification.print(LOG_AREA_DENSFROMF, "Puri->PurificationStart()");
      //Puri->info.print_collected_info();
   }
   catch (...)
   {
      if (Puri != NULL)
      {
         delete Puri;
      }
      throw;
   }

   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Purification finished, %3i multiplications, total time %lf", Puri->info.total_it, Puri->info.total_time);

   {
      std::stringstream ss;
      ss << "Accumulated wall times for writeToFile in PurificationStart()          : " << mat::FileWritable::getStatsTimeWrite();
      do_output(LOG_CAT_TIMINGS, LOG_AREA_DENSFROMF, ss.str().c_str());
   }
   {
      std::stringstream ss;
      ss << "Accumulated wall times for readFromFile in PurificationStart()         : " << mat::FileWritable::getStatsTimeRead();
      do_output(LOG_CAT_TIMINGS, LOG_AREA_DENSFROMF, ss.str().c_str());
   }
   {
      std::stringstream ss;
      ss << "Accumulated wall times for copy and assign in PurificationStart()      : " << mat::FileWritable::getStatsTimeCopyAndAssign();
      do_output(LOG_CAT_TIMINGS, LOG_AREA_DENSFROMF, ss.str().c_str());
   }


   {
      std::stringstream ss;
      ss << "Number of calls to writeToFile in PurificationStart()                  : " << mat::FileWritable::getStatsCountWrite();
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, ss.str().c_str());
   }
   {
      std::stringstream ss;
      ss << "Number of calls to readFromFile in PurificationStart()                 : " << mat::FileWritable::getStatsCountRead();
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, ss.str().c_str());
   }
   {
      std::stringstream ss;
      ss << "Number of calls to FileWritable copy and assign in PurificationStart() : " << mat::FileWritable::getStatsCountCopyAndAssign();
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, ss.str().c_str());
   }

   do_output(LOG_CAT_TIMINGS, LOG_AREA_DENSFROMF,
             "mat::Gblas::time after purification : %12.6f",
             (double)mat::Gblas::time);

   output_current_memory_usage(LOG_AREA_DENSFROMF, "After   Puri->PurificationStart()");

   symmMatrixWrap D_w(Puri->X);
   Puri->clear();    // delete matrices from Puri
   intervalType homoIntervalNew = intervalType(Puri->info.homo_estim_low_F, Puri->info.homo_estim_upp_F);
   intervalType lumoIntervalNew = intervalType(Puri->info.lumo_estim_low_F, Puri->info.lumo_estim_upp_F);


   if (plot_puri_results == SET)
   {
      // plot results
      ostringstream name;
      name << "puri_out_error_" << SCF_step << plot_puri_results_str << ".m";
      Puri->gen_matlab_file_norm_diff(name.str().c_str());
      name.str("");

      name << "puri_out_threshold_" << SCF_step << plot_puri_results_str << ".m";
      Puri->gen_matlab_file_threshold(name.str().c_str());
      name.str("");

      name << "puri_out_nnz_" << SCF_step << plot_puri_results_str << ".m";
      Puri->gen_matlab_file_nnz(name.str().c_str());
      name.str("");

      name << "puri_out_eigs_" << SCF_step << plot_puri_results_str << ".m";
      Puri->gen_matlab_file_eigs(name.str().c_str());
      name.str("");

      name << "puri_out_time_" << SCF_step << plot_puri_results_str << ".m";
      Puri->gen_matlab_file_time(name.str().c_str());
      name.str("");
   }
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Created .m files with results of the purification");

   if (!Puri->info.converged)
   {
      // controlled by the Ergo input flag "purification_ignore_failure"
      if (ignore_purification_failure == UNDEF_VALUE)
      {
         throw "Error in get_dens_from_fock_sparse (GetDensFromFock class) : use_diag_on_error flag is not specified";
      }


      if (ignore_purification_failure == UNSET)
      {
         do_output(LOG_CAT_ERROR, LOG_AREA_DENSFROMF,
                   "Error in purification: Puri->info.converged() "
                   "returned false.");
         Puri->info.print_collected_info();
         return -1;
      }
      else
      {
         do_output(LOG_CAT_WARNING, LOG_AREA_DENSFROMF, "Purification did NOT converged, ignoring.");
      }
   }
   else
   {
      ergo_real acc_error = (ergo_real)Puri->info.accumulated_error_subspace;
      if (acc_error == -1)  // if we do not know that gap
      {
         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,
                   "Purification converged OK");
      }
      else
      {
         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,
                   "Purification converged OK, subspaceError <= %22.11f", acc_error);
      }
   }


   // compute eigenvectors
   if ((eigVecLUMO != NULL) || (eigVecHOMO != NULL))
   {
      triangMatrix invCholFactor_tmp(invCholFactor);
      invCholFactor_tmp.readFromFile();
      /* here: if eigenvector is not computed, it is empty, not NULL */
      if (Puri->info.lumo_eigenvector_is_computed && !eigVecLUMO->is_empty())
      {
         printf("1  %d  %lf  %d  %lf\n",
                Puri->info.lumo_eigenvector_is_computed_in_iter,
                (double)Puri->info.eigValLUMO,
                Puri->info.lumo_eigensolver_iter,
                Puri->info.lumo_eigensolver_time);

         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "LUMO eigenvector is computed.");
         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Perform congruence transformation.");
         // perform congruence transformation
         (*eigVecLUMO) = invCholFactor_tmp * (*eigVecLUMO);
      }


      if (Puri->info.homo_eigenvector_is_computed && !eigVecHOMO->is_empty())
      {
         printf("2  %d  %lf  %d  %lf\n",
                Puri->info.homo_eigenvector_is_computed_in_iter,
                (double)Puri->info.eigValHOMO,
                Puri->info.homo_eigensolver_iter,
                Puri->info.homo_eigensolver_time);

         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "HOMO eigenvector is computed.");
         do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "Perform congruence transformation.");
         // perform congruence transformation
         (*eigVecHOMO) = invCholFactor_tmp * (*eigVecHOMO);
      }
   }    // Note: invCholFactor_tmp goes out of scope


   if (intervalType::intersect(lumoInterval_F_ort_prev_expanded, lumoIntervalNew).empty())
   {
      do_output(LOG_CAT_WARNING, LOG_AREA_DENSFROMF,
                "The intersection of lumoInterval_F_ort_prev_expanded and lumoIntervalNew is empty set!");
   }

   if (intervalType::intersect(homoInterval_F_ort_prev_expanded, homoIntervalNew).empty())
   {
      do_output(LOG_CAT_WARNING, LOG_AREA_DENSFROMF,
                "The intersection of homoInterval_F_ort_prev_expanded and homoIntervalNew is empty set!");
   }

   // Save the improved HOMO/LUMO intervals of F_ort:
   homoInterval_F_ort_prev = homoIntervalNew;
   lumoInterval_F_ort_prev = lumoIntervalNew;

   // Calculate HOMO_LUMO intervals of Finput. We need to expand
   // the F_ort intervals due to the truncation done earlier.
   homoInterval_Finput = homoInterval_F_ort_prev;
   lumoInterval_Finput = lumoInterval_F_ort_prev;

   // Anastasia comment:
   // it may happen that bounds for homo and lumo are very good and intervals are actually empty sets,
   // then the increase() function will throw exception
   if (homoInterval_Finput.low() >= homoInterval_Finput.upp())
   {
      homoInterval_Finput = intervalType(homoInterval_Finput.low(), homoInterval_Finput.low() + 1e-10);
   }
   if (lumoInterval_Finput.low() >= lumoInterval_Finput.upp())
   {
      lumoInterval_Finput = intervalType(lumoInterval_Finput.upp() - 1e-10, lumoInterval_Finput.upp());
   }

   homoInterval_Finput.increase(truncError_1);
   lumoInterval_Finput.increase(truncError_1);

   // Output info about gap.
   ergo_real gapMin = lumoInterval_F_ort_prev.low() - homoInterval_F_ort_prev.upp();
   ergo_real gapMax = lumoInterval_F_ort_prev.upp() - homoInterval_F_ort_prev.low();
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,
             "E(LUMO) - E(HOMO) >= %22.11f = %22.11f eV",
             (double)gapMin, (double)gapMin / UNIT_one_eV);
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,
             "E(LUMO) - E(HOMO) <= %22.11f = %22.11f eV",
             (double)gapMax, (double)gapMax / UNIT_one_eV);
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,
             "HOMO interval : [ %17.12f %17.12f ]",
             (double)homoInterval_F_ort_prev.low(), (double)homoInterval_F_ort_prev.upp());
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,
             "LUMO interval : [ %17.12f %17.12f ]",
             (double)lumoInterval_F_ort_prev.low(), (double)lumoInterval_F_ort_prev.upp());

   puri_stats[stats_prefix + "HOMO_LUMO_gap_lo_eV"] = gapMin / UNIT_one_eV;
   puri_stats[stats_prefix + "HOMO_LUMO_gap_hi_eV"] = gapMax / UNIT_one_eV;


   // we do not need Puri anymore, then delete it
   delete Puri;


   // Check trace of resulting density matrix
   ergo_real trace       = D_w.trace();
   ergo_real wantedTrace = noOfOccupiedOrbs;
   ergo_real traceError  = trace - wantedTrace;
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,
             "Trace of resulting density matrix is %22.11f, error is %18.14f.",
             (double)trace, (double)traceError);
   // Check that relative error in trace is not unreasonably large
   if ((wantedTrace > 0) && (template_blas_fabs(traceError) > wantedTrace / 4))
   {
      throw "Error in get_dens_from_fock_sparse (GetDensFromFock class): traceError is unreasonably large; seems like something went very wrong.";
   }


   // Do truncation to speed up following multiplication operation.
   ergo_real threshold_2 = subspaceThr_2 * (1 - 2 * eigvalueThr_Puri) / (1 + subspaceThr_2);
   // Make sure that eigenvalue movement is not too large:
   threshold_2 = eigvalueThr_2 < threshold_2 ? eigvalueThr_2 : threshold_2;
   size_t nnzD_before_trunc = D_w.nnz();
   double    nnzD_before_trunc_pc = (double)nnzD_before_trunc * 100 / ((double)n * n);
   #ifdef USE_CHUNKS_AND_TASKS
      ergo_real truncError_2 = D_w.thresh_frob(threshold_2);
   #else
      ergo_real truncError_2 = D_w.thresh(threshold_2, truncationNorm);
   #endif
   size_t nnzD_after_trunc = D_w.nnz();
   double    nnzD_after_trunc_pc  = (double)nnzD_after_trunc * 100 / ((double)n * n);
   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,
             "Truncated D_ort ( %s ), selected threshold = %10.6g, "
            "returned error = %10.6g, nnz before = %lu <-> %3.4f %%, nnz after = %lu <-> %3.4f %%",
             mat::getNormTypeString(truncationNorm).c_str(), (double)threshold_2,
             (double)truncError_2, nnzD_before_trunc, nnzD_before_trunc_pc, nnzD_after_trunc, nnzD_after_trunc_pc);
   puri_stats[stats_prefix + "nnz_percentage_D_ort"] = nnzD_after_trunc_pc;

   {
      triangMatrix invCholFactor_tmp(invCholFactor);
      invCholFactor_tmp.readFromFile();
      output_current_memory_usage(LOG_AREA_DENSFROMF, "Before D_w.to_nonnorm_basis");
      Util::TimeMeter timeMeterWriteAndReadAll;
      std::string     sizesStr = mat::FileWritable::writeAndReadAll();
      timeMeterWriteAndReadAll.print(LOG_AREA_DENSFROMF, "FileWritable::writeAndReadAll");
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, ((std::string)"writeAndReadAll sizesStr: '" + sizesStr).c_str());
      Util::TimeMeter timeMeterDortTransf;

      triangMatrixWrap invCholFactor_tmp_w;
      transform_matrix_from_to(invCholFactor_tmp, invCholFactor_tmp_w, params);
      invCholFactor_tmp.clear();
    
      D_w = invCholFactor_tmp_w * D_w * transpose(invCholFactor_tmp_w);

      timeMeterDortTransf.print(LOG_AREA_DENSFROMF, "D_w to non-orthogonal basis");
      output_current_memory_usage(LOG_AREA_DENSFROMF, "After D_w to non-orthogonal basis [D_S = Z * D_ort * ZT]");

#ifndef USE_CHUNKS_AND_TASKS   // eucl_thresh is not implemented for CHT
      // Do truncation again, to reduce memory usage.
      ergo_real threshold_3 = subspaceThr_3 * (1 - 2 * eigvalueThr_Puri - 2 * truncError_2) / (1 + subspaceThr_3);

      //Make sure that eigenvalue movement is not too large:
      threshold_3 = eigvalueThr_3 < threshold_3 ? eigvalueThr_3 : threshold_3;

      //Do truncation, taking into account that we are in 'non-orthogonal basis', passing invCholFactor to thresh
      double    nnzD_S_before_trunc_pc = (double)D_w.nnz() * 100 / ((double)n * n);
      ergo_real truncError_3           = D_w.eucl_thresh(threshold_3, &invCholFactor_tmp_w);
      double    nnzD_S_after_trunc_pc  = (double)D_w.nnz() * 100 / ((double)n * n);
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,
                "Truncated D_S (eucl with Z), selected threshold = %10.6g, returned error = %10.6g, nnz before = %3.4f %%, nnz after = %3.4f %%",
                (double)threshold_3, (double)truncError_3, nnzD_S_before_trunc_pc, nnzD_S_after_trunc_pc);
      puri_stats[stats_prefix + "nnz_percentage_D_S"] = nnzD_S_after_trunc_pc;
#endif
   }

   D_w *= factor;
   transform_matrix_from_to(D_w, resultDens, params);
   D_w.clear();  // clean wrapper

   do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "get_dens_from_fock_sparse ending OK");
   timeMeterTot.print(LOG_AREA_DENSFROMF, "get_dens_from_fock_sparse");

   return 0;
}


/** Function save all needed data to files in order to repeat
 *  recursive expansion in a desired SCF cycle later. The purpose of
 *  the function is mainly testing.
 */
void GetDensFromFock::create_checkpoint(symmMatrix&   Finput,
                                        symmMatrix&   F_ort_prev,
                                        generalVector *eigVecLUMO,
                                        generalVector *eigVecHOMO,
                                        std::string   IDstr)
{
   std::ostringstream checkpoint_ID;
   checkpoint_ID << IDstr << "_" << SCF_step;
   std::cout << "Create checkpoint with ID = " << checkpoint_ID.str() << std::endl;

   /* Save input data */

   std::ostringstream name;
   name << filenameFinput << "_" << checkpoint_ID.str() << ".bin";
   Finput.copyToFile(name.str().c_str());

   name.clear();
   name.str("");
   name << filenameF_ort_prev << "_" << checkpoint_ID.str() << ".bin";
   F_ort_prev.copyToFile(name.str().c_str());


   /* Save class members */

   // matrices

   // Overlap matrix (written to file)
   name.clear();
   name.str("");
   name << filenameOverlap << "_" << checkpoint_ID.str() << ".bin";
   overlapMatrix.copyToFile(name.str().c_str());


   // Inverse Cholesky factor (written to file)
   name.clear();
   name.str("");
   name << filenameinvCholFactor << "_" << checkpoint_ID.str() << ".bin";
   invCholFactor.copyToFile(name.str().c_str());



   // save all numbers of basic arithmetic types
   name.clear();
   name.str("");
   name << file_for_basic_types << "_" << checkpoint_ID.str() << ".data";
   ofstream file_basic(name.str().c_str(), ios::out | ios::app);
   if (!file_basic.is_open())
   {
      throw "GetDensFromFock::create_checkpoint unable open file for writing.";
   }


   file_basic << SCF_step << "\n";
   file_basic << use_diagonalization << "\n";
   file_basic << use_purification << "\n";
   file_basic << store_all_eigenvalues_to_file << "\n";
   file_basic << try_eigv_on_next_iteration_if_fail << "\n";
   file_basic << (double)electronicTemperature << "\n";
   file_basic << std::setprecision(8) << (double)gap_expected_lower_bound << "\n";
   file_basic << std::setprecision(8) << (double)eigvalueErrorLimit << "\n";
   file_basic << std::setprecision(8) << (double)subspaceErrorLimit << "\n";
   file_basic << use_diag_on_error << "\n";
   file_basic << use_diag_on_error_guess << "\n";
   file_basic << create_m_files << "\n";
   file_basic << output_homo_and_lumo_eigenvectors << "\n";
   file_basic << use_prev_vector_as_initial_guess << "\n";
   file_basic << ignore_purification_failure << "\n";
   file_basic << use_rand_perturbation_for_alleigsint << "\n";
   file_basic << stats_prefix << "\n";
   file_basic << plot_puri_results << "\n";
   file_basic << plot_puri_results_str << "\n";
   file_basic << use_acceleration << "\n";
   file_basic << use_new_stopping_criterion << "\n";
   file_basic << eigenvectors_method << "\n";
   file_basic << eigenvectors_iterative_method << "\n";
   file_basic << 1 << "\n"; // old parameter number_of_eigenvalues, remains for the compatibility
   file_basic << std::setprecision(8) << (double)eigensolver_accuracy << "\n";
   file_basic << eigensolver_maxiter << "\n";
   file_basic << n << "\n";
   file_basic << noOfOccupiedOrbs << "\n";
   file_basic << (double)factor << "\n";
   file_basic << std::setprecision(8) << (double)invCholFactor_euclnorm << "\n";
   file_basic << maxMul << "\n";
   file_basic << leavesSizeMax << "\n";
   file_basic << blocksize << "\n";


   if ((output_homo_and_lumo_eigenvectors == SET) && (eigenvectors_method != ""))
   {
      if (!eigVecLUMO->is_empty())
      {
         file_basic << 1 << "\n";
         eigVecLUMO->writeToFile();

         name.clear();
         name.str("");
         name << filenameeigVecLUMO << "_" << checkpoint_ID.str() << ".bin";
         eigVecLUMO->copyToFile(name.str().c_str());

         eigVecLUMO->readFromFile();
      }
      else
      {
         file_basic << 0 << "\n";
      }

      if (!eigVecHOMO->is_empty())
      {
         file_basic << 1 << "\n";
         eigVecHOMO->writeToFile();

         name.clear();
         name.str("");
         name << filenameeigVecHOMO << "_" << checkpoint_ID.str() << ".bin";
         eigVecHOMO->copyToFile(name.str().c_str());

         eigVecHOMO->readFromFile();
      }
      else
      {
         file_basic << 0 << "\n";
      }
    }


   // all non-trivial

   switch (truncationNormPurification)
   {
   case mat::frobNorm:
      file_basic << 1 << "\n";
      break;

   case mat::mixedNorm:
      file_basic << 2 << "\n";
      break;

   case mat::euclNorm:
      file_basic << 3 << "\n";
      break;

   default:
      throw "GetDensFromFock::create_checkpoint unknown truncation norm.";
   }

   switch (stopCriterionNormPurification)
   {
   case mat::frobNorm:
      file_basic << 1 << "\n";
      break;

   case mat::mixedNorm:
      file_basic << 2 << "\n";
      break;

   case mat::euclNorm:
      file_basic << 3 << "\n";
      break;

   default:
      throw "GetDensFromFock::create_checkpoint unknown stopping criterion norm.";
   }

   std::vector<int> blockSizesCopy;
   matrixSizesAndBlocks.getBlockSizeVector(blockSizesCopy);
   file_basic << blockSizesCopy.size() << "\n";
   for (unsigned int i = 0; i < blockSizesCopy.size(); ++i)
   {
      file_basic << blockSizesCopy[i] << "\n";
   }


   file_basic << std::setprecision(16) << (double)homoInterval_Finput.low() << "\n";
   file_basic << std::setprecision(16) << (double)homoInterval_Finput.upp() << "\n";
   file_basic << std::setprecision(16) << (double)lumoInterval_Finput.low() << "\n";
   file_basic << std::setprecision(16) << (double)lumoInterval_Finput.upp() << "\n";

   file_basic << std::setprecision(16) << (double)homoInterval_F_ort_prev.low() << "\n";
   file_basic << std::setprecision(16) << (double)homoInterval_F_ort_prev.upp() << "\n";
   file_basic << std::setprecision(16) << (double)lumoInterval_F_ort_prev.low() << "\n";
   file_basic << std::setprecision(16) << (double)lumoInterval_F_ort_prev.upp() << "\n";

   file_basic.close();
}



inline bool file_exist(const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}


/** Function restores data from files in order to repeat recursive
 *  expansion in a desired SCF cycle. The purpose of the function is
 *  mainly testing.
 */
void GetDensFromFock::restore_from_checkpoint(GetDensFromFock& DensFromFock,
                                              symmMatrix&      Finput,
                                              symmMatrix&      F_ort_prev,
                                              generalVector    *eigVecLUMO,
                                              generalVector    *eigVecHOMO,
                                              std::string      checkpoint_path,
                                              std::string      IDstr,
                                              int              SCF_step)
{
   const char *filenameFinput        = DensFromFock.filenameFinput;
   const char *filenameF_ort_prev    = DensFromFock.filenameF_ort_prev;
   const char *filenameeigVecLUMO    = DensFromFock.filenameeigVecLUMO;
   const char *filenameeigVecHOMO    = DensFromFock.filenameeigVecHOMO;
   const char *filenameOverlap       = DensFromFock.filenameOverlap;
   const char *filenameinvCholFactor = DensFromFock.filenameinvCholFactor;
   const char *file_for_basic_types  = DensFromFock.file_for_basic_types;


   // read all basic data
   std::ostringstream name;
   name << checkpoint_path << "/" << file_for_basic_types << "_" << IDstr << "_" << SCF_step << ".data";
   ifstream file_basic(name.str().c_str(), ios::in);
   if (!file_basic.is_open())
   {
      throw "GetDensFromFock::restore_from_checkpoint unable open file for reading.";
   }

   string tmp;
   double tmp_double;

   file_basic >> DensFromFock.SCF_step;
   file_basic >> DensFromFock.use_diagonalization;
   file_basic >> DensFromFock.use_purification;
   file_basic >> DensFromFock.store_all_eigenvalues_to_file;
   file_basic >> DensFromFock.try_eigv_on_next_iteration_if_fail;
   file_basic >> tmp_double;
   DensFromFock.electronicTemperature = tmp_double;
   file_basic >> tmp_double;
   DensFromFock.gap_expected_lower_bound = tmp_double;
   file_basic >> tmp_double;
   DensFromFock.eigvalueErrorLimit = tmp_double;
   file_basic >> tmp_double;
   DensFromFock.subspaceErrorLimit = tmp_double;
   file_basic >> DensFromFock.use_diag_on_error;
   file_basic >> DensFromFock.use_diag_on_error_guess;
   file_basic >> DensFromFock.create_m_files;
   file_basic >> DensFromFock.output_homo_and_lumo_eigenvectors;
   file_basic >> DensFromFock.use_prev_vector_as_initial_guess;
   file_basic >> DensFromFock.ignore_purification_failure;
   file_basic >> DensFromFock.use_rand_perturbation_for_alleigsint;
   getline(file_basic, tmp);
   getline(file_basic, DensFromFock.stats_prefix);
   //file_basic >> DensFromFock.stats_prefix;
   file_basic >> DensFromFock.plot_puri_results;
   getline(file_basic, tmp);
   getline(file_basic, DensFromFock.plot_puri_results_str);
   //file_basic >> DensFromFock.plot_puri_results_str;
   file_basic >> DensFromFock.use_acceleration;
   file_basic >> DensFromFock.use_new_stopping_criterion;
   getline(file_basic, tmp);
   getline(file_basic, DensFromFock.eigenvectors_method);
   //file_basic >> DensFromFock.eigenvectors_method;
   getline(file_basic, DensFromFock.eigenvectors_iterative_method);
   //file_basic >> DensFromFock.eigenvectors_iterative_method;
   int dummy_number_of_eigenvalues; // removed parameter, remains for the compatibility
   file_basic >> dummy_number_of_eigenvalues;
   file_basic >> tmp_double;
   DensFromFock.eigensolver_accuracy = tmp_double;
   file_basic >> DensFromFock.eigensolver_maxiter;
   file_basic >> DensFromFock.n;
   file_basic >> DensFromFock.noOfOccupiedOrbs;
   file_basic >> tmp_double;
   DensFromFock.factor = tmp_double;
   file_basic >> tmp_double;
   DensFromFock.invCholFactor_euclnorm = tmp_double;
   file_basic >> DensFromFock.maxMul;
   file_basic >> DensFromFock.leavesSizeMax;
   file_basic >> DensFromFock.blocksize;

   int vector_lumo_not_null, vector_homo_not_null;
   file_basic >> vector_lumo_not_null;
   file_basic >> vector_homo_not_null;

   int norm_trunc_ID;
   file_basic >> norm_trunc_ID;
   switch (norm_trunc_ID)
   {
   case 1:
      DensFromFock.truncationNormPurification = mat::frobNorm;
      break;

   case 2:
      DensFromFock.truncationNormPurification = mat::mixedNorm;
      break;

   case 3:
      DensFromFock.truncationNormPurification = mat::euclNorm;
      break;

   default:
      throw "GetDensFromFock::restore_from_checkpoint unknown truncation norm.";
   }


   int norm_st_crit_ID;
   file_basic >> norm_st_crit_ID;
   switch (norm_st_crit_ID)
   {
   case 1:
      DensFromFock.stopCriterionNormPurification = mat::frobNorm;
      break;

   case 2:
      DensFromFock.stopCriterionNormPurification = mat::mixedNorm;
      break;

   case 3:
      DensFromFock.stopCriterionNormPurification = mat::euclNorm;
      break;

   default:
      throw "GetDensFromFock::restore_from_checkpoint unknown stopping criterion norm.";
   }


   int blockSizes_size;
   file_basic >> blockSizes_size;
   std::vector<int> blockSizes(blockSizes_size);
   for (int i = 0; i < blockSizes_size; ++i)
   {
      file_basic >> blockSizes[i];
   }

   DensFromFock.matrixSizesAndBlocks = mat::SizesAndBlocks(blockSizes, DensFromFock.n);

   // eigenvalue bounds
   double lower, upper; // Elias note: changed from ergo_real to double here to make it work when ergo_real is __float128
   file_basic >> lower >> upper;
   DensFromFock.homoInterval_Finput = intervalType(lower, upper);
   file_basic >> lower >> upper;
   DensFromFock.lumoInterval_Finput = intervalType(lower, upper);

   file_basic >> lower >> upper;
   DensFromFock.homoInterval_F_ort_prev = intervalType(lower, upper);
   file_basic >> lower >> upper;
   DensFromFock.lumoInterval_F_ort_prev = intervalType(lower, upper);



   file_basic.close();

   cout << "Finished to read text data. Starting with binary data." << endl;


   // initialize all matrices and vectors with corresponding block dimensions

   name.clear();
   name.str("");
   name << checkpoint_path << "/" << filenameFinput << "_" << IDstr << "_" << SCF_step << ".bin";

   if(!file_exist(name.str())) throw std::runtime_error("File " + name.str() + " does not exist!");

   Finput.resetSizesAndBlocks(DensFromFock.matrixSizesAndBlocks, DensFromFock.matrixSizesAndBlocks); // set data structure
   Finput.copyFromFile(name.str().c_str());

   name.clear();
   name.str("");
   name << checkpoint_path << "/" << filenameF_ort_prev << "_" << IDstr << "_" << SCF_step << ".bin";

   if(!file_exist(name.str())) throw std::runtime_error("File " + name.str() + " does not exist!");

   F_ort_prev.resetSizesAndBlocks(DensFromFock.matrixSizesAndBlocks, DensFromFock.matrixSizesAndBlocks); // set data structure
   F_ort_prev.copyFromFile(name.str().c_str());


   name.clear();
   name.str("");
   name << checkpoint_path << "/" << filenameOverlap << "_" << IDstr << "_" << SCF_step << ".bin";

   if(!file_exist(name.str())) throw std::runtime_error("File " + name.str() + " does not exist!");

   DensFromFock.overlapMatrix.resetSizesAndBlocks(DensFromFock.matrixSizesAndBlocks, DensFromFock.matrixSizesAndBlocks); // set data structure
   DensFromFock.overlapMatrix.copyFromFile(name.str().c_str());

   name.clear();
   name.str("");
   name << checkpoint_path << "/" << filenameinvCholFactor << "_" << IDstr << "_" << SCF_step << ".bin";

   if(!file_exist(name.str())) throw std::runtime_error("File " + name.str() + " does not exist!");

   DensFromFock.invCholFactor.resetSizesAndBlocks(DensFromFock.matrixSizesAndBlocks, DensFromFock.matrixSizesAndBlocks); // set data structure
   DensFromFock.invCholFactor.copyFromFile(name.str().c_str());

   // HOMO and LUMO eigenvectors

   if (vector_lumo_not_null)
   {
      eigVecLUMO->resetSizesAndBlocks(DensFromFock.matrixSizesAndBlocks);// set data structure
      name.clear();
      name.str("");
      name << checkpoint_path << "/" << filenameeigVecLUMO << "_" << IDstr << "_" << SCF_step << ".bin";

      if(!file_exist(name.str())) throw std::runtime_error("File " + name.str() + " does not exist!");

      eigVecLUMO->copyFromFile(name.str().c_str());
      eigVecLUMO->readFromFile();
   }

   if (vector_homo_not_null)
   {
      eigVecHOMO->resetSizesAndBlocks(DensFromFock.matrixSizesAndBlocks);// set data structure
      name.clear();
      name.str("");
      name << checkpoint_path << "/" << filenameeigVecHOMO << "_" << IDstr << "_" << SCF_step << ".bin";

      if(!file_exist(name.str())) throw std::runtime_error("File " + name.str() + " does not exist!");

      eigVecHOMO->copyFromFile(name.str().c_str());
      eigVecHOMO->readFromFile();
   }

}
