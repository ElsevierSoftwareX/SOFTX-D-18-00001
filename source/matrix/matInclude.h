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

/** @file matInclude.h 
 *
 * Copyright(c) Emanuel Rubensson 2006
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date October 2006
 *
 */
#ifndef MAT_MATINCLUDE
#define MAT_MATINCLUDE
#include <iostream>
#include <vector>
#include <fstream>
#include <ios>
#include <cassert>
#include <ctime>
#include <limits>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

/* We need to include config.h to get the PRECISION_QUAD_FLT128 and USE_SSE_INTRINSICS flags. */
#include "config.h"

#include "template_blas_num_limits.h"

#include "Failure.h"
#include "DebugPolicies.h"
#include "SizesAndBlocks.h"
#include "Memory_buffer_thread.h"

#ifdef _OPENMP
#define MAT_OMP_INIT enum omp_failType {noFail = 0, standardFail, runtimeFail, matFail}; \
  volatile omp_failType omp_fail = noFail;				\
  std::exception omp_exce;						\
  std::runtime_error omp_runtime("");					\
  Failure omp_matFail;							\
  omp_set_nested(true);
// if (omp_fail == noFail) {
#define MAT_OMP_START try { 
#define MAT_OMP_END }						\
  catch(Failure & omp_fail_caught) {				\
    omp_fail = matFail;       omp_matFail = omp_fail_caught; }	\
  catch(std::runtime_error & omp_runtime_caught) {		\
    omp_fail = runtimeFail;  omp_runtime = omp_runtime_caught; } \
  catch(std::exception & omp_exce_caught) {			\
    omp_fail = standardFail;  omp_exce = omp_exce_caught; \
} 
#define MAT_OMP_FINALIZE if(omp_fail)					\
    { std::cerr<<"Exception was thrown in OpenMP parallel region\n";	\
      switch (omp_fail) {						\
      case standardFail: throw omp_exce; break;				\
      case  runtimeFail: throw omp_runtime; break;			\
      case      matFail: throw omp_matFail; break;			\
      default: throw Failure("Odd error in omp parallel loop\n");}	\
    }
#else
#define MAT_OMP_INIT
#define MAT_OMP_START
#define MAT_OMP_END
#define MAT_OMP_FINALIZE
#endif

namespace mat{
  class Params {
  protected:
#ifdef _OPENMP
    static unsigned int nProcs;
    static unsigned int matrixParallelLevel; 
#endif
  public:
    static unsigned int getNProcs() {
#ifdef _OPENMP
      if (nProcs == 0)
	throw Failure("mat::Params::getNProcs(): nProcs == 0 Forgot to call setNProcs()?");
      return nProcs;
#else
      return 1;
#endif
    }
    static void setNProcs(unsigned int const nP) {
#ifdef _OPENMP
      nProcs = nP;
#ifdef USE_SSE_INTRINSICS
      Memory_buffer_thread::instance().init_buffers(nProcs);
#endif
#endif
    }
    static unsigned int getMatrixParallelLevel() {
#ifdef _OPENMP
      if (matrixParallelLevel == 0)
	throw Failure("mat::Params::getMatrixParallelLevel(): matrixParallelLevel == 0 Forgot to call setMatrixParallelLevel()?");
      return matrixParallelLevel;
#else
      return 0;
#endif      
    }
    static void setMatrixParallelLevel(unsigned int const mPL) {
#ifdef _OPENMP
      matrixParallelLevel = mPL;
#endif
    }    
  };



  enum property {zero, ful};
  enum normType {frobNorm, euclNorm, mixedNorm};
  normType getNormType(const char* normStr);
  std::string getNormTypeString(normType nType);

  /* getMachineEpsilon(): function for getting the machine epsilon (the
     difference between 1 and the least value greater than 1 that is
     representable) for the given floating-point type.  */
  template<typename Treal>
    inline static Treal getMachineEpsilon() {
    return template_blas_get_machine_epsilon<Treal>();
  }

  template<typename Treal>
    inline static Treal getRelPrecision() {
    return getMachineEpsilon<Treal>();
  }

  class Time {
    static double get_wall_seconds();
    double ticTime;
  public:
    Time();
    void tic();
    float toc();
  };

  class MemUsage {
  private:
    static int getNumberFromBuffer(const char* buffer, const char* s);
  public:
    struct Values {
      float res;
      float virt;
      float peak;
      Values() : res(0), virt(0), peak(0) { }
    };
    static void getMemUsage(Values & values);
  };

} /* end namespace mat */
#endif
