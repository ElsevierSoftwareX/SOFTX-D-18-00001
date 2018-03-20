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

/** @file mat_gblas.h C++ interface to a subset of BLAS and LAPACK
 *
 * This file contains an interface to BLAS and LAPACK routines 
 * which makes it easy to use different precision. Currently single    
 * and double precision is supported. One could also implement      
 * specializations for long double without having to change         
 * any other part in the program that uses the routines below.      
 * It is also possible to use different precision within the same   
 * program without having to recompile the entire library.          
 *                                                                  
 * Copyright(c) Emanuel Rubensson 2005
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date September 2005
 */
#ifndef MAT_GBLAS_HEADER
#define MAT_GBLAS_HEADER
#include <ctime>
#include "Failure.h"

/* We need to include config.h to get USE_LINALG_TEMPLATES and USE_SSE_INTRINSICS flags. */
#include "config.h"

#include "template_lapack_common.h"

#ifdef USE_SSE_INTRINSICS
#include "Memory_buffer_thread.h"
#include "gemm_sse/gemm_sse.h"
#endif

/* LEVEL 3 */
extern "C" void dgemm_(const char *ta,const char *tb,
		       const int *n, const int *k, const int *l,
		       const double *alpha,const double *A,const int *lda,
		       const double *B, const int *ldb,
		       const double *beta, double *C, const int *ldc);
extern "C" void dpptrf_(const char *uplo,const int *n, double* ap, int *info);
extern "C" void dspgst_(const int *itype, const char *uplo,const int *n,
			double* ap,const double *bp,int *info);
extern "C" void dtptri_(const char *uplo,const char *diag,const int *n,
			double* ap,int *info);
/* unit triangular means that a value of 1.0 is assumed  */
/* for the diagonal elements (hence diagonal not stored in packed format) */
extern "C" void dtrmm_(const char *side,const char *uplo,const char *transa,
		       const char *diag,const int *m,const int *n,
		       const double *alpha,const double *A,const int *lda,
		       double *B,const int *ldb);
extern "C" void dsygv_(const int *itype,const char *jobz,
		       const char *uplo,const int *n,
		       double *A,const int *lda,double *B,const int *ldb,
		       double* w,double* work,const int *lwork,int *info);
extern "C" void dggev_(const char *jobbl, const char *jobvr, const int *n,
                       double *A, const int *lda, double *B, const int *ldb,
                       double *alphar, double *alphai, double *beta,
                       double *vl, const int *ldvl,
                       double *vr, const int *ldvr,
                       double *work, const int *lwork, int *info);
extern "C" void dpotrf_(const char *uplo, const int *n, double *A, 
			const int *lda, int *info);
extern "C" void dtrtri_(const char *uplo,const char *diag,const int *n,
			double *A, const int *lda, int *info);
extern "C" void dsyrk_(const char *uplo, const char *trans, const int *n, 
		       const int *k, const double *alpha, const double *A, 
		       const int *lda, const double *beta, 
		       double *C, const int *ldc);
extern "C" void dsymm_(const char *side,const char *uplo,
		       const int *m,const int *n,
		       const double *alpha,const double *A,const int *lda,
		       const double *B,const int *ldb, const double* beta,
		       double *C,const int *ldc);
extern "C" void dpocon_(const char *uplo, const int *n, const double *A,
			const int *lda, const double *anorm, double *rcond,
			double *work, int *iwork, int *info);
extern "C" void dstevx_(const char *jobz, const char *range, const int *n, 
			double *d, double *e, const double *vl, 
			const double *vu, const int *il, const int *iu, 
			const double *abstol, int *m, double *w, double *z, 
			const int *ldz, double *work, int *iwork, int *ifail, 
			int *info);
extern "C" void dstevr_(const char *jobz, const char *range, const int *n, 
			double *d, double *e, const double *vl, 
			const double *vu, const int *il, const int *iu, 
			const double *abstol, int *m, double *w, double *z, 
			const int *ldz, int* isuppz, double *work, int* lwork, 
			int *iwork, int* liwork, int *info);
extern "C" void dsyev_(const char *jobz, const char *uplo, const int *n, 
		       double *a, const int *lda, double *w, double *work, 
		       const int *lwork, int *info);

/* LEVEL 2 */
extern "C" void dgemv_(const char *ta, const int *m, const int *n, 
		       const double *alpha, const double *A, const int *lda, 
		       const double *x, const int *incx, const double *beta, 
		       double *y, const int *incy);
extern "C" void dsymv_(const char *uplo, const int *n, 
		       const double *alpha, const double *A, const int *lda, 
		       const double *x, const int *incx, const double *beta, 
		       double *y, const int *incy);
extern "C" void dtrmv_(const char *uplo, const char *trans, const char *diag,
		       const int *n, const double *A, const int *lda, 
		       double *x, const int *incx);
/* LEVEL 1 */
extern "C" void dscal_(const int* n,const double* da, double* dx,
		       const int* incx);
extern "C" double ddot_(const int* n, const double* dx, const int* incx, 
			const double* dy, const int* incy);
extern "C" void daxpy_(const int* n, const double* da, const double* dx,
		       const int* incx, double* dy,const int* incy);

/* Single precision */
/* LEVEL 3 */
extern "C" void sgemm_(const char *ta,const char *tb,
		       const int *n, const int *k, const int *l,
		       const float *alpha,const float *A,const int *lda,
		       const float *B, const int *ldb,
		       const float *beta, float *C, const int *ldc);
extern "C" void spptrf_(const char *uplo,const int *n, float* ap, int *info);
extern "C" void sspgst_(const int *itype, const char *uplo,const int *n,
			float* ap,const float *bp,int *info);
extern "C" void stptri_(const char *uplo,const char *diag,const int *n,
			float* ap,int *info);
/* unit triangular means that a value of 1.0 is assumed  */
/* for the diagonal elements (hence diagonal not stored in packed format) */
extern "C" void strmm_(const char *side,const char *uplo,const char *transa,
		       const char *diag,const int *m,const int *n,
		       const float *alpha,const float *A,const int *lda,
		       float *B,const int *ldb);
extern "C" void ssygv_(const int *itype,const char *jobz,
		       const char *uplo,const int *n,
		       float *A,const int *lda,float *B,const int *ldb,
		       float* w,float* work,const int *lwork,int *info);
extern "C" void sggev_(const char *jobbl, const char *jobvr, const int *n,
                       float *A, const int *lda, float *B, const int *ldb,
                       float *alphar, float *alphai, float *beta,
                       float *vl, const int *ldvl,
                       float *vr, const int *ldvr,
                       float *work, const int *lwork, int *info);
extern "C" void spotrf_(const char *uplo, const int *n, float *A, 
			const int *lda, int *info);
extern "C" void strtri_(const char *uplo,const char *diag,const int *n,
			float *A, const int *lda, int *info);
extern "C" void ssyrk_(const char *uplo, const char *trans, const int *n, 
		       const int *k, const float *alpha, const float *A, 
		       const int *lda, const float *beta, 
		       float *C, const int *ldc);
extern "C" void ssymm_(const char *side,const char *uplo,
		       const int *m,const int *n,
		       const float *alpha,const float *A,const int *lda,
		       const float *B,const int *ldb, const float* beta,
		       float *C,const int *ldc);
extern "C" void spocon_(const char *uplo, const int *n, const float *A,
			const int *lda, const float *anorm, float *rcond,
			float *work, int *iwork, int *info);
extern "C" void sstevx_(const char *jobz, const char *range, const int *n, 
			float *d, float *e, const float *vl, 
			const float *vu, const int *il, const int *iu, 
			const float *abstol, int *m, float *w, float *z, 
			const int *ldz, float *work, int *iwork, int *ifail, 
			int *info);
extern "C" void sstevr_(const char *jobz, const char *range, const int *n, 
			float *d, float *e, const float *vl, 
			const float *vu, const int *il, const int *iu, 
			const float *abstol, int *m, float *w, float *z, 
			const int *ldz, int* isuppz, float *work, int* lwork, 
			int *iwork, int* liwork, int *info);
extern "C" void ssyev_(const char *jobz, const char *uplo, const int *n, 
		       float *a, const int *lda, float *w, float *work, 
		       const int *lwork, int *info);

/* LEVEL 2 */
extern "C" void sgemv_(const char *ta, const int *m, const int *n, 
		       const float *alpha, const float *A, const int *lda, 
		       const float *x, const int *incx, const float *beta, 
		       float *y, const int *incy);
extern "C" void ssymv_(const char *uplo, const int *n, 
		       const float *alpha, const float *A, const int *lda, 
		       const float *x, const int *incx, const float *beta, 
		       float *y, const int *incy);
extern "C" void strmv_(const char *uplo, const char *trans, const char *diag,
		       const int *n, const float *A, const int *lda, 
		       float *x, const int *incx);
/* LEVEL 1 */
extern "C" void sscal_(const int* n,const float* da, float* dx,
		       const int* incx);
#if 0
// sdot_ is unreliable because of varying return type in different 
// implementations. We therefore always use template dot for single precision 
extern "C" double sdot_(const int* n, const float* dx, const int* incx, 
		       const float* dy, const int* incy);
#endif
extern "C" void saxpy_(const int* n, const float* da, const float* dx,
		       const int* incx, float* dy,const int* incy);

namespace mat
{
  struct Gblas {
    static float time;
    static bool timekeeping;
  };
  
  /*************** Default version throws exception */
  template<class T>
    inline static void gemm(const char *ta,const char *tb,
			    const int *n, const int *k, const int *l,
			    const T *alpha,const T *A,const int *lda,
			    const T *B, const int *ldb,
			    const T *beta,T *C, const int *ldc) {
#ifdef USE_SSE_INTRINSICS
    if (*ta == 'N' && *tb == 'N' && *n == 32 && *k == 32 && *l == 32 && *alpha == 1.0 && *beta == 1) {
      static T * A_packed;
      static T * B_packed;
      static T * C_packed;
      int pack_max_size = 10000;  
      if ( (pack_max_size*sizeof(T))%16 )
	throw std::runtime_error("In gblas gemm: requested buffer size not multiple of 16 bytes");
      static T * buffer;
      Memory_buffer_thread::instance().get_buffer(pack_max_size*3, buffer);
      A_packed = buffer;
      B_packed = buffer + pack_max_size;
      C_packed = buffer + 2*pack_max_size;
      gemm_sse<T>(A,B,C,32,32,32,
		  A_packed, B_packed, C_packed, 
		  pack_max_size, pack_max_size, pack_max_size);  
    }
    else
#endif
      template_blas_gemm(ta,tb,n,k,l,alpha,A,lda,B,ldb,beta,C,ldc);
  }

  
  /*  Computes the Cholesky factorization of a symmetric   *
   *  positive definite matrix in packed storage.          */
  template<class T>
    inline static void pptrf(const char *uplo,const int *n, T* ap, int *info) {
    template_lapack_pptrf(uplo,n,ap,info);
  }
  
  template<class T>
    inline static void spgst(const int *itype, const char *uplo,const int *n,
			     T* ap,const T *bp,int *info) {
    template_lapack_spgst(itype,uplo,n,ap,bp,info);
  }

  /* Computes the inverse of a triangular matrix in packed storage. */
  template<class T>
    inline static void tptri(const char *uplo,const char *diag,const int *n,
			     T* ap,int *info) {
    template_lapack_tptri(uplo,diag,n,ap,info);
  }
  
  template<class T>
    inline static void trmm(const char *side,const char *uplo,
			    const char *transa, const char *diag,
			    const int *m,const int *n,
			    const T *alpha,const T *A,const int *lda,
			    T *B,const int *ldb) {
    template_blas_trmm(side,uplo,transa,diag,m,n,alpha,A,lda,B,ldb);
  }
  
  /* Computes all eigenvalues and the eigenvectors of a generalized  *
   * symmetric-definite generalized eigenproblem,                    *
   * Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x.               */
  template<class T>
    inline static void sygv(const int *itype,const char *jobz,
			    const char *uplo,const int *n,
			    T *A,const int *lda,T *B,const int *ldb,
			    T* w,T* work,const int *lwork,int *info) {
    template_lapack_sygv(itype,jobz,uplo,n,A,lda,B,ldb,w,work,lwork,info);
  }

  template<class T>
    inline static void ggev(const char *jobbl, const char *jobvr, 
			    const int *n, T *A, const int *lda, 
			    T *B, const int *ldb, T *alphar, 
			    T *alphai, T *beta, T *vl, 
			    const int *ldvl, T *vr, const int *ldvr,
			    T *work, const int *lwork, int *info) {
    template_lapack_ggev(jobbl, jobvr, n, A, lda, B, ldb, alphar, alphai, beta, vl, 
			 ldvl, vr, ldvr, work, lwork, info);
  }

  /* Computes the Cholesky factorization of a symmetric *
   * positive definite matrix in packed storage.        */
  template<class T>
    inline static void potrf(const char *uplo, const int *n, T *A, 
			     const int *lda, int *info) {
    template_lapack_potrf(uplo, n, A, lda, info);
  }

  /* Computes the inverse of a triangular matrix.                   */
  template<class T>
    inline static void trtri(const char *uplo,const char *diag,const int *n,
			     T *A, const int *lda, int *info) {
    // Create copies of strings because they cannot be const inside trtri.
    char uploCopy[2];
    char diagCopy[2];
    uploCopy[0] = uplo[0];
    uploCopy[1] = '\0';
    diagCopy[0] = diag[0];
    diagCopy[1] = '\0';
    template_lapack_trtri(uploCopy, diagCopy, n, A, lda, info);
  }

  template<class T>
    inline static void syrk(const char *uplo, const char *trans, const int *n, 
			    const int *k, const T *alpha, const T *A, 
			    const int *lda, const T *beta, 
			    T *C, const int *ldc) {
    template_blas_syrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
  }
  
  template<class T>
    inline static void symm(const char *side,const char *uplo,
			    const int *m,const int *n,
			    const T *alpha,const T *A,const int *lda,
			    const T *B,const int *ldb, const T* beta,
			    T *C,const int *ldc) {
    template_blas_symm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
  }
  
  template<class T>
    inline static void pocon(const char *uplo, const int *n, const T *A,
			     const int *lda, const T *anorm, T *rcond,
			     T *work, int *iwork, int *info) {
    template_lapack_pocon(uplo, n, A, lda, anorm, rcond, work, iwork, info);
  }

  template<class T>
    inline static void stevx(const char *jobz, const char *range, 
			     const int *n, T *d, T *e, const T *vl, 
			     const T *vu, const int *il, const int *iu, 
			     const T *abstol, int *m, T *w, T *z, 
			     const int *ldz, T *work, int *iwork, int *ifail, 
			     int *info) {
    template_lapack_stevx(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, 
			  work, iwork, ifail, info);
  }
  
  template<class T>
    inline static void stevr(const char *jobz, const char *range, const int *n, 
			     T *d, T *e, const T *vl, 
			     const T *vu, const int *il, const int *iu, 
			     const T *abstol, int *m, T *w, T *z, 
			     const int *ldz, int* isuppz, T *work, int* lwork, 
			     int *iwork, int* liwork, int *info) {
    template_lapack_stevr(jobz, range, n, d, e, vl, vu, il, iu, abstol, 
			  m, w, z, ldz, isuppz, 
			  work, lwork, iwork, liwork, info);
  }


  template<class T>
    inline static void syev(const char *jobz, const char *uplo, const int *n, 
			    T *a, const int *lda, T *w, T *work, 
			    const int *lwork, int *info) {
    template_lapack_syev(jobz, uplo, n, a, lda, w, work, lwork, info);
  }


  /* LEVEL 2 */
  template<class T>
    inline static void gemv(const char *ta, const int *m, const int *n, 
			    const T *alpha, const T *A, 
			    const int *lda, 
			    const T *x, const int *incx, 
			    const T *beta, T *y, const int *incy) {
    template_blas_gemv(ta, m, n, alpha, A, lda, x, incx, beta, y, incy);
  }

  template<class T>
    inline static void symv(const char *uplo, const int *n, 
			    const T *alpha, const T *A, 
			    const int *lda, const T *x, 
			    const int *incx, const T *beta, 
			    T *y, const int *incy) {
    template_blas_symv(uplo, n, alpha, A, lda, x, incx, beta, y, incy);
  }

  template<class T>
    inline static void trmv(const char *uplo, const char *trans, 
			    const char *diag, const int *n, 
			    const T *A, const int *lda, 
			    T *x, const int *incx) {
    template_blas_trmv(uplo, trans, diag, n, A, lda, x, incx);    
  }


  /* LEVEL 1 */
  template<class T>
    inline static void scal(const int* n,const T* da, T* dx,
			    const int* incx) {
    template_blas_scal(n, da, dx, incx);
  }

  template<class T>
    inline static T dot(const int* n, const T* dx, const int* incx, 
			const T* dy, const int* incy) {
    return template_blas_dot(n, dx, incx, dy, incy);
  }

  template<class T>
    inline static void axpy(const int* n, const T* da, const T* dx,
			    const int* incx, T* dy,const int* incy) {
    template_blas_axpy(n, da, dx, incx, dy, incy);
  }



  
  /* Below follows specializations for double, single, etc. 
     These specializations are not needed if template_blas and template_lapack are used,
     so in that case we skip this entire section. */
#ifndef USE_LINALG_TEMPLATES


  /*************** Double specialization */
  template<>
    inline void gemm<double>(const char *ta,const char *tb,
			     const int *n, const int *k, const int *l,
			     const double *alpha,
			     const double *A,const int *lda,
			     const double *B, const int *ldb,
			     const double *beta,
			     double *C, const int *ldc) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      dgemm_(ta,tb,n,k,l,alpha,A,lda,B,ldb,beta,C,ldc);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      dgemm_(ta,tb,n,k,l,alpha,A,lda,B,ldb,beta,C,ldc);
    }
  }
  
  template<>
    inline  void pptrf<double>(const char *uplo,const int *n, 
			       double* ap, int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      dpptrf_(uplo,n,ap,info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      dpptrf_(uplo,n,ap,info);
    }
  }  

  template<>
    inline  void spgst<double>(const int *itype, const char *uplo,
			       const int *n,
			       double* ap,const double *bp,int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      dspgst_(itype,uplo,n,ap,bp,info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      dspgst_(itype,uplo,n,ap,bp,info);	
    }
  } 
  
  template<>
    inline  void tptri<double>(const char *uplo,const char *diag,const int *n,
			       double* ap,int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      dtptri_(uplo,diag,n,ap,info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      dtptri_(uplo,diag,n,ap,info);
    }
  }
  
  template<>
    inline  void trmm<double>(const char *side,const char *uplo,
			      const char *transa,
			      const char *diag,const int *m,const int *n,
			      const double *alpha,
			      const double *A,const int *lda,
			      double *B,const int *ldb) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      dtrmm_(side,uplo,transa,diag,m,n,alpha,A,lda,B,ldb);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      dtrmm_(side,uplo,transa,diag,m,n,alpha,A,lda,B,ldb);
    }
  }

  template<>
    inline  void sygv<double>(const int *itype,const char *jobz,
			      const char *uplo,const int *n,
			      double *A,const int *lda,
			      double *B,const int *ldb,
			      double* w,double* work,
			      const int *lwork,int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      dsygv_(itype,jobz,uplo,n,A,lda,B,ldb,w,work,lwork,info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      dsygv_(itype,jobz,uplo,n,A,lda,B,ldb,w,work,lwork,info);	
    }
  }

  template<>
    inline void ggev<double>(const char *jobbl, const char *jobvr, 
			     const int *n, double *A, const int *lda, 
			     double *B, const int *ldb, double *alphar, 
			     double *alphai, double *beta, double *vl, 
			     const int *ldvl, double *vr, const int *ldvr,
			     double *work, const int *lwork, int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      dggev_(jobbl, jobvr, n, A, lda, B, ldb, alphar, alphai, beta, vl, 
	     ldvl, vr, ldvr, work, lwork, info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      dggev_(jobbl, jobvr, n, A, lda, B, ldb, alphar, alphai, beta, vl, 
	     ldvl, vr, ldvr, work, lwork, info);
    }
  }


  template<>
    inline void potrf<double>(const char *uplo, const int *n, double *A, 
			      const int *lda, int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      dpotrf_(uplo, n, A, lda, info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      dpotrf_(uplo, n, A, lda, info);      
    }
  }

  template<>
    inline void trtri<double>(const char *uplo,const char *diag,const int *n,
			      double *A, const int *lda, int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      dtrtri_(uplo, diag, n, A, lda, info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      dtrtri_(uplo, diag, n, A, lda, info);     
    }
  }

  template<>
    inline void syrk<double>(const char *uplo, const char *trans, 
			     const int *n, const int *k, const double *alpha, 
			     const double *A, const int *lda, 
			     const double *beta, double *C, const int *ldc) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      dsyrk_(uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      dsyrk_(uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
    }
  }

  template<>
    inline void symm<double>(const char *side,const char *uplo,
			     const int *m,const int *n, const double *alpha,
			     const double *A,const int *lda,
			     const double *B,const int *ldb, 
			     const double* beta,
			     double *C,const int *ldc) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      dsymm_(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      dsymm_(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
    } 
  }

  template<>
    inline void pocon<double>(const char *uplo, const int *n, 
			      const double *A, const int *lda, 
			      const double *anorm, double *rcond, 
			      double *work, int *iwork, int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      dpocon_(uplo, n, A, lda, anorm, rcond, work, iwork, info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      dpocon_(uplo, n, A, lda, anorm, rcond, work, iwork, info);
    }     
  }

  template<>
    inline void stevx<double>(const char *jobz, const char *range, 
			      const int *n, double *d, double *e, 
			      const double *vl, 
			      const double *vu, const int *il, const int *iu, 
			      const double *abstol, int *m, double *w, 
			      double *z, 
			      const int *ldz, double *work, int *iwork, 
			      int *ifail, int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      dstevx_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, 
	      work, iwork, ifail, info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      dstevx_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, 
	      work, iwork, ifail, info);
    }     
  }

  template<>
    inline void stevr<double>(const char *jobz, const char *range, 
			      const int *n, double *d, double *e, 
			      const double *vl, const double *vu, 
			      const int *il, const int *iu, 
			      const double *abstol, 
			      int *m, double *w, 
			      double *z, const int *ldz, int* isuppz, 
			      double *work, int* lwork, 
			      int *iwork, int* liwork, int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      dstevr_(jobz, range, n, d, e, vl, vu, il, iu, abstol, 
	      m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      dstevr_(jobz, range, n, d, e, vl, vu, il, iu, abstol, 
	      m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
    }     
  }


  
  template<>
    inline void syev<double>(const char *jobz, const char *uplo, const int *n, 
			     double *a, const int *lda, double *w, 
			     double *work, const int *lwork, int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      dsyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      dsyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
    }     
  }


  /* LEVEL 2 */
  template<>
    inline void gemv<double>(const char *ta, const int *m, const int *n, 
			     const double *alpha, const double *A, 
			     const int *lda, 
			     const double *x, const int *incx, 
			     const double *beta, double *y, const int *incy) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      dgemv_(ta, m, n, alpha, A, lda, x, incx, beta, y, incy);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      dgemv_(ta, m, n, alpha, A, lda, x, incx, beta, y, incy);
    }     
  }
  
  template<>
    inline void symv<double>(const char *uplo, const int *n, 
			     const double *alpha, const double *A, 
			     const int *lda, const double *x, 
			     const int *incx, const double *beta, 
			     double *y, const int *incy) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      dsymv_(uplo, n, alpha, A, lda, x, incx, beta, y, incy);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      dsymv_(uplo, n, alpha, A, lda, x, incx, beta, y, incy);
    }     
  }

  template<>
    inline void trmv<double>(const char *uplo, const char *trans, 
			     const char *diag, const int *n, 
			     const double *A, const int *lda, 
			     double *x, const int *incx) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      dtrmv_(uplo, trans, diag, n, A, lda, x, incx);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      dtrmv_(uplo, trans, diag, n, A, lda, x, incx);
    }     
  }


  /* LEVEL 1 */
  template<>
    inline void scal<double>(const int* n,const double* da, double* dx,
			     const int* incx) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      dscal_(n, da, dx, incx);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      dscal_(n, da, dx, incx);
    }     
  }
  
  template<>
    inline double dot<double>(const int* n, const double* dx, const int* incx, 
			      const double* dy, const int* incy) {
    double tmp = 0;
    if (Gblas::timekeeping) {
      clock_t start = clock();
      tmp = ddot_(n, dx, incx, dy, incy);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      tmp = ddot_(n, dx, incx, dy, incy);
    }     
    return tmp;
  }
  
  template<>
    inline void axpy<double>(const int* n, const double* da, const double* dx,
			     const int* incx, double* dy,const int* incy) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      daxpy_(n, da, dx, incx, dy, incy);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      daxpy_(n, da, dx, incx, dy, incy);
    }     
  }

  
  /*************** Single specialization */
  template<>
    inline void gemm<float>(const char *ta,const char *tb,
			    const int *n, const int *k, const int *l,
			    const float *alpha,
			    const float *A,const int *lda,
			    const float *B, const int *ldb,
			    const float *beta,
			    float *C, const int *ldc) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      sgemm_(ta,tb,n,k,l,alpha,A,lda,B,ldb,beta,C,ldc);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      sgemm_(ta,tb,n,k,l,alpha,A,lda,B,ldb,beta,C,ldc);
    }
  }
  
  template<>
    inline void pptrf<float>(const char *uplo,const int *n, 
			     float* ap, int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      spptrf_(uplo,n,ap,info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      spptrf_(uplo,n,ap,info);
    }
  }
  
  template<>
    inline  void spgst<float>(const int *itype, const char *uplo,
			      const int *n,
			      float* ap,const float *bp,int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      sspgst_(itype,uplo,n,ap,bp,info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      sspgst_(itype,uplo,n,ap,bp,info);
    }
  } 

  template<>
    inline  void tptri<float>(const char *uplo,const char *diag,
			      const int *n,
			      float* ap,int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      stptri_(uplo,diag,n,ap,info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      stptri_(uplo,diag,n,ap,info);
    }
  }  

  template<>
    inline  void trmm<float>(const char *side,const char *uplo,
			     const char *transa,
			     const char *diag,const int *m,const int *n,
			     const float *alpha,
			     const float *A,const int *lda,
			     float *B,const int *ldb) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      strmm_(side,uplo,transa,diag,m,n,alpha,A,lda,B,ldb);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      strmm_(side,uplo,transa,diag,m,n,alpha,A,lda,B,ldb);
    }
  }
  
  template<>
    inline  void sygv<float>(const int *itype,const char *jobz,
			     const char *uplo,const int *n,
			     float *A,const int *lda,
			     float *B,const int *ldb,
			     float* w,float* work,
			     const int *lwork,int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      ssygv_(itype,jobz,uplo,n,A,lda,B,ldb,w,work,lwork,info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      ssygv_(itype,jobz,uplo,n,A,lda,B,ldb,w,work,lwork,info);
    }
  }

  template<>
    inline void ggev<float>(const char *jobbl, const char *jobvr, 
			    const int *n, float *A, const int *lda, 
			    float *B, const int *ldb, float *alphar, 
			    float *alphai, float *beta, float *vl, 
			    const int *ldvl, float *vr, const int *ldvr,
			    float *work, const int *lwork, int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      sggev_(jobbl, jobvr, n, A, lda, B, ldb, alphar, alphai, beta, vl, 
	     ldvl, vr, ldvr, work, lwork, info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      sggev_(jobbl, jobvr, n, A, lda, B, ldb, alphar, alphai, beta, vl, 
	     ldvl, vr, ldvr, work, lwork, info);
    }
  }


  template<>
    inline void potrf<float>(const char *uplo, const int *n, float *A, 
			     const int *lda, int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      spotrf_(uplo, n, A, lda, info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      spotrf_(uplo, n, A, lda, info);      
    }
  }

  template<>
    inline void trtri<float>(const char *uplo,const char *diag,const int *n,
			     float *A, const int *lda, int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      strtri_(uplo, diag, n, A, lda, info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      strtri_(uplo, diag, n, A, lda, info);     
    }
  }

  template<>
    inline void syrk<float>(const char *uplo, const char *trans, 
			    const int *n, const int *k, const float *alpha, 
			    const float *A, const int *lda, 
			    const float *beta, float *C, const int *ldc) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      ssyrk_(uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      ssyrk_(uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
    }
  }

  template<>
    inline void symm<float>(const char *side,const char *uplo,
			    const int *m,const int *n, const float *alpha,
			    const float *A,const int *lda,
			    const float *B,const int *ldb, 
			    const float* beta,
			    float *C,const int *ldc) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      ssymm_(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      ssymm_(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
    } 
  }

  template<>
    inline void pocon<float>(const char *uplo, const int *n, 
			     const float *A, const int *lda, 
			     const float *anorm, float *rcond, 
			     float *work, int *iwork, int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      spocon_(uplo, n, A, lda, anorm, rcond, work, iwork, info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      spocon_(uplo, n, A, lda, anorm, rcond, work, iwork, info);
    }     
  }

  template<>
    inline void stevx<float>(const char *jobz, const char *range, 
			     const int *n, float *d, float *e, 
			      const float *vl, 
			     const float *vu, const int *il, const int *iu, 
			     const float *abstol, int *m, float *w, 
			     float *z, 
			     const int *ldz, float *work, int *iwork, 
			     int *ifail, int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      sstevx_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, 
	      work, iwork, ifail, info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      sstevx_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, 
	      work, iwork, ifail, info);
    }     
  }

  template<>
    inline void stevr<float>(const char *jobz, const char *range, 
			     const int *n, float *d, float *e, 
			     const float *vl, const float *vu, 
			     const int *il, const int *iu, 
			     const float *abstol, 
			     int *m, float *w, 
			     float *z, const int *ldz, int* isuppz, 
			     float *work, int* lwork, 
			     int *iwork, int* liwork, int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      sstevr_(jobz, range, n, d, e, vl, vu, il, iu, abstol, 
	      m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      sstevr_(jobz, range, n, d, e, vl, vu, il, iu, abstol, 
	      m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
    }     
  }

  template<>
    inline void syev<float>(const char *jobz, const char *uplo, const int *n, 
			    float *a, const int *lda, float *w, 
			    float *work, const int *lwork, int *info) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      ssyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      ssyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
    }     
  }


  /* LEVEL 2 */
  template<>
    inline void gemv<float>(const char *ta, const int *m, const int *n, 
			    const float *alpha, const float *A, 
			    const int *lda, 
			    const float *x, const int *incx, 
			    const float *beta, float *y, const int *incy) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      sgemv_(ta, m, n, alpha, A, lda, x, incx, beta, y, incy);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      sgemv_(ta, m, n, alpha, A, lda, x, incx, beta, y, incy);
    }     
  }

  template<>
    inline void symv<float>(const char *uplo, const int *n, 
			    const float *alpha, const float *A, 
			    const int *lda, const float *x, 
			    const int *incx, const float *beta, 
			    float *y, const int *incy) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      ssymv_(uplo, n, alpha, A, lda, x, incx, beta, y, incy);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      ssymv_(uplo, n, alpha, A, lda, x, incx, beta, y, incy);
    }     
  }

  template<>
    inline void trmv<float>(const char *uplo, const char *trans, 
			     const char *diag, const int *n, 
			     const float *A, const int *lda, 
			     float *x, const int *incx) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      strmv_(uplo, trans, diag, n, A, lda, x, incx);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      strmv_(uplo, trans, diag, n, A, lda, x, incx);
    }     
  }

  /* LEVEL 1 */
  template<>
    inline void scal<float>(const int* n,const float* da, float* dx,
			    const int* incx) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      sscal_(n, da, dx, incx);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      sscal_(n, da, dx, incx);
    }     
  }
#if 0
  // Sdot has different return type in different BLAS implementations 
  // Therefore the specialization for single precision is removed
  template<>
    inline float dot<float>(const int* n, const float* dx, const int* incx, 
			    const float* dy, const int* incy) {
    float tmp;
    if (Gblas::timekeeping) {
      clock_t start = clock();
      sdot_(n, dx, incx, dy, incy);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      sdot_(n, dx, incx, dy, incy);
    }     
    return tmp;
  }
#endif

  template<>
    inline void axpy<float>(const int* n, const float* da, const float* dx,
			     const int* incx, float* dy,const int* incy) {
    if (Gblas::timekeeping) {
      clock_t start = clock();
      saxpy_(n, da, dx, incx, dy, incy);
      Gblas::time +=	((float)(clock() - start)) / (CLOCKS_PER_SEC);
    }
    else {
      saxpy_(n, da, dx, incx, dy, incy);
    }     
  }

  /* END OF SPECIALIZATIONS */
#endif









  /* Other */
  
  template<class Treal>
    static void fulltopacked(const Treal* full, Treal* packed, const int size){
    int pind=0;
    for (int col=0;col<size;col++)
      {
	for(int row=0;row<=col;row++)
	  {
	    packed[pind]=full[col*size+row];
	    pind++;
	  }
      }
  }

  template<class Treal>
    static void packedtofull(const Treal* packed, Treal* full, const int size){
    int psize=(size+1)*size/2;
    int col=0;
    int row=0;
    for(int pind=0;pind<psize;pind++)
      {
	if (col==row)
	  {
	    full[col*size+row]=packed[pind];
	    col++;
	    row=0;
	  }
	else
	  {
	    full[col*size+row]=packed[pind];
	    full[row*size+col]=packed[pind];
	    row++;
	  }
      }
  }

  template<class Treal>
    static void tripackedtofull(const Treal* packed,Treal* full,
				const int size) {
    int psize=(size+1)*size/2;
    int col=0;
    int row=0;
    for(int pind=0;pind<psize;pind++)
      {
	if (col==row)
	  {
	    full[col*size+row]=packed[pind];
	    col++;
	    row=0;
	  }
	else
	  {
	    full[col*size+row]=packed[pind];
	    full[row*size+col]=0;
	    row++;
	  }
      }
  }

  template<class Treal>
    static void trifulltofull(Treal* trifull, const int size) {
    for(int col = 0; col < size - 1; col++) 
      for(int row = col + 1; row < size; row++) 
	trifull[col * size +  row] = 0;
  }
  

} /* namespace mat */

#endif /* GBLAS */
