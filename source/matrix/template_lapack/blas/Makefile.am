noinst_PROGRAMS = templateblastest templateblastest_threaded

noinst_LIBRARIES = libtemplateblas.a

templateblastest_SOURCES = template_blas_test.cc
templateblastest_LDADD   = libtemplateblas.a $(FLIBS)

templateblastest_threaded_SOURCES = template_blas_test_threaded.cc
templateblastest_threaded_LDADD   = libtemplateblas.a $(FLIBS)

libtemplateblas_a_SOURCES = \
template_blas_basicmath.h	\
template_blas_basicmath.cc	\
template_blas_num_limits.h	\
template_blas_gemm.h		\
template_blas_asum.h		\
template_blas_axpy.h		\
template_blas_common.cc 	\
template_blas_common.h		\
template_blas_copy.h		\
template_blas_dot.h		\
template_blas_gemm.h		\
template_blas_gemv.h		\
template_blas_ger.h		\
template_blas_idamax.h		\
template_blas_nrm2.h		\
template_blas_rot.h		\
template_blas_scal.h		\
template_blas_spmv.h		\
template_blas_spr.h		\
template_blas_spr2.h		\
template_blas_swap.h		\
template_blas_symm.h		\
template_blas_symv.h		\
template_blas_syr2.h		\
template_blas_syr2k.h		\
template_blas_syrk.h		\
template_blas_test.cc		\
template_blas_tpmv.h		\
template_blas_tpsv.h		\
template_blas_trmm.h		\
template_blas_trmv.h		\
template_blas_trsm.h		\
template_blas_trsv.h

EXTRA_DIST =  \
test.sh

TESTS_ENVIRONMENT = top_builddir=$(top_builddir) top_srcdir=$(top_srcdir)

TESTS = \
	test.sh
