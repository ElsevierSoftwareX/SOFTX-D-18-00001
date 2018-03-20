SUBDIRS = \
        template_lapack

noinst_PROGRAMS = mattest mattest_lan test_vector omptest matbench matbenchgemm blastime perturbtest 

noinst_LIBRARIES = libmat.a

mattest_SOURCES = API_test.cc
mattest_LDADD   = libmat.a $(FLIBS)     \
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a   \
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a


test_vector_SOURCES = test_vector.cc
test_vector_LDADD   = libmat.a $(FLIBS)     \
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a   \
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a

mattest_lan_SOURCES = test_LanczosSeveralLargestEig.cc
mattest_lan_LDADD   = libmat.a $(FLIBS)     \
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a   \
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a

matbench_SOURCES = bench.cc
matbench_LDADD   = libmat.a $(FLIBS)     \
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a   \
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a

matbenchgemm_SOURCES = bench_gemm_only.cc
matbenchgemm_LDADD   = libmat.a $(FLIBS)     \
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a   \
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a

omptest_SOURCES = testOmpFail.cc	

blastime_SOURCES = blastime.cc
blastime_LDADD   = libmat.a $(FLIBS)     \
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a   \
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a

perturbtest_SOURCES = Perturb_Test.cc
perturbtest_LDADD   = libmat.a $(FLIBS)     \
        $(top_builddir)/source/matrix/template_lapack/blas/libtemplateblas.a   \
        $(top_builddir)/source/matrix/template_lapack/lapack/libtemplatelapack.a

#bsm_test_SOURCES = bsm-test.c
#bsm_test_LDADD   = libqcsparsemat.a

libmat_a_SOURCES = \
MatrixBase.h            \
Failure.h               \
MatrixGeneral.h         \
Matrix.h                \
MatrixHierarchicBase.h  \
Vector.h                \
VectorHierarchicBase.h  \
VectorGeneral.h         \
SizesAndBlocks.h        \
SizesAndBlocks.cc      \
MatrixSymmetric.h       \
TC2.h			\
Perturbation.h          \
MatrixTriangular.h      \
bisection.h		\
mat_gblas.cc            \
mat_gblas.h             \
general.h               \
matrix_proxy.h          \
FileWritable.h          \
FileWritable.cc        \
sort.h                  \
matInclude.h            \
matInclude.cc          \
Allocator.h		\
AllocatorManager.h	\
AllocatorManager.cc	\
allocate.h		\
allocate.cc		\
ValidPtr.h              \
Lanczos.h               \
LanczosLargestMagnitudeEig.h \
LanczosSeveralLargestEig.h \
DebugPolicies.h         \
Interval.h              \
mat_utils.h             \
MatrixTridiagSymmetric.h \
truncation.h            \
Memory_buffer_thread.h  \
Memory_buffer_thread.cc \
gemm_sse/common.h       \
gemm_sse/g_intrin.h     \
gemm_sse/mm_kernel_inner_sse2_A.h \
gemm_sse/mm_kernel_outer_A.h \
gemm_sse/gemm_sse.h \
gemm_sse/vector_intrin.h

EXTRA_DIST =  \
bootstrap.sh            \
atoms.txt               \
natoms_size.txt         \
matrix_D                \
matrix_F                \
matrix_S                \
test.sh

AM_CPPFLAGS = \
        -I$(top_srcdir)/source/matrix/template_lapack/blas           \
        -I$(top_srcdir)/source/matrix/template_lapack/lapack

TESTS_ENVIRONMENT = top_builddir=$(top_builddir) top_srcdir=$(top_srcdir)

TESTS = \
	test.sh
