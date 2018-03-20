noinst_LIBRARIES = libdft.a

libdft_a_SOURCES = \
	aos.cc		\
	aos.h		\
	cubature_rules.cc	\
	cubature_rules.h	\
	dft.h		\
	dft_common.cc	\
	dft_common.h	\
	fun-becke.c	\
	fun-cam.c	\
	fun-gga.c	\
	fun-kt.c	\
	fun-lb94.c	\
	fun-lyp.c	\
	fun-optx.c	\
	fun-p86c.c	\
	fun-pbec.c	\
	fun-pbex.c	\
	fun-pw86x.c	\
	fun-pw91c.c	\
	fun-pw91x.c	\
	fun-pw92c.c	\
	fun-pz81.c	\
	fun-slater.c	\
	fun-tester.cc	\
	fun-vwn.c	\
	functionals.c	\
	functionals.h	\
	grid_atomic.cc	\
	grid_atomic.h	\
	grid_hicu.cc	\
	grid_hicu.h	\
	grid_interface.h\
	grid_matrix.h	\
	grid_params.h	\
	grid_reader.cc	\
	grid_reader.h	\
	grid_stream.cc	\
	grid_stream.h	\
	integrator.cc	\
	integrator.h	\
	lebedev_laikov.cc\
	lebedev_laikov.h\
	lin_trans.cc	\
	lin_trans.h	\
	rho-mat.cc	\
	rho-mat.h	\
	sparse_matrix.cc	\
	sparse_matrix.h		\
	sparse_pattern.cc	\
	sparse_pattern.h	\
	xc_evaluators.h	\
	xc_matrix.cc	\
	xc_matrix.h	\
	xc_matrix_sparse.cc	\
	xc_matrix_sparse.h

# use matrix for gblas...
AM_CPPFLAGS = \
	-I$(top_srcdir)/source			\
	-I$(top_srcdir)/source/basisset		\
	-I$(top_srcdir)/source/utilities	\
	-I$(top_srcdir)/source/utilities_basic	\
	-I$(top_srcdir)/source/integrals	\
	-I$(top_srcdir)/source/matrix		\
        -I$(top_srcdir)/source/matrix/template_lapack/blas     \
        -I$(top_srcdir)/source/matrix/template_lapack/lapack
