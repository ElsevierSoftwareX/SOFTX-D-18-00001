noinst_LIBRARIES = libelectron_dynamics.a

libelectron_dynamics_a_SOURCES = \
	electron_dynamics.h	\
	tdhf_dynamics.h		\
	tdhf_dynamics.cc

AM_CPPFLAGS = \
	-I$(top_srcdir)/source			\
	-I$(top_srcdir)/source/utilities	\
	-I$(top_srcdir)/source/utilities_basic	\
	-I$(top_srcdir)/source/basisset		\
	-I$(top_srcdir)/source/integrals	\
	-I$(top_srcdir)/source/scf		\
	-I$(top_srcdir)/source/dft		\
        -I$(top_srcdir)/source/matrix     	\
        -I$(top_srcdir)/source/matrix/template_lapack/blas     \
        -I$(top_srcdir)/source/matrix/template_lapack/lapack
