noinst_LIBRARIES = libbasisset.a

libbasisset_a_SOURCES = \
	basisinfo.cc	\
	basisinfo.h	\
	basisset.cc	\
	basisset.h

AM_CPPFLAGS = \
	-I$(top_srcdir)/source			\
	-I$(top_srcdir)/source/utilities	\
	-I$(top_srcdir)/source/utilities_basic	\
	-I$(top_srcdir)/source/integrals	\
	-I$(top_srcdir)/source/matrix		\
        -I$(top_srcdir)/source/matrix/template_lapack/blas     \
        -I$(top_srcdir)/source/matrix/template_lapack/lapack
