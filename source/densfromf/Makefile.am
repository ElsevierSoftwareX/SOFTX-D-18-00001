SUBDIRS = \
	recursive_expansion

noinst_LIBRARIES = libdensfromf.a

libdensfromf_a_SOURCES = \
	densfromf_full.cc	\
	densfromf_full.h	\
	GetDensFromFock.h \
	GetDensFromFock.cc	

AM_CPPFLAGS = \
	-I$(top_srcdir)/source			\
	-I$(top_srcdir)/source/basisset		\
	-I$(top_srcdir)/source/integrals	\
	-I$(top_srcdir)/source/utilities	\
	-I$(top_srcdir)/source/utilities_basic	\
	-I$(top_srcdir)/source/matrix		\
        -I$(top_srcdir)/source/matrix/template_lapack/blas           \
        -I$(top_srcdir)/source/matrix/template_lapack/lapack \
	-I$(top_srcdir)/source/densfromf/recursive_expansion \
	-I$(top_srcdir)/source/densfromf/recursive_expansion/src \
	-I$(top_srcdir)/source/densfromf/recursive_expansion/eigenvectors/src	\
	-I$(top_srcdir)/source/matrix_wrapper \
	-I$(top_srcdir)/source/matrix_wrapper/CHTMatrix/src	\
	-I$(top_srcdir)/source/utilities_basic/file_tools/src	
