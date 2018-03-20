noinst_LIBRARIES = libpurification.a

libpurification_a_SOURCES = \
	constants.h 		\
	purification_general.h	\
	purification_general.cc	\
	purification_sp2acc.h	\
	purification_sp2.h	\
	puri_info.h 	\
	puri_info.cc 	\
	random_matrices.cc\
	random_matrices.h	


AM_CPPFLAGS = \
	-I$(top_srcdir)/source			\
	-I$(top_srcdir)/source/basisset		\
	-I$(top_srcdir)/source/integrals	\
	-I$(top_srcdir)/source/utilities	\
	-I$(top_srcdir)/source/utilities_basic	\
	-I$(top_srcdir)/source/matrix		\
        -I$(top_srcdir)/source/matrix/template_lapack/blas           \
        -I$(top_srcdir)/source/matrix/template_lapack/lapack \
	-I$(top_srcdir)/source/densfromf/recursive_expansion/eigenvectors/src	\
	-I$(top_srcdir)/source/utilities_basic/file_tools/src 
