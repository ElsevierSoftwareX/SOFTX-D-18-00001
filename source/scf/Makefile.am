noinst_LIBRARIES = libscf.a

libscf_a_SOURCES = \
	diis_general.cc	\
	diis_general.h		\
	diis_restricted.cc	\
	diis_restricted.h	\
	diis_unrestricted.cc	\
	diis_unrestricted.h	\
	scf_utils.cc		\
	scf_utils.h		\
	SCF_general.cc		\
	SCF_general.h		\
	SCF_restricted.cc	\
	SCF_restricted.h	\
	SCF_unrestricted.cc	\
	SCF_unrestricted.h	\
	SCF_statistics.h	\
	SCF_statistics.cc	\
	scf.cc			\
	scf.h

AM_CPPFLAGS = \
	-I$(top_srcdir)/source			\
	-I$(top_srcdir)/source/basisset		\
	-I$(top_srcdir)/source/dft		\
	-I$(top_srcdir)/source/matrix		\
	-I$(top_srcdir)/source/matrix_wrapper		\
        -I$(top_srcdir)/source/matrix/template_lapack/blas	\
        -I$(top_srcdir)/source/matrix/template_lapack/lapack	\
	-I$(top_srcdir)/source/densfromf	\
	-I$(top_srcdir)/source/matrix_wrapper/CHTMatrix/src \
	-I$(top_srcdir)/source/utilities	\
	-I$(top_srcdir)/source/utilities_basic	\
	-I$(top_srcdir)/source/integrals
