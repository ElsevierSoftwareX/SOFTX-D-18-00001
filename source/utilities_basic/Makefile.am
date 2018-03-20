SUBDIRS = \
	file_tools 

noinst_LIBRARIES = libutilities_basic.a

libutilities_basic_a_SOURCES = \
	atom_labels.cc			\
	atom_labels.h			\
	barrier.c			\
	barrier.h			\
	csr_matrix.cc			\
	csr_matrix.h			\
	machine_epsilon.cc		\
	machine_epsilon.h		\
	matrix_algebra.cc		\
	matrix_algebra.h		\
	memorymanag.cc			\
	memorymanag.h			\
	molecule.cc			\
	molecule.h			\
	output.cc			\
	output.h			\
	pi.h				\
	realtype.h			\
	polydegree.h			\
	solve_lin_eq_syst.cc		\
	solve_lin_eq_syst.h		\
	units.h				\
	utilities.cc			\
	utilities.h			\
	xyz_file_parser.cc		\
	xyz_file_parser.h

AM_CPPFLAGS = \
	-I$(top_srcdir)/source				\
        -I$(top_srcdir)/source/matrix			\
        -I$(top_srcdir)/source/matrix/template_lapack/blas     \
        -I$(top_srcdir)/source/matrix/template_lapack/lapack
