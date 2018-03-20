#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(ERGO, 3.6, [info@ergoscf.org])
AM_INIT_AUTOMAKE
# Automake-1.11 option for less verbose output
# configure --disable-silent-rules to restore verbose or run "make V=1"
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

AC_CONFIG_SRCDIR([source/ergo_scripted.cc])
AC_CONFIG_HEADER([source/config.h])

if test "x$prefix" = "xNONE"; then
        prefix=/usr/local
fi

dnl #####################################################################
dnl 2. Options
dnl #####################################################################

AC_MSG_CHECKING(which precision to use)
precision_flag=DOUBLE
AC_ARG_WITH(precision,
        AC_HELP_STRING([--with-precision=(single|double|long|quad)],
                [Which precision should be used (default: double)]),
           [ case "$with_precision" in 
             s*) precision_flag=SINGLE;;
             l*) precision_flag=LONG_DOUBLE;;
             q*) precision_flag=QUAD_FLT128;;
             *)  precision_flag=DOUBLE;;
	     esac])
if test $precision_flag = QUAD_FLT128; then
AC_DEFINE(PRECISION_QUAD_FLT128, 1,[which precision should be used])
ergo_cv_precision=4
ergo_cv_quad_flt128=1
fi
if test $precision_flag = LONG_DOUBLE; then
AC_DEFINE(PRECISION_LONG_DOUBLE, 1,[which precision should be used])
ergo_cv_precision=3
fi
if test $precision_flag = SINGLE; then
AC_DEFINE(PRECISION_SINGLE, 1,[which precision should be used])
ergo_cv_precision=1
fi
if test $precision_flag = DOUBLE; then
AC_DEFINE(PRECISION_DOUBLE, 1,[which precision should be used])
ergo_cv_precision=2
fi
AC_MSG_RESULT($precision_flag)
AC_DEFINE_UNQUOTED(FUNC_PRECISION, $ergo_cv_precision,
                   [XC functionals' precision])

AC_ARG_WITH(goto,
        AC_HELP_STRING([--with-goto=PATH_TO_GOTOLIB],
                [Specify the file containing the Goto's BLAS library.]))

ergo_cv_werror=auto
AC_ARG_ENABLE(critical-warnings,
        AC_HELP_STRING([--disable-critical-warnings],
        [Continue compiling even if warnings are emitted (default=auto)]),
        [if test x$enableval = xno -o x$enableval = xyes; then ergo_cv_werror=$enableval; fi])

ergo_cv_performance=no
AC_ARG_ENABLE(performance,
        AC_HELP_STRING([--enable-performance],
        [Compile for performance]),
        [if test x$enableval = xyes; then ergo_cv_performance=yes; fi])

ergo_cv_linalg_templates=no
AC_ARG_ENABLE(linalgebra-templates,
        AC_HELP_STRING([--enable-linalgebra-templates],
        [Use C++ Linear Algebra Template Library (BLAS+LAPACK) (default=no)]),
        [if test x$enableval = xyes; then ergo_cv_linalg_templates=yes; fi])

ergo_cv_coverage=no
AC_ARG_ENABLE(coverage,
        AC_HELP_STRING([--enable-coverage],
        [Enable source code coverage checking using --coverage compiler option (then info is saved that can be used by the gcov and lcov tools) (default=no)]),
        [if test x$enableval = xyes; then ergo_cv_coverage=yes; fi])

AC_MSG_CHECKING(whether to use internal linear algebra template library)
if test x$ergo_cv_linalg_templates = xyes; then
   if test "$with_goto" != ""; then
      AC_ERROR([--enable-linalgebra-templates conflicts with --with-goto])
   fi
   AC_DEFINE(USE_LINALG_TEMPLATES, 1,[Use C++ Linear Algebra Template Library])
fi
AC_MSG_RESULT($ergo_cv_linalg_templates)

ergo_cv_sse_intrinsics=no
AC_ARG_ENABLE(sse-intrinsics,
        AC_HELP_STRING([--enable-sse-intrinsics],
        [Use SSE intrinsics to generate efficient computational kernels (only possible with some compilers/architectures) (default=no)]),
        [if test x$enableval = xyes; then ergo_cv_sse_intrinsics=yes; fi])

AC_MSG_CHECKING(whether to use SSE intrinsics to generate efficient computational kernels)
if test x$ergo_cv_sse_intrinsics = xyes; then
   AC_DEFINE(USE_SSE_INTRINSICS, 1,[Use SSE intrinsics to generate efficient computational kernels])
fi
AC_MSG_RESULT($ergo_cv_sse_intrinsics)

dnl #####################################################################
dnl 3. Programs: compilers and their options.
dnl #####################################################################

# Checks for programs.
AC_PROG_CC([icc pgcc cc gcc])
AC_PROG_CXX([icpc pgcpp c++ g++])
dnl AC_PROG_F77([ifort efc xlf g77])
AC_PROG_RANLIB
AC_PROG_LEX
AC_PROG_YACC

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h memory.h stddef.h stdlib.h string.h strings.h sys/time.h unistd.h])

if test "$ergo_cv_performance" != "no"; then
	for option in -O3 -std=c99 -restrict -fopenmp; do
                SAVE_CFLAGS="$CFLAGS"
                CFLAGS="$CFLAGS $option"
                AC_MSG_CHECKING([whether gcc understands $option])
                AC_TRY_COMPILE([], [],
                        has_option=yes,
                        has_option=no,)
                if test $has_option = no; then
                  CFLAGS="$SAVE_CFLAGS"
                fi
                AC_MSG_RESULT($has_option)
                unset has_option
                unset SAVE_CFLAGS
        done
        unset option
   AC_LANG_PUSH([C++])
	for option in -O3 -fopenmp; do
                SAVE_CXXFLAGS="$CXXFLAGS"
                CXXFLAGS="$CXXFLAGS $option"
                AC_MSG_CHECKING([whether c++ understands $option])
                AC_TRY_COMPILE([], [],
                        has_option=yes,
                        has_option=no,)
                if test $has_option = no; then
                  CXXFLAGS="$SAVE_CXXFLAGS"
                fi
                AC_MSG_RESULT($has_option)
                unset has_option
                unset SAVE_CFLAGS
        done
        unset option
   AC_LANG_POP([C++])
fi

# Enable --coverage compiler option if requested
if test "$ergo_cv_coverage" != "no"; then
	for option in --coverage ; do
                SAVE_CFLAGS="$CFLAGS"
                CFLAGS="$CFLAGS $option"
                AC_MSG_CHECKING([whether gcc understands $option])
                AC_TRY_COMPILE([], [],
                        has_option=yes,
                        has_option=no,)
                if test $has_option = no; then
                  CFLAGS="$SAVE_CFLAGS"
                fi
                AC_MSG_RESULT($has_option)
                unset has_option
                unset SAVE_CFLAGS
        done
        unset option
   AC_LANG_PUSH([C++])
	for option in --coverage ; do
                SAVE_CXXFLAGS="$CXXFLAGS"
                CXXFLAGS="$CXXFLAGS $option"
                AC_MSG_CHECKING([whether c++ understands $option])
                AC_TRY_COMPILE([], [],
                        has_option=yes,
                        has_option=no,)
                if test $has_option = no; then
                  CXXFLAGS="$SAVE_CXXFLAGS"
                fi
                AC_MSG_RESULT($has_option)
                unset has_option
                unset SAVE_CFLAGS
        done
        unset option
   AC_LANG_POP([C++])
fi


# Checks for typedefs, structures, and compiler characteristics.
dnl AC_C_CONST is problematic with gcc, do not quite know why yet.
AC_C_RESTRICT
AC_C_VOLATILE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([atexit floor memset rint strcasecmp strchr strdup strncasecmp])

# Use AC_CHECK_LIB to get math library before doing AC_CHECK_FUNCS for math functions.
AC_CHECK_LIB([m],[exp])

AC_CHECK_FUNCS([pow powf powl])
AC_CHECK_FUNCS([sqrt sqrtf sqrtl])
AC_CHECK_FUNCS([fabs fabsf fabsl])
AC_CHECK_FUNCS([exp expf expl])
AC_CHECK_FUNCS([log logf logl])
AC_CHECK_FUNCS([log10 log10f log10l])
AC_CHECK_FUNCS([erf erff erfl])
AC_CHECK_FUNCS([erfc erfcf erfcl])
AC_CHECK_FUNCS([sin sinf sinl])
AC_CHECK_FUNCS([cos cosf cosl])

# Make sure we have quadmath functions in $ergo_cv_quad_flt128 case
if test "$ergo_cv_quad_flt128" = "1"; then
  AC_CHECK_LIB([quadmath],[sqrtq])
  AC_CHECK_FUNC(sqrtq,,[AC_ERROR(In QUAD_FLT128 case quadmath library required but it was not found)])
fi

# Check for quadmath versions on math functions
AC_CHECK_FUNCS([powq sqrtq fabsq])
AC_CHECK_FUNCS([expq logq log10q])
AC_CHECK_FUNCS([erfq erfcq])
AC_CHECK_FUNCS([sinq cosq])

# Checks for libraries.
AC_CHECK_LIB([pthread], [pthread_create])

# NOTE: Adding the -Werror flag must be done AFTER checking for library functions.
if test "$ergo_cv_werror" = "auto"; then
  if test -d "$srcdir/.svn"; then
     ergo_cv_werror=yes
  else
     ergo_cv_werror=no
  fi
fi

if test "$ergo_cv_werror" != "no"; then
# Extra warnings for C code...
	for option in -Wsign-compare -Wpointer-arith -Wnested-externs \
            -Wchar-subscripts -Wuninitialized -Wmissing-prototypes    \
            -Wreturn-type -Wdeclaration-after-statement -Werror; do
                SAVE_CFLAGS="$CFLAGS"
                CFLAGS="$CFLAGS $option"
                AC_MSG_CHECKING([whether gcc understands $option])
                AC_TRY_COMPILE([], [],
                        has_option=yes,
                        has_option=no,)
                if test $has_option = no; then
                  CFLAGS="$SAVE_CFLAGS"
                fi
                AC_MSG_RESULT($has_option)
                unset has_option
                unset SAVE_CFLAGS
        done
        unset option

   AC_LANG_PUSH([C++]) # consider adding  -Weffc++
	for option in -Wsign-compare -Wpointer-arith \
            -Wchar-subscripts -Wuninitialized -Woverloaded-virtual    \
            -Wreturn-type -Wctor-dtor-privacy -Wnon-virtual-dtor \
            -Wreorder -Wdeprecated -Wunused -Werror; do
                SAVE_CXXFLAGS="$CXXFLAGS"
                CXXFLAGS="$CXXFLAGS $option"
                AC_MSG_CHECKING([whether c++ understands $option])
                AC_TRY_COMPILE([], [],
                        has_option=yes,
                        has_option=no,)
                if test $has_option = no; then
                  CXXFLAGS="$SAVE_CXXFLAGS"
                fi
                AC_MSG_RESULT($has_option)
                unset has_option
                unset SAVE_CFLAGS
        done
        unset option
   AC_LANG_POP([C++])
fi



# Checks for more libraries...
blas=0
if test "$with_goto" != ""; then
   LIBS="$LIBS $with_goto"
   AC_CHECK_FUNC(dgemm_,,[AC_ERROR([GOTO failed to provide dgemm_])])
   AC_CHECK_FUNC([dstevx_],,
     AC_CHECK_LIB([lapack],
                  [dstevx_],
                  [LIBS="-llapack $LIBS"],
                  [AC_ERROR([GOTO needs lapack. Lapack is not found.])]))
   ergo_la_lib=GOTO
   blas=1
elif test x$ergo_cv_linalg_templates = xyes; then
   ergo_la_lib=internal
else
dnl  AC_F77_LIBRARY_LDFLAGS
  # FIXME: Replace `main' with a function in `-llapack':
  AC_CHECK_FUNC(dposv_,,
  AC_CHECK_LIB([acml],[main],[LIBS="$LIBS -lacml";blas=1;ergo_la_lib=ACML],
    AC_CHECK_LIB([mkl],[main],[LIBS="$LIBS -lmkl";blas=1;ergo_la_lib=MKL],
      AC_CHECK_LIB([essl],[main],[LIBS="$LIBS -lessl";blas=1;ergo_la_lib=ESSL],
        AC_CHECK_LIB(atlas,[main],
          AC_CHECK_LIB(f77blas, dswap_,
                       LIBS="$LIBS -lf77blas -latlas";blas=1;ergo_la_lib=ATLAS,,-latlas),
          AC_CHECK_LIB(blas,[dgemm_],
                       LIBS="$LIBS -lblas $FLIBS";ergo_la_lib=generic,,$FLIBS))))
    AC_CHECK_LIB([mkl_lapack],[main],,
      AC_CHECK_LIB([lapack],[dposv_],,
        dnl ATLAS's lapack may require fortran runtime lib.
        AC_CHECK_LIB([lapack],[dgesv_],LIBS="-llapack $LIBS -lcblas -latlas $FLIBS",
        [AC_ERROR([lapack required but not found])],-lcblas $FLIBS))))
  )
fi
AC_DEFINE_UNQUOTED(USE_BLAS_IN_XC, $blas,
                   [Use BLAS for small matrices in XC integration])
dnl AC_F77_DUMMY_MAIN
dnl AC_F77_LIBRARY_LDFLAGS


# We have headers and basic libs.. Checks for bugs in compilers...
AC_MSG_CHECKING([for 10.0 icpc bug])
AC_RUN_IFELSE([AC_LANG_PROGRAM([[#include <math.h>
#include <string.h>
static int f(int n) { return (n == 0) ?  1 : n*f(n-1); }
static void myfunc(int l) { double a = sin(2.2)*f(l), terms;
  memset(&terms, 0, sizeof(terms)); terms += a; }
]],
[[int l,m;for(l=0; l<=2; l++)  for(m=-l; m<=l; m++) myfunc(l); return 0;]])], 
              [AC_MSG_RESULT([ok])],
              [AC_MSG_FAILURE([found. Change compiler to a working one!])],
              [AC_MSG_RESULT([irrelevant])])



dnl #####################################################################
dnl 7. extra functions.
dnl #####################################################################
sprefix=`cd "$srcdir" && pwd`
AC_DEFINE_UNQUOTED(ERGO_SPREFIX, "$sprefix",
                   [The location of the sources])

ERGO_DATA_PREFIX=`eval "echo ${datadir}/ergo"`
dnl autoconf-2.60 leaves a literal "${prefix}" in ERGO_DATA_PREFIX, so
dnl we clean it once more.
ERGO_DATA_PREFIX=`eval "echo ${ERGO_DATA_PREFIX}"`
AC_DEFINE_UNQUOTED(ERGO_DATA_PREFIX, "$ERGO_DATA_PREFIX",
                   [The location of installed ergo data files])
AC_SUBST(ERGO_DATA_PREFIX)

AC_CONFIG_FILES([
Makefile
source/Makefile
source/basisset/Makefile
source/ci/Makefile
source/dft/Makefile
source/electron_dynamics/Makefile
source/integrals/Makefile
source/matrix/Makefile
source/matrix/template_lapack/Makefile
source/matrix/template_lapack/blas/Makefile
source/matrix/template_lapack/lapack/Makefile
source/densfromf/Makefile
source/densfromf/recursive_expansion/Makefile
source/densfromf/recursive_expansion/src/Makefile
source/densfromf/recursive_expansion/eigenvectors/Makefile
source/densfromf/recursive_expansion/eigenvectors/test/Makefile
source/utilities_basic/file_tools/Makefile
source/utilities_basic/file_tools/src/Makefile
source/utilities_basic/file_tools/test/Makefile
source/scf/Makefile
source/test/Makefile
source/utilities/Makefile
source/utilities_basic/Makefile
basis/Makefile
documentation/source/Makefile
params/Makefile
mol/Makefile
test/Makefile
])

AC_OUTPUT

echo ""
echo "Selected precision : $precision_flag"
echo "Critical warnings  : $ergo_cv_werror"
echo "Linear Algebra Lib : $ergo_la_lib"
echo ""
