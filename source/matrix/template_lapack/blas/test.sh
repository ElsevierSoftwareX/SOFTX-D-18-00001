#!/bin/sh

if test "$top_builddir" = ""; then
    top_builddir=..
fi
if test "$top_srcdir" = ""; then
    top_srcdir=..
fi

# Prefer gawk - we know exactly what it can do.
# awk on Sun does not support functions, need to use nawk for this
if gawk '{print 1}'</dev/null > /dev/null 2>&1; then
   AWK=gawk
elif nawk '{print 1}'</dev/null > /dev/null 2>&1; then
   AWK=nawk
else
   AWK=awk
fi


echo
echo Testing template_blas
"$top_builddir"/source/matrix/template_lapack/blas/templateblastest 
if [ $? -eq 0 ]
then
echo template_blas test OK
else
echo ERROR in template_blas test
exit 1
fi

echo
echo Testing template_blas with threads
"$top_builddir"/source/matrix/template_lapack/blas/templateblastest_threaded
if [ $? -eq 0 ]
then
echo template_blas with threads test OK
else
echo ERROR in template_blas with threads test
exit 1
fi

