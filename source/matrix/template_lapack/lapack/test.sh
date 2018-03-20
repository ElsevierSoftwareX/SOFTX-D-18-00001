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
echo Testing template_lapack
"$top_builddir"/source/matrix/template_lapack/lapack/templatelapacktest 
if [ $? -eq 0 ]
then
echo template_lapack test OK
else
echo ERROR in template_lapack test
exit 1
fi

echo
echo Testing template_lapack with threads
"$top_builddir"/source/matrix/template_lapack/lapack/templatelapacktest_threaded
if [ $? -eq 0 ]
then
echo template_lapack with threads test OK
else
echo ERROR in template_lapack with threads test
exit 1
fi

