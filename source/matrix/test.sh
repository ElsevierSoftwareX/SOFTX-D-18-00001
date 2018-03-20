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


echo Testing code sensitive for some compilers when OpenMP is used  
"$top_builddir"/source/matrix/omptest
if [ $? -eq 0 ]
then
echo OpenMP test OK
else
echo ERROR in OpenMP test
exit 1
fi

echo
echo Testing matrix library
"$top_builddir"/source/matrix/mattest "$top_srcdir"/source/matrix
if [ $? -eq 0 ]
then
echo Matrix library tests OK
else
echo ERROR in matrix library test
exit 1
fi

echo
echo Testing matrix library LanczosSeveralLargestEig functionality
"$top_builddir"/source/matrix/mattest_lan
if [ $? -eq 0 ]
then
echo Matrix library LanczosSeveralLargestEig test OK
else
echo ERROR in matrix library LanczosSeveralLargestEig test
exit 1
fi

if [ "$RUN_BENCHMARK" = "1" ]
    then
    echo
    echo Benchmark of matrix library:
    "$top_builddir"/source/matrix/matbench 1000 
    if [ $? -eq 0 ]
	then
	echo Benchmark returned OK
    else
	echo ERROR in matrix library benchmark
	exit 1
    fi
    echo Running BLAS benchmark, result in file blastime.m
    "$top_builddir"/source/matrix/blastime 100 blastime.m 
    if [ $? -eq 0 ]
	then
	echo BLAS benchmark returned OK
    else
	echo ERROR in BLAS benchmark
	exit 1
    fi
else
    echo Skipping matrix library benchmark
    echo To run benchmark, run check as: make check RUN_BENCHMARK=1
fi
