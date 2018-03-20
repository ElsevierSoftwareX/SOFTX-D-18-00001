#!/bin/bash


# Prefer gawk - we know exactly what it can do.
# awk on Sun does not support functions, need to use nawk for this
if gawk '{print 1}'</dev/null > /dev/null 2>&1; then
   AWK=gawk
elif nawk '{print 1}'</dev/null > /dev/null 2>&1; then
   AWK=nawk
else
   AWK=awk
fi

INPUT1=$1
INPUT2=$2
INPUT3=$3
INPUT4=$4
TOL=$5


# check if input files exist (may happen in case eigenvectors are not computed)
if [ ! -e "$INPUT1" ]; then echo "$INPUT1 file does not exist."; exit 1; fi
if [ ! -e "$INPUT2" ]; then echo "$INPUT2 file does not exist."; exit 1; fi
if [ ! -e "$INPUT3" ]; then echo "$INPUT3 file does not exist."; exit 1; fi
if [ ! -e "$INPUT4" ]; then echo "$INPUT4 file does not exist."; exit 1; fi


# we get directory where compare.sh resides, awk_script.awk should be in the same directory
# it is needed of compare.sh called from another directory
BASEDIR=$(dirname "$0") 


$AWK -v tol=$TOL -f $BASEDIR/awk_script.awk $INPUT1 $INPUT2 

exit_status=$?  # exit status from the last command
if [ "$exit_status" -eq 1 ]; then
    echo "HOMO eigenvectors are not equal!"
    exit 1
else
    echo "HOMO OK!"
    # no exit here
fi


$AWK -v tol=$TOL -f $BASEDIR/awk_script.awk $INPUT3 $INPUT4 

exit_status=$?  # exit status from the last command
if [ "$exit_status" -eq 1 ]; then
    echo "LUMO eigenvectors are not equal!"
    exit 1
else
    echo "LUMO OK!"
    exit 0
fi









