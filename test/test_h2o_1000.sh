#!/bin/sh

# Run each test in a separate directory, to allow "make check -j" to work properly.
currtestdirname=tmpdir_test_h2o_1000

rm -rf $currtestdirname
mkdir $currtestdirname
cd $currtestdirname

if [ "$top_builddir" = "" ] || [ "$top_builddir" = ".." ] ; then
    top_builddir=../..
fi
if [ "$top_srcdir" = "" ] || [ "$top_srcdir" = ".." ] ; then
    top_srcdir=../..
fi

if test `"$top_builddir"/source/ergo -e precision` = single; then
    echo SKIPPED
    cd .. ; rm -r $currtestdirname
    exit 0
fi

. "$top_srcdir"/test/functions

errorfilename=ergoscf.out.error.h2o_1000

echo

echo Testing h2o_1000 HF/6-31G**
rm -f ergoscf.out
"$top_builddir"/source/ergo -m "$top_srcdir"/mol/h2o_1000.xyz <<EOINPUT > /dev/null
enable_memory_usage_output = 1
basis = "6-31Gss"
use_simple_starting_guess=1
J_K.fmm_box_size = 44
J_K.exchange_box_size = 44
scf.convergence_threshold = 1e-5
run "HF"
EOINPUT
if 
check_final_energy -76022.6431000 1e-4 ; 
then
echo OK
grep "RESC CONVERGED" ergoscf.out
mv ergoscf.out ergoscf.out.h2o_1000
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi



cd ..
rm -r $currtestdirname

echo
echo h2o_1000 Hartree-Fock tests completed successfully!
echo
