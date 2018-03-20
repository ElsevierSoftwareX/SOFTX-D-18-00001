#!/bin/sh

# Run each test in a separate directory, to allow "make check -j" to work properly.
currtestdirname=tmpdir_test_h2o_216_blyp

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
    exit 0
fi

. "$top_srcdir"/test/functions

errorfilename=ergoscf.out.error.h2o_216_blyp

guess_input='
basis = "STO-2G";
use_simple_starting_guess=1;
J_K.threshold_2el_J = 1e-6;
J_K.threshold_2el_K = 1e-6;
J_K.fmm_box_size = 10.0;
scf.convergence_threshold = 1e-2;
run "HF";
'

echo

echo Testing h2o_216 BLYP-G/6-31G**
rm -f ergoscf.out
echo getting starting guess...
echo $guess_input | "$top_builddir"/source/ergo -m "$top_srcdir"/mol/h2o_216.xyz > /dev/null
rm -f ergoscf.out
echo running BLYP-G/6-31Gss...
"$top_builddir"/source/ergo -m "$top_srcdir"/mol/h2o_216.xyz <<EOINPUT > /dev/null
enable_memory_usage_output = 1
basis = "6-31Gss"
initial_density = "density.bin"
J_K.threshold_2el_J = 1e-10
J_K.threshold_2el_K = 1e-10
J_K.fmm_box_size = 8.8
scf.convergence_threshold = 1e-5
scf.max_no_of_diis_matrices = 4
XC.sparse_mode = 1
XC.radint = 1e-15
XC.angint = 34
run "BLYP"
EOINPUT
if 
check_final_energy -16501.3079147 1e-4 ; 
then
echo OK
grep "RESC CONVERGED" ergoscf.out
mv ergoscf.out ergoscf.out.h2o_216_blyp
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi



cd ..
rm -r $currtestdirname

echo
echo h2o_216 blyp test completed successfully!
echo
