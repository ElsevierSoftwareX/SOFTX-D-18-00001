#!/bin/sh

# Run each test in a separate directory, to allow "make check -j" to work properly.
currtestdirname=tmpdir_test_h2o_216_hf

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

errorfilename=ergoscf.out.error.h2o_216

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

echo Testing h2o_216 HF/6-31G**
rm -f ergoscf.out
echo getting starting guess...
echo $guess_input | "$top_builddir"/source/ergo -m "$top_srcdir"/mol/h2o_216.xyz > /dev/null
rm -f ergoscf.out
echo running HF/6-31Gss...
"$top_builddir"/source/ergo -m "$top_srcdir"/mol/h2o_216.xyz <<EOINPUT > /dev/null
enable_memory_usage_output = 1
basis = "6-31Gss"
initial_density = "density.bin"
J_K.threshold_2el_J = 1e-10
J_K.threshold_2el_K = 1e-10
J_K.fmm_box_size = 8.8
scf.convergence_threshold = 1e-5
scf.max_no_of_diis_matrices = 4
run "HF"
EOINPUT
if 
check_final_energy -16420.8406978 2e-6 ; 
then
echo OK
grep "RESC CONVERGED" ergoscf.out
mv ergoscf.out ergoscf.out.h2o_216
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi



cd ..
rm -r $currtestdirname

echo
echo h2o_216 Hartree-Fock test completed successfully!
echo
