#!/bin/sh

# Run each test in a separate directory, to allow "make check -j" to work properly.
currtestdirname=tmpdir_test_gluala2_hf

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

errorfilename=ergoscf.out.error.gluala2.hf

guess_input='
basis = "STO-2G";
use_simple_starting_guess=1;
J_K.threshold_2el_J = 1e-6;
J_K.threshold_2el_K = 1e-6;
J_K.fmm_box_size = 4.4;
scf.convergence_threshold = 1e-2;
run "HF";
'

echo

echo Testing GluAla2 HF/6-31G**
rm -f ergoscf.out
echo getting starting guess...
echo $guess_input | "$top_builddir"/source/ergo -m "$top_srcdir"/mol/GluAla2.xyz > /dev/null
rm -f ergoscf.out
echo running HF/6-31Gss...
"$top_builddir"/source/ergo -m ../mol/GluAla2.xyz <<EOINPUT > /dev/null
enable_memory_usage_output = 1
basis = "6-31Gss"
initial_density = "density.bin"
J_K.threshold_2el_J = 1e-10
J_K.threshold_2el_K = 1e-10
J_K.fmm_box_size = 3.3
scf.convergence_threshold = 1e-5
run "HF"
EOINPUT
if 
check_final_energy -1437.6933948 2e-7 ; 
then
echo OK
grep CONVERGED ergoscf.out
mv ergoscf.out ergoscf.out.gluala2.hf
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi



cd ..
rm -r $currtestdirname

echo
echo GluAla2 Hartree-Fock test completed successfully!
echo
