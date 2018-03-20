#!/bin/sh

# Run each test in a separate directory, to allow "make check -j" to work properly.
currtestdirname=tmpdir_test_gluala2_b3lyp

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

errorfilename=ergoscf.out.error.gluala2.b3lyp

guess_input_1='
basis = "STO-2G";
use_simple_starting_guess=1;
J_K.threshold_2el_J = 1e-6;
J_K.threshold_2el_K = 1e-6;
J_K.fmm_box_size = 10.0;
scf.convergence_threshold = 1e-2;
run "HF";
'

guess_input_2='
basis = "3-21G";
use_simple_starting_guess=1;
J_K.threshold_2el_J = 1e-6;
J_K.threshold_2el_K = 1e-6;
J_K.fmm_box_size = 10.0;
scf.convergence_threshold = 1e-2;
XC.sparse_mode = 1
XC.radint = 1e-7
XC.angint = 28
run "B3LYP-G";
'

echo

echo Testing GluAla2 B3LYP-G/6-31G**
rm -f ergoscf.out
echo getting starting guess 1...
echo $guess_input_1 | "$top_builddir"/source/ergo -m "$top_srcdir"/mol/GluAla2.xyz > /dev/null
rm -f ergoscf.out
echo getting starting guess 2...
echo $guess_input_2 | "$top_builddir"/source/ergo -m "$top_srcdir"/mol/GluAla2.xyz > /dev/null
rm -f ergoscf.out
echo running B3LYP-G/6-31Gss...
"$top_builddir"/source/ergo -m ../mol/GluAla2.xyz <<EOINPUT > /dev/null
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
run "B3LYP-G"
EOINPUT
if 
check_final_energy -1446.0823389 1e-5 ; 
then
echo OK
grep CONVERGED ergoscf.out
mv ergoscf.out ergoscf.out.gluala2.b3lyp
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi



cd ..
rm -r $currtestdirname

echo
echo GluAla2 b3lyp test completed successfully!
echo
