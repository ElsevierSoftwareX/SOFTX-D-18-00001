#!/bin/sh

currtestdirname=tmpdir_test_gradient

if test "$top_builddir" = ""; then
    top_builddir=..
fi
if test "$top_srcdir" = ""; then
    top_srcdir=..
fi

. "$top_srcdir"/test/functions

# Run each test in a separate directory, to allow "make check -j" to work properly.
currdir=`pwd` ; cd $top_builddir ; top_builddir_pwd=`pwd` ; cd $currdir
rm -rf $currtestdirname ; mkdir $currtestdirname ; cd $currtestdirname
ln -s "$top_builddir_pwd"/source/ergo ./ergo

if test `./ergo -e precision` = 'single'; then
    echo SKIPPED
    cd .. ; rm -r $currtestdirname
    exit 0
fi

errorfilename1=ergoscf.out.error.gradient.1
errorfilename2=ergoscf.out.error.gradient.2
errorfilename3=ergoscf.out.error.gradient.3

epsilon=0.0001
tolerance=0.00001

echo

echo Testing gradient computation for h2o HF/6-31G** using finite differences and fixed basis set, achieved via ghost basis.


coord1=`print_sum_of_values 0.3  $epsilon`
coord2=`print_sum_of_values 0.3 -$epsilon`
echo coord1 is $coord1
echo coord2 is $coord2

rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
ghost_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
basis = "none"
ghost_basis = "6-31Gss"
scf.compute_gradient_fixeddens = 1
run "HF"
EOINPUT
mv ergoscf.out ergoscf.out.gradient

./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.2     $coord1  0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
ghost_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
basis = "none"
ghost_basis = "6-31Gss"
scf.compute_gradient_fixeddens = 1
run "HF"
EOINPUT
mv ergoscf.out ergoscf.out.sav1

./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.2     $coord2  0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
ghost_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
basis = "none"
ghost_basis = "6-31Gss"
run "HF"
EOINPUT
mv ergoscf.out ergoscf.out.sav2

if 
check_gradient_component_by_finite_diff 0 1 $epsilon $tolerance
then
echo OK
rm ergoscf.out.sav1 ergoscf.out.sav2 ergoscf.out.gradient
else
echo ERROR
mv ergoscf.out.sav1 $errorfilename1
mv ergoscf.out.sav2 $errorfilename2
mv ergoscf.out.gradient $errorfilename3
echo output files saved as $errorfilename1 $errorfilename2 $errorfilename3
exit 1
fi


coord1=`print_sum_of_values -1.8  $epsilon`
coord2=`print_sum_of_values -1.8 -$epsilon`
echo coord1 is $coord1
echo coord2 is $coord2

rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
ghost_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
basis = "none"
ghost_basis = "6-31Gss"
scf.compute_gradient_fixeddens = 1
run "HF"
EOINPUT
mv ergoscf.out ergoscf.out.gradient

./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.2     0.3      0.4
H     $coord1 0.1      0.3
H     0.4     1.7      0.2
EOF
ghost_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
basis = "none"
ghost_basis = "6-31Gss"
scf.compute_gradient_fixeddens = 1
run "HF"
EOINPUT
mv ergoscf.out ergoscf.out.sav1

./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.2     0.3      0.4
H     $coord2 0.1      0.3
H     0.4     1.7      0.2
EOF
ghost_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
basis = "none"
ghost_basis = "6-31Gss"
run "HF"
EOINPUT
mv ergoscf.out ergoscf.out.sav2

if 
check_gradient_component_by_finite_diff 1 0 $epsilon $tolerance
then
echo OK
rm ergoscf.out.sav1 ergoscf.out.sav2 ergoscf.out.gradient
else
echo ERROR
mv ergoscf.out.sav1 $errorfilename1
mv ergoscf.out.sav2 $errorfilename2
mv ergoscf.out.gradient $errorfilename3
echo output files saved as $errorfilename1 $errorfilename2 $errorfilename3
exit 1
fi


coord1=`print_sum_of_values 0.2  $epsilon`
coord2=`print_sum_of_values 0.2 -$epsilon`
echo coord1 is $coord1
echo coord2 is $coord2

rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
ghost_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
basis = "none"
ghost_basis = "6-31Gss"
scf.compute_gradient_fixeddens = 1
run "HF"
EOINPUT
mv ergoscf.out ergoscf.out.gradient

./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      $coord1
EOF
ghost_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
basis = "none"
ghost_basis = "6-31Gss"
scf.compute_gradient_fixeddens = 1
run "HF"
EOINPUT
mv ergoscf.out ergoscf.out.sav1

./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      $coord2
EOF
ghost_inline
O     0.2     0.3      0.4
H    -1.8     0.1      0.3
H     0.4     1.7      0.2
EOF
basis = "none"
ghost_basis = "6-31Gss"
run "HF"
EOINPUT
mv ergoscf.out ergoscf.out.sav2

if 
check_gradient_component_by_finite_diff 2 2 $epsilon $tolerance
then
echo OK
rm ergoscf.out.sav1 ergoscf.out.sav2 ergoscf.out.gradient
else
echo ERROR
mv ergoscf.out.sav1 $errorfilename1
mv ergoscf.out.sav2 $errorfilename2
mv ergoscf.out.gradient $errorfilename3
echo output files saved as $errorfilename1 $errorfilename2 $errorfilename3
exit 1
fi


rm -f ergoscf.out
rm -f density.bin

cd ..
rm -r $currtestdirname

echo
echo Gradient tests completed successfully!
echo
