#!/bin/sh

currtestdirname=tmpdir_test_6_d_funcs

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

errorfilename=ergoscf.out.error.6dfuncs

echo


echo Testing h2o HF/6-31G** with 6 d-functions
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.0       0.0       0.0
H    -1.809     0.0       0.0
H     0.453549  1.751221  0.0
EOF
use_6_d_functions = 1
basis = "6-31Gss"
run "HF"
EOINPUT
if 
check_final_energy -76.0231587 1e-7 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing h2o HF/6-31++G** with 6 d-functions
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.0       0.0       0.0
H    -1.809     0.0       0.0
H     0.453549  1.751221  0.0
EOF
use_6_d_functions = 1
basis = "6-31++Gss"
run "HF"
EOINPUT
if 
check_final_energy -76.0307719 1e-7 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi



echo Testing nh3[+] UHF/6-31G** with 6 d-functions
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
N      0.000000     0.000000     0.000000
H      0.000000     0.000000     1.012316
H      0.969771     0.000000    -0.290392
H     -0.390071     0.887881    -0.290336
EOF
use_6_d_functions = 1
basis = "6-31Gss"
charge = 1
spin_polarization = 1
run "HF"
EOINPUT
if
check_final_energy -55.8547924 1e-6 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if
check_final_S2 0.757025 1e-6 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi



cd ..
rm -r $currtestdirname

echo
echo 6-d-function tests completed successfully!
echo
