#!/bin/sh

currtestdirname=tmpdir_test_dft_pure

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

errorfilename=ergoscf.out.error.dftpure

echo

echo Testing cnof BLYP/STO-2G
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
C     0.0       0.0       0.0
N     1.3       0.2       0.4
O     1.5      -0.4       2.0
F    -0.9       2.0       0.5
EOF
basis = "STO-2G"
XC.sparse_mode=1
scf.convergence_threshold = 1e-6
XC.radint=1e-10
run "BLYP"
EOINPUT
if 
check_final_energy -253.0730954 2e-5 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing h2o BLYP/6-31G**
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.0       0.0       0.0
H    -1.809     0.0       0.0
H     0.453549  1.751221  0.0
EOF
basis = "6-31Gss"
XC.type="LMG"
XC.radint = 1e-13
run "BLYP"
EOINPUT
if 
check_final_energy -76.396247 1e-5 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing cnof BLYP/6-31G
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
C     0.0       0.0       0.0
N     1.3       0.2       0.4
O     1.5      -0.4       2.0
F    -0.9       2.0       0.5
EOF
basis = "6-31G"
XC.sparse_mode=1
scf.convergence_threshold = 1e-6
XC.radint=1e-11
run "BLYP"
EOINPUT
if 
check_final_energy -265.352229094 1e-5 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing nh3[+] UBLYP/6-31G**
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
N      0.000000     0.000000     0.000000
H      0.000000     0.000000     1.012316
H      0.969771     0.000000    -0.290392
H     -0.390071     0.887881    -0.290336
EOF
charge = 1
spin_polarization = 1
basis = "6-31Gss"
run "BLYP"
EOINPUT
if
check_final_energy -56.1382272 1e-5 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if
check_final_S2 0.752002 1e-6 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing nh3[+] UBLYP/6-31G** with sparse dft turned off
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
N      0.000000     0.000000     0.000000
H      0.000000     0.000000     1.012316
H      0.969771     0.000000    -0.290392
H     -0.390071     0.887881    -0.290336
EOF
charge = 1
spin_polarization = 1
XC.sparse_mode = 0
basis = "6-31Gss"
run "BLYP"
EOINPUT
if
check_final_energy -56.1382272 1e-5
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if
check_final_S2 0.752002 1e-6
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
echo Pure DFT tests completed successfully!
echo
