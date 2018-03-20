#!/bin/sh

currtestdirname=tmpdir_test_dft_hybrid

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

errorfilename=ergoscf.out.error.dfthybrid

echo

echo Testing h2o B3LYP/6-31G**   using g03-style B3LYP
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.0       0.0       0.0
H    -1.809     0.0       0.0
H     0.453549  1.751221  0.0
EOF
basis = "6-31Gss"
XC.type="LMG"
XC.radint=1e-13
run "B3LYP-G"
EOINPUT
if 
check_final_energy -76.4180487 1e-5 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing cnof B3LYP/6-31G   using g03-style B3LYP
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
C     0.0       0.0       0.0
N     1.3       0.2       0.4
O     1.5      -0.4       2.0
F    -0.9       2.0       0.5
EOF
basis = "6-31G"
XC.type="Turbo"
XC.radint=1e-10
scf.convergence_threshold = 1e-6
run "B3LYP-G"
EOINPUT
if 
check_final_energy -265.391329 1e-5 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing sih4 BHandHLYP/6-31G**
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
Si    0.0000    0.0000    0.0000
H     0.9385    0.9654    0.6265
H     0.7506   -1.2008   -0.4472
H    -1.0315   -0.4013    0.9900
H    -0.6575    0.6367   -1.1694
EOF
basis = "6-31Gss"
run "BHandHLYP"
EOINPUT
if
check_final_energy -291.857715 5e-5 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing sih3[-] UBHandHLYP/6-31G**
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
Si    0.0000    0.0000    0.0000
H     0.9385    0.9654    0.6265
H     0.7506   -1.2008   -0.4472
H    -1.0315   -0.4013    0.9900
EOF
charge = -1
basis = "6-31Gss"
scf.force_unrestricted = 1
run "BHandHLYP"
EOINPUT
if
check_final_energy -291.223592 5e-5 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing sih3 UBHandHLYP/6-31G**
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
Si    0.0000    0.0000    0.0000
H     0.9385    0.9654    0.6265
H     0.7506   -1.2008   -0.4472
H    -1.0315   -0.4013    0.9900
EOF
spin_polarization = 1
basis = "6-31Gss"
run "BHandHLYP"
EOINPUT
if
check_final_energy -291.2083758 5e-5 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing nh3[+] UB3LYP/6-31G**   using g03-style B3LYP
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
run "B3LYP-G"
EOINPUT
if
check_final_energy -56.163575 1e-5 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if
check_final_S2 0.752352 1e-6 ;
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
echo Hybrid DFT tests completed successfully!
echo
