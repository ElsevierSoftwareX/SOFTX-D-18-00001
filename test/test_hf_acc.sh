#!/bin/sh

currtestdirname=tmpdir_test_hf_acc

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

errorfilename=ergoscf.out.error.hf_acc

echo

echo Testing h2 with bond distance 1 au, HF/STO-3G
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
H   0  0 0
H   0  0 1
EOF
basis = "STO-3G"
use_simple_starting_guess=1
scf.purification_with_acceleration = 1
run "HF"
EOINPUT
if 
check_final_energy -1.0659995 1e-7 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing h2 with bond distance 5 au, HF/6-31G**
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
H   -2.5  0  0
H    2.5  0  0
EOF
basis = "6-31Gss"
scf.purification_with_acceleration = 1
run "HF"
EOINPUT
if
check_final_energy -0.8433867 1e-7 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing h2o HF/6-31G**
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.0       0.0       0.0
H    -1.809     0.0       0.0
H     0.453549  1.751221  0.0
EOF
basis = "6-31Gss"
scf.purification_with_acceleration = 1
run "HF"
EOINPUT
if 
check_final_energy -76.0226431 1e-7 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing h2o HF/6-31++G**
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.0       0.0       0.0
H    -1.809     0.0       0.0
H     0.453549  1.751221  0.0
EOF
scf.purification_with_acceleration = 1
basis = "6-31++Gss"
run "HF"
EOINPUT
if 
check_final_energy -76.0304891 1e-7 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing h2o HF/cc-pVTZ
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.0       0.0       0.0
H    -1.809     0.0       0.0
H     0.453549  1.751221  0.0
EOF
scf.purification_with_acceleration = 1
basis = "cc-pVTZ"
run "HF"
EOINPUT
if 
check_final_energy -76.057163 1e-7 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing cnof HF/6-31G
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
C     0.0       0.0       0.0
N     1.3       0.2       0.4
O     1.5      -0.4       2.0
F    -0.9       2.0       0.5
EOF
scf.purification_with_acceleration = 1
basis = "6-31G"
run "HF"
EOINPUT
if 
check_final_energy -264.0542682 1e-6 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing nh3[+] UHF/6-31G**
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
N      0.000000     0.000000     0.000000
H      0.000000     0.000000     1.012316
H      0.969771     0.000000    -0.290392
H     -0.390071     0.887881    -0.290336
EOF
basis = "6-31Gss"
charge = 1
spin_polarization = 1
scf.purification_with_acceleration = 1
run "HF"
EOINPUT
if
check_final_energy -55.8546235 1e-6 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if
check_final_S2 0.757013 1e-6 ;
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
echo Hartree-Fock tests using accelerated recursive expansion completed successfully!
echo
