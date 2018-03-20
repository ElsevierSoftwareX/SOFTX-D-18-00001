#!/bin/sh

currtestdirname=tmpdir_test_ci

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

errorfilename=ergoscf.out.error.ci

echo

echo Testing h2 FCI 6-31G**
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
H     0.0   0.0   0.0
H     0.0   0.0   1.4
EOF
basis = "6-31Gss"
scf.force_unrestricted = 1
do_ci_after_scf = 1
run "HF"
EOINPUT
if 
check_final_ci_corr_energy -0.03387 1e-6 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

echo Testing h2o FCI STO-3G
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.0       0.0       0.0
H    -1.809     0.0       0.0
H     0.453549  1.751221  0.0
EOF
basis = "STO-3G"
scf.force_unrestricted = 1
do_ci_after_scf = 1
run "HF"
EOINPUT
if 
check_final_ci_energy -75.0124258 1e-7 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

echo Testing h2o[+] FCI STO-3G
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.0       0.0       0.0
H    -1.809     0.0       0.0
H     0.453549  1.751221  0.0
EOF
basis = "STO-3G"
scf.force_unrestricted = 1
charge = 1
spin_polarization = 1
do_ci_after_scf = 1
run "HF"
EOINPUT
if
check_final_ci_energy -74.6947713 1e-7 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

echo Testing 2h2 FCI 6-31G
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
H     0.1   0.2   0.3
H     0.3   0.1   1.4
H     2.2   1.1   0.2
H     2.1   1.2   1.5
EOF
basis = "6-31G"
scf.force_unrestricted = 1
do_ci_after_scf = 1
run "HF"
EOINPUT
if
check_final_ci_energy -2.1313212 1e-7 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
#####
if
check_dipole_ci x -0.0307588 1e-5 ;
then
echo Dipole X OK
else
echo ERROR in Dipole X
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
#####
if
check_dipole_ci y -0.016494 1e-5 ;
then
echo Dipole Y OK
else
echo ERROR in Dipole Y
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
#####
if
check_dipole_ci z -0.0020535 1e-5 ;
then
echo Dipole Z OK
else
echo ERROR in Dipole Z
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

echo Testing 2h2 FCI 6-31G with electric field x
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
H     0.1   0.2   0.3
H     0.3   0.1   1.4
H     2.2   1.1   0.2
H     2.1   1.2   1.5
EOF
basis = "6-31G"
scf.force_unrestricted = 1
scf.electric_field_x = -0.03
do_ci_after_scf = 1
run "HF"
EOINPUT
if
check_final_ci_energy -2.1307395 1e-7 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

echo Testing 2h2 FCI 6-31G with electric field y
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
H     0.1   0.2   0.3
H     0.3   0.1   1.4
H     2.2   1.1   0.2
H     2.1   1.2   1.5
EOF
basis = "6-31G"
scf.force_unrestricted = 1
scf.electric_field_y = -0.02
do_ci_after_scf = 1
run "HF"
EOINPUT
if
check_final_ci_energy -2.1310579 1e-7 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

echo Testing 2h2 FCI 6-31G with electric field z
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
H     0.1   0.2   0.3
H     0.3   0.1   1.4
H     2.2   1.1   0.2
H     2.1   1.2   1.5
EOF
basis = "6-31G"
scf.force_unrestricted = 1
scf.electric_field_z = -0.05
do_ci_after_scf = 1
run "HF"
EOINPUT
if
check_final_ci_energy -2.141479 1e-7 ;
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
echo Configuration Interaction [CI] tests completed successfully!
echo
