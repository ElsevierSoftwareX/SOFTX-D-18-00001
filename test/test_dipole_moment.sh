#!/bin/sh

currtestdirname=tmpdir_test_dipole_moment

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

errorfilename=ergoscf.out.error.dipolemoment

echo



echo Testing cnof BHANDHLYP/6-31G
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
C     0.0       0.0       0.0
N     1.3       0.2       0.4
O     1.5      -0.4       2.0
F    -0.9       2.0       0.5
EOF
basis = "6-31G"
XC.radint=1e-12
XC.angint = 35
scf.output_mulliken_pop = 1
run "BHANDHLYP"
EOINPUT
if 
check_final_energy -265.2667394 1e-5 ; 
then
echo Energy OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
#####
if 
check_dipole x 0.0101015 1e-4 ; 
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
check_dipole y -0.0169424 1e-4 ; 
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
check_dipole z -0.5035157 1e-4 ; 
then
echo Dipole Z OK
else
echo ERROR in Dipole Z
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
#####
if 
check_mulliken_charge 0 -0.068157 1e-4 ; 
then
echo Mulliken charge OK
else
echo ERROR in Mulliken
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
#####
if 
check_mulliken_charge 1 0.579839 1e-4 ; 
then
echo Mulliken charge OK
else
echo ERROR in Mulliken charge
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
#####
if 
check_mulliken_charge 2 -0.321502 1e-4 ; 
then
echo Mulliken charge OK
else
echo ERROR in Mulliken charge
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
#####
if 
check_mulliken_charge 3 -0.190181 1e-4 ; 
then
echo Mulliken charge OK
else
echo ERROR in Mulliken charge
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi



echo Testing nh3-triplet BHANDHLYP/6-31G**
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
N      0.000000     0.000000     0.000000
H      0.000000     0.000000     1.012316
H      0.969771     0.000000    -0.290392
H     -0.390071     0.887881    -0.290336
EOF
charge = 0
spin_polarization = 2
basis = "6-31Gss"
scf.output_mulliken_pop = 1
run "BHANDHLYP"
EOINPUT
if
check_final_energy -56.2444854 1e-5 ;
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
check_final_S2 2.001894 1e-6 ;
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
check_dipole x -0.0692112 1e-4 ; 
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
check_dipole y -0.1060111 1e-4 ; 
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
check_dipole z -0.0515257 1e-4 ; 
then
echo Dipole Z OK
else
echo ERROR in Dipole Z
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
#####
if 
check_mulliken_charge 0 0.087443 1e-4 ; 
then
echo Mulliken charge OK
else
echo ERROR in Mulliken
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
#####
if 
check_mulliken_charge 1 -0.029146 1e-4 ; 
then
echo Mulliken charge OK
else
echo ERROR in Mulliken charge
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
#####
if 
check_mulliken_charge 2 -0.029148 1e-4 ; 
then
echo Mulliken charge OK
else
echo ERROR in Mulliken charge
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
#####
if 
check_mulliken_charge 3 -0.029149 1e-4 ; 
then
echo Mulliken charge OK
else
echo ERROR in Mulliken charge
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
#####
if 
check_mulliken_spin 0 0.744599 1e-4 ; 
then
echo Mulliken spin density OK
else
echo ERROR in Mulliken spin density
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
#####
if 
check_mulliken_spin 1 0.418462 1e-4 ; 
then
echo Mulliken spin density OK
else
echo ERROR in Mulliken spin density
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
#####
if 
check_mulliken_spin 2 0.418473 1e-4 ; 
then
echo Mulliken spin density OK
else
echo ERROR in Mulliken spin density
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
#####
if 
check_mulliken_spin 3 0.418465 1e-4 ; 
then
echo Mulliken spin density OK
else
echo ERROR in Mulliken spin density
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


cd ..
rm -r $currtestdirname

echo
echo Dipole moment and Mulliken pop tests completed successfully!
echo
