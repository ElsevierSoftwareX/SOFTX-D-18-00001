#!/bin/sh

currtestdirname=tmpdir_test_ext_elec_field

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

errorfilename=ergoscf.out.error.extelecfield

echo

echo Testing cnof HF/STO-3G, no field
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
C     0.0       0.0       0.0
N     1.3       0.2       0.4
O     1.5      -0.4       2.0
F    -0.9       2.0       0.5
EOF
basis = "STO-3G"
run "HF"
EOINPUT
if 
check_final_energy -260.0484296 1e-6 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

rm ergoscf.out
rm density.bin


echo Testing cnof HF/STO-3G, field in X direction
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
C     0.0       0.0       0.0
N     1.3       0.2       0.4
O     1.5      -0.4       2.0
F    -0.9       2.0       0.5
EOF
basis = "STO-3G"
scf.electric_field_x = 0.02
run "HF"
EOINPUT
if 
check_final_energy -260.0522823 1e-6 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

rm ergoscf.out
rm density.bin


echo Testing cnof HF/STO-3G, field in Y direction
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
C     0.0       0.0       0.0
N     1.3       0.2       0.4
O     1.5      -0.4       2.0
F    -0.9       2.0       0.5
EOF
basis = "STO-3G"
scf.electric_field_y = 0.03
run "HF"
EOINPUT
if 
check_final_energy -260.0445127 1e-6 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

rm ergoscf.out
rm density.bin


echo Testing cnof HF/STO-3G, field in Z direction
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
C     0.0       0.0       0.0
N     1.3       0.2       0.4
O     1.5      -0.4       2.0
F    -0.9       2.0       0.5
EOF
basis = "STO-3G"
scf.electric_field_z = -0.06
run "HF"
EOINPUT
if 
check_final_energy -260.0498294 1e-6 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

rm ergoscf.out
rm density.bin


cd ..
rm -r $currtestdirname

echo
echo External electric field tests completed successfully!
echo
