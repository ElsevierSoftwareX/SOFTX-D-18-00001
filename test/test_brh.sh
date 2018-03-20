#!/bin/sh

currtestdirname=tmpdir_test_brh

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

errorfilename=ergoscf.out.error.brh

echo

echo Testing BrH, HF/STO-3G
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
Br   0.0000    0.0000    0.0000
H    0.6731   -0.2281   -1.2649
EOF
basis = "STO-3G"
run "HF"
EOINPUT
if 
check_final_energy -2545.2272426 1e-6 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

echo Testing BrH, HF/STO-3G, with use_6_d_functions
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
Br   0.0000    0.0000    0.0000
H    0.6731   -0.2281   -1.2649
EOF
basis = "STO-3G"
use_6_d_functions = 1
run "HF"
EOINPUT
if 
check_final_energy -2546.0959772 1e-6 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

echo Testing BrH, HF/3-21G
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
Br   0.0000    0.0000    0.0000
H    0.6731   -0.2281   -1.2649
EOF
basis = "3-21G"
run "HF"
EOINPUT
if 
check_final_energy -2560.6153595 1e-5 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

echo Testing BrH, HF/3-21G, with use_6_d_functions
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
Br   0.0000    0.0000    0.0000
H    0.6731   -0.2281   -1.2649
EOF
basis = "3-21G"
use_6_d_functions = 1
run "HF"
EOINPUT
if 
check_final_energy -2560.6206445 1e-5 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

echo Testing BrH, HF/6-311G**
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
Br   0.0000    0.0000    0.0000
H    0.6731   -0.2281   -1.2649
EOF
basis = "6-311Gss"
run "HF"
EOINPUT
if 
check_final_energy -2572.9546129 1e-7 ;
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
echo BrH \(Hydrogen bromide\) tests completed successfully!
echo
