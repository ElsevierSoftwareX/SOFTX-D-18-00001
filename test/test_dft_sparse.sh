#!/bin/sh

currtestdirname=tmpdir_test_dft_sparse

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

errorfilename=ergoscf.out.error.dftsparse

echo

echo Testing c2h8 LDA/6-31G*
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
C     0.0       0.0       0.0
C     0.0       0.0      10.0
H      0.000000     0.000000     1.084800
H      1.022759     0.000000    -0.361600
H     -0.511380     0.885735    -0.361600
H     -0.511380    -0.885735    -0.361600
H      0.000000     0.000000    11.084800
H      1.022759     0.000000     9.638400
H     -0.511380     0.885735     9.638400
H     -0.511380    -0.885735     9.638400
EOF
basis = "4-31G"
XC.sparse_mode = 1
run "LDA"
EOINPUT
if 
check_final_energy -80.097373423 1e-5 ; 
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
XC.sparse_mode = 1
XC.radint=1e-13
XC.type="LMG"
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

cd ..
rm -r $currtestdirname

echo
echo Sparse DFT tests completed successfully!
echo
