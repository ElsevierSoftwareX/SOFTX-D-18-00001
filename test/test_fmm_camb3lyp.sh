#!/bin/sh

currtestdirname=tmpdir_test_fmm_camb3lyp

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

errorfilename=ergoscf.out.error.fmmcamb3lyp

echo

echo Testing four_h2o CAMB3LYP/6-31G using FMM
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline Angstrom
H                  8.675000    0.000000  -16.640000
H                  9.579936    0.000000  -17.920455
H                  7.203213    0.037757  -17.673336
H                  8.108149    0.037757  -18.953791
H                  9.241330    0.571501   -2.493601
H                 10.146266    0.571501   -3.774055
H                  7.134485    0.659897   -2.982496
H                  8.039421    0.659897   -4.262951
O                  8.675000    0.000000  -17.600000
O                  7.203213    0.037757  -18.633336
O                  9.241330    0.571501   -3.453601
O                  7.134485    0.659897   -3.942496
EOF
basis = "6-31G"
XC.sparse_mode=1
XC.radint=1e-11
XC.angint = 34
run "CAMB3LYP"
EOINPUT
if 
check_final_energy -305.2790607 1e-5 ; 
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


echo Testing cnof CAMB3LYP/6-31G using FMM
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
XC.radint=1e-10
XC.angint = 34
run "CAMB3LYP"
EOINPUT
if 
check_final_energy -265.31129062 1e-5 ; 
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


echo Testing h2o CAMB3LYP/cc-pVDZ using FMM
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.0       0.0       0.0
H    -1.809     0.0       0.0
H     0.453549  1.751221  0.0
EOF
basis = "cc-pVDZ"
XC.sparse_mode=1
XC.radint=1e-10
XC.angint = 34
run "CAMB3LYP"
EOINPUT
if 
check_final_energy -76.391777086 1e-5 ; 
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
echo FMM CAMB3LYP tests completed successfully!
echo
