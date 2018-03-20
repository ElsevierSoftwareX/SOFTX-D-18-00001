#!/bin/sh

currtestdirname=tmpdir_test_dft_hicu

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

errorfilename=ergoscf.out.error.dfthicu

echo



echo Testing twisted h2o BLYP/6-31G** with HiCu grid
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.3       0.11      0.23
H    -1.89      0.18     -0.07
H     0.55      1.66      0.27
EOF
basis = "6-31Gss"
XC.type="HICU"
scf.convergence_threshold = 1e-6
run "BLYP"
EOINPUT
if 
check_final_energy -76.3474656 8e-5 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing H4 BHandHLYP/cc-pVDZ with HiCu grid
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
H   2.0   3.0   12.02
H   2.7   3.4   12.16
H   2.2   3.7   11.85
H   2.1   3.1   12.74
EOF
basis = "cc-pVDZ"
XC.sparse_mode=1
XC.type="HICU"
scf.convergence_threshold = 1e-6
run "BHandHLYP"
EOINPUT
if 
check_final_energy -0.2614294 5e-5 ; 
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing H2 BPW91/cc-pVTZ with HiCu grid
rm -f ergoscf.out
./ergo <<EOINPUT > /dev/null
molecule_inline
H   2.0   3.0   11.12
H   2.5   3.2   11.31
EOF
basis = "cc-pVTZ"
scf.min_number_of_iterations = 2
scf.convergence_threshold = 1e-6
XC.type="HICU"
run "BPW91"
EOINPUT
if
check_final_energy -0.7042593 2e-5 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


echo Testing c2h8 LDA/4-31G with sparse XC and HiCu grid
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
XC.type="HICU"
run "LDA"
EOINPUT
if 
check_final_energy -80.097373423 3e-4 ; 
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
echo DFT HiCu grid tests completed successfully!
echo
