#!/bin/sh

currtestdirname=tmpdir_test_lr_exc_hf

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

errorfilename=ergoscf.out.error.lrhf

echo

echo Testing CNOF HF/STO-2G
rm -f ergoscf.out
./ergo <<EOINPUT 
basis = "STO-2G"
molecule_inline
C     0.0       0.0       0.0
N     1.3       0.2       0.4
O     1.5      -0.4       2.0
F    -0.9       2.0       0.5
EOF
J_K.use_fmm = 0
scf.convergence_threshold = 1e-6
lr.convergence_threshold = 1e-5
get_excited_state "HF" 4
EOINPUT
if 
check_final_energy -252.028041454588873 1e-5 ; 
then
echo Energy OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_lr_eigenvalue 1 0.11364134 1e-5 ;
then
echo Eigenvalue 1 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_lr_eigenvalue 2 0.23956153 1e-5 ;
then
echo Eigenvalue 2 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_lr_eigenvalue 3 0.31319248 1e-5 ;
then
echo Eigenvalue 3 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if [ -n "$KEEP" ]; then
    mv ergoscf.out ergoscf.out_lr_exc_hf_ok
else
 rm ergoscf.out
fi



cd ..
rm -r $currtestdirname

echo
echo TD-HF tests completed successfully!
echo
