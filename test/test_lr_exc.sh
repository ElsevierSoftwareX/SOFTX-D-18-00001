#!/bin/sh

currtestdirname=tmpdir_test_lr_exc

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

errorfilename=ergoscf.out.error.lr_exc

echo

echo Testing CO PBE/6-31Gss
rm -f ergoscf.out
./ergo <<EOINPUT 
basis = "6-31Gss"
molecule_inline
C   0.0 0.0 0.0
O   0   0   2.3
EOF
scf.convergence_threshold = 1e-6
J_K.use_fmm = 0
XC.type="GC2"
lr.convergence_threshold = 1e-5
get_excited_state "PBE" 4
EOINPUT
echo  # this is to get an extra newline
if 
check_final_energy -113.172417358797901 1e-5 ;
then
echo Energy OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_lr_eigenvalue 1 0.28181470 4e-5 ;
then
echo Eigenvalue 1 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_lr_eigenvalue 2 0.28181470 4e-5 ;
then
echo Eigenvalue 2 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_lr_eigenvalue 3 0.309447904 4e-5 ;
then
echo Eigenvalue 3 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if [ -n "$KEEP" ]; then
    mv ergoscf.out ergoscf.out_lr_exc_ok
else
 rm ergoscf.out
fi


cd ..
rm -r $currtestdirname

echo
echo Pure TD-DFT tests completed successfully!
echo
