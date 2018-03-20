#!/bin/sh

currtestdirname=tmpdir_test_exit_status

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

errorfilename=ergoscf.out.error.exitstatus

echo

# Here we put an if-statement around the call to the ergo executable,
# to check that the exit status is OK. In this way we detect if there
# is some problem when the program is funishing, like a seg fault
# error in a destructor. In such cases the output files may seem OK
# but we still want to detect that something went wrong.

echo Testing h2o HF 6-31G**
rm -f ergoscf.out

if
./ergo <<EOINPUT > /dev/null
molecule_inline
O     0.0       0.0       0.0
H    -1.809     0.0       0.0
H     0.453549  1.751221  0.0
EOF
basis = "6-31Gss"
scf.convergence_threshold = 2e-3
scf.purification_subspace_err_limit = 1e-2
scf.use_diag_on_error = 0
run "HF"
EOINPUT
then
    echo Executed ergo, exit status OK.
else
    echo ERROR: executed ergo, exit status indicated error.
    exit 1
fi

if
check_final_energy -76.0226431 1e-4 ;
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


echo Testing Be[-] UHF 6-31G*
rm -f ergoscf.out

if
./ergo <<EOINPUT > /dev/null
molecule_inline
Be   0    0   0
EOF
basis = "6-31Gs"
charge = -1
scf.convergence_threshold = 1e-4
scf.purification_subspace_err_limit = 1e-2
scf.use_diag_on_error = 1
spin_polarization = 1
run "HF"
EOINPUT
then
    echo Executed ergo, exit status OK.
else
    echo ERROR: executed ergo, exit status indicated error.
    exit 1
fi

if
check_final_energy -14.4939039 1e-4 ;
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
echo Exit status tests completed successfully!
echo
