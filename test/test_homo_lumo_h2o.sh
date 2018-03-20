#!/bin/sh          

# Run each test in a separate directory, to allow "make check -j" to work properly.
currtestdirname=tmpdir_test_homo_lumo_h2o

if test "$top_builddir" = ""; then
    top_builddir=..
fi
if test "$top_srcdir" = ""; then
    top_srcdir=..
fi

. "$top_srcdir"/test/functions

# Run each test in a separate directory, to allow "make check -j" to work properly.
currdir=`pwd` ; cd $top_builddir ; top_builddir_pwd=`pwd` ; cd $currdir
currdir=`pwd` ; cd $top_srcdir ; top_srcdir_pwd=`pwd` ; cd $currdir
rm -rf $currtestdirname ; mkdir $currtestdirname ; cd $currtestdirname
ln -s "$top_builddir_pwd"/source/ergo ./ergo
script=./ergo

compare_vectors="$top_srcdir_pwd"/test/compare_homo_lumo.sh


if test `./ergo -e precision` = 'single'; then
    echo SKIPPED
    cd .. ; rm -r $currtestdirname
    exit 0
fi

if test `./ergo -e is_cht_used` = 'chunks_and_tasks_is_used'; then
    echo SKIPPED
    cd .. ; rm -r $currtestdirname
    exit 0
fi


errorfilename=ergoscf.out.error.eigenvectors

filename_eigv_diag_homo=lumo_coefficient_vec_diag.txt
filename_eigv_diag_lumo=homo_coefficient_vec_diag.txt
filename_eigv_homo=homo_coefficient_vec.txt
filename_eigv_lumo=lumo_coefficient_vec.txt

TOL=1e-5 # maximum norm of the difference between 
         #eigenvectors computed using the diagonalization and purification

echo Testing h2o HF/cc-pVTZ - diagonalization
rm -f ergoscf.out
$script <<EOINPUT > /dev/null
molecule_inline
O     0.0       0.0       0.0
H    -1.809     0.0       0.0
H     0.453549  1.751221  0.0
EOF
basis = "cc-pVTZ"
scf.output_homo_and_lumo_eigenvectors = 1 
scf.use_diagonalization = 1
run "HF"
EOINPUT
if 
check_final_energy -76.057163 1e-7 ; 
then
echo Energy OK
else
echo Energy ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

mv $filename_eigv_homo  $filename_eigv_diag_homo
mv $filename_eigv_lumo  $filename_eigv_diag_lumo
rm gabeditfile.gab


echo;

echo Testing h2o HF/cc-pVTZ - non-accelerated purification
rm -f ergoscf.out
$script <<EOINPUT > /dev/null
molecule_inline
O     0.0       0.0       0.0
H    -1.809     0.0       0.0
H     0.453549  1.751221  0.0
EOF
basis = "cc-pVTZ"
scf.output_homo_and_lumo_eigenvectors = 1 
scf.purification_with_acceleration = 0  
scf.output_homo_and_lumo_eigenvectors = 1   
run "HF"
EOINPUT
if 
check_final_energy -76.057163 1e-7 ; 
then
echo Energy OK
else
echo Energy ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

echo "Comparing eigenvectors"
if 
$compare_vectors $filename_eigv_homo $filename_eigv_diag_homo $filename_eigv_lumo $filename_eigv_diag_lumo  $TOL  
then
echo Eigenvectors OK
else
echo Eigenvectors ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


# clean everything
rm $filename_eigv_diag_homo
rm $filename_eigv_diag_lumo
rm $filename_eigv_homo
rm $filename_eigv_lumo
rm gabeditfile.gab



cd ..
rm -r $currtestdirname

echo
echo Eigenvectors test completed successfully!
echo
