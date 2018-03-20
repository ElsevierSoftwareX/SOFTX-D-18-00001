#!/bin/sh          

# Run each test in a separate directory, to allow "make check -j" to work properly.
currtestdirname=tmpdir_test_homo_lumo_h2

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

echo Testing h2 with bond distance 1 au, HF/STO-3G - diagonalization
rm -f ergoscf.out
$script <<EOINPUT
molecule_inline
H   0  0 0
H   0  0 1
EOF
basis = "STO-3G"
use_simple_starting_guess=1              
scf.output_homo_and_lumo_eigenvectors = 1 
scf.use_diagonalization = 1
run "HF"
EOINPUT
if 
check_final_energy -1.0659995 1e-7 ; 
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

echo Testing h2 with bond distance 1 au, HF/STO-3G - non-accelerated purification
rm -f ergoscf.out
$script <<EOINPUT
molecule_inline
H   0  0 0
H   0  0 1
EOF
basis = "STO-3G"
use_simple_starting_guess=1              
scf.purification_with_acceleration = 0  
scf.output_homo_and_lumo_eigenvectors = 1   
scf.min_number_of_iterations = 3
scf.max_no_of_diis_matrices = 1
run "HF"
EOINPUT
if 
check_final_energy -1.0659995 1e-7 ; 
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

rm $filename_eigv_homo
rm $filename_eigv_lumo
rm gabeditfile.gab

echo; 


echo Testing h2 with bond distance 1 au, HF/STO-3G - accelerated purification
rm -f ergoscf.out
$script <<EOINPUT
molecule_inline
H   0  0 0
H   0  0 1
EOF
basis = "STO-3G"
use_simple_starting_guess=1              
scf.purification_with_acceleration = 1  
scf.output_homo_and_lumo_eigenvectors = 1  
scf.min_number_of_iterations = 3
scf.max_no_of_diis_matrices = 1
run "HF"
EOINPUT
if 
check_final_energy -1.0659995 1e-7 ; 
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

rm $filename_eigv_homo
rm $filename_eigv_lumo
rm gabeditfile.gab


echo; 


echo Testing h2 with bond distance 1 au, HF/STO-3G - non-accelerated purification, different initial guess
rm -f ergoscf.out
$script <<EOINPUT
molecule_inline
H   0  0 0
H   0  0 1
EOF
basis = "STO-3G"
use_simple_starting_guess=1              
scf.purification_with_acceleration = 0  
scf.output_homo_and_lumo_eigenvectors = 1   
scf.use_prev_vector_as_initial_guess = 1  
scf.min_number_of_iterations = 4
scf.max_no_of_diis_matrices = 1
run "HF"
EOINPUT
if 
check_final_energy -1.0659995 1e-7 ; 
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

rm $filename_eigv_homo
rm $filename_eigv_lumo
rm gabeditfile.gab

echo; 

# ONLY FOR DEBUGGING
# echo Testing h2 with bond distance 1 au, HF/STO-3G - non-accelerated purification, projection method
# rm -f ergoscf.out
# $script <<EOINPUT
# molecule_inline
# H   0  0 0
# H   0  0 1
# EOF
# basis = "STO-3G"
# use_simple_starting_guess=1              
# scf.purification_with_acceleration = 0  
# scf.output_homo_and_lumo_eigenvectors = 1     
# scf.eigenvectors_method  = "projection"
# scf.use_prev_vector_as_initial_guess = 0  
# scf.eigensolver_accuracy = 1e-7
# scf.min_number_of_iterations = 3
# run "HF"
# EOINPUT
# if 
# check_final_energy -1.0659995 1e-7 ; 
# then
# echo Energy OK
# else
# echo Energy ERROR
# mv ergoscf.out $errorfilename
# echo output file saved as $errorfilename
# exit 1
# fi
# 
# echo "Comparing eigenvectors"
# if 
# $compare_vectors $filename_eigv_homo $filename_eigv_diag_homo $filename_eigv_lumo $filename_eigv_diag_lumo  $TOL  
# then
# echo Eigenvectors OK
# else
# echo Eigenvectors ERROR
# mv ergoscf.out $errorfilename
# echo output file saved as $errorfilename
# exit 1
# fi
# 
# rm $filename_eigv_homo
# rm $filename_eigv_lumo
# rm gabeditfile.gab
# 
# echo; 


# clean everything
rm -f $filename_eigv_diag_homo
rm -f $filename_eigv_diag_lumo
rm -f $filename_eigv_homo
rm -f $filename_eigv_lumo
rm -f gabeditfile.gab



cd ..
rm -r $currtestdirname

echo
echo Eigenvectors test completed successfully!
echo
