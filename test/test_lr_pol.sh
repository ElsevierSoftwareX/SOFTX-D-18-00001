#!/bin/sh

currtestdirname=tmpdir_test_lr_pol

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

errorfilename=ergoscf.out.error.lr_pol

# FIXME: in this test script, if check_pol fails somewhere the script anyway continues, does it really properly detect if the computed polarizability values are wrong?

echo
echo Testing NH3 SVWN5/4-31G polarizability
rm -f ergoscf.out
./ergo <<EOINPUT 
basis = "4-31G"
molecule_inline
N      0.00    -0.0   0.0
H      0.00    -1.9   0.1
H      1.64     0.9   0.0
H     -1.64     0.9   0.0
EOF
scf.convergence_threshold = 1e-6
J_K.use_fmm = 0
lr.convergence_threshold = 1e-5
XC.type="LMG"
XC.radint=1e-9
get_polarisability "SVWN5" all 0.0
initial_density="density.bin"
get_polarisability "SVWN5" all 0.2
EOINPUT
echo  # this is to get an extra newline
if 
check_final_energy -55.985633009 1e-5 ;
then
echo Energy OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

check_pol () {
if 
check_polarisability $1 $2 $3 $4 $5;
then
echo Polarisability $1 $2 at frequency $3 OK
else
echo ERROR Polarisability $1 $2 at frequency $3
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
}

check_pol X X   0.0   -7.877230629 1e-5
check_pol Y X   0.0    0.0 1e-5
check_pol Z X   0.0    0.0 1e-5
check_pol X Y   0.0    0.0 1e-5
check_pol Y Y   0.0   -7.936291744 1e-5
check_pol Z Y   0.0   0.2016103778 1e-5
check_pol X Z   0.0   0.0 1e-5
check_pol Y Z   0.0   0.2016104728 1e-5
check_pol Z Z   0.0    -2.45887058 1e-5
check_pol X X   0.2   -8.774076869 1e-5
check_pol Y X   0.2    0.0 1e-5
check_pol Z X   0.2    0.0 1e-5
check_pol X Y   0.2    0.0 1e-5
check_pol Y Y   0.2   -8.880466698  1e-5
check_pol Z Y   0.2    0.2234719254 1e-5
check_pol X Z   0.2    0.0 1e-5
check_pol Y Z   0.2    0.2234719254 1e-5
check_pol Z Z   0.2   -2.849382309 1e-5

rm density.bin potential.bin

echo
echo Testing HF CAM-B3LYP/aug-cc-pVDZ polarizability
rm -f ergoscf.out
./ergo <<EOINPUT 
basis = "aug-cc-pVDZ"
molecule_inline Angstrom
F   0.0 0.0 0.0
H   0   0   0.924719
EOF
scf.convergence_threshold = 1e-6
J_K.use_fmm = 1
lr.convergence_threshold = 1e-5
XC.type="GC2"
XC.radint=1e-9
# FIXME: XC.sparse_mode= not implemented in response calculations.
XC.sparse_mode = 0
get_polarisability "camb3lyp" all 0.0
initial_density="density.bin"
get_polarisability "camb3lyp" all 0.1
EOINPUT
echo  # this is to get an extra newline


if 
check_final_energy -100.43540064008 1e-5 ;
then
echo Energy OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
check_pol X X   0.000000000    -4.458509918 1e-5
check_pol Y X   0.000000000  0.0  1e-5
check_pol Z X   0.000000000  0.0  1e-5
check_pol X Y   0.000000000  0.0  1e-5
check_pol Y Y   0.000000000    -4.458509918 1e-5
check_pol Z Y   0.000000000  0.0  1e-5
check_pol X Z   0.000000000  0.0  1e-5
check_pol Y Z   0.000000000  0.0  1e-5
check_pol Z Z   0.000000000    -6.444241063 4e-5
check_pol X X   0.100000000    -4.583001179 3e-5
check_pol Y X   0.100000000  0.0  1e-5
check_pol Z X   0.100000000  0.0  1e-5
check_pol X Y   0.100000000  0.0  1e-5
check_pol Y Y   0.100000000    -4.583001187 3e-5
check_pol Z Y   0.100000000  0.0  1e-5
check_pol X Z   0.100000000  0.0  1e-5
check_pol Y Z   0.100000000  0.0  1e-5
check_pol Z Z   0.100000000    -6.584672755 1e-5

rm density.bin potential.bin


# Now compute h2o HF/cc-pVDZ polarizability giving values
# corresponding to the isotropic polarizability value -5.011253
# in [PHYSICAL REVIEW E 92, 063301 (2015)].

echo
echo Testing h2o HF/cc-pVDZ polarizability
rm -f ergoscf.out
./ergo <<EOINPUT 
molecule_inline
O     0.0       0.0       0.0
H    -1.809     0.0       0.0
H     0.453549  1.751221  0.0
EOF
basis = "cc-pVDZ"
get_polarisability "HF" all 0.0
EOINPUT
echo  # this is to get an extra newline
if 
check_final_energy -76.0267949 1e-5 ;
then
echo Energy OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
check_pol X X   0.000000000  -6.225251  1e-5
check_pol Y X   0.000000000  -0.882707  1e-5
check_pol Z X   0.000000000   0.0       1e-5
check_pol X Y   0.000000000  -0.882748  1e-5
check_pol Y Y   0.000000000  -5.768206  1e-5
check_pol Z Y   0.000000000   0.0       1e-5
check_pol X Z   0.000000000   0.0       1e-5
check_pol Y Z   0.000000000   0.0       1e-5
check_pol Z Z   0.000000000  -3.040301  1e-5

rm density.bin potential.bin ergoscf.out

cd ..
rm -r $currtestdirname

echo
echo LR polarizability tests completed successfully!
echo
