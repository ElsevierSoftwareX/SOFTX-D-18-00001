#!/bin/sh

currtestdirname=tmpdir_test_tdhf_dynamics

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

errorfilename=ergoscf.out.error.tdhfdynamics

echo

echo Testing H2 HF/6-311++Gss followed by TDHF electron dynamics
rm -f ergoscf.out
./ergo <<EOINPUT 
use_simple_starting_guess = 1
scf.convergence_threshold = 1e-5
do_electron_dynamics_after_scf = 1
J_K.threshold_1el = 1e-13
J_K.threshold_2el_J = 1e-12
J_K.threshold_2el_K = 1e-10
basis = "6-311++Gss"
molecule_inline Angstrom
H  0.0   0.0   0.7354
H  0.0   0.0   0.0
EOF
run "HF"
EOINPUT
if 
check_final_energy -1.1325074 1e-5 ; 
then
echo Energy OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

# Check instantaneous dipole values. For values that are
# minima or maxima we use a tighter error tolerance since
# those values should be easier to reproduce, less
# dependent on small shifts in time.

if 
check_tdhf_dipole_at_time      0    -0.00000  1e-2  ;
then
echo Dipole value at time 0 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_tdhf_dipole_at_time      14    0.10866  1e-2  ;
then
echo Dipole value at time 14 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_tdhf_dipole_at_time      21    0.14173  1e-2  ;
then
echo Dipole value at time 21 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_tdhf_dipole_at_time      28    0.08504  2e-2  ;
then
echo Dipole value at time 28 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_tdhf_dipole_at_time      45    -0.34016  1e-2  ;
then
echo Dipole value at time 45 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_tdhf_dipole_at_time      49    -0.37323  2e-2  ;
then
echo Dipole value at time 49 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_tdhf_dipole_at_time      54    -0.34016  3e-2  ;
then
echo Dipole value at time 54 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_tdhf_dipole_at_time      75    0.47717  2e-2  ;
then
echo Dipole value at time 75 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_tdhf_dipole_at_time      79    0.50079  1e-2  ;
then
echo Dipole value at time 79 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_tdhf_dipole_at_time      83    0.46299  2e-2  ;
then
echo Dipole value at time 83 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_tdhf_dipole_at_time      105    -0.44409  3e-2  ;
then
echo Dipole value at time 105 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_tdhf_dipole_at_time      107    -0.47717  1e-2  ;
then
echo Dipole value at time 107 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_tdhf_dipole_at_time      110    -0.45354  1e-2  ;
then
echo Dipole value at time 110 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_tdhf_dipole_at_time      114    -0.48189  1e-2  ;
then
echo Dipole value at time 114 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

if 
check_tdhf_dipole_at_time      116    -0.43465  2e-2  ;
then
echo Dipole value at time 116 OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


# More data if needed:
# (this data is supposed to be (almost) from fig 4a 
# in Li et al Phys. Chem. Chem. Phys., 2005, 7, 233–239:
# "A time-dependent Hartree-Fock approach for studying the electronic optical response of molecules in intense fields"
#
# The values were produced by Elias pointing with the mouse
# and looking at pixel indexes in gimp, so they are only rough values.
#
#     0.00000    -0.00000
#    13.50932     0.10866
#    20.96273     0.14173
#    27.95031     0.08504
#    31.67702     0.00472
#    45.18634    -0.34016
#    49.37888    -0.37323
#    53.57143    -0.34016
#    64.75155     0.10394
#    75.46584     0.47717
#    78.72671     0.50079
#    82.91925     0.46299
#    88.97516     0.26457
#    93.63354     0.00472
#    99.68944    -0.20315
#   104.81366    -0.44409
#   106.67702    -0.47717
#   110.40373    -0.45354
#   114.13043    -0.48189
#   116.45963    -0.43465
#   123.44720    -0.05669
#   130.90062     0.20315
#   134.62733     0.35906
#   136.95652     0.39685
#   141.61491     0.30709
#   145.34161     0.35906
#   148.13665     0.29291
#   153.26087     0.04252
#   155.12422    -0.00000
#   158.85093     0.02835
#   160.71429    -0.00945
#   163.50932    -0.14173
#   166.30435    -0.21260
#   172.82609    -0.04252
#   177.95031    -0.11811
#   187.26708     0.04724
#   194.25466    -0.06142
#   200.77640     0.08504
#   207.29814    -0.08504
#   213.35404     0.08976
#   219.40994    -0.06614
#   225.00000     0.04724



cd ..
rm -r $currtestdirname

echo
echo TDHF-dynamics tests completed successfully!
echo
