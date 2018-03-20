#!/bin/sh

currtestdirname=tmpdir_test_mixedbasis

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

errorfilename=ergoscf.out.error.mixedbasis

echo

echo Testing three h2o with HF and mixed basis 6-31G 6-31G* 6-31G**
rm -f ergoscf.out

cat  <<EOM | ./ergo
molecule_inline Angstrom
O        2.634829230    -0.987901199     1.894282198
H        2.634829230    -0.081293319     1.610447471
H        1.874283601    -1.148764019     2.440354030
O       -2.398712033    -0.928257960    79.638548870
H       -2.398712033    -1.863711511    79.804158818
H       -1.553997179    -0.668352659    79.290113803
O        0.126860479     0.972719769   154.235672931
H        0.126860479     1.052575258   155.182310719
H       -0.728382869     0.680994393   153.942490279
EOF
basis = "6-31Gss"
J_K.threshold_1el = 1e-13
J_K.threshold_2el_J = 1e-11
J_K.threshold_2el_K = 1e-11
use_simple_starting_guess = 0
range 1 = 0 3 "6-31G"
range 2 = 3 3 "6-31Gs"
run "HF"
EOM


if 
check_final_energy -228.016771 1e-6 ; 
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
echo Mixed basis set tests completed successfully!
echo
