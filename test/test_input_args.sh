#!/bin/sh

currtestdirname=tmpdir_test_input_args

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

errorfilename=ergoscf.out.error.inputargs

echo


echo Testing nh3[+] UHF/6-31G** with -e arguments and regular input file
rm -f ergoscf.out
echo molecule_inline Angstrom                   > tmp_ergo_arg_test_input_file.ego
echo N      0.000000     0.000000     0.000000 >> tmp_ergo_arg_test_input_file.ego
echo H      0.000000     0.000000     1.012316 >> tmp_ergo_arg_test_input_file.ego
echo H      0.969771     0.000000    -0.290392 >> tmp_ergo_arg_test_input_file.ego
echo H     -0.390071     0.887881    -0.290336 >> tmp_ergo_arg_test_input_file.ego
echo EOF                                       >> tmp_ergo_arg_test_input_file.ego
echo basis = \"6-31Gss\"                       >> tmp_ergo_arg_test_input_file.ego
echo run \"HF\"                                >> tmp_ergo_arg_test_input_file.ego
./ergo -e "charge = 1" -e "spin_polarization = 1" tmp_ergo_arg_test_input_file.ego > /dev/null
if
check_final_energy -55.8546235 1e-6 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
if
check_final_S2 0.757013 1e-6 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

# Do the same test again, but this time with the molecule specified using a -m argument

echo Testing nh3[+] UHF/6-31G** with -e args, -m molecule arg, and regular input file
rm -f ergoscf.out
echo 4                                          > tmp_ergo_arg_test_mol.xyz
echo                                           >> tmp_ergo_arg_test_mol.xyz
echo N      0.000000     0.000000     0.000000 >> tmp_ergo_arg_test_mol.xyz
echo H      0.000000     0.000000     1.012316 >> tmp_ergo_arg_test_mol.xyz
echo H      0.969771     0.000000    -0.290392 >> tmp_ergo_arg_test_mol.xyz
echo H     -0.390071     0.887881    -0.290336 >> tmp_ergo_arg_test_mol.xyz
echo basis = \"6-31Gss\"                        > tmp_ergo_arg_test_input_file.ego
echo run \"HF\"                                >> tmp_ergo_arg_test_input_file.ego
./ergo -m tmp_ergo_arg_test_mol.xyz -e "charge = 1" -e "spin_polarization = 1" tmp_ergo_arg_test_input_file.ego > /dev/null
if
check_final_energy -55.8546235 1e-6 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
if
check_final_S2 0.757013 1e-6 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi

# Do the same test again, but this time with the molecule as a .mol file specified using a -m argument

echo Testing nh3[+] UHF/6-31G** with -e args, -m molecule mol arg, and regular input file
rm -f ergoscf.out
echo BASIS                                      > tmp_ergo_arg_test_mol.mol
echo STO-3G                                    >> tmp_ergo_arg_test_mol.mol
echo Test calculation                          >> tmp_ergo_arg_test_mol.mol
echo -----------------------------------       >> tmp_ergo_arg_test_mol.mol
echo Atomtypes=2 Angstrom                      >> tmp_ergo_arg_test_mol.mol
echo        7.    1                            >> tmp_ergo_arg_test_mol.mol
echo N      0.000000     0.000000     0.000000 >> tmp_ergo_arg_test_mol.mol
echo        1.    3                            >> tmp_ergo_arg_test_mol.mol
echo H      0.000000     0.000000     1.012316 >> tmp_ergo_arg_test_mol.mol
echo H      0.969771     0.000000    -0.290392 >> tmp_ergo_arg_test_mol.mol
echo H     -0.390071     0.887881    -0.290336 >> tmp_ergo_arg_test_mol.mol
echo basis = \"6-31Gss\"                        > tmp_ergo_arg_test_input_file.ego
echo run \"HF\"                                >> tmp_ergo_arg_test_input_file.ego
./ergo -m tmp_ergo_arg_test_mol.mol -e "charge = 1" -e "spin_polarization = 1" tmp_ergo_arg_test_input_file.ego > /dev/null
if
check_final_energy -55.8546235 1e-6 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
if
check_final_S2 0.757013 1e-6 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi


# Do the same test again, but this time with two different input files

echo Testing nh3[+] UHF/6-31G** with -e arguments and two input files
rm -f ergoscf.out
echo molecule_inline Angstrom                   > tmp_ergo_arg_test_input_file_1.ego
echo N      0.000000     0.000000     0.000000 >> tmp_ergo_arg_test_input_file_1.ego
echo H      0.000000     0.000000     1.012316 >> tmp_ergo_arg_test_input_file_1.ego
echo H      0.969771     0.000000    -0.290392 >> tmp_ergo_arg_test_input_file_1.ego
echo H     -0.390071     0.887881    -0.290336 >> tmp_ergo_arg_test_input_file_1.ego
echo EOF                                       >> tmp_ergo_arg_test_input_file_1.ego
echo basis = \"6-31Gss\"                        > tmp_ergo_arg_test_input_file_2.ego
echo run \"HF\"                                >> tmp_ergo_arg_test_input_file_2.ego
./ergo -e "charge = 1" -e "spin_polarization = 1" tmp_ergo_arg_test_input_file_1.ego tmp_ergo_arg_test_input_file_2.ego > /dev/null
if
check_final_energy -55.8546235 1e-6 ;
then
echo OK
else
echo ERROR
mv ergoscf.out $errorfilename
echo output file saved as $errorfilename
exit 1
fi
if
check_final_S2 0.757013 1e-6 ;
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
echo Input arguments tests completed successfully!
echo
