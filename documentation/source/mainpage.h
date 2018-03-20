/* Ergo, version 3.6, a program for linear scaling electronic structure
 * calculations.
 * Copyright (C) 2017 Elias Rudberg, Emanuel H. Rubensson, Pawel Salek,
 * and Anastasia Kruchinina.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Primary academic reference:
 * Kohn−Sham Density Functional Theory Electronic Structure Calculations 
 * with Linearly Scaling Computational Time and Memory Usage,
 * Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek,
 * J. Chem. Theory Comput. 7, 340 (2011),
 * <http://dx.doi.org/10.1021/ct100611z>
 * 
 * For further information about Ergo, see <http://www.ergoscf.org>.
 */
/** @file mainpage.h Main page in documentation
 *  @author Emanuel Rubensson <em>responsible author</em>
 *  @version 0.1
 *  @date    March 2006
 */

/** @mainpage ergo Documentation

@htmlonly <hr> @endhtmlonly

 @section intro_sec Introduction

 @htmlonly

 <p>ErgoSCF.org is the home of Ergo, a quantum chemistry program for
 large-scale self-consistent field calculations.</p>

 Key features of the Ergo program:
 <ul>
 <li>Performs electronic structure calculations using Hartree-Fock and Kohn-Sham density functional theory.</li>
 <li>Written in C++.
 <li>Uses Gaussian basis sets.</li>
 <li>Both core and valence electrons are included in the calculations.</li>
 <li>Both restricted and unrestricted models are implemented for energy calculations.</li>
 <li>Implements a broad range of both pure and hybrid Kohn-Sham density functionals.</li>
 <li>Employs modern linear scaling techniques like fast multipole methods, hierarchic sparse matrix algebra, 
 density matrix purification, and efficient integral screening.</li>
 <li>Linear scaling is achieved not only in terms of CPU usage but also memory utilization.</li>
 <li>The time consuming parts of the code are currently parallelized using the shared-memory paradigm.</li>
 </ul>

 Linear response calculations of polarizabilities and excitation
 energies are possible for the restricted reference density, although
 complete linear scaling is in the current implementation not achieved
 since full dense matrices are still used in parts of the linear
 response implementation.

 <br> <br>

 @endhtmlonly

 @htmlonly <hr> @endhtmlonly

 @section License

 @verbinclude ergo_license_long.txt

 @htmlonly <hr> @endhtmlonly

 @section runergo How to run ergo

The executable "ergo" is the new interface to the Ergo project which
is meant to be scriptable.
<pre>
Usage: ergo [args...]
args can be: input file name
             -e  "input line"
             -d  VARIABLE_NAME variable to describe
             -m  molecule file name
             -c  extra-charges molecule file name
             -h help message that lists also all the available variables.
</pre>
The statements in the input file are divided in two classes: variable
assignments and executable statements (commands). Currently, following
commands are recognized:

run - runs an SCF calculation: <pre>run "HF"</pre>
molecule_inline - defines a molecule
<pre>molecule_inline
C 0 0 0
O 0 0 2.3
EOF
</pre>
molecule - reads the molecule file in a MOLECULE/Dalton or XYZ file format
<pre>
molecule "../nh3.mol"
</pre>

get_excited_state - computes a set of excited states
<pre>get_excited_state "CAMB3LYP" 4</pre>

get_polarisability - computes a polarizability for given frequency.
<pre>get_polarisability "PBE" "X" 0.2</pre>

system - executes a system command
<pre>system "rm density.old; mv density.bin density.old"</pre>

Example of a simplest input file:
<pre>
basis= "6-31Gs"
molecule_inline
O     0 0 0
H     1.2 1.2 0
H    -1.2 1.2 0
EOF
get_polarisability "HF" "Y" 0.01
</pre>

 @htmlonly <hr> @endhtmlonly

 @section gendoc How to generate this documentation
 Run doxygen without arguments in the ergo directory. 
 This will generate html documentation which can be browsed
 via ergo/documentation/html/main.html 
 
@htmlonly <hr> @endhtmlonly

 @section comments How to write comments 
 
 See the Doxygen manual http://www.stack.nl/~dimitri/doxygen/manual.html
 on how to get started writing Doxygen comments. @n
 Comments should be written in JavaDoc style, i.e. using the
 @@ comment prefix rather than @\ and using 
 @verbatim
 /**
  *
  */
 @endverbatim
 rather than 
 @verbatim
 /*!
  *
  */
 @endverbatim
 for doxygen comments.

 We will possibly add latex doc generation in the future so please keep 
 this in mind and use @@htmlonly for all html specific commands.

 
@subsection filecomments File comments
 Each file should begin with a header in the following style:

@verbatim
 /** @file filename.h
     @brief Brief description
     @author Author1 <em>responsible</em>
     @author Author2
  */
@endverbatim
which will result in the following: @n
Brief description
@author Author1   <em>responsible</em>
@author Author2

@htmlonly <hr> @endhtmlonly

 @section install_sec Installation

 @subsection step1 Usual configuration
  
./configure && make
  
Options can be passed - see @ref examples for more information. 

@subsection correctness Correctness test

make check

Verbose variant:
make check VERBOSE=1
  
 @subsection examples Configuration examples 
  
How to install on different machines and libraries (this should probably be updated)

@verbinclude config_examples.txt


@htmlonly <hr> @endhtmlonly

@section knownbugs_sec Known bugs

No known bugs right now.


 * 
 *  
 * 
 *
 *
 *
 *
 * 
 * 
 *
 * 
 */
