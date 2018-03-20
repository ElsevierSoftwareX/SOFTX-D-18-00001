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

/** @file basisset.h

    \brief Code for representing basis set information for Gaussian
    basis sets, and for parsing a text file specifying such a
    basisset.

    @author: Elias Rudberg <em>responsible</em>.
*/

/* -*-mode:c; c-style:k&r; c-basic-offset:4; indent-tabs-mode: nil -*- */
#ifndef BASISSET_HEADER
#define BASISSET_HEADER

#include <vector>
#include "realtype.h"
#include "polydegree.h"

#define MAX_NO_OF_ATOM_TYPES 100

#ifndef BASIS_FUNC_POLY_MAX_DEGREE
#error The constant BASIS_FUNC_POLY_MAX_DEGREE must be defined.
#endif
#if BASIS_FUNC_POLY_MAX_DEGREE<6
#define MAX_NO_OF_SHELLS_PER_ATOM 44
#else
#define MAX_NO_OF_SHELLS_PER_ATOM 88
#endif

#define MAX_NO_OF_CONTR 44

typedef struct
{
  int type;
  int contrCount;
  int shell_ID;
  ergo_real exponentList[MAX_NO_OF_CONTR];
  ergo_real coeffList[MAX_NO_OF_CONTR];
} basisset_shell_struct;

typedef struct
{
  int noOfShells;
  basisset_shell_struct shells[MAX_NO_OF_SHELLS_PER_ATOM];
} basisset_atom_struct;

struct basisset_info {
  std::vector<basisset_atom_struct> atoms;
  basisset_info();
  void clear();
  // Stuff needed for Chunks&Tasks usage
  void write_to_buffer ( char * dataBuffer, size_t const bufferSize ) const;
  size_t get_size() const;
  void assign_from_buffer ( char const * dataBuffer, size_t const bufferSize);
};

int read_basisset_file(basisset_info & result, 
		       const char* fileName,
		       int dirc, 
		       const char *dirv[],
		       int print_raw);


#endif
