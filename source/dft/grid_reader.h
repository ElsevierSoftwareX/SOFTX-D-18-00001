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

/** @file grid_reader.h

    @brief Grid Generator interface. Functions for opening grid file,
    reading chunks from it, and closing the file, are provided.

    @author: Pawel Salek <em>responsible</em>
*/

/* -*-mode:c; c-style:bsd; c-basic-offset:4;indent-tabs-mode:nil; -*- */

#if !defined(_GRID_READER_H_)
#define _GRID_READER_H_ 1

#include "sparse_pattern.h"
#include "grid_stream.h"
#include "grid_interface.h"
#include "grid_params.h"
#include "grid_matrix.h"

struct DftGridReader;

Dft::Matrix* createGridMatrix(const Dft::FullMatrix& mat);
Dft::Matrix* createGridMatrix(const Dft::SparseMatrix& mat);

DftGridReader* grid_open_full(const class GridGenMolInfo *mol_info,
                              const Dft::GridParams& gss,
                              Dft::SparsePattern *pattern,
                              const Dft::Matrix* dmat,
                              const BasisInfoStruct& bis);

bool grid_is_ready();

int grid_getchunk_blocked(DftGridReader* grid_handle, int maxlen,
			  int *nblocks, int *shlblocks, 
			  real (*coor)[3], real *weight);

#define grid_getchunk_plain(r,m,coor,w) \
       (grid_getchunk_blocked((r),(m),NULL,NULL,(coor),(w)))
void grid_close(DftGridReader *rawgrid);
void grid_free_files();
void grid_set_tmpdir(const char *tmpdir);

#endif /* !defined(_GRID_READER_H_) */

 
