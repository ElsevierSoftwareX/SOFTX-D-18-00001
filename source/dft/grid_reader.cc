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

/** @file grid_reader.cc

    @brief Grid Generator interface. Functions for opening grid file,
    reading chunks from it, and closing the file, are provided.

    @author: Pawel Salek <em>responsible</em>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>

#include "dft_common.h"
#include "grid_reader.h"
#include "grid_stream.h"
#include "grid_hicu.h"
#include "grid_matrix.h"

  class FullMatrixWrapper : public Dft::Matrix {
    const Dft::FullMatrix& matrix;
  public:
    explicit FullMatrixWrapper(const Dft::FullMatrix& m) : matrix(m) {}
    virtual ergo_real at(int row, int col) const 
    { return matrix.mat[row + col*matrix.nbast]; }
    virtual bool isSparse() const { return false; }
    virtual const Dft::SparseMatrix* asSparse() const 
    { throw std::runtime_error("Cannot converse full to sparse"); }
    virtual const ergo_real* asFull() const 
    { return matrix.mat; }
  };

  class SparseMatrixWrapper : public Dft::Matrix {
    const Dft::SparseMatrix& matrix;
  public:
    explicit SparseMatrixWrapper(const Dft::SparseMatrix& m) : matrix(m) {}
    virtual ergo_real at(int row, int col) const 
    { return matrix.at(row, col); }
    virtual bool isSparse() const { return true; }
    virtual const Dft::SparseMatrix* asSparse() const
    { return & matrix; }
    virtual const ergo_real* asFull() const 
    { throw std::runtime_error("cannot convert sparse to full"); }
  };


Dft::Matrix*
createGridMatrix(const Dft::FullMatrix& mat)
{
  return new FullMatrixWrapper(mat);
}

Dft::Matrix*
createGridMatrix(const Dft::SparseMatrix& mat)
{
  return new SparseMatrixWrapper(mat);
}


/* ===================================================================
   GENERAL GRID MANIPULATION ROUTINES.
   =================================================================== */
/* grid_get_fname, grid_set_tmpdir: helper routines for creating
 * unique grid file name and setting the temporary directory used for
 * storing the grid file. */

std::string grid_tmpdir;
static char*
grid_get_fname(const char *base, int filenum)
{
    char *res;
    if(grid_tmpdir.empty()) {
        char *tmpdir = getenv("TMPDIR");
        if(tmpdir)
            grid_tmpdir = tmpdir;
        else
            grid_tmpdir = ".";
    }
    res = (char*)malloc(strlen(base) + grid_tmpdir.length() + 15);
    sprintf(res, "%s/%s.%06u.%05d", grid_tmpdir.c_str(),
            base, (unsigned)getpid(), filenum);
    return res;
}

void
grid_set_tmpdir(const char *tmpdir)
{
    grid_tmpdir = tmpdir ? tmpdir : ".";
}

#define GRID_BASE_NAME "ERGO-grid"
#define GRID_PATT_NAME "ERGO-patt"

/* shared grid file */
static FILE *grid_file = NULL;
static int grid_file_open_count = 0;
static char *grid_file_name = NULL;
static char *patt_file_name = NULL;
static pthread_mutex_t grid_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t grdone_mutex = PTHREAD_MUTEX_INITIALIZER;


/** Frees all the cached data if any. */
void
grid_free_files()
{
    if(grid_file_name) {
        unlink(grid_file_name);
        free(grid_file_name);
        grid_file_name = NULL;
    }
    if(patt_file_name) {
        unlink(patt_file_name);
        free(patt_file_name);
        patt_file_name = NULL;
    }
}


static int grid_atexit_registered = 0;
static void
grid_atexit(void)
{
    grid_free_files();
}



bool 
grid_is_ready() {
  if (grid_file_name)
    return true;
  return false;
}

struct DftGridReader {
    FILE *f;
};

static const int MY_MPI_NUM = 0; /* no MPI for now. */
static void
grid_open_stream(const class GridGenMolInfo& molInfo,
		 const Dft::GridParams& gss,
		 Dft::SparsePattern *pattern,
		 DftGridReader *reader)
{
  pthread_mutex_lock(&grdone_mutex);
  if (!grid_file_name) {
    grid_file_name = grid_get_fname(GRID_BASE_NAME, MY_MPI_NUM);
    ErgoGridStream *egStream = grid_stream_new(gss, molInfo);
    if(pattern)
      grid_stream_set_sparse_pattern(egStream, pattern);
    grid_stream_generate(egStream, grid_file_name,
                                   dft_get_num_threads());
    grid_stream_free(egStream);
    if(pattern) {
      patt_file_name = grid_get_fname(GRID_PATT_NAME, MY_MPI_NUM);
      FILE  *f;
      if( (f = fopen(patt_file_name, "wb")) == NULL)
        throw "Cannot open pattern file for writing";
      pattern->save(f);
      if(fclose(f) != 0)
        throw "Cannot close the sparse pattern file";
    }
  }
  
  /* Grid generated if needed, setup time now! */
  if(pattern) {
    FILE  *f;
    if(!patt_file_name)
      throw "Dft::SparsePattern requested but unavailable.";
    if( (f = fopen(patt_file_name, "rb")) == NULL)
      throw "Cannot open pattern file for reading";
    pattern->load(f);
    if(fclose(f) != 0)
      throw "Cannot close the sparse pattern file";
  }
  if(!grid_file)
    grid_file = fopen(grid_file_name, "rb");
  grid_file_open_count++;
  pthread_mutex_unlock(&grdone_mutex);
  reader->f= grid_file;
  if(grid_file == NULL) {
    perror("grid_open_standard: DFT quadrature grid file " GRID_BASE_NAME
           " not found");
    free(reader);
    abort();
  }
}

static void
grid_open_cartesian(const BasisInfoStruct& bis,
		    const Dft::GridParams& gss,
                    const Dft::Matrix* dmat,
                    Dft::SparsePattern *pattern, DftGridReader *reader)
{
  pthread_mutex_lock(&grdone_mutex);
  if (!grid_file_name) {
    grid_file_name = grid_get_fname(GRID_BASE_NAME, MY_MPI_NUM);
    hicu_grid_generate(grid_file_name, bis, 
		       gss.hicuParams.maxError,
		       gss.hicuParams.box_size,
		       gss.hicuParams.start_box_size_debug,
		       gss.hicuParams.use_error_per_volume,
		       gss.hicuParams.do_double_checking,
		       gss.hicuParams.compare_to_refined,
		       gss.hicuParams.use_energy_criterion,
		       gss.hicuParams.use_energy_criterion_only,
		       gss.hicuParams.do_variation_checking,
                       dmat, pattern, dft_get_num_threads(), false);
    if(pattern) {
      patt_file_name = grid_get_fname(GRID_PATT_NAME, MY_MPI_NUM);
      FILE  *f;
      if( (f = fopen(patt_file_name, "wb")) == NULL)
        throw "Cannot open pattern file for writing";
      pattern->save(f);
      if(fclose(f) != 0)
        throw "Cannot close the sparse pattern file";
    }
  }
  /* Grid generated if needed, setup time now! */
  if(pattern) {
    FILE  *f;
    if(!patt_file_name)
      throw "Dft::SparsePattern requested but unavailable.";
    if( (f = fopen(patt_file_name, "rb")) == NULL)
      throw "Cannot open pattern file for reading";
    pattern->load(f);
    if(fclose(f) != 0)
      throw "Cannot close the sparse pattern file";
  }

  if(!grid_file)
    grid_file = fopen(grid_file_name, "rb");
  grid_file_open_count++;
  pthread_mutex_unlock(&grdone_mutex);
  reader->f= grid_file;
  if(grid_file == NULL) {
    perror("grid_open_cartesian: DFT quadrature grid file " GRID_BASE_NAME
           " not found");
    free(reader);
    abort();
  }
}

/** Returns a handle to a grid file. Sets the sparse pattern if
    passed. Observe that sparse pattern must be passed the first time
    to get generated. Otherwise, subsequent calls will not be able to
    set it. */

DftGridReader*
grid_open_full(const class GridGenMolInfo *mol_info,
               const Dft::GridParams& gss,
               Dft::SparsePattern *pattern,
               const Dft::Matrix* dmat,
               const BasisInfoStruct& bis)
{
    DftGridReader *res = dal_new(1,DftGridReader);
    if(!grid_atexit_registered) {
        atexit(grid_atexit);
        grid_atexit_registered = 1;
    } 

    switch(gss.gridType) {
    case Dft::GridParams::TYPE_STANDARD:
      grid_open_stream(*mol_info, gss, pattern, res);
      return res;
    case Dft::GridParams::TYPE_HICU:
      grid_open_cartesian(bis, gss, dmat, pattern, res);
      return res;
    default:
      perror("Error in grid_open: unknown grid type\n");
      free(res);
      abort();
    } /* END SWITCH */
}

/** grid_getchunk_blocked() reads grid data also with screening
    information if only nblocks and shlblocks are provided.

    @param rawgrid shared grid handle.
    @param  maxlen the upper limit on the grid point chunk length.

    @param nBlocks will contain number of active b.f. blocks. May be
                   NULL if uninteresting.

    @param shlBlocks pointer to the shell block range.

    @param    coor array with grid point coordinates.
    @param  weight array with grid point weights.

    @return number of read grid points. -1 on end-of-file.
 */
int
grid_getchunk_blocked(DftGridReader* rawgrid, int maxlen,
                      int *nBlocks, int *shlBlocks, 
                      real (*coor)[3], real *weight)
{
    int points = 0, rc, bl_cnt;
    FILE *f;

    pthread_mutex_lock(&grid_mutex);
    f = rawgrid->f;
    if(fread(&points, sizeof(int), 1, f) <1) {
        points = -1; /* end of file */
        goto unlock_and_quit;
    }
    if(points>maxlen) {
        fprintf(stderr,
                "grid_getchunk: too long vector length in file: %d > %d\n"
                "Calculation will stop.\n", points, maxlen);
        throw "grid_getchunk: too long vector length in file.";
        points = -1; /* stop this! */
        goto unlock_and_quit;
    }

    if(fread(&bl_cnt, sizeof(unsigned), 1, f) <1) {
        puts("OCNT reading error."); points = -1; goto unlock_and_quit;
    }
    if(nBlocks) *nBlocks = bl_cnt;

    if(shlBlocks) {
        rc = fread(shlBlocks, sizeof(int), bl_cnt*2, f);
    } else {
        int buf, cnt;
        rc = 0;
        for(cnt=0; cnt<bl_cnt*2; cnt+=rc) {
            rc = fread(&buf, sizeof(int), 1, f);
            if(rc < 1)
                break;
        }
    }
    if(rc<1) {
        fprintf(stderr,
                "IBLOCKS reading error: failed to read %d blocks.\n",
                *nBlocks);
        points = -1; goto unlock_and_quit;
    }

    if((signed)fread(coor, sizeof(real), 3*points, f) < 3*points) { 
        puts("reading grid point coords failed"); points = -1;
        goto unlock_and_quit;
    }
    if((signed)fread(weight, sizeof(real), points, f) < points) {
        puts("reading grid weights failed."); points = -1;
    }
    unlock_and_quit:
    pthread_mutex_unlock(&grid_mutex);
    return points;
}


/** Closes the shared grid handle that is specifed as the argument. */

void
grid_close(DftGridReader *rawgrid)
{
    pthread_mutex_lock(&grid_mutex);
    if(--grid_file_open_count == 0) {
        fclose(grid_file);
        grid_file = NULL;
    }
    pthread_mutex_unlock(&grid_mutex);
    free(rawgrid);
}
