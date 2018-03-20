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

/** @file sparse_matrix.cc The implementation of sparse matrix optimized
    for XC integration.

    Notes: the atom reordering issues are not discussed. Some
    preliminary experiments suggest that reordering may give 20%
    speedup. The permutation speedup remain to be thoroughly tested.
*/

#include <string.h>

#include <list>
#include <map>
#include <set>
#include <vector>
#include <string.h>

#include "output.h"
#include "dft_common.h"
#include "sparse_matrix.h"

BEGIN_NAMESPACE(Dft)


/** computes a squared distance between two points. */
static inline ergo_real
sqDist(const ergo_real a[], const ergo_real b[])
{
  ergo_real res = 0;
  for(int i=0; i<3; i++) {
    ergo_real d = a[i]-b[i];
    res += d*d;
  }
  return res;
}

class NeighbourList {
  const ShellSpecStruct* shellInfo;
  std::list<int> neighbours;
  ergo_real extent; /**< an approximation for the shell extent. */
public:
  NeighbourList(const ShellSpecStruct* sis, ergo_real thr)
    : shellInfo(sis), extent(0)
  {
    /* We should really check what is left out, an additional term
       that is very relevant for functions with low exponent is
       skipped. */
    for(int i=0; i<sis->noOfContr; i++) {
      ergo_real rad = template_blas_sqrt(-template_blas_log(thr/template_blas_fabs(sis->coeffList[i]))
                                /sis->exponentList[i]);
      if (rad>extent)
        extent = rad;
    }
  }

  void setOverlappingWith(const std::vector<NeighbourList>& list) {
    for(unsigned i=0; i < list.size(); ++i) {
      const ShellSpecStruct& otherShell = *list[i].shellInfo;
      ergo_real sqDistance = sqDist(shellInfo->centerCoords,
                                    otherShell.centerCoords);
      ergo_real minDist = list[i].extent + extent;
      //printf("Comparing to shell # %i  - %f away\n", i, sqrt(sqDistance));
      if(sqDistance<minDist*minDist)
        neighbours.push_back(i);
    }
#if 0
    printf("Radius %f . Found %d neighbours\n",
           radius, neighbours.size());
#endif
  }
  std::list<int>::iterator begin() {
    return neighbours.begin();
  }
  std::list<int>::iterator end() {
    return neighbours.end();
  }
  size_t size() const {
    return neighbours.size();
  }
};


#if 0
static int
findExtremeShell(int noOfShells, const ShellSpecStruct* shellList)
{
  static const ergo_real origin[] = { 0.0, 0.0, 0.0 };

  ergo_real rMax = sqDist(shellList[0].centerCoords, origin);
  int iMax = 0;

  for(int i=1; i<noOfShells; i++) {
    ergo_real r = sqDist(shellList[i].centerCoords, origin);
    if (r>rMax) {
      rMax = r;
      iMax = i;
    }
  }
  /* printf("Extreme Shell: %i\n", iMax); */
  return iMax;
}
#endif


typedef ergo_real *ErgoRealPtr;
SparseMatrix::SparseMatrix(const SparsePattern& pattern_)
  : pattern(pattern_), columns(new ErgoRealPtr[pattern_.size()]),
    n(pattern_.size())
{
  for(int col=0; col<n; col++) {
    int colSize = pattern.getColumnSize(col);
    columns[col] = new ergo_real[colSize];
    for(int row=0; row<colSize; row++)
      columns[col][row] = 0.0;
  }
  createOffsets(pattern);
}


SparseMatrix::SparseMatrix(const SparsePattern& pattern_,
                           const symmMatrix& sMat, const int *aoMap,
                           std::vector<int> const & permutationHML)
  : pattern(pattern_), columns(new ErgoRealPtr[pattern_.size()]),
    n(pattern_.size())
{
  std::vector<int> rowI(n), colI(n);

  for(int col=0; col<n; col++) {
    int colSize = pattern.getColumnSize(col);
    colI.clear();
    rowI.clear();
    columns[col] = new ergo_real[colSize];
    const SparsePattern::Column& column = pattern[col];
    SparsePattern::Column::Iterator colEnd = column.end();
    for(SparsePattern::Column::Iterator row = column.begin();
        row != colEnd; ++row) {
      int pRow = aoMap[*row], pCol = aoMap[col];
      if(col == 0)
      if(pRow<pCol) {
        int t = pRow; pRow = pCol; pCol = t;
      }
      rowI.push_back(pRow);
      colI.push_back(pCol);
    }
    std::vector<ergo_real> columnsTmp;
    sMat.get_values(rowI, colI, columnsTmp, 
                    permutationHML, permutationHML);
    /* FIXME: avoid slow copy somehow */
    std::copy(columnsTmp.begin(), columnsTmp.end(),columns[col]);
#if 0
    printf("get_values() for column %d returned:\n", col);
    for(int row=0; row<colSize; row++)
      printf("%d %d : %f\n", rowI[row], colI[row], columns[col][row]);
#endif
  }
  createOffsets(pattern);
  //print("Densmat");
}

void
SparseMatrix::createOffsets(const SparsePattern& patt)
{
  offsets = new int*[n];
  his = new int*[n];
  cnt = new int[n];
  for(int col=0; col<n; col++) {
    const SparsePattern::IntervalList& intervalList = patt[col].list;
    int numberOfIntervals = intervalList.size();
    int * off = offsets[col] = new int[numberOfIntervals];
    int * hi  = his[col] = new int[numberOfIntervals];
    int offset = 0, last = 0, idx=0;

    for(SparsePattern::IntervalList::const_iterator
          i = intervalList.begin();
        i != intervalList.end(); ++i) {
      offset += i->lo - last; 
      off[idx] = offset;
      hi[idx] = last = i->hi;
      ++idx;
    }
    cnt[col] = numberOfIntervals;
  }
}

#if 1
void
SparseMatrix::print(const char *title) const
{
  puts(title);
  for(int row=0; row<n; row++) {

    const SparsePattern::Column& col2 = pattern[row];
    int prevCol = -1;
    int offset = 0;
    for(SparsePattern::Column::Iterator col = col2.begin();
        col != col2.end(); ++col) {
      /* Take care of skipped columns, if any. */
      for(++prevCol; prevCol < *col; prevCol++)  printf("********* ");
      
      printf("%9.5f ", (double)columns[row][offset]);
      offset++;
    }
    /* Complete tailing columns, if any. */
    for(++prevCol; prevCol < n; prevCol++)  printf("********* ");
    puts("");
  }
}
#endif

void
SparseMatrix::addSymmetrizedTo(symmMatrix& sMat, const int *aoMap, 
                               std::vector<int> const & permutationHML) const
{
  const unsigned BUF_SIZE = 10*n;
  std::vector<int> rowI(BUF_SIZE), colI(BUF_SIZE);
  std::vector<ergo_real> buf(BUF_SIZE);

  colI.clear();
  rowI.clear();
  buf.clear();

  for(int col=0; col<n; col++) {
    const ergo_real *column = columns[col];
    const SparsePattern::Column& patternCol = pattern[col];
    int offset = 0;
    SparsePattern::Column::Iterator colEnd = patternCol.end();
    for(SparsePattern::Column::Iterator row= patternCol.begin();
        row != colEnd; ++row) {
      if(buf.size() == BUF_SIZE) {
        sMat.add_values(rowI, colI, buf, permutationHML, permutationHML);
        colI.clear();
        rowI.clear();
        buf.clear();
      }
      rowI.push_back(aoMap[*row]);
      colI.push_back(aoMap[col]);
      buf.push_back( column[offset]);
      offset++;
    }
  }
  if(!buf.empty())
    sMat.add_values(rowI, colI, buf, permutationHML, permutationHML);
  //print("XC mat");
}

END_NAMESPACE(Dft)

typedef ergo_real real;
static void
zeroorbs(real *tmp, const int *nblocks, const int (*iblocks)[2], int ldaib, int nvclen)
{
  /* DIMENSION TMP(NVCLEN,NBAST),NBLOCKS(NSYM),IBLOCKS(2,LDAIB,NSYM) */
  int ibl, idx, k; 
  for(ibl=0; ibl<nblocks[0]; ibl++)
    for(idx=iblocks[ibl][0]; idx<iblocks[ibl][1]; idx++) {
      real * tmpi = tmp + idx*nvclen;
      for(k=0; k<nvclen; k++) tmpi[k] = 0.0;
    }
}

void
getrho_blocked_lda(int nbast, const Dft::SparseMatrix& dmat,
                   const ergo_real * gao,
                   const int* nblocks, const int (*iblocks)[2],
                   int ldaib, ergo_real *tmp, int nvclen, ergo_real *rho)
{
  /*
    DIMENSION DMAT(NBAST,NBAST), GAO(NVCLEN,NBAST,*)
    DIMENSION NBLOCKS(NSYM), IBLOCKS(2,LDAIB,NSYM), RHO(NVCLEN)
    DIMENSION TMP(NVCLEN,NBAST)
  */
  int colBl, col, rowBl, row, k;
  real * tmpj, d;

  zeroorbs(tmp, nblocks, iblocks, ldaib, nvclen);
    
  for(colBl=0; colBl<nblocks[0]; colBl++)
    for(col=iblocks[colBl][0]; col<iblocks[colBl][1]; col++) {
      const real * gaoi = gao + col*nvclen;
      for(rowBl=0; rowBl<nblocks[0]; rowBl++) {
        int jtop = iblocks[rowBl][1] > col ? col : iblocks[rowBl][1];
        for(row=iblocks[rowBl][0]; row<jtop; row++) {
          d = dmat.at(row, col);
          tmpj = tmp + row*nvclen;
          for(int k=0; k<nvclen; k++)
            tmpj[k] += gaoi[k]*d;
        }
      }
      
      tmpj = tmp + col*nvclen;
      d = dmat.at(col,col)*0.5;
      for(k=0; k<nvclen; k++)
        tmpj[k] += gaoi[k]*d;

    }

  memset(rho, 0, nvclen*sizeof(real));
  for(colBl=0; colBl<nblocks[0]; colBl++)
    for(col=iblocks[colBl][0]; col<iblocks[colBl][1]; col++) {
      const real * gaoi = gao + col*nvclen;
      const real * tmpi = tmp + col*nvclen;
      for(k=0; k<nvclen; k++)
        rho[k] += gaoi[k]*tmpi[k]*2;
    }
}

void
getrho_blocked_gga(int nbast, const Dft::SparseMatrix& dmat,
                   const ergo_real * gao,
                   const int* nblocks, const int (*iblocks)[2],
                   int ldaib, ergo_real *tmp, int nvclen,
                   ergo_real *rho, ergo_real (*grad)[3])
{
  /*
    DIMENSION DMAT(NBAST,NBAST), GAO(NVCLEN,NBAST,*)
    DIMENSION NBLOCKS(NSYM), IBLOCKS(2,LDAIB,NSYM), RHO(NVCLEN)
    DIMENSION TMP(NVCLEN,NBAST*4)
  */
  int colBl, col, rowBl, row, k;
  real * tmpj, d;

  zeroorbs(tmp, nblocks, iblocks, ldaib, nvclen);
    
  for(colBl=0; colBl<nblocks[0]; colBl++)
    for(col=iblocks[colBl][0]; col<iblocks[colBl][1]; col++) {
      const real * gaoi = gao + col*nvclen;
      for(rowBl=0; rowBl<nblocks[0]; rowBl++) {
        for(row=iblocks[rowBl][0]; row<iblocks[rowBl][1]; row++) {
          d = dmat.at(row, col);
          tmpj = tmp + row*nvclen;
          for(int k=0; k<nvclen; k++)
            tmpj[k] += gaoi[k]*d;
        }
      }
    }

  memset(rho,  0,   nvclen*sizeof(real));
  memset(grad, 0, 3*nvclen*sizeof(real));
  for(colBl=0; colBl<nblocks[0]; colBl++)
    for(col=iblocks[colBl][0]; col<iblocks[colBl][1]; col++) {
      const real * gaoi = gao + col*nvclen;
      const real * tmpi = tmp + col*nvclen;
      for(k=0; k<nvclen; k++) {
        rho[k]     += gaoi[k]*tmpi[k];
        grad[k][0] += gaoi[k + nvclen*nbast  ]*tmpi[k]*2;
        grad[k][1] += gaoi[k + nvclen*nbast*2]*tmpi[k]*2;
        grad[k][2] += gaoi[k + nvclen*nbast*3]*tmpi[k]*2;
      }
    }
}
