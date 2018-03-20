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

/** @file sparse_pattern.h

    @brief Class that can be used to store sparse matrix patterns.

    @author: Pawel Salek <em>responsible</em>
*/

#if !defined(_DFT_SPARSE_PATTERN_H_)
#define _DFT_SPARSE_PATTERN_H_ 1

#if !defined(BEGIN_NAMESPACE)
#define BEGIN_NAMESPACE(x) namespace x {
#define END_NAMESPACE(x)   } /* x */
#endif

#include <vector>
#include <stdio.h>

#include "basisinfo.h"

BEGIN_NAMESPACE(Dft)

/** A way to store sparse matrix patterns */
class SparsePattern {
 public:
  /** ranges are upper-exclusive: involve i: lo <= i < hi. */
  struct Interval {
    int lo, hi;
  Interval(int l_, int h_) : lo(l_), hi(h_){}
  };
  typedef std::vector<Interval> IntervalList;
  struct Column {
    IntervalList list;

    void addInterval(int lo, int hi);
    void addIntervals(int nIntervals, int (*intervals)[2]);
    struct Iterator {
      IntervalList::const_iterator current, end;
      int pos;
      Iterator(const IntervalList::const_iterator& beg,
               const IntervalList::const_iterator& end_, int p)
        : current(beg), end(end_), pos(p)
      {}

      Iterator& operator++() {
        ++pos;
#if 0
        if(pos == current->hi)
          printf("Iterator increased to %d current limit %d last? %s %s\n",
                 pos, current->hi,
                 & *current == & *end ? "YES" : "NO",
                 current == end ? "YES" : "NO");
#endif
        if(pos >= current->hi) {
          ++current;
          if(current != end)
            pos = current->lo; 
          else pos = 0;
        }
        return *this;
      }
      bool operator!=(const Iterator& other) const {
        bool res = !(& *current == & *other.current && pos == other.pos);
#if 0
        printf("Iterator::operator!=() compares %p with %p, returns %s \n",
               & *current, & *other.current, res ? "TRUE" : "FALSE");
#endif
        return res;
      }
      int operator*() const {
        //printf("Iterator::operator*() returns %d\n", pos);
        return pos;
      }
      const Interval* operator->() const {
        return &(*current);
      }
      
    };

    Iterator begin() const {
      IntervalList::const_iterator a = list.begin();
      IntervalList::const_iterator b = list.end();
      return Iterator(a, b, a != list.end() ? a->lo : 0);
    }

    Iterator end() const {
      return Iterator(list.end(),list.end(),0);
    }

    int size() const {
      int result = 0;
      for(IntervalList::const_iterator i = list.begin();
          i != list.end(); ++i)
        result += i->hi- i->lo;
      return result;
    }
  };

 private:
  const BasisInfoStruct& bis;
  Column *ranges;
 public:
 explicit SparsePattern(const BasisInfoStruct& bis_)
    : bis(bis_), ranges(new Column[bis_.noOfBasisFuncs])
    { }
    
    ~SparsePattern() {
      delete []ranges;
    }

    /** marks specified ranges as used.
    */
    void add(int nRanges, const int (*range)[2]);

    void save(FILE *f) const;
    void load(FILE *f);
    const Column& operator[](int column) const {
      return ranges[column];
    }

    /** returns the number of stored elements for specified column. */
    int getColumnSize(int col) const {
      return ranges[col].size();
    }

    /** Returns the dimension of the pattern. Auxiliary function. */
    int size() const {
      return bis.noOfBasisFuncs;
    }
    /** returns the total number of nonzero elements. */
    int sizeTotal() const;
};

void setupShellMap(const BasisInfoStruct& bis, int *shellMap, int *aoMap);

END_NAMESPACE(Dft)

#endif /* _DFT_SPARSE_PATTERN_H_ */
