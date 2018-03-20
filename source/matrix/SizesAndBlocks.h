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

/** @file SizesAndBlocks.h

    \brief Class used to keep track of the block sizes used at
    different levels in the hierarchical matrix data structure.
*/

#ifndef MAT_SIZESANDBLOCKS
#define MAT_SIZESANDBLOCKS
#include "matInclude.h"
namespace mat{

  /** Describes dimensions of matrix and its blocks on all levels.
      The key ability is to provide the count and size of blocks, and
      their offset in the entire matrix.  It can generate a
      corresponding object for lower-level blocks. */
  class SizesAndBlocks { 
  public:
    /** Default constructor. */
    SizesAndBlocks()
      :nBlocks(0), nScalars(0), offset(0), nTotalScalars(0) {}
    /** Copy constructor. */
    SizesAndBlocks(SizesAndBlocks const & other);
    /** Constructor used for explicit calls. 
     *  For sizes and blocks at the highest level.
     *  nScalarsInp is the number of total scalar rows/columns in this case.
     */
  SizesAndBlocks(std::vector<int> const & blockSizesInp, 
		 int const nScalarsInp)
    : nBlocks(0), 
      nScalars(nScalarsInp), offset(0), nTotalScalars(nScalarsInp) {
      setup(blockSizesInp);
    }
    /** Assignment operator. */
    SizesAndBlocks& operator=
      (SizesAndBlocks const & other);
    
    bool operator==(SizesAndBlocks const & other) const;
    
    SizesAndBlocks 
      getSizesAndBlocksForLowerLevel(int const blockNumber) const;
    
    inline bool is_empty() const {return blockSizes.empty();}
    inline int const & getNBlocks() const {return nBlocks;}
    inline int const & getNScalars() const {return nScalars;}
    void getBlockSizeVector(std::vector<int> & blockSizesCopy) const;
    /** Returns the blocknumber (between 0 and nBlocks-1) 
     *  that contains elements with the given global index.
     *    
     */
    inline int whichBlock(int const globalIndex) const {
      return (globalIndex - offset) / blockSizes[0]; /* Integer division */
    }
    
    inline int getOffset() const {return offset;}
    inline int getNTotalScalars() const {return nTotalScalars;}
    ~SizesAndBlocks() {}
  protected:
    std::vector<int> blockSizes;
    /**< This is the number of scalars in each block,
     * (not the number of blocks in each block) for 
     * each level starting with the highest level.
     * It should be 1 at the lowest level.
     * Example: [1000 100 10 1]
     * Length is level() + 1
     */
    int nBlocks;          /**< This is the number of blocks in the current 
			   * block.
			   *
			   * == nScalars at lowest level */
      
    int nScalars;      /**< Number of scalars in the current block.  */
    int offset;        /**< Offset in entire system. */
    int nTotalScalars; /**< Total number of scalars in entire system. */
      
  SizesAndBlocks(std::vector<int> const & blockSizesInp, 
		 int const nScalarsInp,
		 int const offsetInp,
		 int const nTotalScalarsInp) 
    : nBlocks(0),
      nScalars(nScalarsInp), offset(offsetInp), 
      nTotalScalars(nTotalScalarsInp) {
	setup(blockSizesInp);
      }
    
    void setup(std::vector<int> const & blockSizesInp);
    
  private:
  }; /* end of class SizesAndBlocks */
  
} /* end namespace mat */
#endif
