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

/** @file SizesAndBlocks.cc

    \brief Class used to keep track of the block sizes used at
    different levels in the hierarchical matrix data structure.
*/

#include "SizesAndBlocks.h"
namespace mat{
  
  SizesAndBlocks::SizesAndBlocks(SizesAndBlocks const & other) 
    : blockSizes(other.blockSizes), 
      nBlocks(other.nBlocks),
      nScalars(other.nScalars), 
      offset(other.offset), 
      nTotalScalars(other.nTotalScalars) {}
  
  SizesAndBlocks& SizesAndBlocks::operator=
  (SizesAndBlocks const & other) {
    nBlocks = other.nBlocks;
    nScalars = other.nScalars; 
    offset = other.offset; 
    nTotalScalars = other.nTotalScalars;
    blockSizes = other.blockSizes;
    return *this;
  }

  bool SizesAndBlocks::operator==(SizesAndBlocks const & other) const {
    bool isEqual = 
      (blockSizes.size() == other.blockSizes.size()) && 
      (nBlocks == other.nBlocks) && 
      (nScalars == other.nScalars) &&
      (offset == other.offset) &&
      (nTotalScalars == other.nTotalScalars); 
    if (isEqual)
      for (unsigned int i = 0; i < blockSizes.size(); i++)
	isEqual = isEqual && 
	  (blockSizes[i] == other.blockSizes[i]);
    return isEqual;
  }

  SizesAndBlocks SizesAndBlocks::
  getSizesAndBlocksForLowerLevel(int const blockNumber) const {
    assert(blockSizes.size() > 1);
    int nScalLowLev;
    if ((blockNumber+1) * blockSizes[0] > nScalars)
      nScalLowLev = nScalars - blockNumber * blockSizes[0];
    else
      nScalLowLev = blockSizes[0];
    std::vector<int> nextBlockSizes(blockSizes.begin() + 1, blockSizes.end());
    assert(offset + blockNumber * blockSizes[0] + nScalLowLev <= nTotalScalars);
    return SizesAndBlocks(nextBlockSizes, 
			  nScalLowLev,
			  offset + blockNumber * blockSizes[0],
			  nTotalScalars);
  }

  void SizesAndBlocks::
  getBlockSizeVector(std::vector<int> & blockSizesCopy) const {
    blockSizesCopy = blockSizes;
  }


  void SizesAndBlocks::setup(std::vector<int> const & blockSizesInp) {
    blockSizes = blockSizesInp;
    assert(!blockSizes.empty());
    assert(blockSizes[blockSizes.size()-1]);
    for (unsigned int ind = 0; ind < blockSizes.size()-1; ind++)
      assert(blockSizes[ind] >= blockSizes[ind+1]);
    nBlocks = nScalars/blockSizes[0]; /* Integer division. */
    if (nScalars%blockSizes[0]) 
      ++nBlocks;
  }

} /* end namespace mat */
