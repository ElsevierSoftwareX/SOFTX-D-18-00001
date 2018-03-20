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

/** @file DebugPolicies.h Classes describing debug policies with different 
 * debug levels. Choice of a higher level gives more tests that the program 
 * executes as expected but at a higher memory and time cost.
 * Normal execution runs at "DebugLevelLow". 
 *
 * Copyright(c) Emanuel Rubensson 2007
 *
 * @author Emanuel Rubensson 
 * @date January 2007
 *
 *
 *
 *
 */
#ifndef MAT_DEBUGPOLICIES
#define MAT_DEBUGPOLICIES

#include <cstdlib>

namespace mat{
#if 0
#define ASSERTALWAYS(x)						\
  this->assertAlways(__FILE__, __LINE__, __DATE__, __TIME__,x)
#define ASSERTDEBUG(x)						\
  this->assertDebug(__FILE__, __LINE__, __DATE__, __TIME__,x)
  /* debugPolicies */

  class DebugLevelHigh {
  public:
    void assertAlways(char const * theFile, int const theLine, 
		      char const * theDate, char const * theTime,
		      bool const statement) const {
      if (!statement) {
	std::cout<<"Assertion failed: "<<theFile<<":"<<theLine
		 <<" Compiled on "<<theDate<<" at "<<theTime<<".\n";
	std::exit(1);
      }
    }
    inline void assertDebug(char const * theFile, int const theLine, 
			    char const * theDate, char const * theTime,
			    bool const statement) const {
      assertAlways(theFile, theLine, theDate, theTime, statement);
    }
  };
  class DebugLevelMedium : public DebugLevelHigh {};
  class DebugLevelLow : public DebugLevelMedium  {
  public:
    inline void assertDebug(char const * theFile, int const theLine, 
			    char const * theDate, char const * theTime,
			    bool const statement) const {}
  };

#else
  

#define ASSERTALWAYS(x)						\
  this->assertAlways(__FILE__, __LINE__, __ID__,x)
#define ASSERTDEBUG(x)						\
  this->assertDebug(__FILE__, __LINE__, __ID__,x)
#endif
  /* debugPolicies */
  class DebugLevelHigh {
  public:
    void assertAlways(char const * theFile, int const theLine, 
		      char const * theId, bool const statement) const {
      if (!statement) {
	std::cout<<"Assertion failed: "<<theFile<<":"<<theLine
		 <<" svn info: "<<theId<<".\n";
	std::exit(1);
      }
    }
    inline void assertDebug(char const * theFile, int const theLine, 
			    char const * theId, bool const statement) const {
      assertAlways(theFile, theLine, theId, statement);
    }
  };
  class DebugLevelMedium : public DebugLevelHigh {};
  class DebugLevelLow : public DebugLevelMedium  {
  public:
    inline void assertDebug(char const * theFile, int const theLine, 
			    char const * theId, bool const statement) const {}
  };


  
} /* end namespace mat */
#endif
