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

/** @file Failure.h

    \brief The Failure class is used for exception handling. It
    inherits std::exception.

    @author Emanuel Rubensson
*/

/* The Failure class is used for exception handling.                *
 * It inherits std::exception which means that an instance of this  *
 * class can be caught by catching std::exception&                  *
 * The "&" must be there, otherwise the failure message             *
 * will be "cut out".                                               *
 * Retrieve the message with a call to the member function "what()" *  
 *                                                                  *
 *                                                                  *
 *     \\\|||///  \    (C) Emanuel Rubensson, August, 2005          *
 *     \ ~   ~ /   \                                                *
 *     | @   @ |    \  mail:  emanuel@theochem.kth.se               *
 * oOo---(_)---oOo---\----------------------------------------------*
 *                                                                  */

#ifndef FAILURE
#define FAILURE
#include <exception>
namespace mat
{
  class Failure : public std::exception {
      const char* message;
    public:
      Failure()
	:message("Failure: No failure information available")
	{}
      explicit Failure (const char* msg)
	:message(msg)
	{}
      virtual ~Failure() throw() {}
      virtual const char* what() const throw()
	{return message;}
  };

  class Acceptable : public Failure {
  public:
    Acceptable()
      :Failure("Acceptable: No acceptable failure information available")
      {}
    explicit Acceptable (const char* msg)
      :Failure(msg)
      {}
  };
  
  class AcceptableMaxIter : public Acceptable {
    int maxiter;
  public:
    AcceptableMaxIter()
      :Acceptable("AcceptableMaxIter: No acceptable failure information available"),
      maxiter(0)
      {}
    explicit AcceptableMaxIter (const char* msg, int maxi = 0)
      :Acceptable(msg), maxiter(maxi)
      {}
    int get_maxiter() const {
      return maxiter;
    }
  };
  
} /* end namespace */
#endif
