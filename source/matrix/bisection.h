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

/** @file bisection.h Bisection method
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date May 7, 2006
 *
 */
#ifndef MAT_BISECTION
#define MAT_BISECTION
#include <cmath>
namespace mat {
  /** Sign function returns the sign of the input. 
   *
   * 1  for positive, 
   * -1 for negative and
   * 0  for zero.
   *
   */
  template<typename Treal>
    inline int sign(Treal value) {
    if (value > 0)
      return 1;
    else if (value < 0)
      return -1;
    else
      return 0;
  }


  /** Bisection algorithm for root finding
   *
   * The bisection method finds the root of a function in the interval
   * [min, max], or more precisely the place where the function changes sign. 
   * It is assumed that the function only changes sign once in the 
   * given interval.
   * 
   * The function is given by a class that has a member function named eval
   * that evaluates the function in the given point.
   *
   */
  template<typename Treal, typename Tfun>
    Treal bisection(Tfun const & fun, Treal min, Treal max, Treal const tol) {
    int sign_min = sign(fun.eval(min));
    int sign_max = sign(fun.eval(max));
    if (sign_min == sign_max)
      throw Failure("bisection(Tfun&, Treal, Treal, Treal): interval "
		    "incorrect");
    Treal middle = (max + min) / 2;
    int sign_middle = sign(fun.eval(middle));
    while (template_blas_fabs(max - min) > tol * 2 && sign_middle != 0) {
      if (sign_middle == sign_min) {
	min = middle;
	sign_min = sign_middle;
      }
      else { /* (sign_middle == sign_max) */
	max = middle;
	sign_max = sign_middle;
      }
      middle = (max + min) / 2;
      sign_middle = sign(fun.eval(middle));
    }
    return middle;
  }

} /* end namespace mat */
#endif
