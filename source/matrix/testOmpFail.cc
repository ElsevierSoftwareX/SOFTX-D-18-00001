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

/** @file testOmpFail.cc

    \brief Test code designed to detect a certain kind of compiler bug
    that appeared when certain compilers are used with OpenMP.
*/

#include <iostream>
#include <cstdlib>

#ifdef _OPENMP
#include <omp.h>
#endif

int main() {
  std::cout<<" Testing code sensitive when certain compilers "
    "are used with OpenMP: \n";
  std::cout<<"  OpenMP is used?   ";
#ifdef _OPENMP
  std::cout<<"YES"<<std::endl;
#else
  std::cout<<"NO"<<std::endl;
#endif
  
  int count = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(2)
#endif
  for (int i = 0; i < 1; i++) {
    try {
      count++;
      for (int j = 0; j < 0; j++) {}
      if (0) throw 0;
    } catch (...) {}
  }
  if (count == 1) {
    std::cout<<"  Test OK" <<std::endl;
    std::exit(0);
  }
  else {
    std::cout<<"  Test Failed" <<std::endl;
    std::exit(1);
  }
}



