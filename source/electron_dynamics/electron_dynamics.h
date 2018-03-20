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

/** @file electron_dynamics.h

    @brief Parameters related to TDHF electron dynamics.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef ELECTRON_DYNAMICS_HEADER
#define ELECTRON_DYNAMICS_HEADER

#include <string>

#include "realtype.h"

namespace ED {

struct Params
{
  ergo_real max_time;
  ergo_real timestep;
  ergo_real dc_pulse_strength;
  ergo_real dc_pulse_time;
  ergo_real ac_pulse_max;
  ergo_real ac_pulse_omega;
  std::string field_type;

  Params() 
  : max_time(125),
    timestep(0.2),
    dc_pulse_strength(0.005),
    dc_pulse_time(10),
    ac_pulse_max(0.07),
    ac_pulse_omega(0.1),
    field_type("none")
  {}
    
};

} /* End of namespace ED.  */

#endif
