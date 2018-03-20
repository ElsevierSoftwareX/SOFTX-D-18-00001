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

/** @file SCF_statistics.cc

    @brief Class for keeping timings and other statistics related to
    self-consistent field (SCF) procedure.

    @author: Emanuel Rubensson <em>responsible</em>.
*/

#include <stdexcept>
#include <fstream>
#include "SCF_statistics.h"

SCF_timer::SCF_timer()
  : stopped_already(false) {
  startTimeWall = Util::TimeMeter::get_wall_seconds();
  Util::TimeMeter::get_current_cpu_times(startTimeCPU_usr, startTimeCPU_sys);
}

void SCF_timer::stop() {
  if (stopped_already)
    throw std::runtime_error("Attempt to stop timer already stopped.");
  elapsedTimeWall = Util::TimeMeter::get_wall_seconds() - startTimeWall;
  double stopTimeCPU_sys, stopTimeCPU_usr;
  Util::TimeMeter::get_current_cpu_times(stopTimeCPU_usr, stopTimeCPU_sys);
  elapsedTimeCPU_sys = stopTimeCPU_sys - startTimeCPU_sys;
  elapsedTimeCPU_usr = stopTimeCPU_usr - startTimeCPU_usr;    
  stopped_already = true;
}

void SCF_statistics::start_timer(std::string identifier) {
  timers[identifier] = SCF_timer();
}
void SCF_statistics::stop_timer(std::string identifier) {
  if ( timers.find(identifier) == timers.end() )
    throw std::runtime_error("Attempt to stop timer not in timer map.");
  timers[identifier].stop();
}

void SCF_statistics::add_value(std::string identifier, double value) {
  if ( values.find(identifier) != values.end() )
    throw std::runtime_error("Attempt to add value already in value map.");
  values[identifier] = value;
}

void  SCF_statistics::add_values( ValueMap & values_to_add) {
  ValueMap::const_iterator it;
  for ( it=values_to_add.begin() ; it != values_to_add.end(); it++ ) 
    values[it->first] = it->second;
}


void SCF_statistics::output_mfile(std::string name) {
  std::string m_name = name + ".m";
  std::ofstream os(m_name.c_str());
  
  // First output the names of all variables as one big comment
  os << "%% SCF_statistics, list of all variables " << std::endl;
  {
    os << "%% Timers: " << std::endl;
    TimerMap::const_iterator it;
    for ( it=timers.begin() ; it != timers.end(); it++ ) {
      std::string s = (*it).first;
      std::string s_wall = s + "_walltime"; 
      std::string s_cpu_sys = s + "_cpu_sys"; 
      std::string s_cpu_usr = s + "_cpu_usr"; 
      os << "% " << s_wall    << std::endl;
      os << "% " << s_cpu_sys << std::endl;
      os << "% " << s_cpu_usr << std::endl;
    }
  }
  {
    os << "%% Other values: " << std::endl;
    ValueMap::const_iterator it;
    for ( it=values.begin() ; it != values.end(); it++ ) {
      std::string s = (*it).first;
      os << "% " << s << std::endl;
    }
  }
  os << "%" << std::endl << std::endl;
  
  // Now output the values
  os << "%% SCF_statistics timers " << std::endl;
  {
    TimerMap::const_iterator it;
    for ( it=timers.begin() ; it != timers.end(); it++ ) {
      std::string s = (*it).first;
      std::string s_wall = s + "_walltime"; 
      std::string s_cpu_sys = s + "_cpu_sys"; 
      std::string s_cpu_usr = s + "_cpu_usr"; 
      double time_cpu_sys = (*it).second.elapsedTimeCPU_sys;
      double time_cpu_usr = (*it).second.elapsedTimeCPU_usr;
      double time_wall = (*it).second.elapsedTimeWall;
      output_value( os, s_wall   , time_wall    );
      output_value( os, s_cpu_sys, time_cpu_sys );
      output_value( os, s_cpu_usr, time_cpu_usr );
    }
  }
  os << "%% SCF_statistics other values " << std::endl;
  {
    ValueMap::const_iterator it;
    for ( it=values.begin() ; it != values.end(); it++ ) {
      std::string s = (*it).first;
      double value = (*it).second;
      output_value( os, s, value );
    }
  }
}

void SCF_statistics::output_value( std::ofstream & os, 
				   std::string id, 
				   double value ) {
  os << "if ( ~exist( '" << id << "' ) )" << std::endl;
  os << "  " << id << " = [];" << std::endl;
  os << "end" << std::endl;
  os << id << " = [" << id << " " << value << "];" 
     << std::endl;  
}
