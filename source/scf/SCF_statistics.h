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

/** @file SCF_statistics.h

    @brief Class for keeping timings and other statistics related to
    self-consistent field (SCF) procedure.

    @author: Emanuel Rubensson <em>responsible</em>.
*/

#ifndef SCF_STATISTICS_HEADER
#define SCF_STATISTICS_HEADER
#include <map>
#include "utilities.h"

struct SCF_timer {
  SCF_timer();
  void stop();
  double elapsedTimeCPU_sys;
  double elapsedTimeCPU_usr;
  double elapsedTimeWall;  
private:
  double startTimeCPU_sys;
  double startTimeCPU_usr;
  double startTimeWall;
  bool stopped_already;
};


class SCF_statistics {
  typedef std::map<std::string, SCF_timer> TimerMap;  
  typedef std::map<std::string, double> ValueMap;  
 public:  
  void start_timer(std::string identifier);
  void stop_timer(std::string identifier);
  void add_value(std::string identifier, double value);
  void add_values( ValueMap & values_to_add);
  void output_mfile(std::string name);
 protected:
  TimerMap timers;
  ValueMap values;  
 private:
  void output_value( std::ofstream & os, std::string id, double value);

};




#endif

