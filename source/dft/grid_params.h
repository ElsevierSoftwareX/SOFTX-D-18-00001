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

/** @file grid_params.h

    @brief Code organizing various settings related to grids used for
    numerical integration.

    @author: Pawel Salek <em>responsible</em>
*/

#if !defined(_GRID_PARAMS_H_)
#define _GRID_PARAMS_H_ 1

#include "grid_atomic.h"

namespace Dft {

/** A structure describing the HiCu grid settings. */
struct HiCuGridParams {
  ergo_real maxError;
  ergo_real box_size;
  ergo_real start_box_size_debug;
  int use_error_per_volume;
  int do_double_checking;
  int compare_to_refined;
  int use_energy_criterion;
  int use_energy_criterion_only;
  int do_variation_checking;
};

/** A structure describing the grid settings. */
struct GridParams {
  /** All the dimensions of the smallest box must be below this
      threshold.  The time goes quickly up as a function of box
      size. Tweak it with an uttermost caution. */
  ergo_real boxSize;
  ergo_real radint;
  int angmin;
  int angmax;
  typedef enum { GC2, LMG, TURBO } RadialScheme;
  typedef enum { TYPE_STANDARD, TYPE_HICU } GridType;
  RadialScheme radialGridScheme;
  GridType  gridType;
  bool cubicBoxes; /**< whether cubic grid boxes should be enforced.
		      Not needed apart from testing. */
  /* The following are HiCu grid parameters. */
  HiCuGridParams hicuParams;
explicit GridParams(ergo_real r_ = 1e-9, int a1 = 6, int a2 = 30,
	   ergo_real bs = 5.0, bool cubic = false,
	   ergo_real hicume = 1e-7, 
	   ergo_real hicubs = 1.5, ergo_real hicusbsd = 0, 
	   int hicuerrpervol = 0,
	   int hicudodoublecheck = 1,
	   int hicuctr = 0, int hicuuec = 0,int hicuueco = 0,
	   int hicudovarcheck = 0)
: boxSize(bs), radint(r_), angmin(a1), angmax(a2), radialGridScheme(LMG),
    gridType(TYPE_STANDARD), cubicBoxes(cubic)
  { 
    hicuParams.maxError = hicume;
    hicuParams.box_size = hicubs;
    hicuParams.start_box_size_debug = hicusbsd;
    hicuParams.use_error_per_volume = hicuerrpervol; 
    hicuParams.do_double_checking = hicudodoublecheck;
    hicuParams.compare_to_refined = hicuctr;
    hicuParams.use_energy_criterion = hicuuec;
    hicuParams.use_energy_criterion_only = hicuueco;
    hicuParams.do_variation_checking = hicudovarcheck;
  }
};

}
#endif /* _GRID_PARAMS_H_ */
