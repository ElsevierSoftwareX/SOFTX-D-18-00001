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

/** @file puri_info.h

    @brief File containing classes IterationInfo and PuriInfo.
    IterationInfo is a class with the information stored for each iteration of the recursive expansion.
    PuriInfo is a class containing general information about the recursive expansion.

    @author Anastasia Kruchinina <em>responsible</em>
    @sa puri_info.cc
*/



#ifndef PURI_INFO_HEADER
#define PURI_INFO_HEADER

#include "output.h"
#include "matrix_typedefs.h"  // definitions of matrix types and interval type
#include "realtype.h"         // definitions of types

//#define CHECK_IF_STOPPED_TOO_LATE_OR_TOO_EARLY


class IterationInfo{
 public:
  typedef ergo_real real;

  int it; //iteration number
  real threshold_X;
  real Xsquare_time;
  real trunc_time;
  real purify_time;
  real total_time;
  real stopping_criterion_time;
  real eucl_diff_time;
  real trace_diff_time ;
  real mixed_diff_time;
  real frob_diff_time;
  real orbital_homo_time;
  real orbital_lumo_time;
  real DX_mult_homo_time;
  real DX_mult_lumo_time;
  real homo_eig_solver_time;
  real lumo_eig_solver_time;
  real XmX2_trace;
  real XmX2_fro_norm;
  real XmX2_mixed_norm;
  real XmX2_eucl;
  real order;
  int poly;
  real gap; // estimated gap
  real NNZ_X;
  real NNZ_X2;
  // bounds for homo and lumo during iterations
  // [lumo_low, lumo_upp] and [1-homo_upp, 1-homo_low]
  real homo_bound_low;
  real homo_bound_upp;
  real lumo_bound_low;
  real lumo_bound_upp;

  real commutation_error;

  real alpha; // for SP2 accelerated

  real constantC;


 IterationInfo():
  it(-1),
    threshold_X(0),
    Xsquare_time(0),
    trunc_time(0),
    purify_time(0),
    total_time(0),
    stopping_criterion_time(0),
    eucl_diff_time(0),
    trace_diff_time(0),
    mixed_diff_time(0),
    frob_diff_time(0),
    orbital_homo_time(0),
    orbital_lumo_time(0),
    DX_mult_homo_time(0),
    DX_mult_lumo_time(0),
    homo_eig_solver_time(0),
    lumo_eig_solver_time(0),
    XmX2_trace(-1),
    XmX2_fro_norm(-1),
    XmX2_mixed_norm(-1),
    XmX2_eucl(-1),
    order(0),
    poly(-1),
    gap(-1),
    NNZ_X(0),
    NNZ_X2(0),
    homo_bound_low(0),
    homo_bound_upp(0),
    lumo_bound_low(0),
    lumo_bound_upp(0),
    commutation_error(0),
    alpha(0),
    constantC(0)
      {};


};


class PuriInfo{
 public:
  typedef ergo_real real;

 PuriInfo() :
  method(0),
    stopping_criterion(0),
    total_it(0),
    converged(0),
    error_subspace(0),
    accumulated_error_subspace(0),
    compute_eigenvectors_in_this_SCF_cycle(false),
    homo_eigenvector_is_computed(false),
    lumo_eigenvector_is_computed(false),
    homo_eigenvector_is_computed_in_iter(-1),
    lumo_eigenvector_is_computed_in_iter(-1),
    homo_eigensolver_iter(-1),
    lumo_eigensolver_iter(-1),
    homo_eigensolver_time(-1),
    lumo_eigensolver_time(-1),
    debug_output(0)
      {};


  void print_collected_info();
  void print_collected_info_printf();

  void get_poly_seq(std::vector<int> &norms);
  void get_vec_frob_norms(std::vector<real> &norms);
  void get_vec_mixed_norms(std::vector<real> &norms);
  void get_vec_traces(std::vector<real> & traces);

  void get_spectrum_bounds(real &lower_spectrum_bound_, real &upper_spectrum_bound_) const;
  void set_spectrum_bounds(const real lower_spectrum_bound_, const real upper_spectrum_bound_);

  int method; // 1 for SP2, 2 for SP2 accelerated

  int stopping_criterion; // 1 if new, 0 if not
  real norm_F_Fprev;
  real total_time;
  int  total_it;
  real time_spectrum_bounds;
  int  estim_total_it;
  int  additional_iterations;

  int converged; // 1 if converged, 0 otherwise

  real error_subspace; // expected maximum error in subspace
  real accumulated_error_subspace; // accumulated error in subspace

  real homo_estim_upp_F;
  real homo_estim_low_F;
  real lumo_estim_upp_F;
  real lumo_estim_low_F;

  bool compute_eigenvectors_in_this_SCF_cycle;
  bool homo_eigenvector_is_computed;
  bool lumo_eigenvector_is_computed;
  int homo_eigenvector_is_computed_in_iter;
  int lumo_eigenvector_is_computed_in_iter;
  int homo_eigensolver_iter;
  int lumo_eigensolver_iter;
  double homo_eigensolver_time;
  double lumo_eigensolver_time;
  real eigValHOMO;
  real eigValLUMO;

  std::vector<IterationInfo> Iterations;
  int debug_output;

  real upper_spectrum_bound;
  real lower_spectrum_bound;

};

#endif
