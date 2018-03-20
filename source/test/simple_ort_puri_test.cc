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

#include <stdio.h>
#include <unistd.h>
#include <memory>
#include <limits>
#include "matrix_typedefs.h"
#include "matrix_utilities.h"
#include "utilities.h"

/** @file simple_ort_puri_test.cc Performs some simple tests of
    density matrix purification in orthogonal basis using artificialy
    generated input matrices. Written by Elias in Nov 2016. */

//const ergo_real penalty_alpha = 22;

static void report_timing(Util::TimeMeter & tm, const char* s) {
  double secondsTaken = tm.get_wall_seconds() - tm.get_start_time_wall_seconds();
  printf("report_timing for '%s': %f wall seconds\n", s, secondsTaken);
}


#if 0
static void get_tridiagonal_matrix_periodic(symmMatrix & F, int n, const std::vector<int> & perm) {
  // Create tridiagonal matrix with just the off-diagonal elements being nonzero
  std::vector<int> rows(n);
  std::vector<int> cols(n);
  std::vector<ergo_real> values(n);
  for(int i = 0; i < n; i++) {
    int row = i;
    int col = i+1;
    if(col == n)
      col = 0;
    rows[i] = row;
    cols[i] = col;
    values[i] = 1;
  }
  F.assign_from_sparse(rows, cols, values, perm, perm);
}

static void get_tridiagonal_matrix_nonperiodic(symmMatrix & F, int n, const std::vector<int> & perm) {
  // Create tridiagonal matrix with just the off-diagonal elements being nonzero
  int nElements = n-1;
  std::vector<int> rows(nElements);
  std::vector<int> cols(nElements);
  std::vector<ergo_real> values(nElements);
  for(int i = 0; i < nElements; i++) {
    rows[i] = i;
    cols[i] = i+1;
    values[i] = 1;
  }
  F.assign_from_sparse(rows, cols, values, perm, perm);
}
#endif


/*****************************************************************
  Here we construct Huckel-style Hamiltonian matrix for a periodic
  molecular system like this:

  H1          H3          H5          H7
  |           |           |           |
  C1 == C2 -- C3 == C4 -- C5 == C6 -- C7 == C8 -- 
        |           |           |           |
        H2          H4          H6          H8

  The basis functions are ordered like this:
  C1 H1 C2 H2 C3 H3 C4 H4 ...

  The Hamiltonian matrix then looks like this:

    C1 H1 C2 H2 C3 H3 C4 H4

C1  a  c  b
H1     d
C2        a  c  b
H2           d
C3              a  c  b
H3                 d
C4                    a  c
H4                       d

... and so on.

The values of the matrix elements a, b, c, d are as follows:
a = alpha (diagonal entry for C atoms)
b = beta (off-diagonal entry for neighboring C-C atoms)
d = alpha + h_A * beta (diagonal entry for H atoms)
c = k_AB * beta (off-diagonal entry for neighboring C-H atoms)

*****************************************************************/
static void get_Huckel_matrix_periodic(symmMatrix & F, int n, const std::vector<int> & perm) {
  assert(n % 2 == 0);
  assert(n >= 4);
  // FIXME: find out what are some reasonable values for the 4 parameters below. They affect gap and sparsity.
  // Setting alpha and h_A to zero still works, then the diagonal becomes zero and the ratio between beta and k_AB determines the properties of the system.
  const ergo_real alpha = 1.3;
  const ergo_real beta = 0.75;
  const ergo_real h_A = 0.95;
  const ergo_real k_AB = 0.15;
  // The number of matrix elements will be n*2 because we have n
  // elements on the diagonal, plus n off-diagonal elements.
  int nElements = 2*n;
  std::vector<int> rows(nElements);
  std::vector<int> cols(nElements);
  std::vector<ergo_real> values(nElements);
  int nCarbons = n/2;
  int count = 0;
  int row, col;
  ergo_real value;
  for(int i = 0; i < nCarbons; i++) {
    // Now do the row of C_i (3 matrix elements)
    row = i*2;
    col = row;
    value = alpha;
    rows[count] = row; cols[count] = col; values[count] = value; count++;
    col = row+1;
    if(col >= n)
      col -= n;
    value = k_AB * beta;
    rows[count] = row; cols[count] = col; values[count] = value; count++;
    col = row+2;
    if(col >= n)
      col -= n;
    value = beta;
    rows[count] = row; cols[count] = col; values[count] = value; count++;
    // Now do the row of H_i (1 matrix element)
    row++;
    col = row;
    value = alpha + h_A * beta;
    rows[count] = row; cols[count] = col; values[count] = value; count++;
  }
  assert(count == nElements);
  F.assign_from_sparse(rows, cols, values, perm, perm);
}

static void assign_from_full_matrix(symmMatrix & A, const ergo_real* A_full, int n, const std::vector<int> & perm) {
  int nElements = n*n;
  std::vector<int> rows(nElements);
  std::vector<int> cols(nElements);
  std::vector<ergo_real> values(nElements);
  ergo_real tolerance = template_blas_sqrt(mat::getMachineEpsilon<ergo_real>());
  int count = 0;
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) {
      ergo_real transpose_absdiff = template_blas_fabs(A_full[i*n+j] - A_full[j*n+i]);
      if(transpose_absdiff > tolerance) {
	std::cerr << "Error in assign_from_full_matrix: (transpose_absdiff > tolerance). transpose_absdiff = " << (double)transpose_absdiff << " tolerance = " << (double)tolerance << std::endl;
	throw std::runtime_error("Error in assign_from_full_matrix: (transpose_absdiff > tolerance).");
      }
      rows[count] = i;
      cols[count] = j;
      values[count] = A_full[i*n+j];
      count++;
    }
  assert(count == nElements);
  normalMatrix A2(A);
  A2.assign_from_sparse(rows, cols, values, perm, perm);
  A = A2;
}

static void print_matrix(const symmMatrix & A, const char* A_name, int n, ergo_real scaleFactor = 1.0) {
  printf("print_matrix for matrix '%s', n = %d\n", A_name, n);
  int nTrunc = 20;
  if(n > nTrunc)
    printf("Only showing part of matrix because (n > nTrunc).\n");
  if(scaleFactor != 1.0)
    printf("NOTE: USING scaleFactor = %f\n", (double)scaleFactor);
  normalMatrix A2(A);
  for(int row = 0; row < n; row++) {
    if(row > nTrunc)
      continue;
    std::vector<int> rows(n);
    std::vector<int> cols(n);
    for(int i = 0; i < n; i++) {
      rows[i] = row;
      cols[i] = i;
    }
    std::vector<ergo_real> values(n);
    A2.get_values(rows, cols, values);
    for(int col = 0; col < n; col++) {
      if(col > nTrunc)
	continue;
      printf(" %7.4f", (double)(values[col]*scaleFactor));
    }
    printf("\n");
  } // end for row
}

static void get_all_matrix_elements_nosymm(const normalMatrix & A, int n, std::vector<ergo_real> & result) {
  result.resize(n*n);
  std::vector<int> rows(n*n);
  std::vector<int> cols(n*n);
  int count = 0;
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) {
      rows[count] = i;
      cols[count] = j;
      count++;
    }
  A.get_values(rows, cols, result);
}

static void get_all_matrix_elements_symm(const symmMatrix & A, int n, std::vector<ergo_real> & result) {
  normalMatrix A2(A);
  get_all_matrix_elements_nosymm(A2, n, result);
}

struct DensMatInfo {
  ergo_real lambda_min;
  ergo_real lambda_max;
  ergo_real lambda_homo;
  ergo_real lambda_lumo;
};

static void get_density_mat_by_diagonalization(const symmMatrix & F,
					       symmMatrix & D,
					       DensMatInfo & info,
					       int n,
					       int n_occ,
					       const std::vector<int> & perm) {
  assert(n >= 1);
  assert(n_occ >= 0 && n_occ <= n);
  // Get full matrix
  std::vector<ergo_real> F_full(n*n);
  get_all_matrix_elements_symm(F, n, F_full);
  // Diagonalize
  int lwork = 3*n*n;
  std::vector<ergo_real> work(lwork);
  std::vector<ergo_real> eigvalList(n);
  std::vector<ergo_real> A(n*n);
  memcpy(&A[0], &F_full[0], n*n*sizeof(ergo_real));
  int syev_info = 0;
  mat::syev("V", "U", &n, &A[0],
	    &n, &eigvalList[0], &work[0], &lwork, 
	    &syev_info);
  if(syev_info != 0)
    throw std::runtime_error("ERROR: mat::syev() failed, gave (syev_info != 0).");
  // Save info about eigenvalues
  info.lambda_min = eigvalList[0];
  info.lambda_max = eigvalList[n-1];
  info.lambda_homo = eigvalList[n_occ-1];
  info.lambda_lumo = eigvalList[n_occ];
  // Put together density matrix using eigenvectors
  std::vector<ergo_real> D_full(n*n);
  for(int i = 0; i < n*n; i++)
    D_full[i] = 0;
  for(int i = 0; i < n_occ; i++) {
    // Add x*xT for cirrent eigenvector x
    const ergo_real* x = &A[i*n];
    for(int a = 0; a < n; a++)
      for(int b = 0; b < n; b++)
	D_full[a*n+b] += x[a]*x[b];
  }
  assign_from_full_matrix(D, &D_full[0], n, perm);
}

static void print_DensMatInfo(const DensMatInfo & info) {
  printf("info.lambda_min = %f\n", (double)info.lambda_min);
  printf("info.lambda_max = %f\n", (double)info.lambda_max);
  printf("info.lambda_homo = %f\n", (double)info.lambda_homo);
  printf("info.lambda_lumo = %f\n", (double)info.lambda_lumo);
}

static ergo_real get_nnz_percentage(int n, const symmMatrix & X) {
  ergo_real nnz = X.nnz();
  ergo_real n2 = (ergo_real)n*n;
  ergo_real nnz_percentage = 100 * nnz / n2;
  return nnz_percentage;
}

static void update_nnz_percentages(int n, const symmMatrix & X, ergo_real & nnz_percentage_min, ergo_real & nnz_percentage_max) {
  ergo_real nnz_percentage = get_nnz_percentage(n, X);
  if(nnz_percentage < nnz_percentage_min)
    nnz_percentage_min = nnz_percentage;
  if(nnz_percentage > nnz_percentage_max)
    nnz_percentage_max = nnz_percentage;
}


#if 0
static void verify_symmetry(int n, const normalMatrix & X) {
  std::vector<ergo_real> fullMat(n*n);
  get_all_matrix_elements_nosymm(X, n, fullMat);
  for(int i = 0; i < n; i++)
    for(int j = i+1; j < n; j++) {
      ergo_real absdiff = template_blas_fabs(fullMat[i*n+j]-fullMat[j*n+i]);
      assert(absdiff < 1e-9);
    }
}

static void get_symm_triple_product_ABA(int n,
					const symmMatrix & A,
					const symmMatrix & B,
					symmMatrix & result_ABA) {
  normalMatrix AB(A);
  AB = A * B;
  normalMatrix ABA(A);
  ABA = AB * A;
  // Check that ABA is really symmetric
  verify_symmetry(n, ABA);
  symmMatrix ABA_symm(ABA);
  result_ABA = ABA_symm;
}

static void get_triple_product_of_symm_mats(int n,
					    const symmMatrix & A,
					    const symmMatrix & B,
					    const symmMatrix & C,
					    normalMatrix & result_ABC) {
  normalMatrix AB(A);
  AB = A * B;
  normalMatrix ABC(A);
  ABC = AB * C;
  result_ABC = ABC;
}

static void get_triple_product_of_normal_mats(int n,
					      const normalMatrix & A,
					      const normalMatrix & B,
					      const normalMatrix & C,
					      normalMatrix & result_ABC) {
  normalMatrix AB(A);
  AB = A * B;
  normalMatrix ABC(A);
  ABC = AB * C;
  result_ABC = ABC;
}

static void get_diagonal_of_matrix(int n,
				   const symmMatrix & A,
				   symmMatrix & result_A_diag) {
  std::vector<int> rowind(n);
  std::vector<int> colind(n);
  for(int i = 0; i < n; i++) {
    rowind[i] = i;
    colind[i] = i;
  }
  std::vector<ergo_real> values(n);
  A.get_values(rowind, colind, values);
  result_A_diag.assign_from_sparse(rowind, colind, values);
}

static void add_to_one_element_symm(int row, int col, ergo_real x, symmMatrix & X, const std::vector<int> & perm) {
  std::vector<int> rowind;
  std::vector<int> colind;
  std::vector<ergo_real> values;
  if(row == col) {
    rowind.resize(1);
    colind.resize(1);
    values.resize(1);
    rowind[0] = row;
    colind[0] = col;
    values[0] = x;
  }
  else {
    // Two off-diagonal elements
    rowind.resize(2);
    colind.resize(2);
    values.resize(2);
    rowind[0] = row;
    colind[0] = col;
    values[0] = x;
    rowind[1] = col;
    colind[1] = row;
    values[1] = x;
  }
  normalMatrix tmp(X);
  tmp.assign_from_sparse(rowind, colind, values, perm, perm);
  symmMatrix tmp2(tmp);
  X += tmp2;
}

static void add_to_one_element_nosymm(int row, int col, ergo_real x, normalMatrix & X, const std::vector<int> & perm) {
  std::vector<int> rowind;
  std::vector<int> colind;
  std::vector<ergo_real> values;
  rowind.resize(1);
  colind.resize(1);
  values.resize(1);
  rowind[0] = row;
  colind[0] = col;
  values[0] = x;
  normalMatrix tmp(X);
  tmp.assign_from_sparse(rowind, colind, values, perm, perm);
  X += tmp;
}

static ergo_real evaluate_f1_symm(int n,
				  const symmMatrix & D,
				  const symmMatrix & E) {
  symmMatrix Po(D);
  symmMatrix Pv(D);
  Pv.add_identity(-1);
  Pv *= -1.0;
  // Get X = Po * E * Pv
  normalMatrix X(D);
  get_triple_product_of_symm_mats(n, Po, E, Pv, X);
  normalMatrix XTX(D);
  XTX = transpose(X) * X;
  return 2*XTX.trace();
}

static ergo_real evaluate_f1_nosymm(int n,
				    const symmMatrix & D,
				    const normalMatrix & E) {
  normalMatrix Po(D);
  normalMatrix Pv(D);
  Pv.add_identity(-1);
  Pv *= -1.0;
  normalMatrix ET(E);
  ET = transpose(E);
  // Get X = Po * E * Pv
  normalMatrix Xov(D);
  get_triple_product_of_normal_mats(n, Po, E, Pv, Xov);
  normalMatrix XovTXov(D);
  XovTXov = transpose(Xov) * Xov;
  normalMatrix Xvo(D);
  get_triple_product_of_normal_mats(n, Pv, E, Po, Xvo); // NOTE: it must be E here, not ET!
  normalMatrix XvoTXvo(D);
  XvoTXvo = transpose(Xvo) * Xvo;
#if 0
  // Check if Xvo is transpose of Xov (it is not!)
  normalMatrix T;
  T = transpose(Xov);
  T += (-1.0)*Xvo;
  std::cout << "T.frob() = " << T.frob() << "    Xvo.frob() = " << Xvo.frob() << "  Xov.frob() = " << Xov.frob() << std::endl;
#endif
  return XovTXov.trace() + XvoTXvo.trace();
}

static void zeroize_elements_inside_pattern(const symmMatrix & patternMat,
					    const normalMatrix & E_in,
					    normalMatrix & E_result,
					    const std::vector<int> & perm) {
  // Get all elements of patternMat in order to get the sparsity pattern
  normalMatrix patternMat_nosymm(patternMat);
  std::vector<int> rowind;
  std::vector<int> colind;
  std::vector<ergo_real> values;
  patternMat_nosymm.get_all_values(rowind, colind, values);
  // Get corresponding elements of E_in
  E_in.get_values(rowind, colind, values);
  normalMatrix tmp(E_in);
  tmp.assign_from_sparse(rowind, colind, values, perm, perm); // Now tmp contains the elements we want to remove.
  E_result = E_in;
  E_result -= tmp;
}

static void get_g_of_E_minus_E0(int n,
				const symmMatrix & E0,
				const symmMatrix & patternMat,
				const symmMatrix & E,
				const std::vector<int> & perm,
				symmMatrix & result_g_of_E_minus_E0) {
  symmMatrix E_minus_E0(E);
  E_minus_E0 -= E0;
  normalMatrix E_minus_E0_normal(E_minus_E0);
  normalMatrix g_of_E_minus_E0(E);
  zeroize_elements_inside_pattern(patternMat, E_minus_E0_normal, g_of_E_minus_E0, perm);
  verify_symmetry(n, g_of_E_minus_E0);
  symmMatrix g_of_E_minus_E0_symm(g_of_E_minus_E0);
  result_g_of_E_minus_E0 = g_of_E_minus_E0_symm;
}

static ergo_real evaluate_f2_nosymm(int n,
				    const symmMatrix & patternMat, // matrix defining sparsity pattern
				    const symmMatrix & E0,
				    const normalMatrix & E,
				    const std::vector<int> & perm) {
  normalMatrix E_minus_E0(E);
  normalMatrix E0_nosymm(E0);
  E_minus_E0 -= E0_nosymm;
  normalMatrix g_of_E_minus_E0(E);
  zeroize_elements_inside_pattern(patternMat, E_minus_E0, g_of_E_minus_E0, perm);
  normalMatrix Y(g_of_E_minus_E0);
  normalMatrix YTY(E);
  YTY = transpose(Y) * Y;
  ergo_real factor = 1;
  ergo_real result = factor * YTY.trace();
  //  printf("evaluate_f2_nosymm returning %12.6g  (Y.frob() = %12.6g)\n", result, Y.frob());
  return result;
}

static ergo_real evaluate_func_nosymm(int n,
				      const symmMatrix & D,
				      const symmMatrix & E0,
				      const symmMatrix & patternMat, // matrix defining sparsity pattern
				      const normalMatrix & E,
				      const std::vector<int> & perm,
				      const std::string & funcName) {
  if(funcName == "f1")
    return evaluate_f1_nosymm(n, D, E);
  else if(funcName == "f2")
    return evaluate_f2_nosymm(n, patternMat, E0, E, perm);
  else
    throw std::runtime_error("ERROR in evaluate_func_nosymm: unknown funcName.");
}

static ergo_real evaluate_f_complete_symm(int n,
					  const symmMatrix & D,
					  const symmMatrix & E0,
					  const symmMatrix & patternMat, // matrix defining sparsity pattern
					  const symmMatrix & E,
					  const std::vector<int> & perm) {
  normalMatrix E_nosymm(E);
  ergo_real f1Value = evaluate_f1_nosymm(n, D, E_nosymm);
  ergo_real f2Value = evaluate_f2_nosymm(n, patternMat, E0, E_nosymm, perm);
  return f1Value + penalty_alpha * f2Value;
}

static void get_numerical_gradient_of_func(int n,
					   const symmMatrix & D,
					   const symmMatrix & E0,
					   const symmMatrix & patternMat, // matrix defining sparsity pattern
					   const symmMatrix & E_in,
					   symmMatrix & result_numerical_gradient_of_f1,
					   const std::vector<int> & perm,
					   const std::string & funcName) {
  normalMatrix E_nosymm(E_in);
  // Get result numerically, as full matrix
  const ergo_real h = 1e-6;
  std::vector<ergo_real> M(n*n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) {
      normalMatrix E1(E_nosymm);
      add_to_one_element_nosymm(i, j, h, E1, perm);
      normalMatrix E2(E_nosymm);
      add_to_one_element_nosymm(i, j, -h, E2, perm);
      ergo_real f1_1 = evaluate_func_nosymm(n, D, E0, patternMat, E1, perm, funcName);
      ergo_real f1_2 = evaluate_func_nosymm(n, D, E0, patternMat, E2, perm, funcName);
      ergo_real gradientValue = (f1_1 - f1_2) / (2*h);
      if(funcName == "f2")
	printf("i j = %2d %2d  f2 gradientValue = %22.11f\n", i, j, (double)gradientValue);
      M[i*n+j] = gradientValue;
    }
  assign_from_full_matrix(result_numerical_gradient_of_f1, &M[0], n, perm);
}
#endif


static void do_truncation(int n,
			  symmMatrix & X,
			  ergo_real truncation_threshold,
			  const symmMatrix & D_in,
			  bool use_alt_trunc,
			  const std::vector<int> & perm) {
  if(use_alt_trunc == false) {
    X.frob_thresh(truncation_threshold);
    return;
  }
  // OK, use alternative truncation approach.

  symmMatrix A(X);
  A.frob_thresh(truncation_threshold);
  if(A.nnz() == X.nnz())
    return;
  symmMatrix zeroMat(A);
  zeroMat.clear();

  bool do_truncation_internally = true;

  // Use truncated D-matrix to make calculations faster.
  // We can truncate it a lot and still get reasonable results, but using only the diagonal will not work since then it commutes with any E-matrix.
  symmMatrix D(D_in);
  if(do_truncation_internally)
    D.frob_thresh(2.5);//truncation_threshold*100); // FIXME change amount of truncation of D here, seems like lots o truncation can be used while still giving OK results.

  // Now A contains the truncated X, that defines the sparsity pattern we want to use.
  // Get all elements of A in order to get the sparsity pattern
  std::vector<int> rowind;
  std::vector<int> colind;
  std::vector<ergo_real> values;
  A.get_all_values(rowind, colind, values);
  symmMatrix E(zeroMat);

  symmMatrix YY(zeroMat);
  YY.clear();

  normalMatrix ED(zeroMat);
  normalMatrix DE(zeroMat);
  normalMatrix E_new_3(zeroMat);
  normalMatrix DED(zeroMat);

  const int nOptSteps = 2; // NOTE: CHANGE nOptSteps HERE!
  for(int optStep = 0; optStep < nOptSteps; optStep++) {
    E = X - A; // this gives E such that X = A + E
    if(do_truncation_internally)
      E.frob_thresh(truncation_threshold*0.01); // FIXME change amount of truncation of E here

    // Compute gradient (?) matrix ED+DE-2*DED
    ED = E * D;
    if(do_truncation_internally)
      ED.frob_thresh(truncation_threshold*0.01); // FIXME choose amount of truncation here
    DED = D * ED;
    symmMatrix DED_symm(DED);
    DE = transpose(ED); // Same as D * E
    E_new_3 = ED + DE;
    symmMatrix E_new_symm(E_new_3);
    E_new_symm += (ergo_real)(-2.0) * DED_symm;
    // Now E_new_symm contains ED+DE-2*DED
    if(0) {
      symmMatrix analytical_gradient_of_f1(zeroMat);
      analytical_gradient_of_f1 = (ergo_real)2.0 * E_new_symm; // NOTE: multiply (ED+DE-2*DED) by factor 2 to get gradient
    }

    values.clear();
    E_new_symm.get_values(rowind, colind, values);
    YY.assign_from_sparse(rowind, colind, values);
    A += (ergo_real)2.0*YY; // FIXME: HOW TO SET FACTOR HERE? ARBITATY 0.3 VALUE!?
  } // end for optStep
  X = A;
  return;
}



#if 0
if(0) {
  symmMatrix analytical_gradient_of_f1(zeroMat);
  analytical_gradient_of_f1 = 2.0 * E_new_symm; // NOTE: multiply (ED+DE-2*DED) by factor 2 to get gradient
  symmMatrix numerical_gradient_of_f1(zeroMat);
  get_numerical_gradient_of_func(n, D, E0, A, E, numerical_gradient_of_f1, perm, "f1");
  print_matrix(numerical_gradient_of_f1, "numerical_gradient_of_f1", n, 1e2);
  print_matrix(analytical_gradient_of_f1, "analytical_gradient_of_f1", n, 1e2);
  ergo_real diff_f1 = symmMatrix::frob_diff(numerical_gradient_of_f1, analytical_gradient_of_f1);
  std::cout << "frobdiff between numerical and analytical gradients of f1: " << diff_f1 << "  (numerical_gradient_of_f1.frob() = " << numerical_gradient_of_f1.frob() << ")." << std::endl;
  symmMatrix g_of_E_minus_E0(zeroMat);
  get_g_of_E_minus_E0(n, E0, A, E, perm, g_of_E_minus_E0);
  symmMatrix analytical_gradient_of_f2(zeroMat);
  analytical_gradient_of_f2 = 2.0 * g_of_E_minus_E0; // NOTE multiply by 2 to get gradient
  symmMatrix numerical_gradient_of_f2(zeroMat);
  get_numerical_gradient_of_func(n, D, E0, A, E, numerical_gradient_of_f2, perm, "f2");
  print_matrix(numerical_gradient_of_f2, "numerical_gradient_of_f2", n, 1e5);
  print_matrix(analytical_gradient_of_f2, "analytical_gradient_of_f2", n, 1e5);
  ergo_real diff_f2 = symmMatrix::frob_diff(numerical_gradient_of_f2, analytical_gradient_of_f2);
  std::cout << "frobdiff between numerical and analytical gradients of f2: " << diff_f2 << "  (numerical_gradient_of_f2.frob() = " << numerical_gradient_of_f2.frob() << ")." << std::endl;
  if(analytical_gradient_of_f2.frob() > 1e-9)
    exit(0);
 }
symmMatrix analytical_gradient_of_f1(zeroMat);
analytical_gradient_of_f1 = 2.0 * E_new_symm; // NOTE: multiply (ED+DE-2*DED) by factor 2 to get gradient
symmMatrix g_of_E_minus_E0(zeroMat);
get_g_of_E_minus_E0(n, E0, A, E, perm, g_of_E_minus_E0);
symmMatrix analytical_gradient_of_f2(zeroMat);
analytical_gradient_of_f2 = 2.0 * g_of_E_minus_E0; // NOTE multiply by 2 to get gradient
symmMatrix gradient(analytical_gradient_of_f1);
gradient += penalty_alpha * analytical_gradient_of_f2;
E -= 0.004 * gradient;
#endif



static void report_subspace_error(int currIterCount, int n, const symmMatrix & X, const symmMatrix & D) {
  normalMatrix XD(X);
  XD = X * D;
  normalMatrix DX(X);
  DX = transpose(XD);
  normalMatrix diff(XD);
  diff -= DX;
  printf("report_subspace_error for currIterCount = %2d: diff.frob() = %g\n", currIterCount, (double)diff.frob());
}


static void report_subspace_error_via_diagonalization(int currIterCount, int n, int n_occ, const symmMatrix & X, const symmMatrix & D_ref, const std::vector<int> & perm) {
  symmMatrix D(X);
  D.clear();
  symmMatrix minusX(X);
  minusX *= -1.0;
  DensMatInfo info;
  get_density_mat_by_diagonalization(minusX, D, info, n, n_occ, perm);
  ergo_real diff = symmMatrix::frob_diff(D, D_ref);
  printf("report_subspace_error_via_diagonalization for currIterCount = %2d: diff = %g\n", currIterCount, (double)diff);
}


static void get_density_mat_by_purification(const symmMatrix & F,
					    symmMatrix & result_D,
					    int n,
					    int n_occ,
					    const std::vector<int> & perm,
					    const DensMatInfo & info,
					    ergo_real truncation_threshold,
					    const symmMatrix & D_ref,
					    bool use_alt_trunc,
					    bool verify_each_step) {
  printf("Starting get_density_mat_by_purification() now, truncation_threshold = %g, use_alt_trunc = %d.\n", (double)truncation_threshold, (int)use_alt_trunc);
  // FIXME: Truncate D_ref here?
  //  symmMatrix D_ref(D_ref_in);
  //  D_ref.frob_thresh(truncation_threshold);
  // Start by transformation to get all eigenvalues in [0,1] in reverse order
  symmMatrix X(F);
  X.add_identity(-info.lambda_max);
  X *= ((ergo_real)(-1.0) / (info.lambda_max - info.lambda_min));
  // Prepare X2, needed later
  symmMatrix X2(F);
  X2.clear();
  // Loop until converged
  const int maxIter = 222;
  int iterCount = 0;
  int extraFinalStepsCounter = 0;
  ergo_real X2_time_total = 0;
  ergo_real truncation_time_total = 0;
  ergo_real X_nnz_percentage_min = 100;
  ergo_real X_nnz_percentage_max = 0;
  ergo_real X2_nnz_percentage_min = 100;
  ergo_real X2_nnz_percentage_max = 0;
  while(true) {
    Util::TimeMeter tm_X2;
    X2 = (ergo_real)1.0 * X * X;
    X2_time_total += tm_X2.get_wall_seconds() - tm_X2.get_start_time_wall_seconds();
    update_nnz_percentages(n, X, X_nnz_percentage_min, X_nnz_percentage_max);
    update_nnz_percentages(n, X2, X2_nnz_percentage_min, X2_nnz_percentage_max);
    printf("Iteration %2d : nnz for X = %5.1f %% nnz for X2 = %5.1f %%\n", iterCount, (double)get_nnz_percentage(n, X), (double)get_nnz_percentage(n, X2));
    ergo_real X_X2_diff = symmMatrix::frob_diff(X, X2);
    if(X_X2_diff < 1e-1)
      extraFinalStepsCounter++;
    if(extraFinalStepsCounter == 7)
      break;
    ergo_real trace = X.trace();
    if(trace > n_occ) {
      // Choose polynomial X2
      X = X2;
    }
    else {
      // Choose polynomial 2*X - X2
      X *= 2.0;
      X += ((ergo_real)-1.0) * X2;
    }
    // Turn off use_alt_trunc in final iterations to allow eigenvalues to converge
    bool use_alt_trunc_effective = use_alt_trunc;
    if(iterCount > 12) // extraFinalStepsCounter > 4) // FIXME CHOOSE WHEN TO STOP USING ALT TRUNC HERE
      use_alt_trunc_effective = false;
    Util::TimeMeter tm_truncation;
    do_truncation(n, X, truncation_threshold, D_ref, use_alt_trunc_effective, perm);
    truncation_time_total += tm_truncation.get_wall_seconds() - tm_truncation.get_start_time_wall_seconds();
    if(verify_each_step) {
      report_subspace_error(iterCount, n, X, D_ref);
      report_subspace_error_via_diagonalization(iterCount, n, n_occ, X, D_ref, perm);
    }
    iterCount++;
    if(iterCount > maxIter)
      throw std::runtime_error("ERROR in get_density_mat_by_purification: (iterCount > maxIter).");
  }
  printf("Purification done, iterCount = %d, X2_time_total = %f seconds, truncation_time_total = %f seconds\n", iterCount, (double)X2_time_total, (double)truncation_time_total);
  printf("NNZ values for X : min %5.1f %% max %5.1f %%\n", (double)X_nnz_percentage_min, (double)X_nnz_percentage_max);
  printf("NNZ values for X2: min %5.1f %% max %5.1f %%\n", (double)X2_nnz_percentage_min, (double)X2_nnz_percentage_max);
  printf("Calling report_subspace_error_via_diagonalization for final result X matrix at end of get_density_mat_by_purification:\n");
  if(verify_each_step)
    report_subspace_error_via_diagonalization(iterCount, n, n_occ, X, D_ref, perm);
  result_D = X;
}

static void verify_idempotency(const symmMatrix & X) {
  symmMatrix X2(X);
  X2 = (ergo_real)1.0 * X * X;
  ergo_real diff = symmMatrix::frob_diff(X, X2);
  printf("verify_idempotency(): diff = %g\n", (double)diff);
  ergo_real tolerance = template_blas_sqrt(mat::getMachineEpsilon<ergo_real>());
  assert(diff < tolerance);
}

static void verify_gap(const DensMatInfo & info) {
  ergo_real gap_abs = info.lambda_lumo - info.lambda_homo;
  ergo_real spectrum_width = info.lambda_max - info.lambda_min;
  ergo_real gap_rel = gap_abs / spectrum_width;
  printf("verify_gap(): gap_abs = %f, gap_rel = %f\n", (double)gap_abs, (double)gap_rel);
  assert(gap_rel > 1e-6);
}


#if 0
static void check_decay_internal(int n,
				 const std::vector<ergo_real> & D_full,
				 int dist_step,
				 int dist_max) {
  // Consider distances from 0 up to dist_max
  for(int dist = 0; dist < dist_max; dist += dist_step) {
    ergo_real maxAbsValue = 0;
    for(int i = 0; i < n; i++)
      for(int j = i; j < n; j++) {
	int ijdiff  = j - i;
	if(ijdiff >= dist && ijdiff < n/2) {
	  ergo_real absValue = template_blas_fabs(D_full[i*n+j]);
	  if(absValue > maxAbsValue)
	    maxAbsValue = absValue;
	}
      } // end for i j
    printf("dist %5d : maxAbsValue = %12.8f\n", dist, (double)maxAbsValue);
  } // end for dist
}

static void check_decay(const symmMatrix & D, const char* name, int n) {
  printf("Starting check_decay() for '%s', n = %d\n", name, n);
  std::vector<ergo_real> D_full(n*n);
  get_all_matrix_elements_symm(D, n, D_full);
  // First part with step 1
  int dist_step = 1;
  int dist_max = 20;
  if(dist_max > n/2)
    dist_max = n/2;
  printf("==========  first part with dist_step = 1  ===============\n");
  check_decay_internal(n, D_full, dist_step, dist_max);
  dist_step *= 2;
  dist_max *= 2;
  printf("==========  first part with dist_step = 2  ===============\n");
  check_decay_internal(n, D_full, dist_step, dist_max);
  dist_step *= 2;
  dist_max *= 2;
  printf("==========  first part with dist_step = 4  ===============\n");
  check_decay_internal(n, D_full, dist_step, dist_max);
  dist_step = n / 40;
  if(dist_step == 0)
    dist_step = 1;
  dist_max = n/2;
  printf("==========  now with dist_step = %d  ===============\n", dist_step);
  check_decay_internal(n, D_full, dist_step, dist_max);
  printf("check_decay() done.\n");
}
#endif


/*
ELIAS NOTE 2016-11-04: the following input seems to work for blockSize = 1 now:
./simple_ort_puri_test 14 1e-2 1e-1 1
*/

int main(int argc, char *argv[])
{
  try {
    Util::TimeMeter tm_everything;

    // OK n values: 10 14 18 22 26 30 50 70 90 110 150 210 290 410 610 1010 1410 2010 2410 3010
    int n = 8;
    ergo_real truncation_threshold = 1e-5;
    ergo_real result_diff_tolerance = 1e-4;
    int use_alt_trunc = 0;
    int verify_each_step = 0;

#ifdef PRECISION_SINGLE
    result_diff_tolerance = 1e-3;
#endif

    if(argc >= 2)
      n = atoi(argv[1]);
    if(argc >= 3)
      truncation_threshold = atof(argv[2]);
    if(argc >= 4)
      result_diff_tolerance = atof(argv[3]);
    if(argc >= 5)
      use_alt_trunc = atoi(argv[4]);
    if(argc >= 6)
      verify_each_step = atoi(argv[5]);

    int n_occ = n/2;
    int blockSize = 8;

    std::cout << "n = " << n << std::endl;
    std::cout << "truncation_threshold = " << (double)truncation_threshold << std::endl;
    std::cout << "result_diff_tolerance = " << (double)result_diff_tolerance << std::endl;
    std::cout << "n_occ = " << n_occ << std::endl;
    std::cout << "blockSize = " << blockSize << std::endl;

    if(n <= 1)
      throw std::runtime_error("ERROR: (n <= 1)");

#ifdef _OPENMP
    int defThreads;
    const char *env = getenv("OMP_NUM_THREADS");
    if ( !(env && (defThreads=atoi(env)) > 0) )
      defThreads = 1;
    mat::Params::setNProcs(defThreads);
    mat::Params::setMatrixParallelLevel(1);
    std::cout<<"OpenMP is used, number of threads set to "
	     <<mat::Params::getNProcs()<<". Matrix parallel level: "
	     <<mat::Params::getMatrixParallelLevel()<<"."<<std::endl;
#endif

    // Prepare stuff needed to use matrix library
    mat::SizesAndBlocks sizeBlockInfo;
    static const int sparseMatrixBlockFactor = 2;
    sizeBlockInfo =
      prepareMatrixSizesAndBlocks(n,
				  blockSize,
				  sparseMatrixBlockFactor,
				  sparseMatrixBlockFactor,
				  sparseMatrixBlockFactor);
    std::vector<int> perm(n);
    for(int i = 0; i < n; i++)
      perm[i] = i;

    // Create artificial effective Hamiltonian matrix F
    symmMatrix F;
    F.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
    // get_tridiagonal_matrix_periodic(F, n, perm);
    get_Huckel_matrix_periodic(F, n, perm);
    print_matrix(F, "F", n);

    // Get reference density matrix using diagonalization
    symmMatrix D_ref;
    D_ref.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
    DensMatInfo densMatInfo;
    Util::TimeMeter tm_get_density_mat_by_diagonalization;
    get_density_mat_by_diagonalization(F, D_ref, densMatInfo, n, n_occ, perm);
    report_timing(tm_get_density_mat_by_diagonalization, "tm_get_density_mat_by_diagonalization");
    print_matrix(D_ref, "D_ref", n);
    print_DensMatInfo(densMatInfo);

    // Verify that D_ref is idempotent
    verify_idempotency(D_ref);

    // Verify that we have significant gap
    verify_gap(densMatInfo);

    // Get density matrix using purification
    symmMatrix D_puri;
    D_puri.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
    Util::TimeMeter tm_get_density_mat_by_purification;
    get_density_mat_by_purification(F, D_puri, n, n_occ, perm, densMatInfo, truncation_threshold, D_ref, use_alt_trunc, verify_each_step);
    report_timing(tm_get_density_mat_by_purification, "tm_get_density_mat_by_purification");

    // Check accuracy of result
    ergo_real diff_D_puri_vs_D_ref = symmMatrix::frob_diff(D_puri, D_ref);
    std::cout << "diff_D_puri_vs_D_ref = " << (double)diff_D_puri_vs_D_ref << std::endl;
    symmMatrix D_ref_truncated(D_ref);
    D_ref_truncated.frob_thresh(truncation_threshold);
    ergo_real diff_D_ref_vs_D_ref_truncated = symmMatrix::frob_diff(D_ref, D_ref_truncated);
    std::cout << "diff_D_ref_vs_D_ref_truncated = " << (double)diff_D_ref_vs_D_ref_truncated << std::endl;
    assert(diff_D_puri_vs_D_ref < result_diff_tolerance);

    // Investigate decay of matrix elements in D
    //    check_decay(D_ref, "D_ref", n);
    //    check_decay(D_puri, "D_puri", n);

    report_timing(tm_everything, "tm_everything");
  }
  catch(std::runtime_error & e) {
    std::cout << "Error: std::runtime_error caught: " << e.what() << std::endl;
    return -1;
  }
  catch(std::exception & e) {
    std::cout << "Error: std::exception caught: " << e.what() << std::endl;
    return -1;
  }
  catch(...) {
    printf("ERROR: exception caught in simple_ort_puri_test main().\n");
    return -1;
  }

  puts("simple_ort_puri_test finished OK.");
  return 0;
}
