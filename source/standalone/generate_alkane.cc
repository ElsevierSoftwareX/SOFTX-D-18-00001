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

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <string>
#include <cstring>
#include <vector>
#include <cassert>

typedef struct {
  std::string label;
  double coord[3];
} atom_struct;

struct vec3d_struct {
  double coord[3];
  void set(double x, double y, double z) {
    coord[0] = x;
    coord[1] = y;
    coord[2] = z;
  }
  vec3d_struct(double x, double y, double z) {
    set(x, y, z);
  }
};

int write_xyz_file(int nAtoms, atom_struct atomList[], const char* outFileName, const char* comment) {
  const int maxLabelLen = 4;
  char paddingString[888];
  int k, len;
  FILE* f = fopen(outFileName, "wt");
  if(f == NULL) {
    printf("error opening file '%s' for writing\n", outFileName);
    return -1;
  }
  // First line should cantain the number of atoms
  fprintf(f, "%i\n", nAtoms);
  // Second line is a comment line
  fprintf(f, "%s\n", comment);
  for(int i = 0; i < nAtoms; i++) {
    len = strlen(atomList[i].label.c_str());
    if(len > maxLabelLen) {
      printf("error: label too long: '%s'\n", atomList[i].label.c_str());
      return -1;
    }
    for(k = 0; k < maxLabelLen-len; k++)
      paddingString[k] = ' ';
    paddingString[k] = '\0';
    fprintf(f, "%s%s %15.9f %15.9f %15.9f\n",
	    atomList[i].label.c_str(),
	    paddingString,
	    atomList[i].coord[0],
	    atomList[i].coord[1],
	    atomList[i].coord[2]);
  } // END FOR i
  fclose(f);
  printf("File '%s' created OK\n", outFileName);
  return 0;
}

static double get_dist(vec3d_struct const & pos1, vec3d_struct const & pos2) {
  double sum_of_squares = 0;
  for(int i = 0; i < 3; i++) {
    double d = pos1.coord[i] - pos2.coord[i];
    sum_of_squares += d*d;
  }
  return sqrt(sum_of_squares);
}

static double get_vec_length(vec3d_struct const & A) {
  double sum = 0;
  for(int i = 0; i < 3; i++)
    sum += A.coord[i]*A.coord[i];
  return sqrt(sum);
}

static double get_dot_product(vec3d_struct const & A, vec3d_struct const & B) {
  double sum = 0;
  for(int i = 0; i < 3; i++)
    sum += A.coord[i]*B.coord[i];
  return sum;
}

static void normalize_vector(vec3d_struct & A) {
  double d = get_vec_length(A);
  for(int i = 0; i < 3; i++)
    A.coord[i] /= d;
}

static void add_to_vec(vec3d_struct & A, vec3d_struct const & B, double factor) {
  for(int i = 0; i < 3; i++)
    A.coord[i] += factor * B.coord[i];
}

int main(int argc, char* argv[])
{
  printf("generate_alkane version 1.0\n");
  if(argc != 3) {
    printf("Usage:\n");
    printf("generate_alkane N filename\n");
    printf("where\n");
    printf("   N        : number of C atoms.\n");
    printf("   filename : filename of new xyz file to create.\n");
    return -1;
  }
  int N = atoi(argv[1]);
  if(N < 1) {
    printf("Error: (N < 1)\n");
    return -1;
  }
  const char* filename = argv[2];
  printf("N = %d\n", N);
  printf("filename = '%s'\n", filename);

  const double C_C_dist = 1.54;
  const double C_H_dist = 1.09;

  // To get the angles right, start by setting up unit vectors pointing in the directions we need.
  // Create perfect tetrahedron A-B-C-D
  double h = 1 / sqrt(8);
  vec3d_struct A(0, 1, -h);
  vec3d_struct B(sqrt(3)/2, -0.5, -h);
  vec3d_struct C(-sqrt(3)/2, -0.5, -h);
  vec3d_struct D(0, 0, sqrt(1+h*h));
  // Check that tetrahedron ABCD is correct; distance to origin should be same for all of the points, and the distances between them should also be the same.
  double lengthA = get_vec_length(A);
  double lengthB = get_vec_length(B);
  double lengthC = get_vec_length(C);
  double lengthD = get_vec_length(D);
  double tol = 1e-8;
  assert(fabs(lengthA-lengthB) < tol);
  assert(fabs(lengthA-lengthC) < tol);
  assert(fabs(lengthA-lengthD) < tol);
  double distAB = get_dist(A, B);
  double distAC = get_dist(A, C);
  double distAD = get_dist(A, D);
  double distBC = get_dist(B, C);
  double distBD = get_dist(B, D);
  double distCD = get_dist(C, D);
  assert(fabs(distAB-distAC) < tol);
  assert(fabs(distAB-distAD) < tol);
  assert(fabs(distAB-distBC) < tol);
  assert(fabs(distAB-distBD) < tol);
  assert(fabs(distAB-distCD) < tol);
  // Normalize vectors
  normalize_vector(A);
  normalize_vector(B);
  normalize_vector(C);
  normalize_vector(D);

  // Check angle between vectors
  double dot_product_AB = get_dot_product(A, B);
  double dot_product_AC = get_dot_product(A, C);
  double dot_product_AD = get_dot_product(A, D);
  double dot_product_BC = get_dot_product(B, C);
  double dot_product_BD = get_dot_product(B, D);
  double dot_product_CD = get_dot_product(C, D);
  assert(fabs(dot_product_AB-dot_product_AC) < tol);
  assert(fabs(dot_product_AB-dot_product_AD) < tol);
  assert(fabs(dot_product_AB-dot_product_BC) < tol);
  assert(fabs(dot_product_AB-dot_product_BD) < tol);
  assert(fabs(dot_product_AB-dot_product_CD) < tol);
  double angle = acos(dot_product_AB);
  double pi = 2 * asin(1);
  double angle_in_degrees = angle * 180 / pi;
  printf("Angle between vectors: %15.10f <--> %15.10f degrees\n", angle, angle_in_degrees);

  int NC = N;

  // Place C atoms
  std::vector<vec3d_struct> C_pts;
  C_pts.reserve(NC);

  vec3d_struct currPos(0, 0, 0);
  for(int i = 0; i < NC; i++) {
    C_pts.push_back(currPos);
    if(i % 2 == 0) {
      // Move in direction D
      add_to_vec(currPos, D, C_C_dist);
    }
    else {
      // Move in direction -A
      add_to_vec(currPos, A, -C_C_dist);
    }
  }

  // Now add H atoms for each C atom
  std::vector<vec3d_struct> H_pts;
  H_pts.reserve(2*NC+2);
  for(int i = 0; i < NC; i++) {
    vec3d_struct Cpos = C_pts[i];
    if(i == 0) {
      // Special case: first C atom, should have 3 H atoms.
      vec3d_struct Hpos1 = Cpos;
      add_to_vec(Hpos1, A, C_H_dist);
      H_pts.push_back(Hpos1);
      vec3d_struct Hpos2 = Cpos;
      add_to_vec(Hpos2, B, C_H_dist);
      H_pts.push_back(Hpos2);
      vec3d_struct Hpos3 = Cpos;
      add_to_vec(Hpos3, C, C_H_dist);
      H_pts.push_back(Hpos3);
      if(NC == 1) {
	// Add 4th H atom in this case
	vec3d_struct Hpos4 = Cpos;
	add_to_vec(Hpos4, D, C_H_dist);
	H_pts.push_back(Hpos4);
      }
    }
    else if(i == NC-1) {
      // Special case: last C atom, should have 3 H atoms.
      if(i % 2 == 0) {
	vec3d_struct Hpos1 = Cpos;
	add_to_vec(Hpos1, D, C_H_dist);
	H_pts.push_back(Hpos1);
	vec3d_struct Hpos2 = Cpos;
	add_to_vec(Hpos2, B, C_H_dist);
	H_pts.push_back(Hpos2);
	vec3d_struct Hpos3 = Cpos;
	add_to_vec(Hpos3, C, C_H_dist);
	H_pts.push_back(Hpos3);
      }
      else {
	vec3d_struct Hpos1 = Cpos;
	add_to_vec(Hpos1, A, -C_H_dist);
	H_pts.push_back(Hpos1);
	vec3d_struct Hpos2 = Cpos;
	add_to_vec(Hpos2, B, -C_H_dist);
	H_pts.push_back(Hpos2);
	vec3d_struct Hpos3 = Cpos;
	add_to_vec(Hpos3, C, -C_H_dist);
	H_pts.push_back(Hpos3);
      }
    }
    else {
      // Normal case, 2 H atoms.
      if(i % 2 == 0) {
	vec3d_struct Hpos1 = Cpos;
	add_to_vec(Hpos1, B, C_H_dist);
	H_pts.push_back(Hpos1);
	vec3d_struct Hpos2 = Cpos;
	add_to_vec(Hpos2, C, C_H_dist);
	H_pts.push_back(Hpos2);
      }
      else {
	vec3d_struct Hpos1 = Cpos;
	add_to_vec(Hpos1, B, -C_H_dist);
	H_pts.push_back(Hpos1);
	vec3d_struct Hpos2 = Cpos;
	add_to_vec(Hpos2, C, -C_H_dist);
	H_pts.push_back(Hpos2);
      }
    }
  }
  int NH = H_pts.size();
  assert(NH == 2*NC+2);

  // Now create final list of atoms
  int nAtomsTot = NC + NH;
  std::vector<atom_struct> atoms;
  atoms.reserve(nAtomsTot);
  // Add C atoms
  for(int i = 0; i < NC; i++) {
    atom_struct newAtom;
    newAtom.label = "C";
    for(int k = 0; k < 3; k++)
      newAtom.coord[k] = C_pts[i].coord[k];
    atoms.push_back(newAtom);
  }
  // Add H atoms
  for(int i = 0; i < NH; i++) {
    atom_struct newAtom;
    newAtom.label = "H";
    for(int k = 0; k < 3; k++)
      newAtom.coord[k] = H_pts[i].coord[k];
    atoms.push_back(newAtom);
  }

  char comment[888];
  sprintf(comment, "C%dH%d xyz file created by the generate_alkane program, NC = %d, NH = %d, C_C_dist = %f, C_H_dist = %f", NC, NH, NC, NH, C_C_dist, C_H_dist);
  int nAtoms = atoms.size();
  if(write_xyz_file(nAtoms, &atoms[0], filename, comment) != 0) {
    printf("Error in write_xyz_file().\n");
    return -1;
  }

  printf("Done. Alkane xyz file '%s' created, %d C atoms and %d H atoms, %d atoms in total.\n", filename, NC, NH, nAtoms);

  return 0;
}
