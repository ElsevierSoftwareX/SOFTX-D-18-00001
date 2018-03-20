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
/**
 * @file transform.h
 *
 * @author Anastasia Kruchinina
 * \date
 *
 * Contains specializations of template function transform_matrix_from_to
 * for various matrix types. If Chunks and Tasks are not used, than
 * no transformation is needed.
 */

#ifndef HEADER_TRANSFORM
#define HEADER_TRANSFORM

#include <iostream>
#include "matrix_typedefs.h"
#include "matrix_typedefs_chtml.h"

typedef ergo_real real;



template<typename TYPE1, typename TYPE2>
inline void transform_matrix_from_to(const TYPE1& A, TYPE2& B, const ParamsType& P)
{
   throw std::runtime_error("Error in transform_matrix_from_to : it is not implemented for given template parameters.");
}


/*  "FAKE" TRANSFORMATIONS  */

template<>
inline void transform_matrix_from_to<symmMatrix, symmMatrix>
   (const symmMatrix& A, symmMatrix& B, const ParamsType& P)
{
   B = A;
}


template<>
inline void transform_matrix_from_to<normalMatrix, normalMatrix>
   (const normalMatrix& A, normalMatrix& B, const ParamsType& P)
{
   B = A;
}


template<>
inline void transform_matrix_from_to<triangMatrix, triangMatrix>
   (const triangMatrix& A, triangMatrix& B, const ParamsType& P)
{
   B = A;
}


#ifdef USE_CHUNKS_AND_TASKS



/* FROM ERGO TO CHT */


template<typename MatrixType, typename MatrixTypeWrapper>
inline void get_sparse_matrix_data(const MatrixType& X,
                                   vector<int>&      rows,
                                   vector<int>&      cols,
                                   vector<real>&     vals)
{
   throw std::runtime_error("Error in transform.h : get_sparse_matrix_data is not implemented for a given template parameters.");
}


template<>
inline void get_sparse_matrix_data<symmMatrix, chtml::CHTSymmMatrix<real, ParamsType> >(const symmMatrix& X,
                                                                                        vector<int>&      rows,
                                                                                        vector<int>&      cols,
                                                                                        vector<real>&     vals)
{
   rows.clear();
   cols.clear();
   vals.clear();
   X.get_all_values(rows, cols, vals);

   size_t count = 0;
   for (size_t i = 0; i < rows.size(); ++i)
   {
      if (vals[i] == 0)
      {
         continue;
      }
      rows[count] = rows[i];
      cols[count] = cols[i];
      vals[count] = vals[i];
      count++;
   }

   rows.resize(count);
   cols.resize(count);
   vals.resize(count);
}


template<>
inline void get_sparse_matrix_data<symmMatrix, chtml::CHTGeneralMatrix<real, ParamsType> >(const symmMatrix& X,
                                                                                           vector<int>&      rows,
                                                                                           vector<int>&      cols,
                                                                                           vector<real>&     vals)
{
   rows.clear();
   cols.clear();
   vals.clear();
   X.get_all_values(rows, cols, vals);

   size_t count = 0;
   for (size_t i = 0; i < rows.size(); ++i)
   {
      if (vals[i] == 0)
      {
         continue;
      }
      rows[count] = rows[i];
      cols[count] = cols[i];
      vals[count] = vals[i];
      count++;
   }


   // here we have just upper triangle
   // now set the lower triangle

   rows.resize(count);
   cols.resize(count);
   vals.resize(count);
   rows.reserve(count * 2);
   cols.reserve(count * 2);
   vals.reserve(count * 2);

   size_t N = rows.size();
   for (size_t i = 0; i < N; ++i)
   {
      if (rows[i] != cols[i])
      {
         rows.push_back(cols[i]);
         cols.push_back(rows[i]);
         vals.push_back(vals[i]);
         count++;
      }
   }
}


template<>
inline void get_sparse_matrix_data<triangMatrix, chtml::CHTTriangMatrix<real, ParamsType> >(const triangMatrix& X,
                                                                                            vector<int>&        rows,
                                                                                            vector<int>&        cols,
                                                                                            vector<real>&       vals)
{
   rows.clear();
   cols.clear();
   vals.clear();
   X.get_all_values(rows, cols, vals);

   size_t count = 0;
   for (size_t i = 0; i < rows.size(); ++i)
   {
      if (vals[i] == 0)
      {
         continue;
      }
      rows[count] = rows[i];
      cols[count] = cols[i];
      vals[count] = vals[i];
      count++;
   }

   rows.resize(count);
   cols.resize(count);
   vals.resize(count);
}


template<>
inline void get_sparse_matrix_data<normalMatrix, chtml::CHTGeneralMatrix<real, ParamsType> >(const normalMatrix& X,
                                                                                             vector<int>&        rows,
                                                                                             vector<int>&        cols,
                                                                                             vector<real>&       vals)
{
   rows.clear();
   cols.clear();
   vals.clear();
   X.get_all_values(rows, cols, vals);

   size_t count = 0;
   for (size_t i = 0; i < rows.size(); ++i)
   {
      if (vals[i] == 0)
      {
         continue;
      }
      rows[count] = rows[i];
      cols[count] = cols[i];
      vals[count] = vals[i];
      count++;
   }

   rows.resize(count);
   cols.resize(count);
   vals.resize(count);
}


template<>
inline void transform_matrix_from_to<symmMatrix, chtml::CHTSymmMatrix<real, ParamsType> >
   (const symmMatrix& A, chtml::CHTSymmMatrix<real, ParamsType>& B, const ParamsType& P)
{
   B.set_matrix_params(P);

   int n = A.get_nrows();
   int m = A.get_ncols();
   // check dim of B

   vector<int>  rows;
   vector<int>  cols;
   vector<real> vals;
   get_sparse_matrix_data<symmMatrix, chtml::CHTSymmMatrix<real, ParamsType> >(A, rows, cols, vals);
   B.create_CHT_matrix_from_sparse(rows, cols, vals);
}


template<>
inline void transform_matrix_from_to<symmMatrix, chtml::CHTGeneralMatrix<real, ParamsType> >
   (const symmMatrix& A, chtml::CHTGeneralMatrix<real, ParamsType>& B, const ParamsType& P)
{
   B.set_matrix_params(P);

   int n = A.get_nrows();
   int m = A.get_ncols();
   // check dim of B

   vector<int>  rows;
   vector<int>  cols;
   vector<real> vals;
   get_sparse_matrix_data<symmMatrix, chtml::CHTGeneralMatrix<real, ParamsType> >(A, rows, cols, vals);
   B.create_CHT_matrix_from_sparse(rows, cols, vals);
}


template<>
inline void transform_matrix_from_to<normalMatrix, chtml::CHTGeneralMatrix<real, ParamsType> >
   (const normalMatrix& A, chtml::CHTGeneralMatrix<real, ParamsType>& B, const ParamsType& P)
{
   B.set_matrix_params(P);

   int n = A.get_nrows();
   int m = A.get_ncols();
   // check dim of B

   vector<int>  rows;
   vector<int>  cols;
   vector<real> vals;
   get_sparse_matrix_data<normalMatrix, chtml::CHTGeneralMatrix<real, ParamsType> >(A, rows, cols, vals);
   B.create_CHT_matrix_from_sparse(rows, cols, vals);
}


template<>
inline void transform_matrix_from_to<triangMatrix, chtml::CHTTriangMatrix<real, ParamsType> >
   (const triangMatrix& A, chtml::CHTTriangMatrix<real, ParamsType>& B, const ParamsType& P)
{
   B.set_matrix_params(P);

   int n = A.get_nrows();
   int m = A.get_ncols();
   // check dim of B

   vector<int>  rows;
   vector<int>  cols;
   vector<real> vals;
   get_sparse_matrix_data<triangMatrix, chtml::CHTTriangMatrix<real, ParamsType> >(A, rows, cols, vals);
   B.create_CHT_matrix_from_sparse(rows, cols, vals);
}


/* FROM CHT TO ERGO */


template<typename MatrixTypeCHT, typename MatrixType>
inline void set_sparse_matrix_from_data(MatrixType&         X,
                                        const vector<int>&  rows,
                                        const vector<int>&  cols,
                                        const vector<real>& vals)
{
   throw std::runtime_error("Error in transform_matrix_from_to : set_sparse_matrix_from_data is not implemented for a given template parameters.");
}


template<>
inline void set_sparse_matrix_from_data<chtml::CHTGeneralMatrix<real, ParamsType>, symmMatrix>(symmMatrix&         A,
                                                                                               const vector<int>&  rows,
                                                                                               const vector<int>&  cols,
                                                                                               const vector<real>& vals)
{
   // we need just upper triangle
   vector<int> rows1;
   rows1.resize(rows.size());
   vector<int> cols1;
   cols1.resize(cols.size());
   vector<real> vals1;
   vals1.resize(vals.size());
   size_t count = 0;
   for (size_t i = 0; i < rows.size(); ++i)
   {
      if (rows[i] > cols[i])
      {
         continue;
      }
      rows1[count] = rows[i];
      cols1[count] = cols[i];
      vals1[count] = vals[i];
      count++;
   }
   rows1.resize(count);
   cols1.resize(count);
   vals1.resize(count);
   A.assign_from_sparse(rows1, cols1, vals1);
}


template<>
inline void set_sparse_matrix_from_data<chtml::CHTSymmMatrix<real, ParamsType>, symmMatrix>(symmMatrix&         A,
                                                                                            const vector<int>&  rows,
                                                                                            const vector<int>&  cols,
                                                                                            const vector<real>& vals)
{
   A.assign_from_sparse(rows, cols, vals);
}


template<>
inline void transform_matrix_from_to<chtml::CHTSymmMatrix<real, ParamsType>, symmMatrix>
   (const chtml::CHTSymmMatrix<real, ParamsType>& A, symmMatrix& B, const ParamsType& P)
{
   int n, m, nB, mB;

   vector<int>  rows;
   vector<int>  cols;
   vector<real> vals;

   n  = A.get_nrows();
   nB = B.get_nrows();
   assert(nB == n);
   m  = A.get_ncols();
   mB = B.get_ncols();
   assert(mB == m);
   A.get_matrix(rows, cols, vals);
   set_sparse_matrix_from_data<chtml::CHTSymmMatrix<real, ParamsType>, symmMatrix>(B, rows, cols, vals);
}


template<>
inline void transform_matrix_from_to<chtml::CHTGeneralMatrix<real, ParamsType>, symmMatrix>
   (const chtml::CHTGeneralMatrix<real, ParamsType>& A, symmMatrix& B, const ParamsType& P)
{
   int n, m, nB, mB;

   vector<int>  rows;
   vector<int>  cols;
   vector<real> vals;

   n  = A.get_nrows();
   nB = B.get_nrows();
   assert(nB == n);
   m  = A.get_ncols();
   mB = B.get_ncols();
   assert(mB == m);
   A.get_matrix(rows, cols, vals);
   set_sparse_matrix_from_data<chtml::CHTGeneralMatrix<real, ParamsType>, symmMatrix>(B, rows, cols, vals);
}


#endif



#endif // HEADER_TRANSFORM
