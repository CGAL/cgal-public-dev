// Copyright (c) 2019 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Jean-Philippe Bauchet, Florent Lafarge, Gennadii Sytov, Dmitry Anisimov
//

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_OSQP_SOLVER
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_OSQP_SOLVER

#include <CGAL/license/Levels_of_detail.h>

#include <osqp/osqp.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <map>
#include <utility>
#include <set>
#include <vector>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
    typename GeomTraits>
  class OSQP_solver {

  public:
    using Traits = GeomTraits;
    using FT = typename GeomTraits::FT;
    using Sparse_matrix_FT = typename Eigen::SparseMatrix<FT, Eigen::ColMajor>;
    using Sparse_matrix_FT_iterator = typename Sparse_matrix_FT::InnerIterator;
    using Dense_vector_FT = typename Eigen::Matrix<FT, Eigen::Dynamic, 1>;

    void solve(const std::size_t number_of_segments,
               const std::size_t number_of_edges,
               const Sparse_matrix_FT & P_mat,
               const Sparse_matrix_FT & A_mat,
               const Dense_vector_FT & q_v,
               const Dense_vector_FT & l_v,
               const Dense_vector_FT & u_v,
               std::vector<FT> & result) {

      /*
      const c_int n = number_of_segments + number_of_edges; // number of variables
      const c_int m = 2 * number_of_edges + n; // number of constraints
      const c_int P_nnz = n;
      const c_int A_nnz = 6 * number_of_edges + n;

      CGAL_precondition(P_mat.nonZeros() == n);
      CGAL_precondition(A_mat.nonZeros() == A_nnz);
      CGAL_precondition(q_v.nonZeros() == n);
      CGAL_precondition(l_v.nonZeros() == m);
      CGAL_precondition(u_v.nonZeros() == m);

      c_float P_x[n];
      c_int   P_i[n];
      c_int   P_p[n+1];
      build_P_data(n, P_mat, P_x, P_i, P_p);

      c_float A_x[A_nnz];
      c_int   A_i[A_nnz];
      c_int   A_p[n+1];
      build_A_data(A_mat, A_x, A_i, A_p);

      c_float q[n];
      c_float l[m];
      c_float u[m];
      build_vectors(n, m, q_v, l_v, u_v, q, l, u);

      // Problem settings.
      OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

      // Structures
      OSQPWorkspace *work; // Workspace
      OSQPData *data;      // OSQPData

      // Populate data.
      data    = (OSQPData *)c_malloc(sizeof(OSQPData));
      data->n = n;
      data->m = m;
      data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
      data->q = q;
      data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
      data->l = l;
      data->u = u;

      // Define Solver settings as default.
      osqp_set_default_settings(settings);
      settings->eps_rel = 1.0e-15;
      settings->verbose = false;

      // Setup workspace.
      work = osqp_setup(data, settings);

      // Solve Problem.
      osqp_solve(work);

      result.clear();
      c_float *i = work->solution->x;
      for(int j = 0; j < n; j++) {
        result.push_back(FT(i[j]));
      }

      // Clean workspace.
      osqp_cleanup(work);
      c_free(data->A);
      c_free(data->P);
      c_free(data);
      c_free(settings); */
    }

  private:
    void build_P_data(const c_int n, const Sparse_matrix_FT & P_mat,
                      c_float * P_x, c_int * P_i, c_int * P_p) {

      std::size_t it = 0;
      for (int i = 0; i < P_mat.outerSize(); ++i) {
        for (Sparse_matrix_FT_iterator m_i(P_mat, i); m_i; ++m_i) {
          const double val = CGAL::to_double(m_i.value());
          P_x[it] = val;
          ++it;
        }
      }
      P_p[0] = 0;
      for (int i = 0; i < n; ++i) {
        P_i[i] = i;
        P_p[i] = i;
      }
      P_p[n] = n;
    }

    void build_A_data(const Sparse_matrix_FT & A_mat,
                      c_float * A_x, c_int * A_i, c_int * A_p) {

      std::size_t it = 0;
      for (int i = 0; i < A_mat.outerSize(); ++i) {
        for (Sparse_matrix_FT_iterator m_i(A_mat, i); m_i; ++m_i) {
          const double val = CGAL::to_double(m_i.value());
          const std::size_t idx = m_i.row();
          A_x[it] = val;
          A_i[it] = idx;
          ++it;
        }
      }
      A_p[0] = 0;
      for (int i = 1; i <= A_mat.outerSize(); ++i) {
        const std::size_t coln = A_mat.innerVector(i-1).nonZeros();
        A_p[i] = A_p[i-1] + coln;
      }
    }

    void build_vectors(const c_int n, const c_int m,
                       const Dense_vector_FT & q_v,
                       const Dense_vector_FT & l_v,
                       const Dense_vector_FT & u_v,
                       c_float * q, c_float * l, c_float * u) {

      for (int i = 0; i < m; ++i) {
        if (i < n)
          q[i] = CGAL::to_double(q_v[i]);
        l[i] = CGAL::to_double(l_v[i]);
        u[i] = CGAL::to_double(u_v[i]);
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_OSQP_SOLVER
