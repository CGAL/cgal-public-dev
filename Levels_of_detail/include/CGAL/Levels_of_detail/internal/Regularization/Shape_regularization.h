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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_REGULARIZATION
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_REGULARIZATION

#include <CGAL/license/Levels_of_detail.h>

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <utility>
#include <set>
#include <map>

#include <CGAL/Levels_of_detail/internal/Regularization/OSQP_solver.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
    typename GeomTraits,
    typename InputRange,
    typename NeighborQuery,
    typename RegularizationType>
  class Shape_regularization {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Neighbor_query = NeighborQuery;
    using Regularization_type = RegularizationType;
    using FT = typename GeomTraits::FT;
    using FT_triplet = Eigen::Triplet<FT>;
    using QP_solver = internal::OSQP_solver<Traits>;
    using Sparse_matrix_FT = typename Eigen::SparseMatrix<FT, Eigen::ColMajor>;
    using Dense_vector_FT = typename Eigen::Matrix<FT, Eigen::Dynamic, 1>;

    Shape_regularization(
      InputRange& input_range,
      NeighborQuery& neighbor_query,
      RegularizationType& regularization_type) :
    m_input_range(input_range),
    m_neighbor_query(neighbor_query),
    m_regularization_type(regularization_type),
    m_qp_solver(QP_solver()),
    m_parameters(Parameters()) { }

    void regularize() {
      if(m_input_range.size() < 2) return;

      m_graph.clear();
      build_graph_of_neighbors();
      if(m_graph.size() == 0)
        return;
      CGAL_postcondition(m_graph.size() > 0);

      m_bounds.clear();
      m_bounds.reserve(m_input_range.size());
      m_max_bound = FT(0);

      obtain_bounds();

      CGAL_postcondition(m_max_bound > 0);
      CGAL_postcondition(m_bounds.size() == m_input_range.size());

      m_targets.clear();

      obtain_targets();
      if(m_targets.size() == 0) return;
      CGAL_postcondition(m_targets.size() > 0);

      build_OSQP_solver_data();

      std::vector<FT> result_qp;
      std::size_t n = m_input_range.size() + m_targets.size();
      result_qp.reserve(n);

      m_qp_solver.solve(m_input_range.size(), m_targets.size(), m_P_mat, m_A_mat, m_q, m_l, m_u, result_qp);
      CGAL_postcondition(result_qp.size() == n);

      m_regularization_type.update(result_qp);
    }

  private:
    class Parameters {
      public:
        const FT m_weight;
        const FT m_lambda;
        const FT m_neg_inf;
        const FT m_pos_inf;
        const FT m_val_pos;
        const FT m_val_neg;

        Parameters():
        m_weight(FT(100000)),
        m_lambda(FT(4)/FT(5)),
        m_neg_inf(FT(-10000000000)),
        m_pos_inf(FT(10000000000)),
        m_val_pos(FT(2) * m_lambda),
        m_val_neg(FT(-2) * m_lambda) { }
    };

  private:
    Input_range& m_input_range;
    Neighbor_query& m_neighbor_query;
    Regularization_type& m_regularization_type;
    QP_solver m_qp_solver;
    std::set <std::pair<std::size_t, std::size_t>> m_graph;
    std::map <std::pair<std::size_t, std::size_t>, FT> m_targets;
    const Parameters m_parameters;

    // Variables for the OSQP solver:
    Sparse_matrix_FT m_P_mat;
    Sparse_matrix_FT m_A_mat;
    Dense_vector_FT m_q;
    Dense_vector_FT m_l;
    Dense_vector_FT m_u;
    FT m_max_bound;
    std::vector<FT> m_bounds;

    void build_graph_of_neighbors() {
      std::vector <std::size_t> neighbors;
      std::pair<std::size_t, std::size_t> p;

      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        neighbors.clear();
        m_neighbor_query(i, neighbors);
        for (const std::size_t index : neighbors) {
          i < index ? p = std::make_pair(i, index) : p = std::make_pair(index, i);
          m_graph.insert(p);
        }
      }
    }

    void obtain_bounds() {
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const FT bound = m_regularization_type.bound(i);
        CGAL_postcondition(bound > 0);

        m_bounds.push_back(bound);
        if (m_max_bound < bound)
          m_max_bound = bound;
      }
    }

    void obtain_targets() {
      for(const auto &gi : m_graph) {
        const std::size_t i = gi.first;
        const std::size_t j = gi.second;

        FT tar_val = m_regularization_type.target_value(i, j);
        if (CGAL::abs(tar_val) < m_regularization_type.bound(i) + m_regularization_type.bound(j))
          m_targets[gi] = tar_val;
      }
    }

    void build_quadratic_matrix(const std::size_t n, const std::size_t k) {
      std::vector<FT_triplet> vec;
      vec.reserve(k);

      for(std::size_t i = 0; i < n; ++i) {
        FT val = FT(0);
        if (i < k)
          val = FT(2) * m_parameters.m_weight * (FT(1) - m_parameters.m_lambda) / (m_bounds[i] * m_bounds[i] * FT(k));

        vec.push_back(FT_triplet(i, i, val));
      }
      CGAL_postcondition(vec.size() == n);

      m_P_mat.resize(n, n);
      m_P_mat.setFromTriplets(vec.begin(), vec.end());
      m_P_mat.makeCompressed();
    }

    void build_linear_part_vactor(const std::size_t n, const std::size_t k) {
      m_q.resize(n);
      for (std::size_t i = 0; i < n; ++i) {
        i < k ? m_q[i] = FT(0) : m_q[i] = m_parameters.m_lambda * m_parameters.m_weight / (FT(4) * m_max_bound * FT(n - k));
      }
    }

    void build_linear_constraints_matrix(const std::size_t n,
                                         const std::size_t m,
                                         const std::size_t k,
                                         const std::size_t e,
                                         const std::size_t A_nnz) {
      std::vector<FT_triplet> vec;
      vec.reserve(A_nnz);

      std::size_t it = 0;
      std::size_t ij = k;

      for (const auto& ti : m_targets) {
        const std::size_t i = ti.first.first;
        const std::size_t j = ti.first.second;

        vec.push_back(FT_triplet(it, i, m_parameters.m_val_neg));
        vec.push_back(FT_triplet(it, j, m_parameters.m_val_pos));
        vec.push_back(FT_triplet(it, ij, -1));
        ++it;

        vec.push_back(FT_triplet(it, i, m_parameters.m_val_pos));
        vec.push_back(FT_triplet(it, j, m_parameters.m_val_neg));
        vec.push_back(FT_triplet(it, ij, -1));
        ++it;
        ++ij;
      }

      for (std::size_t i = e * 2; i < m; i++) {
        for (std::size_t j = 0; j < n; j++) {
          if(j == i - e * 2)
            vec.push_back(FT_triplet(i, j, 1));
        }
      }
      CGAL_postcondition(vec.size() == A_nnz);

      m_A_mat.resize(m, n);
      m_A_mat.setFromTriplets(vec.begin(), vec.end());
      m_A_mat.makeCompressed();
    }

    void build_bounds_vectors(const std::size_t m, const std::size_t k, const std::size_t e) {
      m_u.resize(m);
      m_l.resize(m);
      auto ti = m_targets.begin();

      for(std::size_t i = 0; i < m; ++i) {

        if (i < 2 * e) {
          const FT val = ti->second;
          if (i % 2 == 0)
            m_u[i] =  m_parameters.m_val_neg * val;
          else {
            m_u[i] =  m_parameters.m_val_pos * val;
            ++ti;
          }
          m_l[i] = m_parameters.m_neg_inf;
        }
        else if (i < 2 * e + k) {
          m_l[i] = -1 * m_max_bound;
          m_u[i] = m_max_bound;
        }
        else {
          m_l[i] =  m_parameters.m_neg_inf;
          m_u[i] =  m_parameters.m_pos_inf;
        }
      }
    }

    void build_OSQP_solver_data() {
      const std::size_t k = m_input_range.size(); // k segments
      const std::size_t e = m_targets.size(); // e edges
      const std::size_t n = k + e; // number of variables
      const std::size_t m = 2 * e + n; // number of constraints
      const std::size_t A_nnz = 6 * e + n;  // number of entries in the constraint matrix

      build_quadratic_matrix(n, k);
      build_linear_part_vactor(n, k);
      build_linear_constraints_matrix(n, m, k, e, A_nnz);
      build_bounds_vectors(m, k, e);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_SHAPE_REGULARIZATION
