
#ifndef CGAL_SHAPE_REGULARIZATION
#define CGAL_SHAPE_REGULARIZATION

// #include <CGAL/license/Shape_regularization.h>

#include <Eigen/Sparse>
#include <Eigen/Dense> // for vectors of bounds
#include <vector>
#include <utility> // for pairs
#include <set>
#include <map>

#include <CGAL/Shape_regularization/internal/OSQP_solver.h>

// CGAL includes.
// #include <CGAL/assertions.h>
// #include <CGAL/property_map.h>

namespace CGAL {
namespace Regularization {

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
    weight(FT(100000)), 
    lambda(FT(4)/FT(5)), 
    neg_inf(FT(-10000000000)),
    pos_inf(FT(10000000000)),
    val_pos(FT(2) * lambda),
    val_neg(FT(-2) * lambda) {
      
      CGAL_precondition(input_range.size() > 0);

    }

    /* takes instances neighbor_query, RegularizationType and solver.
    Algorithm implementation:
    1) Build neighbor graph from input range, use std::set which contains 
    std::pair e.g (0,2) for segment 1
    2) build data for QP solver
    3) call QP solver, send the matrices
    4) call update() from Rotated_segments_regularization_2 class */
    void regularize() { 

      std::vector<std::size_t> neighbors;
      std::pair<std::size_t, std::size_t> p;
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        neighbors.clear();
        m_neighbor_query(i, neighbors);
        for (const std::size_t index : neighbors) {
          i < index ? p = std::make_pair(i, index) : p = std::make_pair(index, i);
          m_graph.insert(p);
        }
      }

      //get bounds
      m_bounds.reserve(m_input_range.size());
      m_theta_max = FT(0);
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const FT theta = m_regularization_type.bound(i);
        CGAL_postcondition(theta > 0);
        m_bounds.push_back(theta);
        if (m_theta_max < theta) {
          m_theta_max = theta;
        }
      }
      CGAL_postcondition(m_theta_max > 0);
      CGAL_postcondition(m_bounds.size() == m_input_range.size());

      //calculate m_t_ijs
      for(const auto &gi : m_graph) {
        FT t_ij = m_regularization_type.target_value(gi.first, gi.second);
        m_t_ijs[gi] = t_ij;
      }

      build_OSQP_solver_data(); 
      // print_OSQP_solver_data_debug();
      std::vector<FT> result_qp;
      std::size_t n = m_input_range.size() + m_graph.size();
      result_qp.reserve(n);
      m_qp_solver.solve(m_input_range.size(), m_graph.size(), m_P_mat, m_A_mat, m_q, m_l, m_u, result_qp);
      CGAL_postcondition(result_qp.size() == n);

      m_regularization_type.update(result_qp);

    }
    
  private:
    // Fields.
    Input_range& m_input_range;
    Neighbor_query& m_neighbor_query;
    Regularization_type& m_regularization_type;
    QP_solver m_qp_solver;
    std::set <std::pair<std::size_t, std::size_t>> m_graph;
    std::map <std::pair<std::size_t, std::size_t>, FT> m_t_ijs;
    
    // variables for the OSQP solver:
    Sparse_matrix_FT m_P_mat;
    Sparse_matrix_FT m_A_mat;
    Dense_vector_FT m_q;
    Dense_vector_FT m_l;
    Dense_vector_FT m_u;
    const FT weight;
    const FT lambda;
    const FT neg_inf;
    const FT pos_inf;
    const FT val_pos;
    const FT val_neg;
    FT m_theta_max;
    std::vector<FT> m_bounds;


    void print_OSQP_solver_data_debug() {
      std::cout << std::endl << "m_P_mat: " << std::endl << m_P_mat;
      std::cout << std::endl << "m_A_mat: " << std::endl << m_A_mat;
      std::cout << std::endl << "m_q: " << std::endl << m_q;
      std::cout << std::endl << "m_l: " << std::endl << m_l;
      std::cout << std::endl << "m_u: " << std::endl << m_u << std::endl;
    }

    void build_quadratic_matrix(const std::size_t n, const std::size_t k) {

      std::vector<FT_triplet> vec;
      vec.reserve(k);
      FT val;
      for(std::size_t i = 0; i < n; ++i) {
        // i < k ? val = FT(2) * weight * (FT(1) - lambda) / (theta_max * theta_max * k) : val = FT(0);
        i < k ? val = FT(2) * weight * (FT(1) - lambda) / (m_bounds[i] * m_bounds[i] * k) : val = FT(0);
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
        i < k ? m_q[i] = FT(0) : m_q[i] = lambda * weight / (4 * m_theta_max * (n - k));
      }

    }

    void build_linear_constraints_matrix(const std::size_t n, 
                                         const std::size_t m, 
                                         const std::size_t k,
                                         const std::size_t e,
                                         const std::size_t A_nnz) {

      std::vector<FT_triplet> vec;
      vec.reserve(A_nnz);

      std::size_t it = 0, ij = k;
      for (const auto& gi : m_graph) {
        vec.push_back(FT_triplet(it, gi.first, val_neg));
        vec.push_back(FT_triplet(it, gi.second, val_pos));
        vec.push_back(FT_triplet(it, ij, -1));
        ++it;
        vec.push_back(FT_triplet(it, gi.first, val_pos));
        vec.push_back(FT_triplet(it, gi.second, val_neg));
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
      std::set<std::pair<std::size_t, std::size_t>>::iterator gi = m_graph.begin();
      for(std::size_t i = 0; i < m; ++i) {
        if (i < 2 * e) {
          if (i % 2 == 0) 
            m_u[i] = val_neg * m_t_ijs[*gi];
          else {
            m_u[i] = val_pos * m_t_ijs[*gi];
            gi++;
          }
          m_l[i] = neg_inf;
        }
        else if (i < 2 * e + k) {
          m_l[i] = -1 * m_theta_max;
          m_u[i] = m_theta_max;
        }
        else {
          m_l[i] = neg_inf;
          m_u[i] = pos_inf;
        }
      }

    }


    void build_OSQP_solver_data() {

      const std::size_t k = m_input_range.size(); // k segments
      const std::size_t e = m_graph.size(); // e edges
      const std::size_t n = k + e; // number of variables
      const std::size_t m = 2 * e + n; // number of constraints
      const std::size_t A_nnz = 6 * e + n;  // number of entries in the constraint matrix

      build_quadratic_matrix(n, k);
      build_linear_part_vactor(n, k);
      build_linear_constraints_matrix(n, m, k, e, A_nnz);
      build_bounds_vectors(m, k, e);

    }

  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION
