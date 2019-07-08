
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
    using QP_solver = internal::OSQP_solver<Traits, Input_range>;

    Shape_regularization(
      InputRange& input_range, 
      NeighborQuery& neighbor_query, 
      RegularizationType& regularization_type) :
    m_input_range(input_range),
    m_neighbor_query(neighbor_query),
    m_regularization_type(regularization_type),
    m_qp_solver(QP_solver(input_range)) {
      
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

      //calculate m_t_ijs
      for(auto const &gi : m_graph) {
        FT t_ij = m_regularization_type.target_value(gi.first, gi.second);
        m_t_ijs[gi] = t_ij;
      }

      build_OSQP_solver_data(); 
      // print_OSQP_solver_data_debug();
      std::vector<FT> result_qp;
      m_qp_solver.solve(m_graph, m_t_ijs, result_qp);

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
    Eigen::SparseMatrix<FT, Eigen::ColMajor> m_P_mat;
    Eigen::SparseMatrix<FT, Eigen::ColMajor> m_A_mat;
    Eigen::SparseVector<FT> m_q;
    std::vector <FT> m_l;
    std::vector <FT> m_u;

    void print_OSQP_solver_data_debug() {
      std::cout << std::endl << "m_P_mat: " << std::endl << m_P_mat;
      std::cout << std::endl << "m_A_mat: " << std::endl << m_A_mat;
      std::cout << std::endl << "m_q: " << std::endl << m_q;
      std::cout << std::endl << "m_l: " << std::endl;
      for (std::size_t i = 0; i < m_l.size(); ++i) {
        std::cout << m_l[i] << " ";
      }
      std::cout << std::endl << "m_u: " << std::endl;
      for (std::size_t i = 0; i < m_u.size(); ++i) {
        std::cout << m_u[i] << " ";
      }
      std::cout << std::endl;
    }

    void build_OSQP_solver_data() {
      const std::size_t k = m_input_range.size(); // k segments
      const std::size_t e = m_graph.size();
      const std::size_t n = k + e; // number of variables
      const std::size_t m = 2 * e + n; // number of constraints
      //csc *P - quadratic part of the cost P in csc format (size n x n)
      // std::size_t P_nnz = n;
      // helper variables for calculation
      const std::size_t weight = 100000;
      const FT lambda = FT(4)/FT(5);
      const FT theta_max = FT(10.0);
      const FT neg_inf = -10000000.0;
      const FT pos_inf =  10000000.0;

      std::vector<FT_triplet> vec;
      vec.reserve(k);
      for(std::size_t i = 0; i < k; ++i) {
        FT val = 2 * weight * (1 - lambda) / (theta_max * theta_max * k);
        vec.push_back(FT_triplet(i, i, val));
      }

      m_P_mat.resize(n, n);
      m_P_mat.setFromTriplets(vec.begin(), vec.end());
      m_P_mat.makeCompressed();

      
      m_q.resize(n);
      for (std::size_t i = k; i < n; ++i) {
        m_q.coeffRef(i) = lambda * weight / (4 * theta_max * (n - k));
      }

      const FT val_pos = 2 * lambda;
      const FT val_neg = -2 * lambda;
      //csc *A - linear constraints matrix A in csc format (size m x n) 

      std::size_t A_nnz = 6 * e + n;  
      
      vec.clear();
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

      m_A_mat.resize(m, n);
      m_A_mat.setFromTriplets(vec.begin(), vec.end());
      m_A_mat.makeCompressed();
      
      m_u.clear();
      m_l.clear();
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
          m_l[i] = -1 * theta_max;
          m_u[i] = theta_max;
        }
        else {
          m_l[i] = neg_inf;
          m_u[i] = pos_inf;
        }
      }

    }

  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION
