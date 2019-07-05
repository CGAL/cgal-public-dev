
#ifndef CGAL_SHAPE_REGULARIZATION
#define CGAL_SHAPE_REGULARIZATION

// #include <CGAL/license/Shape_regularization.h>

#include <vector>
#include <utility> // for pairs
#include <set>

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
    using QP_solver = internal::OSQP_solver<Traits, Input_range>;
    // using Point = typename GeomTraits::Point_2;
    using Segment = typename GeomTraits::Segment_2;

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

    void regularize() { 
      // takes instances neighbor_query, RegularizationType and solver.
      //Algorithm implementation:
      //1) Build neighbor graph from input range, use std::set which contains 
      //std::pair e.g (0,2) for segment 1
      //2) build data for QP solver
      //3) call QP solver, send the matrices
      //4) call update() from Rotated_segments_regularization_2 class

      std::vector<std::size_t> neighbors;
      for(std::size_t i = 0; i < m_input_range.size(); ++i) {
        neighbors.clear();
        m_neighbor_query(i, neighbors);
        for (const std::size_t index : neighbors) {
          std::pair<std::size_t, std::size_t> p;
          i < index ? p = std::make_pair(i, index) : p = std::make_pair(index, i);
          m_graph.insert(p);
        }
      }

      //calculate m_t_ijs
      for(auto const &gi : m_graph) {
        m_t_ijs.push_back(m_regularization_type.target_value(gi.first, gi.second));
      }

      std::vector<FT> result_qp;
      m_qp_solver.solve(m_graph, m_regularization_type.get_t_ijs_map(), result_qp);

      m_regularization_type.update(result_qp);

    }
    
  private:
    // Fields.
    Input_range& m_input_range;
    Neighbor_query& m_neighbor_query;
    Regularization_type& m_regularization_type;
    QP_solver m_qp_solver;
    std::set<std::pair<std::size_t, std::size_t>> m_graph;
    std::vector<FT> m_t_ijs;

  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION
