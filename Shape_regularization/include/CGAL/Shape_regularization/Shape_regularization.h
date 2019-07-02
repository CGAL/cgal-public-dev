
#ifndef CGAL_SHAPE_REGULARIZATION
#define CGAL_SHAPE_REGULARIZATION

// #include <CGAL/license/Shape_regularization.h>

#include <vector>
#include <utility> // for pairs
#include <set> // for sets (for graph implementation) 

// CGAL includes.
// #include <CGAL/assertions.h>
// #include <CGAL/property_map.h>

namespace CGAL {
namespace Regularization {

  template<
    typename GeomTraits,
    typename InputRange,
    typename NeighborQuery, 
    typename RegularizationType,
    typename QPSolver >
  class Shape_regularization {

  public:

    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Neighbor_query = NeighborQuery;
    using Regularization_type = RegularizationType;
    using QP_solver = QPSolver;
    using FT = typename GeomTraits::FT;
    // using Point = typename GeomTraits::Point_2;
    using Segment = typename GeomTraits::Segment_2;

    Shape_regularization(
      InputRange& input_range, 
      NeighborQuery& neighbor_query, 
      RegularizationType& regularization_type,
      const QPSolver qp_solver = QPSolver()) :
    m_input_range(input_range),
    m_neighbor_query(neighbor_query),
    m_regularization_type(regularization_type),
    m_qp_solver(qp_solver) {

      CGAL_precondition(input_range.size() > 0);
      // clear();
    }

    void regularize() { 
      // takes instances neighbor_query, RegularizationType and solver.
      //Algorithm implementation:
      //1) Build neighbor graph from input range, use std::set which contains 
      //std::pair e.g (0,2) for segment 1
      //2) build data for QP solver
      //3) call QP solver, send the matrices
      //4) call update() from Rotated_segments_regularization_2 class

      for(int i = 0; i < m_input_range.size(); ++i) {
        std::vector<int> result;
        m_neighbor_query(i, result);
        for(int j = 0; j < result.size(); ++j) {
          std::pair<int, int> p;
          if(i < result[j]) { 
            p = std::make_pair(i, result[j]);
          }
          else {
            p = std::make_pair(result[j], i);
          }
          m_graph.insert(p);
        }
      }

      for(auto const &gi : m_graph) {
        // std::cout << "(" << gi.first << ", " << gi.second << ")" << std::endl;
        //calculate m_t_ijs
        m_t_ijs.push_back(m_regularization_type.target_value(gi.first, gi.second));
      }

      std::cout << std::endl;
      for(int i = 0; i < m_t_ijs.size(); ++i) {
        std::cout << m_t_ijs[i] << std::endl;
      } 

      // m_regularization_type.debug_trmu_ijs();

      // m_qp_solver.test();

      std::vector<FT> result_qp;
      m_qp_solver.solve(m_graph, m_regularization_type.get_t_ijs_map(),  m_regularization_type.get_r_ijs_map(), result_qp);
     /* std::cout << std::endl;
      for (int i = 0; i < result_qp.size(); ++i) {
        std::cout << result_qp[i] << " " << std::endl;
      } */

      m_regularization_type.update(result_qp);

    }
    
  private:
    // Fields.
    Input_range& m_input_range;
    Neighbor_query& m_neighbor_query;
    Regularization_type& m_regularization_type;
    QP_solver m_qp_solver;
    std::set<std::pair<int, int>> m_graph;
    std::vector<FT> m_t_ijs;

  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION