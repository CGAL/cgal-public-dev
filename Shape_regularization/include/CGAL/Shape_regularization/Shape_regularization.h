
#ifndef CGAL_SHAPE_REGULARIZATION
#define CGAL_SHAPE_REGULARIZATION

#include <vector>
#include <utility> // for pairs
#include <set> // for sets (for graph implementation) 

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
    using Neighbor_query = NeighborQuery;
    using Regularization_type = RegularizationType;
    using Input_range = InputRange;
    using QP_solver = QPSolver;

    Shape_regularization(
      const InputRange& input_range, 
      NeighborQuery& neighbor_query, 
      RegularizationType& regularization_type,
      const GeomTraits traits = GeomTraits(),
      const QPSolver = QPSolver()) :
      m_input_range(input_range),
      m_neighbor_query(neighbor_query),
      m_regularization_type(regularization_type) {
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
      for(std::size_t i = 0; i < m_input_range.size(); ++i) {
          
      }
    }
  private:
  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION