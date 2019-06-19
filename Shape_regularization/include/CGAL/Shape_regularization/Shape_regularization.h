
#ifndef CGAL_SHAPE_REGULARIZATION
#define CGAL_SHAPE_REGULARIZATION

#include <vector>
#include <utility> // for pairs
#include <set> // for sets (for graph implementation)

namespace CGAL {
namespace Shape_detection {

template<
  typename GeomTraits,
  typename NeighborQuery, 
  typename RegularizationType,
  typename InputRange,
  typename QPSolver 
>
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
            const QPSolver = QPSolver(),
            FT max_angle_in_degrees = FT(10.0),
            FT max_difference_in_meters = FT(0.00001)) :
          m_input_range(input_range),
          m_neighbor_query(neighbor_query),
          m_regularization_type(regularization_type) {
                CGAL_precondition(input_range.size() > 0);
        }

        void regularize() { // takes instances neighbor_query, RegularizationType and solver.
            //Algorithm implementation:
            //1) Build neighbor graph from input range, use std::set which contains std::pair e.g (0,2) for segment 1
            //2) build data for QP solver
            //3) call QP solver, send the matrices
            //4) call update() from Rotated_segments_regularization_2 class
            for(std::size_t i = 0; i < m_input_range.size(); ++i) {
                
            }
        }
    private:
};

} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION