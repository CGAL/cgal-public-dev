#ifndef CGAL_LEVELS_OF_DETAIL_BOUNDARY_EXTRACTOR_2_H
#define CGAL_LEVELS_OF_DETAIL_BOUNDARY_EXTRACTOR_2_H

// STL includes.
#include <vector>

// Internal includes.
#include <CGAL/Levels_of_detail/enumerations.h>
#include <CGAL/Levels_of_detail/internal/structures.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Boundary_extractor_2 {

  public:
    using Traits = GeomTraits;
    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;

    using Polygon_face_2 = Polygon_face_2<Traits>;

    void create_segments(
      const std::vector<Polygon_face_2>& polygon_faces,
      const std::vector<std::size_t>& indices,
      std::vector<Segment_2>& segments) const {

      segments.clear();
      
      for (std::size_t i = 0; i < indices.size(); ++i) {
        const Polygon_face_2& polygon_face = polygon_faces[indices[i]];

        const auto& neighbors = polygon_face.neighbors;
        for (std::size_t j = 0; j < neighbors.size(); ++j) {
          
          const Polygon_face_2& neighbor_face = polygon_faces[neighbors[j]];
          if (neighbor_face.visibility == Visibility_label::OUTSIDE) 
            segments.push_back(polygon_face.edges[j]);
        }
      }
    }
    
  }; // Boundary_extractor_2

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_BOUNDARY_EXTRACTOR_2_H
