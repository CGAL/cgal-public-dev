#ifndef CGAL_LEVELS_OF_DETAIL_STRUCTURES_H
#define CGAL_LEVELS_OF_DETAIL_STRUCTURES_H

// STL includes.
#include <vector>

// Internal includes.
#include <CGAL/Levels_of_detail/enumerations.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  struct Polygon_face_2 {

  public:
    using Traits = GeomTraits;
    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;

    std::vector<Point_2> vertices;
    Visibility_label visibility;
    std::vector<std::size_t> neighbors;
    std::vector<Segment_2> edges;

    Polygon_face_2() :
    visibility(Visibility_label::OUTSIDE) 
    { }
  };

  template<typename GeomTraits>
  struct Face_info {

  public:
    using Traits = GeomTraits;
  };

  template<typename GeomTraits>
  struct Vertex_info {

  public:
    using Traits = GeomTraits;
  };

  template<typename GeomTraits>
  struct Building {

  public:
    using Traits = GeomTraits;
    using Segment_2 = typename Traits::Segment_2;
    using Triangle_2 = typename Traits::Triangle_2;

    std::vector<Triangle_2> footprint;
    std::vector<Segment_2> boundaries;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_STRUCTURES_H
