#ifndef CGAL_LEVELS_OF_DETAIL_POLYGON_FACES_2_STORED_CONNECTIVITY_H
#define CGAL_LEVELS_OF_DETAIL_POLYGON_FACES_2_STORED_CONNECTIVITY_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/structures.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Polygon_faces_2_stored_connectivity {

  public:
    using Traits = GeomTraits;
    using Polygon_face_2 = Polygon_face_2<Traits>;

    Polygon_faces_2_stored_connectivity(
      const std::vector<Polygon_face_2>& polygon_faces) :
    m_polygon_faces(polygon_faces)
    { }

    void get_neighbors(
      const std::size_t query_index,
      std::vector<std::size_t>& neighbors) const {

      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_polygon_faces.size());

      neighbors.clear();
      neighbors = m_polygon_faces[query_index].neighbors;
    }

  private:
    const std::vector<Polygon_face_2>& m_polygon_faces;

  }; // Polygon_faces_2_stored_connectivity

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_POLYGON_FACES_2_STORED_CONNECTIVITY_H
