#ifndef CGAL_LEVELS_OF_DETAIL_POLYGON_FACES_2_VISIBILITY_CONDITIONS_H
#define CGAL_LEVELS_OF_DETAIL_POLYGON_FACES_2_VISIBILITY_CONDITIONS_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enumerations.h>
#include <CGAL/Levels_of_detail/internal/structures.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Polygon_faces_2_visibility_conditions {

  public:
    using Traits = GeomTraits;
    using Polygon_face_2 = Polygon_face_2<Traits>;

    Polygon_faces_2_visibility_conditions(
      const std::vector<Polygon_face_2>& polygon_faces) :
    m_polygon_faces(polygon_faces)
    { }

    bool belongs_to_region(
      const std::size_t query_index, 
      const std::vector<std::size_t>& region) const {

      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_polygon_faces.size());
      
      if (region.size() == 1)
        if (m_polygon_faces[region[0]].visibility == Visibility_label::OUTSIDE)
          return false;
      return m_polygon_faces[query_index].visibility == Visibility_label::INSIDE;
    }

    bool is_valid_region(const std::vector<std::size_t>& region) const {
      
      if (region.size() == 1)
        if (m_polygon_faces[region[0]].visibility == Visibility_label::OUTSIDE)
          return false;
      return region.size() >= 1;
    }

    void update(const std::vector<std::size_t>&) {
      // skipped!
    }

  private:
    const std::vector<Polygon_face_2>& m_polygon_faces;

  }; // Polygon_faces_2_visibility_conditions

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_POLYGON_FACES_2_VISIBILITY_CONDITIONS_H
