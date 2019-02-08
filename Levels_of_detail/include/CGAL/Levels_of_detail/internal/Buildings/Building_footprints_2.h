#ifndef CGAL_LEVELS_OF_DETAIL_BUILDING_FOOTPRINTS_2_H
#define CGAL_LEVELS_OF_DETAIL_BUILDING_FOOTPRINTS_2_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/structures.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Building_footprints_2 {

  public:
    using Traits = GeomTraits;
    using Point_2 = typename Traits::Point_2;
    using Triangle_2 = typename Traits::Triangle_2;

    using Polygon_face_2 = Polygon_face_2<Traits>;

    using Face_info = Face_info<Traits>;
    using Vertex_info = Vertex_info<Traits>;

    using VB = CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, Traits>;
    using FB = CGAL::Triangulation_face_base_with_info_2<Face_info, Traits>;
    
    using CFB = CGAL::Constrained_triangulation_face_base_2<Traits, FB>;
    using TAG = CGAL::Exact_predicates_tag;
    using TDS = CGAL::Triangulation_data_structure_2<VB, CFB>;

    using Triangulation = 
    CGAL::Constrained_Delaunay_triangulation_2<Traits, TDS, TAG>;
    using Vertex_handle = typename Triangulation::Vertex_handle;

    void create_footprint_triangles(
      const std::vector<Polygon_face_2>& polygon_faces,
      const std::vector<std::size_t>& indices,
      std::vector<Triangle_2>& triangles) const {

      triangles.clear();
      
      Triangulation triangulation;
      std::vector<Vertex_handle> vhs;
      
      for (std::size_t i = 0; i < indices.size(); ++i) {
        const Polygon_face_2& polygon_face = polygon_faces[indices[i]];

        insert_vertices(polygon_face, vhs, triangulation);
        insert_constraints(vhs, triangulation);

        for (auto fit = triangulation.finite_faces_begin(); 
        fit != triangulation.finite_faces_end(); ++fit) {

          const Point_2& p1 = fit->vertex(0)->point();
          const Point_2& p2 = fit->vertex(1)->point();
          const Point_2& p3 = fit->vertex(2)->point();

          triangles.push_back(Triangle_2(p1, p2, p3));
        }
      }
    }

  private:

    void insert_vertices(
      const Polygon_face_2& polygon_face,
      std::vector<Vertex_handle>& vhs,
      Triangulation& triangulation) const {

      triangulation.clear();
      const auto& vertices = polygon_face.vertices;

      vhs.clear();
      vhs.resize(vertices.size());

      for (std::size_t i = 0; i < vertices.size(); ++i)
        vhs[i] = triangulation.insert(vertices[i]);
    }

    void insert_constraints(
      const std::vector<Vertex_handle>& vhs,
      Triangulation& triangulation) const {

      for (std::size_t i = 0; i < vhs.size(); ++i) {
        const std::size_t ip = (i + 1) % vhs.size();

        if (vhs[i] != vhs[ip])
          triangulation.insert_constraint(vhs[i], vhs[ip]);
      }
    }

  }; // Building_footprints_2

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_BUILDING_FOOTPRINTS_2_H
