#ifndef CGAL_LEVELS_OF_DETAIL_PLANAR_GROUND_ESTIMATOR_H
#define CGAL_LEVELS_OF_DETAIL_PLANAR_GROUND_ESTIMATOR_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/structures.h>
#include <CGAL/Levels_of_detail/internal/utilities.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Planar_ground_estimator {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;
    using Triangle_2 = typename Traits::Triangle_2;
    using Plane_3 = typename Traits::Plane_3;

    using Vertex_info = Vertex_info<Traits>;
    using Face_info = Face_info<Traits>;

    using VB = CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, Traits>;
    using FB = CGAL::Triangulation_face_base_with_info_2<Face_info, Traits>;
    
    using CFB = CGAL::Constrained_triangulation_face_base_2<Traits, FB>;
    using TAG = CGAL::Exact_predicates_tag;
    using TDS = CGAL::Triangulation_data_structure_2<VB, CFB>;

    using Triangulation = 
    CGAL::Constrained_Delaunay_triangulation_2<Traits, TDS, TAG>;

    using Vertex_handle = typename Triangulation::Vertex_handle;
    using Face_handle = typename Triangulation::Finite_faces_iterator;

    typename Traits::Orientation_2 orientation_2;

    Planar_ground_estimator(const Plane_3& ground_plane) :
    m_ground_plane(ground_plane)
    { }

    void initialize(const std::vector<Point_2>& bbox) {

      clear();
      std::vector<Vertex_handle> vhs;
      vhs.reserve(4);

      // Add bounding box vertices.
      for (const Point_2& p : bbox) {
        const Vertex_handle vh = m_triangulation.insert(p);
        vh->info().z = internal::position_on_plane_3(p, m_ground_plane).z();
        vhs.push_back(vh);
      }

      // Add bounding box edges as constraints.
      for (std::size_t i = 0; i < vhs.size(); ++i) {
        const std::size_t ip = (i + 1) % vhs.size();

        if (vhs[i] != vhs[ip])
          m_triangulation.insert_constraint(vhs[i], vhs[ip]);
      }
    }

    template<typename Urban_object>
    void add_urban_object(const Urban_object& object) {

      // Add object boundaries as constraints.
      for (const Segment_2& segment : object.boundaries()) {

        const Point_2& source = segment.source();
        const Point_2& target = segment.target();

        const Vertex_handle svh = m_triangulation.insert(source);
        const Vertex_handle tvh = m_triangulation.insert(target);

        svh->info().z = 
          internal::position_on_plane_3(source, m_ground_plane).z();
        tvh->info().z = 
          internal::position_on_plane_3(target, m_ground_plane).z();

        if (svh != tvh)
          m_triangulation.insert_constraint(svh, tvh);
      }
    }

    template<typename Urban_object>
    void tag_faces(const Urban_object& object) {

      // Tag all faces that belong to this object. 
      for (auto fh = m_triangulation.finite_faces_begin();
      fh != m_triangulation.finite_faces_end(); ++fh) {
        if (is_valid(fh, object.footprint())) {

          fh->info().urban_tag = object.urban_tag();
          fh->info().object_index = object.index();
          fh->info().tagged = true;
        }
      }
    }

    void clear() {
      m_triangulation.clear();
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_as_triangle_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) {
      
      internal::Indexer<Point_2> indexer;
      std::size_t num_vertices = 0;

      for (auto fit = m_triangulation.finite_faces_begin(); 
      fit != m_triangulation.finite_faces_end(); ++fit) {

        cpp11::array<std::size_t, 3> face;
        for (std::size_t k = 0; k < 3; ++k) {
          const Point_2& p = fit->vertex(k)->point();

          const std::size_t idx = indexer(p);
          if (idx == num_vertices) {

            *(output_vertices++) = 
              internal::position_on_plane_3(p, m_ground_plane);
            ++num_vertices;
          }
          face[k] = idx;
        }
        *(output_faces++) = std::make_pair(face, fit->info().urban_tag);
      }
    }

  private:
    const Plane_3& m_ground_plane;
    Triangulation m_triangulation;

    bool is_valid(
      const Face_handle& fh,
      const std::vector<Triangle_2>& footprint) const {
      
      if (fh->info().tagged)
        return false;

      const Point_2& p1 = fh->vertex(0)->point();
      const Point_2& p2 = fh->vertex(1)->point();
      const Point_2& p3 = fh->vertex(2)->point();

      const Triangle_2 triangle = Triangle_2(p1, p2, p3);
      const Point_2 b = internal::triangle_barycenter_2(triangle);

      for (const Triangle_2& tri : footprint) {
        if (tri.has_on_bounded_side(b) || 
            tri.has_on_boundary(b)) return true;
      }
      return false;
    }

  }; // Planar_ground_estimator

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_PLANAR_GROUND_ESTIMATOR_H
