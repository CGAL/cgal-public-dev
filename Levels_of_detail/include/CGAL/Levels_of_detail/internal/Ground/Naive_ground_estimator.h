#ifndef CGAL_LEVELS_OF_DETAIL_NAIVE_GROUND_ESTIMATOR_H
#define CGAL_LEVELS_OF_DETAIL_NAIVE_GROUND_ESTIMATOR_H

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

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Naive_ground_estimator {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Triangle_3 = typename Traits::Triangle_3;

    using Face_info = Face_info<Traits>;
    using Vertex_info = Vertex_info<Traits>;

    using VB = CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, Traits>;
    using FB = CGAL::Triangulation_face_base_with_info_2<Face_info, Traits>;
    
    using CFB = CGAL::Constrained_triangulation_face_base_2<Traits, FB>;
    using TAG = CGAL::Exact_predicates_tag;
    using TDS = CGAL::Triangulation_data_structure_2<VB, CFB>;

    using Triangulation = 
    CGAL::Constrained_Delaunay_triangulation_2<Traits, TDS, TAG>;

    Naive_ground_estimator(
      const Input_range& input_range,
      const Point_map point_map) :
    m_input_range(input_range),
    m_point_map(point_map)
    { }

    void create_triangles(std::vector<Triangle_3>& triangles) const {

      Triangulation triangulation;
      insert_vertices(triangulation);
      
      triangles.clear();
      triangles.reserve(triangulation.number_of_faces());

      for (auto fit = triangulation.finite_faces_begin(); 
      fit != triangulation.finite_faces_end(); ++fit) {

        const Point_2& p1 = fit->vertex(0)->point();
        const Point_2& p2 = fit->vertex(1)->point();
        const Point_2& p3 = fit->vertex(2)->point();

        const FT h1 = fit->vertex(0)->info().height;
        const FT h2 = fit->vertex(1)->info().height;
        const FT h3 = fit->vertex(2)->info().height;

        const Point_3 a = Point_3(p1.x(), p1.y(), h1);
        const Point_3 b = Point_3(p2.x(), p2.y(), h2);
        const Point_3 c = Point_3(p3.x(), p3.y(), h3);

        triangles.push_back(Triangle_3(a, b, c));
      }
    }

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;

    void insert_vertices(Triangulation& triangulation) const {

      std::size_t i = 0;
      triangulation.clear();

      for (auto it = m_input_range.begin(); it != m_input_range.end(); ++it, ++i) {
        if (i % 10 != 0) 
          continue;

        const Point_3& point = get(m_point_map, *it);
        
        const auto vh = triangulation.insert(Point_2(point.x(), point.y()));
        vh->info().height = point.z();
      }
    }

  }; // Naive_ground_estimator

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_NAIVE_GROUND_ESTIMATOR_H
