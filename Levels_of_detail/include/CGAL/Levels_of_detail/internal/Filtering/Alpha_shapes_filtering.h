#ifndef CGAL_LEVELS_OF_DETAIL_ALPHA_SHAPES_FILTERING_H
#define CGAL_LEVELS_OF_DETAIL_ALPHA_SHAPES_FILTERING_H

// CGAL includes.
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

// LOD includes.
#include <CGAL/Levels_of_detail/internal/Utilities.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<class GeometricTraits>
class Alpha_shapes_filtering {

public:
  using Traits = GeometricTraits;
            
  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;

  using Vb = CGAL::Alpha_shape_vertex_base_2<Traits>;
  using Fb = CGAL::Alpha_shape_face_base_2<Traits>;
  using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
  using Triangulation_2 = CGAL::Delaunay_triangulation_2<Traits, Tds>;
  using Alpha_shape_2 = CGAL::Alpha_shape_2<Triangulation_2>;

  Alpha_shapes_filtering(const FT alpha) : 
  m_alpha(alpha) 
  { }

  template<class Input_range, class Point_map>
  void add_points(const Input_range &input_range, Point_map point_map) {
    insert_in_triangulation(input_range, point_map);
  }

  template<typename Output>
  void get_filtered_points(Output &output, FT sampling) {
    
    CGAL_precondition(m_alpha > FT(0));
    CGAL_precondition(m_triangulation.number_of_vertices() >= 3);

    Alpha_shape_2 alpha_shape(m_triangulation, m_alpha, Alpha_shape_2::GENERAL);
    for (
      auto it = alpha_shape.alpha_shape_edges_begin(); 
      it != alpha_shape.alpha_shape_edges_end(); 
      ++it) {

      const Point_2 &source = it->first->vertex((it->second + 1) % 3)->point();
      const Point_2 &target = it->first->vertex((it->second + 2) % 3)->point();

      const FT dist = static_cast<FT>(
        CGAL::sqrt(
          CGAL::to_double(
            CGAL::squared_distance(source, target))));

      CGAL_precondition(sampling > FT(0));
      const std::size_t nb_pts = static_cast<std::size_t>(dist / sampling) + 1;
      
      CGAL_precondition(nb_pts > 0);
      for (std::size_t i = 0; i <= nb_pts; ++i) {

        const FT ratio = static_cast<FT>(i) / static_cast<FT>(nb_pts);
        output.push_back(
          Point_2(
            source.x() * (FT(1) - ratio) + target.x() * ratio,
            source.y() * (FT(1) - ratio) + target.y() * ratio));
      }
    }
  }

private:
  const FT m_alpha;
  Triangulation_2 m_triangulation;

  template<class Input_range, class Point_map>
  void insert_in_triangulation(
    const Input_range &input_range, 
    Point_map point_map) {
                
    for (auto it = input_range.begin(); it != input_range.end(); ++it)
      m_triangulation.insert(
        point_2_from_point_3(get(point_map, *it)));
  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_ALPHA_SHAPES_FILTERING_H
