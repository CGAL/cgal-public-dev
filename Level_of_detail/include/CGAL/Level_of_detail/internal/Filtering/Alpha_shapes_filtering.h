#ifndef CGAL_LEVEL_OF_DETAIL_ALPHA_SHAPES_FILTERING_H
#define CGAL_LEVEL_OF_DETAIL_ALPHA_SHAPES_FILTERING_H

// CGAL includes.
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Level_of_detail/internal/utils.h>

namespace CGAL {

namespace Level_of_detail {

template<class GeomTraits>
class Alpha_shapes_filtering {

public:
  using Kernel           = GeomTraits;
            
  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;

  using Vb = CGAL::Alpha_shape_vertex_base_2<Kernel>;
  using Fb = CGAL::Alpha_shape_face_base_2<Kernel>;
  using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
  using Triangulation_2 = CGAL::Delaunay_triangulation_2<Kernel, Tds>;
  using Alpha_shape_2   = CGAL::Alpha_shape_2<Triangulation_2>;

  using Alpha_vertex_iterator = typename Alpha_shape_2::Alpha_shape_vertices_iterator;
  using Vertex_handle         = typename Triangulation_2::Vertex_handle;

  Alpha_shapes_filtering(const FT alpha) : 
    m_alpha(alpha) 
  {
    
  }

  template<class InputRange, class PointMap>
  void add_points(const InputRange& range, PointMap point_map) {
    CGAL_precondition(range.size() > 2);

    insert_in_triangulation(range, point_map);
  }

  template <typename Output>
  void get_filtered_points (Output &output, FT sampling)
  {
    CGAL_precondition(m_alpha > FT(0));
    Alpha_shape_2 alpha_shape(m_triangulation, m_alpha, Alpha_shape_2::GENERAL);

#if 0
    for (Alpha_vertex_iterator av_it = alpha_shape.alpha_shape_vertices_begin();
         av_it != alpha_shape.alpha_shape_vertices_end(); ++av_it)
      output.push_back((*av_it)->point());
#else
    for (typename Alpha_shape_2::Alpha_shape_edges_iterator ae_it = alpha_shape.alpha_shape_edges_begin();
         ae_it != alpha_shape.alpha_shape_edges_end(); ++ae_it)
    {
      const Point_2& source = ae_it->first->vertex ((ae_it->second + 1) % 3)->point();
      const Point_2& target = ae_it->first->vertex ((ae_it->second + 2) % 3)->point();

      FT dist = std::sqrt (CGAL::squared_distance (source, target));

      std::size_t nb_pts = std::size_t(dist / sampling) + 1;

      for (std::size_t i = 0; i <= nb_pts; ++ i)
      {
        FT ratio = i / FT(nb_pts);
        output.push_back (Point_2 (source.x() * (1. - ratio) + target.x() * ratio,
                                   source.y() * (1. - ratio) + target.y() * ratio));
      }
    }
#endif
  }

private:

  Triangulation_2 m_triangulation;
  const FT m_alpha;

  template<class InputRange, class PointMap>
  void insert_in_triangulation(const InputRange& range, PointMap point_map) {
                
    for (typename InputRange::const_iterator ce_it = range.begin(); ce_it != range.end(); ++ce_it)
      m_triangulation.insert (internal::point_2_from_point_3 (get(point_map, *ce_it)));
  }
};

} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_ALPHA_SHAPES_FILTERING_H
