#ifndef CGAL_LEVELS_OF_DETAIL_GROUND_H
#define CGAL_LEVELS_OF_DETAIL_GROUND_H

// STL includes.
#include <algorithm>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>

// Ground.
#include <CGAL/Levels_of_detail/internal/Ground/Smooth_ground_estimator.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class Ground {

  public:
    using Data_structure = DataStructure;
    
    using Traits = typename DataStructure::Traits;
    using Filtered_range = typename Data_structure::Filtered_range;
    using Point_map = typename Data_structure::Point_map;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;

    using Smooth_ground_estimator = 
    Smooth_ground_estimator<Traits, Filtered_range, Point_map>;

    Ground(Data_structure& data_structure) :
    m_data(data_structure)
    { }

    void make_planar() {

      if (m_data.verbose) 
        std::cout << std::endl << "- Computing planar ground"
        << std::endl;

      if (m_data.verbose) 
        std::cout << "* fitting plane"
      << std::endl;

      const auto& points = m_data.ground_points();
      CGAL_precondition(points.size() >= 3);

      auto& planar_ground = m_data.planar_ground;

      internal::plane_from_points_3(
        points, 
        m_data.point_map, 
        planar_ground.plane);

      std::vector<Point_3> tmp;
      internal::bounding_box_on_plane_3(
        points, 
        m_data.point_map, 
        planar_ground.plane, 
        tmp);

      planar_ground.bounding_box.clear();
      planar_ground.bounding_box.reserve(4);

      const Point_2 p1 = Point_2(tmp[0].x(), tmp[0].y());
      const Point_2 p2 = Point_2(tmp[1].x(), tmp[1].y());
      const Point_2 p3 = Point_2(tmp[2].x(), tmp[2].y());
      const Point_2 p4 = Point_2(tmp[3].x(), tmp[3].y());

      const FT distance = internal::distance_2(p1, p3) / FT(4);

      const Point_2 a = Point_2(p1.x() - distance, p1.y() - distance);
      const Point_2 b = Point_2(p2.x() + distance, p2.y() - distance);
      const Point_2 c = Point_2(p3.x() + distance, p3.y() + distance);
      const Point_2 d = Point_2(p4.x() - distance, p4.y() + distance);

      planar_ground.bounding_box.push_back(a);
      planar_ground.bounding_box.push_back(b);
      planar_ground.bounding_box.push_back(c);
      planar_ground.bounding_box.push_back(d);
    }

    void make_smooth() {

      if (m_data.verbose) 
        std::cout << std::endl << "- Computing smooth ground"
        << std::endl;

      if (m_data.verbose) 
        std::cout << "* creating triangulation"
      << std::endl;
      
      const auto& points = m_data.ground_points();
      CGAL_precondition(points.size() >= 3);

      Smooth_ground_estimator estimator(
        points, 
        m_data.point_map);

      estimator.create_triangles(
        m_data.smooth_ground.triangles);
    }

    template<typename OutputIterator>
    void return_as_polygon(OutputIterator output) const {

      CGAL_precondition(!m_data.planar_ground.bounding_box.empty());
      for (const Point_2& p : m_data.planar_ground.bounding_box)
        *(output++) = Point_3(p.x(), p.y(), 
        internal::position_on_plane_3(p, m_data.planar_ground.plane).z());
    }

    template<
    typename VerticesOutputIterator, 
    typename FacesOutputIterator>
    void return_as_triangle_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const {
      
      const auto& triangles = m_data.smooth_ground.triangles;
      
      internal::Indexer<Point_3> indexer;
      std::size_t num_vertices = 0;

      for (std::size_t i = 0; i < triangles.size(); ++i) {
        cpp11::array<std::size_t, 3> face;
          
        for (std::size_t j = 0; j < 3; ++j) {
          const auto& point = triangles[i][j];

          const std::size_t idx = indexer(point);
          if (idx == num_vertices) {

            *(output_vertices++) = point;
            ++num_vertices;
          }
          face[j] = idx;
        }
        *(output_faces++) = face;
      }
    }

  private:
    Data_structure& m_data;
    
  }; // Ground

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_GROUND_H
