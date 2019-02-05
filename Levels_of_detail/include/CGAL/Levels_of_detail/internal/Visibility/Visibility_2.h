#ifndef CGAL_LEVELS_OF_DETAIL_VISIBILITY_2_H
#define CGAL_LEVELS_OF_DETAIL_VISIBILITY_2_H

// STL includes.
#include <vector>
#include <memory>
#include <cstdlib>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enumerations.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename Connectivity,
  typename VisibilityMap,
  typename PolygonFace>
  class Visibility_2 {
			
  public:
    using Traits = GeomTraits;
    using Visibility_map = VisibilityMap;
    using Polygon_face = PolygonFace;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Vector_2 = typename Traits::Vector_2;
    using Triangle_2 = typename Traits::Triangle_2;

    Visibility_2(
      const Connectivity& connectivity,
      Visibility_map visibility_map) :
    m_connectivity(connectivity),
    m_visibility_map(visibility_map),
    m_num_probes(100)
    { }

    void compute(std::vector<Polygon_face>& polygon_faces) const {

      for (std::size_t i = 0; i < polygon_faces.size(); ++i)
        compute_monte_carlo_label(polygon_faces[i]);
    }
    
  private:
    const Connectivity& m_connectivity;
    Visibility_map m_visibility_map;
    const std::size_t m_num_probes;

    void compute_monte_carlo_label(Polygon_face& polygon_face) const {

      FT area = FT(0);
      std::vector< std::pair<Triangle_2, FT> > probe;
      
      const auto& vertices = polygon_face.vertices;
      probe.reserve(vertices.size() - 2);

      for (std::size_t i = 1; i < vertices.size() - 1; ++i) {
        
        Triangle_2 triangle(vertices[0], vertices[i], vertices[i + 1]);
        probe.push_back(std::make_pair(triangle, area));
        area += CGAL::abs(triangle.area());
      }
              
      probe.push_back(std::make_pair(Triangle_2(), area));
              
      FT mean_value = FT(0);
      for (std::size_t i = 0; i < m_num_probes; ++i) {

        std::vector<std::size_t> neighbors;
        m_connectivity.get_neighbors(
          random_point_in_triangles(probe), neighbors);

        CGAL_precondition(neighbors.size() == 1);
        mean_value += get_function_value(neighbors[0]);
      }
      mean_value /= static_cast<FT>(m_num_probes);
              
      if (mean_value >= FT(1) / FT(2))
        polygon_face.visibility = Visibility_label::BUILDING;
      else
        polygon_face.visibility = Visibility_label::OUTSIDE;

      std::cout << int(polygon_face.visibility) << std::endl;
    }

    Point_2 random_point_in_triangles(
      const std::vector< std::pair<Triangle_2, double> >& probe) const {

      const FT key = static_cast<FT>(
      CGAL::to_double(probe.back().second) * 
      (rand() / static_cast<double>(RAND_MAX)));

      for (std::size_t i = 0; i < probe.size() - 1; ++i)
        if (probe[i].second < key && key <= probe[i + 1].second)
          return random_point_in_triangle(probe[i].first);
      
      std::cerr << "Error: probability is out of range" << std::endl;
      return Point_2(FT(0), FT(0));
    }
      
    Point_2 random_point_in_triangle(const Triangle_2& triangle) const {
      
      const Vector_2 v(triangle[0], triangle[1]);
      const Vector_2 w(triangle[0], triangle[2]);

      Point_2 out = triangle[0];

      double rv = rand() / static_cast<double>(RAND_MAX);
      double rw = rand() / static_cast<double>(RAND_MAX);

      if (rv + rw > 1.0) {

        rv = 1.0 - rv;
        rw = 1.0 - rw;
      }

      const FT bv = static_cast<FT>(rv);
      const FT bw = static_cast<FT>(rw);

      out += bv * v;
      out += bw * w;

      return out;
    }
              
    FT get_function_value(const std::size_t point_id) const {
      return get(m_visibility_map, point_id);
    }

  }; // Visibility_2

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_VISIBILITY_2_H
