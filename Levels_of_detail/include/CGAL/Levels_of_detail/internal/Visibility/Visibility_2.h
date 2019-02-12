#ifndef CGAL_LEVELS_OF_DETAIL_VISIBILITY_2_H
#define CGAL_LEVELS_OF_DETAIL_VISIBILITY_2_H

// STL includes.
#include <vector>
#include <memory>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enumerations.h>
#include <CGAL/Levels_of_detail/internal/structures.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename InputConnectivity,
  typename VisibilityMap>
  class Visibility_2 {
			
  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Connectivity = InputConnectivity;
    using Visibility_map = VisibilityMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Vector_2 = typename Traits::Vector_2;
    using Triangle_2 = typename Traits::Triangle_2;

    using Polygon_face_2 = Polygon_face_2<Traits>;

    Visibility_2(
      const Input_range& input_range,
      const Connectivity& connectivity,
      Visibility_map visibility_map) :
    m_input_range(input_range),
    m_connectivity(connectivity),
    m_visibility_map(visibility_map),
    m_num_probes(100)
    { }

    void compute(std::vector<Polygon_face_2>& polygon_faces) const {

      for (std::size_t i = 0; i < polygon_faces.size(); ++i)
        compute_monte_carlo_label(polygon_faces[i]);
    }
    
  private:
    const Input_range& m_input_range;
    const Connectivity& m_connectivity;
    Visibility_map m_visibility_map;
    const std::size_t m_num_probes;

    void compute_monte_carlo_label(Polygon_face_2& polygon_face) const {

      std::vector< std::pair<Triangle_2, FT> > probability;
      create_probability(polygon_face, probability);
              
      const FT mean_value = 
      compute_mean_value(probability);
              
      if (mean_value >= FT(1) / FT(2))
        polygon_face.visibility = Visibility_label::INSIDE;
      else
        polygon_face.visibility = Visibility_label::OUTSIDE;
    }

    void create_probability(
      const Polygon_face_2& polygon_face, 
      std::vector< std::pair<Triangle_2, FT> >& probability) const {

      const auto& vertices = polygon_face.vertices;
      
      probability.clear();
      probability.reserve(vertices.size() - 2);

      FT area = FT(0);
      for (std::size_t i = 1; i < vertices.size() - 1; ++i) {
        
        Triangle_2 triangle(vertices[0], vertices[i], vertices[i + 1]);
        probability.push_back(std::make_pair(triangle, area));

        area += CGAL::abs(triangle.area());
      }
      probability.push_back(std::make_pair(Triangle_2(), area));
    }

    FT compute_mean_value(
      const std::vector< std::pair<Triangle_2, FT> >& probability) const {

      Point_2 point;
      FT mean_value = FT(0);

      std::vector<std::size_t> neighbors;
      for (std::size_t i = 0; i < m_num_probes; ++i) {

        compute_random_point_in_triangles(probability, point);
        m_connectivity.get_neighbors(point, neighbors);

        CGAL_precondition(neighbors.size() == 1);
        mean_value += get_function_value(neighbors[0]);
      }
      
      mean_value /= static_cast<FT>(m_num_probes);
      return mean_value;
    }

    void compute_random_point_in_triangles(
      const std::vector< std::pair<Triangle_2, FT> >& probability,
      Point_2& point) const {

      const FT key = static_cast<FT>(
      CGAL::to_double(probability.back().second) * 
      (rand() / static_cast<double>(RAND_MAX)));

      for (std::size_t i = 0; i < probability.size() - 1; ++i) {
        if (probability[i].second < key && key <= probability[i + 1].second) {
         
          compute_random_point_in_triangle(probability[i].first, point);
          return;
        }
      }
      
      std::cerr << 
        "Error (compute_random_point_in_triangles()): probability is out of range!" 
      << std::endl;
      point = Point_2(FT(0), FT(0));
    }
    
    void compute_random_point_in_triangle(
      const Triangle_2& triangle,
      Point_2& point) const {
      
      const Vector_2 v(triangle[0], triangle[1]);
      const Vector_2 w(triangle[0], triangle[2]);

      point = triangle[0];

      double rv = rand() / static_cast<double>(RAND_MAX);
      double rw = rand() / static_cast<double>(RAND_MAX);

      if (rv + rw > 1.0) {
        
        rv = 1.0 - rv;
        rw = 1.0 - rw;
      }

      const FT bv = static_cast<FT>(rv);
      const FT bw = static_cast<FT>(rw);

      point += bv * v;
      point += bw * w;
    }
              
    FT get_function_value(const std::size_t index) const {
      return get(m_visibility_map, *(m_input_range.begin() + index));
    }

  }; // Visibility_2

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_VISIBILITY_2_H
