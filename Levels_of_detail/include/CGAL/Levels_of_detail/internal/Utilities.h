#ifndef CGAL_LEVELS_OF_DETAIL_UTILITIES_H
#define CGAL_LEVELS_OF_DETAIL_UTILITIES_H

// STL includes.
#include <map>
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Eigen_diagonalize_traits.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename FT>
  static FT max_value() {
    return FT(1000000000000);
  }

  template<typename Point>
  class Indexer {
  
  public:
    std::size_t operator()(const Point& point) {
      const auto pair = m_indices.insert(
        std::make_pair(
          point, 
          m_indices.size()));
      const auto& item = pair.first;
      const std::size_t idx = item->second;
      return idx;
    }
    void clear() { m_indices.clear(); }

  private:
    std::map<Point, std::size_t> m_indices;
  };

  template<typename Segment_2>
  typename Kernel_traits<Segment_2>::Kernel::FT
  average_spacing_2(
    const std::vector<Segment_2>& segments,
    const std::size_t num_neighbors) {

    using Traits = 
    typename Kernel_traits<Segment_2>::Kernel;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;

    using Local_traits 
    = CGAL::Exact_predicates_inexact_constructions_kernel;
		using Local_FT = typename Local_traits::FT;
    using Local_point_3 = typename Local_traits::Point_3;

    CGAL_precondition(segments.size() > 0);
    std::vector<Local_point_3> points;
    points.reserve(segments.size() * 2);

    for (std::size_t i = 0; i < segments.size(); ++i) {
      const Segment_2& segment = segments[i];

      const Point_2& source = segment.source();
      const Point_2& target = segment.target();

      const Local_FT sx = static_cast<Local_FT>(CGAL::to_double(source.x()));
      const Local_FT sy = static_cast<Local_FT>(CGAL::to_double(source.y()));

      const Local_FT tx = static_cast<Local_FT>(CGAL::to_double(target.x()));
      const Local_FT ty = static_cast<Local_FT>(CGAL::to_double(target.y()));

      points.push_back(Local_point_3(sx, sy, Local_FT(0)));
      points.push_back(Local_point_3(tx, ty, Local_FT(0)));
    }

    const Local_FT average_spacing = 
    CGAL::compute_average_spacing<CGAL::Sequential_tag>(
      points, num_neighbors, CGAL::parameters::point_map(
        CGAL::Identity_property_map<Local_point_3>()).
        geom_traits(Local_traits()));
                
    return static_cast<FT>(average_spacing);
  }

  template<
  typename Point_2,
  typename Line_2>
  typename Kernel_traits<Point_2>::Kernel::FT
  boundary_points_on_line_2(
    const std::vector<Point_2> &points,
    const std::vector<std::size_t> &indices,
    const Line_2 &line,
    Point_2 &p,
    Point_2 &q) {

    using Traits = 
    typename Kernel_traits<Point_2>::Kernel;
    using FT = typename Traits::FT;
    using Vector_2 = typename Traits::Vector_2;

    typename Traits::Compute_scalar_product_2 scalar_product_2;

    FT min_proj_value = max_value<FT>();
    FT max_proj_value = -max_value<FT>();

    const Vector_2 ref_vector = line.to_vector();
    const Point_2& ref_point = points[indices[0]];
    
    for (std::size_t i = 0; i < indices.size(); ++i) {
      const Point_2& point = points[indices[i]];
      
      Vector_2 curr_vector(ref_point, point);
      const FT value = scalar_product_2(curr_vector, ref_vector);
      
      if (value < min_proj_value) {
        min_proj_value = value;
        p = point;
      }
      
      if (value > max_proj_value) {
        max_proj_value = value;
        q = point;
      }
    }
  }

  template<
  typename Point_2,
  typename Line_2>
  typename Kernel_traits<Point_2>::Kernel::FT
  points_squared_length_2(
    const std::vector<Point_2> &points,
    const std::vector<std::size_t> &indices,
    const Line_2 &line) {

    using Traits = 
    typename Kernel_traits<Point_2>::Kernel;
    using FT = typename Traits::FT;
    
    typename Traits::Compute_squared_distance_2 squared_distance_2;

    Point_2 p, q;
    boundary_points_on_line_2(points, indices, line, p, q);

    const FT squared_length = 
    squared_distance_2(line.projection(p), line.projection(q));

    return squared_length;
  }

	template<
  typename Point_2, 
  typename Line_2>
  typename Kernel_traits<Point_2>::Kernel::FT
	line_from_points_2(
    const std::vector<Point_2>& points, 
    const std::vector<std::size_t>& indices,
    Line_2& line) {

    using Traits = 
    typename Kernel_traits<Point_2>::Kernel;    
    using FT = typename Traits::FT;

    using Local_traits 
    = CGAL::Exact_predicates_inexact_constructions_kernel;
		using Local_FT = typename Local_traits::FT;
		using Local_line_2 = typename Local_traits::Line_2;
    using Local_point_2 = typename Local_traits::Point_2;

		using Diagonalize_traits = CGAL::Eigen_diagonalize_traits<Local_FT, 2>;

    CGAL_precondition(indices.size() >= 2);
    
    std::vector<Local_point_2> local_points(indices.size());
    for (std::size_t i = 0; i < indices.size(); ++i) {

      CGAL_precondition(indices[i] >= 0 && indices[i] < points.size());
      const Point_2& point = points[indices[i]];

      const Local_FT x = static_cast<Local_FT>(CGAL::to_double(point.x()));
      const Local_FT y = static_cast<Local_FT>(CGAL::to_double(point.y()));

      local_points[i] = Local_point_2(x, y);
    }

    Local_line_2 fitted_line;
    Local_point_2 fitted_centroid;

    const FT quality = static_cast<FT>(
      CGAL::linear_least_squares_fitting_2(
        local_points.begin(), local_points.end(), 
        fitted_line, fitted_centroid, CGAL::Dimension_tag<0>(),
        Local_traits(), Diagonalize_traits()));

    line = Line_2(
      static_cast<FT>(fitted_line.a()), 
      static_cast<FT>(fitted_line.b()), 
      static_cast<FT>(fitted_line.c()));

    return quality;
  }

	template<
  typename Items, 
  typename Point_map, 
  typename Plane_3>
	typename Kernel_traits<
  typename boost::property_traits<Point_map>::value_type>::Kernel::FT
  plane_from_points_3(
    const Items& items, 
    const Point_map point_map, 
    Plane_3& plane) {

    using Traits = 
    typename Kernel_traits<
    typename boost::property_traits<Point_map>::value_type>::Kernel;
    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;

    using Local_traits 
    = CGAL::Exact_predicates_inexact_constructions_kernel;
		using Local_FT = typename Local_traits::FT;
		using Local_point_3 = typename Local_traits::Point_3;
		using Local_plane_3 = typename Local_traits::Plane_3;

		using Diagonalize_traits = CGAL::Eigen_diagonalize_traits<Local_FT, 3>;

		CGAL_precondition(items.size() >= 3);
		std::vector<Local_point_3> local_points(items.size());
				
    std::size_t i = 0;
		for (auto it = items.begin(); it != items.end(); ++it, ++i) {
			const Point_3& point = get(point_map, *it);

			const Local_FT x = static_cast<Local_FT>(CGAL::to_double(point.x()));
			const Local_FT y = static_cast<Local_FT>(CGAL::to_double(point.y()));
			const Local_FT z = static_cast<Local_FT>(CGAL::to_double(point.z()));

			local_points[i] = Local_point_3(x, y, z);
		}

		Local_plane_3 fitted_plane;
    Local_point_3 fitted_centroid;

		const FT quality = static_cast<FT>(
      CGAL::linear_least_squares_fitting_3(
        local_points.begin(), local_points.end(), 
        fitted_plane, fitted_centroid, CGAL::Dimension_tag<0>(),
        Local_traits(), Diagonalize_traits()));

		plane = Plane_3(
      static_cast<FT>(fitted_plane.a()), 
      static_cast<FT>(fitted_plane.b()), 
      static_cast<FT>(fitted_plane.c()), 
      static_cast<FT>(fitted_plane.d()));

    return quality;
	}

  template<
  typename Segment_2, 
  typename Point_2>
  void bounding_box_2(
    const std::vector<Segment_2>& segments,
    std::vector<Point_2>& bounding_box) {
    
    using Traits = typename Kernel_traits<Segment_2>::Kernel;
    using FT = typename Traits::FT;

    CGAL_precondition(segments.size() > 0);

    FT minx = max_value<FT>();
    FT miny = max_value<FT>();
    FT maxx = -max_value<FT>();
    FT maxy = -max_value<FT>();

    for (std::size_t i = 0; i < segments.size(); ++i) {
      const Segment_2& segment = segments[i];
                      
      const Point_2& source = segment.source();
      const Point_2& target = segment.target();

      minx = CGAL::min(minx, source.x()); minx = CGAL::min(minx, target.x());
      miny = CGAL::min(miny, source.y()); miny = CGAL::min(miny, target.y());

      maxx = CGAL::max(maxx, source.x()); maxx = CGAL::max(maxx, target.x());
      maxy = CGAL::max(maxy, source.y()); maxy = CGAL::max(maxy, target.y());
    }

    bounding_box.clear();
    bounding_box.resize(4);

    bounding_box[0] = Point_2(minx, miny);
    bounding_box[1] = Point_2(maxx, miny);
    bounding_box[2] = Point_2(maxx, maxy);
    bounding_box[3] = Point_2(minx, maxy);
  }

  template<
  typename Items, 
  typename Point_map, 
  typename Plane_3, 
  typename Point_3>
  void bounding_box_on_plane_3(
    const Items& items, 
    const Point_map& point_map,
    const Plane_3& plane, 
    std::vector<Point_3>& bounding_box) {

    using Traits = 
    typename Kernel_traits<
    typename boost::property_traits<Point_map>::value_type>::Kernel;
    using FT = typename Traits::FT;
                
    CGAL_precondition(items.size() > 0);

    FT minx = max_value<FT>();
    FT miny = max_value<FT>();
    FT maxx = -max_value<FT>();
    FT maxy = -max_value<FT>();
                
    FT z = FT(0), size = FT(0);
    for (auto it = items.begin(); it != items.end(); ++it, size += FT(1)) {
					
      const Point_3& point = get(point_map, *it);
      const Point_3 projected = plane.projection(point);

      minx = CGAL::min(minx, projected.x());
      miny = CGAL::min(miny, projected.y());

      maxx = CGAL::max(maxx, projected.x());
      maxy = CGAL::max(maxy, projected.y());

      z += projected.z();
    }
    z /= size;

    bounding_box.clear();
    bounding_box.resize(4);

    bounding_box[0] = Point_3(minx, miny, z);
    bounding_box[1] = Point_3(maxx, miny, z);
    bounding_box[2] = Point_3(maxx, maxy, z);
    bounding_box[3] = Point_3(minx, maxy, z);
  }

  template<typename Point_3>
  typename Kernel_traits<Point_3>::Kernel::Point_2
  point_2_from_point_3(const Point_3& point_3) {

    return typename Kernel_traits<Point_3>::Kernel::Point_2(
      point_3.x(), point_3.y());
  }

  template<
  typename Point_2, 
  typename Plane_3>
  typename Kernel_traits<Point_2>::Kernel::Point_3
  position_on_plane_3(const Point_2& point, const Plane_3& plane) {

    using Traits = typename Kernel_traits<Point_2>::Kernel;
    using FT = typename Traits::FT;
    using Line_3 = typename Traits::Line_3;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    using Intersect_3 = typename Traits::Intersect_3;

    static Vector_3 vertical(FT(0), FT(0), FT(1));
    const Line_3 line(Point_3(point.x(), point.y(), FT(0)), vertical);
  
    typename CGAL::cpp11::result_of<Intersect_3(Line_3, Plane_3)>::type
      inter = CGAL::intersection(line, plane);
      
    if (inter)
      if (const Point_3* p = boost::get<Point_3>(&*inter))
        return *p;

    std::cerr << 
      "Error (position_on_plane): cannot compute the 3D position!" 
    << std::endl;

    return Point_3(FT(0), FT(0), FT(0));
  }

  template<typename GeomTraits>
  struct Point_3_from_point_2_and_plane {

  public:
    using Traits = GeomTraits;

    using argument_type = typename Traits::Point_2;
    using result_type = typename Traits::Point_3;

    using Plane_3 = typename Traits::Plane_3;
    const Plane_3& m_plane;

    Point_3_from_point_2_and_plane(const Plane_3& plane) : 
    m_plane(plane) 
    { }

    result_type operator()(const argument_type& point_2) const {
      return position_on_plane_3(point_2, m_plane);
    }
  };

  template<typename Point_2>
  typename Kernel_traits<Point_2>::Kernel::FT
  distance_2(
    const Point_2& p, 
    const Point_2& q) {
      
    using Traits = typename Kernel_traits<Point_2>::Kernel;
    using FT = typename Traits::FT;
    
    typename Traits::Compute_squared_distance_2 squared_distance_2;

    return static_cast<FT>(
        CGAL::sqrt(
          CGAL::to_double(
            squared_distance_2(p, q))));
  }

  template<typename Point_2, typename Line_2>
  typename Kernel_traits<Point_2>::Kernel::FT
  distance_2(
    const Point_2& p, 
    const Line_2& line) {
      
    using Traits = typename Kernel_traits<Point_2>::Kernel;
    using FT = typename Traits::FT;
    
    typename Traits::Compute_squared_distance_2 squared_distance_2;

    return static_cast<FT>(
        CGAL::sqrt(
          CGAL::to_double(
            squared_distance_2(p, line))));
  }

  template<typename Vector_2>
  typename Kernel_traits<Vector_2>::Kernel::FT
  cos_angle_2(
    const Vector_2& p,
    const Vector_2& q) {
    
    using Traits = typename Kernel_traits<Vector_2>::Kernel;
    typename Traits::Compute_scalar_product_2 scalar_product_2;

    return CGAL::abs(
      scalar_product_2(p, q));
  }

  template<typename GeomTraits>
  struct Segment_3_from_points_and_plane {

  public:
    using Traits = GeomTraits;
    using Point_2 = typename Traits::Point_2;
    using Line_2 = typename Traits::Line_2;
    using Plane_3 = typename Traits::Plane_3;

    using Points_2 = std::vector<Point_2>;

    using argument_type = std::vector<std::size_t>;
    using result_type = typename Traits::Segment_3;

    const Points_2& m_points;
    const Plane_3& m_plane;

    Segment_3_from_points_and_plane(
      const Points_2& points,
      const Plane_3& plane) : 
    m_points(points),
    m_plane(plane) 
    { }

    result_type operator()(const argument_type& indices) const {

      Line_2 line;
      line_from_points_2(m_points, indices, line);

      Point_2 p, q;
      boundary_points_on_line_2(m_points, indices, line, p, q);

      return result_type(
        position_on_plane_3(p, m_plane),
        position_on_plane_3(q, m_plane));
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_UTILITIES_H
