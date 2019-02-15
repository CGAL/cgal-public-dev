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
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename FT>
  static FT max_value() {
    return FT(1000000000000);
  }

  template<
  typename GeomTraits,
  typename InputRange,  
  typename PointMap>
  class Point_2_from_index_map {

  public:    
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using Point_2 = typename Traits::Point_2;

    using key_type = std::size_t;
    using value_type = Point_2;
    using reference = const value_type&;
    using category = boost::lvalue_property_map_tag;

    const Input_range& m_input_range;
    const Point_map m_point_map;

    Point_2_from_index_map(
      const Input_range& input_range,
      const Point_map point_map) : 
    m_input_range(input_range),
    m_point_map(point_map)
    { }

    friend reference get(
      const Point_2_from_index_map& point_2_from_index_map, 
      const key_type index) {

      const auto& point_3 = get(
        point_2_from_index_map.m_point_map, 
        *(point_2_from_index_map.m_input_range.begin() + index));

      // Hack to satisfy Spatial_searching classes, which require lvalue
      // property map. Here, as a CGAL::Point_2 is basically a
      // CGAL::Point_3 if we forget the third coordinate (both are based
      // on a cpp11::array<FT, 3/2>), we can just reinterpret the
      // reference and it works well.
      return reinterpret_cast<const value_type&>(point_3);
    }
  };

  template<
  typename GeomTraits,
  typename InputRange,  
  typename PointMap>
  class Point_3_from_index_map {

  public:    
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using Point_3 = typename Traits::Point_3;

    using key_type = std::size_t;
    using value_type = Point_3;
    using reference = const value_type&;
    using category = boost::lvalue_property_map_tag;

    const Input_range& m_input_range;
    const Point_map m_point_map;

    Point_3_from_index_map(
      const Input_range& input_range,
      const Point_map point_map) : 
    m_input_range(input_range),
    m_point_map(point_map)
    { }

    friend reference get(
      const Point_3_from_index_map& point_3_from_index_map, 
      const key_type index) {

      return get(point_3_from_index_map.m_point_map, 
        *(point_3_from_index_map.m_input_range.begin() + index));
    }
  };

  template<typename Point>
  class Index_to_point_map {
                        
  public:
    using value_type = Point;
    using reference = const value_type&;
    using key_type = std::size_t;
    using category = boost::lvalue_property_map_tag;

    Index_to_point_map(const std::vector<Point>& points) : 
    m_points(points) { 

      CGAL_precondition(m_points.size() > 0);
    }

    reference operator[](key_type index) const { 
                
      CGAL_precondition(index >= 0);
      CGAL_precondition(index < m_points.size());
      return m_points[index];
    }

    friend inline reference get(
      const Index_to_point_map& index_to_point_map, 
      key_type key) { 
      return index_to_point_map[key];
    }
                
  private:
    const std::vector<Point>& m_points;
  };

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
  void boundary_points_on_line_2(
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
    const std::vector<Point_2>& points,
    const std::vector<std::size_t>& indices,
    const Line_2& line) {

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

  template<typename Vector_3>
  void average_vector_3(
    const std::vector<Vector_3>& vectors,
    const std::vector<std::size_t>& indices,
    Vector_3& res) {

    using Traits = 
    typename Kernel_traits<Vector_3>::Kernel;
    using FT = typename Traits::FT;

    FT x = FT(0), y = FT(0), z = FT(0);
    for (std::size_t i = 0; i < indices.size(); ++i) {
                    
      const std::size_t index = indices[i];
      const Vector_3& v = vectors[index];

      x += v.x();
      y += v.y();
      z += v.z();
    }

    x /= static_cast<FT>(indices.size());
    y /= static_cast<FT>(indices.size());
    z /= static_cast<FT>(indices.size());

    res = Vector_3(x, y, z);
  }

  template<typename Vector_3>
  typename Kernel_traits<Vector_3>::Kernel::FT
  angle_deg_3(
    const Vector_3& m, 
    const Vector_3& n) {
				
    using Traits = 
    typename Kernel_traits<Vector_3>::Kernel;
    using FT = typename Traits::FT;

    typename Traits::Construct_cross_product_vector_3 cross_product_3;
    typename Traits::Compute_squared_length_3 squared_length_3;
    typename Traits::Compute_scalar_product_3 dot_product_3;

		const auto cross = cross_product_3(m, n);
			
    const FT length = 
      static_cast<FT>(
        CGAL::sqrt(
          CGAL::to_double(
            squared_length_3(cross))));

		const FT dot = dot_product_3(m, n);

		FT angle_rad = 
      static_cast<FT>(
        std::atan2(
          CGAL::to_double(length), 
          CGAL::to_double(dot)));
                
    const FT half_pi = static_cast<FT>(CGAL_PI) / FT(2);
      
    if (angle_rad > half_pi) 
      angle_rad = static_cast<FT>(CGAL_PI) - angle_rad;

		const FT angle_deg = 
      angle_rad * FT(180) / static_cast<FT>(CGAL_PI);
      
    return angle_deg;
	}

  template<
  typename Input_range, 
  typename Point_map,
  typename Point_3>
  void barycenter_3(
    const Input_range& input_range,
    const Point_map point_map,
    const std::vector<std::size_t>& indices, 
    Point_3& b) {

    using Traits = 
    typename Kernel_traits<Point_3>::Kernel;
    using FT = typename Traits::FT;

    CGAL_precondition(indices.size() > 0);

    FT x = FT(0), y = FT(0), z = FT(0);
    for (std::size_t i = 0; i < indices.size(); ++i) {    
      const Point_3& p = get(point_map, *(input_range.begin() + indices[i]));

      x += p.x();
      y += p.y();
      z += p.z();
    }

    x /= static_cast<FT>(indices.size());
    y /= static_cast<FT>(indices.size());
    z /= static_cast<FT>(indices.size());

    b = Point_3(x, y, z);
  }

  template<
  typename Input_range, 
  typename Point_map,
  typename Line_3>
  typename Kernel_traits<Line_3>::Kernel::FT
  line_from_points_3(
    const Input_range& input_range,
    const Point_map point_map,
    const std::vector<std::size_t>& indices, 
    Line_3& line) {

    using Traits = 
    typename Kernel_traits<Line_3>::Kernel;
    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;

    using Local_traits 
    = CGAL::Exact_predicates_inexact_constructions_kernel;
		using Local_FT = typename Local_traits::FT;
		using Local_line_3 = typename Local_traits::Line_3;
    using Local_point_3 = typename Local_traits::Point_3;

		using Diagonalize_traits = CGAL::Eigen_diagonalize_traits<Local_FT, 3>;

    CGAL_precondition(indices.size() >= 2);

		std::vector<Local_point_3> local_points(indices.size());
		for (std::size_t i = 0; i < indices.size(); ++i) {
				
      CGAL_precondition(indices[i] >= 0 && indices[i] < input_range.size());
			const Point_3& p = get(point_map, *(input_range.begin() + indices[i]));

      const Local_FT x = static_cast<Local_FT>(CGAL::to_double(p.x()));
      const Local_FT y = static_cast<Local_FT>(CGAL::to_double(p.y()));
      const Local_FT z = static_cast<Local_FT>(CGAL::to_double(p.z()));

			local_points[i] = Local_point_3(x, y, z);
		}

    Local_line_3 fitted_line;
    Local_point_3 fitted_centroid;

    const FT quality = static_cast<FT>(
      CGAL::linear_least_squares_fitting_3(
        local_points.begin(), local_points.end(), 
        fitted_line, fitted_centroid, CGAL::Dimension_tag<0>(),
        Local_traits(), Diagonalize_traits()));

    const Local_point_3& a = fitted_line.point(0);
    const Local_point_3& b = fitted_line.point(1);

		line = Line_3(
      Point_3(static_cast<FT>(a.x()), static_cast<FT>(a.y()), static_cast<FT>(a.z())), 
      Point_3(static_cast<FT>(b.x()), static_cast<FT>(b.y()), static_cast<FT>(b.z())));
  
    return quality;
  }

  template<
  typename Input_range, 
  typename Point_map,
  typename Line_3>
  typename Kernel_traits<Line_3>::Kernel::FT
  average_distance_to_line_3(
    const Input_range& input_range,
    const Point_map point_map,
    const std::vector<std::size_t>& indices, 
    const Line_3 &line) {

    using Traits = 
    typename Kernel_traits<Line_3>::Kernel;
    using FT = typename Traits::FT;

    typename Traits::Compute_squared_distance_3 squared_distance_3;

    CGAL_precondition(indices.size() > 0);

    FT average_distance = FT(0);
    for (std::size_t i = 0; i < indices.size(); ++i) {

      CGAL_precondition(indices[i] >= 0 && indices[i] < input_range.size());
			const auto& p = get(point_map, *(input_range.begin() + indices[i]));
      
      const auto q = line.projection(p);
      average_distance += 
      static_cast<FT>(
        CGAL::sqrt(
          CGAL::to_double(
            squared_distance_3(p, q))));
    }
                
    average_distance /= static_cast<FT>(indices.size());
    return average_distance;
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

		CGAL_precondition(items.size() >= 2);
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
  typename Items, 
  typename Point_map, 
  typename Plane_3>
	typename Kernel_traits<
  typename boost::property_traits<Point_map>::value_type>::Kernel::FT
  plane_from_points_3(
    const Items& items, 
    const Point_map point_map, 
    const std::vector<std::size_t>& indices,
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

		CGAL_precondition(indices.size() >= 2);
		std::vector<Local_point_3> local_points(indices.size());

		for (std::size_t i = 0; i < indices.size(); ++i) {
			const Point_3& point = get(point_map, *(items.begin() + indices[i]));

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
  typename Items,
  typename Point_map,
  typename Plane_3,
  typename FT>
  FT points_area_3(
    const Items& items,
    const Point_map point_map,
    const std::vector<std::size_t>& indices,
    const Plane_3& plane,
    const FT alpha) {

    using Traits = 
    typename Kernel_traits<Plane_3>::Kernel;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Triangle_2 = typename Traits::Triangle_2;
    
    using Vb = CGAL::Alpha_shape_vertex_base_2<Traits>;
    using Fb = CGAL::Alpha_shape_face_base_2<Traits>;
    using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
    using Triangulation_2 = CGAL::Delaunay_triangulation_2<Traits, Tds>;
    using Alpha_shape_2 = CGAL::Alpha_shape_2<Triangulation_2>;

    Triangulation_2 triangulation;
    for (std::size_t i = 0; i < indices.size(); ++i) {
      
      const Point_3& point = get(point_map, *(items.begin() + indices[i]));
      triangulation.insert(point_2_from_point_3(point));
    }

    FT total_area = FT(0);
    Alpha_shape_2 alpha_shape(triangulation, alpha, Alpha_shape_2::GENERAL);

    for (auto fit = alpha_shape.finite_faces_begin(); 
      fit != alpha_shape.finite_faces_end(); 
      ++fit) {

      const auto type = alpha_shape.classify(fit);
      if (type == Alpha_shape_2::INTERIOR) {

        const Point_2& p1 = fit->vertex(0)->point();
        const Point_2& p2 = fit->vertex(1)->point();
        const Point_2& p3 = fit->vertex(2)->point();

        const Triangle_2 triangle = Triangle_2(p1, p2, p3);
        total_area += triangle.area();
      }
    }

    return total_area;
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

  template<typename Point_3, typename Plane_3>
  typename Kernel_traits<Point_3>::Kernel::FT
  distance_3(
    const Point_3& p, 
    const Plane_3& plane) {
      
    using Traits = typename Kernel_traits<Point_3>::Kernel;
    using FT = typename Traits::FT;
    
    typename Traits::Compute_squared_distance_3 squared_distance_3;

    return static_cast<FT>(
        CGAL::sqrt(
          CGAL::to_double(
            squared_distance_3(p, plane))));
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

  template<typename Vector_3>
  typename Kernel_traits<Vector_3>::Kernel::FT
  cos_angle_3(
    const Vector_3& p,
    const Vector_3& q) {
    
    using Traits = typename Kernel_traits<Vector_3>::Kernel;
    typename Traits::Compute_scalar_product_3 scalar_product_3;

    return CGAL::abs(
      scalar_product_3(p, q));
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

  template<
  typename Segment_2, 
  typename Plane_3>
  typename Kernel_traits<Segment_2>::Kernel::Segment_3
  segment_3_from_segment_2_and_plane(
    const Segment_2& segment, 
    const Plane_3& plane) {

    using Traits = 
    typename Kernel_traits<Segment_2>::Kernel;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Segment_3 = typename Traits::Segment_3;

    const Point_2& source = segment.source();
    const Point_2& target = segment.target();

    const Point_3 p1 = position_on_plane_3(source, plane);
    const Point_3 p2 = position_on_plane_3(target, plane);

    return Segment_3(p1, p2);
  }

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_UTILITIES_H
