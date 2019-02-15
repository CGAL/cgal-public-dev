#ifndef CGAL_LEVELS_OF_DETAIL_ROOF_ESTIMATOR_H
#define CGAL_LEVELS_OF_DETAIL_ROOF_ESTIMATOR_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/number_utils.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits, 
  typename InputRange,
  typename PointMap>
  class Roof_estimator {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Vector_2 = typename Traits::Vector_2;
    using Vector_3 = typename Traits::Vector_3;
    using Plane_3 = typename Traits::Plane_3;
    using Line_2 = typename Traits::Line_2;

		typename Traits::Compute_squared_length_3 squared_length_3;
			
		typename Traits::Compute_scalar_product_2 dot_product_2;
		typename Traits::Compute_scalar_product_3 dot_product_3;

		typename Traits::Compute_determinant_2 cross_product_2;
		typename Traits::Construct_cross_product_vector_3 cross_product_3;

    Roof_estimator(
      const Input_range& input_range, 
      const Point_map point_map,
      const std::vector< std::vector<std::size_t> >& roof_indices) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_roof_indices(roof_indices)
    { }

    void estimate(std::vector< std::vector<Point_3> >& roofs) const {
      
      roofs.clear();
      for (std::size_t i = 0; i < m_roof_indices.size(); ++i)
        estimate_roof(m_roof_indices[i], roofs);
    }

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const std::vector< std::vector<std::size_t> >& m_roof_indices;

    void estimate_roof(
      const std::vector<std::size_t>& indices,
      std::vector< std::vector<Point_3> >& roofs) const {

      Plane_3 plane;
      internal::plane_from_points_3(
        m_input_range, m_point_map, indices, plane);

      std::vector<Point_3> points;
      project_points_on_plane(indices, plane, points);

      CGAL_precondition(points.size() >= 2);

      const Vector_3 m = plane.orthogonal_vector();
      const Vector_3 n = Vector_3(FT(0), FT(0), FT(1));

      FT angle_3d; Vector_3 axis;
		  const bool success = 
        compute_angle_and_axis(m, n, angle_3d, axis);

      if (!success) 
        return;
			
      const FT angle_3d_deg = 
        angle_3d * FT(180) / static_cast<FT>(CGAL_PI);

			if (angle_3d_deg != FT(0) && angle_3d_deg != FT(180))
				rotate_points(angle_3d, axis, points);

      Vector_2 roof_direction;
      estimate_roof_direction(points, roof_direction);

      const Vector_2 y_direction = Vector_2(FT(0), FT(1));

      FT angle_2d;
			compute_angle(roof_direction, y_direction, angle_2d);

			Point_2 barycenter;
      compute_barycenter(points, barycenter);

      rotate_points(angle_2d, barycenter, points);

			std::vector<Point_3> boundary;
			compute_bounding_box(points, boundary);

      rotate_points(-angle_2d, barycenter, boundary);

      if (angle_3d_deg != FT(0) && angle_3d_deg != FT(180))
				rotate_points(-angle_3d, axis, boundary);
				
			if (!is_valid_boundary(boundary)) {
				
				boundary.clear();
				return;
			}
			roofs.push_back(boundary);
    }

    void project_points_on_plane(
      const std::vector<std::size_t>& indices, 
      const Plane_3& plane, 
      std::vector<Point_3>& points) const {
      
      CGAL_precondition(indices.size() > 0);
      points.reserve(indices.size());

      for (std::size_t i = 0; i < indices.size(); ++i) {			
        const Point_3& p = 
          get(m_point_map, *(m_input_range.begin() + indices[i]));

				points.push_back(plane.projection(p));
      }
    }

    bool compute_angle_and_axis(
      const Vector_3& m, 
      const Vector_3& n, 
      FT& angle, 
      Vector_3& axis) const {
				
			const auto cross = cross_product_3(m, n);
			const FT length = 
        static_cast<FT>(
          CGAL::sqrt(
            CGAL::to_double(
              squared_length_3(cross))));

			const FT dot = dot_product_3(m, n);

			angle = 
        static_cast<FT>(
          std::atan2(
            CGAL::to_double(length), 
            CGAL::to_double(dot)));

			const FT angle_deg = 
        angle * FT(180) / static_cast<FT>(CGAL_PI);

			if (angle_deg == FT(0) || angle_deg == FT(180)) 
				return true;

			if (length == FT(0)) {           
        std::cerr << "error weight quality estimator: length = 0" << std::endl;
        exit(EXIT_FAILURE);
      }
                
			CGAL_precondition(length > FT(0));
			axis = cross / length;

      const FT half_pi = static_cast<FT>(CGAL_PI) / FT(2);
      if (angle > half_pi) {
                    
        angle = static_cast<FT>(CGAL_PI) - angle;
        axis = -axis;
      }
			return true;
		}

		void rotate_points(
      const FT angle, 
      const Vector_3& axis, 
      std::vector<Point_3>& points) const {

			for (std::size_t i = 0; i < points.size(); ++i) {
					
				Point_3& p = points[i];
				rotate_point(angle, axis, p);
			}
		}

		void rotate_point(
      const FT angle, 
      const Vector_3& axis, 
      Point_3& p) const {

			const double tmp_angle = CGAL::to_double(angle);

			const FT c = static_cast<FT>(std::cos(tmp_angle));
			const FT s = static_cast<FT>(std::sin(tmp_angle));

			const FT C = FT(1) - c;

			const FT x = axis.x();
			const FT y = axis.y();
			const FT z = axis.z();

			p = Point_3(
        (x * x * C + c)     * p.x() + (x * y * C - z * s) * p.y() + (x * z * C + y * s) * p.z(),
				(y * x * C + z * s) * p.x() + (y * y * C + c)     * p.y() + (y * z * C - x * s) * p.z(),
				(z * x * C - y * s) * p.x() + (z * y * C + x * s) * p.y() + (z * z * C + c)     * p.z());
		}

    void estimate_roof_direction(
      const std::vector<Point_3>& points, 
      Vector_2& direction) const {

      using Local_traits = CGAL::Exact_predicates_inexact_constructions_kernel;

			using Local_FT = typename Local_traits::FT;
      using Local_point_2 = typename Local_traits::Point_2;
      using Local_line_2 = typename Local_traits::Line_2;

      using Diagonalize_traits = CGAL::Eigen_diagonalize_traits<Local_FT, 2>;

      CGAL_precondition(points.size() > 1);

			std::vector<Local_point_2> tmp_points(points.size());
			for (std::size_t i = 0; i < points.size(); ++i) {
				const Point_3& p = points[i];

				const Local_FT x = static_cast<Local_FT>(CGAL::to_double(p.x()));
				const Local_FT y = static_cast<Local_FT>(CGAL::to_double(p.y()));

				tmp_points[i] = Local_point_2(x, y);
			}

      Local_line_2 tmp_line;
			Local_point_2 centroid;
          
			CGAL::linear_least_squares_fitting_2(
        tmp_points.begin(), tmp_points.end(), 
        tmp_line, centroid, CGAL::Dimension_tag<0>(), 
        Local_traits(), Diagonalize_traits());

			const Line_2 line = Line_2(
        Point_2(static_cast<FT>(tmp_line.point(0).x()), static_cast<FT>(tmp_line.point(0).y())), 
				Point_2(static_cast<FT>(tmp_line.point(1).x()), static_cast<FT>(tmp_line.point(1).y())));

      direction = line.to_vector();
    }

    void compute_angle(
      const Vector_2& m, 
      const Vector_2& n, 
      FT& angle) const {
				
			const FT cross = cross_product_2(m, n);
			const FT dot = dot_product_2(m, n);

			angle = 
        static_cast<FT>(
          std::atan2(
            CGAL::to_double(cross), 
            CGAL::to_double(dot)));
		}

		void compute_barycenter(
      const std::vector<Point_3>& points, 
      Point_2& b) const {

			FT bx = FT(0), by = FT(0);
			for (std::size_t i = 0; i < points.size(); ++i) {
							
				const Point_3& p = points[i];

				bx += p.x();
				by += p.y();
			}

			bx /= static_cast<FT>(points.size());
      by /= static_cast<FT>(points.size());

      b = Point_2(bx, by);
    }

		void rotate_points(
      const FT angle, 
      const Point_2& barycenter, 
      std::vector<Point_3>& points) const {

			for (std::size_t i = 0; i < points.size(); ++i) {
					
				Point_3& p = points[i];
				rotate_point(angle, barycenter, p);
			}
		}

    void rotate_point(
      const FT angle, 
      const Point_2& barycenter, 
      Point_3& p) const {

			FT x = p.x();
			FT y = p.y();

			x -= barycenter.x();
			y -= barycenter.y();

			p = Point_3(x, y, p.z());

      const double tmp_angle = CGAL::to_double(angle);

      const FT c = static_cast<FT>(std::cos(tmp_angle));
			const FT s = static_cast<FT>(std::sin(tmp_angle));

			x = p.x() * c - p.y() * s;
			y = p.y() * c + p.x() * s;

			x += barycenter.x();
			y += barycenter.y();

			p = Point_3(x, y, p.z());
		} 

		void compute_bounding_box(
      const std::vector<Point_3>& points, 
      std::vector<Point_3>& bbox) const {

      const FT big_value = FT(100000000000000);

			FT minx =  big_value, miny =  big_value;
			FT maxx = -big_value, maxy = -big_value;

			FT z = FT(0);
			for (std::size_t i = 0; i < points.size(); ++i) {
				const Point_3& p = points[i];

				minx = CGAL::min(minx, p.x());
				miny = CGAL::min(miny, p.y());

				maxx = CGAL::max(maxx, p.x());
			  maxy = CGAL::max(maxy, p.y());

				z += p.z();
			}
			z /= static_cast<FT>(points.size());

      bbox.clear();
			bbox.resize(4);

      bbox[0] = Point_3(minx, miny, z);
			bbox[1] = Point_3(maxx, miny, z);
			bbox[2] = Point_3(maxx, maxy, z);
			bbox[3] = Point_3(minx, maxy, z);
		}

		bool is_valid_boundary(const std::vector<Point_3>& boundary) const {

			if (std::isnan(CGAL::to_double(boundary[0].x())) ||
					std::isnan(CGAL::to_double(boundary[0].y())) ||
					std::isnan(CGAL::to_double(boundary[0].z()))  ) return false;

      if (std::isnan(CGAL::to_double(boundary[1].x())) ||
					std::isnan(CGAL::to_double(boundary[1].y())) ||
					std::isnan(CGAL::to_double(boundary[1].z()))  ) return false;

      if (std::isnan(CGAL::to_double(boundary[2].x())) ||
					std::isnan(CGAL::to_double(boundary[2].y())) ||
					std::isnan(CGAL::to_double(boundary[2].z()))  ) return false;

      if (std::isnan(CGAL::to_double(boundary[3].x())) ||
					std::isnan(CGAL::to_double(boundary[3].y())) ||
					std::isnan(CGAL::to_double(boundary[3].z()))  ) return false;

			if (boundary.size() < 3) 
        return false;
			
      return true;
		}

  }; // Roof_estimator

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_ROOF_ESTIMATOR_H
