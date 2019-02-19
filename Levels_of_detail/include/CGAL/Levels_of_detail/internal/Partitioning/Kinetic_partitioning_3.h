#ifndef CGAL_LEVELS_OF_DETAIL_KINETIC_PARTITIONING_3_H
#define CGAL_LEVELS_OF_DETAIL_KINETIC_PARTITIONING_3_H

// STL includes.
#include <vector>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

// Kinetic includes.
#include "kinetic3/defs_cgal.h"

// CGAL includes.
#include <CGAL/point_generators_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/structures.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<class GeomTraits>
  class Kinetic_partitioning_3 {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;

    using Polyhedron_facet_3 = Polyhedron_facet_3<Traits>;

    using Exact_traits = CGAL::Exact_predicates_exact_constructions_kernel;
    
    using JP_polygon  = std::vector<typename Exact_traits::Point_3>;
    using JP_polygons = std::vector<JP_polygon>;

    using JP_FT = JPTD::FT;
    using JP_point_3 = JPTD::CGAL_Point_3;

    typename Traits::Construct_cross_product_vector_3 cross_product_3;
    typename Traits::Compute_squared_length_3 squared_length_3;
    typename Traits::Compute_scalar_product_3 dot_product_3;

    using Local_traits = CGAL::Exact_predicates_inexact_constructions_kernel;
		using Local_FT = Local_traits::FT;
    using Local_point_2 = Local_traits::Point_2;

    using Point_creator_2 = Creator_uniform_2<Local_FT, Local_point_2>;

    Kinetic_partitioning_3(
      const std::vector< std::vector<Point_3> >& walls,
      const std::vector< std::vector<Point_3> >& roofs,
      const std::vector<Point_3>& ground,
      const std::size_t max_intersections) :
    m_walls(walls),
    m_roofs(roofs),
    m_ground(ground),
    m_max_intersections(max_intersections),
    m_up_scale(FT(3)),
    m_down_scale(FT(1) / FT(2)),
    m_z_scale(FT(10)),
    m_tolerance(FT(1) / FT(100000)),
    m_fixed_disc_radius(FT(1) / FT(1000)),
    m_num_points_in_disc(25) { 
      
      srand(0);
    }

    void compute(std::vector<Polyhedron_facet_3>& polyhedrons) const {
      
      JP_polygons jp_polygons;

      set_ground(jp_polygons);
      set_walls(jp_polygons);
      set_roofs(jp_polygons);
    }

  private:
    const std::vector< std::vector<Point_3> >& m_walls;
    const std::vector< std::vector<Point_3> >& m_roofs;
    const std::vector<Point_3>& m_ground;

    // External parameters.
    const std::size_t m_max_intersections;

    // Internal parameters.
    const FT m_up_scale;
    const FT m_down_scale;
    const FT m_z_scale;

    const FT m_tolerance;
    const FT m_fixed_disc_radius;

    const std::size_t m_num_points_in_disc;

    void set_ground(JP_polygons& jp_polygons) const {
      
      std::vector<Point_3> polygon = m_ground;
      process_polygon(polygon, jp_polygons, m_up_scale, FT(1));
    }

    void set_walls(JP_polygons& jp_polygons) const {
      
      std::vector<Point_3> polygon;
      for (std::size_t i = 0; i < m_walls.size(); ++i) {
        const auto& wall = m_walls[i];

        polygon = wall;
        process_polygon(polygon, jp_polygons, m_down_scale, m_z_scale);
      }
    }

    void set_roofs(JP_polygons& jp_polygons) const {
      
      std::vector<Point_3> polygon;
      for (std::size_t i = 0; i < m_roofs.size(); ++i) {
        const auto& roof = m_roofs[i];

        polygon = roof;
        process_polygon(polygon, jp_polygons, m_up_scale, FT(1));
      }
    }

    void process_polygon(
      std::vector<Point_3> &polygon, 
      JP_polygons& jp_polygons, 
      const FT scale, 
      const FT z_extender) const {

      if (polygon.size() == 0) 
        return;

      scale_polygon(scale, z_extender, polygon);
      perturb_polygon_vertices(polygon);
                
      JP_polygon jp_polygon(polygon.size());

      for (std::size_t i = 0; i < polygon.size(); ++i)
        jp_polygon[i] = JP_point_3(
          JP_FT(polygon[i].x()), JP_FT(polygon[i].y()), JP_FT(polygon[i].z()));

      jp_polygons.push_back(jp_polygon);
    }

    void scale_polygon(
      const FT scale, 
      const FT z_extender, 
      std::vector<Point_3>& polygon) const {

      Point_3 b;
      compute_barycenter(polygon, b);

      for (std::size_t i = 0; i < polygon.size(); ++i) {
        Point_3& p = polygon[i];

        const FT x = (p.x() - b.x()) * scale + b.x();
        const FT y = (p.y() - b.y()) * scale + b.y();
        const FT z = (p.z() - b.z()) * scale * z_extender + b.z();

        p = Point_3(x, y, z);
      }
    }

    void compute_barycenter(
      const std::vector<Point_3>& polygon, 
      Point_3& b) const {

      CGAL_precondition(polygon.size() > 0);
      FT x = FT(0), y = FT(0), z = FT(0);

      for (std::size_t i = 0; i < polygon.size(); ++i) {
        const Point_3& p = polygon[i];

        x += p.x();
        y += p.y();
        z += p.z();
      }

      x /= static_cast<FT>(polygon.size());
      y /= static_cast<FT>(polygon.size());
      z /= static_cast<FT>(polygon.size());

      b = Point_3(x, y, z);
    }

    void perturb_polygon_vertices(std::vector<Point_3>& polygon) const {
      
      Vector_3 source_normal;
      bool success = compute_source_normal(polygon, source_normal);

      if (!success) 
        return;

      const Vector_3 target_normal = 
        Vector_3(FT(0), FT(0), FT(1));;

      FT angle; Vector_3 axis;
      success = compute_angle_and_axis(source_normal, target_normal, angle, axis);

      if (!success) 
        return;

      Point_3 b;
      compute_barycenter(polygon, b);

      const FT angle_deg = angle * FT(180) / static_cast<FT>(CGAL_PI);
                
      if (angle_deg != FT(0) && angle_deg != FT(180))
        rotate_polygon(b, angle, axis, polygon);
                
      perturb_horizontal_polygon(polygon);

      if (angle_deg != FT(0) && angle_deg != FT(180))
        rotate_polygon(b, -angle, axis, polygon);
    }

    bool compute_source_normal(
      const std::vector<Point_3>& polygon, 
      Vector_3& normal) const {

      CGAL_precondition(polygon.size() >= 3);

      const Vector_3 zero = Vector_3(FT(0), FT(0), FT(0));
      for (std::size_t i = 0; i < polygon.size(); ++i) {

        const std::size_t ip  = (i + 1) % polygon.size();
        const std::size_t ipp = (i + 2) % polygon.size();
                    
        const Point_3& p1 = polygon[i];
        const Point_3& p2 = polygon[ip];
        const Point_3& p3 = polygon[ipp];

        const Vector_3 v1 = Vector_3(p2, p1);
        const Vector_3 v2 = Vector_3(p2, p3);

        normal = cross_product_3(v1, v2);
        if (!are_equal_points(normal, zero)) {
                     
          normalize(normal);
          return true;
        }
      }
      return false;
    }

    template<class Point>
    bool are_equal_points(
      const Point& p, 
      const Point& q) const {

      const FT eps = m_tolerance;
      return 
      (CGAL::abs(p.x() - q.x()) < eps) && 
      (CGAL::abs(p.y() - q.y()) < eps) && 
      (CGAL::abs(p.z() - q.z()) < eps);
    }

    void normalize(Vector_3& v) const {
      v /= static_cast<FT>(
        CGAL::sqrt(
          CGAL::to_double(v.squared_length())));
    }

    bool compute_angle_and_axis(
      const Vector_3& m, 
      const Vector_3& n, 
      FT& angle, 
      Vector_3& axis) const {
				
			const auto cross = cross_product_3(m, n);
			
      const FT length = static_cast<FT>(
        CGAL::sqrt(
          CGAL::to_double(
            squared_length_3(cross))));

			const FT dot = dot_product_3(m, n);

			angle = static_cast<FT>(
        std::atan2(
          CGAL::to_double(length), 
          CGAL::to_double(dot)));
			
      const FT angle_deg = 
        angle * FT(180) / static_cast<FT>(CGAL_PI);

			if (angle_deg == FT(0) || angle_deg == FT(180)) 
				return true;

			if (length == FT(0)) {
        std::cout << "Error: length = 0" << std::endl;
        exit(0);
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

    void rotate_polygon(
      const Point_3& b, 
      const FT angle, 
      const Vector_3& axis, 
      std::vector<Point_3>& polygon) const {

      Point_3 q;
      for (std::size_t i = 0; i < polygon.size(); ++i) {   
        Point_3& p = polygon[i];

        q = Point_3(p.x() - b.x(), p.y() - b.y(), p.z() - b.z());
        rotate_point(angle, axis, q);
        p = Point_3(q.x() + b.x(), q.y() + b.y(), q.z() + b.z());
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

    void perturb_horizontal_polygon(std::vector<Point_3>& polygon) const {

      const FT disc_radius = m_fixed_disc_radius;
      for (std::size_t i = 0; i < polygon.size(); ++i) {

        Point_3& p = polygon[i];
        perturb_point_inside_disc(disc_radius, p);
      }
    }

    void perturb_point_inside_disc(
      const FT disc_radius, 
      Point_3& p) const {

      std::vector<Local_point_2> points_2;
      points_2.reserve(m_num_points_in_disc);

      CGAL::Random_points_in_disc_2<Local_point_2, Point_creator_2> 
        random_points_2(CGAL::to_double(disc_radius));
      
      CGAL::cpp11::copy_n(
        random_points_2, 
        m_num_points_in_disc, 
        std::back_inserter(points_2));

      const std::size_t rand_index = size_t_rand(m_num_points_in_disc - 1);
      const Local_point_2& q = points_2[rand_index];

      const FT qx = static_cast<FT>(q.x());
      const FT qy = static_cast<FT>(q.y());

      p = Point_3(p.x() + qx, p.y() + qy, p.z());
    }

    std::size_t size_t_rand(const std::size_t maxv) const {
      return static_cast<std::size_t>(rand() % maxv);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_KINETIC_PARTITIONING_3_H
