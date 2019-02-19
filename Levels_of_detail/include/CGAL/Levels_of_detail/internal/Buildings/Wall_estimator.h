#ifndef CGAL_LEVELS_OF_DETAIL_WALL_ESTIMATOR_H
#define CGAL_LEVELS_OF_DETAIL_WALL_ESTIMATOR_H

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/number_utils.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>
#include <CGAL/Levels_of_detail/internal/structures.h>

// Shape detection.
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Coplanar_faces_connectivity.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Coplanar_faces_conditions.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Wall_estimator {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Segment_2 = typename Traits::Segment_2;
    using Plane_3 = typename Traits::Plane_3;
    using Vector_3 = typename Traits::Vector_3;

    using Connectivity = Coplanar_faces_connectivity<Traits>;
    using Conditions = Coplanar_faces_conditions<Traits>;
    using Region_growing = Region_growing<Connectivity, Conditions>;

    typename Traits::Construct_cross_product_vector_3 cross_product_3;
    typename Traits::Compute_squared_length_3 squared_length_3;
    typename Traits::Compute_scalar_product_3 dot_product_3;

    using Face_info = Face_info<Traits>;
    using Vertex_info = Vertex_info<Traits>;

    using VB = CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, Traits>;
    using FB = CGAL::Triangulation_face_base_with_info_2<Face_info, Traits>;
    
    using CFB = CGAL::Constrained_triangulation_face_base_2<Traits, FB>;
    using TAG = CGAL::Exact_predicates_tag;
    using TDS = CGAL::Triangulation_data_structure_2<VB, CFB>;

    using CDT = CGAL::Constrained_Delaunay_triangulation_2<Traits, TDS, TAG>;
    
    using Face_handle = typename CDT::Face_handle;
    using Vertex_handle = typename CDT::Vertex_handle;
    using Edge = typename CDT::Edge;

    Wall_estimator(
      const std::vector<Segment_2>& boundaries,
      const Plane_3& ground_plane,
      const FT building_height) :
    m_boundaries(boundaries),
    m_ground_plane(ground_plane),
    m_building_height(building_height),
    m_tolerance(FT(1) / FT(100000)),
    m_default_height(-FT(100000000000000)),
    m_max_num_iters(100)
    { }

    void estimate(std::vector< std::vector<Point_3> >& result) const {
      
      std::vector< std::vector<Point_3> > walls;

      walls.clear();
      walls.reserve(m_boundaries.size());

      for (std::size_t i = 0; i < m_boundaries.size(); ++i)
        estimate_wall(m_boundaries[i], walls);

      std::vector< std::vector<std::size_t> > regions;
      
      detect_coplanar_walls(walls, regions);
      merge_coplanar_walls(walls, regions, result);
    }

  private:
    const std::vector<Segment_2>& m_boundaries;
    const Plane_3& m_ground_plane;
    const FT m_building_height;

    const FT m_tolerance;
    const FT m_default_height;

    const std::size_t m_max_num_iters;

    void estimate_wall(
      const Segment_2& boundary,
      std::vector< std::vector<Point_3> >& walls) const {

      const Point_2& source = boundary.source();
      const Point_2& target = boundary.target();
      
      std::vector<Point_3> wall(4);

      wall[0] = internal::position_on_plane_3(source, m_ground_plane);
      wall[1] = internal::position_on_plane_3(target, m_ground_plane); 

      const FT height1 = m_building_height - wall[1].z();
      const FT height0 = m_building_height - wall[0].z();

      wall[2] = Point_3(wall[1].x(), wall[1].y(), wall[1].z() + height1);
      wall[3] = Point_3(wall[0].x(), wall[0].y(), wall[0].z() + height0);

      walls.push_back(wall);
    }

    void detect_coplanar_walls(
      const std::vector< std::vector<Point_3> >& walls,
      std::vector< std::vector<std::size_t> >& regions) const {

      Connectivity connectivity(walls);
      Conditions conditions(walls);

      std::vector<std::size_t> indices(walls.size());
      for (std::size_t i = 0; i < walls.size(); ++i)
        indices[i] = i;

      Region_growing region_growing(
        indices,
        connectivity,
        conditions);

      region_growing.detect(regions);
    }

    void merge_coplanar_walls(
      const std::vector< std::vector<Point_3> >& walls,
      const std::vector< std::vector<std::size_t> >& regions,
      std::vector< std::vector<Point_3> >& result) const {

      result.clear();
      result.reserve(regions.size());

      for (std::size_t i = 0; i < regions.size(); ++i)
        merge(walls, regions[i], result);
    }

    void merge(
      const std::vector< std::vector<Point_3> >& walls,
      const std::vector<std::size_t>& region,
      std::vector< std::vector<Point_3> >& result) const {

      Vector_3 source_normal;
      bool success = compute_source_normal(walls, region, source_normal);

      if (!success) {
        std::cerr << "Error: source normal!" << std::endl;
        exit(EXIT_FAILURE);
      }

      const Vector_3 target_normal = Vector_3(FT(0), FT(0), FT(1));

      if (source_normal == -target_normal) 
        source_normal = target_normal;

      FT angle; Vector_3 axis;
      success = compute_angle_and_axis(
        source_normal, target_normal, angle, axis);

      if (!success) {
        std::cerr << "Error: angle and axis!" << std::endl;   
        exit(EXIT_FAILURE); 
      }

      Point_3 b;
      compute_barycenter(walls, region, b);

      const FT angle_deg = 
        angle * FT(180) / static_cast<FT>(CGAL_PI);
                
      std::vector< std::vector<Point_3> > rotated;
      if (angle_deg != FT(0) && angle_deg != FT(180))
        rotate_walls(walls, region, angle, axis, b, rotated);
      
      CDT cdt;
      triangulate_region_facets(rotated, cdt);

      std::vector<Point_3> final_face;
      if (cdt.number_of_faces() != 0) 
        get_back_region_facets(cdt, final_face);

      fix_orientation(final_face);

      std::vector< std::vector<Point_3> > region_facets = { final_face };
      std::vector<std::size_t> indices = { 0 };

      rotated.clear();
      if (angle_deg != FT(0) && angle_deg != FT(180))
        rotate_walls(region_facets, indices, -angle, axis, b, rotated);

      result.push_back(rotated[0]);
    }

    bool compute_source_normal(
      const std::vector< std::vector<Point_3> >& walls,
      const std::vector<std::size_t>& region, 
      Vector_3& normal) const {
                
      if (region.size() == 0) 
        return false;

      Vector_3 tmp_normal; 
      FT x = FT(0), y = FT(0), z = FT(0);

      // all walls should be equally oriented!
      for (std::size_t i = 0; i < 1; ++i) {
        const bool success = compute_wall_normal(walls[region[i]], tmp_normal);

        if (!success)
          return false;

        x += tmp_normal.x();
        y += tmp_normal.y();
        z += tmp_normal.z();
      }

      // x /= static_cast<FT>(region.size());
      // y /= static_cast<FT>(region.size());
      // z /= static_cast<FT>(region.size());

      normal = Vector_3(x, y, z);
      const Vector_3 zero = Vector_3(FT(0), FT(0), FT(0));

      if (are_equal_points(normal, zero))
        return false;

      normalize(normal);
      return true;
    }

    bool compute_wall_normal(
      const std::vector<Point_3>& wall, 
      Vector_3& normal) const {
                
      CGAL_precondition(wall.size() >= 3);
      if (wall.size() < 3)
        return false;

      const bool success = compute_cross_product(wall, normal);
      if (success) {
                    
        normalize(normal);
        return true;
      }
      return false;
    }

    bool compute_cross_product(
      const std::vector<Point_3>& wall, 
      Vector_3& normal) const {

      const Vector_3 zero = Vector_3(FT(0), FT(0), FT(0));
      for (std::size_t i = 0; i < wall.size(); ++i) {

        const std::size_t ip = (i + 1) % wall.size();
        const std::size_t ipp = (i + 2) % wall.size();

        const Point_3& p1 = wall[i];
        const Point_3& p2 = wall[ip];
        const Point_3& p3 = wall[ipp];

        const Vector_3 v1 = Vector_3(p2, p1);
        const Vector_3 v2 = Vector_3(p2, p3);

        normal = cross_product_3(v1, v2);
        if (!are_equal_points(normal, zero)) 
          return true;
      }
      return false;
    }

    template<class Point>
    bool are_equal_points(const Point& p, const Point& q) const {

      const FT eps = m_tolerance;
      return 
      (CGAL::abs(p.x() - q.x()) < eps) && 
      (CGAL::abs(p.y() - q.y()) < eps) && 
      (CGAL::abs(p.z() - q.z()) < eps);
    }

    void normalize(Vector_3& v) const {
      v /= static_cast<FT>(
        CGAL::sqrt(
          CGAL::to_double(
            v.squared_length())));
    }

    bool compute_angle_and_axis(
      const Vector_3& m, const Vector_3& n, 
      FT& angle, Vector_3& axis) const {

			const auto cross = cross_product_3(m, n);
			const FT length = 
        static_cast<FT>(
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
                 
        std::cerr << "Error: length = 0" << std::endl;
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

    void compute_barycenter(
      const std::vector< std::vector<Point_3> >& walls, 
      const std::vector<std::size_t>& region,
      Point_3& b) const {

      CGAL_precondition(region.size() > 0);

      Point_3 tmp_b; FT x = FT(0), y = FT(0), z = FT(0);
      for (std::size_t i = 0; i < region.size(); ++i) {
                    
        const std::vector<Point_3>& wall = walls[region[i]];
        compute_wall_barycenter(wall, tmp_b);

        x += tmp_b.x();
        y += tmp_b.y();
        z += tmp_b.z();
      }

      x /= static_cast<FT>(region.size());
      y /= static_cast<FT>(region.size());
      z /= static_cast<FT>(region.size());

      b = Point_3(x, y, z);
    }

    void compute_wall_barycenter(
      const std::vector<Point_3>& wall, 
      Point_3& b) const {
                
      CGAL_precondition(wall.size() != 0);
      FT x = FT(0), y = FT(0), z = FT(0);

      for (std::size_t i = 0; i < wall.size(); ++i) {
        const Point_3& p = wall[i];

        x += p.x();
        y += p.y();
        z += p.z();
      }

      x /= static_cast<FT>(wall.size());
      y /= static_cast<FT>(wall.size());
      z /= static_cast<FT>(wall.size());

      b = Point_3(x, y, z);
    }

    void rotate_walls(
      const std::vector< std::vector<Point_3> >& walls,
      const std::vector<std::size_t>& region, 
      const FT angle, 
      const Vector_3& axis, 
      const Point_3& b,
      std::vector< std::vector<Point_3> >& rotated) const {

      rotated.resize(region.size());
      for (std::size_t i = 0; i < region.size(); ++i)
        rotate_wall(walls[region[i]], angle, axis, b, rotated[i]);
    }

    void rotate_wall(
      const std::vector<Point_3>& wall, 
      const FT angle, 
      const Vector_3& axis, 
      const Point_3& b,
      std::vector<Point_3>& rotated) const {

      if (angle == FT(0)) {
        
        rotated = wall;
        return;
      }

      rotated.clear();
      rotated.resize(wall.size());

      Point_3 q;
      for (std::size_t i = 0; i < wall.size(); ++i) {
        const Point_3& p = wall[i];

        q = Point_3(p.x() - b.x(), p.y() - b.y(), p.z() - b.z());
        rotate_point(angle, axis, q);
        rotated[i] = Point_3(q.x() + b.x(), q.y() + b.y(), q.z() + b.z());
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

    void triangulate_region_facets(
      const std::vector< std::vector<Point_3> >& region_facets, 
      CDT& cdt) const {

			std::vector< std::vector<Vertex_handle> > vhs;
      insert_points(region_facets, cdt, vhs);
                
      std::vector< std::pair<Vertex_handle, Vertex_handle> > final_vhs;
      update_constraints(region_facets, vhs, final_vhs);

      insert_constraints(final_vhs, cdt);
    }

    void insert_points(
      const std::vector< std::vector<Point_3> >& region_facets, 
      CDT& cdt, 
      std::vector< std::vector<Vertex_handle> >& vhs) const {
                
      cdt.clear();
      vhs.clear();

      vhs.resize(region_facets.size());
			for (std::size_t i = 0; i < region_facets.size(); ++i) {
				const auto& region_facet = region_facets[i];

				vhs[i].resize(region_facet.size());
				for (std::size_t j = 0; j < region_facet.size(); ++j) {
          const Point_3& p = region_facet[j];

					vhs[i][j] = cdt.insert(Point_2(p.x(), p.y()));
					vhs[i][j]->info().height = p.z();
				}
			}
    }

    void update_constraints(
      const std::vector< std::vector<Point_3> >& region_facets, 
      const std::vector< std::vector<Vertex_handle> >& vhs, 
      std::vector< std::pair<Vertex_handle, Vertex_handle> >& final_vhs) const {

      for (std::size_t i = 0; i < region_facets.size(); ++i) {

        for (std::size_t j = 0; j < region_facets[i].size(); ++j) {
          const std::size_t jp = (j + 1) % region_facets[i].size();

          if (is_boundary_edge(region_facets[i][j], region_facets[i][jp], i, region_facets)) {

            const auto final_constraint = std::make_pair(vhs[i][j], vhs[i][jp]);
            final_vhs.push_back(final_constraint);
          }
        }
      }
    }

    bool is_boundary_edge(
      const Point_3& p1, const Point_3& p2, 
      const std::size_t facet_index, 
      const std::vector< std::vector<Point_3> >& region_facets) const {

      for (std::size_t i = 0; i < region_facets.size(); ++i) {
        if (i == facet_index) 
          continue;

        for (std::size_t j = 0; j < region_facets[i].size(); ++j) {
          const std::size_t jp = (j + 1) % region_facets[i].size();

          if (are_equal_edges(p1, p2, region_facets[i][j], region_facets[i][jp])) 
            return false;
        }
      }
      return true;
    }

    bool are_equal_edges(
      const Point_3& p1, const Point_3& p2, 
      const Point_3& q1, const Point_3& q2) const {
      
      return 
      (are_equal_points(p1, q1) && are_equal_points(p2, q2)) || 
      (are_equal_points(p1, q2) && are_equal_points(p2, q1));
    }

    void insert_constraints(
      const std::vector< std::pair<Vertex_handle, Vertex_handle> >& final_vhs, 
      CDT& cdt) const {
                
      for (std::size_t i = 0; i < final_vhs.size(); ++i) {
        const auto& final_constraint = final_vhs[i];
                    
        if (final_constraint.first != final_constraint.second)
          cdt.insert_constraint(final_constraint.first, final_constraint.second);
      }
    }

    void get_back_region_facets(
      const CDT& cdt, 
      std::vector<Point_3>& region_facet) const {

      region_facet.clear();
      if (cdt.number_of_faces() == 0) return;

      Face_handle fh;
      bool success = find_first_face_handle(cdt, fh);
      if (!success) 
        return;

      success = traverse_cdt(fh, cdt, region_facet);
      if (!success) 
        return;

      if (region_facet.size() < 3) 
        return;
    }

    bool find_first_face_handle(
      const CDT& cdt, 
      Face_handle& fh) const {

      for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
        fh = static_cast<Face_handle>(fit);

        const Vertex_handle& vh1 = fh->vertex(0);
        const Vertex_handle& vh2 = fh->vertex(1);
        const Vertex_handle& vh3 = fh->vertex(2);

        const Point_2& p1 = vh1->point();
        const Point_2& p2 = vh2->point();
        const Point_2& p3 = vh3->point();

        for (std::size_t i = 0; i < 3; ++i) {
                        
          const Edge edge = std::make_pair(fh, i);
          if (cdt.is_constrained(edge)) 
            return true;
        }
      }
      return false;
    }

    bool traverse_cdt(
      const Face_handle& fh, 
      const CDT& cdt, 
      std::vector<Point_3>& region_facet) const {
                
      Edge edge;
      region_facet.clear();

      const bool success = find_first_edge(fh, cdt, edge);
      if (!success) 
        return false;

      CGAL_precondition(edge.second >= 0 && edge.second <= 2);

      Vertex_handle vh = edge.first->vertex((edge.second + 2) % 3);
      Vertex_handle end = vh;
                
      if (vh->info().height == m_default_height) 
        return false;

      const Point_2& p = vh->point();
      region_facet.push_back(Point_3(p.x(), p.y(), vh->info().height));

      std::size_t num_iters = 0; 
      do {
        
        get_next_vertex_handle(cdt, vh, edge);
        const Point_2& q = vh->point();

        if (vh->info().height == m_default_height) 
          return false;
        
        if (vh == end) 
          break;

        region_facet.push_back(Point_3(q.x(), q.y(), vh->info().height));
                    
        if (num_iters == m_max_num_iters) 
          return false;
        ++num_iters;

      } while (vh != end);

      return is_valid_traversal(region_facet);
    }

    bool find_first_edge(
      const Face_handle& fh, 
      const CDT& cdt, 
      Edge& edge) const {

      for (int i = 0; i < 3; ++i) {              
        edge = std::make_pair(fh, i);

        if (cdt.is_constrained(edge)) 
          return true;
      }
      return false;
    }

    void get_next_vertex_handle(
      const CDT& cdt, 
      Vertex_handle& vh, 
      Edge& edge) const {

      const int index = edge.first->index(vh);
      Edge next = std::make_pair(edge.first, (index + 2) % 3);

      while (!cdt.is_constrained(next)) {

        const Face_handle fh = next.first->neighbor(next.second);
        const Vertex_handle tmp = next.first->vertex((next.second + 1) % 3);
                    
        const std::size_t tmp_index = fh->index(tmp);
        next = std::make_pair(fh, (tmp_index + 2) % 3);
      }

      vh = next.first->vertex((next.second + 2) % 3);
      edge = next;
    }

    bool is_valid_traversal(
      const std::vector<Point_3>& region_facet) const {

      if (region_facet.size() < 3) 
        return false;

      for (std::size_t i = 0; i < region_facet.size(); ++i) {
        const Point_3& p = region_facet[i];

        for (std::size_t j = 0; j < region_facet.size(); ++j) {
          const Point_3& q = region_facet[j];

          if (i == j) 
            continue;
          
          if (are_equal_points(p, q)) 
            return false;
        }
      }
      return true;
    }

    void fix_orientation(
      std::vector<Point_3>& region_facet) const {

      std::vector<Point_2> polygon(region_facet.size());
      for (std::size_t i = 0; i < region_facet.size(); ++i) {
                        
        const Point_3& p = region_facet[i];
        polygon[i] = Point_2(p.x(), p.y());
      }

      if (CGAL::orientation_2(
        polygon.begin(), polygon.end()) == CGAL::CLOCKWISE) 
        std::reverse(region_facet.begin(), region_facet.end());
    }

  }; // Wall_estimator

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_WALL_ESTIMATOR_H
