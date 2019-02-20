#ifndef CGAL_LEVELS_OF_DETAIL_WEIGHT_QUALITY_ESTIMATOR_3_H
#define CGAL_LEVELS_OF_DETAIL_WEIGHT_QUALITY_ESTIMATOR_3_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/Polygon_2.h>
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enumerations.h>
#include <CGAL/Levels_of_detail/internal/structures.h>
#include <CGAL/Levels_of_detail/internal/utilities.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Weight_quality_estimator_3 {
			
  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;

    using Polyhedron_facet_3 = Polyhedron_facet_3<Traits>;
    using Graphcut_face_3= Graphcut_face_3<Traits>;

    using Polygon = CGAL::Polygon_2<Traits>;

    typename Traits::Compute_squared_length_3 squared_length_3;
    typename Traits::Compute_scalar_product_3 dot_product_3;
		typename Traits::Construct_cross_product_vector_3 cross_product_3;

    Weight_quality_estimator_3(
      const std::vector<Polyhedron_facet_3>& polyhedrons) :
    m_polyhedrons(polyhedrons),
    m_tolerance(FT(1) / FT(100000))
    { }

    void estimate(std::vector<Graphcut_face_3>& gc_faces) const {

      gc_faces.clear();
      for (int i = 0; i < m_polyhedrons.size(); ++i)
        add_graphcut_faces(i, gc_faces);
      normalize_weights(gc_faces);
    }
    
  private:
    const std::vector<Polyhedron_facet_3>& m_polyhedrons;
    const FT m_tolerance;

    void add_graphcut_faces(
      const int poly_index, std::vector<Graphcut_face_3>& gc_faces) const {

      for (int i = 0; i < m_polyhedrons[poly_index].faces.size(); ++i)
        process_polyhedron_face(poly_index, i, gc_faces);
    }

    void process_polyhedron_face(
      const int poly_index, 
      const int face_index, 
      std::vector<Graphcut_face_3>& gc_faces) const {

      if (was_already_added(poly_index, face_index, gc_faces))
        return;

      const auto& polyhedron = m_polyhedrons[poly_index];
                    
      const auto& faces = polyhedron.faces;
      const auto& vertices = polyhedron.vertices;

      Graphcut_face_3 gc_face;
      gc_face.data = std::make_pair(poly_index, face_index);

      find_neighbors(
        faces, vertices, 
        poly_index, face_index, 
        gc_face.neighbors);
                    
      compute_weight(faces[face_index], vertices, gc_face);
      compute_quality(gc_face);

      if (gc_face.quality > FT(1) / FT(2)) 
        gc_face.is_valid = true;
      else 
        gc_face.is_valid = false;

      gc_faces.push_back(gc_face);
    }

    bool was_already_added(
      const int poly_index, 
      const int face_index, 
      const std::vector<Graphcut_face_3>& gc_faces) const {

      for (std::size_t i = 0; i < gc_faces.size(); ++i) {
        const auto& neighbors = gc_faces[i].neighbors;

        const auto& neigh_1 = neighbors.first;
        const auto& neigh_2 = neighbors.second;

        if (neigh_1.first == poly_index && neigh_1.second == face_index) 
          return true;
        
        if (neigh_2.first == poly_index && neigh_2.second == face_index) 
          return true;
      }
      return false;
    }

    void find_neighbors(
      const std::vector< std::vector<std::size_t> >& faces, 
      const std::vector<Point_3>& vertices, 
      const int poly_index, 
      const int face_index, 
      std::pair< std::pair<int, int>, std::pair<int, int> >& neighbors) const {

      auto& neigh_1 = neighbors.first;
      auto& neigh_2 = neighbors.second;

      neigh_1.first = poly_index;
      neigh_1.second = face_index;

      find_neighbor(poly_index, faces[face_index], vertices, neigh_2);
    }

    void find_neighbor(
      const int poly_index, 
      const std::vector<std::size_t>& face, 
      const std::vector<Point_3>& vertices, 
      std::pair<int, int>& neigh) const {

      neigh.first = -1;
      neigh.second = -1;

      for (std::size_t i = 0; i < m_polyhedrons.size(); ++i) {
        if (i == poly_index) continue;

        for (std::size_t j = 0; j < m_polyhedrons[i].faces.size(); ++j) {
          if (are_equal_faces(
            face, vertices, 
            m_polyhedrons[i].faces[j], 
            m_polyhedrons[i].vertices)) {

            neigh.first = i;
            neigh.second = j;

            return;
          }
        }
      }
    }

    bool are_equal_faces(
      const std::vector<std::size_t>& f1, 
      const std::vector<Point_3>& v1, 
      const std::vector<std::size_t>& f2, 
      const std::vector<Point_3>& v2) const {

      if (f1.size() != f2.size()) 
        return false;

      std::size_t count = 0;
      for (std::size_t i = 0; i < f1.size(); ++i) {
        for (std::size_t j = 0; j < f2.size(); ++j) {

          if (CGAL::abs(v1[f1[i]].x() - v2[f2[j]].x()) < m_tolerance && 
              CGAL::abs(v1[f1[i]].y() - v2[f2[j]].y()) < m_tolerance && 
              CGAL::abs(v1[f1[i]].z() - v2[f2[j]].z()) < m_tolerance) {
            
            ++count; break;
          }
        }
      }
      return count == f1.size();
    }

    void compute_weight(
      const std::vector<std::size_t>& face, 
      const std::vector<Point_3>& vertices, 
      Graphcut_face_3& gc_face) const {

      gc_face.weight = compute_face_area(face, vertices);
    }

    FT compute_face_area(
      const std::vector<std::size_t>& face, 
      const std::vector<Point_3>& vertices) const {

      Vector_3 source_normal;
      bool success = compute_source_normal(face, vertices, source_normal);

      if (!success) 
        return get_default_weight();

      const Vector_3 target_normal = 
        Vector_3(FT(0), FT(0), FT(1));

      if (source_normal == -target_normal) 
        source_normal = target_normal;

      FT angle; Vector_3 axis;
      success = compute_angle_and_axis(source_normal, target_normal, angle, axis);

      if (!success) 
        return get_default_weight();
                
      Point_3 b;
      compute_3d_polygon_barycenter(face, vertices, b);
                    
      Polygon polygon;
      success = create_polygon(face, vertices, b, angle, axis, polygon);

      if (!success) 
        return get_default_weight();
      
      return CGAL::abs(polygon.area());
    }

    FT get_default_weight() const {
      return FT(0);
    }

    bool compute_source_normal(
      const std::vector<std::size_t>& face, 
      const std::vector<Point_3>& vertices, 
      Vector_3& normal) const {
                
      CGAL_precondition(face.size() >= 3);

      const Vector_3 zero = Vector_3(FT(0), FT(0), FT(0));
      for (std::size_t i = 0; i < face.size(); ++i) {

        const std::size_t ip = (i + 1) % face.size();
        const std::size_t ipp = (i + 2) % face.size();
                    
        const Point_3& p1 = vertices[face[i]];
        const Point_3& p2 = vertices[face[ip]];
        const Point_3& p3 = vertices[face[ipp]];

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

    template<typename Point>
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
          CGAL::to_double(
            v.squared_length())));
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

    void compute_3d_polygon_barycenter(
      const std::vector<std::size_t>& face, 
      const std::vector<Point_3>& vertices, 
      Point_3& b) const {
                
      CGAL_precondition(face.size() != 0);
      FT x = FT(0), y = FT(0), z = FT(0);

      for (std::size_t i = 0; i < face.size(); ++i) {
        const Point_3& p = vertices[face[i]];

        x += p.x();
        y += p.y();
        z += p.z();
      }

      x /= static_cast<FT>(face.size());
      y /= static_cast<FT>(face.size());
      z /= static_cast<FT>(face.size());

      b = Point_3(x, y, z);
    }

    bool create_polygon(
      const std::vector<std::size_t>& face, 
      const std::vector<Point_3>& vertices, 
      const Point_3& b, 
      const FT angle, 
      const Vector_3& axis, 
      Polygon& polygon) const {

      polygon.clear();

      for (std::size_t i = 0; i < face.size(); ++i) {
        const Point_3& p = vertices[face[i]];

				const FT angle_deg = angle * FT(180) / static_cast<FT>(CGAL_PI);
        if (angle_deg != FT(0) && angle_deg != FT(180)) {

          Point_3 q = Point_3(p.x() - b.x(), p.y() - b.y(), p.z() - b.z());
          rotate_point(angle, axis, q);
          polygon.push_back(Point_2(q.x() + b.x(), q.y() + b.y()));

        } else 
          polygon.push_back(Point_2(p.x(), p.y()));
      }

      if (!polygon.is_simple()) 
        return false;

      if (polygon.is_clockwise_oriented()) 
        polygon.reverse_orientation();

      return polygon.size() >= 3;
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

    void compute_quality(Graphcut_face_3& gc_face) const {
      gc_face.quality = FT(1);
    }

    void normalize_weights(std::vector<Graphcut_face_3>& gc_faces) const {
                
			FT total_weight = FT(0);
			for (std::size_t i = 0; i < gc_faces.size(); ++i) {
					
				const auto& gc_face = gc_faces[i];
				total_weight += gc_face.weight;
			}

			for (std::size_t i = 0; i < gc_faces.size(); ++i)
				gc_faces[i].weight /= total_weight;
    }

  }; // Weight_quality_estimator_3

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_WEIGHT_QUALITY_ESTIMATOR_3_H
