#ifndef CGAL_LEVELS_OF_DETAIL_ROOF_WALL_EXTRACTOR_H
#define CGAL_LEVELS_OF_DETAIL_ROOF_WALL_EXTRACTOR_H

// STL includes.
#include <vector>
#include <algorithm>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>
#include <CGAL/Levels_of_detail/internal/structures.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Roof_wall_extractor {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    using Plane_3 = typename Traits::Plane_3;

    using Roof = Roof<Traits>;
    using Wall = Wall<Traits>;

    using Polyhedron_facet_3 = Polyhedron_facet_3<Traits>;

    typename Traits::Construct_cross_product_vector_3 cross_product_3;
    typename Traits::Compute_squared_length_3 squared_length_3;
    typename Traits::Compute_scalar_product_3 dot_product_3;

    Roof_wall_extractor(
      const std::vector<Polyhedron_facet_3>& polyhedrons,
      const Plane_3& ground_plane) :
    m_polyhedrons(polyhedrons),
    m_ground_plane(ground_plane),
    m_tolerance(FT(1) / FT(100000)),
    m_angle_threshold(FT(10)),
    m_distance_threshold(FT(1) / FT(2))
    { }

    void extract(
      std::vector<Roof>& roofs,
      std::vector<Wall>& walls) const {

      roofs.clear();
      walls.clear();

      std::vector< std::vector<Point_3> > faces;
      for (std::size_t i = 0; i < m_polyhedrons.size(); ++i)
				process_polyhedron(i, faces);
      
      create_roofs_and_walls(faces, roofs, walls);
    }

  private:
    const std::vector<Polyhedron_facet_3>& m_polyhedrons;
    const Plane_3& m_ground_plane;

    const FT m_tolerance;
    const FT m_angle_threshold;
    const FT m_distance_threshold;

    void process_polyhedron(
      const std::size_t poly_index, 
      std::vector< std::vector<Point_3> >& result) const {

      const Polyhedron_facet_3& polyhedron = m_polyhedrons[poly_index];
      std::vector<Point_3> face;

      if (polyhedron.visibility == Visibility_label::OUTSIDE) 
        return;

			const auto& faces = polyhedron.faces;
			const auto& vertices = polyhedron.vertices;

      for (const auto& f : faces) {

        face.clear();
        face.reserve(f.size());

        for (const auto i : f)
          face.push_back(vertices[i]);

        if (!is_interior_face(face, poly_index)) 
          result.push_back(face);
      }
    }

    bool is_interior_face(
      const std::vector<Point_3>& face, 
      const std::size_t poly_index) const {

      for (std::size_t i = 0; i < m_polyhedrons.size(); ++i) {
        const auto& polyhedron = m_polyhedrons[i];
                    
        if (
          polyhedron.visibility == Visibility_label::OUTSIDE || 
          i == poly_index) 
          continue;

        const auto& faces = polyhedron.faces;
				const auto& vertices = polyhedron.vertices;

        for (std::size_t j = 0; j < faces.size(); ++j)
          if (are_equal_faces(face, faces[j], vertices)) 
            return true;
      }
      return false;
    }

    bool are_equal_faces(
      const std::vector<Point_3>& f1, 
      const std::vector<std::size_t>& face, 
      const std::vector<Point_3>& vertices) const {

      if (f1.size() != face.size()) 
        return false;

      std::size_t count = 0;
      for (std::size_t i = 0; i < f1.size(); ++i) {
        for (std::size_t j = 0; j < face.size(); ++j) {
                        
          if (are_equal_points(f1[i], vertices[face[j]])) {
                            
            ++count;
            break;
          }
        }
      }
      return count == f1.size();
    }

    template<typename Point>
    bool are_equal_points(
      const Point& p, 
      const Point& q) const {
      
      return 
        CGAL::abs(p.x() - q.x()) < m_tolerance && 
        CGAL::abs(p.y() - q.y()) < m_tolerance && 
        CGAL::abs(p.z() - q.z()) < m_tolerance;
    }

    void create_roofs_and_walls(
      const std::vector< std::vector<Point_3> >& faces,
      std::vector<Roof>& roofs,
      std::vector<Wall>& walls) const {

      for (const auto& face : faces) {
        if (is_ground_face(face))
          continue;

        if (is_vertical_face(face))
          walls.push_back(Wall(face));
        else
          roofs.push_back(Roof(face));
      }
    }

    bool is_ground_face(const std::vector<Point_3>& face) const {

      for (const auto& p : face) {
        const Point_3 q = m_ground_plane.projection(p);
        if (internal::distance_3(p, q) > m_distance_threshold)
          return false;
      }
      return true;
    }

    bool is_vertical_face(const std::vector<Point_3>& face) const {

      Vector_3 face_normal;
			const bool success = set_face_normal(face, face_normal);

      if (!success) 
        return true;

			const Vector_3 ground_normal = Vector_3(FT(0), FT(0), FT(1));

      const FT angle = compute_angle(face_normal, ground_normal);
      const FT angle_diff = CGAL::abs(FT(90) - CGAL::abs(angle));

      if (angle_diff < m_angle_threshold) 
        return true;
    
      return false;
    }

    bool set_face_normal(
      const std::vector<Point_3>& face,  
      Vector_3& face_normal) const {

      CGAL_precondition(face.size() >= 3);

      const Vector_3 zero = Vector_3(FT(0), FT(0), FT(0));
      for (std::size_t i = 0; i < face.size(); ++i) {

        const std::size_t ip = (i + 1) % face.size();
        const std::size_t ipp = (i + 2) % face.size();
                    
        const Point_3& p1 = face[i];
        const Point_3& p2 = face[ip];
        const Point_3& p3 = face[ipp];

        const Vector_3 v1 = Vector_3(p2, p1);
        const Vector_3 v2 = Vector_3(p2, p3);

        face_normal = cross_product_3(v1, v2);
        if (!are_equal_points(face_normal, zero)) {
                     
          normalize(face_normal);
          return true;
        }
      }
      return false;
		}

    void normalize(Vector_3& v) const {
      v /= static_cast<FT>(
        CGAL::sqrt(
          CGAL::to_double(
            v.squared_length())));
    }

    FT compute_angle(const Vector_3& m, const Vector_3& n) const {
				
			const auto cross = cross_product_3(m, n);
			const FT length = static_cast<FT>(
        CGAL::sqrt(
          CGAL::to_double(
            squared_length_3(cross))));

			const FT dot = dot_product_3(m, n);

			FT angle_rad = static_cast<FT>(
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

  }; // Roof_wall_extractor

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_ROOF_WALL_EXTRACTOR_H
