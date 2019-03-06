#ifndef CGAL_LEVELS_OF_DETAIL_SMOOTH_GROUND_ESTIMATOR_H
#define CGAL_LEVELS_OF_DETAIL_SMOOTH_GROUND_ESTIMATOR_H

// STL includes.
#include <map>
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enumerations.h>
#include <CGAL/Levels_of_detail/internal/structures.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename KnnSearch,
  typename PointRange,
  typename PointMap>
  class Smooth_ground_estimator : public Planar_ground_estimator<GeomTraits> {

  public:
    using Traits = GeomTraits;
    using Knn_search = KnnSearch;
    using Point_range = PointRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Plane_3 = typename Traits::Plane_3;
    using Triangle_3 = typename Traits::Triangle_3;

    using Base = Planar_ground_estimator<Traits>;
    using Vertex_handle = typename Base::Triangulation::Vertex_handle;
    using Face_handle = typename Base::Triangulation::Face_handle;
    using Point_iterator = typename Point_range::const_iterator;

    struct Candidate_face {
      Face_handle face;
      FT max_error;
      std::vector<Point_iterator> inliers;

      Candidate_face(Face_handle face = Face_handle()) : 
      face(face), 
      max_error(FT(0))
      { }
    };

    using Candidate_face_ptr = std::shared_ptr<Candidate_face>;

    struct Compare_candidates {
      
      bool operator()(
        const Candidate_face_ptr& a, 
        const Candidate_face_ptr& b) const {
        
        if (a->max_error != b->max_error)
          return a->max_error > b->max_error;
        return a->face < b->face;
      }
    };

    using Face_map = std::map<Face_handle, Candidate_face_ptr>;
    using Face_queue = std::set<Candidate_face_ptr, Compare_candidates>;

    Smooth_ground_estimator(
      const Plane_3& ground_plane,
      const Knn_search& knn_search,
      const Point_range& points,
      const Point_map point_map,
      const FT ground_precision) :
    Base(ground_plane),
    m_ground_plane(ground_plane),
    m_knn_search(knn_search),
    m_points(points),
    m_point_map(point_map),
    m_ground_precision(ground_precision)
    { }

    template<typename Urban_object>
    void tag_faces(const Urban_object& object) {

      // Tag all faces that belong to this object. 
      for (auto fh = Base::m_triangulation.finite_faces_begin();
      fh != Base::m_triangulation.finite_faces_end(); ++fh) {
        if (Base::is_valid(fh, object.footprint())) {

          fh->info().urban_tag = object.urban_tag();
          fh->info().object_index = object.index();
          fh->info().tagged = true;
          set_fixed_ground_heights(fh);
        }
      }
    }

    void finilize() {

      // Set real heights of all ground faces.
      for (auto fh = Base::m_triangulation.finite_faces_begin();
      fh != Base::m_triangulation.finite_faces_end(); ++fh) {
        
        if(fh->info().urban_tag == Urban_object_type::GROUND)
          set_real_ground_heights(fh);
      }

      // Refine mesh until we reach the m_ground_precision.
      refine();
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_as_triangle_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces,
      std::size_t& num_vertices,
      internal::Indexer<Point_3>& indexer) {

      for (auto fit = Base::m_triangulation.finite_faces_begin(); 
      fit != Base::m_triangulation.finite_faces_end(); ++fit) {

        cpp11::array<std::size_t, 3> face;
        for (std::size_t k = 0; k < 3; ++k) {
          const Point_3 p = get_point_3(fit->vertex(k));

          const std::size_t idx = indexer(p);
          if (idx == num_vertices) {

            *(output_vertices++) = p;
            ++num_vertices;
          }
          face[k] = idx;
        }
        *(output_faces++) = 
        std::make_pair(face, fit->info().urban_tag);
      }
    }

  private:
    const Plane_3& m_ground_plane;
    const Knn_search& m_knn_search;
    const Point_range& m_points;
    const Point_map m_point_map;
    const FT m_ground_precision;

    Point_3 get_point_3(const Vertex_handle& vh) const {

      auto fh = Base::m_triangulation.incident_faces(vh);
      const auto start = fh;

      FT sum = FT(0); FT num_faces = FT(0);
      do {
        if (Base::m_triangulation.is_infinite(fh)) {
          ++fh; continue;
        }

        sum += fh->info().z[fh->index(vh)];
        num_faces += FT(1);

        ++fh;
      } while (fh != start);

      CGAL_precondition(num_faces > FT(0));
      const FT z = sum / num_faces;

      const Point_2& p = vh->point();
      return Point_3(p.x(), p.y(), z);
    }

    template<typename Face>
    void set_fixed_ground_heights(Face& face) const {
      
      Point_2 b;
      Base::barycenter(face, b);
      const FT z = internal::position_on_plane_3(b, m_ground_plane).z();
      for (std::size_t i = 0; i < 3; ++i)
        face->info().z[i] = z;
    }

    template<typename Face>
    void set_real_ground_heights(Face& face) const {

      typename Knn_search::Neighbors neighbors;
      for (std::size_t i = 0; i < 3; ++i) {
        
        const Point_2& p = face->vertex(i)->point();
        m_knn_search.get_neighbors(p, neighbors);

        FT z = FT(0);
        for (const auto& it : neighbors)
          z += get(m_point_map, *it).z();
        z /= static_cast<FT>(neighbors.size());

        face->info().z[i] = z;
      }
    }

    void refine() {

      const FT tolerance = m_ground_precision;

      Face_map face_map;
      Face_queue todo;
      init_queue(face_map, todo, m_ground_precision);

      Face_handle hint;
      while (!todo.empty()) {
        
        const Candidate_face_ptr candidate = *(todo.begin());
        todo.erase(todo.begin());

        const bool out_of_tolerance = (candidate->max_error > tolerance * tolerance);
        const bool badly_shaped = (!well_shaped(candidate->face));

        if (!out_of_tolerance && !badly_shaped)
          continue;

        if (is_too_small(candidate->face, FT(3) * tolerance))
          continue;

        // Get circumcenter and conflict zone.
        Point_2 center = CGAL::circumcenter(
          candidate->face->vertex(0)->point(),
          candidate->face->vertex(1)->point(),
          candidate->face->vertex(2)->point());

        typename Base::Triangulation::Locate_type lt; int li;
        hint = Base::m_triangulation.locate(center, lt, li, hint);
        
        if (lt == Base::Triangulation::VERTEX ||
            Base::m_triangulation.is_infinite(hint) ||
            hint->info().urban_tag == Urban_object_type::BUILDING ||
            hint->info().urban_tag == Urban_object_type::TREE) {
          
          if (out_of_tolerance) {
            center = CGAL::barycenter(
              candidate->face->vertex(0)->point(), FT(1),
              candidate->face->vertex(1)->point(), FT(1),
              candidate->face->vertex(2)->point(), FT(1));

            hint = Base::m_triangulation.locate(center, lt, li, hint);
        
            if (lt == Base::Triangulation::VERTEX ||
                Base::m_triangulation.is_infinite(hint) ||
                hint->info().urban_tag == Urban_object_type::BUILDING ||
                hint->info().urban_tag == Urban_object_type::TREE)
              continue;
          } else {
            continue;
          }
        }

        std::vector<Face_handle> conflict;
        Base::m_triangulation.get_conflicts(
          center, std::back_inserter(conflict));

        // Recover points and remove conflict cells from local structures.
        std::vector<Point_iterator> points;
        for (std::size_t i = 0; i < conflict.size(); ++i) {
          
          if (Base::m_triangulation.is_infinite(conflict[i]) || 
              conflict[i]->info().urban_tag == Urban_object_type::BUILDING || 
              conflict[i]->info().urban_tag == Urban_object_type::TREE)
            continue;

          const auto filter = face_map.find(conflict[i]);
          if (filter == face_map.end())
            continue;

          const Candidate_face_ptr cface = filter->second;
          if (!(cface->inliers.empty()))
            std::copy(
              cface->inliers.begin(), cface->inliers.end(),
              std::back_inserter(points));

          face_map.erase(filter);
          todo.erase(cface);
        }

        // Insert new vertex.
        const Vertex_handle v = Base::m_triangulation.insert(center, hint);
        std::vector<Candidate_face_ptr> new_faces;
        
        auto circ = Base::m_triangulation.incident_faces(v);
        const auto start = circ;
        do {
          
          circ->info().urban_tag = Urban_object_type::GROUND;
          set_real_ground_heights(circ);
          const Candidate_face_ptr cface = std::make_shared<Candidate_face>(circ);
          face_map.insert(std::make_pair(circ, cface));
          new_faces.push_back(cface);
          ++circ;
        }
        while (circ != start);

        // Redistribute points.
        for (std::size_t i = 0; i < points.size(); ++i) {
          
          const Point_3& point_3 = get(m_point_map, *(points[i]));
          const Point_2 point_2 = internal::point_2_from_point_3(point_3);
        
          hint = Base::m_triangulation.locate(point_2, hint);

          const auto filter = face_map.find(hint);
          CGAL_assertion(filter != face_map.end());
          Candidate_face_ptr cface = filter->second;
          
          const Triangle_3 triangle = internal::triangle_3<Triangle_3>(hint);
          const FT sq_dist = CGAL::squared_distance(point_3, triangle);

          cface->inliers.push_back(points[i]);
          cface->max_error = CGAL::max(cface->max_error, sq_dist);
        }

        // Insert new faces.
        for (std::size_t i = 0; i < new_faces.size(); ++i)
          todo.insert(new_faces[i]);
      }
    }

    void init_queue(
      Face_map& face_map,
      Face_queue& todo,
      FT tolerance) {

      for(auto fit = Base::m_triangulation.finite_faces_begin();
      fit != Base::m_triangulation.finite_faces_end(); ++fit) {

        if (Base::m_triangulation.is_infinite(fit) || 
            fit->info().urban_tag == Urban_object_type::BUILDING ||
            fit->info().urban_tag == Urban_object_type::TREE)
          continue;
        
        face_map.insert(
          std::make_pair(
            fit, 
            std::make_shared<Candidate_face>(fit)));
      }

      Face_handle hint;
      for (auto pit = m_points.begin(); pit != m_points.end(); ++pit) {
        
        const Point_3& point_3 = get(m_point_map, *pit);
        const Point_2 point_2 = internal::point_2_from_point_3(point_3);
        
        hint = Base::m_triangulation.locate(point_2, hint);
        if (Base::m_triangulation.is_infinite(hint) || 
            hint->info().urban_tag == Urban_object_type::BUILDING || 
            hint->info().urban_tag == Urban_object_type::TREE)
          continue;

        const auto filter = face_map.find(hint);
        CGAL_assertion(filter != face_map.end());
        Candidate_face_ptr candidate = filter->second;
        
        const Triangle_3 triangle = internal::triangle_3<Triangle_3>(hint);
        const FT sq_dist = CGAL::squared_distance(point_3, triangle);

        candidate->inliers.push_back(pit);
        candidate->max_error = CGAL::max(candidate->max_error, sq_dist);
      }

      for (auto fit = face_map.begin(); fit != face_map.end(); ++fit)
        todo.insert(fit->second);
    }

    bool well_shaped(const Face_handle& fh) const {
      
      const Point_2& pa = fh->vertex(0)->point();
      const Point_2& pb = fh->vertex(1)->point();
      const Point_2& pc = fh->vertex(2)->point();

      double area = 2.0 * CGAL::to_double(CGAL::area(pa, pb, pc));
      area = area * area;

      const double a = CGAL::to_double(CGAL::squared_distance(pb, pc));
      const double b = CGAL::to_double(CGAL::squared_distance(pc, pa));
      const double c = CGAL::to_double(CGAL::squared_distance(pa, pb));

      double q;
      if(a < b) {
        if(a < c)
          q = area / (b * c);
        else
          q = area / (a * b);
      } else {
        if(b < c)
          q = area / (a * c);
        else
          q = area / (a * b);
      }
      return q > 0.125;
    }

    bool is_too_small(
      const Face_handle& fh, 
      const FT min_size) const {
      
      for (std::size_t i = 0; i < 3; ++i) {
        
        const Point_2& a = fh->vertex(i)->point();
        const Point_2& b = fh->vertex((i + 1) % 3)->point();
        
        if (CGAL::squared_distance(a, b) > min_size * min_size)
          return false;
      }
      return true;
    }

  }; // Smooth_ground_estimator

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_SMOOTH_GROUND_ESTIMATOR_H
