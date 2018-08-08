#ifndef CGAL_LEVEL_OF_DETAIL_TRIANGULATION_GROUND_REFINER_H
#define CGAL_LEVEL_OF_DETAIL_TRIANGULATION_GROUND_REFINER_H

// STL includes.
#include <vector>

// LOD includes.
#include <CGAL/Level_of_detail/internal/Data/Kd_tree_with_data_creator.h>
#include <CGAL/Level_of_detail/internal/Property_maps/Point_2_from_iterator_map.h>
#include <CGAL/Level_of_detail/internal/utils.h>
#include <CGAL/number_utils.h>

namespace CGAL {

namespace Level_of_detail {

template<class GeomTraits, class Triangulation,
         typename PointRange, typename PointMap,
         typename Tree>
class Triangulation_ground_refiner
{

public:
  using Kernel        = GeomTraits;

  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  using Triangle_3 = typename Kernel::Triangle_3;
  using FT = typename Kernel::FT;

  using Point_iterator = typename PointRange::const_iterator;
  using Face_handle = typename Triangulation::Face_handle;
  using Vertex_handle = typename Triangulation::Vertex_handle;

private:
  
  Triangulation& m_triangulation;
  PointRange m_points;
  PointMap m_point_map;
  const Tree& m_tree;

  struct Candidate_face
  {
    Face_handle face;
    FT max_error;
    std::vector<Point_iterator> inliers;

    Candidate_face (Face_handle face = Face_handle())
      : face (face), max_error (FT(0))
    { }
    
  };

  using Candidate_face_ptr = std::shared_ptr<Candidate_face>;

  struct Compare_candidates
  {
    bool operator() (const Candidate_face_ptr& a, const Candidate_face_ptr& b) const
    {
      if (a->max_error != b->max_error)
        return a->max_error > b->max_error;
      // else
      return a->face < b->face;
    }
  };

  using Face_map = std::map<Face_handle, Candidate_face_ptr>;
  using Face_queue = std::set<Candidate_face_ptr, Compare_candidates>;

public:

  Triangulation_ground_refiner (Triangulation& triangulation,
                                const PointRange& points,
                                PointMap point_map,
                                const Tree& tree)
    : m_triangulation (triangulation)
    , m_points (points)
    , m_point_map (point_map)
    , m_tree (tree)
  {
  }

  void output_face (Candidate_face_ptr candidate)
  {
    static int nb = 0;

    if (candidate->inliers.empty())
      return;
    
    ++ nb;

    std::string fname = "face_";
    if (nb < 10)
      fname = fname + "0" + std::to_string(nb) + ".xyz";
    else
      fname = fname + std::to_string(nb) + ".xyz";
    
    std::ofstream f (fname);
    f.precision(18);

    for (std::size_t i = 0; i < candidate->inliers.size(); ++ i)
      f << get(m_point_map, *(candidate->inliers[i])) << std::endl;
  }

  void refine (FT tolerance)
  {
    Face_map face_map;
    Face_queue todo;
    init_queue (face_map, todo, tolerance);

    Face_handle hint;

    std::ofstream file ("circumcenters.xyz");
    std::ofstream file2 ("centers.xyz");

    std::ofstream file3 ("not_found.polylines.txt");
    file.precision(18);
    file2.precision(18);
    file3.precision(18);

    
    while (!todo.empty())
    {
      Candidate_face_ptr candidate = *(todo.begin());
      todo.erase (todo.begin());
//      face_map.erase (candidate->face);

      bool out_of_tolerance =  (candidate->max_error > tolerance * tolerance);
      bool badly_shaped = (!well_shaped (candidate->face));

      if (!out_of_tolerance && !badly_shaped)
        continue;

      if (too_small (candidate->face, 3. * tolerance))
        continue;
      
      // Get circumcenter and conflict zone
      Point_2 center = CGAL::circumcenter (candidate->face->vertex(0)->point(),
                                           candidate->face->vertex(1)->point(),
                                           candidate->face->vertex(2)->point());
      if (out_of_tolerance)
      {
        center = CGAL::barycenter (candidate->face->vertex(0)->point(), 1.,
                                   candidate->face->vertex(1)->point(), 1.,
                                   candidate->face->vertex(2)->point(), 1.);
        file2 << center << " 0" << std::endl;
      }
      else
        file << center << " 0" << std::endl;

      typename Triangulation::Locate_type lt;
      int li;
      hint = m_triangulation.locate (center, lt, li, hint);
      
      if (lt == Triangulation::VERTEX ||
          m_triangulation.is_infinite (hint) ||
          hint->info().visibility_label() != Visibility_label::OUTSIDE)
      {
        continue;
      }


      std::vector<Face_handle> conflict;
      m_triangulation.get_conflicts (center, std::back_inserter (conflict));

      // Recover points and remove conflict cells from local structures
      std::vector<Point_iterator> points;
      
      for (std::size_t i = 0; i < conflict.size(); ++ i)
      {
        if (m_triangulation.is_infinite(conflict[i])
            || conflict[i]->info().visibility_label() != Visibility_label::OUTSIDE)
          continue;

        typename Face_map::iterator fiter = face_map.find(conflict[i]);
        if (fiter == face_map.end())
          continue;

        Candidate_face_ptr cface = fiter->second;
        
        if (!(cface->inliers.empty()))
          std::copy (cface->inliers.begin(), cface->inliers.end(),
                     std::back_inserter (points));

        face_map.erase (fiter);
        todo.erase(cface);
      }

      // Insert new vertex
      Vertex_handle v = m_triangulation.insert (center, hint);

      std::vector<Candidate_face_ptr> new_faces;
      
      typename Triangulation::Face_circulator circ = m_triangulation.incident_faces(v);
      typename Triangulation::Face_circulator start = circ;
      do
      {
        circ->info().visibility_label() = Visibility_label::OUTSIDE;
        compute_heights (circ);
        Candidate_face_ptr cface = std::make_shared<Candidate_face>(circ);
        face_map.insert (std::make_pair (circ, cface));
        new_faces.push_back (cface);
        ++ circ;
      }
      while (circ != start);

      // Redistribute points
      for (std::size_t i = 0; i < points.size(); ++ i)
      {
        const Point_3& point_3 = get(m_point_map, *(points[i]));
        Point_2 point_2 = internal::point_2_from_point_3(point_3);
      
        hint = m_triangulation.locate (point_2, hint);

        typename Face_map::iterator fiter = face_map.find(hint);
        CGAL_assertion (fiter != face_map.end());
        Candidate_face_ptr cface = fiter->second;
        
        Triangle_3 triangle = internal::triangle_3<Triangle_3>(hint);
      
        FT sq_dist = CGAL::squared_distance (point_3, triangle);

        // if (sq_dist > tolerance * tolerance)
        {
          cface->inliers.push_back (points[i]);
          cface->max_error = (std::max)(cface->max_error, sq_dist);
        }
      }

      // insert new faces
      for (std::size_t i = 0; i < new_faces.size(); ++ i)
        todo.insert (new_faces[i]);
    }
    
  }

private:

  void init_queue (Face_map& face_map,
                   Face_queue& todo,
                   FT tolerance)
  {

    for (typename Triangulation::Finite_faces_iterator
           it = m_triangulation.finite_faces_begin();
         it != m_triangulation.finite_faces_end(); ++ it)
    {
      if (m_triangulation.is_infinite(it)
          || it->info().visibility_label() != Visibility_label::OUTSIDE)
        continue;
      
      face_map.insert (std::make_pair (it, std::make_shared<Candidate_face>(it)));
    }

    Face_handle hint;
    for (Point_iterator it = m_points.begin(); it != m_points.end(); ++ it)
    {
      const Point_3& point_3 = get(m_point_map, *it);
      Point_2 point_2 = internal::point_2_from_point_3(point_3);
      
      hint = m_triangulation.locate (point_2, hint);

      if (m_triangulation.is_infinite(hint)
          || hint->info().visibility_label() != Visibility_label::OUTSIDE)
        continue;

      typename Face_map::iterator fiter = face_map.find(hint);
      CGAL_assertion (fiter != face_map.end());
      Candidate_face_ptr candidate = fiter->second;
      
      Triangle_3 triangle = internal::triangle_3<Triangle_3>(hint);
      
      FT sq_dist = CGAL::squared_distance (point_3, triangle);

//      if (sq_dist > tolerance * tolerance)
      {
        candidate->inliers.push_back (it);
        candidate->max_error = (std::max)(candidate->max_error, sq_dist);
      }
    }

    for (typename Face_map::iterator it = face_map.begin();
         it != face_map.end(); ++ it)
      todo.insert (it->second);
  }

  bool well_shaped (Face_handle face) const
  {
    const Point_2& pa = face->vertex(0)->point();
    const Point_2& pb = face->vertex(1)->point();
    const Point_2& pc = face->vertex(2)->point();

    double area = 2 * CGAL::to_double (CGAL::area(pa, pb, pc));
    area = area * area;

    double a = CGAL::to_double (CGAL::squared_distance(pb, pc));
    double b = CGAL::to_double (CGAL::squared_distance(pc, pa));
    double c = CGAL::to_double (CGAL::squared_distance(pa, pb));

    double q;
    if(a<b)
    {
      if(a<c)
        q = area/(b*c);
      else
        q = area/(a*b);
    }
    else
    {
      if(b<c)
        q = area/(a*c);
      else
        q = area/(a*b);
    }
    return q > 0.125;
  }

  void compute_heights (Face_handle face) const
  {
    for (std::size_t j = 0; j < 3; ++ j)
    {
      const Point_2& p2 = face->vertex(j)->point();
      typename Tree::Neighbors neighbors;
      m_tree.search_knn_2 (p2, neighbors);
      double h_mean = 0.;
      for (std::size_t i = 0; i < neighbors.size(); ++ i)
        h_mean += get (m_point_map, *(neighbors[i])).z();
      face->info().height(j) = h_mean / neighbors.size();
    }
  }

  bool too_small (Face_handle face, FT minimum_size) const
  {
    for (std::size_t i = 0; i < 3; ++ i)
    {
      const Point_2 a = face->vertex(i)->point();
      const Point_2 b = face->vertex((i+1)%3)->point();
      if (CGAL::squared_distance(a,b) > minimum_size * minimum_size)
        return false;
    }
    return true;
  }

};

} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_TRIANGULATION_GROUND_REFINER_H
