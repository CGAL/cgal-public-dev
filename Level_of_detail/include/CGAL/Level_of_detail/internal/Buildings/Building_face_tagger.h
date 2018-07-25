#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_FACE_TAGGER_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_FACE_TAGGER_H

// STL includes.
#include <utility>

// CGAL includes.
#include <CGAL/Barycentric_coordinates_2/Segment_coordinates_2.h>

// LOD includes.
#include <CGAL/Level_of_detail/Enumerations.h>

namespace CGAL {

namespace Level_of_detail {

namespace BC  = CGAL::Barycentric_coordinates;

template<class GeomTraits, class Triangulation>
class Building_face_tagger {
			
public:
  using Kernel        = GeomTraits;

  typename Kernel::Compute_squared_distance_2 squared_distance_2;

  using FT        = typename Kernel::FT;
  using Line_2    = typename Kernel::Line_2;
  using Point_2   = typename Kernel::Point_2;
  using Segment_2 = typename Kernel::Segment_2;
            
  using Triangulation_face_handle    = typename Triangulation::Face_handle; 
  using Triangulation_faces_iterator = typename Triangulation::Finite_faces_iterator;
  using All_faces_iterator = typename Triangulation::All_faces_iterator;
  using Triangulation_vertex_handle  = typename Triangulation::Vertex_handle;

  using Visibility_label     = LOD::Visibility_label;

  using Segment_coordinates = std::pair<FT, FT>;

  Building_face_tagger(
    const Triangulation &triangulation, 
    const FT constraints_threshold) :
    m_triangulation(triangulation),
    m_constraints_threshold(constraints_threshold),
    m_tolerance(-FT(1) / FT(10)),
    m_building_index(-1),
    m_start_new_building(true) { 

    set_default_building_indices();
  }

  void tag_according_to_constraints(const std::vector<Segment_2>& segments) {

    for (Triangulation_faces_iterator tf_it = m_triangulation.finite_faces_begin();
         tf_it != m_triangulation.finite_faces_end(); ++tf_it) {
      flood(segments, tf_it);
      m_start_new_building = true;
    }
  }

private:
  const Triangulation &m_triangulation;
            
  const FT m_constraints_threshold;
  const FT m_tolerance;

  int  m_building_index;
  bool m_start_new_building;

  void set_default_building_indices() {
    for (All_faces_iterator tf_it = m_triangulation.all_faces_begin();
         tf_it != m_triangulation.all_faces_end(); ++tf_it)
      tf_it->info().group_number() = -1;
  }

  void flood(const std::vector<Segment_2>& segments, Triangulation_face_handle face_handle) {
				
    // If this face is not valid due to some criteria, we do not handle this face and skip it.
    if (is_not_valid_face(face_handle)) return;

    // Check if we should increase the buildings counter.
    set_building_index_to_face(face_handle);

    // Continue propagating.
    for (size_t i = 0; i < 3; ++i) {
      Triangulation_face_handle face_handle_neighbour = face_handle->neighbor(i);

      if (!is_constrained_edge(face_handle, i, segments)) 
        flood(segments, face_handle_neighbour);
    }
  }

  bool is_not_valid_face(Triangulation_face_handle face_handle) const {

    // This face was already handled before.
    const int building_number = face_handle->info().group_number();
    if (building_number >= 0) return true;

    // This is infinite face.
    if (m_triangulation.is_infinite(face_handle)) return true;

    // This face is outside of any building.
    const Visibility_label visibility_label = face_handle->info().visibility_label();
    if (visibility_label == Visibility_label::OUTSIDE) return true;

    // The face is valid.
    return false;
  }

  void set_building_index_to_face(Triangulation_face_handle face_handle) {
                
    if (m_start_new_building) {
      m_start_new_building = false; ++m_building_index;
    }
    face_handle->info().group_number() = m_building_index;
  }

  bool is_constrained_edge(Triangulation_face_handle face_handle,
                           const size_t vertex_index,
                           const std::vector<Segment_2>& segments) const {

    const Point_2 &p1 = face_handle->vertex((vertex_index + 1) % 3)->point();
    const Point_2 &p2 = face_handle->vertex((vertex_index + 2) % 3)->point();

    return check_constraint(p1, p2, segments);
  }

  bool check_constraint(const Point_2 &p1, const Point_2 &p2,
                        const std::vector<Segment_2>& segments) const {
				
    CGAL_precondition(segments.size() > 0);

    for (std::size_t i = 0; i < segments.size(); ++ i)
      if (is_constrained(p1, p2, segments[i]))
        return true;
    return false;
  }

  bool is_constrained(const Point_2 &p1, const Point_2 &p2, const Segment_2 &segment) const {

    const Point_2 &source = segment.source();
    const Point_2 &target = segment.target();

    Line_2 line(source, target);

    const Point_2 pr1 = line.projection(p1);
    const Point_2 pr2 = line.projection(p2);

    const FT squared_threshold = m_constraints_threshold * m_constraints_threshold;

    if (squared_distance_2(p1, pr1) > squared_threshold) return false;
    if (squared_distance_2(p2, pr2) > squared_threshold) return false;
				
    Segment_coordinates bc = CGAL::make_pair(BC::compute_segment_coordinates_2(source, target, p1, Kernel()));
    const bool state1 = bc.first > m_tolerance && bc.second > m_tolerance;

    bc = CGAL::make_pair(BC::compute_segment_coordinates_2(source, target, p2, Kernel()));
    const bool state2 = bc.first > m_tolerance && bc.second > m_tolerance;

    bc = CGAL::make_pair(BC::compute_segment_coordinates_2(p1, p2, source, Kernel()));
    const bool state3 = bc.first > m_tolerance && bc.second > m_tolerance;

    bc = CGAL::make_pair(BC::compute_segment_coordinates_2(p1, p2, target, Kernel()));
    const bool state4 = bc.first > m_tolerance && bc.second > m_tolerance;

    if ( (state1 && state2) || (state3 && state4) ) return true;
    return false;
  }
};

} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_FACE_TAGGER_H
