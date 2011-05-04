// Offset_construction_2.h

#ifndef Offset_construction_2_h
#define Offset_construction_2_h

#include <CGAL/Offset_statistics_2.h>
#include <CGAL/approximated_offset_2.h>
#include <CGAL/offset_polygon_2.h>

#include <CGAL/create_straight_skeleton_2.h>
//#include <CGAL/create_offset_polygons_2.h>
#include <boost/shared_ptr.hpp>

namespace CGAL {


/*! \classf
 * A class implementing construction of polygon offsets 
 * via existing CGAL methods with govering statistics. 
 * The methods include exact offset construction and approximate
 * (Ron Wein) offset construction. Will add constructions via 
 * skeletons in the future.
 * The inputs are defined by caller application,
 * the class keeps them as _references_. The outputs are defined 
 * by the class and caller application may keep a reference to them.
 */

template <typename CPolygon_2>
class Offset_construction_2 : public Offset_statistics_2<CPolygon_2>
{
public:

  typedef Offset_statistics_2<CPolygon_2> Base;
  typedef typename Base::Types Types;
  typedef typename Types::Polygon_2 Polygon_2;
  typedef typename Types::Polygon_traits_2 Polygon_traits_2;  
  typedef typename Types::Input_rational Input_rational;
  
  typedef typename Base::Polygon_statistics Polygon_statistics;
  typedef typename Base::Offset_statistics Offset_statistics;

  typedef typename Types::Approximate_offset_polygon_2 Approximate_offset_polygon_2; 
  typedef typename Types::Approximate_polygon_2 Approximate_polygon_2;
  
  typedef typename Types::Conic_traits_2 Conic_traits_2;
  typedef typename Types::Rat_polygon_2 Rat_polygon_2;
  typedef typename Types::Exact_rational Exact_rational;
  typedef typename Types::Exact_offset_polygon_2 Exact_offset_polygon_2; 
  typedef typename Types::Exact_polygon_2 Exact_polygon_2;

  typedef typename CGAL::Straight_skeleton_2<Polygon_traits_2> Straight_skeleton_2;
  typedef typename boost::shared_ptr<Straight_skeleton_2> Straight_skeleton_ptr_2;

  // typedef typename boost::shared_ptr<Polygon_2>  Polygon_ptr_2;
  // typedef typename std::vector<Polygon_ptr_2>    Polygon_sequence_2 ;

  struct Statistics
  {
      // exact offset statistics
      Polygon_statistics exact_input;
      Offset_statistics exact_stats;

      // approximate offset statistics
      Offset_statistics approximate_stats;

      // skeleton statistics
      Polygon_statistics inner_skeleton;
      Polygon_statistics outer_skeleton;
  };

  /// input parameters stored _as references_
  Offset_construction_2(
    const Polygon_2& i_polygon,
    const Input_rational& i_offset,
    const Input_rational& i_eps) :
    m_polygon(i_polygon),
    m_offset(i_offset),
    m_epsilon(i_eps)
  {
  }

  ~Offset_construction_2(void)
  {
  }

  const Statistics& statistics() const { return m_statistics; }

  const Exact_offset_polygon_2& exact_offset_polygon() const 
  { 
    return m_exact_offset_polygon; 
  }
  const Exact_polygon_2& exact_outer_boundary() const 
  { 
    return m_exact_outer_boundary; 
  }

  const Approximate_offset_polygon_2& approximate_offset_polygon() const 
  { 
    return m_approximate_offset_polygon; 
  }
  const Approximate_polygon_2& approximate_outer_boundary() const 
  { 
    return m_approximate_outer_boundary; 
  }

  const Straight_skeleton_ptr_2& inner_skeleton_ptr() const 
  { 
    return mp_inner_skeleton; 
  }

  const Straight_skeleton_ptr_2& outer_skeleton_ptr() const 
  { 
    return mp_outer_skeleton; 
  }

  void clear();

  // exact offset construction
  void compute_exact_offset_polygon();

  // approximate offset construction
  void compute_approximate_offset_polygon();

  // inner skeleton construction
  void compute_inner_skeleton();

  // outer skeleton construction
  void compute_outer_skeleton();

private:

  // input parameters
  const Polygon_2& m_polygon;
  const Input_rational& m_offset;
  const Input_rational& m_epsilon;

  Rat_polygon_2 m_rat_polygon;
  Exact_offset_polygon_2 m_exact_offset_polygon;
  Exact_polygon_2 m_exact_outer_boundary;

  Approximate_offset_polygon_2 m_approximate_offset_polygon;
  Approximate_polygon_2 m_approximate_outer_boundary;

  Straight_skeleton_ptr_2 mp_inner_skeleton;
  Straight_skeleton_ptr_2 mp_outer_skeleton;
  
  // segment dalaunay graph / voronoi diagram 
  /*
  Segment_Delaunay_graph_traits_without_intersections_2<K,MTag> Segment_Delaunay_graph_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>
  */

  Statistics m_statistics;

};


template <typename CPolygon_2>
void Offset_construction_2<CPolygon_2>::compute_inner_skeleton()
{
  mp_inner_skeleton.reset();
  
  CGAL::Timer timer;
  timer.start();
  if(!m_polygon.is_empty() )
  {
    assert(m_polygon.is_simple());
    assert(m_polygon.is_counterclockwise_oriented());

    Polygon_traits_2 traits;
    //Straight_skeleton_ptr_2 is_ptr = 
    mp_inner_skeleton = 
      CGAL::create_interior_straight_skeleton_2(
      m_polygon.vertices_begin(), 
      m_polygon.vertices_end(),
      traits);
  }
  timer.stop();

  m_statistics.inner_skeleton.time = timer.time();
  m_statistics.inner_skeleton.size = 0;
  
  if(mp_inner_skeleton)
  {
      const Straight_skeleton_2 & inner_skeleton = *mp_inner_skeleton;
      m_statistics.inner_skeleton.size = inner_skeleton.size_of_halfedges();
      
      /* 
      // Get not rounded offset contours

       Polygon_traits_2 traits;
      // Get offset contours
      Polygon_sequence_2 offset_contours = 
        CGAL::create_offset_polygons_2(
        m_offset - m_epsilon, 
        *mp_inner_skeleton, 
        traits);
    
      // Dump the generated offset polygons
      OUT_DEBUG << offset_contours.size() << " offset contours obtained" << std::endl;
    
      for (typename Polygon_sequence_2::const_iterator i = offset_contours.begin(); 
        i != offset_contours.end(); ++ i )
      {
        // Each element in the offset_contours sequence is a shared pointer 
        // to a Polygon_2 instance.
        OUT_DEBUG << (*i)->size() << " vertices in offset contour\n" ;
    
        for (typename Polygon_2::Vertex_const_iterator j = (*i)->vertices_begin(); 
          j != (*i)->vertices_end(); ++ j )
           OUT_DEBUG << "(" << j->x() << "," << j->y() << ")" << std::endl ;
      }
      */
  
  }

  //if (verbose())
  {
    OUT_DEBUG << "Inner skeleton of size " << m_statistics.inner_skeleton.size
    << " computed in " << m_statistics.inner_skeleton.time << std::endl;
  }

}

template <typename CPolygon_2>
void Offset_construction_2<CPolygon_2>::compute_outer_skeleton()
{
  mp_outer_skeleton.reset();

  CGAL::Timer timer;
  timer.start();
  if(!m_polygon.is_empty() )
  {
    assert(m_polygon.is_simple());
    assert(m_polygon.is_counterclockwise_oriented());

    Input_rational max_offset = 0;
    // compute max offset
    {
      Input_rational x_diff = m_polygon.left_vertex()->x() - m_polygon.right_vertex()->x();
      Input_rational y_diff = m_polygon.top_vertex()->x() - m_polygon.bottom_vertex()->x();
      max_offset = (x_diff > y_diff? x_diff: y_diff);
    }

    LOG_DEBUG << "max_offset = " << max_offset << ", m_offset = " << m_offset << std::endl;

    Polygon_traits_2 traits;
    //Straight_skeleton_ptr_2 is_ptr = 
    mp_outer_skeleton = 
      CGAL::create_exterior_straight_skeleton_2(
      max_offset + m_offset,
      m_polygon, 
      traits);
  }
  timer.stop();

  m_statistics.outer_skeleton.time = timer.time();
  m_statistics.outer_skeleton.size = 0;
  
  if(mp_outer_skeleton)
  {
      const Straight_skeleton_2 & outer_skeleton = *mp_outer_skeleton;
      m_statistics.outer_skeleton.size = outer_skeleton.size_of_halfedges();
  }

  //if (verbose())
  {
    OUT_DEBUG << "Inner skeleton of size " << m_statistics.outer_skeleton.size
      << " computed in " << m_statistics.outer_skeleton.time << std::endl;
  }
}

template <typename CPolygon_2>
void Offset_construction_2<CPolygon_2>::clear()
{
  m_exact_offset_polygon.clear();
  m_approximate_offset_polygon.clear();
}

template <typename CPolygon_2>
void Offset_construction_2<CPolygon_2>::compute_exact_offset_polygon()
{
  // Compute the offset polygon.
  Exact_rational radius = type_cast<Exact_rational>(m_offset);
//  LOG_DEBUG << "exact radius " << radius << std::endl;
  {
      CGAL::Timer timer;
      timer.start();

      copy(m_polygon, m_rat_polygon);
      timer.stop();

      m_statistics.exact_input.time = timer.time();
      m_statistics.exact_input.size = m_rat_polygon.size();

      //if (verbose())
      {
          OUT_DEBUG << "Exact polygon of size " << m_statistics.exact_input.size
            << " computed in " << m_statistics.exact_input.time << std::endl;
      }
  }
//  LOG_DEBUG << "exact source polygon size " << m_rat_polygon.size() << std::endl;

  {
      CGAL::Timer timer;
      if(! m_rat_polygon.is_empty() ) {
        timer.start();
        Conic_traits_2 conic_traits;
        m_exact_offset_polygon = offset_polygon_2 (m_rat_polygon, radius, conic_traits);
        timer.stop();
        
        m_exact_outer_boundary = m_exact_offset_polygon.outer_boundary();
        m_statistics.exact_stats.time = timer.time();
        m_statistics.exact_stats.size = m_exact_outer_boundary.size();
        m_statistics.exact_stats.holes = m_exact_offset_polygon.number_of_holes();
      }
  }
  
  // TODO write operator<<
  //if (verbose())
  {
    OUT_DEBUG << "The exact offset polygon has "
              << m_statistics.exact_stats.size << " vertices, "
              << m_statistics.exact_stats.holes << " holes." << std::endl;
    OUT_DEBUG << "Exact offset computation took "
              << m_statistics.exact_stats.time << " seconds." << std::endl << std::endl;
  }
}

template <typename CPolygon_2>
void Offset_construction_2<CPolygon_2>::compute_approximate_offset_polygon()
{
// Approximate the offset polygon.
  CGAL::Timer                  timer;

  timer.start();
  m_approximate_offset_polygon = approximated_offset_2 (m_polygon, m_offset - m_epsilon, to_double(m_epsilon));
  timer.stop();

  m_approximate_outer_boundary = m_approximate_offset_polygon.outer_boundary();

  m_statistics.approximate_stats.time = timer.time();
  m_statistics.approximate_stats.size = m_approximate_outer_boundary.size();
  m_statistics.approximate_stats.holes = m_approximate_offset_polygon.number_of_holes();

  // TODO write operator<<
  //if (verbose())
  {
    OUT_DEBUG << "The approximate offset polygon has "
              << m_statistics.approximate_stats.size << " vertices, "
              << m_statistics.approximate_stats.holes << " holes." << std::endl;
    OUT_DEBUG << "Approximate offset computation took "
              << m_statistics.approximate_stats.time << " seconds." << std::endl << std::endl;
  }

}

} // namespace CGAL


#endif // Offset_construction_2_h
