// Minkowski_sum_construction_2.h

#ifndef Minkowski_sum_construction_2_h
#define Minkowski_sum_construction_2_h

#include <CGAL/Offset_statistics_2.h>
#include <CGAL/Circle_approximation_2.h>
#include <CGAL/Self_intersection_handler_2.h>
#include <CGAL/General_polygon_set_2.h>

#include <CGAL/minkowski_sum_2.h>

namespace CGAL {


/*! \classf
 * A class implementing construction of approximate
 * polygon offsets via minkowski sum with kgon (that 
 * represents approximate circle). 
 * Input parameters stored _as references_, computed
 * data recalculated on demand (via compute_*() functions).
 */

template <typename CPolygon_2>
class Minkowski_sum_construction_2 : public Offset_statistics_2<CPolygon_2>
{
public:

  typedef Offset_statistics_2<CPolygon_2> Base;
  typedef typename Base::Types Types;

  typedef typename Types::Polygon_2 Polygon_2;
  typedef typename Types::Polygon_traits_2 Polygon_traits_2;
  typedef typename Types::Input_rational Input_rational;
  typedef typename Base::Offset_statistics Offset_statistics;
  typedef typename Types::Rat_traits_2 Rat_traits_2;
  typedef typename Types::Approximate_offset_polygon_2 Approximate_offset_polygon_2; 
  typedef typename Types::Approximate_polygon_2 Approximate_polygon_2;
  typedef typename Types::Approximate_polygon_list_2 Approximate_polygon_list_2;
  
  struct Statistics
  {
      Offset_statistics kgon_sum_stats;
      Offset_statistics kgon_offset_stats;

      Offset_statistics kgon_diff_stats;
      
      double kgon_circles_time;
  };

  typedef typename CGAL::Circle_approximation_2<CPolygon_2>    Circle_app_2;
  typedef typename Circle_app_2::Directed_slope_set Directed_slope_set;
  typedef typename Circle_app_2::Equal_directed_slope_2 Equal_directed_slope_2;

  typedef typename Types::Polygon_with_holes_2  Kgon_sum_polygon_2;
  
  typedef typename CGAL::General_polygon_set_2<Rat_traits_2> Polygon_with_hole_set_2;
  typedef typename std::list<Approximate_offset_polygon_2> Approximate_offset_polygon_list_2;
  typedef typename std::list<Kgon_sum_polygon_2> Kgon_sum_polygon_list_2;
  
  typedef typename CGAL::Arrangement_2<Rat_traits_2>             Arrangement_2;
  typedef typename CGAL::Self_intersection_handler_2<CPolygon_2> Intersection_handler_2;

  typedef typename Rat_traits_2::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Rat_traits_2::Point_2                      Arc_point_2;

  typedef typename Polygon_traits_2::Point_2  Point_2;
  typedef typename Polygon_traits_2::Circle_2 Circle_2;
  typedef typename Polygon_traits_2::Segment_2  Segment_2;
  typedef typename Polygon_traits_2::Direction_2  Direction_2;
  typedef typename Polygon_traits_2::Vector_2  Vector_2;

  typedef typename std::list<Circle_2> Circle_list_2;

  typedef Approximate_polygon_2  Contour_2;
  typedef typename Contour_2::Curve_iterator Contour_2_iterator;
  struct Arc_data_2
  {
    
    Arc_data_2(
      const Point_2& i_source,
      const Point_2& i_target,
      const Circle_2& i_circle,
      Contour_2& i_contour,
      Contour_2_iterator i_cit
      ) :
      circle(i_circle),
      source(i_source),
      target(i_target),
      contour(i_contour),
      cit(i_cit)
      {
//        CGAL::Comparison_result start_center = CGAL::compare(start_point.y(), arc_center.y());
//        CGAL::Comparison_result end_center = CGAL::compare(end_point.y(), arc_center.y());
        
//        is_upper = (start_center == LARGER) || (end_center == LARGER);
      //    is_upper = true;
      }
//    bool is_upper_half() { return is_upper; }
    Point_2 source;
    Point_2 target;
    Circle_2 circle;
    Contour_2& contour;
    Contour_2_iterator cit;
    //bool is_upper;
  };

  struct Arc_detection_data_2
  {
    Arc_detection_data_2(typename Polygon_2::Edge_const_circulator ecirc) :
      first_ecirc(ecirc),
      last_ecirc(ecirc),
      known_arc_position(false)
      {}
    typename Polygon_2::Edge_const_circulator first_ecirc;
    typename Polygon_2::Edge_const_circulator last_ecirc;
    Segment_2 first_segment;
    Segment_2 last_segment;
    Point_2 arc_center;
    bool known_arc_position;
    CGAL::Orientation orientation;
  };
  
  typedef typename std::multimap<Point_2, Arc_data_2> Arc_center_multimap_2;
  typedef typename Arc_center_multimap_2::value_type Arc_center_val_2;

  /// input parameters stored _as references_
  Minkowski_sum_construction_2(
    const Polygon_2& i_polygon,
    const Input_rational& i_offset,
    const Input_rational& i_eps,
    Circle_app_2& i_circle_app) :
    m_polygon(i_polygon),
    m_offset(i_offset),
    m_epsilon(i_eps),
    m_circle_app(i_circle_app),
    m_kgon(m_circle_app.kgon()) // also stored as reference
  {
    //m_circle_app.verbose(false);
  }

  ~Minkowski_sum_construction_2(void)
  {
  }

  virtual void verbose(const bool& be_verbose) { Base::verbose(be_verbose); m_circle_app.verbose(be_verbose); }

  const Statistics& statistics() const { return m_statistics; }

  const Kgon_sum_polygon_2& kgon_sum() const { return m_kgon_sum; }
  const Polygon_2& kgon_sum_outer_boundary() const { return m_kgon_sum_outer_boundary; }
  const Polygon_with_hole_set_2& kgon_diff() const { return m_kgon_diff; }
  const Polygon_2& kgon_diff_first_boundary() const { return m_kgon_diff_first_boundary; }
  const Circle_list_2& kgon_circles() const { return m_kgon_circles; }
  const Circle_list_2& kgon_skipped_circles() const { return m_kgon_skipped_circles; }
  const Circle_list_2& kgon_intersecting_circles() const { return m_kgon_intersecting_circles; }
  const Approximate_offset_polygon_2& kgon_offset_polygon() const { return m_kgon_offset_polygon; }
  const Approximate_polygon_2& kgon_offset_outer_boundary() const { return m_kgon_offset_outer_boundary; }

  void clear();

  /// compute minkowski sum with kgon
  void compute_kgon_sum();

  /// compute complement of minkowski sum of complement of the input with kgon, that is minkowski "difference"
  void compute_kgon_diff();
  
  /// detect arcs by kgon preimage sequences
  void compute_kgon_induced_arcs();
  /// detect self-intersections caused by arcs and modify result to avoid them
  void compute_kgon_offset_polygon();

private:

  /// detect arcs by kgon preimage sequences in single contour
  void compute_polygon_arcs(const Polygon_2& polygon, const Directed_slope_set& kgon_slopes);

  // TODO: move arc detection to subclass
  void complete_and_add_circle(const Segment_2& first_segment, const Segment_2& last_segment,
    const Point_2& arc_center, const bool& known_arc_position, 
    const CGAL::Orientation& orientation, const Directed_slope_set& kgon_slopes,
    Approximate_polygon_2& contour, int& total_not_known);
  void add_x_monotone_arcs(const Point_2& start_point, const Point_2& end_point,
    const Circle_2& circle, Approximate_polygon_2& contour);

  void overwrite(const Polygon_2& p1, Approximate_polygon_2& p2)
  {
    p2.clear();
    if(p1.is_empty()) return;

    typename Polygon_2::Edge_const_circulator ecirc = p1.edges_circulator();
    typename Polygon_2::Edge_const_circulator ecirc_end = ecirc;
    do
    {
      p2.push_back(X_monotone_curve_2 ((*ecirc).source(), (*ecirc).target()));
      ++ecirc;
    }while(ecirc != ecirc_end);
  }

//   void overwrite(const Approximate_polygon_2& p1, Polygon_2& p2)
//   {
//     // pints of p1 are rational
//     p2.clear();
//     if(p1.is_empty()) return;
// 
//     typename Approximate_polygon_2::Curve_const_iterator cit = p1.curves_begin();
//     typename Approximate_polygon_2::Curve_const_iterator cit_end = p1.curves_end();
//     do
//     {
//       p2.push_back(Point_2((*cit).source().x().alpha(), (*cit).source().y().alpha()));
//       ++cit;
//     }while(cit != cit_end);
//   }

  static bool is_simple(const Approximate_polygon_2& polygon) { 
    Arrangement_2 arr;
    insert (arr, polygon.curves_begin(), polygon.curves_end());
    LOG_DEBUG << "no of faces " << arr.number_of_faces() << std::endl;
    return (arr.number_of_faces() <= 2); // supposed to be one closed contour
  }

  // compute bounding box of at least radius away from polygon
  // place a polygon as a hole in it
  // cut the polygon with a single hole to two simple polygons (left and right)
  // by cutting through hole via upper and lower vertex
  static void mold_polygon(const Polygon_2& i_polygon,
                          const Input_rational& radius,
                          Polygon_2& o_polygon_left,
                          Polygon_2& o_polygon_right);


  // input parameters
  const Polygon_2& m_polygon;
  const Input_rational& m_offset;
  const Input_rational& m_epsilon;
  Circle_app_2& m_circle_app;
  const Polygon_2& m_kgon;

  // offset results
  Kgon_sum_polygon_2 m_kgon_sum;
  Polygon_2 m_kgon_sum_outer_boundary;
  Polygon_with_hole_set_2 m_kgon_diff;
  Polygon_2 m_kgon_diff_first_boundary;
  
  Circle_list_2 m_kgon_circles;
  Circle_list_2 m_kgon_skipped_circles;
  Circle_list_2 m_kgon_intersecting_circles;
  Approximate_offset_polygon_2 m_kgon_offset_polygon;
  Approximate_polygon_2 m_kgon_offset_outer_boundary;
  Approximate_polygon_list_2 m_kgon_sum_with_arcs;

  Arc_center_multimap_2 m_arc_center_multimap;
  // Kgon_sum_polygon_2 m_kgon_sum_arcs;
  // Kgon_sum_polygon_2 m_kgon_sum_with_arcs;

  Statistics m_statistics;
};

template <typename CPolygon_2>
void Minkowski_sum_construction_2<CPolygon_2>::clear()
{
  m_kgon_sum.clear();
  m_kgon_sum_outer_boundary.clear();

  m_kgon_diff.clear();
  m_kgon_diff_first_boundary.clear();

  m_kgon_offset_polygon.clear();
  m_kgon_offset_outer_boundary.clear();

  m_kgon_circles.clear();
  m_kgon_skipped_circles.clear();
  m_kgon_intersecting_circles.clear();
  
  m_kgon_sum_with_arcs.clear();
}


template <typename CPolygon_2>
void Minkowski_sum_construction_2<CPolygon_2>::compute_kgon_offset_polygon()
{
  // assumes that kgon_sum and kgon_circles (induced arcs) are already computed
  m_kgon_offset_polygon.clear();
  m_kgon_offset_outer_boundary.clear();

  CGAL::Timer timer;
  timer.start();

  // Detect self-intersections and replace self-intersecting
  // contours with the minkowski sum results
  if(!m_kgon_sum_with_arcs.empty())
  {
    // initialize handler data
    // TODO: clean-up global data by reorganizing classes
    Intersection_handler_2::smp_offset = &m_offset;
    Intersection_handler_2::smp_circles = &m_kgon_intersecting_circles;
    Intersection_handler_2::smp_ac_multimap = &m_arc_center_multimap;
    //Intersection_handler_2::smp_ac_freeze = 
    
    Intersection_handler_2 handler;
    handler.handle_self_intersections(m_kgon_sum_with_arcs, m_kgon_sum, m_kgon_offset_polygon);
  }
  else
  {
    LOG_DEBUG << "empty outer_contour" << std::endl;
  }

  timer.stop();

  m_kgon_offset_outer_boundary = m_kgon_offset_polygon.outer_boundary();

  m_statistics.kgon_offset_stats.time = timer.time();
  m_statistics.kgon_offset_stats.size = m_kgon_offset_outer_boundary.size();
  m_statistics.kgon_offset_stats.holes = m_kgon_offset_polygon.number_of_holes();

  // TODO write operator<<
  //if (verbose())
  {
    OUT_DEBUG << "The kgon offset polygon has "
              << m_statistics.kgon_offset_stats.size << " vertices, "
              << m_statistics.kgon_offset_stats.holes << " holes."<< std::endl;
    OUT_DEBUG << "The kgon offset polygon has " << m_kgon_intersecting_circles.size()
              << " frozen intersecting circles." << std::endl;
    OUT_DEBUG << "Kgon offset computation took "
              << m_statistics.kgon_offset_stats.time << " seconds." << std::endl;
  }

}

template <typename CPolygon_2>
void Minkowski_sum_construction_2<CPolygon_2>::compute_kgon_sum()
{
  m_kgon_sum.clear();
  m_kgon_sum_outer_boundary.clear();

  CGAL::Timer timer;

//  LOG_DEBUG << "counterclockwise:" 
//    << " polygon("  << m_polygon.is_counterclockwise_oriented()
//    << ") kgon(" << m_kgon.is_counterclockwise_oriented() << ")" << std::endl;

  timer.start();
  if(! m_polygon.is_empty() )
    m_kgon_sum = minkowski_sum_2 (m_polygon, m_kgon);
  timer.stop();

  m_kgon_sum_outer_boundary = m_kgon_sum.outer_boundary();

  m_statistics.kgon_sum_stats.time = timer.time();
  m_statistics.kgon_sum_stats.size = m_kgon_sum_outer_boundary.size();
  m_statistics.kgon_sum_stats.holes = m_kgon_sum.number_of_holes();

  // TODO write operator
  //if (verbose())
  {
    OUT_DEBUG << "The Kgon Sum polygon has "
              << m_statistics.kgon_sum_stats.size << " vertices in ob, and "
              << m_statistics.kgon_sum_stats.holes << " holes." << std::endl;
    OUT_DEBUG << "Kgon Sum computation took "
              << m_statistics.kgon_sum_stats.time << " seconds." << std::endl;
  }

}



template <typename CPolygon_2>
void Minkowski_sum_construction_2<CPolygon_2>::compute_kgon_induced_arcs()
{
  m_kgon_circles.clear();
  m_kgon_skipped_circles.clear();
  m_kgon_intersecting_circles.clear();
  m_arc_center_multimap.clear(); // helper structure of computed arcs
  m_kgon_sum_with_arcs.clear();

  CGAL::Timer timer;
  timer.start();

  Directed_slope_set kgon_slopes(m_kgon.edges_begin(), m_kgon.edges_end());
 
  // detect arcs with "freezed" waggly segments
  m_kgon_sum_outer_boundary = m_kgon_sum.outer_boundary();
  compute_polygon_arcs(m_kgon_sum_outer_boundary, kgon_slopes);
  
  LOG_DEBUG << m_kgon_sum.number_of_holes()
            << " HOLES" << std::endl;
   for(typename Kgon_sum_polygon_2::Hole_const_iterator hit = m_kgon_sum.holes_begin();
      hit != m_kgon_sum.holes_end();
      ++hit)
  {
    LOG_DEBUG << "HOLE:" << std::endl;
    compute_polygon_arcs(*hit, kgon_slopes);
  }
  
  timer.stop();

  m_statistics.kgon_circles_time = timer.time();
  m_statistics.kgon_sum_stats.circles = m_kgon_circles.size() + m_kgon_skipped_circles.size();

  //if (verbose())
  {
    OUT_DEBUG << "Kgon offset " << m_statistics.kgon_sum_stats.circles << " circles detection took "
              << m_statistics.kgon_circles_time << " seconds." << std::endl << std::endl;
  }

}

template <typename CPolygon_2>
void Minkowski_sum_construction_2<CPolygon_2>::compute_polygon_arcs(
  const Polygon_2& polygon, const Directed_slope_set& kgon_slopes)
{
 // int total_removed_vertices = 0;

  CGAL::set_pretty_mode(std::clog);
  LOG_DEBUG << "compute arcs of :" << std::endl;
  LOG_DEBUG << polygon << std::endl;
  LOG_DEBUG << std::endl;

  bool is_ccw_pgn = polygon.is_counterclockwise_oriented();
  bool is_ccw_kgon = m_kgon.is_counterclockwise_oriented();
  LOG_DEBUG << "counterclockwise:" 
    << " polygon("  << is_ccw_pgn
    << ") kgon(" << is_ccw_kgon << ")" << std::endl;

  CGAL::Orientation orientation = CGAL::COUNTERCLOCKWISE;
  Polygon_2 ccw_polygon = polygon;
  if(!is_ccw_pgn)
  {
    orientation = CGAL::CLOCKWISE;
//    LOG_DEBUG << "clockwise->counterclockwise polygon" << std::endl;
//    ccw_polygon.reverse_orientation(); 
  }

  m_circle_app.log_input_slopes();
  typename Directed_slope_set::const_iterator  sit;
  for (sit = kgon_slopes.begin(); sit != kgon_slopes.end(); ++sit)
  {
    Circle_app_2::log_slope("k___gon__slope", *sit);
    //LOG_DEBUG << "k___gon__slope (" << *sit << ") \t"  << sit->direction() << std::endl;
  }
  m_circle_app.log_circle_edges(ccw_polygon);

  int total_not_known = 0;

  Approximate_polygon_2 contour;

  typename Polygon_2::Edge_const_circulator ecirc = ccw_polygon.edges_circulator();

  // find first non-circular slope or first in the seq. of circ slopes
  // to avoid connecting the last circ slope to the first
  bool is_circle_slope = m_circle_app.is_circle_slope(*ecirc);

  typename Directed_slope_set::const_iterator slope_iter_end = kgon_slopes.end();
  typename Directed_slope_set::const_iterator slope_iter = slope_iter_end;

  if(is_circle_slope)
  { 
    // slope_iter = kgon_slopes.find(*ecirc);
    slope_iter = find_if(kgon_slopes.begin(), kgon_slopes.end(),
        Equal_directed_slope_2(*ecirc));
    // BUG Roza: comb.dat crash with dependent type
    LOG_DEBUG << "assert1(circle_slope in kgon): " << *ecirc << std::endl;
    assert(slope_iter != slope_iter_end);
    
    bool is_next_kgon_slope = false;
    do {
      ++ecirc;
      ++slope_iter; // circulate
      if(slope_iter == slope_iter_end)
      {
        slope_iter = kgon_slopes.begin();
      }
      is_circle_slope = m_circle_app.is_circle_slope(*ecirc);
      if(is_circle_slope)
      {
        is_next_kgon_slope = (slope_iter->direction() == (*ecirc).direction());
        Circle_app_2::log_slope("next_slope", *slope_iter);
        Circle_app_2::log_slope("next__edge", *ecirc);
      }
    } while( is_circle_slope && is_next_kgon_slope );
  }

  LOG_DEBUG << std::endl << "FIRST slope :" << *ecirc << std::endl;

  // iterate over all slopes
  // find sequences of circular ones 
  // with the same slope sequence as in the kgon
  bool start_new_circle = true; // should start a new circle
  bool computing_arc = false; // in the process of circle creation
  bool known_arc_position = false;

  Segment_2 first_segment;
  Segment_2 last_segment;
  Point_2 arc_center;

  slope_iter = slope_iter_end;
  Input_rational sqr_offset = m_offset * m_offset;

  typename Polygon_2::Edge_const_circulator ecirc_end = ecirc;
  do{

#define SAVE_MAP_LOOKUPS
#ifndef SAVE_MAP_LOOKUPS
    // if circle slope
    //    if new circle
    //      begin new circle                  (B)
    //    else (old circle)
    //      if !next slope
    //         complete and add old circle    (A)
    //         begin new circle               (B)
    //    compute arc center (if needed)      (C)
    //    get next preimage                   (D)
    // else (not circle slope)
    //     if old circle
    //        complete and add old circle     (A)
    //     add segment                        (E)

    bool is_kgon_slope = m_circle_app.is_circle_slope(*ecirc);
    if(is_kgon_slope)
    {
      // just to make sure
      typename Directed_slope_set::const_iterator sit = 
        find_if(kgon_slopes.begin(), kgon_slopes.end(),
        Equal_directed_slope_2(*ecirc));
    //  assert(sit != slope_iter_end);
      is_kgon_slope = (sit != slope_iter_end);
      if(!is_kgon_slope) {
        LOG_DEBUG << "BUG: not polygon ";
        Circle_app_2::log_slope("but not kgon slope", *ecirc);
      }
      // init slope iter 
      if(!computing_arc) slope_iter = sit;
    }

    // is current kgon slope goes in the expected direction,
    // that is the same direction as the next preimage slope
    bool is_next_kgon_slope = false;
    if(is_kgon_slope) 
    {
      assert(slope_iter != slope_iter_end);
      is_next_kgon_slope = (slope_iter->direction() == (*ecirc).direction());
      Circle_app_2::log_slope("next_slope", *slope_iter);
      Circle_app_2::log_slope("next__edge", *ecirc);
    }

    LOG_DEBUG << "computing_arc " << computing_arc 
      << " is_kgon_slope " << is_kgon_slope 
      << " is_next_kgon_slope " << is_next_kgon_slope << std::endl;

    // A. complete and add old circle
    if( computing_arc && 
        (!is_kgon_slope ||
        ( is_kgon_slope && !is_next_kgon_slope) ) )
    {
      LOG_DEBUG << "a) complete and add old circle " << *ecirc << std::endl;

      complete_and_add_circle(first_segment, last_segment, arc_center,
        known_arc_position, orientation, kgon_slopes, contour, total_not_known);

      start_new_circle = true;
      computing_arc = false;
      known_arc_position = false;
    }

    last_segment = *ecirc;

    LOG_DEBUG << "is_kgon_slope " << is_kgon_slope 
      << " start_new_circle " << start_new_circle 
      << " is_next_kgon_slope " << is_next_kgon_slope << std::endl;

    // B. begin new circle
    if( is_kgon_slope && 
       ( start_new_circle ||
         (!start_new_circle && !is_next_kgon_slope) ) )
    {
      LOG_DEBUG << "b) begin new circle " << *ecirc << std::endl;

      first_segment = last_segment;
      start_new_circle = false;
      computing_arc = true;

      // find and compare preimage to first segment
      //slope_iter = kgon_slopes.find(first_segment);
      slope_iter = find_if(kgon_slopes.begin(), kgon_slopes.end(),
        Equal_directed_slope_2(first_segment));

      /*
      if(slope_iter == slope_iter_end)
      {
        LOG_DEBUG << "!!! skipping " << last_segment << std::endl;

        computing_arc = false;
        continue;
      }

      LOG_DEBUG << "preimage " << *slope_iter << std::endl;
      */
      LOG_DEBUG << "assert2(circle_slope in kgon)" << first_segment << std::endl;
      assert(slope_iter != slope_iter_end);

      known_arc_position = false;
    }

    LOG_DEBUG << "is_kgon_slope " << is_kgon_slope 
//      << " compute_arc_center " << compute_arc_center 
      << " known_arc_position " << known_arc_position << std::endl;

    // C. compute arc center
    if(is_kgon_slope && !known_arc_position)
    {
      bool compute_arc_center = (slope_iter->to_vector() == last_segment.to_vector());
      if(compute_arc_center)
      {
        LOG_DEBUG << "c) compute arc center " << *ecirc << std::endl;

        known_arc_position = true;
        arc_center = last_segment.source() + Vector_2(slope_iter->source(), CGAL::ORIGIN);
        LOG_DEBUG << "arc_center " << arc_center << std::endl;
      }
    }

    // D. get next preimage
    if(is_kgon_slope)
    {
      LOG_DEBUG << "d) get next preimage " << *ecirc << std::endl;
      // move to the next slope on kgon
      ++slope_iter; // circulate
      if(slope_iter == slope_iter_end)
        slope_iter = kgon_slopes.begin();       
    }
    else
   // E. add segment
    {
      LOG_DEBUG << "e) add segment " << *ecirc << std::endl;
      // solid edge - add to circularized offset
      contour.push_back(X_monotone_curve_2 ((*ecirc).source(), (*ecirc).target()));

      // TODO: if convex angle with previous segment - add skipped circle

      computing_arc = false;
      known_arc_position = false;
    }

#else // SAVE_MAP_LOOKUPS

    // a version with less map look-ups:
    // if not creating circle
    //    if circle slope
    //       begin new circle                 (B)
    //       compute arc center (if needed)   (C)
    //       get next preimage                (D)
    //    else (polygon slope)
    //       add segment                      (E)
    // else (creating circle)
    //   if next slope
    //      compute arc center (if needed)    (C)
    //      get next preimage                 (D)
    //   else (not next slope)
    //      complete and add old circle       (A)
    //      if circle slope
    //         begin new circle                 (B)
    //         compute arc center (if needed)   (C)
    //         get next preimage                (D)
    //      else (polygon slope)
    //         add segment                      (E)
    //

    // does current edge have a kgon slope
    // Note: this can be an expensive map look-up, 
    // so it has to be computed only when absolutely necessary
    bool is_kgon_slope = false;

    // is current edge goes in the expected kgon slope direction,
    // that is the same direction as the next preimage slope
    bool is_next_kgon_slope = false;

    Circle_app_2::log_slope("nxt__edge", *ecirc);

    // check if new edge corresponds to next kgon slope
    if(slope_iter != slope_iter_end) // => computing_arc
    {
      Circle_app_2::log_slope("exp_slope", *slope_iter);
      is_next_kgon_slope = (slope_iter->direction() == (*ecirc).direction());
    }

    // check if new edge has a kgon slope (look it up if needed)
    if(is_next_kgon_slope)
    {
      is_kgon_slope = true;
    }
    else // !is_next_kgon_slope
    {
        // the expensive map look-up => init slope iter
        slope_iter = find_if(kgon_slopes.begin(), kgon_slopes.end(),
          Equal_directed_slope_2(*ecirc));
        is_kgon_slope = (slope_iter != slope_iter_end);

        if(slope_iter != slope_iter_end)
        {
          Circle_app_2::log_slope("nxt_slope", *slope_iter);
        }
    }

    LOG_DEBUG << "is_kgon_slope " << is_kgon_slope
      << " is_next_kgon_slope " << is_next_kgon_slope << std::endl;

    // A. complete and add old circle
    if( computing_arc && !is_next_kgon_slope )
    {
      LOG_DEBUG << "a) complete and add old circle " << *ecirc << std::endl;

      complete_and_add_circle(first_segment, last_segment, arc_center,
        known_arc_position, orientation, kgon_slopes, contour, total_not_known);

      computing_arc = false;
      known_arc_position = false;
    }

    last_segment = *ecirc;

    // B. begin new circle
    if( is_kgon_slope && (!computing_arc || (computing_arc && !is_next_kgon_slope) ) )
    {
      LOG_DEBUG << "b) begin new circle " << *ecirc << std::endl;

      first_segment = last_segment;
      computing_arc = true;
      known_arc_position = false;
    }

    // C. compute arc center
    if(is_kgon_slope && !known_arc_position)
    {
      bool compute_arc_center = (slope_iter->to_vector() == last_segment.to_vector());
      if(compute_arc_center)
      {
        LOG_DEBUG << "c) compute arc center " << *ecirc << std::endl;

        known_arc_position = true;
        arc_center = last_segment.source() + Vector_2(slope_iter->source(), CGAL::ORIGIN);
        LOG_DEBUG << "arc_center " << arc_center << std::endl;
      }
    }

    // D. get next preimage
    if(is_kgon_slope)
    {
      LOG_DEBUG << "d) get next preimage " << *ecirc << std::endl;

      // move to the next slope on kgon
      ++slope_iter; // circulate
      if(slope_iter == slope_iter_end)
        slope_iter = kgon_slopes.begin();       
    }
    else
   // E. add segment
    {
      LOG_DEBUG << "e) add segment " << *ecirc << std::endl;

      // solid edge - add to circularized offset
      contour.push_back(X_monotone_curve_2 ((*ecirc).source(), (*ecirc).target()));

      // TODO: if convex angle with previous segment - add skipped circle
      // (because there should have been arc between 2 consequent segs

      computing_arc = false;
      known_arc_position = false;
    }


    LOG_DEBUG << "is_kgon_slope " << is_kgon_slope 
//      << " compute_arc_center " << compute_arc_center 
      << " known_arc_position " << known_arc_position << std::endl;

      
#endif
  
  }while(++ecirc != ecirc_end);

  // complete the arc if needed
  if(computing_arc)
  {
    LOG_DEBUG << "a)) complete and add old circle " << *ecirc << std::endl;

    complete_and_add_circle(first_segment, last_segment, arc_center,
      known_arc_position, orientation, kgon_slopes, contour, total_not_known);
 }

  m_kgon_sum_with_arcs.push_back(contour);

  LOG_DEBUG << "Kgon offset with " << m_kgon_circles.size() << " replaced circles "
    << ", with " << m_kgon_skipped_circles.size() << " skipped circles" 
    << " and with " << total_not_known << " undetected circle positions." << std::endl;

//  LOG_DEBUG << "Removed "<< total_removed_vertices << " vertices." << std::endl;

}

template <typename CPolygon_2>
void Minkowski_sum_construction_2<CPolygon_2>::complete_and_add_circle(
  const Segment_2& first_segment,
  const Segment_2& last_segment,
  const Point_2& arc_center, 
  const bool& known_arc_position,
  const CGAL::Orientation& orientation,
  const Directed_slope_set& kgon_slopes,
  Approximate_polygon_2& contour,
  int& total_not_known)
{
  Input_rational sqr_offset = m_offset * m_offset;

  Point_2 start_point = first_segment.target();
  Point_2 end_point = last_segment.source();

  LOG_DEBUG << "->ARC: " << first_segment << std::endl;

  //static bool freeze_end_segments = true;
  static bool freeze_end_segments = false;
  // freeze_end_segments - by default is true
  // set to false if small non-convexities at arc ends are not important
  // TODO: choose if segment needs to be frozen according to actual geometry

  bool include_first = false; // do not incorporate first segment in the arc
  bool include_last = false; // do not incorporate last segment in the arc
  if(!freeze_end_segments)
  {
    // check if first segment is equal to preimage and should be included
    include_first = false;
    {
      typename Directed_slope_set::const_iterator slope_it = 
        find_if(kgon_slopes.begin(), kgon_slopes.end(),
        Equal_directed_slope_2(first_segment));
      assert(slope_it != kgon_slopes.end());
      include_first = (slope_it->to_vector() == first_segment.to_vector());
    }
    if(include_first)
    {
       start_point = first_segment.source();
    }

    // check if last segment is equal to preimage and should be included
    include_last = (first_segment == last_segment)? include_first: false;
    if(first_segment != last_segment)
    {
      typename Directed_slope_set::const_iterator slope_it = 
        find_if(kgon_slopes.begin(), kgon_slopes.end(),
        Equal_directed_slope_2(last_segment));
      assert(slope_it != kgon_slopes.end());
      include_last = (slope_it->to_vector() == last_segment.to_vector());
    }
    if(include_last)
    {
       end_point = last_segment.target();
    }
  }

  LOG_DEBUG << "include_first: " << include_first << std::endl;
  LOG_DEBUG << "include_last: " << include_last << std::endl;
  
  // add first (frosen) segment
  if(!include_first)
  {
    contour.push_back(X_monotone_curve_2 (first_segment.source(), first_segment.target()));
  }


  if(((first_segment == last_segment) && !include_first) || start_point == end_point)
  // not enough segments - skip the arc, add frosen segments
  {
    if(known_arc_position)
    {
      Circle_2 circle(arc_center, sqr_offset);
      m_kgon_skipped_circles.push_back(circle);
      LOG_DEBUG << "_ARC_: " << circle << std::endl;
    }
    else
    {
      ++total_not_known;
      LOG_DEBUG << "ARC: unknown" << std::endl;
    }
  }
  else
  // create an arc from start to end point, add frosen segments
  {
    LOG_DEBUG << "assert3(known_arc_position)" << known_arc_position << std::endl;
    assert(known_arc_position);
    // arc of the offset will always be counterclockwise
    Circle_2 circle(arc_center, sqr_offset, CGAL::COUNTERCLOCKWISE);//orientation);
    m_kgon_circles.push_back(circle);
    LOG_DEBUG << "ARC: " << circle << " from " << start_point << " to " << end_point << std::endl;

    add_x_monotone_arcs(start_point, end_point, circle, contour);
  }

  if(first_segment != last_segment)
  {
    LOG_DEBUG << "<-ARC: " << last_segment << std::endl;
    if(!known_arc_position)
    {
      // verify there are no skipped segments between the frozen segments
      LOG_DEBUG << "assert4(no skipped segments) " << start_point << " == " << end_point << std::endl;
      assert(start_point == end_point);
    }
    // add last (frosen) segment
    if(!include_last)
    {
      contour.push_back(X_monotone_curve_2 (last_segment.source(), last_segment.target()));
    }
  }
}


template <typename CPolygon_2>
void Minkowski_sum_construction_2<CPolygon_2>::add_x_monotone_arcs(
  const Point_2& start_point, 
  const Point_2& end_point,
  const Circle_2& circle, 
  Approximate_polygon_2& contour)
{

  Point_2 arc_center = circle.center();
  CGAL::Orientation orientation = circle.orientation();
  
  // if arc is not x-monotone split it in two
  // (in case of an offset it is not possible that arc is longer then 2 quadrants
  // therefore it is enough to check if y-coordinates of start and end are from
  // the same side of the arc center. if not - the arc should be split according 
  // to its orientation)
  bool is_start_below_center = start_point.y() < arc_center.y();
  bool is_start_above_center = start_point.y() > arc_center.y();
  bool is_end_below_center = end_point.y() < arc_center.y();
  bool is_end_above_center = end_point.y() > arc_center.y();
  bool from_below_to_above = is_start_below_center && is_end_above_center;
  bool from_above_to_below = is_start_above_center && is_end_below_center;

  if(from_below_to_above || from_above_to_below)
  {
    // have to split
    // in (center.x - radius, center.y) or in (center.x + radius, center.y)
    Point_2 mid_point = arc_center;
    Vector_2 offset_vector(m_offset, 0);
    if((orientation == CGAL::COUNTERCLOCKWISE && from_above_to_below) ||
      (orientation == CGAL::CLOCKWISE && from_below_to_above) )
    {
      mid_point = (arc_center - offset_vector);
    }
    else
    {
      mid_point = (arc_center + offset_vector);   
    }
    LOG_DEBUG << "ARC: " << circle << " split in " << mid_point << std::endl;

    Arc_point_2 arc_start(start_point.x(), start_point.y());
    Arc_point_2 arc_mid(mid_point.x(), mid_point.y());    
    Arc_point_2 arc_end(end_point.x(), end_point.y());

    Contour_2_iterator cit;
    cit = contour.insert(X_monotone_curve_2 (circle, arc_start, arc_mid, orientation));
    m_arc_center_multimap.insert(Arc_center_val_2(arc_center,
      Arc_data_2(start_point, mid_point, circle, contour, cit)) );
    cit = contour.insert(X_monotone_curve_2 (circle, arc_mid, arc_end, orientation));
    m_arc_center_multimap.insert(Arc_center_val_2(arc_center,
      Arc_data_2(mid_point, end_point, circle, contour, cit)) );
  }
  else
  {
    // no need to split
    Arc_point_2 arc_start(start_point.x(), start_point.y());   
    Arc_point_2 arc_end(end_point.x(), end_point.y());
    Contour_2_iterator cit;
    cit = contour.insert(X_monotone_curve_2 (circle, arc_start, arc_end, orientation));
    m_arc_center_multimap.insert(Arc_center_val_2(arc_center,
      Arc_data_2(start_point, end_point, circle, contour, cit)) );
  }
}

template <typename CPolygon_2>
void Minkowski_sum_construction_2<CPolygon_2>::compute_kgon_diff()
{
  m_kgon_diff.clear();
  m_kgon_diff_first_boundary.clear();
  /*  CGAL::set_pretty_mode(std::clog);
  LOG_DEBUG << "molding P:" << m_polygon.is_counterclockwise_oriented()  << std::endl;
  LOG_DEBUG << m_polygon << std::endl;
  LOG_DEBUG << std::endl;*/

  CGAL::Timer                  timer;

  timer.start();
  
  // compute the complement of the polygon and cut it in two parts
  Polygon_2 left_mold, right_mold;
  mold_polygon(m_polygon, m_offset, left_mold, right_mold);

/*  LOG_DEBUG << "left: " << left_mold.is_counterclockwise_oriented()  << std::endl;
  LOG_DEBUG << left_mold << std::endl;
  LOG_DEBUG << std::endl;
  LOG_DEBUG << "right:" << right_mold.is_counterclockwise_oriented()  << std::endl;
  LOG_DEBUG << right_mold << std::endl;
  LOG_DEBUG << std::endl;*/

  // offset the complement of the polyon by offseting and then joining the two parts
  Kgon_sum_polygon_2 inset_r_left = minkowski_sum_2 (left_mold, m_kgon);
  Kgon_sum_polygon_2 inset_r_right = minkowski_sum_2 (right_mold, m_kgon);

/*  LOG_DEBUG << "left+:"  << std::endl;
  LOG_DEBUG << inset_r_left << std::endl;
  LOG_DEBUG << std::endl;
  LOG_DEBUG << "right+:" << std::endl;
  LOG_DEBUG << inset_r_right << std::endl;
  LOG_DEBUG << std::endl;*/

  Kgon_sum_polygon_2 inset_r;
  bool success = join(inset_r_left, inset_r_right, inset_r);

  timer.stop();

  m_statistics.kgon_diff_stats.time = timer.time();
  m_statistics.kgon_diff_stats.size = inset_r.outer_boundary().size();
  m_statistics.kgon_diff_stats.holes = inset_r.number_of_holes();

  LOG_DEBUG << "reverse " << m_statistics.kgon_diff_stats.holes << " holes" << std::endl;
  
  // TODO: in case of holes in input we need to compute their offset
  // and join them here before the second complement
  
//  assert(inset_r.number_of_holes() != 0);
//  P = *inset_r.holes_begin();
//  P.reverse_orientation();

//  Polygon_2 empty_boundary;
//  Kgon_sum_polygon_2 inset_r_unbounded(empty_boundary, inset_r.holes_begin(), inset_r.holes_end());
  // the inverse of the holes is the resulting inset

  m_kgon_diff_first_boundary.clear();
  m_kgon_diff.clear();
  typename Kgon_sum_polygon_2::Hole_iterator hit = inset_r.holes_begin();
  if( hit != inset_r.holes_end())
  {
      m_kgon_diff_first_boundary = *hit;
      m_kgon_diff_first_boundary.reverse_orientation();
      m_statistics.kgon_diff_stats.size = m_kgon_diff_first_boundary.size();

      do
      {
        Approximate_polygon_2 pgn;
        overwrite(*hit, pgn);
        pgn.reverse_orientation();
        m_kgon_diff.insert(pgn);
        ++hit;
      }
      while(hit != inset_r.holes_end());
  }

  //if (verbose())
  {
    OUT_DEBUG << "The complement Kgon Sum"
      << " has " << (m_statistics.kgon_diff_stats.holes == 0? "ob": "first hole")
      << " of size " << m_statistics.kgon_diff_stats.size
      << ", with " << m_statistics.kgon_diff_stats.holes << " holes." << std::endl;
    OUT_DEBUG << "The complement Kgon Sum"
      << " was computed in " << m_statistics.kgon_diff_stats.time << " seconds." << std::endl;
  }
  
  if(m_statistics.kgon_diff_stats.holes > 1)
  {
    LIM_DEBUG << "!!! inset multiple contours not handled yet !!!" << std::endl;
  }
  else
  if(m_statistics.kgon_diff_stats.holes == 0)
  {
    LIM_DEBUG << "!!! inset is empty - not handled yet !!!" << std::endl;
  }

//   if(!m_kgon_diff.is_empty())
//   {
//     Approximate_offset_polygon_list_2 pwh_list;
//     m_kgon_diff.polygons_with_holes(std::back_inserter(pwh_list));
//     Approximate_offset_polygon_2 first_polygon_with_hole = *(pwh_list.begin());
//     Approximate_polygon_2 first_boundary = first_polygon_with_hole.outer_boundary();
//     overwrite(first_boundary, m_kgon_diff_first_boundary);
//   }

}

// compute bounding box of at least radius away from polygon
// place a polygon as a hole in it
// cut the polygon with a single hole to two simple polygons (left and right)
// by cutting through hole via upper and lower vertex
template <typename CPolygon_2>
void Minkowski_sum_construction_2<CPolygon_2>::mold_polygon(const Polygon_2& i_polygon,
                            const Input_rational& radius,
                            Polygon_2& o_polygon_left,
                            Polygon_2& o_polygon_right)
{
  // get outer box
  Bbox_2 polygon_box = i_polygon.bbox();
  double d_radius_ceil = ceil(to_double(radius));
  Input_rational x_min = floor(polygon_box.xmin()) - d_radius_ceil;
  Input_rational y_min = floor(polygon_box.ymin()) - d_radius_ceil;
  Input_rational x_max = ceil(polygon_box.xmax()) + d_radius_ceil;
  Input_rational y_max = ceil(polygon_box.ymax()) + d_radius_ceil;

//  LOG_DEBUG << "box2 " << " (" << x_min << ", " << y_min << ") - (" << x_max << ", " << y_max << ")" << std::endl;

  Polygon_2 hole = i_polygon;
  hole.reverse_orientation();

  typedef typename Polygon_2::Vertex_iterator Vertex_iterator;

  // get top and bottom vertex
  Vertex_iterator top_vertex = hole.top_vertex();
  Vertex_iterator bottom_vertex = hole.bottom_vertex();

  Input_rational y_top = top_vertex->y();
  Input_rational y_bottom = bottom_vertex->y();
  Input_rational x_top = top_vertex->x();
  Input_rational x_bottom = bottom_vertex->x();

//  LOG_DEBUG << "pgn" << " top(" << x_top << ", " << y_top << ") - bottom(" << x_bottom << ", " << y_bottom << ")" << std::endl;

  // box points in counterclockwise order
  /*  two half-boxes layout:
      ymax: xmin ymax               xtop ymax          xmax ymax

                                    xtop ytop

                           xbottom ybottom

      ymin: xmin ymin      xbottom ymin                xmax ymin
  */


  static int HALF_BOX_SIZE = 5;
  // separate the box to left and right parts
  Point_2 box_points_l[] = {
    Point_2(x_top, y_top),
    Point_2(x_top, y_max),
    Point_2(x_min, y_max),
    Point_2(x_min, y_min),
    Point_2(x_bottom, y_min)};

  Point_2 box_points_r[] = {
    Point_2(x_bottom, y_bottom),
    Point_2(x_bottom, y_min),
    Point_2(x_max, y_min),
    Point_2(x_max, y_max),
    Point_2(x_top, y_max)};

  // add box vertices to the output
  o_polygon_left = Polygon_2(box_points_l, box_points_l + HALF_BOX_SIZE);
  o_polygon_right = Polygon_2(box_points_r, box_points_r + HALF_BOX_SIZE);

  // add polygon vertices to the output
  if(top_vertex < bottom_vertex)
  {
    // bottom to top (left) needs 2 insertions
    o_polygon_left.insert(o_polygon_left.vertices_begin() , hole.vertices_begin(), top_vertex);
    o_polygon_right.insert(o_polygon_right.vertices_begin() , top_vertex, bottom_vertex);
    o_polygon_left.insert(o_polygon_left.vertices_begin() , bottom_vertex, hole.vertices_end());
  }
  else
  {
    // top to bottom (right) needs 2 insertions
    assert(top_vertex > bottom_vertex);
    o_polygon_right.insert(o_polygon_right.vertices_begin() , hole.vertices_begin(), bottom_vertex);
    o_polygon_left.insert(o_polygon_left.vertices_begin() , bottom_vertex, top_vertex);
    o_polygon_right.insert(o_polygon_right.vertices_begin() , top_vertex, hole.vertices_end());
  }

//  LOG_DEBUG << "left mold size " << o_polygon_left.size()<< std::endl;
//  LOG_DEBUG << "right mold size " << o_polygon_right.size()<< std::endl;

}


} // namespace CGAL


#endif // Minkowski_sum_construction_2_h
