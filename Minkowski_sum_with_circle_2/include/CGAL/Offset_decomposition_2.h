// Offset_decomposition_2.h

#ifndef Offset_decomposition_2_h
#define Offset_decomposition_2_h

#include <CGAL/Offset_statistics_2.h>
#include <CGAL/Circle_approximation_2.h>
#include <CGAL/Minkowski_sum_construction_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <list>

namespace CGAL {


/*! \classf
 * A class implementing decomposition of polygon offsets 
 * via offset/inset constructions with govering statistics. 
 * The inputs are defined by caller application,
 * the class keeps them as _references_. The outputs are defined 
 * by the class and caller application may keep a reference to them.
 */

template <typename CPolygon_2>
class Offset_decomposition_2 : public Offset_statistics_2<CPolygon_2>
{
public:

  typedef Offset_statistics_2<CPolygon_2> Base;
  typedef typename Base::Types Types;
  typedef typename Types::Polygon_2 Polygon_2;
  typedef typename Types::Polygon_traits_2 Polygon_traits_2;  
  typedef typename Types::Rat_traits_2 Rat_traits_2;
  typedef typename Types::Input_rational Input_rational;
  typedef Minkowski_sum_construction_2<CPolygon_2> Min_sum_construction_2;
  typedef Offset_construction_2<CPolygon_2> Off_construction_2;

  typedef typename Polygon_traits_2::Point_2   Point_2;
//  typedef typename Polygon_traits_2::Segment_2   Segment_2;
//  typedef typename Polygon_traits_2::Direction_2   Direction_2;
//  typedef typename CGAL::Aff_transformation_2<Polygon_traits_2>  Transformation_2;
  
  
  typedef typename Types::Exact_alg_kernel Exact_alg_kernel;
  typedef typename Exact_alg_kernel::Segment_2 Exact_segment_2;
  typedef typename std::list<Exact_segment_2> Exact_segment_list_2; 

  typedef typename Base::Offset_statistics Offset_statistics;
  typedef typename Base::Polygon_statistics Polygon_statistics;

  typedef typename Types::Conic_traits_2 Conic_traits_2;

  typedef typename Types::Exact_rational Exact_rational;
  typedef typename Types::Exact_algebraic Exact_algebraic;
  typedef typename Types::Rat_polygon_2 Rat_polygon_2;
  typedef typename Types::Rat_polygon_with_holes_2 Rat_polygon_with_holes_2;
  typedef typename Types::Exact_offset_polygon_2 Exact_offset_polygon_2; 
  typedef typename Types::Exact_polygon_2 Exact_polygon_2;
  typedef typename Types::Exact_polygon_list_2 Exact_polygon_list_2;
  
  typedef typename CGAL::Circle_approximation_2<CPolygon_2>    Circle_app_2;
  typedef typename CGAL::Polygon_with_holes_2<Polygon_traits_2>  Kgon_sum_polygon_2;
  typedef typename CGAL::General_polygon_set_2<Rat_traits_2> Polygon_with_hole_set_2;
  typedef typename std::list<Kgon_sum_polygon_2>                 Kgon_sum_polygon_list_2;
  
  struct Statistics
  {
      // kgon sum statistics
      Polygon_statistics outer_core_input;
      Polygon_statistics inner_core_input;

      Offset_statistics outer_core_stats;
      Offset_statistics inner_core_stats;

      Polygon_statistics outer_eps_app_stats;
      Polygon_statistics inner_eps_app_stats;
      
      Polygon_statistics outer_core_app_stats;
      Polygon_statistics inner_core_app_stats;
      
      Polygon_statistics outer_core_ms_stats;
      Polygon_statistics inner_core_ms_stats;

      Polygon_statistics outer_core_msr_stats;
      Polygon_statistics inner_core_msr_stats;
      
      Polygon_statistics outer_diff_stats;
      Polygon_statistics inner_diff_stats;

      Polygon_statistics outer_approximability_stats;
      Polygon_statistics inner_approximability_stats;
      
      double approximability_time;
  };

  /// input parameters stored _as references_
  Offset_decomposition_2(
    const Polygon_2& i_polygon,
    const Input_rational& i_offset,
    const Input_rational& i_eps,
    const Input_rational& i_delta) :
    m_polygon(i_polygon),
    m_offset(i_offset),
    m_epsilon(i_eps),
    m_delta(i_delta),
    // minkowski sum construction parameters and classes:
    m_ms_offset(i_eps),
    m_ms_epsilon(i_delta),
    m_ms_inner(true),
    m_circle_app(m_polygon, m_ms_offset, m_ms_epsilon, m_ms_inner),
    m_ms_kgon(m_circle_app.kgon()),
    m_ms_polygon(m_polygon),
    m_ms_constructor(m_ms_polygon, m_ms_offset, m_ms_epsilon, m_circle_app),
    m_offset_constructor(m_inner_core_app, m_offset, m_epsilon)
  {
    m_circle_app.verbose(false);
    m_ms_constructor.verbose(false);
    m_offset_constructor.verbose(false);
  }

  ~Offset_decomposition_2(void)
  {
  }

  const Statistics& statistics() const { return m_statistics; }

  const Exact_polygon_2& outer_core() const { return m_outer_core; }

  const Exact_polygon_2& inner_core() const { return m_inner_core; }

  const bool& approximability() const { return m_approximability; }

  const Exact_segment_list_2& approximable_data() const { return m_approximable_data; }
  const Exact_segment_list_2& non_approximable_data() const { return m_non_approximable_data; }

  const Polygon_2& outer_core_app() const { return m_outer_core_app; }
  const Polygon_2& inner_core_app() const { return m_inner_core_app; }

  const Polygon_2& ms_kgon() const { return m_ms_kgon; }

  const Kgon_sum_polygon_2& outer_eps_app() const { return m_outer_eps_app; }
  const Kgon_sum_polygon_2& inner_eps_app() const { return m_inner_eps_app; }

  const Kgon_sum_polygon_2& outer_core_ms() const { return m_outer_core_ms; }
  const Kgon_sum_polygon_2& inner_core_ms() const { return m_inner_core_ms; }

  const Kgon_sum_polygon_2& outer_core_msr() const { return m_outer_core_msr; }
  const Kgon_sum_polygon_2& inner_core_msr() const { return m_inner_core_msr; }

  const Kgon_sum_polygon_2& outer_core_diff() const { return m_outer_core_diff; }
  const Kgon_sum_polygon_2& inner_core_diff() const { return m_inner_core_diff; }

  const bool& inner_approximability() const { return m_inner_app; }
  const bool& outer_approximability() const { return m_outer_app; }

  const bool& approximability_app() const { return m_approximability_app; }

  // exact offset of the computed core (approximate source polygon)
  const Exact_polygon_2& exact_outer_boundary() const { return m_offset_constructor.exact_outer_boundary(); }
  const Exact_offset_polygon_2& exact_offset_polygon() const { return m_offset_constructor.exact_offset_polygon(); }

  void clear();

  // inset/ P-tilda r - eps
  void compute_outer_core();
  // inset/ P-tag-tilda r + eps (not implemented yet)
  void compute_inner_core();

  /// decide if (convex) input polygon close up to epsilon to a given offset of some polygon
  void compute_approximability();

//  // outer: outer eps-offset delta-approximation of input polygon (Q)
//  void compute_outer_core_off() { compute_core_off(false); }
//  // inner: inner eps-offset delta-approximation of input polygon (Q)
//  void compute_inner_core_off() { compute_core_off(true); }

  // outer: outer eps-offset delta-approximation of input polygon (Q)
  void compute_outer_eps_app() { compute_eps_app(false); }
  // inner: inner eps-offset delta-approximation of input polygon (Q)
  void compute_inner_eps_app() { compute_eps_app(true); }

  // outer: inner r-inset delta-approximation of outer eps-offset delta-approximation of input polygon (Q)
  void compute_outer_core_app() { compute_core_app(false); }
  // inner: outer r-inset delta-approximation of inner eps-offset delta-approximation of input polygon (Q)
  void compute_inner_core_app() { compute_core_app(true); }

  // outer: outer r+eps-offset delta-approximation of outer core approximation
  void compute_outer_core_ms() { compute_core_ms(false); }
  // inner: inner r+eps-offset delta-approximation of inner core approximation
  void compute_inner_core_ms() { compute_core_ms(true); }

  // outer: outer r-offset delta-approximation of outer core approximation
  void compute_outer_core_msr() { compute_core_msr(false); }
  // inner: inner r-offset delta-approximation of inner core approximation
  void compute_inner_core_msr() { compute_core_msr(true); }

  // outer: outer r+eps-offset delta-approximation diff with input
  void compute_outer_core_diff() { compute_core_diff(false); }
  // inner: inner r+eps-offset delta-approximation diff with input
  void compute_inner_core_diff() { compute_core_diff(true); }

  // outer: is outer diff empty?
  void decide_outer_approximability() { m_outer_app = decide_approximability_app(false); }
  // inner: is inner diff empty?
  void decide_inner_approximability() { m_inner_app = decide_approximability_app(true); }

  /// decide with delta-precision if input polygon close up to epsilon to a given offset of some polygon.
  void compute_approximability_app();

  void save_inner_core()
  {
    std::ofstream ofs("inner_core.dat");
    ofs << m_inner_core_app;
    //std::ofstream ofs("inner_core_2.dat");
    // add point-point printing
  }

  void compute_exact_offset()
  {
    m_offset_constructor.compute_exact_offset_polygon();
  }
  
  void update_slopes() {  m_circle_app.update_slopes(); }

private:

  // eps-offset delta-approximation of input polygon (Q)
  void compute_eps_app(const bool& is_inner_app);
  // r-inset delta-approximation of eps-offset of input polygon (Q)
  void compute_core_app(const bool& is_inner_app);
  // r+eps-offset delta-approximation of core approximation
  void compute_core_ms(const bool& is_inner_app);
  // r-offset delta-approximation of core approximation
  void compute_core_msr(const bool& is_inner_app);  
  // diff of input polygon (Q) and core approximation offset
  void compute_core_diff(const bool& is_inner_app);  

  bool decide_approximability_app(const bool& is_inner_app);
  
  static bool decide_approximability_convex(
    const Input_rational& epsilon, 
    const Input_rational& radius,
    const Rat_polygon_2& Q, 
    const Exact_polygon_2& P/*_tilda*/, 
    Exact_segment_list_2& o_yes_data, Exact_segment_list_2& o_no_data);


    const Polygon_2& m_polygon;
    const Input_rational& m_offset;
    const Input_rational& m_epsilon;
    const Input_rational& m_delta;

    Polygon_2 m_ms_polygon;
    Input_rational m_ms_offset;
    Input_rational m_ms_epsilon;
    const Polygon_2& m_ms_kgon;
    bool m_ms_inner;
    Circle_app_2 m_circle_app;
    Min_sum_construction_2 m_ms_constructor;
    Off_construction_2 m_offset_constructor;
     
    Exact_polygon_list_2 m_outer_cores;
    Exact_polygon_2 m_outer_core;

    Exact_polygon_list_2 m_inner_cores;
    Exact_polygon_2 m_inner_core;

    //Rat_polygon_with_holes_2 m_rat_polygon_with_holes;
    Rat_polygon_2 m_rat_oc_polygon;
    Rat_polygon_2 m_rat_ic_polygon;
    
    bool m_approximability;
    Exact_segment_list_2 m_approximable_data;
    Exact_segment_list_2 m_non_approximable_data;

    Kgon_sum_polygon_2 m_inner_eps_app;
    Kgon_sum_polygon_2 m_outer_eps_app;
    
    Polygon_2 m_inner_core_app;
    Polygon_2 m_outer_core_app;
    
    Kgon_sum_polygon_2 m_inner_core_ms;
    Kgon_sum_polygon_2 m_outer_core_ms;

    Kgon_sum_polygon_2 m_inner_core_msr;
    Kgon_sum_polygon_2 m_outer_core_msr;  
    
    Kgon_sum_polygon_2 m_inner_core_diff;
    Kgon_sum_polygon_2 m_outer_core_diff;

    bool m_inner_app;
    bool m_outer_app;    
    bool m_approximability_app;

    Statistics m_statistics;
};


template <typename CPolygon_2>
void Offset_decomposition_2<CPolygon_2>::clear()
{
  m_outer_cores.clear();
  m_outer_core.clear();

  m_inner_cores.clear();
  m_inner_core.clear();

  m_rat_oc_polygon.clear();
  m_rat_ic_polygon.clear();

  m_approximable_data.clear();
  m_non_approximable_data.clear();

  m_inner_eps_app.clear();
  m_outer_eps_app.clear();
  m_inner_core_app.clear();
  m_outer_core_app.clear();

  m_inner_core_ms.clear();
  m_outer_core_ms.clear();
  m_inner_core_msr.clear();
  m_outer_core_msr.clear();

  m_inner_core_diff.clear();
  m_outer_core_diff.clear();

}


// inset/ P-tilda r - eps
template <typename CPolygon_2>
void Offset_decomposition_2<CPolygon_2>::compute_outer_core()
{
  m_outer_cores.clear();
  m_outer_core.clear();

  m_rat_oc_polygon.clear();

  if(m_polygon.is_empty()) return;

  /// in case of the convex Q, P-tilda is equivalent to Inset(r-eps)

  // in case of non-convex Q, Inset(r-eps) is fully
  // contained in P-tilda (the difference is in the
  // arcs created by concave vertices of Q)

  Exact_rational radius_minus_eps = type_cast<Exact_rational>(m_offset - m_epsilon);
//  LOG_DEBUG << "exact inset " << radius_minus_eps << std::endl;
  {
      CGAL::Timer timer;
      timer.start();

      copy(m_polygon, m_rat_oc_polygon);
      timer.stop();

      m_statistics.outer_core_input.time = timer.time();
      m_statistics.outer_core_input.size = m_rat_oc_polygon.size();

      //if (verbose())
      {
        OUT_DEBUG << "Conversion of input polygon of size " << m_statistics.outer_core_input.size
          << " computed in " << m_statistics.outer_core_input.time << std::endl;
      }
  }

  CGAL::Timer timer;

  //LOG_DEBUG << "before inset" << std::endl;
  timer.start();
  Conic_traits_2 conic_traits;
  inset_polygon_2 (m_rat_oc_polygon, radius_minus_eps, conic_traits,
                      std::back_inserter(m_outer_cores));
  //LOG_DEBUG << "after inset" << std::endl;
  timer.stop();

  m_statistics.outer_core_stats.time = timer.time();
  m_statistics.outer_core_stats.size = 0;

  //if (verbose())
  {
    OUT_DEBUG << "The reconstruction inset has "
            << m_outer_cores.size() << " polygons." << std::endl;
  }

  if(!m_outer_cores.empty())
  {
      m_outer_core = *(m_outer_cores.begin());

      m_statistics.outer_core_stats.size = m_outer_core.size();

      //if (verbose())
      {
        OUT_DEBUG << "the first with "
          << m_statistics.outer_core_stats.size << " vertices."<< std::endl;
      }
  }

  //if (verbose())
  {
    OUT_DEBUG << "Reconstruction inset computation took "
      << m_statistics.outer_core_stats.time << " seconds." << std::endl << std::endl;
  }
}

// approximability
template <typename CPolygon_2>
void Offset_decomposition_2<CPolygon_2>::compute_approximability()
{
  LOG_DEBUG << "compute_approximability()" << std::endl;
  m_approximable_data.clear();
  m_non_approximable_data.clear();

  // assumes outer core has been computed  
  if(m_rat_oc_polygon.is_empty()) return;

  if(!m_rat_oc_polygon.is_convex())
  {
    LIM_DEBUG << "!!! Approximability for non convex input is not implemented yet !!!" << std::endl;
    return;
  }
  else
  {
    m_approximability = decide_approximability_convex(
      m_epsilon, m_offset, m_rat_oc_polygon, m_outer_core,
      m_approximable_data, m_non_approximable_data);
  }

/*
  // TODO: move to the caller
  Input_rational eps_10 = m_epsilon / 10;
  Input_rational eps_100 = m_epsilon / 100;
  bool eps_10_app, eps_100_app;

  if(m_forward_construction)
  {
    eps_10_app = decide_approximability_convex(
        eps_10, m_offset, m_rat_polygon, m_outer_core);
    eps_100_app = decide_approximability_convex(
        eps_100, m_offset, m_rat_polygon, m_outer_core);
  }
*/
  //if (verbose())
  {
    OUT_DEBUG << "Polygon is " << (m_approximability? "": "not ") << "approximable "
    << "with offset " << m_offset << " and epsilon " << m_epsilon << std::endl;
  }
}

// approximate approximability
template <typename CPolygon_2>
void Offset_decomposition_2<CPolygon_2>::compute_approximability_app()
{
  LOG_DEBUG << "compute_approximability_app()" << std::endl;

  m_approximability_app = false;
  if (m_inner_app)
  {
    m_approximability_app = m_inner_app;
  }
  else
  {
    if (!m_outer_app)
    {
      m_approximability_app = m_outer_app;
    }
    else
    {
      LOG_DEBUG << "!!! Approximability is not decided for offset " << m_offset
        << " and epsilon " << m_epsilon
        << " with delta " << m_delta  << " !!!"<< std::endl;
    }
  }

  //if (verbose())
  {
    OUT_DEBUG << "Polygon is " << ((m_inner_app || !m_outer_app)?
      (m_approximability_app? "": "not ") : "undecided if ")
      << "approximable " << "with offset " << m_offset << " and epsilon " << m_epsilon
      << " with delta " << m_delta << std::endl;

    OUT_DEBUG << "Approximability " << (m_inner_app? "one-sided":"two-sided") << " decision procedure took "
      << m_statistics.inner_approximability_stats.time + (m_inner_app? 0: m_statistics.outer_approximability_stats.time)
      << " seconds." << std::endl << std::endl;
  }
}

// inset/ P-tag-tilda r + eps
template <typename CPolygon_2>
void Offset_decomposition_2<CPolygon_2>::compute_inner_core()
{
  m_inner_cores.clear();
  m_inner_core.clear();

  m_rat_ic_polygon.clear();

  if(m_polygon.is_empty()) return;

  /// Inset(r+eps)

  Exact_rational radius_plus_eps = type_cast<Exact_rational>(m_offset + m_epsilon);
//  LOG_DEBUG << "exact inset " << radius_minus_eps << std::endl;
  {
      CGAL::Timer timer;
      timer.start();

      copy(m_polygon, m_rat_ic_polygon);
      timer.stop();

      m_statistics.inner_core_input.time = timer.time();
      m_statistics.inner_core_input.size = m_rat_ic_polygon.size();

      //if (verbose())
      {
        OUT_DEBUG << "Conversion of input polygon of size " << m_statistics.inner_core_input.size
          << " computed in " << m_statistics.inner_core_input.time << std::endl;
      }
  }

  CGAL::Timer timer;

  timer.start();
  Conic_traits_2 conic_traits;
  inset_polygon_2 (m_rat_ic_polygon, radius_plus_eps, conic_traits,
                      std::back_inserter(m_inner_cores));
  timer.stop();

  m_statistics.inner_core_stats.time = timer.time();
  m_statistics.inner_core_stats.size = 0;

  //if (verbose())
  {

    OUT_DEBUG << "The reconstruction inset has "
            << m_inner_cores.size() << " polygons." << std::endl;
  }

  if(!m_inner_cores.empty())
  {
      m_inner_core = *(m_inner_cores.begin());

      //if (verbose())
      {
        OUT_DEBUG << "the first with "
            << m_inner_core.size() << " vertices."<< std::endl;
      }
  }

  //if (verbose())
  {
    OUT_DEBUG << "Reconstruction inset computation took "
            << timer.time() << " seconds." << std::endl << std::endl;
  }
}

/*!
 * Check for the convex Q if it is approximable up to epsilon
 * by offset radius from P (where P is P-tilda, i.e. r-eps inset of Q).
 * It is enough to check that any vertex of (convex) Q is in an r+eps
 * neighbourhood of some vertex of P (see the paper for the proof).
 * Since vertices of Q and P go are supposed to go almost "in parallel",
 * the test can be done in linear time, by walking on both polygons and
 * testing vertex distances.
 */
template <typename CPolygon_2>
bool Offset_decomposition_2<CPolygon_2>::decide_approximability_convex(
    const Input_rational& epsilon, const Input_rational& radius,
    const Rat_polygon_2& Q, const Exact_polygon_2& P/*_tilda*/,
    Exact_segment_list_2& o_yes_data, Exact_segment_list_2& o_no_data)
{
    typedef typename Rat_polygon_2::Vertex_iterator QVertex_iterator;
    typedef typename Rat_polygon_2::Point_2 QPoint_2;

//   typedef typename Exact_polygon_2::General_polygon_2::Vertex_iterator PVertex_iterator;
//   typedef typename Exact_polygon_2::General_polygon_2 Exact_general_polygon_2;
#ifdef _MSC_VER
    typedef typename Exact_polygon_2::Curve_const_iterator PCurve_iterator;
    typedef typename Exact_polygon_2::Point_2 PPoint_2;
#else
    typedef typename Exact_polygon_2::General_polygon_2::Curve_const_iterator PCurve_iterator;
    typedef typename Exact_polygon_2::General_polygon_2::Point_2 PPoint_2;
#endif

//    typedef typename PPoint_2::Kernel::FT PDist_2;
    typedef Exact_algebraic PDist_2;
    // find closest vertex of p to the first vertex of q
    // check that distance is less then r+eps
    bool approximability = false;

    // assert(!Q.empty() && !P.empty());

    if(Q.is_empty()) return false;
    if(P.is_empty()) return false;

    CGAL::Orientation p_orientation = P.orientation();
    CGAL::Orientation q_orientation = Q.orientation();
    Rat_polygon_2 oriented_Q = Q;
    if(p_orientation != q_orientation)
    {
      LOG_DEBUG << "matching Q orientation to P " << (p_orientation == CGAL::CLOCKWISE? "clockwise": "counterclockwise") << std::endl;
      oriented_Q.reverse_orientation();
    }
      
    
    // assumption: both polygons have the same orientation
    // that is are counterclockwise oriented

    // we have a polygon Q (with rational coordinates)
    // and a polygon P (with possible sqrt in coordinates)

    // we take an arbitrary vertex of (convex) Q and
    // find a closest vertex of (convex) P.
    // 
    PDist_2 in_dist = type_cast<PDist_2>(radius + epsilon);
    in_dist *= in_dist;   // (r + eps)^2
    LOG_DEBUG << "in_dist: " << in_dist << std::endl;

    QVertex_iterator q_start = oriented_Q.vertices_begin();
    QVertex_iterator qit = q_start;

    QPoint_2 pnt = *qit;
    PPoint_2 q_pnt = PPoint_2(
                     type_cast<PDist_2>(pnt.x()),
                     type_cast<PDist_2>(pnt.y()));

    PCurve_iterator p_start = P.curves_begin();
    PCurve_iterator pit = p_start;
    PCurve_iterator p_curr = pit;

    PPoint_2 p_pnt = pit->source();
    PDist_2 min_dist = CGAL::squared_distance(q_pnt, p_pnt);
    ++pit;

    for(; pit != P.curves_end(); ++pit)
    {
        p_pnt = pit->source();

        PDist_2 p_dist = CGAL::squared_distance(q_pnt, p_pnt);
        if (p_dist < min_dist)
        {
            min_dist = p_dist;
            p_curr = pit;
        }
    }

    Exact_segment_2 curr_segment = Exact_segment_2(q_pnt, p_curr->source());
    bool in_range = (min_dist <= in_dist);

    if (!in_range)
    {
       LOG_DEBUG << "no app: did not find a match for " << *qit << std::endl;
       o_no_data.push_back(curr_segment);
//       return approximability;
    }
    else
    { 
       LOG_DEBUG << "app: found a match for " << *qit << std::endl;
       o_yes_data.push_back(curr_segment);
    }

    PDist_2 max_dist = min_dist;
    // go "in parallel" on vertices of p and q
    // and match q to a locally closest p neighbour
    // Q approximable if hausdorff distance
    // (maximum of these min distances)
    // is in r+e range
    for(++qit; qit != oriented_Q.vertices_end(); ++qit)
    {
        pit = p_curr;

        QPoint_2 pnt = *qit;
        PPoint_2 q_pnt = PPoint_2(
                         type_cast<PDist_2>(pnt.x()),
                         type_cast<PDist_2>(pnt.y()));

        PPoint_2 p_pnt = pit->source();
//        LOG_DEBUG << "q: " << q_pnt << std::endl;
//        LOG_DEBUG << "p: " << p_pnt << std::endl;
        min_dist = CGAL::squared_distance(q_pnt, p_pnt);
        PCurve_iterator min_pit = pit;
//        LOG_DEBUG << "dist: " << min_dist << std::endl;
        
        ++pit;
        if(pit == P.curves_end()) pit = p_start;

        do
        {
            PPoint_2 p_pnt = pit->source();
//            LOG_DEBUG << "p: " << p_pnt << std::endl;

            PDist_2 dist = CGAL::squared_distance(q_pnt, p_pnt);
//            LOG_DEBUG << "dist: " << dist << std::endl;
            if (dist > min_dist)
            {
//              LOG_DEBUG << ">" << min_dist << std::endl;
              break;
            }
            else
            {
//              LOG_DEBUG << "<=" << min_dist << std::endl;
              min_dist = dist;
              min_pit = pit;
            }

            ++pit;
            if(pit == P.curves_end()) pit = p_start;            
        }
        while(pit != p_curr);

        curr_segment = Exact_segment_2(q_pnt, min_pit->source());
        in_range = (min_dist <= in_dist);

        if(!in_range)
        {
            LOG_DEBUG << "no app: did not find a match for " << *qit << std::endl;
            o_no_data.push_back(curr_segment);
        }
        else
        {
            LOG_DEBUG << "app: found a match for " << *qit << std::endl;
            o_yes_data.push_back(curr_segment);
        }

        // continue matching from current P vertex
        p_curr = min_pit;

        // update hausdorff distance (to max of observed minimums)
        if(min_dist > max_dist)
          max_dist = min_dist;
    }

    LOG_DEBUG << "max_dist: " << max_dist << std::endl;
    approximability = (max_dist <= in_dist);

    return approximability;
}


// eps-offset delta-approximation of input polygon (Q)
template <typename CPolygon_2>
void Offset_decomposition_2<CPolygon_2>::compute_eps_app(const bool& is_inner_app)
{
  Kgon_sum_polygon_2& eps_app = (is_inner_app? m_inner_eps_app: m_outer_eps_app);
  eps_app.clear();

  Polygon_statistics& eps_app_stats = (is_inner_app? m_statistics.inner_eps_app_stats: m_statistics.outer_eps_app_stats);

  //////////////////////////////
  // offset(eps) inner/outer
//  Input_rational delta_hat1 = delta / epsilon;    // delta / epsilon;
//  Input_rational delta1 = delta;   // epsilon * ratio
//  LOG_DEBUG << "delta1 " << delta1 << " (" << delta_hat1 << ")" << std::endl;

  m_ms_offset = m_epsilon;
  m_ms_epsilon = m_delta;
  m_ms_inner = is_inner_app;

  m_circle_app.compute_kgon();
//  Polygon_2 kgon_eps = circ_app.kgon();

  m_ms_polygon = m_polygon;
  m_ms_constructor.compute_kgon_sum();

  eps_app = m_ms_constructor.kgon_sum();

  eps_app_stats.time = m_circle_app.statistics().kgon_stats.time + m_ms_constructor.statistics().kgon_sum_stats.time;

  //if (verbose())
  {
    OUT_DEBUG << "Approximability 1st step computed in " << eps_app_stats.time << " seconds." << std::endl;
  }
}

// r-inset delta-approximation of eps-offset of input polygon (Q)
template <typename CPolygon_2>
void Offset_decomposition_2<CPolygon_2>::compute_core_app(const bool& is_inner_app)
{
  Kgon_sum_polygon_2& eps_app = (is_inner_app? m_inner_eps_app: m_outer_eps_app);
  Polygon_2& core_app = (is_inner_app? m_inner_core_app: m_outer_core_app);
  core_app.clear();

  Polygon_statistics& core_app_stats = (is_inner_app? m_statistics.inner_core_app_stats: m_statistics.outer_core_app_stats);
  
  Polygon_2 eps_app_outer_boundary = eps_app.outer_boundary();
  if(eps_app.number_of_holes() > 0)
  {
    LIM_DEBUG << "!!! offset eps holes not handled yet !!!" << std::endl;
  }

  //////////////////////////////
  // inset(r) outer/inner
//  Input_rational delta2 = delta1;  // epsilon * ratio
//  Input_rational delta_hat2 = delta2 / radius; // eps_hat * ratio
//  LOG_DEBUG << "delta2 " << delta2 << " (" << delta_hat2 << "), kgon r size " << kgon_r.size() << std::endl;

  m_ms_offset = m_offset;
  m_ms_epsilon = m_delta;
  m_ms_inner = !is_inner_app;

  m_circle_app.compute_kgon();
//  Polygon_2 kgon_r = circ_app.kgon();

  m_ms_polygon = eps_app_outer_boundary;
  m_ms_constructor.compute_kgon_diff();

  core_app_stats.time = m_circle_app.statistics().kgon_stats.time + m_ms_constructor.statistics().kgon_diff_stats.time;

  //if (verbose())
  {
    OUT_DEBUG << "Approximability 2nd step computed in " << core_app_stats.time << " seconds." << std::endl;
  }
  
  Polygon_with_hole_set_2 inset_r = m_ms_constructor.kgon_diff();
  core_app = m_ms_constructor.kgon_diff_first_boundary();

  CGAL::set_pretty_mode(std::clog);
  LOG_DEBUG << "created the polygon P:" << std::endl;
  LOG_DEBUG << core_app << std::endl;
  LOG_DEBUG << std::endl;

}

// r-offset delta-approximation of core approximation
template <typename CPolygon_2>
void Offset_decomposition_2<CPolygon_2>::compute_core_msr(const bool& is_inner_app)
{
  Polygon_2& core_app = (is_inner_app? m_inner_core_app: m_outer_core_app);
  Kgon_sum_polygon_2& core_msr = (is_inner_app? m_inner_core_msr: m_outer_core_msr);
  core_msr.clear();

  Polygon_statistics& core_msr_stats = (is_inner_app? m_statistics.inner_core_msr_stats: m_statistics.outer_core_msr_stats);

  //////////////////////////////
  // offset(r) inner/outer
  m_ms_offset = m_offset;
  m_ms_epsilon = m_delta;
  m_ms_inner = is_inner_app;

  m_circle_app.compute_kgon();
//  Polygon_2 kgon_r = m_circle_app.kgon();

  m_ms_polygon = core_app;
  m_ms_constructor.compute_kgon_sum();

  core_msr = m_ms_constructor.kgon_sum();

  core_msr_stats.time = m_circle_app.statistics().kgon_stats.time + m_ms_constructor.statistics().kgon_sum_stats.time;
  OUT_DEBUG << "Core r-offset approximation computed in " << core_msr_stats.time << " seconds." << std::endl;

}

// r+eps-offset delta-approximation of core approximation
template <typename CPolygon_2>
void Offset_decomposition_2<CPolygon_2>::compute_core_ms(const bool& is_inner_app)
{
  Polygon_2& core_app = (is_inner_app? m_inner_core_app: m_outer_core_app);
  Kgon_sum_polygon_2& core_ms = (is_inner_app? m_inner_core_ms: m_outer_core_ms);
  core_ms.clear();

  Polygon_statistics& core_ms_stats = (is_inner_app? m_statistics.inner_core_ms_stats: m_statistics.outer_core_ms_stats);

  //////////////////////////////
  // offset(r+eps) inner/outer

//  Input_rational delta3 = delta2;  // epsilon * ratio
//  Input_rational delta_hat3 = delta3 / (radius + epsilon); // eps_hat * ratio / (1 + eps_hat)
// LOG_DEBUG << "delta3 " << delta3 << " (" << delta_hat3 << "), kgon r+eps size " << kgon_r_plus_eps.size() << std::endl;

  m_ms_offset = m_offset + m_epsilon;
  m_ms_epsilon = m_delta;
  m_ms_inner = is_inner_app;

  m_circle_app.compute_kgon();
//  Polygon_2 kgon_r_plus_eps = m_circle_app.kgon();

  m_ms_polygon = core_app;
  m_ms_constructor.compute_kgon_sum();

  Kgon_sum_polygon_2 offset_r_plus_eps = m_ms_constructor.kgon_sum();
//  Polygon_2 offset_r_plus_eps_outer_boundary = offset_r_plus_eps.outer_boundary();

  core_ms = offset_r_plus_eps;

  core_ms_stats.time = m_circle_app.statistics().kgon_stats.time + m_ms_constructor.statistics().kgon_sum_stats.time;
  OUT_DEBUG << "Approximability 3rd step computed in " << core_ms_stats.time << " seconds." << std::endl;
  
  // TODO: in case of holes in inset_r we need to compute their inset
  // and cut them out here before the difference
    
}

// diff of input polygon (Q) and core approximation offset
template <typename CPolygon_2>
void Offset_decomposition_2<CPolygon_2>::compute_core_diff(const bool& is_inner_app)
{
  Kgon_sum_polygon_2& core_ms = (is_inner_app? m_inner_core_ms: m_outer_core_ms);
  Kgon_sum_polygon_2& core_diff = (is_inner_app? m_inner_core_diff: m_outer_core_diff);
  Polygon_statistics& diff_stats = (is_inner_app? m_statistics.inner_diff_stats: m_statistics.outer_diff_stats);
  core_diff.clear();
  
  //////////////////////////////
  // difference
  Kgon_sum_polygon_list_2 diff_list;
  typename Kgon_sum_polygon_list_2::const_iterator pit;
  
  CGAL::Timer                  timer;

  timer.start();
  CGAL::difference (m_polygon, core_ms, std::back_inserter(diff_list));
  timer.stop();

  diff_stats.time = timer.time();
  diff_stats.size = diff_list.size();

  OUT_DEBUG << "The difference (Q - Q'): " << diff_stats.size;
  OUT_DEBUG << " computed in " << diff_stats.time << " seconds." << std::endl;
  
  CGAL::set_pretty_mode(std::clog);
  for (pit = diff_list.begin(); pit != diff_list.end(); ++pit) {
    LOG_DEBUG << "--> " << std::endl;
    LOG_DEBUG << *pit << std::endl;
  }
  LOG_DEBUG << std::endl;

  if(!diff_list.empty())
  {
    core_diff = *(diff_list.begin());
  }
}


/*!
  * Construct approximate Q_eps = offset(eps), P_tilde = inset(r), Q' = offset(r+eps)
  * and check if Q is completely inside Q'
  * delta is the approximation quality of these constructions
  * core_app (P_tilde) is the constructed core approximation
  * inner core approximation: inner, outer, inner - if yes then also yes for exact
  * outer core approximation: outer, inner, outer - if no then also no for exact
 */
template <typename CPolygon_2>
bool Offset_decomposition_2<CPolygon_2>::decide_approximability_app(const bool& is_inner_app)
{
  LOG_DEBUG << "decide_approximability_app " << (is_inner_app? "inner":"outer" ) << std::endl;

  Kgon_sum_polygon_2& eps_app = (is_inner_app? m_inner_eps_app: m_outer_eps_app);
  Polygon_2& core_app = (is_inner_app? m_inner_core_app: m_outer_core_app);
  Kgon_sum_polygon_2& core_ms = (is_inner_app? m_inner_core_ms: m_outer_core_ms);
  Kgon_sum_polygon_2& core_diff = (is_inner_app? m_inner_core_diff: m_outer_core_diff);

  Polygon_statistics& approximability_stats = (is_inner_app? 
    m_statistics.inner_approximability_stats: m_statistics.outer_approximability_stats);
  const Polygon_statistics& eps_app_stats = (is_inner_app?
    m_statistics.inner_eps_app_stats: m_statistics.outer_eps_app_stats);
  const Polygon_statistics& core_app_stats = (is_inner_app? 
    m_statistics.inner_core_app_stats: m_statistics.outer_core_app_stats);
  const Polygon_statistics& core_ms_stats = (is_inner_app? 
    m_statistics.inner_core_ms_stats: m_statistics.outer_core_ms_stats);
  const Polygon_statistics& diff_stats = (is_inner_app? 
    m_statistics.inner_diff_stats: m_statistics.outer_diff_stats);

  LOG_DEBUG << "decide_approximability_app for offset " << m_offset
    << " and eps " << m_epsilon
    << " with delta " << m_delta << std::endl;
  
  bool approximable = !is_inner_app;

  LOG_DEBUG << "compute_eps_app" << std::endl;
  if(eps_app.is_unbounded()) // is empty
    compute_eps_app(is_inner_app);
  
  LOG_DEBUG << "compute_core_app" << std::endl;
  if(core_app.is_empty())
    compute_core_app(is_inner_app);

  LOG_DEBUG << "compute_core_ms" << std::endl;  
  if(core_ms.is_unbounded()) // is empty
    compute_core_ms(is_inner_app);

  LOG_DEBUG << "compute_core_diff" << std::endl;  
  if(core_diff.is_unbounded()) // is empty
    compute_core_diff(is_inner_app);

  approximable = (core_diff.is_unbounded()); // is empty

  approximability_stats.time = eps_app_stats.time + core_app_stats.time + core_ms_stats.time + diff_stats.time;

  return approximable;
}



} // namespace CGAL


#endif // Offset_decomposition_2_h
