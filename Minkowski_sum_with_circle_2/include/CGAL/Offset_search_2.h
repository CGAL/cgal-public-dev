// Offset_search_2.h

#ifndef Offset_search_2_h
#define Offset_search_2_h

#include <CGAL/Offset_statistics_2.h>
#include <CGAL/Offset_decomposition_2.h>

namespace CGAL {


/*! \classf
 * A class implementing search for r and eps parameters
 * of the polygonal offset approximation 
 * via offset decomposition algorithms with govering statistics. 
 * The inputs are defined by caller application,
 * the class keeps them as _references_. The outputs are defined 
 * by the class and a caller application may keep a reference to them.
 */

template <typename CPolygon_2>
class Offset_search_2 : public Offset_statistics_2<CPolygon_2>
{
public:

  typedef Offset_statistics_2<CPolygon_2> Base;
  typedef typename Base::Types Types;
  typedef typename Types::Polygon_2 Polygon_2;
  typedef typename Types::Polygon_traits_2 Polygon_traits_2;  
  typedef typename Types::Rat_traits_2 Rat_traits_2;
  typedef typename Types::Input_rational Input_rational;
  typedef Offset_decomposition_2<CPolygon_2> Off_decomposition_2;

  typedef typename Polygon_traits_2::Point_2   Point_2;
//  typedef typename Polygon_traits_2::Segment_2   Segment_2;
//  typedef typename Polygon_traits_2::Direction_2   Direction_2;
//  typedef typename CGAL::Aff_transformation_2<Polygon_traits_2>  Transformation_2;
  
  
  //typedef typename Types::Exact_alg_kernel Exact_alg_kernel;
  //typedef typename Exact_alg_kernel::Segment_2 Exact_segment_2;
  //typedef typename std::list<Exact_segment_2> Exact_segment_list_2; 

  typedef typename Base::Offset_statistics Offset_statistics;
  typedef typename Base::Polygon_statistics Polygon_statistics;

  //typedef typename Types::Conic_traits_2 Conic_traits_2;

  //typedef typename Types::Exact_rational Exact_rational;
  //typedef typename Types::Exact_algebraic Exact_algebraic;
  //typedef typename Types::Rat_polygon_2 Rat_polygon_2;
  //typedef typename Types::Rat_polygon_with_holes_2 Rat_polygon_with_holes_2;
  //typedef typename Types::Exact_offset_polygon_2 Exact_offset_polygon_2; 
  //typedef typename Types::Exact_polygon_2 Exact_polygon_2;
  //typedef typename Types::Exact_polygon_list_2 Exact_polygon_list_2;
  //
  //typedef typename CGAL::Circle_approximation_2<CPolygon_2>    Circle_app_2;
  //typedef typename CGAL::Polygon_with_holes_2<Polygon_traits_2>  Kgon_sum_polygon_2;
  //typedef typename CGAL::General_polygon_set_2<Rat_traits_2> Polygon_with_hole_set_2;
  //typedef typename std::list<Kgon_sum_polygon_2>                 Kgon_sum_polygon_list_2;
  
  struct Statistics
  {
      // kgon sum statistics
      Polygon_statistics outer_approximability_stats;
      Polygon_statistics inner_approximability_stats;
      
      double eps_search_time;
      int eps_search_decisions;
      int eps_search_iterations;

      double r_search_time;
      int r_search_decisions;
      int r_search_iterations;

  };

  /// input parameters stored _as references_
  Offset_search_2(
    const Polygon_2& i_polygon,
    Input_rational& i_offset,
    Input_rational& i_eps,
    Input_rational& i_delta) :
    m_polygon(i_polygon),
    m_offset(i_offset),
    m_epsilon(i_eps),
    m_delta(i_delta),
    // decomposer 
    m_decomposer(i_polygon, i_offset, i_eps, i_delta),
    m_inner_decomposer(i_polygon, i_offset, i_eps, i_delta),
    m_outer_decomposer(i_polygon, i_offset, i_eps, i_delta)
  {
    m_decomposer.verbose(false);
    m_inner_decomposer.verbose(false);
    m_outer_decomposer.verbose(false);
  }

  ~Offset_search_2(void)
  {
  }

  const Statistics& statistics() const { return m_statistics; }

  void compute_critical_epsilon();
  void compute_maximal_radius();

private:

  const Polygon_2& m_polygon;
  Input_rational& m_offset;
  Input_rational& m_epsilon;
  Input_rational& m_delta;

  Off_decomposition_2 m_decomposer;
  Off_decomposition_2 m_inner_decomposer;
  Off_decomposition_2 m_outer_decomposer;

  Statistics m_statistics;
};
  

template <typename CPolygon_2>
void Offset_search_2<CPolygon_2>::compute_critical_epsilon()
{
  LOG_DEBUG << "compute_critical_epsilon for r=" << m_offset << " and delta=" << m_delta  << " (eps diff limit)" << std::endl;

  m_statistics.eps_search_time = 0;
  m_statistics.eps_search_decisions = 0;
  m_statistics.eps_search_iterations = 0;

  // we assume eps = r has an answer YES (approximable with core = input polygon)
  // and eps = 0 has an answer NO (not approximable)

  // we perform a binary search for eps in [0 ... r] with initial delta = I/8
  // if decision algorithms return YES or NO we update the known eps limits accordingly
  // if limits got closer then required delta_limit we stop and return yes limit
  
  // if both sides of the decision algorithm are UNDECIDED
  // it means that we are at most 2xdelta from the critical epsilon
  // therefore we halve the delta and rerun the decision with the same epsilon

  // if the result is YES or NO we update the limits (and keep the delta) accordingly
  // if the result is UNDECIDED again we know for sure that yes and no limits can
  // be updated to eps_curr +-4xdelta
  // however we reduce the interval even further since at least one of the
  // eps_curr +- 2xdelta limits has a definitive answer

  Input_rational delta_limit = m_delta;
  Input_rational eps_yes = m_offset;
  Input_rational eps_no = 0;

  Input_rational eighth(1, 8);
  Input_rational eps_curr = (eps_yes + eps_no)/2;
  Input_rational interval = (eps_yes - eps_no);
  Input_rational delta_ratio = eighth;
  Input_rational delta_curr = interval * delta_ratio;
  Input_rational delta_yes = delta_curr;
  
  for(; interval > delta_limit;
        interval = (eps_yes - eps_no),
        eps_curr = (eps_yes + eps_no)/2,
        delta_curr = interval * delta_ratio)
  {
    m_epsilon = eps_curr;
    m_delta = delta_curr;

    LOG_DEBUG << "eps_curr = " << eps_curr << std::endl;
    LOG_DEBUG << "interval = " << interval << std::endl;
    LOG_DEBUG << "delta_ratio = " << delta_ratio << std::endl;
    LOG_DEBUG << "delta_curr = " << delta_curr << std::endl;

    ++m_statistics.eps_search_iterations;
    
    // check if YES
    m_decomposer.clear();
    m_decomposer.decide_inner_approximability();
    bool is_yes = m_decomposer.inner_approximability();
    ++m_statistics.eps_search_decisions;
    m_statistics.eps_search_time += m_decomposer.statistics().inner_approximability_stats.time;

    if(is_yes)
    {
      OUT_DEBUG << "YES: eps = " << m_epsilon << ", delta = " << m_delta << std::endl;
      eps_yes = m_epsilon;
      delta_yes = m_delta;
      delta_ratio = eighth;      
    }
    else
    {
      // check if NO
      m_decomposer.decide_outer_approximability();
      bool is_no = !m_decomposer.outer_approximability();
      ++m_statistics.eps_search_decisions;
      m_statistics.eps_search_time += m_decomposer.statistics().outer_approximability_stats.time;

      if(is_no)
      {
        OUT_DEBUG << "NO_: eps = " << m_epsilon << ", delta = " << m_delta << std::endl;
        eps_no = m_epsilon;
        delta_ratio = eighth;
      }
      else
      {
        // is UNDECIDED
        OUT_DEBUG << "UND: eps = " << m_epsilon << ", delta = " << m_delta << std::endl;

        // halve the delta
        m_delta = delta_curr/2;
         
        // check if YES
        m_decomposer.clear();
        m_decomposer.decide_inner_approximability();
        bool is_yes = m_decomposer.inner_approximability();
        ++m_statistics.eps_search_decisions;
        m_statistics.eps_search_time += m_decomposer.statistics().inner_approximability_stats.time;

        if(is_yes)
        {
          OUT_DEBUG << "YES2: eps = " << m_epsilon << ", delta = " << m_delta << std::endl;
          eps_yes = m_epsilon;
          delta_yes = m_delta;
          delta_ratio = delta_ratio/2;          
        }
        else
        {
          // check if NO
          m_decomposer.decide_outer_approximability();
          bool is_no = !m_decomposer.outer_approximability();
          ++m_statistics.eps_search_decisions;
          m_statistics.eps_search_time += m_decomposer.statistics().outer_approximability_stats.time;

          if(is_no)
          {
            OUT_DEBUG << "NO_2: eps = " << m_epsilon << ", delta = " << m_delta << std::endl;
            eps_no = m_epsilon;
            delta_ratio = delta_ratio/2;
          }
          else
          {
            // is UNDECIDED
            OUT_DEBUG << "UND2: eps = " << m_epsilon << ", delta = " << m_delta << std::endl;

            // at least one of the +-2*delta limits should give definitive result
            // both +-4*delta have definitive result
            m_epsilon = eps_curr + 2 * m_delta;
            // check if YES
            m_inner_decomposer.clear();
            m_inner_decomposer.decide_inner_approximability();
            bool is_yes = m_inner_decomposer.inner_approximability();
            ++m_statistics.eps_search_decisions;
            m_statistics.eps_search_time += m_inner_decomposer.statistics().inner_approximability_stats.time;

            if(is_yes)
            {
              OUT_DEBUG << "YES2: eps = " << m_epsilon << ", delta = " << m_delta << std::endl;
              eps_yes = m_epsilon;
              delta_yes = m_delta;
            }
            else
            {
              m_epsilon = eps_curr + 4 * m_delta;
              OUT_DEBUG << "YES3: eps = " << m_epsilon << ", delta = " << m_delta << std::endl;
              eps_yes = m_epsilon;
              delta_yes = m_delta;
            }

            m_epsilon = eps_curr - 2 * m_delta;
            // check if NO
            m_outer_decomposer.clear();
            m_outer_decomposer.decide_outer_approximability();
            bool is_no = !m_outer_decomposer.outer_approximability();
            ++m_statistics.eps_search_decisions;
            m_statistics.eps_search_time += m_outer_decomposer.statistics().outer_approximability_stats.time;
            if(is_no)
            {
              OUT_DEBUG << "NO_2: eps = " << m_epsilon << ", delta = " << m_delta << std::endl;
              eps_no = m_epsilon;
            }
            else
            {
              m_epsilon = eps_curr - 4 * m_delta;
              OUT_DEBUG << "NO_3: eps = " << m_epsilon << ", delta = " << m_delta << std::endl;
              eps_no = m_epsilon;
            }
            
            if(!is_yes && !is_no)
            {
              m_epsilon = eps_curr;
              LIM_DEBUG << "BUG!: eps = " << m_epsilon << ", delta = " << m_delta << std::endl;
            }
          }
        } 
      }      
    }
  }

  OUT_DEBUG << "Critical epsilon is within " << interval << " from eps_yes = " << eps_yes
            << " (computed with delta = " << delta_yes << ") final delta = " << m_delta << std::endl;
  OUT_DEBUG << "Critical epsilon search has taken " << m_statistics.eps_search_time << " time for " << m_statistics.eps_search_decisions
            << " decision procedures in " << m_statistics.eps_search_iterations << " iterations." << std::endl;

}

template <typename CPolygon_2>
void Offset_search_2<CPolygon_2>::compute_maximal_radius()
{
  LOG_DEBUG << "compute_maximal_radius for eps=" << m_epsilon << " and delta=" << m_delta << " (r diff limit)" << std::endl;

  // we assume r = eps has an answer YES (approximable with core = input polygon)
  // and r = bounding_box_long_side has an answer NO (not approximable)

  // we perform a binary search for r in [eps ... bbs] with initial delta = eps/8
  // if decision algorithms return YES or NO we update the known r limits accordingly
  // if limits got closer then required Delta_limit we stop and return yes limit

  // if both sides of the decision algorithm are UNDECIDED on current r
  // it means that we are at most 2xdelta from the critical epsilon of this r
  // therefore we halve the delta and rerun the decision with the same r

  // if the result is YES or NO we update the limits (and keep the delta) accordingly
  // if the result is UNDECIDED again we know for sure that yes and no limits can
  // be updated to eps_curr +-4xdelta
  // however we reduce the interval even further since at least one of the
  // eps_curr +- 2xdelta limits has a definitive answer


  // we assume r = eps has an answer YES (approximable with core = input polygon)
  // we double r with delta = 1/4 eps until the answer is not YES
  // after that we subdivide the resulting yes/not yes interval
  // reducing the delta when necessary


  Input_rational delta_limit = m_delta;
  Input_rational r_yes = m_epsilon;
  Input_rational r_no = m_epsilon + m_delta;

  double bboxSize = 1;
  if(m_polygon.size() != 0)
  {
    CGAL::Bbox_2 bbox = m_polygon.bbox();
    double bboxXDim = (bbox.xmax() - bbox.xmin());
    double bboxYDim = (bbox.ymax() - bbox.ymin());
    bboxSize = bboxXDim > bboxYDim? bboxXDim: bboxYDim;
    LOG_DEBUG << "bbox of size " << bboxSize << std::endl;

    r_no = bboxSize;
  }


  Input_rational eighth(1, 8);
  Input_rational r_curr = (r_no + r_yes)/2;
  Input_rational interval = (r_no - r_yes);

  Input_rational delta_ratio = eighth;
  Input_rational delta_curr = m_epsilon * delta_ratio;
  Input_rational delta_yes = delta_curr;

  for(; interval > delta_limit;
        interval = (r_no - r_yes),
        r_curr = (r_no + r_yes)/2,
        delta_curr = delta_curr / 2)
  {
    m_offset = r_curr;
    m_delta = delta_curr;

    ++m_statistics.r_search_iterations;
    
    // check if YES
    m_decomposer.clear();
    m_decomposer.decide_inner_approximability();
    bool is_yes = m_decomposer.inner_approximability();
    ++m_statistics.r_search_decisions;
    m_statistics.r_search_time += m_decomposer.statistics().inner_approximability_stats.time;

    if(is_yes)
    {
      OUT_DEBUG << "YES: r = " << m_offset << ", delta = " << m_delta << std::endl;
      r_yes = m_offset;
      delta_yes = m_delta;
    }
    else
    {
      // check if NO
      m_decomposer.decide_outer_approximability();
      bool is_no = !m_decomposer.outer_approximability();
      ++m_statistics.r_search_decisions;
      m_statistics.r_search_time += m_decomposer.statistics().outer_approximability_stats.time;

      if(is_no)
      {
        OUT_DEBUG << "NO_: r = " << m_offset << ", delta = " << m_delta << std::endl;
        r_no = m_offset;
      }
      else
      {
        // is UNDECIDED
        OUT_DEBUG << "UND: r = " << m_offset << ", delta = " << m_delta << std::endl;
        //delta_ratio = delta_ratio/2;
      }
    }
  }

  OUT_DEBUG << "Maximal radius is within " << interval << " from r_yes = " << r_yes
            << " (computed with delta = " << delta_yes << ") final delta = " << m_delta << std::endl;
  OUT_DEBUG << "Maximal radius search has taken " << m_statistics.r_search_time << " time for " << m_statistics.r_search_decisions
            << " decision procedures in " << m_statistics.r_search_iterations << " iterations." << std::endl;

}

} // namespace CGAL


#endif // Offset_search_2_h