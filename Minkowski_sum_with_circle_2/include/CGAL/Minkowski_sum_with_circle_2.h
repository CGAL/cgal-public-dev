// Minkowski_sum_with_circle_2.h

#ifndef Minkowski_sum_with_circle_2_h
#define Minkowski_sum_with_circle_2_h

#include "colors.h"
#include <CGAL/Offset_statistics_2.h>
#include <CGAL/Minkowski_sum_construction_2.h>
#include <CGAL/Offset_construction_2.h>
#include <CGAL/Offset_decomposition_2.h>
#include <CGAL/Offset_search_2.h>
/*
#include <CGAL/approximated_offset_2.h>
#include <CGAL/offset_polygon_2.h>

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/Timer.h>
#include <iostream>
*/
namespace CGAL {


/*! \classf
 * A class serving as data holder for the GUI application, 
 * provides "on demand" update to the data, that is change 
 * in the input parameters is recorded as "needs_update" 
 * satus for the influenced data. 
 * The data are various constructions of polygon approximate 
 * and exect offsets and their decompositions performed with 
 * various comparisons. 
 */

template <typename CPolygon_2>
class Minkowski_sum_with_circle_2
{
public:
  typedef CPolygon_2 Polygon_2;

  // actual data constructors
  typedef typename Offset_statistics_2<CPolygon_2>::Types Types;
  typedef Minkowski_sum_construction_2<CPolygon_2> Min_sum_construction_2;
  typedef typename Min_sum_construction_2::Circle_app_2 Circle_app_2;
  typedef Offset_construction_2<CPolygon_2> Construction_2;
  typedef Offset_decomposition_2<CPolygon_2> Decomposition_2;
  typedef Offset_search_2<CPolygon_2> Search_2;

  // circle approximation types
  typedef typename Circle_app_2::Kgon_type_2 Kgon_type_2;
  typedef typename Circle_app_2::Polygon_size Polygon_size;

  // minkowski sum types
  typedef typename Types::Input_rational  Input_rational;
  typedef typename Min_sum_construction_2::Kgon_sum_polygon_2 Kgon_sum_polygon_2;
  typedef typename Min_sum_construction_2::Circle_list_2 Circle_list_2;
  
  // approximate (rational) offset types
  typedef typename Types::Approximate_polygon_2 Approximate_polygon_2;
  typedef typename Types::Approximate_offset_polygon_2  Approximate_offset_polygon_2;

  typedef typename Types::Approximate_polygon_list_2 Approximate_polygon_list_2;

  // exact offset types
  typedef typename Types::Rat_polygon_2 Rat_polygon_2;
  typedef typename Types::Rat_polygon_with_holes_2  Rat_polygon_with_holes_2;

  typedef typename Types::Conic_traits_2  Conic_traits_2;
  typedef typename Types::Exact_traits_2 Exact_traits_2;
  typedef typename Types::Exact_polygon_2 Exact_polygon_2;
  typedef typename Types::Exact_offset_polygon_2  Exact_offset_polygon_2;

  typedef typename Types::Exact_polygon_list_2  Exact_polygon_list_2;
  typedef typename Types::Exact_offset_polygon_list_2 Exact_offset_polygon_list_2;

  typedef typename Construction_2::Straight_skeleton_ptr_2 Straight_skeleton_ptr_2;
  typedef typename Construction_2::Straight_skeleton_2 Straight_skeleton_2;

  typedef typename Decomposition_2::Exact_alg_kernel Exact_alg_kernel;
  typedef typename Decomposition_2::Exact_segment_list_2 Exact_segment_list_2;

  struct Input
  {
      std::string polygon_name;
      unsigned int polygon_size;

      std::string offset;
      std::string eps;

//      Kgon_type_2 kgon_type;
//      unsigned int kgon_size;

      Input() :
        polygon_name(), polygon_size(0),
            offset("1/1"), eps("1/10")//,
//            kgon_type(Kgon_regular), kgon_size(8) 
        {}
  };

  struct Statistics
  {
      // input statistics
      Input input;
  };

  /// available data types 
  enum Data_type { 
    Data_polygon = 0, 
    // forward
    Data_kgon,
    Data_kgon_sum,
    Data_kgon_induced_circles,
    Data_kgon_offset_polygon,
    Data_exact_offset_polygon,
    Data_approximate_offset_polygon, 
    // forward & backward
    Data_inner_skeleton,
    Data_outer_skeleton,
    Data_inner_core,
    Data_outer_core,
    Data_approximability,
    Data_approximate_inner_core,
    Data_approximate_outer_core,
    Data_approximate_inner_kgon_sum,
    Data_approximate_outer_kgon_sum,
    Data_approximate_inner_decision,
    Data_approximate_outer_decision,   
    Data_approximate_approximability,
//    Data_search_epsilon,
//    Data_search_radius,
    Data_type_size
  };

    Minkowski_sum_with_circle_2();

    // input polygon and offset parameters
    const Polygon_2& polygon() const { return m_polygon; }
    void polygon(const Polygon_2& polygon) { 
      m_polygon = polygon;
      
      // make sure polygon has correct orientation
      if(!m_polygon.is_empty()) 
      {
        // keep directions according to polygon orientation
        if(!m_polygon.is_counterclockwise_oriented())
        {
          LOG_DEBUG << "clockwise->counterclockwise polygon" << std::endl;
          m_polygon.reverse_orientation();
        }
      }

      CGAL::set_pretty_mode(std::clog);
      LOG_DEBUG << "changed the input polygon:" << std::endl;
      LOG_DEBUG << m_polygon << std::endl;
      LOG_DEBUG << std::endl;
  
      polygon_changed(); 
    }

    const Input_rational& offset() const { return m_offset; }
    void offset(const Input_rational& offset) { 
      if(m_offset != offset) offset_changed();
      m_offset = offset; 
    }

    const Input_rational& epsilon() const { return m_epsilon; }
    void epsilon(const Input_rational& epsilon, const bool& update_size = true)
    {
      if(m_epsilon != epsilon)
      {
        m_epsilon = epsilon;
        epsilon_changed();
      }

      //Input_rational eps_from_size = m_epsilon;
      //eps_from_kgon_size(m_kgon_size, m_offset, eps_from_size);

      //bool should_update_kgon_size = (eps_from_size > m_epsilon);
      if(update_size)
      {
          bool changed = m_circle_app.update_kgon_size();
          if(changed) kgon_size_changed();
      }
    }

    const Input_rational& delta() const { return m_delta; }
    void delta(const Input_rational& delta) {
      if(m_delta != delta) delta_changed();
      m_delta = delta;
    }
    
    // circle approximation by kgon and its parameters (offset is circle radius)
    const Circle_app_2& circle_app() const { return m_circle_app; }

    const Polygon_size& kgon_size() const { return m_circle_app.kgon_size(); }
    const Kgon_type_2& kgon_type() const { return m_circle_app.kgon_type(); }
    void kgon_size(const Polygon_size& i_kgon_size, const bool& update_eps = false) 
    { 
      if(kgon_size() != i_kgon_size) kgon_size_changed();
      m_circle_app.kgon_size(i_kgon_size); 

      if(update_eps)
      {
          bool changed = m_circle_app.update_epsilon(m_epsilon);
          if(changed) epsilon_changed();
      }
    }
    void kgon_type(const Kgon_type_2& i_kgon_type) { 
      if(kgon_type() != i_kgon_type) kgon_type_changed();
      m_circle_app.kgon_type(i_kgon_type); 
    }

    // circle approximation result
    const Polygon_2& kgon() const { return m_circle_app.kgon(); }

    // minkowski sum constructions
    const Kgon_sum_polygon_2& kgon_sum() const 
    { 
      return m_ms_constructor.kgon_sum(); 
    }
    const Polygon_2& kgon_sum_outer_boundary() const 
    { 
      return m_ms_constructor.kgon_sum_outer_boundary(); 
    }
    const Circle_list_2& kgon_circles() const 
    { 
      return m_ms_constructor.kgon_circles(); 
    }
    const Circle_list_2& kgon_skipped_circles() const 
    { 
      return m_ms_constructor.kgon_skipped_circles(); 
    }
    const Circle_list_2& kgon_intersecting_circles() const
    {
      return m_ms_constructor.kgon_intersecting_circles();
    }
    const Approximate_offset_polygon_2& kgon_offset_polygon() const
    { 
      return m_ms_constructor.kgon_offset_polygon(); 
    }
    const Approximate_polygon_2& kgon_offset_outer_boundary() const 
    { 
      return m_ms_constructor.kgon_offset_outer_boundary(); 
    }

    // offset constructions
    const Exact_offset_polygon_2& exact_offset_polygon() const
    { 
      return m_polygon_offset_constructor.exact_offset_polygon(); 
    }
    const Exact_polygon_2& exact_outer_boundary() const 
    { 
      return m_polygon_offset_constructor.exact_outer_boundary();
    }

    const Approximate_offset_polygon_2& approximate_offset_polygon() const 
    { 
      return m_polygon_offset_constructor.approximate_offset_polygon();
    }
    const Approximate_polygon_2& approximate_outer_boundary() const 
    { 
      return m_polygon_offset_constructor.approximate_outer_boundary();
    }

    const Straight_skeleton_ptr_2& inner_skeleton_ptr(const bool& is_fwd) const 
    {
      const Construction_2& constructor = (is_fwd? m_kgon_sum_offset_constructor: m_polygon_offset_constructor);
      return constructor.inner_skeleton_ptr();
    }

    const Straight_skeleton_ptr_2& outer_skeleton_ptr(const bool& is_fwd) const 
    {
      const Construction_2& constructor = (is_fwd? m_kgon_sum_offset_constructor: m_polygon_offset_constructor);
      return constructor.outer_skeleton_ptr();
    }

    // decomposition constructions
    const Exact_polygon_2& outer_core(const bool& is_fwd) const
    {
      const Decomposition_2& decomposer = (is_fwd? m_kgon_sum_decomposer: m_polygon_decomposer);
      return decomposer.outer_core();
    }
    const Exact_polygon_2& inner_core(const bool& is_fwd) const
    {
      const Decomposition_2& decomposer = (is_fwd? m_kgon_sum_decomposer: m_polygon_decomposer);
      return decomposer.inner_core();
    }
    
    const bool& approximability(const bool& is_fwd) const
    {
      const Decomposition_2& decomposer = (is_fwd? m_kgon_sum_decomposer: m_polygon_decomposer);
      return decomposer.approximability();
    }
    const bool& approximability() const
    {
      return approximability(m_forward_construction);
    }
    const Exact_segment_list_2& approximable_data(const bool& is_fwd) const
    {
      const Decomposition_2& decomposer = (is_fwd? m_kgon_sum_decomposer: m_polygon_decomposer);
      return decomposer.approximable_data();
    }
    const Exact_segment_list_2& non_approximable_data(const bool& is_fwd) const
    {
      const Decomposition_2& decomposer = (is_fwd? m_kgon_sum_decomposer: m_polygon_decomposer);
      return decomposer.non_approximable_data();
    }    

    const Polygon_2& outer_core_app(const bool& is_fwd) const
    {
      const Decomposition_2& decomposer = (is_fwd? m_kgon_sum_decomposer: m_polygon_decomposer);
      return decomposer.outer_core_app();
    }
    const Polygon_2& inner_core_app(const bool& is_fwd) const
    {
      const Decomposition_2& decomposer = (is_fwd? m_kgon_sum_decomposer: m_polygon_decomposer);
      return decomposer.inner_core_app();
    }
    const Kgon_sum_polygon_2& outer_core_ms(const bool& is_fwd) const
    {
      const Decomposition_2& decomposer = (is_fwd? m_kgon_sum_decomposer: m_polygon_decomposer);
      return decomposer.outer_core_ms();
    }
    const Kgon_sum_polygon_2& inner_core_ms(const bool& is_fwd) const
    {
      const Decomposition_2& decomposer = (is_fwd? m_kgon_sum_decomposer: m_polygon_decomposer);
      return decomposer.inner_core_ms();
    }    
    const Kgon_sum_polygon_2& outer_core_diff(const bool& is_fwd) const
    {
      const Decomposition_2& decomposer = (is_fwd? m_kgon_sum_decomposer: m_polygon_decomposer);
      return decomposer.outer_core_diff();
    }
    const Kgon_sum_polygon_2& inner_core_diff(const bool& is_fwd) const
    {
      const Decomposition_2& decomposer = (is_fwd? m_kgon_sum_decomposer: m_polygon_decomposer);
      return decomposer.inner_core_diff();
    }    
    const bool& inner_approximability(const bool& is_fwd) const
    {
      const Decomposition_2& decomposer = (is_fwd? m_kgon_sum_decomposer: m_polygon_decomposer);
      return decomposer.inner_approximability();
    }
    const bool& outer_approximability(const bool& is_fwd) const
    {
      const Decomposition_2& decomposer = (is_fwd? m_kgon_sum_decomposer: m_polygon_decomposer);
      return decomposer.outer_approximability();
    }

    const bool& approximability_app(const bool& is_fwd) const
    {
      const Decomposition_2& decomposer = (is_fwd? m_kgon_sum_decomposer: m_polygon_decomposer);
      return decomposer.approximability_app();
    }

    const QColor& approximability_indicator() const
    {
      return m_app_indicator;
    }

    const Exact_offset_polygon_2& core_exact_offset_polygon(const bool& is_fwd) const
    {
      const Decomposition_2& decomposer = (is_fwd? m_kgon_sum_decomposer: m_polygon_decomposer);
      return decomposer.exact_offset_polygon();
    }
    const Exact_polygon_2& core_exact_outer_boundary(const bool& is_fwd) const
    {
      const Decomposition_2& decomposer = (is_fwd? m_kgon_sum_decomposer: m_polygon_decomposer);
      return decomposer.exact_outer_boundary();
    }

    // search constructions

    
   
    // update on-demand functions
    void needs_update();
    void update();

    void needs_update(const Data_type& i_data);
    bool update(const Data_type& i_data); // returns true if model is changed

    void polygon_changed();
    void offset_changed();
    void epsilon_changed();
    void delta_changed();
    void kgon_type_changed();
    void kgon_size_changed();

    void clear();

    /// construction type parameter (forward - construction, backward - decomposition)
    const bool& forward_construction() const { return m_forward_construction; }
    void forward_construction(const bool& status) { m_forward_construction = status; }

    /// statistics 
    void get_times(double& exact, double& approximate, double& kgon, unsigned int& kgon_size, 
      double& kgon_sum, double& arcs, double& intersects);

private:

    const QColor& get_approximability_color()
    {    
      // color according to results: yes - green, don't know - yellow, no - red
     if(inner_approximability(m_forward_construction)) // is "yes" - approximable
      {
        return polygon_color.green;
      }
      else
      if(!outer_approximability(m_forward_construction)) // is "no" - not approximable
      {
        return polygon_color.red;
      }

      return polygon_color.yellow;
    }

    // input polygon and offset and decomposition parameters
    Polygon_2 m_polygon;
    Input_rational m_offset;
    Input_rational m_epsilon;
    Input_rational m_delta;

    // search offset parameters
    Input_rational m_search_offset;
    Input_rational m_search_epsilon;
    Input_rational m_search_delta;

    // circle approximation by kgon
    Circle_app_2 m_circle_app;
    // offset via minkowski sum
    Min_sum_construction_2 m_ms_constructor;
    // standard offsets of the input polygon
    Construction_2 m_polygon_offset_constructor;
    // standard offsets of the minkowski sum polygon
    Construction_2 m_kgon_sum_offset_constructor;
    // decomposition of the input polygon
    Decomposition_2 m_polygon_decomposer;
    // decomposition of the minkowski sum polygon
    Decomposition_2 m_kgon_sum_decomposer;

    // decomposition of the input polygon
    Search_2 m_polygon_search;
    // decomposition of the minkowski sum polygon
    Search_2 m_kgon_sum_search;

    Statistics m_statistics;
    bool m_needs_update[Data_type_size];
    bool m_forward_construction;
    bool m_is_inner;
    QColor m_app_indicator;
};


template <typename CPolygon_2>
Minkowski_sum_with_circle_2<CPolygon_2>::Minkowski_sum_with_circle_2() :
m_forward_construction(true),
m_is_inner(true),
m_app_indicator(polygon_color.yellow),
m_offset(1, 1), // 1, 1
m_epsilon(1, 10), // 1, 10
m_delta(m_epsilon/4), // 1, 40
m_circle_app(m_polygon, m_offset, m_epsilon, m_is_inner),
m_ms_constructor(m_polygon, m_offset, m_epsilon, m_circle_app),
m_polygon_offset_constructor(m_polygon, m_offset, m_epsilon),
m_kgon_sum_offset_constructor(kgon_sum_outer_boundary(), m_offset, m_epsilon),
m_polygon_decomposer(m_polygon, m_offset, m_epsilon, m_delta),
m_kgon_sum_decomposer(kgon_sum_outer_boundary(), m_offset, m_epsilon, m_delta),
m_search_offset(m_offset), // 1, 1
m_search_epsilon(m_epsilon), // 1, 10
m_search_delta(m_delta), // 1, 40
m_polygon_search(m_polygon, m_search_offset, m_search_epsilon, m_search_delta),
m_kgon_sum_search(kgon_sum_outer_boundary(), m_search_offset, m_search_epsilon, m_search_delta)

{
  needs_update();
}

template <typename CPolygon_2>
void Minkowski_sum_with_circle_2<CPolygon_2>::get_times(double& exact, 
                                                   double& approximate, 
                                                   double& kgon, 
                                                   unsigned int& kgon_size, 
                                                   double& kgon_sum, 
                                                   double& arcs,
                                                   double& intersects)
{
  exact = m_polygon_offset_constructor.statistics().exact_stats.time;
  approximate = m_polygon_offset_constructor.statistics().approximate_stats.time;
  kgon = m_circle_app.statistics().kgon_stats.time;
  kgon_size = m_circle_app.statistics().kgon_stats.size;
  kgon_sum = m_ms_constructor.statistics().kgon_sum_stats.time;
  arcs = m_ms_constructor.statistics().kgon_circles_time;
  intersects = m_ms_constructor.statistics().kgon_offset_stats.time;
}


template <typename CPolygon_2>
void Minkowski_sum_with_circle_2<CPolygon_2>::update()
{
  // update kgon
  update(Data_kgon);

  // update sum
  update(Data_kgon_sum);
  // reconstruct arcs
  update(Data_kgon_induced_circles);
  update(Data_kgon_offset_polygon);

  update(Data_exact_offset_polygon);
  update(Data_approximate_offset_polygon);

  update(Data_outer_skeleton);
  update(Data_inner_skeleton);

  update(Data_outer_core);
  update(Data_inner_core);

  update(Data_approximability);

  update(Data_approximate_inner_core);
  update(Data_approximate_outer_core);
  update(Data_approximate_inner_kgon_sum);
  update(Data_approximate_outer_kgon_sum);
  update(Data_approximate_inner_decision);
  update(Data_approximate_outer_decision);

  update(Data_approximate_approximability);

}

template <typename CPolygon_2>
void Minkowski_sum_with_circle_2<CPolygon_2>::needs_update()
{
  for(int data_type = Data_polygon; data_type != Data_type_size; ++data_type)
  {
    needs_update(static_cast<Data_type>(data_type));
  }
}

template <typename CPolygon_2>
void Minkowski_sum_with_circle_2<CPolygon_2>::needs_update(const Data_type& i_data)
{
  m_needs_update[i_data] = true;
}

template <typename CPolygon_2>
void Minkowski_sum_with_circle_2<CPolygon_2>::polygon_changed()
{
  needs_update(Data_polygon);
  //if(kgon_type == Kgon_dependent) - other types too, due to slopes
  needs_update(Data_kgon);
  needs_update(Data_kgon_sum);
  needs_update(Data_kgon_induced_circles);
  needs_update(Data_kgon_offset_polygon);

  needs_update(Data_exact_offset_polygon);
  needs_update(Data_approximate_offset_polygon);

  needs_update(Data_inner_skeleton);
  needs_update(Data_outer_skeleton);

  needs_update(Data_inner_core);
  needs_update(Data_outer_core);

  needs_update(Data_approximability);

  needs_update(Data_approximate_inner_core);
  needs_update(Data_approximate_outer_core);
  needs_update(Data_approximate_inner_kgon_sum);
  needs_update(Data_approximate_outer_kgon_sum);
  needs_update(Data_approximate_inner_decision);
  needs_update(Data_approximate_outer_decision);

  needs_update(Data_approximate_approximability);
  
}

template <typename CPolygon_2>
void Minkowski_sum_with_circle_2<CPolygon_2>::offset_changed()
{
  needs_update(Data_kgon);
  needs_update(Data_kgon_sum);
  needs_update(Data_kgon_induced_circles);
  needs_update(Data_kgon_offset_polygon);

  needs_update(Data_exact_offset_polygon);
  needs_update(Data_approximate_offset_polygon);

  if(m_forward_construction)
  {
    needs_update(Data_inner_skeleton);
  }
  needs_update(Data_outer_skeleton); // boundary depends on offset
    
  needs_update(Data_inner_core);
  needs_update(Data_outer_core);

  needs_update(Data_approximability);

  needs_update(Data_approximate_inner_core);
  needs_update(Data_approximate_outer_core);
  needs_update(Data_approximate_inner_kgon_sum);
  needs_update(Data_approximate_outer_kgon_sum);
  needs_update(Data_approximate_inner_decision);
  needs_update(Data_approximate_outer_decision);

  needs_update(Data_approximate_approximability);
  
}

template <typename CPolygon_2>
void Minkowski_sum_with_circle_2<CPolygon_2>::epsilon_changed()
{
  m_delta = m_epsilon/4;
  delta_changed();
  
  needs_update(Data_approximate_offset_polygon);

  // epsilon can change kgon and all that depends on it
  // TODO: tune it more presicely to avoid unnecessary updates
  needs_update(Data_kgon);
  needs_update(Data_kgon_sum);
  needs_update(Data_kgon_induced_circles);
  needs_update(Data_kgon_offset_polygon);

  needs_update(Data_inner_skeleton);
  needs_update(Data_outer_skeleton);
  
  needs_update(Data_inner_core);
  needs_update(Data_outer_core);

  needs_update(Data_approximability);

  needs_update(Data_approximate_inner_core);
  needs_update(Data_approximate_outer_core);
  needs_update(Data_approximate_inner_kgon_sum);
  needs_update(Data_approximate_outer_kgon_sum);
  needs_update(Data_approximate_inner_decision);
  needs_update(Data_approximate_outer_decision);

  needs_update(Data_approximate_approximability);
  
}

template <typename CPolygon_2>
void Minkowski_sum_with_circle_2<CPolygon_2>::delta_changed()
{

  needs_update(Data_approximate_inner_core);
  needs_update(Data_approximate_outer_core);
  needs_update(Data_approximate_inner_kgon_sum);
  needs_update(Data_approximate_outer_kgon_sum);
  needs_update(Data_approximate_inner_decision);
  needs_update(Data_approximate_outer_decision);

  needs_update(Data_approximate_approximability);

}

template <typename CPolygon_2>
void Minkowski_sum_with_circle_2<CPolygon_2>::kgon_type_changed()
{
  needs_update(Data_kgon);
  needs_update(Data_kgon_sum);
  needs_update(Data_kgon_induced_circles);
  needs_update(Data_kgon_offset_polygon);

  if(m_forward_construction)
  {
    needs_update(Data_inner_skeleton);
    needs_update(Data_outer_skeleton);
    
    needs_update(Data_inner_core);
    needs_update(Data_outer_core);
    
    needs_update(Data_approximability);

    needs_update(Data_approximate_inner_core);
    needs_update(Data_approximate_outer_core);
    needs_update(Data_approximate_inner_kgon_sum);
    needs_update(Data_approximate_outer_kgon_sum);
    needs_update(Data_approximate_inner_decision);
    needs_update(Data_approximate_outer_decision);

    needs_update(Data_approximate_approximability);    
  }
}

template <typename CPolygon_2>
void Minkowski_sum_with_circle_2<CPolygon_2>::kgon_size_changed()
{
  needs_update(Data_kgon);
  needs_update(Data_kgon_sum);
  needs_update(Data_kgon_induced_circles);
  needs_update(Data_kgon_offset_polygon);

  if(m_forward_construction)
  {
    needs_update(Data_inner_skeleton);
    needs_update(Data_outer_skeleton);
    
    needs_update(Data_inner_core);
    needs_update(Data_outer_core);
    
    needs_update(Data_approximability);

    needs_update(Data_approximate_inner_core);
    needs_update(Data_approximate_outer_core);
    needs_update(Data_approximate_inner_kgon_sum);
    needs_update(Data_approximate_outer_kgon_sum);
    needs_update(Data_approximate_inner_decision);
    needs_update(Data_approximate_outer_decision);

    needs_update(Data_approximate_approximability);
  }
}


template <typename CPolygon_2>
bool Minkowski_sum_with_circle_2<CPolygon_2>::update(const Data_type& i_data)
{
  if(m_needs_update[i_data])
  {
    // for dual (fwd/bwd) modes:
    // kgon_sum used as input in forward (construction) mode
    // poligon used as input in backward (decomposition) mode
    const Data_type& input_data = m_forward_construction? Data_kgon_sum: Data_polygon;
    Decomposition_2& decomposer = (m_forward_construction? m_kgon_sum_decomposer: m_polygon_decomposer);
    Construction_2& skeleton_constructor = (m_forward_construction? m_kgon_sum_offset_constructor: m_polygon_offset_constructor);
    Search_2& searcher = (m_forward_construction? m_kgon_sum_search: m_polygon_search);
   
    switch(i_data)
    {
      case Data_polygon:
        m_circle_app.update_slopes();
        m_polygon_decomposer.update_slopes();
        break;
      case Data_kgon:
        update(Data_polygon);
        m_circle_app.compute_kgon();
        break;
      case Data_kgon_sum:
        update(Data_polygon);
        update(Data_kgon);
        m_ms_constructor.compute_kgon_sum();
        m_kgon_sum_decomposer.update_slopes();
        break;
      case Data_kgon_induced_circles:
        update(Data_kgon_sum);
        m_ms_constructor.compute_kgon_induced_arcs();
        break;
      case Data_kgon_offset_polygon:
        update(Data_kgon_induced_circles);
        m_ms_constructor.compute_kgon_offset_polygon();
        break;
      case Data_exact_offset_polygon:
        update(Data_polygon);
        m_polygon_offset_constructor.compute_exact_offset_polygon();
        break;
      case Data_approximate_offset_polygon: 
        update(Data_polygon);
        m_polygon_offset_constructor.compute_approximate_offset_polygon();
        break;
        
      case Data_inner_skeleton:
        update(input_data);
        skeleton_constructor.compute_inner_skeleton();
        break;
      case Data_outer_skeleton:
        update(input_data);
        skeleton_constructor.compute_outer_skeleton();
        break;
        
      case Data_inner_core:
        update(input_data);
        decomposer.compute_inner_core();
        break;
      case Data_outer_core:
        update(input_data);
        decomposer.compute_outer_core();      
        break;
      case Data_approximability:
        update(Data_outer_core);
        decomposer.compute_approximability();
        break;
        
      case Data_approximate_inner_core:
        update(input_data);
        decomposer.compute_inner_core_app();
        //decomposer.save_inner_core();
        break;
      case Data_approximate_outer_core:
        update(input_data);
        decomposer.compute_outer_core_app();
        break;
      case Data_approximate_inner_kgon_sum:
        update(Data_approximate_inner_core);
        decomposer.compute_inner_core_ms();
        break;
      case Data_approximate_outer_kgon_sum:
        update(Data_approximate_outer_core);
        decomposer.compute_outer_core_ms();
        break;
      case Data_approximate_inner_decision:
        update(Data_approximate_inner_kgon_sum);
        decomposer.compute_inner_core_diff();
        decomposer.decide_inner_approximability();
        break;
      case Data_approximate_outer_decision:
        update(Data_approximate_outer_kgon_sum);
        decomposer.compute_outer_core_diff();
        decomposer.decide_outer_approximability();
        break;
      case Data_approximate_approximability:
        update(Data_approximate_inner_decision);
        if(!decomposer.inner_approximability())
          // outer core data should be computed only if needed
          update(Data_approximate_outer_decision); 
        decomposer.compute_approximability_app();
        m_app_indicator = get_approximability_color();

        // TODO: incorporate search properly in UI
        m_search_offset = m_offset;
        m_search_epsilon = m_epsilon;
        m_search_delta = m_delta;
        
        //searcher.compute_critical_epsilon();
        //searcher.compute_maximal_radius();

        break;

//      case Data_core_offset:
//        update(Data_approximate_inner_core);
//        decomposer.compute_exact_offset();
//        break;

    };

    m_needs_update[i_data] = false;

    return true;
  }

  return false;
}

template <typename CPolygon_2>
void Minkowski_sum_with_circle_2<CPolygon_2>::clear()
{
  m_polygon.clear();

  m_circle_app.clear();
  m_ms_constructor.clear();
  m_polygon_offset_constructor.clear();
  m_kgon_sum_offset_constructor.clear();

  m_polygon_decomposer.clear();
  m_kgon_sum_decomposer.clear();

  needs_update();
}

} // namespace CGAL


#endif // Minkowski_sum_with_circle_2_h
