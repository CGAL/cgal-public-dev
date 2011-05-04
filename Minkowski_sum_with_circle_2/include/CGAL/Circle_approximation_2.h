// Circle_approximation_2

#ifndef CGAL_CIRCLE_APPROXIMATION_2_H
#define CGAL_CIRCLE_APPROXIMATION_2_H

#include <set>
#include <map>
#include <iterator>
#include <algorithm>

#include <CGAL/Offset_statistics_2.h>
#include <CGAL/Rational_points_on_circle_2.h>
#include <CGAL/enum.h>
#include <CGAL/Gmpq.h>
#include <CGAL/CORE_BigRat.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Fraction_traits.h>

#include <CGAL/Timer.h>
#include "output_debug.h"

namespace CGAL {

/*! \classf
 * A class implementing various circle approximations 
 * that can be used in a Minkowski Sum to create approximate
 * offset of a polygon.
 * The circle approximation is a kgon with a rational
 * coordinate vertices and unique edge slopes, that differ 
 * from the slopes of the polygon to offset. The "slopes" here
 * are directed, i.e. segments with opposite directions considered 
 * to have different slopes.
 * Currently, only "inner" circle approximations are implemented,
 * that is the ones with all vertices on the circle boundary, with
 * the resulting kgon completely inscribed in the circle.
 */

template <typename CPolygon_2>
class Circle_approximation_2 : public Offset_statistics_2<CPolygon_2>
{
public:

  typedef Offset_statistics_2<CPolygon_2> Base;
  typedef typename Base::Types Types;
  typedef typename Types::Polygon_2 Polygon_2;
  typedef typename Types::Polygon_traits_2 Polygon_traits_2;
  typedef typename Types::Input_rational Input_rational;
  typedef typename Base::Polygon_statistics Polygon_statistics;
  typedef typename CGAL::Rational_points_on_circle_2<Polygon_traits_2> Rational_points_2;


  enum Kgon_type_2 { 
    Kgon_regular = 0, 
    Kgon_random, 
    Kgon_dependent,
    Kgon_const
    };
    
  typedef typename Polygon_2::Container Container;
  typedef typename Container::size_type Polygon_size;

  typedef typename Types::Polygon_with_holes_2 Polygon_with_holes_2;

  typedef typename Polygon_traits_2::Point_2   Point_2;
  typedef typename Polygon_traits_2::Segment_2   Segment_2;
  typedef typename Polygon_traits_2::Direction_2   Direction_2;
  typedef typename CGAL::Aff_transformation_2<Polygon_traits_2>  Transformation_2;

  typedef typename Container::size_type   size_type;

  typedef typename Rational_points_2::Points_on_circle_2 Points_on_circle_2;

  // Note(!): 
  // CGAL::compare_slopes does not take direction into account
  // therefore not used
  /// operand for comparing directed segment slopes in containers
  struct Less_directed_slope_2
  {
    bool operator()(const Segment_2& s1, const Segment_2& s2) const
    {
      return s1.direction() < s2.direction();
    }
  };

  /// operand for comparing directed segment slopes in algorithms 
  struct Equal_directed_slope_2
  {
    Equal_directed_slope_2(const Segment_2& s)
    {
      dir = s.direction();
    }

    bool operator()(const Segment_2& s) const
    {
      return dir == s.direction();
    }
    Direction_2 dir;
  };

  typedef typename std::set<Segment_2, Less_directed_slope_2> Directed_slope_set;

  struct Statistics
  {
      // kgon statistics
      Polygon_statistics kgon_stats;
  };

  const Statistics& statistics() const { return m_statistics; }

  /*!
   * Constructor with an input polygon (the one that is being offseted).
   * Sorts i_polygon segment slopes according to segment _direction_, 
   * that is the same slope used in two different directions treated as
   * different slopes.
   * \param i_polygon The input polygon.
   * \param i_offset The circle radius.
   * \param i_epsilon The max allowed distance from disk of the same radius.
   * \param m_kgon The output polygon.
   * \param m_kgon_size The minimal number of kgon vertices.
   */
  Circle_approximation_2(
  const Polygon_2& i_polygon,
  const Input_rational& i_offset,
  const Input_rational& i_epsilon,
  const bool& i_is_inner) :
    m_polygon(i_polygon),
    m_radius(i_offset),
    m_epsilon(i_epsilon),
    m_is_inner(i_is_inner),
    // outer approximation simulated via the inner approximation
    // of the circle with radius = offset + epsilon
    m_scaling(m_is_inner? m_radius: m_radius + m_epsilon),
    m_rat_points(m_scaling, i_epsilon),
    m_kgon_type(Kgon_regular),
    m_kgon_size(8) // 8
  {
    update_slopes();
    //update_kgon_size();
    Rational_points_2::precompute_short_rationals();
  }

  void update_slopes();

  const Polygon_2& kgon() const { return m_kgon; }
//  void kgon(const Polygon_2& kgon) { m_kgon = kgon; }

  const Polygon_size& kgon_size() const { return m_kgon_size; }
  const Kgon_type_2& kgon_type() const { return m_kgon_type; }
  void kgon_size(const Polygon_size& kgon_size) { 
    m_kgon_size = kgon_size; 
  }
  void kgon_type(const Kgon_type_2& kgon_type) { 
    m_kgon_type = kgon_type; 
  }

  bool update_kgon_size()
  {
    Polygon_size size = m_kgon_size;
    kgon_size_from_eps(m_epsilon, m_radius, size);
    LOG_DEBUG << "update kgon size from " << m_kgon_size
              << " to " << size << std::endl;
    bool size_has_changed = (m_kgon_size != size);
    kgon_size(size);
    return size_has_changed;
  }

  bool update_epsilon(Input_rational& o_epsilon) const
  {
    Input_rational epsilon = m_epsilon;
    eps_from_kgon_size(m_kgon_size, m_radius, o_epsilon);
    LOG_DEBUG << "update epsilon from " << epsilon
              << " to " << o_epsilon << std::endl;
    return (epsilon != o_epsilon);
  }

  void clear()
  {
    LOG_DEBUG << "Circle_approximation_2::clear()" << std::endl;
    update_slopes();
    clear_kgon();
  }

  void clear_kgon()
  {
    m_kgon.clear();
  }
  
 /*!
   * Construct a circle approximation with a kgon tailored
   * for the input polygon according to type and epsilon (or size)
   * parameters.
   * \param m_polygon The input polygon.
   * \param m_epsilon The max allowed distance from disk of the same radius.
   * \param m_radius The circle radius.
   * \param m_kgon The output polygon.
   */
  void compute_kgon();

  /// nicely formatted output of segment slope according to direction
  static void log_slope(const std::string& prefix, const Segment_2& seg) 
  { 
    std::ostringstream slope_str;
 //   LOG_DEBUG << "circle__offset (" << seg << ") \t"  << seg.direction() << std::endl;

    Direction_2 dir = seg.direction();
    Input_rational x = dir.dx(), y = dir.dy();
    Input_rational slope = 0;
    if(x != 0 && y != 0)
    {
      slope = y/x;
      slope_str << slope;
 //     LOG_DEBUG << " \t" << slope << " \t" << seg.direction();
    }
    else
    if(y == 0)
    {
      slope_str << ((x > 0)? "+0": "-0");

//      LOG_DEBUG << " \t" << ((x > 0)? "+0": "-0") << " \t" << seg.direction();
    }
    else
    {
      slope_str << ((y > 0)? "+inf": "-inf");

//      LOG_DEBUG << " \t" << ((y > 0)? "+inf": "-inf") << " \t" << seg.direction();
    }

    LOG_DEBUG << prefix << " (" << slope_str.str() << ") \t"  
      << seg.direction() << " \t" << seg << std::endl;

  };


  /// does segment in question has a kgon slope (that is not a polygon slope)
  bool is_circle_slope(const Segment_2& i_segment) const
  {
    return is_unique_slope(i_segment);
  }

  void log_input_slopes() const;

  void log_circle_edges(const Polygon_with_holes_2&
                        i_polygon_with_holes) const;
  
  void log_circle_edges(const Polygon_2& i_polygon) const;

  size_type number_of_unique_slopes() {
    return m_input_slopes.size();
  }

  static void kgon_size_from_eps(const Input_rational& epsilon, 
    const Input_rational& radius, Polygon_size& kgon_size);
  static void eps_from_kgon_size(const Polygon_size& kgon_size, 
    const Input_rational& radius, Input_rational& epsilon);

 
private:

    /*!
   * Construct a circle approximation with an almost regular kgon
   * (vertices are close to regulary distributed on a circle and 
   * there are at least m_kgon_size of them).
   * Note: Rational coordinates of vertices are "naive", converted 
   * from double, that is unnecessarily long.
   * \param m_kgon The output polygon.
   * \param m_kgon_size The minimal number of kgon vertices.
   * \param m_radius The circle radius.
   */
  void regular_kgon_size();

 /*!
   * Construct a circle approximation with an almost regular kgon
   * (vertices are close to regulary distributed on a circle and 
   * ensure distance of at most m_epsilon between kgon and circle).
   * Note: Rational coordinates
   * \param m_kgon The output polygon.
   * \param m_epsilon The max dist from disk of the same radius.
   * \param m_radius The circle radius.
   */
  void regular_kgon_eps();

 /*!
   * Construct a circle approximation with a kgon tailored
   * for the input polygon (vertices are chosen in a way that
   * ensures distance of at most i_eps between approximate
   * minkowski sum and exact offset segments).
   * \param m_kgon The output polygon.
   * \param m_epsilon The max dist from disk of the same radius.
   * \param m_radius The circle radius.
   */
  void dependent_kgon();

  /*!
   * Construct a circle approximation with a random kgon
   * (randomly chosen vertices on the circle according 
   * to m_kgon_size).
   * \param m_kgon The output polygon.
   * \param m_radius The circle radius.
   * \param m_kgon_size The minimal number of kgon vertices.
   */
  void random_kgon();
  
  /*!
   * Construct a circle approximation with a predefined kgon
   * (currently it is a rombic square - useful for debugging).
   * \param m_kgon The output polygon.
   * \param m_radius The circle radius.
   */
  void const_kgon();
 
  bool has_non_unique_slopes(const Polygon_2& i_polygon,
    Directed_slope_set& o_non_unique_slopes) const;

  /// does segment in question has unique (not in an input polygon) slope
  bool is_unique_slope(const Segment_2& i_segment) const;

  /*!
   * Replace polygon edges with non-unique slopes (slopes that a polygon
   * to be offset has) by "splitting" the edge to two almost equal parts
   * (create mid-point with the rational coordinates on the circle).
   * \param i_non_unique_slopes The problematic slopes.
   * \param io_polygon The polygon that needs fixing.
   */
  void fix_non_unique_slopes_by_segment_splitting(
      const Directed_slope_set& i_non_unique_slopes, Polygon_2& io_polygon) const;
  

  static Transformation_2 get_scaling(const Input_rational& i_radius)
  {
    typedef typename CGAL::Fraction_traits<Input_rational> FracTraits;
    typename FracTraits::Numerator_type num;
    typename FracTraits::Denominator_type denom(1);
    typename FracTraits::Decompose decomposer;

    decomposer(i_radius, num, denom);
    Transformation_2 scale(CGAL::SCALING, num, denom); 
   // LOG_DEBUG << "scaling by: " << num <<" / " << denom << std::endl;
   // LOG_DEBUG << "scaling: " << scale << std::endl;
    return scale;
  }

  const Polygon_2& m_polygon;
  const Input_rational& m_radius;
  const Input_rational& m_epsilon;
  const bool& m_is_inner;
  Input_rational m_scaling;
  
  Directed_slope_set m_input_slopes;

  Polygon_2 m_kgon;
  Polygon_size m_kgon_size;
  Kgon_type_2 m_kgon_type;

  Statistics m_statistics;
  
  Rational_points_2 m_rat_points;
};


template <typename CPolygon_2>
void Circle_approximation_2<CPolygon_2>::update_slopes()
{
  LOG_DEBUG << "Polygon edges: " << m_polygon.size() << std::endl;

  if(!m_polygon.is_empty()) 
  {
    assert(m_polygon.is_counterclockwise_oriented());

    // sort and store the input polygon slopes
    m_input_slopes = Directed_slope_set(m_polygon.edges_begin(), m_polygon.edges_end());
  }
  else
  {
    m_input_slopes.clear();
  }

  LOG_DEBUG << "Unique directed slopes: " << m_input_slopes.size() << std::endl;
}

template <typename CPolygon_2>
void Circle_approximation_2<CPolygon_2>::compute_kgon()
{
  CGAL::Timer timer;
  timer.start();

  switch(m_kgon_type)
  {
  case Kgon_regular:
    {
      //regular_kgon_size(m_kgon_size);
      regular_kgon_eps();
      //kgon_size(m_kgon.size());
    }
    break;
  case Kgon_random:
    {
      random_kgon();
    }
    break;
  case Kgon_dependent:
    {
      //regular_kgon_eps();
      //kgon_size(m_kgon.size());

      dependent_kgon();
      //kgon_size(m_kgon.size());
    }
    break;
  default:
    {
      const_kgon();
    }
    break;
  };

  timer.stop();

  m_statistics.kgon_stats.time = timer.time();
  m_statistics.kgon_stats.size = kgon().size();

  //if (verbose())
  {
    OUT_DEBUG << (m_is_inner? "Inner ": "Outer ") << "Kgon of size " << m_statistics.kgon_stats.size
        << " computed in " << m_statistics.kgon_stats.time << std::endl;
  }

}

template <typename CPolygon_2>
bool Circle_approximation_2<CPolygon_2>::has_non_unique_slopes(
    const Polygon_2& i_polygon, Directed_slope_set& o_non_unique_slopes) const
{
  o_non_unique_slopes.clear();
  
  // compare the slopes to the slopes of the original polygon
  Directed_slope_set outputSlopes(i_polygon.edges_begin(), i_polygon.edges_end());
  
  //typename Directed_slope_set::iterator result_iter = o_non_unique_slopes.begin();
  set_intersection(m_input_slopes.begin(), m_input_slopes.end(),
                   outputSlopes.begin(), outputSlopes.end(),
                   typename std::insert_iterator<Directed_slope_set>(
                       o_non_unique_slopes,
                       o_non_unique_slopes.begin()
                   ),
                   Less_directed_slope_2());
  
  LOG_DEBUG << o_non_unique_slopes.size() << " non unique slopes" <<std::endl;
  return !o_non_unique_slopes.empty();
}

template <typename CPolygon_2>
bool Circle_approximation_2<CPolygon_2>::is_unique_slope(
    const Segment_2& i_segment) const
{
  typename Directed_slope_set::const_iterator slope_iter =
    find_if(m_input_slopes.begin(), m_input_slopes.end(),
            Equal_directed_slope_2(i_segment));
  
  return (slope_iter == m_input_slopes.end());
}

template <typename CPolygon_2>
void Circle_approximation_2<CPolygon_2>::
fix_non_unique_slopes_by_segment_splitting(
    const Directed_slope_set& i_non_unique_slopes, Polygon_2& io_polygon) const
{
  // if there is a non-unique slope get rid of it
  // by splitting a problematic segment to 2 segments
  // (with rational mid-point)
  std::cerr
    << "fix_non_unique_slopes_by_segment_splitting is not implemented yet"
    << std::endl;
/*
  typename Directed_slope_set::const_iterator slope_iter = i_non_unique_slopes.begin();
  typename Directed_slope_set::const_iterator slope_iter_end = i_non_unique_slopes.end();

  while(slope_iter != slope_iter_end)
  {


    ++ slope_iter;
  }
*/
}

template <typename CPolygon_2>
void Circle_approximation_2<CPolygon_2>::const_kgon()
{
  m_kgon.clear();

  m_scaling = m_is_inner? m_radius: m_radius + m_epsilon;
  
  // creating by default square polygon
//  Point_2 points[] = {
//    Point_2(0.0, m_scaling),
//    Point_2(-m_scaling, 0.0),
//    Point_2(0.0, -m_scaling),
//    Point_2(m_scaling, 0.0)
//  };

//  m_kgon = Polygon_2(points, points+4);

  // counterclockwise orientation
  m_kgon.push_back(Point_2(0.0, m_scaling));
  m_kgon.push_back(Point_2(-m_scaling, 0.0));
  m_kgon.push_back(Point_2(0.0, -m_scaling));
  m_kgon.push_back(Point_2(m_scaling, 0.0));
}


template <typename CPolygon_2>
void Circle_approximation_2<CPolygon_2>::regular_kgon_size()
{
  m_kgon.clear();

  m_scaling = m_is_inner? m_radius: m_radius + m_epsilon;
  Transformation_2 scale = get_scaling(m_scaling);

  const Input_rational pi = CGAL_PI; // not exact
  const Input_rational angle_fraction = 2 * pi / m_kgon_size;
  
  // create regular k-gon
  size_type count;
  for(count=0; count < m_kgon_size; ++count)
  {
    Input_rational angle = angle_fraction * count;
    double double_angle = CGAL::to_double(angle);
    double x = std::sin(double_angle);
    double y = std::cos(double_angle);
    
    Point_2 exactPoint =
      Rational_points_2::convert_to_exact_point_on_unit_circle(Point_2(x, y));
//    LOG_DEBUG << "point: " << exactPoint << std::endl;
//    LOG_DEBUG << "scaled: " << scale(exactPoint) << std::endl;
    m_kgon.push_back(scale(exactPoint));
  }
  
  // make all its coordinates rational, but still on a circle boundary
  
  // compare the slopes to the slopes of the original polygon
  Directed_slope_set non_unique_slopes;
  bool need_to_fix = has_non_unique_slopes(m_kgon, non_unique_slopes);
  
  LOG_DEBUG << (need_to_fix ? "" : "no ") << "need to fix" << std::endl;
  
  // if there is a non-unique slope get rid of it
  // by splitting a problematic segment to 2 segments
  // (with rational mid-point)
  if (need_to_fix) {
    fix_non_unique_slopes_by_segment_splitting(non_unique_slopes, m_kgon);
  }
  
  CGAL::set_pretty_mode(std::clog);
  LOG_DEBUG << "created the polygon Br:" << std::endl;
  LOG_DEBUG << m_kgon << std::endl;
  LOG_DEBUG << std::endl;
}

template <typename CPolygon_2>
void Circle_approximation_2<CPolygon_2>::regular_kgon_eps()
{
  m_kgon.clear();

  m_scaling = m_is_inner? m_radius: m_radius + m_epsilon;
  Transformation_2 scale = get_scaling(m_scaling);
  // compute k from eps
  
  Polygon_2 first_quadrant_polygon;
  m_rat_points.get_first_quadrant_rationals(first_quadrant_polygon);

  Polygon_2 unit_polygon;
  Rational_points_2::mirror_first_quadrant(first_quadrant_polygon, unit_polygon);

  LOG_DEBUG << "found " << unit_polygon.size() << " rational points" << std::endl;

  // compare the slopes to the slopes of the original polygon
  Directed_slope_set non_unique_slopes;
  bool need_to_fix = has_non_unique_slopes(unit_polygon, non_unique_slopes);
  
  LOG_DEBUG << (need_to_fix ? "" : "no ") << "need to fix" << std::endl;
  
  // if there is a non-unique slope get rid of it
  // by splitting a problematic segment to 2 segments
  // (with rational mid-point)
  if (need_to_fix) {
    fix_non_unique_slopes_by_segment_splitting(non_unique_slopes, unit_polygon);
  }

  typename Polygon_2::Vertex_const_iterator vit;
  for (vit = unit_polygon.vertices_begin(); vit != unit_polygon.vertices_end(); ++vit) 
  {
    m_kgon.push_back(scale(Point_2(vit->x(), vit->y())));
  }

  CGAL::set_pretty_mode(std::clog);
  LOG_DEBUG << "created the polygon Br:" << std::endl;
  if(m_kgon.size() < size_type(10000))
    LOG_DEBUG << m_kgon << std::endl;
  else
    LOG_DEBUG << "\t of huge size " << m_kgon.size() << std::endl;
  LOG_DEBUG << std::endl;
}

template <typename CPolygon_2>
void Circle_approximation_2<CPolygon_2>::dependent_kgon()
{
  m_kgon.clear();

  Transformation_2 rotate(CGAL::ROTATION, Input_rational(-1), Input_rational(0)); // rotate by 3/2 PI

  m_scaling = m_is_inner? m_radius: m_radius + m_epsilon;
  Transformation_2 scale = get_scaling(m_scaling);

  Transformation_2 rotate2(CGAL::ROTATION, Input_rational(0), Input_rational(-1)); // rotate by PI
//  Transformation_2 scale2(CGAL::SCALING, Input_rational(1), Input_rational(-1)); // inverse y

//  Transformation_2 scale2(CGAL::SCALING, Input_rational(-1), Input_rational(1)); // inverse

  Polygon_2 unit_polygon;
  log_input_slopes();

  bool first = true;
  Input_rational first_param(0, 1), last_param(0, 1);
  // compute dir normals (outward from tangent circle)
  typename Directed_slope_set::const_iterator dit;
  for(dit = m_input_slopes.begin(); dit != m_input_slopes.end(); ++dit)
  {
    Direction_2 dir = dit->direction();

    LOG_DEBUG << "dir: " << dir << std::endl;

    Direction_2 normal_dir = rotate(dir);
    LOG_DEBUG << "normal_dir: " << normal_dir << std::endl;

    Point_2 below, above;
    bool is_single = false;
    Rational_points_2::approximate_to_point_on_unit_circle(normal_dir, m_epsilon,
      is_single, below, above);

    // TODO: verify eps?
    unit_polygon.push_back(below);
    // TODO: verify "above" points for proper order with next "below"
    if(!is_single)
      unit_polygon.push_back(above);
  }

  LOG_DEBUG << "found " << unit_polygon.size() << " rational points" << std::endl;

  // compare the slopes to the slopes of the original polygon
  Directed_slope_set non_unique_slopes;
  bool need_to_fix = has_non_unique_slopes(unit_polygon, non_unique_slopes);
  
  LOG_DEBUG << (need_to_fix ? "" : "no ") << "need to fix" << std::endl;
  
  // if there is a non-unique slope get rid of it
  // by splitting a problematic segment to 2 segments
  // (with rational mid-point)
  if (need_to_fix) {
    fix_non_unique_slopes_by_segment_splitting(non_unique_slopes, unit_polygon);
  }

  typename Polygon_2::Vertex_const_iterator vit;
  for (vit = unit_polygon.vertices_begin(); vit != unit_polygon.vertices_end(); ++vit) 
  {
    m_kgon.push_back(scale(Point_2(vit->x(), vit->y())));
  }

  CGAL::set_pretty_mode(std::clog);
  LOG_DEBUG << "created the polygon Br:" << std::endl;
  LOG_DEBUG << m_kgon << std::endl;
  LOG_DEBUG << std::endl;
}


template <typename CPolygon_2>
void Circle_approximation_2<CPolygon_2>::random_kgon() 
{
  m_kgon.clear();

  typedef typename CGAL::Random_points_on_circle_2<Point_2> Random_points_on_circle_2;
  typedef typename std::vector<Point_2> PointVector_2;
  
  double radius_one = 1.0;
  Random_points_on_circle_2 point_generator(radius_one);
  LOG_DEBUG << "radius: " << m_radius <<std::endl;
  
  m_scaling = m_is_inner? m_radius: m_radius + m_epsilon;
  Transformation_2 scale = get_scaling(m_scaling);
  LOG_DEBUG << "scaling: " << scale << std::endl;
  
  Directed_slope_set point_sorter;
  //PointVector_2 points;
  size_type count;
  for (count=0; count < m_kgon_size; ++count)
  {
    Point_2 generated_point = *point_generator++;
    Input_rational rad = generated_point.x()*generated_point.x() + generated_point.y()*generated_point.y();
    LOG_DEBUG << "gen point: " << " (" << rad << ")" << generated_point.x() <<" " << generated_point.y() << std::endl;
    Point_2 exact_point =
      Rational_points_2::convert_to_exact_point_on_unit_circle(generated_point);
    rad = exact_point.x()*exact_point.x() + exact_point.y()*exact_point.y(); 
    LOG_DEBUG << "xct point: " << " (" << rad << ")" << exact_point.x() <<" " << exact_point.y() << std::endl;

//	double x = to_double(exact_point.x());
//	double y = to_double(exact_point.y());
//	LOG_DEBUG << "point: " << x <<" " <<y << std::endl;
    
//    LOG_DEBUG << "point: " << exactPoint << std::endl;
//    LOG_DEBUG << "scaled: " << scale(exactPoint) << std::endl;
    point_sorter.insert(Segment_2(CGAL::ORIGIN, scale(exact_point)));
//    points.push_back(scale(exact_point));
  }
  
  // BUG: Slope_set slopes appear to be unsorted
  // (receiving self-intersecting polygon)
  // PATCH: compute convex_hull of points
  /*
  PointVector_2 sorted_points;
  CGAL::convex_hull_2(points.begin(), points.end(), 
                      std::back_inserter(sorted_points));
  
  m_kgon.insert(m_kgon.vertices_begin(), sorted_points.begin(), sorted_points.end());
  */
    typename Directed_slope_set::const_iterator set_iter = point_sorter.begin();
    typename Directed_slope_set::const_iterator set_iter_end = point_sorter.end();
    for(; set_iter != set_iter_end; ++set_iter)
    {
      m_kgon.push_back(set_iter->target());
    }
  
  CGAL::set_pretty_mode(std::clog);
  LOG_DEBUG << "created the polygon Br:" << std::endl;
  LOG_DEBUG << m_kgon << std::endl;
  LOG_DEBUG << std::endl;
  
  Directed_slope_set non_unique_slopes;
  bool need_to_fix = has_non_unique_slopes(m_kgon, non_unique_slopes);

  LOG_DEBUG <<(need_to_fix?"":"no ") <<"need to fix" << std::endl;

  // if need to fix - retry the random points?

}

template <typename CPolygon_2>
void Circle_approximation_2<CPolygon_2>::log_input_slopes() const
{
  typename Directed_slope_set::const_iterator sit;
  for (sit = m_input_slopes.begin(); sit != m_input_slopes.end(); ++sit)
  {
    log_slope("polygon__slope", *sit);
      //LOG_DEBUG << "polygon__slope (" << *sit << ") \t"  << sit->direction() << std::endl;
  }
}


template <typename CPolygon_2>
void Circle_approximation_2<CPolygon_2>::log_circle_edges(
    const Polygon_2& i_polygon) const
{
  typename Polygon_2::Edge_const_iterator  eit;
  
  LOG_DEBUG << " " << i_polygon.size() << " edges:" <<std::endl;
  for (eit = i_polygon.edges_begin(); eit != i_polygon.edges_end(); ++eit)
  {
    Segment_2 seg = *eit;
    bool isCircleSlope = is_circle_slope(seg);
    if(isCircleSlope)
    {
       log_slope("circle__offset", seg);

//      LOG_DEBUG << "circle__offset (" << seg << ") \t"  << seg.direction() << std::endl;
    }
    else
    {
       log_slope("polygon_offset", seg);
//      LOG_DEBUG << "polygon_offset (" << seg << ") \t"  << seg.direction() << std::endl;
    }
  }
  LOG_DEBUG << " " << std::endl;
}

template <typename CPolygon_2>
void Circle_approximation_2<CPolygon_2>::log_circle_edges(
    const Polygon_with_holes_2& i_polygon_with_holes) const
{
  if (!i_polygon_with_holes.is_unbounded())
  {
    LOG_DEBUG << "Outer boundary = " << std::endl;
    log_circle_edges(i_polygon_with_holes.outer_boundary());
  }
  else
    LOG_DEBUG << "Unbounded polygon." << std::endl;
  
  typename Polygon_with_holes_2::Hole_const_iterator  hit;
  unsigned int k = 1;
  
  LOG_DEBUG << "  " << i_polygon_with_holes.number_of_holes()
            << " holes:" << std::endl;
  for (hit = i_polygon_with_holes.holes_begin();
       hit != i_polygon_with_holes.holes_end(); ++hit, ++k)
  {
    LOG_DEBUG << "    Hole #" << k << " = ";
    log_circle_edges(*hit);
  }
  LOG_DEBUG << std::endl;
}

template <typename CPolygon_2>
void Circle_approximation_2<CPolygon_2>::kgon_size_from_eps(
  const Input_rational& epsilon, 
  const Input_rational& radius, 
  Polygon_size& kgon_size)
{
  LOG_DEBUG << "kgon_size_from_eps(" << epsilon << ", " << radius << ", " << kgon_size << ")";
  // k >= PI / arccos(1 - eps/r)
  kgon_size = (CGAL_PI / std::acos(to_double(1 - epsilon / radius))) + 1;
  LOG_DEBUG << "eps " << epsilon << " => kgon_size " << kgon_size << std::endl;
}

template <typename CPolygon_2>
void Circle_approximation_2<CPolygon_2>::eps_from_kgon_size(
  const Polygon_size& kgon_size, 
  const Input_rational& radius, 
  Input_rational& epsilon)
{
  LOG_DEBUG << "eps_from_kgon_size(" << kgon_size << ", " << radius << ", " << kgon_size << ")" << std::endl;
  
    // eps >= r * (1 - cos (PI/k))
    epsilon = radius * (1 - std::cos(CGAL_PI / kgon_size) );
  LOG_DEBUG << " epsilon = radius * (1 - std::cos(CGAL_PI / kgon_size) ) = " << epsilon << std::endl;
    
//    epsilon = radius * std::sin(CGAL_PI * ((kgon_size - 2) / (2 * kgon_size)));
    int kgon_size_root = int( ceil(sqrt(double(kgon_size))) );
  LOG_DEBUG << " kgon_size_root = int( ceil(sqrt(double(kgon_size))) ) = " << kgon_size_root << std::endl;

    if(epsilon <= 0) epsilon = radius / kgon_size_root; // tight enough eps = radius / sqrt(kgon_size)
    LOG_DEBUG << "kgon_size " << kgon_size << "  => eps " << epsilon << std::endl;
}

} // namespace CGAL

#endif // CGAL_CIRCLE_APPROXIMATION_2_H
