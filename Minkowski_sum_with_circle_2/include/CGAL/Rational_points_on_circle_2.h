// Rational_points_on_circle_2

#ifndef Rational_points_on_circle_2_H
#define Rational_points_on_circle_2_H

#include <CGAL/Fraction_traits.h>
#include <CGAL/Polygon_2.h>
//#include <ext/hash_map>
//namespace std { using namespace __gnu_cxx; }

#include "output_debug.h"

namespace CGAL {

/*! \classf
 * A class implementing selection of a set of rational
 * points or single rational point on a 2D circle 
 * according to given parameters.
 * The points are selected on a unit circle and then 
 * are scaled according to required circle radius. 
 */

#define DEC 10
#define BIN 2
#define BASE DEC

//#define DEC_DEGREE 6    // 1 000 000
//#define BIN_DEGREE 20   // 1 048 576
#define DEC_DEGREE 5    // 100 000
#define BIN_DEGREE 17   // 131 072

#define BASE_DEGREE DEC_DEGREE

template <typename Traits>
class Rational_points_on_circle_2
{
public:

  typedef Traits Traits_2;
  typedef typename CGAL::Polygon_2<Traits_2> Polygon_2;
  typedef Rational_points_on_circle_2<Traits>  Rational_points_2;
  
  typedef typename Traits_2::FT Input_rational;
  typedef typename Traits_2::Point_2   Point_2; 
  typedef typename Traits_2::Direction_2   Direction_2;

  typedef typename std::map<Input_rational, Point_2> Points_on_circle_2;
//  typedef typename std::hash_map<Input_rational, Point_2> Points_on_circle_2;

  typedef typename std::map<int, Points_on_circle_2> Degree_to_points_on_circle_2;
  
  struct Quadrant_2
  {
    bool neg_x;
    bool neg_y;
    Quadrant_2(): neg_x(false), neg_y(false) {}
    Quadrant_2(bool x_is_neg, bool y_is_neg): neg_x(x_is_neg), neg_y(y_is_neg) {}
  };

  struct Rational_point_2
  {
    Input_rational x;
    Input_rational y;

    Quadrant_2 quadrant;

    Rational_point_2(const Point_2& i_point) :
      x(i_point.x()),
      y(i_point.y()),
      quadrant(x < 0, y < 0)
    {}

    Rational_point_2(const Input_rational& i_x, const Input_rational& i_y) :
      x(i_x),
      y(i_y),
      quadrant(x < 0, y < 0)
    {}

    void into_first_quadrant() 
    { 
      quadrant.neg_x = (x < 0);
      quadrant.neg_y = (y < 0);

      flip_x_if_needed();
      flip_y_if_needed();
      swap_if_needed();
    }

    void from_first_quadrant() 
    {
      swap_if_needed();
      flip_y_if_needed();
      flip_x_if_needed();
    }

    void swap_if_needed() 
    { 
      if(quadrant.neg_x^quadrant.neg_y) 
      {
        Input_rational temp = x; x = y; y = temp;
      } 
    }

    void flip_x_if_needed() { if(quadrant.neg_x) x = -x; }
    void flip_y_if_needed() { if(quadrant.neg_y) y = -y; }
  };
  
  Rational_points_on_circle_2(const Input_rational& i_radius,
    const Input_rational& i_epsilon) :
    m_radius(i_radius),
    m_epsilon(i_epsilon),
    m_eps_hat(i_epsilon / i_radius)
  {
  }

  /*!
   * Convert mashine-precision point "on circle" to the rational point
   * lying exactly on the circle. The exact point will be on a line trough
   * (-1, 0) an the inexact point (the bit-size of received coordinates
   * will be relatively big).
   * \param i_point Not exact circle point.
   */
  static Point_2 convert_to_exact_point_on_unit_circle(const Point_2& i_point);
  
  /*!
   * Convert given direction to the rational point(s)
   * lying exactly on the circle. The exact point(s) will be close up to
   * epsilon to the exact point on the circle in a given direction.
   * The method used is the one from Ron's article: todo - describe.
   * \param i_dir The direction of the exact circle point.
   * \param i_eps The maximum distance between the exact point and the chord connecting two approximations.
   * \param is_exact true if single point (exactly in the given direction) was returned
   * \param o_below the output point below the input
   * \param o_above the output point above the input
   */  
  static void approximate_to_point_on_unit_circle(const Direction_2& i_dir, const Input_rational& i_eps,
    bool& is_exact, Point_2& o_below, Point_2& o_above);

  /*!
   * Mirror rational points found in the first quadrant to other quadrants.
   */
  static void mirror_first_quadrant(const Polygon_2& i_polygon, Polygon_2& o_polygon);  
  
  void update_ratio();
 
  /// get squared distance between points
  static Input_rational get_sqr_dist(const Point_2& point1, const Point_2& point2)
  {
    Input_rational x_diff = point1.x() - point2.x();
    Input_rational y_diff = point1.y() - point2.y();
    Input_rational sqr_dist = x_diff*x_diff + y_diff*y_diff;
    return sqr_dist;
  }

  /// precompute all rational points on circle for various epsilon requirements
  static void precompute_short_rationals() 
  {
    if(!m_degree_to_points.empty()) return;

    static const int MAX_DEGREE = BASE_DEGREE; //5; //2;
    int denom = BASE;
    for(int degree = 1; degree <= MAX_DEGREE; ++degree)
    {
      Points_on_circle_2 points;
      m_degree_to_points[degree] = points;
      precompute_short_rationals(denom, &m_degree_to_points[degree]);
      denom *= BASE;
    }
  }
 
  /*! 
   * find short rationals in the first quadrant with at most delta distance 
   * between each two points (i.e. hold the epsilon)
  */
  void get_first_quadrant_rationals(Polygon_2& o_first_quadrant_polygon);

//  static void map_into_first_quadrant(Point_2& io_point, Quadrant_2& o_quadrant);
//  static void map_from_first_quadrant(Point_2& io_point, const Quadrant_2& i_quadrant);  
  
private:

  // remove common divisors of two numbers
  static void remove_common_factor(Input_rational& num, Input_rational& denom);

  // converging series for sqrt approximation
  static void next_in_series(const int& index, const Input_rational& sqr_target, 
    const Input_rational& prev_in_series, Input_rational& next_in_series);

  // approximate sqrt from below/above up to eps
  static void approximate_sqrt(const Input_rational& sqr_target, 
    const Input_rational& eps, Input_rational& app_sqrt, 
    const bool& from_below = true, const Input_rational& sqrt_begin = 1);

  /// convert slope param to point on circle
  static Point_2 get_param_point(const Input_rational& param)
  {
    Input_rational sqr_param = param*param;
    Input_rational one_plus_sqr_param = 1 + sqr_param;
    Point_2 param_point((1 - sqr_param)/one_plus_sqr_param, 2*param/one_plus_sqr_param);
    return param_point;
  }

  /// convert point on circle to slope param
  static Input_rational get_point_param(const Point_2& point)
  {
    Input_rational param = point.y() / (point.x() + 1);
    return param;
  }

  /// convert degree of rationals (necessary for the epsilon) to denominator of param t
  static int deg_to_denom(const int& deg)
  {
    int denom = BASE;
    for(int degree = 1; degree < deg; ++degree)
    {
      denom *= BASE;
    }
    return denom;
  }

  /// precompute all rational points on circle for the given epsilon requirements
  static void precompute_short_rationals(const int& target, Points_on_circle_2* p_map);  

  /// get all short rationals with param t degree (t = BASE^-degree)
  static void get_short_rationals(const int& degree, Points_on_circle_2** pp_map);
  /// find degree of rationals necessary for the epsilon
  int get_rationals_degree();
  /// get maximal (delta) distance between rational points that does not violate the epsilon
  Input_rational get_sqr_rationals_dist();

  // split the interval of given param
  void split_short_rational(const Input_rational& param, Input_rational& point_t, Point_2& point)
  {
    // TODO: split param with its predecessor
  }  

  const Input_rational& m_radius;
  const Input_rational& m_epsilon;    
  Input_rational m_eps_hat;  
  
  static Degree_to_points_on_circle_2 m_degree_to_points;  
};

template <typename Traits>
typename Rational_points_on_circle_2<Traits>::Degree_to_points_on_circle_2 
Rational_points_on_circle_2<Traits>::m_degree_to_points;

template <typename Traits>
void Rational_points_on_circle_2<Traits>::update_ratio()
{
  m_eps_hat = m_epsilon / m_radius;
}

template <typename Traits>
typename Rational_points_on_circle_2<Traits>::Point_2
Rational_points_on_circle_2<Traits>::convert_to_exact_point_on_unit_circle(
  const Point_2& i_point)
{
  Input_rational inexact_x = i_point.x();
  Input_rational inexact_y = i_point.y();
  
  /// Parametrize rational points on a unit circle by a slope t of the line
  /// passing thru (-1, 0) and a rational circle point (i.e. y = t * (1 + x))
  /// This way the rational point will have next coordinates:
  /// ( (1 - t^2) / (1 + t^2) , 2t / (1 + t^2) )

  Point_2 exact_point;
  // take care of not dividing by 0 and not missing (-1, 0) point
  if (inexact_x == -1 || (inexact_y == 0 && inexact_x < 0))
  {
    exact_point = Point_2(-1, 0);
  }
  else
  {
    Input_rational slope = inexact_y/(inexact_x + 1);
    Input_rational slope_sq = slope*slope;

    //approximate_from_below()
    
    Input_rational exact_x = (1 - slope_sq)/(1 + slope_sq);
    Input_rational exact_y = (2 * slope)/(1 + slope_sq);
    
    exact_point = Point_2(exact_x, exact_y);
  }
  
  return exact_point;
}
/*
template <typename Traits>
void Rational_points_on_circle_2<Traits>::
map_into_first_quadrant(Point_2& io_point, Quadrant_2& o_quadrant)
{
  Rational_point_2 p(io_point);
  p.into_first_quadrant();
  o_quadrant = p.quadrant;
  io_point = Point_2(p.x, p.y);
}

template <typename Traits>
void Rational_points_on_circle_2<Traits>::
map_from_first_quadrant(Point_2& io_point, const Quadrant_2& i_quadrant)
{
  Rational_point_2 p(io_point);
  p.quadrant = i_quadrant;
  p.from_first_quadrant();
  io_point = Point_2(p.x, p.y);
}
*/


template <typename Traits>
void Rational_points_on_circle_2<Traits>::
remove_common_factor(Input_rational& num, Input_rational& denom)
{
  /*
  Input_rational ratio = x/y;
  x = numerator(ratio);
  y = denominator(ratio);
  */
  typedef typename CGAL::Fraction_traits<Input_rational> FracTraits;
  typename FracTraits::Numerator_type numerator;
  typename FracTraits::Denominator_type denominator(1);
  typename FracTraits::Decompose decomposer;

  decomposer(num/denom, numerator, denominator);
  num = numerator;
  denom = denominator;
}

// converging series for sqrt approximation
template <typename Traits>
void Rational_points_on_circle_2<Traits>::
next_in_series(const int& index, const Input_rational& sqr_target, 
    const Input_rational& prev_in_series, Input_rational& next_in_series)
{
  next_in_series = (prev_in_series + sqr_target/prev_in_series)/2;
  //LOG_DEBUG << "\t[" << index << "]: " << next_in_series  
  //  << "(" << to_double(next_in_series) << ")" << std::endl;
}

// approximate sqrt from below/above up to eps
template <typename Traits>
void Rational_points_on_circle_2<Traits>::
approximate_sqrt(const Input_rational& sqr_target, 
    const Input_rational& eps, Input_rational& app_sqrt, 
    const bool& from_below, const Input_rational& sqrt_begin)
{
    app_sqrt = 0;
    const Input_rational sqr_sqr_target = sqr_target * sqr_target;

    const int MAX_INDEX = 10;
    Input_rational next = sqrt_begin;
    int index = 0;
    bool val = false;
    bool inverse_val = false;

    LOG_DEBUG << "sqr target: " << sqr_target << " eps: " << eps << std::endl;
    do
    {
      Input_rational prev_in_series = next;
      next_in_series(++index, sqr_target, prev_in_series, next);

      Input_rational sqr_next = next * next;
      //LOG_DEBUG << "\t" << to_double(sqr_next);

      Input_rational val_diff = sqr_target - sqr_next;
      if(!from_below)
        val_diff = -val_diff;

      //LOG_DEBUG << "\t" << to_double(val_diff);

      if(val_diff >= 0)
        val = (val_diff <= eps);

      if(!val)
      {
        Input_rational inv_val_diff = (sqr_target - sqr_sqr_target/sqr_next);
        if(!from_below)
          inv_val_diff = -inv_val_diff;

        //LOG_DEBUG << "\t" << to_double(inv_val_diff);

        if(inv_val_diff >= 0)
          inverse_val = (inv_val_diff <= eps);
      }

      //LOG_DEBUG << std::endl;

    }while((!val && !inverse_val) && (index <= MAX_INDEX));

    if(val)
    {
      app_sqrt = next;
    }
    else
    //if(inverse_val)
    {
      app_sqrt = sqr_target/next;
    }
}


template <typename Traits>
void Rational_points_on_circle_2<Traits>::
approximate_to_point_on_unit_circle(const Direction_2& i_dir, 
                                    const Input_rational& i_eps,
                                    bool& is_exact, 
                                    Point_2& o_below, 
                                    Point_2& o_above)
{

  // point on circle
  Rational_point_2 poc(i_dir.dx(), i_dir.dy());

  // handle special cases
  if(poc.y == 0)
  {
    Input_rational unit_y = 0;
    Input_rational unit_x = (poc.quadrant.neg_x)? -1: 1;

    is_exact = true;
    o_below = Point_2(unit_x, unit_y);
    o_above = o_below;

    return;
  }
  else
  if(poc.x == 0)
  {
    Input_rational unit_x = 0;
    Input_rational unit_y = (poc.quadrant.neg_y)? -1: 1;

    is_exact = true;
    o_below = Point_2(unit_x, unit_y);
    o_above = o_below;

    return;
  }

  assert(poc.y != 0);
  poc.into_first_quadrant();

  /*   
  double d_x = to_double(x);
  double d_y = to_double(y);
  double d_len = std::sqrt(d_x*d_x + d_y*d_y);
  double d_slope = (d_len - d_x)/d_y;

  // eps and rp max dist for the unit circle 
  Input_rational eps_hat = m_epsilon / m_radius;
  Input_rational eps_hat_inverse = 1 / eps_hat;

  // find eps <= BASE^(-degree*2) 
  int denom_degree = 1;
  int denom = BASE;
  Input_rational denom_sqr = denom * denom;
  while(denom_sqr < eps_hat_inverse)
  {
    ++denom_degree;
    denom *= BASE;
    denom_sqr = denom * denom;
  }

  int param_nom = int(d_slope * denom);
  Input_rational param1(param_nom, denom);
  Input_rational param2(param_nom + 1, denom);
  */
    
  remove_common_factor(poc.x, poc.y);

  // rational points on unit circle parametrization:
  // t = y/(1+x) = (1-x)/y => for vector of length d = sqrt(x^2 + y^2)
  // t = (d-x/y) can be not rational. to get proper delta between our 
  // rationals we need to compute t with accuracy y*sqrt(eps-hat).
  // therefore we compute two rational points with half of this accuracy 
  // from above and from below to get dist at most delta between them
  Input_rational& x = poc.x;
  Input_rational& y = poc.y;

  Input_rational target = x*x + y*y;
  Input_rational app_val = x + y;
  Input_rational sqr_eps = (i_eps*y)/2;

  Input_rational below = 0;
  Input_rational above = 0;
  approximate_sqrt(target, sqr_eps, below, true, app_val);
  LOG_DEBUG << "below app: " << below << " with diff " << (target - below*below) << std::endl;
  approximate_sqrt(target, sqr_eps, above, false, app_val);
  LOG_DEBUG << "above app: " << above << " with diff " << (above*above - target) << std::endl;

  if(below == above)
  {
    is_exact = true;

    Input_rational param = (below - x)/y;
    LOG_DEBUG << "above = below => " << param << std::endl;
    Input_rational sqr = param*param;
    Input_rational unit_x = (1 - sqr)/(1 + sqr);
    Input_rational unit_y = 2*param/(1 + sqr);

    poc.x = unit_x;
    poc.y = unit_y;
    poc.from_first_quadrant();

    o_below = Point_2(poc.x, poc.y);
    o_above = o_below;
  }
  else      
  {
    is_exact = false;

    Input_rational param1 = (below - x)/y;
    Input_rational param2 = (above - x)/y;
    LOG_DEBUG << "above != below => " << param2 <<" - " <<param1 
      << " = " << to_double(param2 - param1) << std::endl;

    Input_rational param = param1;
    Input_rational sqr = param*param;
    Input_rational unit_x = (1 - sqr)/(1 + sqr);
    Input_rational unit_y = 2*param/(1 + sqr);

    poc.x = unit_x;
    poc.y = unit_y;
    poc.from_first_quadrant();

    o_below = Point_2(poc.x, poc.y);

    param = param2;
    sqr = param*param;
    unit_x = (1 - sqr)/(1 + sqr);
    unit_y = 2*param/(1 + sqr);

    poc.x = unit_x;
    poc.y = unit_y;
    poc.from_first_quadrant();

    o_above = Point_2(poc.x, poc.y);
  }
}

template <typename Traits>
void Rational_points_on_circle_2<Traits>::
mirror_first_quadrant(
  const Polygon_2& i_polygon, 
  Polygon_2& o_polygon)
{
  o_polygon = i_polygon;
  typename Polygon_2::Vertex_const_iterator vit;
  for (vit = i_polygon.vertices_begin(); vit != i_polygon.vertices_end(); ++vit) 
  {
    o_polygon.push_back(Point_2(- vit->y(), vit->x()));
  }
  for (vit = i_polygon.vertices_begin(); vit != i_polygon.vertices_end(); ++vit) 
  {
    o_polygon.push_back(Point_2(- vit->x(), - vit->y()));
  }
  for (vit = i_polygon.vertices_begin(); vit != i_polygon.vertices_end(); ++vit) 
  {
    o_polygon.push_back(Point_2(vit->y(), - vit->x()));
  }
}

template <typename Traits>
void Rational_points_on_circle_2<Traits>::
precompute_short_rationals(
  const int& target, 
  Points_on_circle_2* p_map)
{
  LOG_DEBUG << "precompute denom " << target << std::endl;
  assert(p_map);

  Points_on_circle_2& poc_map = *p_map;
  poc_map.clear();

  CGAL::Timer timer;
  timer.start();

  Point_2 first_point = Point_2(1, 0);
  Point_2 last_point = Point_2(0, 1);
  poc_map[0] = first_point;
  poc_map[1] = last_point;

  // Iterate over all nominators of target denominator 
  // and compute rational points on circle
  int denom = target;

  // Note: Probably no point in iterating over all smaller
  // denominators for computing shorter rationals, because
  // it will create an overhead of too many points (squared 
  // number of what we compute for one numerator). May be 
  // worth testing.
 
//  int start_target = 2;
// //  if(target >= 1000) 
//    start_target = target; // -2;
//  for(int denom = start_target; denom <= target; ++denom)
//  {
    for(int nom = 1; nom < denom; ++nom)
    {
      Input_rational param(nom, denom);
      // no reason to use "find" here:
      // point computation is faster then search in the map
      //typename Points_on_circle_2::const_iterator it = poc_map.find(param);
      //if(it == poc_map.end())
      {
        poc_map[param] = get_param_point(param);
      }
    }
//  }
  timer.stop();

//    LOG_DEBUG << "Number of points = " << poc_map.size() << std::endl;

/*    // write map to file
  std::ostringstream fileNameStream;
  fileNameStream << "rationals_" << target << ".txt";
  std::string rationalsFile = fileNameStream.str();
  std::ofstream out_file(rationalsFile.c_str());

  Points_on_circle_2::const_iterator it;
  Point_2 prev_point = first_point;
/ *
  Input_rational prev_point_t = 0;
  Input_rational min_sqr_dist = 2;
  Input_rational min_point = 1;
  Input_rational min_point2 = 0;
  Input_rational max_sqr_dist = 0;
  Input_rational max_point = 0;
  Input_rational max_point2 = 0;
* /
//    bool first = true;
  for ( it=poc_map.begin() ; it != poc_map.end(); ++it)
  {
    Point_2 curr_point = (*it).second;
    Input_rational curr_point_t = (*it).first;

    Input_rational sqr_dist = get_sqr_dist(curr_point, prev_point);

//      out_file << (*it).first << "\t" << (*it).second << "\t" << sqr_dist << endl;
/ *
    if(first)
    {
      first = false;
      continue;
    }

    if(min_sqr_dist > sqr_dist)
    {
      min_sqr_dist = sqr_dist;
      min_point = curr_point_t;
      min_point2 = prev_point_t;
    }

    if(max_sqr_dist < sqr_dist)
    {
      max_sqr_dist = sqr_dist;
      max_point = curr_point_t;
      max_point2 = prev_point_t;
   }
* /
    prev_point = curr_point;
//      prev_point_t = curr_point_t;
  }

//    LOG_DEBUG << "Min dist = " << min_sqr_dist << ", Min t_diff = ("
//      << min_point << "," <<min_point2 <<")" << std::endl;
//    LOG_DEBUG << "Max dist = " << max_sqr_dist << ", Max t_diff = ("
//      << max_point << "," <<max_point2 <<")"<< std::endl;
*/
   LOG_DEBUG << "RPOC of size " << poc_map.size()
     << " computed in " << timer.time() << std::endl;
/*
#ifdef SHOULD_LOG
  // log distances of small rat point maps
  if(poc_map.size() < 100)
  {
    typename Points_on_circle_2::const_iterator it;
    Point_2 prev_point = first_point;
    for ( it=poc_map.begin() ; it != poc_map.end(); ++it)
    {
      Point_2 curr_point = (*it).second;
      Input_rational sqr_dist = get_sqr_dist(curr_point, prev_point);

      LOG_DEBUG << "\t" << (*it).first << "\t" << (*it).second << "\t" << sqr_dist << std::endl;
      prev_point = curr_point;
    }
  }
#endif
*/
}


template <typename Traits>
void Rational_points_on_circle_2<Traits>::
get_short_rationals(
  const int& degree, 
  Points_on_circle_2** pp_map)
{
  typename Degree_to_points_on_circle_2::iterator it = 
    m_degree_to_points.find(degree);
  if(it == m_degree_to_points.end())
  {
    Points_on_circle_2 points;
    m_degree_to_points[degree] = points;

    precompute_short_rationals(deg_to_denom(degree), 
      &m_degree_to_points[degree]);
    it = m_degree_to_points.find(degree);
  }
 
  *pp_map = &(it->second); //&m_degree_to_points[degree];
  LOG_DEBUG << *pp_map << std::endl;
}

/// find degree of rationals necessary for the epsilon
template <typename Traits>
int Rational_points_on_circle_2<Traits>::
get_rationals_degree()
{
  // eps and rp max dist for the unit circle 
  Input_rational eps_hat = m_epsilon / m_radius;
  Input_rational eps_hat_inverse = 1 / eps_hat;

  LOG_DEBUG << "RPOC: 1/ eps_hat = " << eps_hat_inverse << std::endl;

  // find eps <= BASE^(-degree*2)
  int denom_degree = 1;
  Input_rational denom = BASE;
  Input_rational denom_sqr = denom * denom;
  while(denom_sqr < eps_hat_inverse)
  {
    ++denom_degree;
    denom *= BASE;
    denom_sqr = denom * denom;
  }

  LOG_DEBUG << "RPOC: lambda = " << denom << "(= " << BASE << "^" << denom_degree <<")" << std::endl << std::endl;
  return denom_degree;
}

/// get maximal (delta) distance between rational points on circle that does not violate the epsilon
template <typename Traits>
typename Rational_points_on_circle_2<Traits>::Input_rational Rational_points_on_circle_2<Traits>::
get_sqr_rationals_dist()
{
  // unit circle equations:
  // 1) x + eps_hat = 1 => x = 1 - eps_hat
  // 2) y = sqrt(1 - x^2) => y^2 = 1 - (1 - eps_hat)^2 = 2*eps_hat - eps_hat^2
  // 3) delta = 2 * y => delta^2 = 4 * eps_hat * (2 - eps_hat)
  
  Input_rational eps_hat = m_epsilon / m_radius;
  return 4 * eps_hat * (2 - eps_hat);
}


/*! 
 * Find short rationals in the first quadrant with the 
 * epsilon (at most delta distance between each two).
 *
 * Iterate over precomputed set of rational points selecting 
 * each time the farthest point with distance below delta.
 * 
 * Iterate with next_point computing distance to the previous 
 * inserted point prev_point. When distance to next_point gets 
 * bigger then delta curr_point is supposed to be the farthest 
 * point with distance smaller then delta.
 * 
 * If distances in the initial set are not guaranteed to be 
 * smaller then delta, it is necessary to introduce splitting
 * of greater then delta intervals. In such case curr_point
 * will be equal to prev_point.
 *
 * TODO: use the fact that in our sets distances between points 
 * grow smaller and smaller to iterate more efficiently. 
 * For example jump at least number of previously skipped points 
 * when looking for the next point.
 */
template <typename Traits>
void Rational_points_on_circle_2<Traits>::
get_first_quadrant_rationals(Polygon_2& o_first_quadrant_polygon)
{
  // eps and rp max dist for the unit circle 
  Input_rational eps_hat = m_epsilon / m_radius;

  // get a set of points we choose from (guarantees eps via given denominator size)
  Points_on_circle_2* p_points = NULL;
  int denom_degree = get_rationals_degree();
  get_short_rationals(denom_degree, &p_points);
  //LOG_DEBUG << "short_rationals = " << p_points << std::endl;
  Points_on_circle_2& m_rat_points = *p_points;

  CGAL_assertion(!m_rat_points.empty());
 
  // get maximal allowed distance between points
  Input_rational delta_sqr = get_sqr_rationals_dist();
  LOG_DEBUG << "sqr_delta = " << delta_sqr << std::endl;

  // insert the first point (1, 0)
  Point_2 curr_point = m_rat_points[0];
  typename Points_on_circle_2::size_type rat_count = 1;
  o_first_quadrant_polygon.push_back(curr_point);
  Point_2 prev_point = curr_point;

  // find farthest point within delta distance
  typename Points_on_circle_2::const_iterator it;
  int log_point_counter = 0;
  for ( it = m_rat_points.begin() ; it != m_rat_points.end(); ++it)
  {
    Point_2 next_point = (*it).second;
    Input_rational next_dist_sqr = get_sqr_dist(next_point, prev_point);

    if(log_point_counter < 100)
      LOG_DEBUG << next_dist_sqr << " ";
    
    if(next_dist_sqr > delta_sqr)
    {
      if(log_point_counter < 100)
      {
        log_point_counter++;
        LOG_DEBUG << " => "<< curr_point << std::endl;
      }
      bool need_to_split = (curr_point == prev_point);
      if(!need_to_split)
      // prev is our guy
      {
        o_first_quadrant_polygon.push_back(curr_point);
        ++rat_count;
        prev_point = curr_point;
 
        next_dist_sqr = get_sqr_dist(next_point, prev_point);
        need_to_split = (next_dist_sqr > delta_sqr);
      }

      // TODO: implement splitting for sparse sets 
      // (with more then delta distance between points) 
      if(need_to_split)
      // look for (next_t + prev_t)/2 till in delta dist
      {
        // curr_point = mid(next, prev);
        ERR_DEBUG << "!!! split not implemented yet !!!" << std::endl;
        //exit 1;
      }
    }

    curr_point = next_point;
  }

  LOG_DEBUG << "first quadrant points = " << o_first_quadrant_polygon.size() << " of " << m_rat_points.size() << std::endl;

}

} // namespace CGAL

#endif // Rational_points_on_circle_2_H