#ifndef CGAL_BEZIER_TRAITS_CACHE_H
#define CGAL_BEZIER_TRAITS_CACHE_H

/*! \file
 * Header file for the _Bezier_cache class.
 */

#include <set>
#include <map>
#include <ostream>
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_d_1.h>
#include <CGAL/Polynomial.h>


namespace CGAL {

/*! \class _Bezier_cache
 * Stores all cached intersection points and vertical tangency points.
 */

template <class AlgebraicKernel_d_1_> 
class _Bezier_cache
{
public:

  // Setting the typedefs for the various type names
  typedef AlgebraicKernel_d_1_								Algebraic_kernel_d_1;

  typedef typename Algebraic_kernel_d_1::Coefficient					Coefficient;

  // TODO Maybe just use the Coefficient type for Integer later
  typedef typename Get_arithmetic_kernel<Coefficient>::Arithmetic_kernel::Integer	Integer;

  typedef typename Algebraic_kernel_d_1::Polynomial_1					Polynomial_1;
  typedef typename Algebraic_kernel_d_1::Algebraic_real_1				Algebraic_real_1;

  typedef typename CGAL::Polynomial_traits_d<Polynomial_1>				Polynomial_traits_1;


  typedef _Bezier_cache<Algebraic_kernel_d_1>                	Self;

  /// \name Type definitions for the vertical tangency-point mapping.
  //@{
  typedef size_t                                     		Curve_id;
  typedef std::list<Algebraic_real_1>		     		Vertical_tangency_list;
  typedef typename Vertical_tangency_list::const_iterator  	Vertical_tangency_iter;
  //@}

  /// \name Type definitions for the intersection-point mapping.
  //@{

  /*! \struct Intersection_point_2
   * Representation of an intersection point (in both parameter and physical
   * spaces).
   */
  struct Intersection_point_2
  {
    Algebraic_real_1          s;      // The parameter for the first curve.
    Algebraic_real_1          t;      // The parameter for the second curve.
    Algebraic_real_1          x;      // The x-coordinate.
    Algebraic_real_1          y;      // The y-coordinate.

    /*! Constructor. */
    Intersection_point_2 (const Algebraic_real_1& _s, const Algebraic_real_1& _t,
                          const Algebraic_real_1& _x, const Algebraic_real_1& _y) :
      s(_s), t(_t),
      x(_x), y(_y)
    {}
  };

  typedef std::pair<Curve_id, Curve_id>              		   Curve_pair;
  typedef std::pair<Algebraic_real_1, Algebraic_real_1>            Parameter_pair;
  typedef std::list<Intersection_point_2>            		   Intersection_list;
  typedef typename Intersection_list::const_iterator		   Intersection_iter;

private:

  typedef std::map<Curve_id, Vertical_tangency_list> Vert_tang_map;
  typedef typename Vert_tang_map::value_type         Vert_tang_map_entry;
  typedef typename Vert_tang_map::iterator           Vert_tang_map_iterator;

  /*! \struct Less_curve_pair
   * An auxiliary functor for comparing pair of curve IDs.
   */
  struct Less_curve_pair
  {
    bool operator() (const Curve_pair& cp1, const Curve_pair& cp2) const
    {
      // Compare the pairs of IDs lexicographically.
      return (cp1.first < cp2.first ||
              (cp1.first == cp2.first && cp1.second < cp2.second));
    }
  };

  typedef std::list<Algebraic_real_1>                Parameter_list;
  typedef std::pair<Intersection_list, bool>         Intersection_info;
  typedef std::map<Curve_pair,
                   Intersection_info,
                   Less_curve_pair>                  Intersect_map;
  typedef typename Intersect_map::value_type         Intersect_map_entry;
  typedef typename Intersect_map::iterator           Intersect_map_iterator;

  /*! \struct My_point_2
   * Representation of exact and approximate point coordinates.
   */
  struct My_point_2
  {
    typename Parameter_list::const_iterator  prm_it;
    Algebraic_real_1                         x;
    Algebraic_real_1                         y;
    double                                   app_x;
    double                                   app_y;

    /*! Default constructor. */
    My_point_2 () :
      app_x (0),
      app_y (0)
    {}

    /*! Constructor. */
    My_point_2 (typename Parameter_list::const_iterator it,
                const Algebraic_real_1& _x, const Algebraic_real_1& _y) :
      prm_it (it),
      x (_x),
      y (_y),
      //app_x (_x.to_double()),
      //app_y (_y.to_double())
      app_x (CGAL::to_double(_x)),
      app_y (CGAL::to_double(_y))
    {}

    /*! Get the parameter value. */
    const Algebraic_real_1& parameter () const
    {
      return (*prm_it);
    }

    /*! Check if the given two points are equal. */
    bool equals (const My_point_2& other) const
    {
	typename Algebraic_kernel_d_1::Compare_1 compare_1 = algebraic_kernel_d_1.compare_1_object();

	return (compare_1 (x, other.x) == EQUAL && compare_1 (y, other.y) == EQUAL);
        //return (CGAL::compare (x, other.x) == EQUAL &&
        //        CGAL::compare (y, other.y) == EQUAL);
    }
  };

  /*! \struct Less_distance_iter
   * An auxiliary functor for comparing approximate distances.
   */
  typedef std::list<My_point_2>                      Point_list;
  typedef typename Point_list::iterator              Point_iter;
  typedef std::pair<double, Point_iter>              Distance_iter;

  struct Less_distance_iter
  {
    bool operator() (const Distance_iter& dit1,
                     const Distance_iter& dit2) const
    {
      return (dit1.first < dit2.first);
    }
  };


  // Data members:
  static Algebraic_kernel_d_1	algebraic_kernel_d_1;
  // Nt_traits         	nt_traits;        /*! Number-type traits. */

  Vert_tang_map     	vert_tang_map;    /*! Maps curves to their vertical
                                          tangency parameters. */
  Intersect_map     	intersect_map;    /*! Maps curve pairs to their
                                          intersection parameter pairs. */

  // Copy constructor and assignment operator - not supported.
  _Bezier_cache (const Self& );
  Self& operator= (const Self& );

public:

  /*! Constructor. */
  _Bezier_cache ()
  {}

  /*!
   * Get the vertical tangency parameters of a given curve (X(t), Y(t)).
   * \param id The curve ID.
   * \param polyX The coefficients of X(t).
   * \param normX The normalizing factor of X(t).
   * \return A list of parameters 0 < t_1 < ...< t_k < 1 such that X'(t_i) = 0.
   */
  const Vertical_tangency_list&
  get_vertical_tangencies (const Curve_id& id,
                           const Polynomial_1& polyX, const Integer& normX);

  /*!
   * Get the intersection parameters of two given curve (X_1(s), Y_1(s))
   * and (X_2(t), Y_2(t)).
   * \param id1 The ID of the first curve.
   * \param polyX_1 The coefficients of X_1.
   * \param normX_1 The normalizing factor of X_1.
   * \param polyY_1 The coefficients of Y_1.
   * \param normY_1 The normalizing factor of Y_1.
   * \param id2 The ID of the second curve.
   * \param polyX_2 The coefficients of X_2.
   * \param normX_2 The normalizing factor of X_2.
   * \param polyY_2 The coefficients of Y_2.
   * \param normY_2 The normalizing factor of Y_2.
   * \param do_ovlp Output: Do the two curves overlap.
   * \pre id1 < id2 (swap their order if this is not the case).
   * \return A list of parameter pairs (s_1, t_1), ..., (s_k, t_k) such that
   *         X_1(s_i) = X_2(t_i) and Y_1(s_i) = Y_2(t_i).
   *         Each pair is also associated with the physical point coordinates.
   */
  const Intersection_list&
  get_intersections (const Curve_id& id1,
                     const Polynomial_1& polyX_1, const Integer& normX_1,
                     const Polynomial_1& polyY_1, const Integer& normY_1,
                     const Curve_id& id2,
                     const Polynomial_1& polyX_2, const Integer& normX_2,
                     const Polynomial_1& polyY_2, const Integer& normY_2,
                     bool& do_ovlp);

  /*!
   * Mark two curves as overlapping.
   * \param id1 The ID of the first curve.
   * \param id2 The ID of the second curve.
   */
  void mark_as_overlapping (const Curve_id& id1,
                            const Curve_id& id2);

private:

  /*!
   * Compute all s-parameter values such that the curve (X_1(s), Y_1(s))
   * intersects with (X_2(t), Y_2(t)) for some t.
   * \param polyX_1 The coefficients of X_1.
   * \param normX_1 The normalizing factor of X_1.
   * \param polyY_1 The coefficients of Y_1.
   * \param normY_1 The normalizing factor of Y_1.
   * \param polyX_2 The coefficients of X_2.
   * \param normX_2 The normalizing factor of X_2.
   * \param polyY_2 The coefficients of Y_2.
   * \param normY_2 The normalizing factor of Y_2.
   * \param s_vals Output: A list of s-values such that (X_1(s), Y_1(s))
   *                       is an intersection point with (X_2(t), Y_2(t)).
   *                       The list is given in ascending order.
   *                       Note that the values are not bounded to [0,1].
   * \return Do the two curves overlap (if they do, s_vals is empty).
   */
  bool _intersection_params (const Polynomial_1& polyX_1, const Integer& normX_1,
                             const Polynomial_1& polyY_1, const Integer& normY_1,
                             const Polynomial_1& polyX_2, const Integer& normX_2,
                             const Polynomial_1& polyY_2, const Integer& normY_2,
                             Parameter_list& s_vals) const;

  /*!
   * Compute all s-parameter values of the self intersection of (X(s), Y(s))
   * with (X(t), Y(t)) for some t.
   * \param polyX The coefficients of X.
   * \param polyY The coefficients of Y.
   * \param s_vals Output: A list of s-values of the self-intersection points.
   *                       Note that the values are not bounded to [0,1].
   */
  void _self_intersection_params (const Polynomial_1& polyX,
                                  const Polynomial_1& polyY,
                                  Parameter_list& s_vals) const;

  /*!
   * Compute the resultant of two bivariate polynomials in x and y with
   * respect to y. The bivariate polynomials are given as vectors of
   * polynomials, where bp1[i] is a coefficient of y^i, which is in turn a
   * polynomial in x.
   * \param bp1 The first bivariate polynomial.
   * \param bp2 The second bivariate polynomial.
   * \return The resultant polynomial (a polynomial in x).
   */

  Polynomial_1 _compute_resultant (const std::vector<Polynomial_1>& bp1,
                                 const std::vector<Polynomial_1>& bp2) const;
};

// ---------------------------------------------------------------------------
// Get the vertical tangency parameters of a given curve (X(t), Y(t)).
//
template<class AlgebraicKernel_d_1_>
const typename _Bezier_cache<AlgebraicKernel_d_1_>::Vertical_tangency_list&
_Bezier_cache<AlgebraicKernel_d_1_>::get_vertical_tangencies
        (const Curve_id& id,
         const Polynomial_1& polyX, const Integer& /* normX */)
{
  // Try to find the curve ID in the map.
  Vert_tang_map_iterator    map_iter = vert_tang_map.find (id);

  if (map_iter != vert_tang_map.end())
  {
    // Found in the map: return the cached parameters' list.
    return (map_iter->second);
  }

  // We need to compute the vertical tangency parameters: find all t-values
  // such that X'(t) = 0, and store them in the cache.
  Vertical_tangency_list&  vert_tang_list = vert_tang_map[id];


  typename Polynomial_traits_1::Differentiate	differentiate;
  const Polynomial_1&		polyX_der = differentiate (polyX);
  //const Polynomial_1&        polyX_der = nt_traits.derive (polyX);
  
  typename Algebraic_kernel_d_1::Solve_1 solve_1 = algebraic_kernel_d_1.solve_1_object();

  /* The solve_1 returns a list of std::pair
   * of type Parameter_pair. However, we need only 
   * the roots without the mults.
   */
  typename std::list<Parameter_pair> root_pair;
  solve_1 (polyX_der, std::back_inserter(root_pair));
  //nt_traits.compute_polynomial_roots (polyX_der, 0, 1, std::back_inserter(vert_tang_list));


  // Copy the roots back
  typename std::list<Parameter_pair>::iterator rpit = root_pair.begin();
  typename std::back_insert_iterator<Vertical_tangency_list> oit = std::back_inserter(vert_tang_list);

  while(rpit != root_pair.end())
  {
	*oit++ = rpit->first;
	rpit++;
  }

  return (vert_tang_list);
}

// ---------------------------------------------------------------------------
// Get the intersection parameters of two given curve (X_1(s), Y_1(s))
// and (X_2(t), Y_2(t)).
//
template<class AlgebraicKernel_d_1_>
const typename _Bezier_cache<AlgebraicKernel_d_1_>::Intersection_list&
_Bezier_cache<AlgebraicKernel_d_1_>::get_intersections
        (const Curve_id& id1,
         const Polynomial_1& polyX_1, const Integer& normX_1,
         const Polynomial_1& polyY_1, const Integer& normY_1,
         const Curve_id& id2,
         const Polynomial_1& polyX_2, const Integer& normX_2,
         const Polynomial_1& polyY_2, const Integer& normY_2,
         bool& do_ovlp)
{
  CGAL_precondition (id1 <= id2);

  // Construct the pair of curve IDs, and try to find it in the map.
  Curve_pair                curve_pair (id1, id2);
  Intersect_map_iterator    map_iter = intersect_map.find (curve_pair);

  if (map_iter != intersect_map.end())
  {
    // Found in the map: return the cached information.
    do_ovlp = map_iter->second.second;
    return (map_iter->second.first);
  }

  // We need to compute the intersection-parameter pairs.
  Intersection_info&      info = intersect_map[curve_pair];


  typename Algebraic_kernel_d_1::Compare_1 compare_1 = algebraic_kernel_d_1.compare_1_object();
  typename Algebraic_kernel_d_1::Construct_algebraic_real_1 construct_algebraic_real_1 = algebraic_kernel_d_1.construct_algebraic_real_1_object();

  // Check if we have to compute a self intersection (a special case),
  // or a regular intersection between two curves.
  if (id1 == id2)
  {
    // Compute all parameter values that lead to a self intersection.
    Parameter_list          s_vals;

    _self_intersection_params (polyX_1, polyY_1,
                               s_vals);

    // Match pairs of parameter values that correspond to the same point.
    // Note that we make sure that both parameter pairs are in the range [0, 1]
    // (if not, the self-intersection point is imaginary).

    typename Parameter_list::iterator  s_it;
    typename Parameter_list::iterator  t_it;

    const Algebraic_real_1&            one = construct_algebraic_real_1 (1);
    // const Algebraic_real_1&		one(1);

    // Use the construct functor for algebraic_real_1
    const Algebraic_real_1&            denX = construct_algebraic_real_1 (normX_1);
    const Algebraic_real_1&            denY = construct_algebraic_real_1 (normY_1);
    //const Algebraic_real_1&            denX = nt_traits.convert (normX_1);
    //const Algebraic_real_1&            denY = nt_traits.convert (normY_1);
    Point_list                         pts1;

    for (s_it = s_vals.begin(); s_it != s_vals.end(); ++s_it)
    {
      if (CGAL::sign (*s_it) == NEGATIVE)
        continue;

      // if (CGAL::compare (*s_it, one) == LARGER)
      if (compare_1 (*s_it, one) == LARGER)
        break;

      typename Polynomial_traits_1::Evaluate	evaluate;

// TODO Eliminate the CGAL::to_double for the "evaluate" function and the division later
      // Evaluate returns the Polynomial_traits_1::Coefficient
      const Algebraic_real_1&  x = construct_algebraic_real_1( evaluate (polyX_1, CGAL::to_double(*s_it)) );
      const Algebraic_real_1&  y = construct_algebraic_real_1( evaluate (polyY_1, CGAL::to_double(*s_it)) );
      //const Algebraic_real_1&  x = nt_traits.evaluate_at (polyX_1, *s_it);
      //const Algebraic_real_1&  y = nt_traits.evaluate_at (polyY_1, *s_it);

      for (t_it = s_it; t_it != s_vals.end(); ++t_it)
      {
        //if (CGAL::compare (*t_it, one) == LARGER)
        if (compare_1 (*t_it, one) == LARGER)
          break;

        //if (CGAL::compare (nt_traits.evaluate_at (polyX_1, *t_it),
        //                   x) == EQUAL &&
        //    CGAL::compare (nt_traits.evaluate_at (polyY_1, *t_it),
        //                   y) == EQUAL)

	const Algebraic_real_1& tempX = construct_algebraic_real_1( evaluate (polyX_1, CGAL::to_double(*t_it)) );
	const Algebraic_real_1& tempY = construct_algebraic_real_1( evaluate (polyY_1, CGAL::to_double(*t_it)) );

        if (compare_1 (tempX, x) == EQUAL &&
            compare_1 (tempY, y) == EQUAL)
        {
          info.first.push_back (Intersection_point_2 (*s_it, *t_it, CGAL::to_double(x) / CGAL::to_double(denX), 
								CGAL::to_double(y) / CGAL::to_double(denY)));
        }
      }
    }

    info.second = false;
    return (info.first);
  }

  // Compute s-values and t-values such that (X_1(s), Y_1(s)) and
  // (X_2(t), Y_2(t)) are the intersection points.
  Parameter_list          s_vals;

  do_ovlp = _intersection_params (polyX_1, normX_1, polyY_1, normY_1,
                                  polyX_2, normX_2, polyY_2, normY_2,
                                  s_vals);

  if (do_ovlp)
  {
    // Update the cache and return an empty list of intersection parameters.
    info.second = true;
    return (info.first);
  }

  Parameter_list          t_vals;

  do_ovlp = _intersection_params (polyX_2, normX_2, polyY_2, normY_2,
                                  polyX_1, normX_1, polyY_1, normY_1,
                                  t_vals);

  CGAL_assertion (! do_ovlp);

  // s-values and t-values reported lie in the range [0,1]. We have
  // to pair them together. Note that the numbers of s-values and t-values
  //can be different as a s-value in [0,1] may correspond to a t-values
  //outside this range. To make the pairing correct, we choose the one with 
  //less values and look amongst other values the correct one (thus the swapt below).

  // Compute all points on (X_1, Y_1) that match an s-value from the list.
  typename Parameter_list::iterator  s_it;


  const Algebraic_real_1&                   denX_1 = construct_algebraic_real_1 (normX_1);
  const Algebraic_real_1&                   denY_1 = construct_algebraic_real_1 (normY_1);
  //const Algebraic_real_1&                   denX_1 = nt_traits.convert (normX_1);
  //const Algebraic_real_1&                   denY_1 = nt_traits.convert (normY_1);

  Point_list                         pts1;


  for (s_it = s_vals.begin(); s_it != s_vals.end(); ++s_it)
  {
//TODO Need to remove the calls to CGAL::to_double
    const Algebraic_real_1&  x = CGAL::to_double(evaluate (polyX_1, CGAL::to_double(*s_it))) / CGAL::to_double(denX_1);
    const Algebraic_real_1&  y = CGAL::to_double(evaluate (polyY_1, CGAL::to_double(*s_it))) / CGAL::to_double(denY_1);
    //const Algebraic_real_1&  x = nt_traits.evaluate_at (polyX_1, *s_it) / denX_1;
    //const Algebraic_real_1&  y = nt_traits.evaluate_at (polyY_1, *s_it) / denY_1;

    pts1.push_back (My_point_2 (s_it, x, y));
  }

  // Compute all points on (X_2, Y_2) that match a t-value from the list.
  typename Parameter_list::iterator  t_it;
  const Algebraic_real_1&                   denX_2 = construct_algebraic_real_1 (normX_2);
  const Algebraic_real_1&                   denY_2 = construct_algebraic_real_1 (normY_2);
  //const Algebraic_real_1&                   denX_2 = nt_traits.convert (normX_2);
  //const Algebraic_real_1&                   denY_2 = nt_traits.convert (normY_2);

  Point_list                         pts2;

  for (t_it = t_vals.begin(); t_it != t_vals.end(); ++t_it)
  {
//CHANGED
    const Algebraic_real_1&  x = CGAL::to_double(evaluate (polyX_2, CGAL::to_double(*t_it))) / CGAL::to_double(denX_2);
    const Algebraic_real_1&  y = CGAL::to_double(evaluate (polyY_2, CGAL::to_double(*t_it))) / CGAL::to_double(denY_2);
    //const Algebraic_real_1&  x = nt_traits.evaluate_at (polyX_2, *t_it) / denX_2;
    //const Algebraic_real_1&  y = nt_traits.evaluate_at (polyY_2, *t_it) / denY_2;

    pts2.push_back (My_point_2 (t_it, x, y));
  }

  typename Polynomial_traits_1::Degree	degree;
  // Go over the points in the pts1 list.
  const bool                x2_simpler = degree(polyX_2) < degree(polyX_1);
  const bool                y2_simpler = degree(polyY_2) < degree(polyY_1);
  //const bool                x2_simpler = nt_traits.degree(polyX_2) <
  //                                       nt_traits.degree(polyX_1);
  //const bool                y2_simpler = nt_traits.degree(polyY_2) <
  //                                       nt_traits.degree(polyY_1);

  Point_iter                pit1;
  Point_iter                pit2;
  double                    dx, dy;
  Algebraic_real_1          s, t;

  const Algebraic_real_1	one = construct_algebraic_real_1 (1);
  //const Algebraic_real_1    one (1);
  
  unsigned int              k;

  //pointers are used to set the list pts1_ptr as the one with the less values
  Point_list* pts1_ptr=&pts1;
  Point_list* pts2_ptr=&pts2;
  bool swapt=pts1.size() > pts2.size();
  if (swapt)
    std::swap(pts1_ptr,pts2_ptr);
  
  for (pit1 = pts1_ptr->begin(); pit1 != pts1_ptr->end(); ++pit1)
  {
    // Construct a vector of distances from the current point to all other
    // points in the pts2 list.
    const int                     n_pts2 = pts2_ptr->size();
    std::vector<Distance_iter>    dist_vec (n_pts2);

    for (k = 0, pit2 = pts2_ptr->begin(); pit2 != pts2_ptr->end(); k++, ++pit2)
    {
      // Compute the approximate distance between the two current points.
      dx = pit1->app_x - pit2->app_x;
      dy = pit1->app_y - pit2->app_y;

      dist_vec[k] = Distance_iter (dx*dx + dy*dy, pit2);
    }
      
    // Sort the vector according to the distances from *pit1.
    std::sort (dist_vec.begin(), dist_vec.end(), 
               Less_distance_iter());

    // Go over the vector entries, starting from the most distant from *pit1
    // to the closest and eliminate pairs of points (we expect that
    // eliminating the distant points is done easily). We stop when we find
    // a pair for *pit1 or when we are left with a single point.
    bool                    found = false;

    // Returns an Algebraic_real_1
    const Algebraic_real_1& s = pit1->parameter();

    Algebraic_real_1        t;

    for (k = n_pts2 - 1; !found && k > 0; k--)
    {
      pit2 = dist_vec[k].second;

      if (pit1->equals (*pit2))
      {
        // Obtain the parameter value, and try to simplify the representation
        // of the intersection point.
        t = pit2->parameter();

        if (x2_simpler)
          pit1->x = pit2->x;
        if (y2_simpler)
          pit1->y = pit2->y;

        // Remove this point from pts2, as we found a match for it.
        pts2_ptr->erase (pit2);
        found = true;
      }
    }

    if (! found)
    {
      // We are left with a single point - pair it with *pit1.
      pit2 = dist_vec[0].second;

      // Obtain the parameter value, and try to simplify the representation
      // of the intersection point.
      t = pit2->parameter();

      if (x2_simpler)
        pit1->x = pit2->x;
      if (y2_simpler)
        pit1->y = pit2->y;

      // Remove this point from pts2, as we found a match for it.
      pts2_ptr->erase (pit2);
    }


    // Check that s- and t-values both lie in the legal range of [0,1].
    CGAL_assertion(CGAL::sign (s) != NEGATIVE && compare_1 (s, one) != LARGER &&
                   CGAL::sign (t) != NEGATIVE && compare_1 (t, one) != LARGER);
    //CGAL_assertion(CGAL::sign (s) != NEGATIVE && CGAL::compare (s, one) != LARGER &&
    //               CGAL::sign (t) != NEGATIVE && CGAL::compare (t, one) != LARGER);
    
    if (!swapt)
      info.first.push_back (Intersection_point_2 (s, t,pit1->x, pit1->y));
    else
      info.first.push_back (Intersection_point_2 (t, s,pit1->x, pit1->y));
  }

  info.second = false;
  return (info.first);
}

// ---------------------------------------------------------------------------
// Mark two curves as overlapping.
//
template<class Algebraic_kernel_d_1_>
void _Bezier_cache<Algebraic_kernel_d_1_>::mark_as_overlapping (const Curve_id& id1,
                                                   const Curve_id& id2)
{
  CGAL_precondition (id1 < id2);

  // Construct the pair of curve IDs, and try to find it in the map.
  Curve_pair                curve_pair (id1, id2);
  Intersect_map_iterator    map_iter = intersect_map.find (curve_pair);

  if (map_iter != intersect_map.end())
  {
    // Found in the map: Make sure the curves are marked as overlapping.
    CGAL_assertion (map_iter->second.second);
    return;
  }

  // Add a new entry and mark the curves as overlapping.
  Intersection_info&      info = intersect_map[curve_pair];

  info.second = true;
  return;
}

// ---------------------------------------------------------------------------
// Compute all s-parameter values of the intersection of (X_1(s), Y_1(s))
// with (X_2(t), Y_2(t)) for some t.
//
template<class Algebraic_kernel_d_1>
bool _Bezier_cache<Algebraic_kernel_d_1>::_intersection_params
        (const Polynomial_1& polyX_1, const Integer& normX_1,
         const Polynomial_1& polyY_1, const Integer& normY_1,
         const Polynomial_1& polyX_2, const Integer& normX_2,
         const Polynomial_1& polyY_2, const Integer& normY_2,
         Parameter_list& s_vals) const
{
  // Clear the output parameter list.
  if (! s_vals.empty())
    s_vals.clear();

  // We are looking for the s-values s_0 such that (X_1(s_0), Y_1(s_0)) is
  // an intersection with some point on (X_2(t), Y_2(t)). We therefore have
  // to solve the system of bivariate polynomials in s and t:
  //    I: X_2(t) - X_1(s) = 0
  //   II: Y_2(t) - Y_1(s) = 0
  //
  Integer                  coeff;
  int                      k;

  // Define types for the various functions required of the polynomials
  //typedef Polynomial_1::Construct_polynomial	construct_polynomial;
  typename Polynomial_traits_1::Construct_polynomial	construct_polynomial;
  typename Polynomial_traits_1::Degree			degree;
  typename Polynomial_traits_1::Scale			scale;
  typename Polynomial_traits_1::Get_coefficient		get_coefficient;


  // Construct the bivariate polynomial that corresponds to Equation I.
  // Note that we represent a bivariate polynomial as a vector of univariate
  // polynomials, whose i'th entry corresponds to the coefficient of t^i,
  // which is in turn a polynomial it s.

  const int	degX_2 = degree (polyX_2);
  //const int                degX_2 = nt_traits.degree (polyX_2);
  std::vector<Polynomial_1>  coeffsX_st (degX_2 < 0 ? 1 : (degX_2 + 1));

  for (k = degX_2; k >= 0; k--)
  {
    coeff = get_coefficient (polyX_2, k) * normX_1;
    coeffsX_st[k] = construct_polynomial (coeff);
    //coeff = nt_traits.get_coefficient (polyX_2, k) * normX_1;
    //coeffsX_st[k] = nt_traits.construct_polynomial (&coeff, 0);
  }
  coeffsX_st[0] = coeffsX_st[0] - scale (polyX_1, normX_2);
  //coeffsX_st[0] = coeffsX_st[0] - nt_traits.scale (polyX_1, normX_2);


  // Construct the bivariate polynomial that corresponds to Equation II.
  const int                degY_2 = degree (polyY_2);
  //const int                degY_2 = nt_traits.degree (polyY_2);
  std::vector<Polynomial_1>  coeffsY_st (degY_2 < 0 ? 1 : (degY_2 + 1));
    
  for (k = degY_2; k >= 0; k--)
  {
    coeff = get_coefficient (polyY_2, k) * normY_1;
    coeffsY_st[k] = construct_polynomial (coeff);
    //coeff = nt_traits.get_coefficient (polyY_2, k) * normY_1;
    //coeffsY_st[k] = nt_traits.construct_polynomial (&coeff, 0);
  }
  coeffsY_st[0] = coeffsY_st[0] - scale (polyY_1, normY_2);
  //coeffsY_st[0] = coeffsY_st[0] - nt_traits.scale (polyY_1, normY_2);

  // Compute the resultant of the two bivariate polynomials and obtain
  // a polynomial in s.
  Polynomial_1            res = _compute_resultant (coeffsX_st, coeffsY_st);

  //if (nt_traits.degree (res) < 0)
  if (degree (res) < 0)
  {
    // If the resultant is identically zero, then the two curves overlap.
    return (true);
  }

  // Compute the roots of the resultant polynomial and mark that the curves do
  // not overlap. The roots we are interested in must be in the interval [0,1].
  typename Algebraic_kernel_d_1::Solve_1 solve_1 = algebraic_kernel_d_1.solve_1_object();


  /* The solve_1 returns a list of std::pair
   * of type Parameter_pair. However, we need only 
   * the roots without the mults.
   */
  typename std::list<Parameter_pair> root_pair;
  solve_1 (res, std::back_inserter(root_pair));
  //nt_traits.compute_polynomial_roots (polyX_der, 0, 1, std::back_inserter(vert_tang_list));


  // Copy the roots back
  typename std::list<Parameter_pair>::iterator rpit = root_pair.begin();
  typename std::back_insert_iterator<Parameter_list> oit = std::back_inserter(s_vals);

  while(rpit != root_pair.end())
  {
	*oit++ = rpit->first;
	rpit++;
  }

  //solve_1 (res, std::back_inserter(s_vals));
  //nt_traits.compute_polynomial_roots (res,0,1,std::back_inserter (s_vals));
  return (false);

}

// ---------------------------------------------------------------------------
// Compute all s-parameter values of the self intersection of (X(s), Y(s))
// with (X(t), Y(t)) for some t.
//
template<class Algebraic_kernel_d_1_>
void _Bezier_cache<Algebraic_kernel_d_1_>::_self_intersection_params
        (const Polynomial_1& polyX,
         const Polynomial_1& polyY,
         Parameter_list& s_vals) const
{
  // Clear the output parameter list.
  if (! s_vals.empty())
    s_vals.clear();

  // We are looking for the solutions of the following system of bivariate
  // polynomials in s and t:
  //    I: X(t) - X(s) / (t - s) = 0
  //   II: Y(t) - Y(s) / (t - s) = 0
  //
  Integer                 *coeffs;
  int                      i, k;

  // Construct the bivariate polynomial that corresponds to Equation I.
  // Note that we represent a bivariate polynomial as a vector of univariate
  // polynomials, whose i'th entry corresponds to the coefficient of t^i,
  // which is in turn a polynomial it s.

  // Types for the polynomial functions
  typename Polynomial_traits_1::Construct_polynomial	construct_polynomial;
  typename Polynomial_traits_1::Degree			degree;
  typename Polynomial_traits_1::Get_coefficient		get_coefficient;


  const int                degX = degree (polyX);
  //const int                degX = nt_traits.degree (polyX);

  CGAL_assertion(degX > 0);
  if (degX <= 0) return; //no self intersection if X is constant
  std::vector<Polynomial_1>  coeffsX_st (degX);

  coeffs = new Integer [degX];

  for (i = 0; i < degX; i++)
  {
    for (k = i + 1; k < degX; k++)
      coeffs[k - i - 1] = get_coefficient (polyX, k);
      //coeffs[k - i - 1] = nt_traits.get_coefficient (polyX, k);

    //coeffsX_st[i] = nt_traits.construct_polynomial (coeffs, degX - i - 1);
    coeffsX_st[i] = construct_polynomial (coeffs, coeffs + degX - i - 1);
  }

  delete[] coeffs;

  // Construct the bivariate polynomial that corresponds to Equation II.
  const int                degY = degree (polyY);
  //const int                degY = nt_traits.degree (polyY);

  CGAL_assertion(degY > 0);
  if (degY <= 0) return; //no self intersection if Y is constant  
  std::vector<Polynomial_1>  coeffsY_st (degY);
    
  coeffs = new Integer [degY];

  for (i = 0; i < degY; i++)
  {
    for (k = i + 1; k < degY; k++)
      coeffs[k - i - 1] = get_coefficient (polyY, k);
      //coeffs[k - i - 1] = nt_traits.get_coefficient (polyY, k);

    coeffsY_st[i] = construct_polynomial (coeffs, coeffs + degY - i - 1);
    //coeffsY_st[i] = nt_traits.construct_polynomial (coeffs, degY - i - 1);
  }

  delete[] coeffs;

  // Compute the resultant of the two bivariate polynomials and obtain
  // a polynomial in s.
  Polynomial_1            res = _compute_resultant (coeffsX_st, coeffsY_st);

  //if (nt_traits.degree (res) < 0)
  if (degree (res) < 0)
      return;

  // Compute the roots of the resultant polynomial.
  //Algebraic_kernel_d_1		ak;
  typename Algebraic_kernel_d_1::Solve_1 solve_1 = algebraic_kernel_d_1.solve_1_object();


  /* The solve_1 returns a list of std::pair
   * of type Parameter_pair. However, we need only 
   * the roots without the mults.
   */
  typename std::list<Parameter_pair> root_pair;
  solve_1 (res, std::back_inserter(root_pair));


  // Copy the roots back
  typename std::list<Parameter_pair>::iterator rpit = root_pair.begin();
  typename std::back_insert_iterator<Parameter_list> oit = std::back_inserter(s_vals);

  while(rpit != root_pair.end())
  {
	*oit++ = rpit->first;
	rpit++;
  }

  //solve_1 (res, std::back_inserter(s_vals));
  //nt_traits.compute_polynomial_roots (res, std::back_inserter (s_vals));
  return;
}

// ---------------------------------------------------------------------------
// Compute the resultant of two bivariate polynomials.
//
template<class Algebraic_kernel_d_1_>
typename _Bezier_cache<Algebraic_kernel_d_1_>::Polynomial_1
_Bezier_cache<Algebraic_kernel_d_1_>::_compute_resultant
        (const std::vector<Polynomial_1>& bp1,
         const std::vector<Polynomial_1>& bp2) const
{

/*
  // Create the a bivariate polynomial from the input vector polynomials 
  // and subsequently the Sylvester matrix of polynomial coefficients. We 
  // compute the resultant from the Polynomial_traits_d::Resultant
  //
*/

  typedef typename CGAL::Polynomial_type_generator<Coefficient, 2>::Type 	Polynomial_2;
  typedef typename CGAL::Polynomial_traits_d<Polynomial_2> 			Polynomial_traits_2;

  typename Polynomial_traits_2::Construct_polynomial 				construct_polynomial_2;
  typename Polynomial_traits_2::Resultant 					resultant;

  return resultant ( construct_polynomial_2(bp1.begin(), bp1.end()), construct_polynomial_2(bp2.begin(), bp2.end()) );

}

} //namespace CGAL

#endif
