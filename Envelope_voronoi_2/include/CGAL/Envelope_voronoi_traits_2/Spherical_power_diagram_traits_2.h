// Copyright (c) 2005-2007 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id:  $
//
//
// Author(s): Ophir Setter          <ophir.setter@post.tau.ac.il>
//

/*! \file Spherical_power_diagram_traits.h
  This file contains an implementation of traits class for a power
  diagram of circles on the unit spheres. The circles are created by
  intersection the unit sphere with a rational plane.
*/

#ifndef SPHERICAL_POWER_DIAGRAM_2_H
#define SPHERICAL_POWER_DIAGRAM_2_H

#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Envelope_voronoi_2/envelope_voronoi_assertions.h>

namespace CGAL {

template <class T_Kernel>
class Spherical_power_diagram_traits_2 :
    public Arr_geodesic_arc_on_sphere_traits_2<T_Kernel>
{
public:
  typedef T_Kernel                                             Kernel;
  typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel>          Base;
  typedef Spherical_power_diagram_traits_2<Kernel>             Self;

  typedef typename Kernel::Vector_3                            Vector_3;
  typedef typename Kernel::Direction_3                         Direction_3;
  typedef typename Kernel::Plane_3                             Plane_3;
  typedef typename Kernel::Line_3                              Line_3;
  typedef typename Kernel::Point_3                             Point_3;
  typedef typename Kernel::Ray_3                               Ray_3;
  typedef typename Kernel::Direction_2                         Direction_2;

  typedef typename Base::Point_2                               Point_2;
  typedef typename Base::X_monotone_curve_2
    X_monotone_curve_2;
  typedef typename Base::Curve_2                               Curve_2;

  /*! \todo Get rid of the local var. 'ker'. Pass it through the constructor
   * instead!
   */
  class Site_2 : public Plane_3 {
  public:
    Site_2() {}

    /*! The site represents a circle on the sphere that is contined in the
     * given plane.
     */
    Site_2(const Plane_3& p) : Plane_3(p) {
      CGAL_envelope_voronoi_precondition_code(Kernel ker;);
      CGAL_envelope_voronoi_precondition_code
        (Point_3 o = ker.construct_point_3_object()(CGAL::ORIGIN););
      CGAL_envelope_voronoi_precondition
        (ker.compute_squared_distance_3_object()(o, p) <= 1);
    }

    /*! The site represents a point on the sphere (power circle of radius 0).
     */
    Site_2(const typename Kernel::Point_3& p) {
      Kernel ker;
      Point_3 o = ker.construct_point_3_object()(CGAL::ORIGIN);

      CGAL_envelope_voronoi_precondition
        (ker.compute_squared_distance_3_object()(o, p) <= 1);

      Vector_3 v = ker.construct_vector_3_object()(o, p);
      Direction_3 d = ker.construct_direction_3_object()(v);
      *this = ker.construct_plane_3_object()(p, d);
    }
  };

protected:

  //! The function return the bisector plane of the two sites.
  /*! The function returns the plance that creates the great-circle which
      is the bisector of the two power circle sites.
      The function assumes that the sites are not equal. If the two planes
      are parallel, then the bisector is another parallel plane that passes
      through the origin.
      \param s1 The first site.
      \param s2 The second site.
      \return The plane that creates the bisector between the two sites.
  */
  static Plane_3 bisector_plane(const Site_2& s1, const Site_2& s2) {
    Kernel k;
    CGAL_envelope_voronoi_precondition(k.equal_3_object() (s1, s2) == false);
    typename Kernel::Intersect_3 intersect = k.intersect_3_object();
    const Plane_3& pl1 = s1;
    const Plane_3& pl2 = s2;
    // Object obj = intersect(pl1, pl2);
    typename boost::result_of<typename Kernel::Intersect_3(Plane_3,
                                                           Plane_3)>::type
      result = intersect(pl1, pl2);

    // const Line_3* pl = object_cast<Line_3>(&obj);
    // if (pl != NULL)
    if (result) {
      const Line_3* pl =  std::get_if<Line_3>(&*result);

      // intersection should be empty - don't support s1 and s2 that are same.
      CGAL_envelope_voronoi_assertion(pl);

      // The intersection of the two generating plane is a 3D line.
      Point_3 o = k.construct_point_3_object()(ORIGIN);

      CGAL_envelope_voronoi_assertion_code
        (const Point_3& p1 = k.construct_point_on_3_object()(*pl, 0));
      CGAL_envelope_voronoi_assertion_code
        (const Point_3& p2 = k.construct_point_on_3_object()(*pl, 1));
      CGAL_envelope_voronoi_assertion
        (k.coplanar_orientation_3_object()(o, p1, p2) != COLLINEAR);
      return k.construct_plane_3_object()(*pl, o);
    }

    // intersection should be empty - don't support s1 and s2 that are same.
    // CGAL_envelope_voronoi_assertion(obj.is_empty());

    // planes are parallel, return a parallel plane through the origin.
    Point_3 o = k.construct_point_3_object()(ORIGIN);
    Vector_3 orth_v = k.construct_orthogonal_vector_3_object()(s1);
    return k.construct_plane_3_object()(o, orth_v);
  }

public:
  class Compare_distance_at_point_2
  {

    //! The function compare the distance between a point and two sites.
    /*! The function compares the distance between a point on the sphere and
      two circle sites. Currently, the point in the arrangement traits is
      represented as an extended Kernel::Direction_3.
      \return
    */
  public:
    Comparison_result operator()(const Site_2& h1,
                                 const Site_2& h2,
                                 const Point_2& p) const
    {
      // TODO: this function should intersect rays and planes.
      //       This way we avoid the are_strictly along a line
      //       (which is probably more expensive).

      typedef typename Kernel::FT                          FT;

      Kernel kernel;
      const Plane_3& pl1 = h1;
      const Plane_3& pl2 = h2;
      if (kernel.equal_3_object() (pl1, pl2) == true) return EQUAL;

      Point_3 o = kernel.construct_point_3_object() (ORIGIN);
      Line_3 line = kernel.construct_line_3_object() (o, p);
      Ray_3 ray = kernel.construct_ray_3_object()(o, p);
      Point_3 point = kernel.construct_point_on_3_object()(ray, 1);

      typename Kernel::Intersect_3 intersect_3 = kernel.intersect_3_object();
      typename Kernel::Compute_squared_distance_3 dist_3 =
        kernel.compute_squared_distance_3_object();
      typename Kernel::Are_strictly_ordered_along_line_3 are_along_line =
        kernel.are_strictly_ordered_along_line_3_object();

      // Object obj1 = intersect_3(line, pl1);
      // const Point_3* p1 = object_cast<Point_3>(&obj1);
      typename
        boost::result_of<typename Kernel::Intersect_3(Line_3, Plane_3)>::type
        result1 = intersect_3(line, pl1);

      // Object obj2 = intersect_3(line, pl2);
      // const Point_3* p2 = object_cast<Point_3>(&obj2);
      typename
        boost::result_of<typename Kernel::Intersect_3(Line_3, Plane_3)>::type
        result2 = intersect_3(line, pl2);

      const Point_3* p1 = nullptr;
      if (result1) {
        p1 = std::get_if<Point_3>(&*result1);
        CGAL_envelope_voronoi_assertion(p1);
      }

      const Point_3* p2 = nullptr;
      if (result2) {
        p2 =  std::get_if<Point_3>(&*result2);
        CGAL_envelope_voronoi_assertion(p2);
      }

      if (!p1 && !p2) return EQUAL;
      if (!p1) return LARGER;
      if (!p2) return SMALLER;

      // is_opp_n tells if point pn is on the other side of the origin.
      bool is_opp_1 = are_along_line(point, o, *p1);
      bool is_opp_2 = are_along_line(point, o, *p2);
      if (is_opp_1 && !is_opp_2) return LARGER;
      if (!is_opp_1 && is_opp_2) return SMALLER;

      CGAL_envelope_voronoi_assertion(is_opp_1 == is_opp_2);
      FT dist1 = dist_3 (o, *p1);
      FT dist2 = dist_3 (o, *p2);

      return is_opp_1 ? opposite(compare(dist1, dist2)) :
        compare(dist1, dist2);
    }
  };

  //!
  Compare_distance_at_point_2 compare_distance_at_point_2_object() const
  { return Compare_distance_at_point_2(); }

  //!
  class Construct_point_on_x_monotone_2 {
    /*! The function constructs a point on an x-monotone curve.
      \todo Maybe move this to the geodesic traits.
    */
  private:
    typedef Self Traits;

    //! The base traits (in case it has state).
    const Traits* m_traits;

  public:
    /*! Constructor
     * \param traits the traits instance
     */
    Construct_point_on_x_monotone_2(const Traits* traits) : m_traits(traits) {}

    //!
    Point_2 operator() (const X_monotone_curve_2& xcurve) const {
      const Kernel& kernel = *m_traits;
      auto cons_point_3 = kernel.construct_point_3_object();
      auto cons_ray_3 = kernel.construct_ray_3_object();
      auto cons_vector_3 = kernel.construct_vector_3_object();
      auto cons_cross = kernel.construct_cross_product_vector_3_object();
      auto cons_dir_3 = kernel.construct_direction_3_object();
      auto cons_opp_dir = kernel.construct_opposite_direction_3_object();
      auto cons_opp_vec = kernel.construct_opposite_vector_3_object();
      auto eq = kernel.equal_3_object();

      auto cons_point_2 = m_traits->construct_point_2_object();

      Point_2 s = xcurve.source();
      Point_2 t = xcurve.target();

      Point_3 o = cons_point_3(ORIGIN);

      Ray_3 ray_s = cons_ray_3(o, s);
      Vector_3 vec_s = cons_vector_3(ray_s);
      Direction_3 dir_s = cons_dir_3(vec_s);

      Ray_3 ray_t = cons_ray_3(o, t);
      Vector_3 vec_t = cons_vector_3(ray_t);
      Direction_3 dir_t = cons_dir_3(vec_t);

      Ray_3 ray_norm = cons_ray_3(o, xcurve.normal());
      Vector_3 vec_norm = cons_vector_3(ray_norm);

      if (eq(dir_s, cons_opp_dir(dir_t))) {
        Vector_3 vec_res = cons_cross(vec_t, vec_norm);
        return cons_point_2(cons_dir_3(vec_res));
      }

      Vector_3 vec_res = vec_s + vec_t;
      Vector_3 cross = cons_cross(vec_s, vec_t);
      Direction_3 dir_cross = cons_dir_3(cross);

      if (eq(dir_cross, xcurve.normal()) == false)
        vec_res = cons_opp_vec(vec_res);

      return cons_point_2(cons_dir_3(vec_res));
    }
  };

  //!
  Construct_point_on_x_monotone_2 construct_point_on_x_monotone_2_object() const
  { return Construct_point_on_x_monotone_2(this); }

  //!
  class Compare_dominance_2 {
  public:
    Comparison_result operator()(const Site_2& h1,
                                 const Site_2& h2) const
    {
      // One site dominates the other only if they are the same cite
      CGAL_envelope_voronoi_assertion_code(Kernel k;);
      CGAL_envelope_voronoi_assertion(k.equal_3_object() (h1, h2) == true);
      return EQUAL;
    }
  };

  //!
  Compare_dominance_2 compare_dominance_2_object() const
  { return Compare_dominance_2(); }

  //!
  class Compare_distance_above_2 {
  private:
    typedef Self Traits;

    //! The base traits (in case it has state).
    const Traits* m_traits;

  public:
    /*! Constructor
     * \param traits the traits instance
     */
    Compare_distance_above_2(const Traits* traits) : m_traits(traits) {}

    //!
    Comparison_result operator()(const Site_2& h1,
                                 const Site_2& h2,
                                 const X_monotone_curve_2& cv) const
    {
      typedef typename Kernel::FT                         FT;

      const Kernel& kernel = *m_traits;

      auto ortho = kernel.construct_orthogonal_vector_3_object();
      Vector_3 v1 = ortho(h1);
      Vector_3 v2 = ortho(h2);

      auto oppo = kernel.construct_opposite_vector_3_object();
      v1 = oppo(v1);
      v2 = oppo(v2);

      auto dire = kernel.construct_direction_3_object();
      auto ctp = m_traits->construct_point_2_object();
      Point_2 d1 = Point_2(ctp(dire(v1)));
      Point_2 d2 = Point_2(ctp(dire(v2)));

      Comparison_result res = m_traits->compare_y(d1, d2);
      if (res != EQUAL)
        return CGAL::opposite(res);

      // The center of the circle with the larger radius is contained inside
      // the area of dominance.
      // We where is the center of the larger circle with respect to the arc.
      // (Do not forget to check if the curve is directed right.)
      // \todo confirm that the following code is true. Voronoi diagram of
      // points on the sphere returns opposite...

      CGAL_envelope_voronoi_assertion(cv.is_vertical() == true);

      Point_3 o = kernel.construct_point_3_object()(ORIGIN);
      Plane_3 cv_plane = kernel.construct_plane_3_object()(o, cv.normal());

      FT h1_dist = kernel.compute_squared_distance_3_object()(o, h1);
      FT h2_dist = kernel.compute_squared_distance_3_object()(o, h2);
      bool h1_larger = (CGAL::compare(h1_dist, h2_dist) == SMALLER);

      typename Kernel::Construct_translated_point_3   translate =
        kernel.construct_translated_point_3_object();
      if (h1_larger)
        res = kernel.oriented_side_3_object()(cv_plane, translate(o, v1));
      else
        res = kernel.oriented_side_3_object()(cv_plane, translate(o, v2));

      if (cv.is_directed_right()) return res;
      return CGAL::opposite(res);
    }
  };

  //!
  Compare_distance_above_2 compare_distance_above_2_object() const
  { return Compare_distance_above_2(this); }

  /*! class Given to circles on a sphere (represented as planed intersecting
    the sphere) create their bisector on the unit sphere. A precondition is
    that both of the planes intersect the unit sphere.
  */
  class Construct_bisector_2 {
  private:
    typedef Self Traits;

    //! The base traits (in case it has state)
    const Traits* m_traits;

  public:
    /*! Constructor
     * \param traits the traits instance
     */
    Construct_bisector_2(const Traits* traits) : m_traits(traits) {}

    //!
    template <class OutputIterator>
      OutputIterator operator()(const Site_2& s1,
                                const Site_2& s2,
                                OutputIterator o) const
    {
      const Kernel& kernel = *m_traits;

      // check that both planes intersect the unit sphere.
      CGAL_envelope_voronoi_precondition_code (Point_3 origin = \
                              kernel.construct_point_3_object() (ORIGIN));
      CGAL_envelope_voronoi_precondition (kernel.compute_squared_distance_3_object()     \
                         (origin, s1) <= 1);
      CGAL_envelope_voronoi_precondition (kernel.compute_squared_distance_3_object() \
                         (origin, s2) <= 1);

      // check that both plane are oriented correctly.
      CGAL_envelope_voronoi_precondition_msg (                          \
        kernel.oriented_side_3_object() (s1, origin) == ON_POSITIVE_SIDE, \
        "the intersecting plane should be oriented the other way.");
      CGAL_envelope_voronoi_precondition_msg (                          \
        kernel.oriented_side_3_object() (s2, origin) == ON_POSITIVE_SIDE, \
        "the intersecting plane should be oriented the other way.");

      Plane_3 bis = bisector_plane(s1, s2);

#if defined(CGAL_FULL_X_MONOTONE_GEODESIC_ARC_ON_SPHERE_IS_SUPPORTED)
      Curve_2 c(bis);
#else
      // we use 4 different directions to make sure that there is no
      // arc bigger then 180 degrees.
      typedef typename Kernel::Vector_3                Vector_3;
      typedef typename Kernel::Direction_3             Direction_3;

      Vector_3 v1 = kernel.construct_base_vector_3_object() (bis, 1);
      Vector_3 v2 = kernel.construct_base_vector_3_object() (bis, 2);

      Direction_3 d1 = kernel.construct_direction_3_object() (v1);
      Direction_3 d2 = kernel.construct_direction_3_object() (v2);
      Direction_3 d3 = kernel.construct_opposite_direction_3_object() (d1);
      Direction_3 d4 = kernel.construct_opposite_direction_3_object() (d2);
      auto p1 = m_traits->construct_point_2_object()(d1);
      auto p2 = m_traits->construct_point_2_object()(d2);
      auto p3 = m_traits->construct_point_2_object()(d3);
      auto p4 = m_traits->construct_point_2_object()(d4);

      Curve_2 c1 = m_traits->construct_curve_2_object()(p1, p2);
      Curve_2 c2 = m_traits->construct_curve_2_object()(p2, p3);
      Curve_2 c3 = m_traits->construct_curve_2_object()(p3, p4);
      Curve_2 c4 = m_traits->construct_curve_2_object()(p4, p1);
#endif

      // Temporary!!! Initiating new traits locally.
      typedef std::list<CGAL::Object>                  Object_list;
      Object_list x_monotones;

#if defined(CGAL_FULL_X_MONOTONE_GEODESIC_ARC_ON_SPHERE_IS_SUPPORTED)
      m_traits->make_x_monotone_2_object()(c, std::back_inserter(x_monotones));
#else
      m_traits->make_x_monotone_2_object()(c1, std::back_inserter(x_monotones));
      m_traits->make_x_monotone_2_object()(c2, std::back_inserter(x_monotones));
      m_traits->make_x_monotone_2_object()(c3, std::back_inserter(x_monotones));
      m_traits->make_x_monotone_2_object()(c4, std::back_inserter(x_monotones));
#endif

      for (auto it = x_monotones.begin(); it != x_monotones.end(); ++it) {
        X_monotone_curve_2 curve;
        if (assign(curve, *it)) {
          *o++ = curve;
          continue;
        }
        // all have to be x-mono curves.
        CGAL_error_msg("All bisectors sub-curves should be x-mono curves");
      }

      return o;
    }
  };

  Construct_bisector_2 construct_bisector_2_object() const
  { return Construct_bisector_2(this); }
};

} //namespace CGAL

#endif
