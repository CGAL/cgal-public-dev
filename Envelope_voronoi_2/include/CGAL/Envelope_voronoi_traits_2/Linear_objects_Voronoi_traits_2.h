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


/*! \file Linear_objects_Voronoi_traits_2 
  A traits for constructing Voronoi diagram of linear objects.
  The traits class supports interior-disjoint points, segments, rays, and lines.
  The implementation is based on the Curved_kernel_via_analysis_2 of CGAL
*/

#ifndef CGAL_LINEAR_OBJECTS_VORONOI_TRAITS_2_H
#define CGAL_LINEAR_OBJECTS_VORONOI_TRAITS_2_H

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Algebraic_structure_traits.h>

#include <CGAL/Envelope_3/Envelope_base.h>
#include <CGAL/Envelope_voronoi_2/envelope_voronoi_assertions.h>
#include <CGAL/Envelope_voronoi_traits_2/Algebraic_apollonius_traits_base_2.h>

namespace CGAL {

//! \todo Move to Algebraic_structure_traits
/*@{*/
template <typename T>
bool is_zero(const boost::numeric::interval<T> &num)
{
  return (num == boost::numeric::interval<T>(0));
}

template <typename T>
bool is_one(const boost::numeric::interval<T> &num)
{
  return (num == boost::numeric::interval<T>(1));
}

/*@}*/

namespace Envelope_voronoi_traits_2
{
  /*!
   * \class Linear_object_site
   *        The class represents a segment voronoi site. To make life easier
   *        in computing the segment voronoi diagram we treat the end points 
   *        of the segment as sites on their own. This means that the class 
   *        should support both a segment site and a point site. The segment
   *        sites are open segment actually.
   */
  template <typename NT_>
  class Linear_object_site
  {
    //! \name Types
    /*@{*/
    typedef NT_                                                  NT;

  public:
    typedef CGAL::Simple_cartesian<NT>                           Kernel;
    typedef Linear_object_site<NT>                               Self;

    typedef typename Kernel::Point_2                             Point_2;
    typedef typename Kernel::Segment_2                           Segment_2;
    typedef typename Kernel::Ray_2                               Ray_2;
    typedef typename Kernel::Line_2                              Line_2;

  protected:
    //! The Site_type enum specifies the type of the possible sites.
    enum Site_type {Point_site, Segment_site, Ray_site, Line_site};
    /*@}*/

    //! \name Member Variables.
    /*@{*/ 
  
    //! Kernel point to hold the point in case that the site is a point.
    Point_2 m_point;

    //! Kernel segment to hold the point in case that the site is a segment.
    Segment_2 m_segment;
 
    //! Kernel ray to hold the point in case that the site is a ray.
    Ray_2 m_ray;

    //! Kernel line to hold the point in case that the site is a line.
    Line_2 m_line;
  
    //! The type of the object (point/segment/ray/line).
    Site_type m_type;
    /*@}*/

  public:
    Linear_object_site(const Point_2& point)
      : m_point(point), m_type(Point_site)
    {}

    Linear_object_site(const Segment_2& segment)
      : m_segment(segment), m_type(Segment_site)
    {}

    Linear_object_site(const Ray_2& ray)
      : m_ray(ray), m_type(Ray_site)
    {}

    Linear_object_site(const Line_2& line)
      : m_line(line), m_type(Line_site)
    {}

    bool is_point() const
    {
      return m_type == Point_site;
    }

    bool is_segment() const
    {
      return m_type == Segment_site;
    }
      
    bool is_ray() const
    {
      return m_type == Ray_site;
    }

    bool is_line() const
    {
      return m_type == Line_site;
    }
      
    const Point_2& point() const
    {
      CGAL_envelope_voronoi_precondition (is_point());
      return m_point;
    }

    Point_2& point()
    {
      CGAL_envelope_voronoi_precondition (is_point());
      return m_point;
    }

    const Segment_2& segment() const
    {
      CGAL_envelope_voronoi_precondition (is_segment());
      return m_segment;
    }

    Segment_2& segment()
    {
      CGAL_envelope_voronoi_precondition (is_segment());
      return m_segment;
    }

    const Ray_2& ray() const
    {
      CGAL_envelope_voronoi_precondition (is_ray());
      return m_ray;
    }

    Ray_2& ray()
    {
      CGAL_envelope_voronoi_precondition (is_ray());
      return m_ray;
    }
      
    const Line_2& line() const
    {
      CGAL_envelope_voronoi_precondition (is_line());
      return m_line;
    }

    Line_2& line()
    {
      CGAL_envelope_voronoi_precondition (is_line());
      return m_line;
    }

    //! Returns the supporting line of the underlying object.
    /*! Returns the supporting line of the underlying object.
      \pre The object is not a point site.  
      \return The supporting line of the site.
    */
    Line_2 supporting_line() const
    {
      CGAL_envelope_voronoi_precondition (is_point() == false);
      
      switch (type())
      {
       case Segment_site: 
        return segment().supporting_line();
       case Ray_site: 
        return ray().supporting_line();
       case Line_site: 
        return line();
       default: 
        CGAL_error();
        return Line_2();
      }
    }
  
    //! Returns a point on the site.
    /*! The function is returns a point on the site in the same way as the 
      Kernel classes work.
      For each i we return a different point. If the site has a source
      (segment/ray) then 0 returns the source point. If a site has a target
      (segment) then 1 returns the target vertex. If the site is a point then
      we return a point. In any case point_on(1) - point_on(0) creates a vector
      in the same direction as the underlying object.
      \pre i should be 0 or 1. If the site is a point, only 0 is legal.
      \post See above.
      \param i An index for the point to be returned.
      \return A point on site.
    */
    Point_2 point_on(int i) const
    {
      CGAL_envelope_voronoi_precondition (i == 0 || i == 1);
      CGAL_envelope_voronoi_precondition (is_point() == false || i == 0);
      
      switch (type())
      {
       case Point_site:
        return point();
       case Segment_site: 
        return segment().point(i);
       case Ray_site: 
        return ray().point(i);
       case Line_site: 
        return line().point(i);
       default: 
        CGAL_error();
        return Point_2();
      }
    }

    //! equality operator
    bool operator== (const Self& s) const
    {
      if (type() != s.type())
        return false;

      switch (type())
      {
       case Point_site:
        return (m_point == s.m_point);
       case Segment_site: 
        return (m_segment == s.m_segment);
       case Ray_site: 
        return (m_ray == s.m_ray);
       case Line_site: 
        return (m_line == s.m_line);
       default: 
        CGAL_error();
        return false;
      }
    }

  protected:
    Site_type type() const
    {
      return m_type;
    }
  };
}

/*!
 * \class Linear_objects_Voronoi_traits_2 
 * A traits class to calculate the Voronoi diagram of set of linear segments
 * based on the CKvA.
 * \todo Do not inherit from apllonius. There is no point.
 */  
template <class CurvedKernel_>
class Linear_objects_Voronoi_traits_2
  : public Algebraic_apollonius_traits_base_2< 
  CurvedKernel_, 
  Linear_objects_Voronoi_traits_2< CurvedKernel_ >,
  typename Envelope_voronoi_traits_2::Linear_object_site< 
    typename CurvedKernel_::Curve_kernel_2::         
    Coordinate_1::Rational
    >,
  CGAL::Field_tag
  >
{  
  typedef CurvedKernel_                                     Curved_kernel_2;

  typedef typename Curved_kernel_2::Curve_kernel_2          Algebraic_kernel_2;
  typedef typename Algebraic_kernel_2::Coefficient          Coefficient;
  typedef typename Algebraic_kernel_2::
  Coordinate_1::Rational                                    Rational;
  
  typedef typename CGAL::
  Get_arithmetic_kernel<Coefficient>::
  Arithmetic_kernel::Field_with_sqrt                        Field_with_sqrt;

  typedef boost::numeric::interval<Field_with_sqrt>         Interval;

public:

  typedef Envelope_voronoi_traits_2::
  Linear_object_site<Rational>                            Site_2;

  typedef Linear_objects_Voronoi_traits_2<Curved_kernel_2>  Self;

  typedef Algebraic_apollonius_traits_base_2< 
    Curved_kernel_2, Self, Site_2, CGAL::Field_tag >        Base;

  typedef typename Base::Point_2                            Point_2;
  typedef typename Base::X_monotone_curve_2                 X_monotone_curve_2;
  typedef typename Base::Curve_2                            Curve_2;
  typedef typename Base::Multiplicity                       Multiplicity;

  typedef typename Base::Surface_3                          Surface_3;
  typedef typename Base::Xy_monotone_surface_3            Xy_monotone_surface_3;
  

protected:
  //! Compute the distance from a point to a site.
  /*! Compute the distance from a point to a site - we actually compute the
    SQUARED distance.
    \param px The x-coordinate of the point.
    \param py The y-coordinate of the point.
    \param site The site we compute the distance to.
    \param out_excepted_at_endpoint Pointer to a result boolean. The boolean is
    True if the result is excepted at an enpoint of a ray/segment. This is 
    needed by the proximity functors (see below).
    \return The distance from the point to this site. (+ output argument).
  */
  template <typename NT>
  static NT distance(const NT &px,
                     const NT &py, const Site_2& site, 
                     bool *out_excepted_at_endpoint = NULL)
  {
    if (out_excepted_at_endpoint)
      *out_excepted_at_endpoint = false;
      
    // (qx, qy) will be the closest point to (px, py) on site.
    NT qx, qy;
    if (site.is_point() == false)
    {
      NT x1, y1, x2, y2;
      x1 = site.point_on(0).x();
      y1 = site.point_on(0).y();
      x2 = site.point_on(1).x();
      y2 = site.point_on(1).y();

      // we compute closest point on the linear object to the point (px, py).
      // u is a rational which tells us where we are on the segment: 
      // u = 0 - we are on (x1, y1), u = 1 - we are on (x2, y2), every other
      // u value we are on (x1, y1) + u *(x2 - x1, y2 - y1)
      // http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
      NT dx = x2 - x1;
      NT dy = y2 - y1;

      NT u = ((px - x1)*dx + (py - y1)*dy) / (dx*dx + dy*dy);

      // if the site is a ray we cannot have a negative u.
      // if the site is a segment we have to have a value between [0, 1].
      if (site.is_ray())
        u = CGAL::max<NT>(u, 0);
      if (site.is_segment())
      {
        u = CGAL::max<NT>(u, 0);
        u = CGAL::min<NT>(u, 1);
      }
      if (out_excepted_at_endpoint && 
          ((bool)(CGAL::is_zero(u)) || (bool)(CGAL::is_one(u))))
        *out_excepted_at_endpoint = true;
      
      qx = x1 + u*dx;
      qy = y1 + u*dy;
    }
    else
    {
      CGAL_envelope_voronoi_assertion(site.is_point());

      qx = site.point_on(0).x(); 
      qy = site.point_on(0).y();
    }

    // distance between p and q.
    NT x = qx - px;
    NT y = qy - py;
    return x*x + y*y;
  }

public:

  //! Make_xy_monotone_3 functor required by the EnvelopeTraits_3 concept.
  /*! Make_xy_monotone_3
    The traits class treats each segment as a set of 3~different Voronoi sites:
    the two vertices of the segment and the open set that is the interior of
    the segment. (A ray is treated as 2~different sites.)
    The handling of all predicates, specifically predicates related to 
    bisector construction, are made much simpler with the cost of an enlarged
    set of sites.
    The split of a site into several sites is implemented in this functor.
  */
  class Make_xy_monotone_3
  {
  public:
    template <class OutputIterator>
    OutputIterator operator() (Surface_3 S,
                               bool is_lower, 
                               OutputIterator oi) const
    {
      // In any case, S is a surface.
      *oi++ = S;

      // The endpoints (2 - in case of a segment and 1 - in the case of a ray).
      // are the other "surfaces" to consider.
      if (S.is_segment())
      {
        typename Site_2::Segment_2 &s = S.segment();
        *oi++ = Xy_monotone_surface_3(s.source());
        *oi++ = Xy_monotone_surface_3(s.target());
      }
      
      if (S.is_ray())
      {
        *oi++ = Xy_monotone_surface_3(S.ray().source());
      }

      return oi;
    }
  };

  Make_xy_monotone_3
  make_xy_monotone_3_object () const
  {
    return Make_xy_monotone_3();
  }

  //! Construct_projected_intersections_2 class constructs the bisector curve 
  //  of two sites.
  /*! After the splitting we only need to consider three types of bisectors,
    namely, the bisector between two points, the bisector between two
    lines (where a line can also relate to the interior of a segment or a
    ray) that is a two-line bisector, and the bisector between a point and
    a line that is a full parabola. 
  */
  class Construct_projected_intersections_2
  {
    typedef typename Site_2::Point_2                       KPoint_2;
    
  protected:
    const Self * m_traits;
    
  public:
    Construct_projected_intersections_2(const Self * traits)
      : m_traits(traits)
    {}
    
    template <class OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s1,
                              const Xy_monotone_surface_3& s2,
                              OutputIterator o) const
    {
      //! \todo Assertion that the sites are interior-disjoint.
        
        
      // if one of the sites is a point and the other is a linear
      // object I prefer that s1 will be the point and s2 the linear obj.
      if (s1.is_point() == false && s2.is_point())
        return operator() (s2, s1, o);
        
      // Two points: we use the kernel to get the line bisector and then 
      // construct the polynomial that represents it.
      if (s1.is_point() && s2.is_point())
      {
        if (CGAL::compare_xy(s1.point(), s2.point()) == EQUAL)
          return o;
          
        typename Site_2::Line_2 line = CGAL::bisector(s1.point(), s2.point());
          
        return construct_bisector_quadratic(0, 0, 0, 
                                            line.a(), line.b(), line.c(), o);
      }
        
      // The bisector of two linear objects is a degenerate hyperbola in
      // the general case. If the two objects are collinear then we take
      // the line bisector of their endpoints.
      if (s1.is_point() == false && s2.is_point() == false)
      {
        typename Site_2::Line_2 line1 = s1.supporting_line();
        typename Site_2::Line_2 line2 = s2.supporting_line();
          
          
        if ((line1 == line2) || (line1 == line2.opposite()))
        {
          // The sites are interior-disjoint, so they cannot overlap
          // (and more specifically, cannot be lines).
          CGAL_envelope_voronoi_assertion_msg(                        \
                                              s1.is_segment() || s1.is_ray(),                           \
                                              "The sites are not interior-disjoint.");
          CGAL_envelope_voronoi_assertion_msg(                        \
                                              s2.is_segment() || s2.is_ray(),                           \
                                              "The sites are not interior-disjoint.");
            
          typename Site_2::Point_2 s1_point = s1.point_on(0);
          typename Site_2::Point_2 s2_point = s2.point_on(0);

          if (CGAL::compare_xy(s1_point, s2_point) == SMALLER)
          {
            typename Site_2::Point_2 p = 
              s1.is_segment() ? s1.segment().max() : s1.ray().source();
              
            return
              bisector_of_collinear_point_and_linear_obj(p, s2, o);
          }
          else
          {
            typename Site_2::Point_2 p = 
              s1.is_segment() ? s1.segment().min() : s1.ray().source();
            return 
              bisector_of_collinear_point_and_linear_obj(p, s2, o);
          }
        }
        else
        {
          // The segments have different supporting lines. This is the general
          // case of a degenerate hyperbola.

          // The bisector of the two segments is like the bisector of their 2 
          // supporting lines - a degenerate hyperbola:
          // Let the equations of the lines be
          //
          // L1 : A1x + B1y + C1 = 0 and L2 : A2x + B2y + C2 = 0
          //
          // Then the equation of the hyperbola is:
          //
          // rx^2 + sy^2 + txy + ux + vy + w = 0 
          //
          // where:
          //
          // r = a1^2*b2^2 - b1^2*a2^2         s = b1^2*a2^2 - a1^2*b2^2
          //
          // t = 2*a1*b1 (a2^2 + b2^2) - 2*a2*b2 (a1^2 + b1^2)
          //
          // u = 2*c1*a1 (a2^2 + b2^2) - 2*c2*a2 (a1^2 + b1^2)
          //
          // v = 2*b1*c1 (a2^2 + b2^2) - 2*b2*c2 (a1^2 + b1^2)
          //
          // w = c1^2 (a2^2 + b2^2) - c2^2 (a1^2 + b1^2)
          //
          typename Site_2::Line_2 line1 = s1.supporting_line();
          typename Site_2::Line_2 line2 = s2.supporting_line();
            
          Rational a1 = line1.a();
          Rational b1 = line1.b();
          Rational c1 = line1.c();
            
          Rational a2 = line2.a();
          Rational b2 = line2.b();
          Rational c2 = line2.c();
            
          Rational a1_2 = a1 * a1; Rational b1_2 = b1 * b1; 
          Rational a2_2 = a2 * a2; Rational b2_2 = b2 * b2; 
            
          Rational a1_b1 = a1_2 + b1_2; Rational a2_b2 = a2_2 + b2_2;
            
          Rational r = a1_2 * b2_2 - b1_2 * a2_2;
          Rational s = b1_2 * a2_2 - a1_2 * b2_2;
          Rational t = 2*a1*b1*a2_b2 - 2*a2*b2*a1_b1;
          Rational u = 2*c1*a1*a2_b2 - 2*c2*a2*a1_b1;
          Rational v = 2*b1*c1*a2_b2 - 2*b2*c2*a1_b1;
          Rational w = c1*c1*a2_b2 - c2*c2*a1_b1;

          // construct the degenerate hyperbola.
          return construct_bisector_quadratic(r, s, t, u, v, w, o);
        }
            
      }

      // we took care of point-point and line-line cases.
      // this case is point-line which results in a parabola in the general
      // case, or just a rational line in case that the point is collinear 
      // with the line.
      CGAL_envelope_voronoi_assertion(                \
                                      s1.is_point() && s2.is_point() == false);

      if (s2.supporting_line().has_on(s1.point()))
      {
        return 
          bisector_of_collinear_point_and_linear_obj(s1.point(), 
                                                     s2, o);
      }
      else
      {
        // the point and the segment aren't collinear - the bisector is a
        // parabolic arc. If the line equation is: Ax + By + C = 0 and the
        // point is (p, q) then the geometric place of all the points with
        // equal distance from the point and the line is:
        //
        // rx^2 + sy^2 + txy + ux + vy + w = 0
        //
        // where:
        //
        // r = B^2, s = A^2, t = -2AB
        // u = -2*[p(A^2 + B^2) + CA]
        // v = -2*[q(A^2 + B^2) + CB]
        // w = (p^2 + q^2)(A^2 + B^2) - C^2
        //
        Rational a = s2.supporting_line().a();
        Rational b = s2.supporting_line().b();
        Rational c = s2.supporting_line().c();
        
        Rational p = s1.point().x(); Rational q = s1.point().y();
        
        Rational r = b*b;
        Rational s = a*a;
        Rational t = -2*a*b;
        
        Rational a_2_b_2 = a*a + b*b;
        Rational u = -2*(p*a_2_b_2 + c*a);
        Rational v = -2*(q*a_2_b_2 + c*b);
        Rational w = (p*p + q*q)*a_2_b_2 - c*c;
        
        return construct_bisector_quadratic(r, s, t, u, v, w, o);
      }
    }
    
  protected:
    
    //! construct_bisector_quadratic function construct a quadratic curve
    // to represent a bisector of two Voronoi sites.
    /*! The function constructs the quadratic of the form:
      $rx^2 + sy^2 + txy + ux + vy + w = 0$.
      It uses the traits class function to construct the conic, and then 
      converts the types to the correct type required by the concept.
      \return An output iterator whose value type is CGAL::Object which contain
      pair<X_monotone_curve_2,Multiplicity> or Point_2.
    */
    template <class OutputIterator>
    OutputIterator construct_bisector_quadratic (
                                                 const Rational& r,
                                                 const Rational& s,
                                                 const Rational& t,
                                                 const Rational& u,
                                                 const Rational& v,
                                                 const Rational& w,
                                                 OutputIterator o) const
    {
      typedef std::list<X_monotone_curve_2> Xcurve_list;

      std::list<X_monotone_curve_2> x_conix;
      m_traits->construct_conic(r, s, t, u, v, w, 
                                std::back_inserter(x_conix));

      convert_to_bisector(x_conix.begin(), x_conix.end(), o);

      return o;
    }

    //! Construct the bisector of a point and a linear object (segment, ray).
    /*! Construct the bisector of a point and a linear object (segment, ray).
      \param point The point
      \param site The linear object (segment or ray).
      \param o Output iterator
      \pre As in this whole class, the point should not lay in the interior 
      of the ray/segment.
      \return An output iterator whose value type is CGAL::Object which contain
      pair<X_monotone_curve_2,Multiplicity> or Point_2.
    */
    template <class OutputIterator>
    OutputIterator bisector_of_collinear_point_and_linear_obj (
                                                               const KPoint_2& point, 
                                                               const Site_2& site,
                                                               OutputIterator o) const
    {
      typedef std::list<X_monotone_curve_2> Xcurve_list;
        
      std::list<X_monotone_curve_2> x_conix;
      m_traits->bisector_of_collinear_point_and_linear_obj
        (point, site, std::back_inserter(x_conix));
        
      convert_to_bisector(x_conix.begin(), x_conix.end(), o);
        
      return o;
    }
    
    /*! Converts a sequence of X_monotone_curve_2 to a sequence of CGAL::Object
      which contain pair<X_monotone_curve_2,Multiplicity>.
      The multiplicity used is always 1 (this is the multiplicity of the 
      bisector in this function).
      \return An output iterator containing the sequence.
    */
    template <class InputIterator, class OutputIterator>
    OutputIterator convert_to_bisector (InputIterator begin, InputIterator end,
                                        OutputIterator o) const
    {
      while (begin != end)
      {
        *o++ = make_object(std::make_pair(*begin++,
                                          Multiplicity(1)));
      }
        
      return o;
    }

  };

  Construct_projected_intersections_2 
  construct_projected_intersections_2_object() const
  {
    return Construct_projected_intersections_2(this);
  }

  
  class Construct_projected_boundary_2
  {
  protected:
    const Self * m_traits;
    
  public:
    
    Construct_projected_boundary_2(const Self * traits)
      : m_traits(traits)
    {}
    
    //! Constructs the boundary of the distance surfaces.
    /*! The interior of a segment behaves like a full straight line inside the 
      region between the two perpendicular lines at its vertices. Moreover, 
      it has no effect on any other point outside that region as the vertices 
      of the segment dominate the rest of the Euclidean plane.
      The domain of the distance function to the interior of the segment is, 
      therefore, defined to be the aforementioned region. The domain of a ray 
      is the half-plane that contains it and bounded by the perpendicular 
      line at its source. The domain of a distance function to a point site 
      and a line site remains the whole plane.
      
      \return a sequence of CGAL::Object that wraps 
      pair<X_monotone_curve_2, Oriented_side>
    */
    template <class OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s,
                              OutputIterator o) const
    {
      // the domain of the distance functions from points and lines is
      // the whole plane.
      if (s.is_point() || s.is_line())
        return o;

      // if this is a segment/ray then the boundary of the surface is 
      // the vertical line/s from the endpoint/s.
      // We need to consturct the two lines/one line and specify which 
      // areas are above and below.
      CGAL_envelope_voronoi_assertion(s.is_segment() || s.is_ray());
      typename Site_2::Point_2 source = s.point_on(0);
      typename Site_2::Point_2 target = s.point_on(1);
        
      // The direction of the source with respect to the target determine
      // the slope of the lines which determines the value of the 
      // Oriented_side (see below, and compare_above func).
      CGAL::Oriented_side source_side = ON_NEGATIVE_SIDE;
      if (compare_above(source, target) != LARGER)
      {
        source_side = ON_POSITIVE_SIDE;
      }
      // source is not equal target
      CGAL_envelope_voronoi_assertion(source_side != ON_ORIENTED_BOUNDARY);
        
      typedef std::list<X_monotone_curve_2> Xcurve_list;
      std::list<X_monotone_curve_2> x_conix;
        
      // Both ray and segement have the boundary at the source.
      m_traits->bisector_of_collinear_point_and_linear_obj
        (source, s, std::back_inserter(x_conix));
      o = convert_to_boundary(x_conix.begin(), x_conix.end(), 
                              source_side, o);
        
      // Segment has boundary at the target too.
      if (s.is_segment())
      {
        x_conix.clear();
        m_traits->bisector_of_collinear_point_and_linear_obj
          (target, s, std::back_inserter(x_conix));
        o = convert_to_boundary(x_conix.begin(), x_conix.end(), 
                                CGAL::opposite(source_side), o);
      }
      return o;
    }

  protected:
    /*!
     * Checks which of the points is above the bisector of the two and
     * which is below the bisector.
     * \param p1 The first point.
     * \param p2 The second point.
     */
    static Comparison_result compare_above(
                                           const typename Site_2::Point_2& p1,
                                           const typename Site_2::Point_2& p2)
    {
      Comparison_result res = CGAL::compare_y(p1, p2);
      if (res != EQUAL)
        return res;
      return CGAL::opposite(CGAL::compare_x(p1, p2));
    }

    /*! Converts a sequence of X_monotone_curve_2 to a sequence of 
      CGAL::Object that wraps pair<X_monotone_curve_2, Oriented_side>.
      \param begin begin of the sequence.
      \param end end of the sequence.
      \param side the assigned Oriented_side.
      \param o The output iterator.
      \return a sequence of CGAL::Object that wraps 
      pair<X_monotone_curve_2, Oriented_side>.
    */
    template <class InputIterator, class OutputIterator>
    OutputIterator convert_to_boundary (InputIterator begin, InputIterator end,
                                        CGAL::Oriented_side side,
                                        OutputIterator o) const
    {
      while (begin != end)
      {
        *o++ = make_object(std::make_pair(*begin++, side));
      }
        
      return o;
    }

  };

  Construct_projected_boundary_2 
  construct_projected_boundary_2_object() const
  {
    return Construct_projected_boundary_2(this);
  }
  
  
  //! Proximity predicates section.
  /*! We need all the below functor becuase we need to perform a symbolic
    pertubation on some of the site. In other words, comparing the
    distances of a point in the plane to two points or two linear objects 
    is trivial (we just compare the ``regular'' squared distance from 
    that point to the sites).
    When comparing the distances to a point and to a linear object we have
    to be more careful, as the interiors of segments and rays are not a
    closed sets. We compare the squared distances to the point and to the 
    linear object. If an ``equal'' result is the outcome of considering an 
    endpoint of the linear object, we decide that the point site is closer,
    as the endpoint of a linear object (a segment or a ray) is not contained
    in its interior.
  */
  /*@{*/      
  class Compare_z_at_xy_3
  {
  protected:
    const Base * m_base;
      
  public:
        
    Compare_z_at_xy_3(const Base* base)
      : m_base(base) {}

    Comparison_result operator()(const Point_2& p,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const
    {
      Rational xl = m_base->kernel().lower_boundary_x_2_object() (p.xy());
      Rational xu = m_base->kernel().upper_boundary_x_2_object() (p.xy());
        
      Rational yl = m_base->kernel().lower_boundary_y_2_object() (p.xy());
      Rational yu = m_base->kernel().upper_boundary_y_2_object() (p.xy());
        
      Interval xi(xl, xu);
      Interval yi(yl, yu);

      // I have no idea why this try-catch is needed, but for some reason
      // on my gcc 4.2.3 an "boost::interval: uncertain comparison" exception
      // is thrown even though we picked the tribool comparison.
      // My only reaction is WTF?!?!
      try
      {
        Interval d1 = distance(xi, yi, h1);
        Interval d2 = distance(xi, yi, h2);
        
        using namespace boost::numeric::interval_lib::compare::tribool;
        std::cout << "comparing using tribool" << std::endl;
        boost::tribool res = (d1 < d2);
        std::cout << "END comparing using tribool" << std::endl;
        
        if (res)
          return SMALLER;
        if (!res)
          return LARGER;
      }
      catch(const boost::numeric::interval_lib::comparison_error &)
      {}

      Field_with_sqrt px, py;
      bool except_end1, except_end2;
      convert_to_CORE(p, px, py);
      Field_with_sqrt de1 = distance(px, py, h1, &except_end1);
      Field_with_sqrt de2 = distance(px, py, h2, &except_end2);
      CGAL::Comparison_result res = CGAL::compare(de1, de2);

      // If equal, prefer the point to be closer.
      if (res == EQUAL && except_end1 != except_end2)
      {
        if (except_end1)
          return CGAL::LARGER;
        else
          return CGAL::SMALLER;
      }
      return res;
    }
      
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const
    {
      // Take a representative from the interior of the point 
      // compare surfaces at that point 
      Point_2 p = m_base->construct_interior_vertex_2_object() (cv);
      return (*this)(p, h1, h2);
    }
      
      
    Comparison_result operator()(const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const
        
    {
      // We support only interior-disjoint objects, so the only sites that
      // are allowed to overlap are two points.
      CGAL_envelope_voronoi_assertion (h1.is_point() && h2.is_point());
      return EQUAL;
    }
      
  };
  
  Compare_z_at_xy_3 compare_z_at_xy_3_object() const
  {
    return Compare_z_at_xy_3(this);
  }

  class Compare_z_at_xy_above_3
  {
  public:
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const
    {
      std::pair<Rational, Rational> p_above = 
        Voronoi_diagram_of_lines_3::
        rational_point_above_or_below_arc<Base>(cv, true);
        
      bool except_end1, except_end2;
      Rational de1 = distance<Rational>(p_above.first, p_above.second, h1, 
                                        &except_end1);
      Rational  de2 = distance<Rational>(p_above.first, p_above.second, h2,
                                         &except_end2);
      CGAL::Comparison_result res = CGAL::compare(de1, de2);
      
      // if equal, prefer the point to be closer.
      if (res == EQUAL && except_end1 != except_end2)
      {
        if (except_end1)
          return CGAL::LARGER;
        else
          return CGAL::SMALLER;
      }
      
      return res;
    }
      
  };

  Compare_z_at_xy_above_3 compare_z_at_xy_above_3_object() const
  {
    return Compare_z_at_xy_above_3();
  }

  class Compare_z_at_xy_below_3
  {
  protected: 
    const Self *m_self;
  public:

    Compare_z_at_xy_below_3(const Self *self)
      : m_self(self) {}
      
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const
    {
      Comparison_result res = CGAL::opposite(m_self->
                                             compare_z_at_xy_above_3_object() 
                                             (cv, h1, h2));
#ifndef CGAL_NDEBUG 
      std::pair<Rational, Rational> p_below = 
        Voronoi_diagram_of_lines_3::
        rational_point_above_or_below_arc<Base>(cv, false);

      Rational de1 = distance<Rational>(p_below.first, p_below.second, h1);
      Rational de2 = distance<Rational>(p_below.first, p_below.second, h2);

      // a segment and its end point have equal on one hand and 
      // something else on the other. 
      // \todo make the assertion work in case of equal.
      //       maybe it has to do with the proximity predicates todos.
      Comparison_result res2 = CGAL::compare(de1, de2);
      if (res != EQUAL && res2 != EQUAL)
        assert(res == res2);
#endif
        
      return res;
    }

  };

  Compare_z_at_xy_below_3 compare_z_at_xy_below_3_object() const
  {
    return Compare_z_at_xy_below_3(this);
  }
  /*@}*/
  
protected:
  //! construct_conic
  /*! The function constructs X_monotone_curve_2 objects of CKvA that 
    repesent the conic: r*x^2 + s*y^2 + t*x*y + u*x + v*y + w = 0
    The output is written to given output iterator.
    \return The updated output iterator - value_type of X_monotone_curve_2
  */
  template <class OutputIterator>
  OutputIterator construct_conic (
                                  const Rational& r,
                                  const Rational& s,
                                  const Rational& t,
                                  const Rational& u,
                                  const Rational& v,
                                  const Rational& w,
                                  OutputIterator o) const
  {
    typedef CGAL::Polynomial< Coefficient >           Polynomial_1;
    typedef CGAL::Polynomial< Polynomial_1 >          Polynomial_2;
    typedef CGAL::Polynomial_traits_d< Polynomial_1 > Polynomial_traits_1;
    typedef CGAL::Polynomial_traits_d< Polynomial_2 > Polynomial_traits_2;
        
    typedef std::list<CGAL::Object> Object_list;
    typedef CGAL::Exponent_vector EV;
    typedef std::pair<EV, Coefficient> Monomial;
    typedef CGAL::Fraction_traits< Rational > FT;

    typename FT::Decompose decompose;

    // getting integer from the rational coefficients.
    Coefficient rn, rd;
    Coefficient sn, sd;
    Coefficient tn, td;
    Coefficient un, ud;
    Coefficient vn, vd;
    Coefficient wn, wd;
    decompose(r, rn, rd);
    decompose(s, sn, sd);
    decompose(t, tn, td);
    decompose(u, un, ud);
    decompose(v, vn, vd);
    decompose(w, wn, wd);
        
    Coefficient d = rd*sd*td*ud*vd*wd;
    // need to multiply each numenator with the multiplication of the
    // denominator (except the current denom).
    rn *= sd*td*ud*vd*wd;
    sn *= rd*td*ud*vd*wd;
    tn *= rd*sd*ud*vd*wd;
    un *= rd*sd*td*vd*wd;
    vn *= rd*sd*td*ud*wd;
    wn *= rd*sd*td*ud*vd;
        
    std::vector<Monomial> monomials;
        
    monomials.push_back(Monomial(EV(2, 0), rn));
    monomials.push_back(Monomial(EV(0, 2), sn));
    monomials.push_back(Monomial(EV(1, 1), tn));
    monomials.push_back(Monomial(EV(1, 0), un));
    monomials.push_back(Monomial(EV(0, 1), vn));
    monomials.push_back(Monomial(EV(0, 0), wn));
        
    typename Polynomial_traits_2::Construct_polynomial 
      construct_polynomial;
    Polynomial_2 result = construct_polynomial(monomials.begin(),
                                               monomials.end());
        
    Object_list x_objs;
        
    Curve_2 cur = 
      Base::instance().kernel().construct_curve_2_object()(result);
        
    this->make_x_monotone_2_object()
      (cur, std::back_inserter(x_objs));
        
        
    for (Object_list::iterator it = x_objs.begin();
         it != x_objs.end(); ++it)
    {
      const X_monotone_curve_2 *xcurve;
      if ((xcurve = object_cast<X_monotone_curve_2> (&(*it))) != NULL)
      {
        *o++ = *xcurve;

      }
      else
      {
        CGAL_error_msg("This function supports only X_monotone_curve_2");
      }

    }

    return o;
  }

  
  /*!
   * Returns the line bisector of a point and a segment that are collinear.
   * The bisector is output to the output iterator X_monotone_curve_2.
   * \param point The point.
   * \param seg The segment.
   * \param o The output iterator we output the result to.
   * \pre The point and the segment are collinear and the point is does not
   *      intersect the open segment.
   * \return The update output iterator
   */
  template <class OutputIterator>
  OutputIterator bisector_of_collinear_point_and_linear_obj(
                                                            const typename Site_2::Point_2& point, 
                                                            const Site_2& site,
                                                            OutputIterator o) const
  {
    CGAL_envelope_voronoi_precondition (                              \
                                        CGAL::are_ordered_along_line(point, site.point_on(0),           \
                                                                     site.point_on(1)) ||               \
                                        CGAL::are_ordered_along_line(point, site.point_on(1),           \
                                                                     site.point_on(0)));
      
    // The bisector is the perpendicular line to the supporting line at the
    // midpoint of point and one of the endpoints of site.
    typename Site_2::Point_2 p;
    if (site.is_ray())
    {
      p = site.ray().source();
    }
    else
    {
      CGAL_envelope_voronoi_assertion(site.is_segment());
      typename Site_2::Segment_2 seg = site.segment();
      if (CGAL::compare_xy(seg.min(), point) == SMALLER)
      {
        p = seg.max();
      }
      else
      {
        CGAL_envelope_voronoi_assertion (CGAL::compare_xy(seg.max(), point) \
                                         == LARGER);
        p = seg.min();              
      }
    }
      
    // The bisector is... see above.
    typename Site_2::Line_2 line = 
      site.supporting_line().perpendicular(CGAL::midpoint(p, point));

    return construct_conic(0, 0, 0, 
                           line.a(), line.b(), line.c(), o);
  }
    
};


} //namespace CGAL

#endif  // CGAL_LINEAR_OBJECTS_VORONOI_TRAITS_2_H
