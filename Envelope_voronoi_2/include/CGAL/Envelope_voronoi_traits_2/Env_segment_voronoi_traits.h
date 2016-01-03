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


/*! \file Env_segment_voronoi_traits.h
  Contains an implementation for the segment Voronoi diagram based on our
  extension for the conic traits. THIS IS A VERY EXPERIMENTAL VERSION AND
  PROBABLY DOES NOT WORK WELL.
*/

#ifndef CGAL_ENV_VORONOI_SEGMENT_TRAITS_H
#define CGAL_ENV_VORONOI_SEGMENT_TRAITS_H

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/representation_tags.h>
#include <CGAL/functions_on_enums.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/number_utils.h> 
#include <CGAL/Cartesian_converter.h>

#include <CGAL/Arr_trace_traits.h>

#include <CGAL/Envelope_3/Envelope_base.h>
#include <CGAL/Envelope_voronoi_2/envelope_voronoi_assertions.h>


namespace CGAL {

/*!
 * \class The class represents a segment voronoi site. To make life easier
 *        in computing the segment voronoi diagram we treat the end points 
 *        of the segment as sites on their own. This means that the class 
 *        should support both a segment site and a point site. The segment
 *        sites are open segment actually.
 */
template <class Env_segment_voronoi_traits_2_>
class _Segment_voronoi_site
{
 public:
  typedef Env_segment_voronoi_traits_2_                        Env_segment_voronoi_traits_2;
  typedef _Segment_voronoi_site<Env_segment_voronoi_traits_2> Self;
  typedef typename Env_segment_voronoi_traits_2::Rat_kernel    Rat_kernel;
  typedef typename Env_segment_voronoi_traits_2::Alg_kernel    Alg_kernel;
  typedef typename Rat_kernel::FT              Rational;
  typedef typename Alg_kernel::FT              Algebraic;

  typedef typename Rat_kernel::Point_2         Rat_point_2;
  typedef typename Rat_kernel::Segment_2       Rat_segment_2;

  typedef typename Alg_kernel::Point_2         Alg_point_2;
  typedef typename Alg_kernel::Segment_2       Alg_segment_2;

  typedef typename Env_segment_voronoi_traits_2::Point_2       Point_2;
  
 protected:
  Rat_segment_2 _rat_segment;
  Rat_point_2 _rat_point;
  
  Alg_segment_2 _alg_segment;
  Alg_point_2 _alg_point;

  bool _is_segment;
        
 public:
  _Segment_voronoi_site(const Rat_point_2& point)
    : _rat_point(point), _is_segment(false)
    {
      Cartesian_converter<Rat_kernel, Alg_kernel> conv;
      _alg_point = conv(_rat_point);
    }

  _Segment_voronoi_site(const Rat_segment_2& segment)
    : _rat_segment(segment), _is_segment(true)
    {
      Cartesian_converter<Rat_kernel, Alg_kernel> conv;
      _alg_segment = conv(_rat_segment);
    }

  const bool is_segment() const
    {
      return _is_segment;
    }
      
  const bool is_point() const
    {
      return !is_segment();
    }

  const Rat_segment_2 segment() const
    {
      CGAL_envelope_voronoi_precondition (is_segment());
      return _rat_segment;
    }
      
  const Alg_segment_2 alg_segment() const
    {
      CGAL_envelope_voronoi_precondition (is_segment());
      return _alg_segment;
    }

  const Rat_point_2 point() const
    {
      CGAL_envelope_voronoi_precondition (is_point());
      return _rat_point;
    }
        
  const Algebraic squared_distance(const Point_2& p) const
    {
      if (is_segment())
        return CGAL::squared_distance(p, _alg_segment);
      return CGAL::squared_distance(p, _alg_point);
    }

  bool operator== (const Self& s) const
    {
      if (is_segment() != s.is_segment())
        return false;

      if (is_segment())
        return (_rat_segment == s._rat_segment);

      return (_rat_point == s._rat_point);
    }
};


/*!
 * \class A traits class for computing segment voronoi diagram using lower envelope
 *        algorightm. The template parameters are needed because that the class is
 *        based on the Arr_conic_traits_2.
 *
 * The class is templated with two parameters: 
 * Rat_kernel A kernel that provides the input objects or coefficients.
 *            Rat_kernel::FT should be an integral or a rational type.
 * Alg_kernel A geometric kernel, where Alg_kernel::FT is the number type
 *            for the coordinates of arrangement vertices, which are algebraic
 *            numbers of degree up to 4 (preferably it is CORE::Expr).
 * Nt_traits A traits class for performing various operations on the integer,
 *           rational and algebraic types. 
 */
#ifdef DEBUG_CONIC_TRAITS
template <class Rat_kernel_, class Alg_kernel_, class Nt_traits_>
class Env_segment_voronoi_traits_2 : 
public CGAL::Arr_trace_traits_2< 
    Arr_conic_traits_2<Rat_kernel_, Alg_kernel_, Nt_traits_> >
#else
template <class Rat_kernel_, class Alg_kernel_, class Nt_traits_>
class Env_segment_voronoi_traits_2 : public Arr_conic_traits_2<Rat_kernel_, Alg_kernel_, Nt_traits_>
#endif
{
 public:
  typedef Rat_kernel_                          Rat_kernel;
  typedef Alg_kernel_                          Alg_kernel;
  typedef Nt_traits_                           Nt_traits;
  typedef typename Rat_kernel::FT              Rational;
  typedef typename Alg_kernel::FT              Algebraic;
  typedef Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits> 
    Base;
  typedef Env_segment_voronoi_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
    Self;
  typedef unsigned int                         Multiplicity;

  typedef typename Rat_kernel::Point_2         Rat_point_2;
  typedef typename Rat_kernel::Line_2          Rat_line_2;
  typedef typename Rat_kernel::Segment_2       Rat_segment_2;

  typedef typename Alg_kernel::Point_2         Alg_point_2;
  typedef typename Alg_kernel::Line_2          Alg_line_2;
  typedef typename Alg_kernel::Segment_2       Alg_segment_2;

  typedef typename Base::Point_2               Point_2;
  typedef typename Base::Curve_2               Curve_2;
  typedef typename Base::X_monotone_curve_2    X_monotone_curve_2;
  typedef std::pair<X_monotone_curve_2, 
    Multiplicity>                              Intersection_curve;
  typedef typename Base::Has_boundary_category Has_boundary_category;
    
  typedef _Segment_voronoi_site<Self>          Surface_3;
  typedef Surface_3                            Xy_monotone_surface_3;
    
  class Make_xy_monotone_3
    {
    public:
      
      template <class OutputIterator>
        OutputIterator operator()(const Surface_3& s,
                                  bool is_lower,
                                  OutputIterator o) 
        {
          *o++ = s;
          return o;
        }
    };
  
  Make_xy_monotone_3 make_xy_monotone_3_object()
    {
      return Make_xy_monotone_3();
    }
  
  class Compare_z_at_xy_3
    {
    public:
      
      Comparison_result operator()(const Point_2& p,
                                   const Xy_monotone_surface_3& h1,
                                   const Xy_monotone_surface_3& h2)
        {
          Comparison_result res = 
            CGAL::compare(h1.squared_distance(p), h2.squared_distance(p));
          if (res != EQUAL)
            return res;
          
          // If both sites are of the same type, then they are equal.
          // Otherwise, we preffer the point over the segment.
          if (h1.is_point() == h2.is_point())
            return res;

          if (h1.is_point())
            return SMALLER;
          else
            return LARGER;
        }

      Comparison_result operator()(const X_monotone_curve_2& cv,
                                   const Xy_monotone_surface_3& h1,
                                   const Xy_monotone_surface_3& h2)
        {
          Comparison_result res = operator ()(cv.source(), h1, h2);
          CGAL_envelope_voronoi_assertion (res == operator ()(cv.target(), h1, h2));
          return res;
        }
      
      
      Comparison_result operator()(const Xy_monotone_surface_3& h1,
                                   const Xy_monotone_surface_3& h2)

        {
          // If we got here then the site should be equal.
          CGAL_envelope_voronoi_assertion (h1 == h2);
          return EQUAL;
        }
    
    };
  
  Compare_z_at_xy_3 compare_z_at_xy_3_object()
    {
      return Compare_z_at_xy_3();
    }

  class Compare_z_at_xy_above_3
    {
    protected:     
    public:
      Comparison_result operator()(const X_monotone_curve_2& cv,
                                   const Xy_monotone_surface_3& h1,
                                   const Xy_monotone_surface_3& h2)
        {
          if (h1.is_segment() && h2.is_point())
            return CGAL::opposite(operator() (cv, h2, h1));

          if (h1.is_point() && h2.is_point())
          {
            return CGAL::opposite(_compare_above(h1.point(), h2.point()));
          }

          if (h1.is_segment() && h2.is_segment())
          {
            Rat_line_2 line1 = h1.segment().supporting_line();
            Rat_line_2 line2 = h2.segment().supporting_line();
            
            if ((line1 == line2) ||
                (line1 == line2.opposite()))
            {
              // the segments should not overlap.
              CGAL_envelope_voronoi_assertion_code 
                (
                  CGAL::Object obj = CGAL::intersection(h1.segment(), h2.segment());
                );
              CGAL_envelope_voronoi_assertion (object_cast<Rat_segment_2>(&obj) == NULL);
              
              // the two segments are on the same line. same as points.
              // we need to make sure that we don't take the same point
              // in case that the the two segment share an end point.
              if (CGAL::compare_xy(h1.segment().source(), 
                                   h2.segment().source()) == EQUAL)
                return CGAL::opposite(_compare_above(h1.segment().source(), 
                                                     h2.segment().target()));
              return CGAL::opposite(_compare_above(h1.segment().source(), 
                                                   h2.segment().source()));
            }
            
            // The segments should not intersect (except at their end points).
            CGAL_envelope_voronoi_assertion_code 
              (
                const Rat_point_2& source1 = h1.segment().source();
                const Rat_point_2& target1 = h1.segment().target();
                const Rat_point_2& source2 = h2.segment().source();
                const Rat_point_2& target2 = h2.segment().target();
              );
            CGAL_envelope_voronoi_assertion ((CGAL::compare_xy(source1, source2) == EQUAL) ||
                            (CGAL::compare_xy(source1, target2) == EQUAL) ||
                            (CGAL::compare_xy(target1, source2) == EQUAL) ||
                            (CGAL::compare_xy(target1, target2) == EQUAL) ||
                            (CGAL::do_intersect(h1.segment(), h2.segment()) == false));
            
            Alg_segment_2 segment1 = h1.alg_segment();

            Comparison_result res = CGAL::opposite(CGAL::orientation(cv.left(), 
                                                                     cv.right(), segment1.source()));
            if (res == EQUAL)
              res = CGAL::opposite(CGAL::orientation(cv.left(), cv.right(), segment1.target()));
              
            return res;
          }
          
          // h1 is a point, h2 is a segment
          CGAL_envelope_voronoi_assertion (h1.is_point() && h2.is_segment());

          // If the point is on the line of the segment then we treat it as
          // two point sites.
          if (CGAL::collinear(h1.point(), h2.segment().source(), h2.segment().target()))
          {
            // make sure that we don't take the same point in case that the point is 
            // on a segment end point.
            Comparison_result res = CGAL::opposite(_compare_above(h1.point(), h2.segment().source()));
            if (res == EQUAL)
              res = CGAL::opposite(_compare_above(h1.point(), h2.segment().target()));

            return res;
          }

          // point isn't collinear with the segment. the bisector is a parabolic
          // arc.
          if (cv.is_facing_up())
            return LARGER;
          else
            return SMALLER;
        }

    };

  Compare_z_at_xy_above_3 compare_z_at_xy_above_3_object()
    {
      return Compare_z_at_xy_above_3();
    }

  class Compare_z_at_xy_below_3
    {
    public:
      Comparison_result operator()(const X_monotone_curve_2& cv,
                                   const Xy_monotone_surface_3& h1,
                                   const Xy_monotone_surface_3& h2)
        {
          Compare_z_at_xy_above_3 cmp_above;
          return CGAL::opposite(cmp_above(cv, h1, h2));
        }

    };

  Compare_z_at_xy_below_3 compare_z_at_xy_below_3_object()
    {
      return Compare_z_at_xy_below_3();
    }


  class Construct_projected_boundary_2
    {
    public:

      template <class OutputIterator>
        OutputIterator operator()(const Xy_monotone_surface_3& s,
                                  OutputIterator o) const
        {
          typedef std::pair<X_monotone_curve_2, CGAL::Oriented_side> Boundary_curve;

          if (s.is_point())
            return o;

          // if this is a segment then the boundary of the surface is the vertical
          // lines from the endpoints of the segment.
          
          Curve_2 cv_above, cv_below;
          if (_compare_above(s.segment().source(), s.segment().target()) == LARGER)
          {
            cv_above = _bisector_of_collinear_point_and_segment(s.segment().source(), s.segment());
            cv_below = _bisector_of_collinear_point_and_segment(s.segment().target(), s.segment());
          }
          else
          {
            cv_above = _bisector_of_collinear_point_and_segment(s.segment().target(), s.segment());
            cv_below = _bisector_of_collinear_point_and_segment(s.segment().source(), s.segment());          
          }

          X_monotone_curve_2 xcurve;
          Self self;
          std::list<CGAL::Object> x_mono_list;
          self.make_x_monotone_2_object() (cv_above, std::back_inserter(x_mono_list));
          CGAL_envelope_voronoi_assertion (x_mono_list.size() == 1);
          if (assign(xcurve, x_mono_list.front()))
          {
            *o++ = make_object(Boundary_curve(xcurve, ON_NEGATIVE_SIDE));
          }

          x_mono_list.clear();
          self.make_x_monotone_2_object() (cv_below, std::back_inserter(x_mono_list));
          CGAL_envelope_voronoi_assertion (x_mono_list.size() == 1);
          if (assign(xcurve, x_mono_list.front()))
          {
            *o++ = make_object(Boundary_curve(xcurve, ON_POSITIVE_SIDE));
          }

          return o;
        }
    };

  Construct_projected_boundary_2 
    construct_projected_boundary_2_object()
    {
      return Construct_projected_boundary_2();
    }


  class Construct_projected_intersections_2
    {
      
    protected:

     
    public:

      template <class OutputIterator>
        OutputIterator operator()(const Xy_monotone_surface_3& s1,
                                  const Xy_monotone_surface_3& s2,
                                  OutputIterator o) const
        {
          Curve_2 bisector;
          if (s1.is_segment() && s2.is_point())
            return operator() (s2, s1, o);

          if (s1.is_point() && s2.is_point())
          {
            if (CGAL::compare_xy(s1.point(), s2.point()) == EQUAL)
              return o;

            Rat_line_2 line = CGAL::bisector(s1.point(), s2.point());
            bisector = Curve_2(0, 0, 0, line.a(), line.b(), line.c());
          }

          if (s1.is_segment() && s2.is_segment())
          {
            Rat_line_2 line1 = s1.segment().supporting_line();
            Rat_line_2 line2 = s2.segment().supporting_line();
            
            if ((line1 == line2) ||
                (line1 == line2.opposite()))
            {
              // the segments should not overlap.
              CGAL_envelope_voronoi_assertion_code 
                (
                  CGAL::Object obj = CGAL::intersection(s1.segment(), s2.segment());
                );
              CGAL_envelope_voronoi_assertion (object_cast<Rat_segment_2>(&obj) == NULL);

              if (CGAL::compare_xy(s1.segment().min(), 
                                   s2.segment().max()) == SMALLER)
              {
                bisector = _bisector_of_collinear_point_and_segment(s1.segment().max(),
                                                                    s2.segment());
              }
              else
              {
                bisector = _bisector_of_collinear_point_and_segment(s1.segment().min(),
                                                                    s2.segment());
              }
            }
            else
            {
              // The segments have different supporting lines.
              
              // The segments should not intersect (except at their end points).
              CGAL_envelope_voronoi_assertion_code 
                (
                  const Rat_point_2& source1 = s1.segment().source();
                  const Rat_point_2& target1 = s1.segment().target();
                  const Rat_point_2& source2 = s2.segment().source();
                  const Rat_point_2& target2 = s2.segment().target();
                  );
              CGAL_envelope_voronoi_assertion ((CGAL::compare_xy(source1, source2) == EQUAL) ||
                              (CGAL::compare_xy(source1, target2) == EQUAL) ||
                              (CGAL::compare_xy(target1, source2) == EQUAL) ||
                              (CGAL::compare_xy(target1, target2) == EQUAL) ||
                              (CGAL::do_intersect(s1.segment(), s2.segment()) == false));


              // The bisector of the two segments is like the bisector of their 2 supporting
              // lines - a degenerate hyperbola:
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
              Rat_line_2 line1 = s1.segment().supporting_line();
              Rat_line_2 line2 = s2.segment().supporting_line();
              
              Rational a1 = line1.a(); Rational b1 = line1.b(); Rational c1 = line1.c();
              Rational a2 = line2.a(); Rational b2 = line2.b(); Rational c2 = line2.c();
              
              Rational a1_2 = a1 * a1; Rational b1_2 = b1 * b1; 
              Rational a2_2 = a2 * a2; Rational b2_2 = b2 * b2; 

              Rational a1_b1 = a1_2 + b1_2; Rational a2_b2 = a2_2 + b2_2;

              Rational r = a1_2 * b2_2 - b1_2 * a2_2;
              Rational s = b1_2 * a2_2 - a1_2 * b2_2;
              Rational t = 2*a1*b1*a2_b2 - 2*a2*b2*a1_b1;
              Rational u = 2*c1*a1*a2_b2 - 2*c2*a2*a1_b1;
              Rational v = 2*b1*c1*a2_b2 - 2*b2*c2*a1_b1;
              Rational w = c1*c1*a2_b2 - c2*c2*a1_b1;

              // We also need two intersection points that will be used as the source
              // and target points. We use the algebraic kernel to get 2 points on the
              // bisector of these two line segments.
              
              // The segments need to be in the "correct" direction. Meaning, they both
              // have to face "outward" from their relative intersection point. 
              
              Rat_segment_2 seg1 = s1.segment();
              Rat_segment_2 seg2 = s2.segment();
              CGAL::Orientation o1 = CGAL::orientation(seg1.source(),seg1.target(), seg2.source());
              if (o1 == COLLINEAR)
                o1 = CGAL::orientation(seg1.source(),seg1.target(), seg2.target());
              
              CGAL::Orientation o2 = CGAL::orientation(seg2.source(), seg2.target(), seg1.source());
              if (o2 == COLLINEAR)
                o2 = CGAL::orientation(seg2.source(), seg2.target(), seg1.target());

              Alg_segment_2 alg_seg1 = s1.alg_segment();
              Alg_segment_2 alg_seg2 = s2.alg_segment();
              if (o1 == o2)
                alg_seg2 = alg_seg2.opposite();
              
              Alg_line_2 l = CGAL::bisector(alg_seg1.supporting_line(), alg_seg2.supporting_line());
              Alg_point_2 p1 = l.point(0);
              Alg_point_2 p2 = l.point(1);

              // creating the conic
              bisector = Curve_2 (r, s, t, u, v, w, COLLINEAR, 
                                  p1, p2, true, true);
            }
            
          }

          if (s1.is_point() && s2.is_segment())
          {
            if (CGAL::orientation(s1.point(), s2.segment().source(), 
                                  s2.segment().target()) == COLLINEAR)
            {
              bisector = _bisector_of_collinear_point_and_segment(s1.point(), s2.segment());
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
              Rational a = s2.segment().supporting_line().a();
              Rational b = s2.segment().supporting_line().b();
              Rational c = s2.segment().supporting_line().c();
              
              Rational p = s1.point().x(); Rational q = s1.point().y();

              Rational r = b*b;
              Rational s = a*a;
              Rational t = -2*a*b;

              Rational a_2_b_2 = a*a + b*b;
              Rational u = -2*(p*a_2_b_2 + c*a);
              Rational v = -2*(q*a_2_b_2 + c*b);
              Rational w = (p*p + q*q)*a_2_b_2 - c*c;

              bisector = Curve_2(r, s, t, u, v, w);
            }
          }

          // split into x-monotone
          Self self;
          std::list<Object> x_mono_list;
          self.make_x_monotone_2_object() (bisector, std::back_inserter(x_mono_list));
          
          typename list<Object>::iterator it;
          for (it = x_mono_list.begin(); it != x_mono_list.end(); ++it)
          {
            X_monotone_curve_2 x;
            if (assign(x, *it))
              *o++ = make_object(Intersection_curve(x, 1));
          }
          
          return o;
        }
    };

  Construct_projected_intersections_2 construct_projected_intersections_2_object()
    {
      return Construct_projected_intersections_2();
    }
  
 protected:
  /*!
   * Returns the line bisector of a point and a segment that are collinear.
   * \param point The point.
   * \param seg The segment.
   * \pre The point and the segment are collinear and the point is does not
   *      intersect the open segment.
   * \return The bisector between the point and the segment.
   */
  static Curve_2 _bisector_of_collinear_point_and_segment(const Rat_point_2& point, 
                                                          const Rat_segment_2& seg)
    {
      CGAL_envelope_voronoi_precondition (CGAL::orientation(point, seg.source(), 
                                           seg.target()) == COLLINEAR);
      CGAL_envelope_voronoi_precondition ((point == seg.source()) || (point == seg.target()) ||
                         (CGAL::do_intersect(point, seg) == false));
      
      // The bisector (except for the special case later) is the same bisector
      // of the two point: point and the closer end-point of seg.
      const Rat_point_2 *p;
      if (CGAL::compare_xy(seg.min(), point) == SMALLER)
      {
        p = &seg.max();
      }
      else
      {
        CGAL_envelope_voronoi_assertion (CGAL::compare_xy(seg.max(), point) == LARGER);
        p = &seg.min();              
      }
      
      // if the segment and point are touching then the bisector 
      // is the perpendicular line to the segment in that point.
      Rat_line_2 line;
      if (CGAL::compare_xy(*p, point) == EQUAL)
      {
        line = seg.supporting_line().perpendicular(*p);
      }
      else
      {
        line = CGAL::bisector(*p, point);
      }

      return Curve_2(0, 0, 0, line.a(), line.b(), line.c());
    }

  /*!
   * Checks which of the points is above the bisector of the two and
   * which is below the bisector.
   * \param p1 The first point.
   * \param p2 The second point.
   */
  static Comparison_result _compare_above(const Rat_point_2& p1,
                                     const Rat_point_2& p2)
    {
      Comparison_result res = CGAL::compare_y(p1, p2);
      if (res != EQUAL)
        return res;
      return CGAL::opposite(CGAL::compare_x(p1, p2));
    }


};

template <class T>
std::ostream& 
operator<< (std::ostream& os, 
            const CGAL::_Segment_voronoi_site<T>& site)
{
  if (site.is_point())  
    os << "Point: " << site.point();
  else
    os << "Segment: " << site.segment();

  return os;
}

} //namespace CGAL

#endif
