// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $URL: svn+ssh://ophirset@scm.gforge.inria.fr/svn/cgal/trunk/Arrangement_2/include/CGAL/Arr_conic_traits_2.h $
// $Id: Arr_conic_traits_2.h 5446 2007-11-27 17:16:54Z ophirset $
//
//
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>

#ifndef CGAL_ARR_CONIC_TRAITS_2_H
#define CGAL_ARR_CONIC_TRAITS_2_H

/*! \file
 * The conic traits-class for the arrangement package.
 */

#include <CGAL/tags.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_geometry_traits/Conic_arc_2.h>
#include <CGAL/Arr_geometry_traits/Conic_x_monotone_arc_2.h>
#include <CGAL/Arr_geometry_traits/Conic_point_2.h>

#include <fstream>

namespace CGAL {

#define DCT_TRACE(x)

/*!
 * \class A traits class for maintaining an arrangement of conic arcs (bounded
 * segments of algebraic curves of degree 2 at most).
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
template <typename Rat_kernel_, typename Alg_kernel_, typename Nt_traits_>
class Arr_conic_traits_2 {
public:
  using Rat_kernel = Rat_kernel_;
  using Alg_kernel = Alg_kernel_;
  using Nt_traits = Nt_traits_;
  using Rational = typename Rat_kernel::FT;
  using Rat_point_2 = typename Rat_kernel::Point_2;
  using Rat_segment_2 = typename Rat_kernel::Segment_2;
  using Rat_line_2 = typename Rat_kernel::Line_2;
  using Rat_circle_2 = typename Rat_kernel::Circle_2;
  using Algebraic = typename Alg_kernel::FT;
  using Integer = typename Nt_traits::Integer;
  using Self = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;
  using Has_left_category = Tag_true;
  using Has_merge_category = Tag_true;
  using Has_boundary_category = Tag_true;
  using Boundary_category = Arr_unbounded_boundary_tag;
  using Curve_2 = _Conic_arc_2<Rat_kernel, Alg_kernel, Nt_traits>;
  using X_monotone_curve_2 = _Conic_x_monotone_arc_2<Curve_2>;
  using Point_2 = _Conic_point_2<Alg_kernel>;

private:
  // Type definition for the intersection points mapping.
  using Conic_id = typename X_monotone_curve_2::Conic_id;
  using Intersection_point_2 = typename X_monotone_curve_2::Intersection_point_2;
  using Intersection_map = typename X_monotone_curve_2::Intersection_map;

  mutable Intersection_map  inter_map;   // Mapping conic pairs to their intersection
                                         // points.

public:
  /*! Default constructor.
   */
  Arr_conic_traits_2() {}

  /*! Get the next conic index. */
  static unsigned int get_index() {
    static unsigned int index = 0;
    return (++index);
  }

  /// \name Basic functor definitions.
  //@{

  class Compare_x_2 {
  public:
    /*! Compare the x-coordinates of two points.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2);
     *         SMALLER if x(p1) < x(p2);
     *         EQUAL if x(p1) = x(p2).
     */
    Comparison_result operator() (const Point_2 & p1, const Point_2 & p2) const
    {
      Alg_kernel   ker;
      return (ker.compare_x_2_object() (p1, p2));
    }

    /*! Compare the x-ordering of a vertical line passing through the point p
     * and an unbounded end of the curve xc.
     * \param p The first point.
     * \param xc The unbounded curve.
     * \param ce The end of the curve to check.
     * \pre The curve end has a bounded x-coordinate and an unbounded
     * y-coordinate.
     * \return LARGER if p > xc.
     *         SMALLER if p < xc.
     *         EQUAL if p = xc - a vertical line.
     */
    Comparison_result operator()(const Point_2& p, const X_monotone_curve_2& xc,
                                 Arr_curve_end ce) const
    {
      DCT_TRACE("Compare_x_2:" << endl);
      DCT_TRACE("Point: " << p << endl);
      DCT_TRACE("Curve: " << xc << endl);
      DCT_TRACE("End: " << ce << endl);

      CGAL_precondition (xc.infinite_in_x(ce) == ZERO);
      CGAL_precondition (xc.infinite_in_y(ce) != ZERO);
      CGAL_precondition (xc.is_vertical() || xc.is_hyperbola());

      Alg_kernel ker;
      if (xc.is_vertical()) {
        DCT_TRACE("Result: " << ker.compare_x_2_object()(p, xc.source()) << endl);
        return ker.compare_x_2_object() (p, xc.source());
      }

      CGAL_assertion(xc.is_hyperbola() && xc.has_vertical_asymptote());
      typename Alg_kernel::Line_2 asymptote = xc.get_algebraic_vertical_asymptote();
      Comparison_result res = ker.compare_x_at_y_2_object()(p, asymptote);
      if (res != EQUAL) {
        DCT_TRACE("Result: " << res << endl);
        return res;
      }

      DCT_TRACE("Result: " << ker.compare_x_2_object()(p, xc.source()) << endl);
      return ker.compare_x_2_object() (p, xc.source());
    }

    /*! Compare the x-ordering of the unbounded curve ends of xc1 and xc2.
     * \param xc1 The first unbounded curve.
     * \param ce1 The end of the first curve to check.
     * \param xc2 The second unbounded curve.
     * \param ce2 The end of the second curve to check.
     * \pre The xc1 curve end has a bounded x-coordinate and an unbounded y-coordinate.
     * \pre The xc2 curve end has a bounded x-coordinate and an unbounded y-coordinate.
     * \return LARGER if xc1 > xc2.
     *         SMALLER if xc1 < xc2.
     *         EQUAL if p = xc - they are the same curve.
     */
    Comparison_result operator()(const X_monotone_curve_2& xc1,
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2& xc2,
                                 Arr_curve_end ce2) const {
      DCT_TRACE("Compare_x_2:" << endl);
      DCT_TRACE("Curve1: " << xc1 << endl);
      DCT_TRACE("End1: " << ce1 << endl);
      DCT_TRACE("Curve2: " << xc2 << endl);
      DCT_TRACE("End2: " << ce2 << endl);
      DCT_TRACE("Result: Unknown" << endl);

      CGAL_precondition (xc1.infinite_in_x(ce1) == ZERO);
      CGAL_precondition (xc1.infinite_in_y(ce1) != ZERO);
      CGAL_precondition (xc2.infinite_in_x(ce2) == ZERO);
      CGAL_precondition (xc2.infinite_in_y(ce2) != ZERO);

      if (xc1.is_vertical())
        return operator()(xc1.source(), xc2, ce2);
      if (xc2.is_vertical())
        return opposite(operator()(xc2.source(), xc1, ce1));

      // both are not vertical
      CGAL_assertion(xc1.is_hyperbola());
      CGAL_assertion(xc2.is_hyperbola());

      // if they don't have the same asymptote that's easy.
      typename Rat_kernel::Line_2 asymptote1 = xc1.get_vertical_asymptote();
      typename Rat_kernel::Line_2 asymptote2 = xc2.get_vertical_asymptote();

      Rat_kernel rat_ker;
      typename Rat_kernel::Point_2 p = rat_ker.construct_point_2_object()(ORIGIN);
      Comparison_result res =
        rat_ker.compare_x_at_y_2_object() (p, asymptote1, asymptote2);
      if (res != EQUAL) return res;

      // if they are not in the same side of the asymptote then it is also easy.
      Alg_kernel ker;
      Comparison_result side_1 = ker.compare_x_at_y_2_object()
        (xc1.source(), xc1.get_algebraic_vertical_asymptote());
      Comparison_result side_2 = ker.compare_x_at_y_2_object()
        (xc2.source(), xc2.get_algebraic_vertical_asymptote());
      if (side_1 != side_2) return side_1;

      // They are in the same side - we compare the y ordering of some point
      // close to the hyperbola.
      CGAL_assertion_msg((side_1 != EQUAL) && (side_2 != EQUAL),
                         "point on asymptote");

      // Each of the hyperbolas are of the form: (x + m) (y - a*x - b) = c
      // (normalized by t) where m is the same for both hyperbolas.
      // First we look at the ratio of c(s). If they are equal we look at the
      // a*x + b expression (see below) and if they are still equal, we look at
      // the "a" factor.
      Nt_traits nt_traits;
      const Rational r1 = nt_traits.convert_to_rational (xc1.r());
      const Rational t1 = nt_traits.convert_to_rational (xc1.t());
      const Rational v1 = nt_traits.convert_to_rational (xc1.v());
      const Rational u1 = nt_traits.convert_to_rational (xc1.u());
      const Rational w1 = nt_traits.convert_to_rational (xc1.w());

      const Rational r2 = nt_traits.convert_to_rational (xc2.r());
      const Rational t2 = nt_traits.convert_to_rational (xc2.t());
      const Rational v2 = nt_traits.convert_to_rational (xc2.v());
      const Rational u2 = nt_traits.convert_to_rational (xc2.u());
      const Rational w2 = nt_traits.convert_to_rational (xc2.w());

      CGAL_assertion ((CGAL::sign(t1)!=ZERO) && (CGAL::sign(t2))!=ZERO);

      const Rational m = v1 / t1;
      CGAL_assertion_msg (m == v2 / t2, "hyperbolas don't have the same asymptote");

      const Rational a1 = -r1 / t1;
      const Rational a2 = -r2 / t2;

      const Rational b1 = -(u1/t1 + a1*m);
      const Rational b2 = -(u2/t2 + a2*m);

      const Rational c1 = -(w1/t1 + b1*m);
      const Rational c2 = -(w2/t2 + b2*m);

      Sign inf_dir = xc1.infinite_in_x(ce1);
      Sign s = CGAL::sign (c1 - c2);
      if (s != ZERO) return s * side_1 * inf_dir;

      // if they are equal then it is detemined by comparing the two -m*a+b expressions:
      s = CGAL::sign (-m*(a1-a2) + (b1-b2));
      if (s != ZERO) return s;

      // compare the a(s)
      s = CGAL::sign (a1 - a2);
      if (s != ZERO) return s * side_1 * inf_dir;

      return EQUAL;
    }
  };

  /*! Get a Compare_x_2 functor object. */
  Compare_x_2 compare_x_2_object() const { return Compare_x_2(); }

  Compare_x_2 compare_x_near_boundary_2_object() const { return Compare_x_2(); }

  class Compare_xy_2 {
  public:
    /*! Compares two points lexigoraphically: by x, then by y.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2), or if x(p1) = x(p2) and y(p1) > y(p2);
     *         SMALLER if x(p1) < x(p2), or if x(p1) = x(p2) and y(p1) < y(p2);
     *         EQUAL if the two points are equal.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const {
      Alg_kernel ker;
      return ker.compare_xy_2_object()(p1, p2);
    }
  };

  /*! Get a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object() const { return Compare_xy_2(); }

  /*! A function object that determines whether a curve end is bounded.
   */
  class Is_bounded_2 {
  public:
    /*! Is the end of an x-monotone curve bounded?
     * \param xcv The x-monotone curve.
     * \param ce The end of xcv identifier.
     * \return true is the curve end is bounded, and false otherwise
     */
    bool operator() (const X_monotone_curve_2 & xcv, Arr_curve_end ce) const {
      Sign sx = xcv.infinite_in_x(ce);
      Sign sy = xcv.infinite_in_y(ce);
      return (sx == ZERO && sy == ZERO) ? true : false;
    }
  };

  /*! Obtain a Is_bounded_2 function object. */
  Is_bounded_2 is_bounded_2_object() const { return Is_bounded_2(); }

  class Parameter_space_in_x_2 {
  public:
    /*! Determines if the x-coordinate of the curve xc is bounded.
     * \param xc The curve.
     * \param ce The curve end.
     * \return ARR_LEFT_BOUNDARY - if the curve grows to minus infinity.
     *         ARR_RIGHT_BOUNDARY - if the curve grows to plus infinity.
     *         ARR_INTERIOR - if it is finite.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2& xc,
                                   Arr_curve_end ce) {
      DCT_TRACE("Parameter_space_in_x_2:" << endl);
      DCT_TRACE("Curve: " << xc << endl);
      DCT_TRACE("End: " << ce << endl);

      DCT_TRACE("Result: " << (xc.infinite_in_x(ce)) << endl);
      Sign s = xc.infinite_in_x(ce);
      switch (s) {
       case POSITIVE: return ARR_RIGHT_BOUNDARY; break;
       case NEGATIVE: return ARR_LEFT_BOUNDARY; break;
       default: return ARR_INTERIOR;
      }
    }
  };

  /*! Get a Parameter_space_in_x_2 functor object. */
  Parameter_space_in_x_2 parameter_space_in_x_2_object () const
  { return Parameter_space_in_x_2(); }

  class Parameter_space_in_y_2 {
  public:
    /*! Determines if the y-coordinate of the curve xc is bounded.
     * \param xc The curve.
     * \param ce The curve end.
     * \return ARR_LEFT_BOUNDARY - if the curve grows to minus infinity.
     *         ARR_RIGHT_BOUNDARY - if the curve grows to plus infinity.
     *         ARR_INTERIOR - if it is finite.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2& xc,
                                   Arr_curve_end ce) {
      DCT_TRACE("Parameter_space_in_y_2:" << endl);
      DCT_TRACE("Curve: " << xc << endl);
      DCT_TRACE("End: " << ce << endl);
      DCT_TRACE("Result: " <<  xc.infinite_in_y(ce) << endl);
      Sign s = xc.infinite_in_y(ce);
      switch (s) {
       case POSITIVE: return ARR_RIGHT_BOUNDARY; break;
       case NEGATIVE: return ARR_LEFT_BOUNDARY; break;
       default: return ARR_INTERIOR;
      }
    }
  };

  /*! Get a Parameter_space_in_y_2 functor object. */
  Parameter_space_in_y_2 parameter_space_in_y_2_object () const
  { return Parameter_space_in_y_2(); }

  class Construct_min_vertex_2 {
  public:
    /*! Get the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \pre The curve is not unbounded to the minimum direction.
     * \return The left endpoint.
     */
    const Point_2& operator() (const X_monotone_curve_2 & cv) const {
      DCT_TRACE("Construct_min_vertex_2: " << cv.left() << endl);
      CGAL_precondition (!cv.is_left_unbounded());
      return (cv.left());
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object () const
  { return Construct_min_vertex_2(); }

  class Construct_max_vertex_2 {
  public:
    /*! Get the right endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \pre cv is not unbounded to the maximum direction.
     * \return The right endpoint.
     */
    const Point_2& operator() (const X_monotone_curve_2 & cv) const {
      DCT_TRACE("Construct_max_vertex_2: " << cv.right() << endl);
      CGAL_precondition (!cv.is_right_unbounded());
      return (cv.right());
    }
  };

  /*! Get a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object () const
  { return Construct_max_vertex_2(); }

  class Is_vertical_2 {
  public:
    /*! Check whether the given x-monotone curve is a vertical segment.
     * \param cv The curve.
     * \return (true) if the curve is a vertical segment; (false) otherwise.
     */
    bool operator() (const X_monotone_curve_2& cv) const {
      DCT_TRACE("Is_vertical_2:" << endl);
      DCT_TRACE("Result: " << cv.is_vertical() << endl);
      return (cv.is_vertical());
    }
  };

  /*! Get an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object () const { return Is_vertical_2(); }

  class Compare_y_at_x_2 {
  private:
    mutable Intersection_map& m_inter_map;   // The map of intersection points.
  public:
    /*! Constructor. */
    // TODO: check with Ron const + mutable.
    Compare_y_at_x_2 (Intersection_map& map) : m_inter_map(map) {}

    /*! Return the location of the given point with respect to the input curve.
     * \param cv The curve.
     * \param p The point.
     * \pre p is in the x-range of cv.
     * \return SMALLER if y(p) < cv(x(p)), i.e. the point is below the curve;
     *         LARGER if y(p) > cv(x(p)), i.e. the point is above the curve;
     *         EQUAL if p lies on the curve.
     */
    Comparison_result operator()(const Point_2 & p,
                                 const X_monotone_curve_2 & cv) const {
      DCT_TRACE("Compare_y_at_x_2:" << endl);
      DCT_TRACE("Curve: " << cv << endl);
      DCT_TRACE("Point: " << p << endl);

      Alg_kernel ker;
      if (cv.is_vertical()) {
        // A special treatment for vertical segments:
        // In case p has the same x c-ordinate of the vertical segment, compare
        // it to the segment endpoints to determine its position.
        Comparison_result res1 = ker.compare_y_2_object() (p, cv.left());
        Comparison_result res2 = ker.compare_y_2_object() (p, cv.right());

        res1 = (cv.is_left_unbounded()) ? SMALLER : res1;
        res2 = (cv.is_right_unbounded()) ? LARGER : res2;

        if (res1 == res2) return res1;
        else return EQUAL;
      }

      // Check whether the point is exactly on the curve.
      if (cv.contains_point(p)) return EQUAL;

      // Get a point q on the x-monotone arc with the same x coordinate as p.
      Comparison_result x_res = ker.compare_x_2_object() (p, cv.left());
      Point_2 q;

      if (! cv.is_left_unbounded() && x_res == EQUAL) q = cv.left();
      else {
        CGAL_precondition (cv.is_left_unbounded() || x_res != SMALLER);

        x_res = ker.compare_x_2_object() (p, cv.right());
        if (! cv.is_right_unbounded() && x_res == EQUAL)
          q = cv.right();
        else {
          CGAL_precondition (cv.is_right_unbounded() || x_res != LARGER);
          q = cv.get_point_at_x (p);
        }
      }

      // Compare p with the a point of the curve with the same x coordinate.
      DCT_TRACE("Result: " << (ker.compare_y_2_object() (p, q)) << endl);
      return (ker.compare_y_2_object() (p, q));
    }

    /*! Returns SMALLER, EQUAL or LARGER according to the y-ordering of the
     * two curves xc1 and xc2 at x = +-oo.
     * \param xc1 First curve.
     * \param xc2 Second curve.
     * \param ce What end to check.
     * \pre Both xc1 and xc2 are unbounded in ce direction.
     * \return SMALLER if xc1 < xc2.
     *         LARGER if xc1 > xc2.
     *         EQUAL if xc1 = xc2 - are the same curve.
     */
    Comparison_result operator()(const X_monotone_curve_2 & xc1,
                                 const X_monotone_curve_2 & xc2,
                                 Arr_curve_end ce) const {
      CGAL_precondition(Parameter_space_in_x_2() (xc1, ce) != ARR_INTERIOR);
      CGAL_precondition(Parameter_space_in_x_2() (xc2, ce) != ARR_INTERIOR);

      Comparison_result ret = xc1.compare_y_at_infinity(xc2, ce, m_inter_map);
      CGAL_expensive_postcondition
        (ret == xc1.compare_y_at_infinity(xc2, ce, m_inter_map, false));
      return ret;
    }
  };

  /*! Get a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object() const {
    // TODO: check with Ron about mutable. What is the meaning
    // of const in traits classes.
    // TODO: if mutable - change to const function
    return Compare_y_at_x_2(inter_map);
  }

  Compare_y_at_x_2 compare_y_near_boundary_2_object () const
  { return Compare_y_at_x_2(inter_map); }

  class Compare_y_at_x_left_2 {
  public:
    /*! Compares the y value of two x-monotone curves immediately to the left
     * of their intersection point.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its left.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the left of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator()(const X_monotone_curve_2& cv1,
                                 const X_monotone_curve_2& cv2,
                                 const Point_2& p) const {
      DCT_TRACE("Compare_y_at_x_left_2:" << endl);
      DCT_TRACE("Curve1: " << cv1 << endl);
      DCT_TRACE("Curve2: " << cv2 << endl);
      DCT_TRACE("Point: " << p << endl);

      // Make sure that p lies on both curves, and that both are defined to its
      // left (so their left endpoint is lexicographically smaller than p).
      CGAL_precondition(cv1.contains_point(p) && cv2.contains_point(p));

      CGAL_precondition_code(Alg_kernel ker;);
      CGAL_precondition(cv1.is_left_unbounded() ||
                        ker.compare_xy_2_object()(p, cv1.left()) == LARGER);
      CGAL_precondition(cv2.is_left_unbounded() ||
                        ker.compare_xy_2_object()(p, cv2.left()) == LARGER);

      // If one of the curves is vertical, it is below the other one.
      if (cv1.is_vertical()) {
        if (cv2.is_vertical()) {
          // Both are vertical:
          DCT_TRACE("Res: " << EQUAL << endl);
          return EQUAL;
        }
        else {
          DCT_TRACE("Res: " << SMALLER << endl);
          return SMALLER;
        }
      }
      else if (cv2.is_vertical()) {
        DCT_TRACE("Res: " << LARGER << endl);
        return LARGER;
      }

      // Compare the two curves immediately to the left of p:
      DCT_TRACE("Res: " << cv1.compare_to_left(cv2, p) << endl);
      return cv1.compare_to_left(cv2, p);
    }
  };

  /*! Get a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(); }

  class Compare_y_at_x_right_2 {
  public:
    /*! Compares the y value of two x-monotone curves immediately to the right
     * of their intersection point.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its right.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the right of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator()(const X_monotone_curve_2& cv1,
                                 const X_monotone_curve_2& cv2,
                                 const Point_2& p) const {
      DCT_TRACE("Compare_y_at_x_right_2:" << endl);
      DCT_TRACE("Curve1: " << cv1 << endl);
      DCT_TRACE("Curve2: " << cv2 << endl);
      DCT_TRACE("Point: " << p << endl);

      // Make sure that p lies on both curves, and that both are defined to its
      // left (so their left endpoint is lexicographically smaller than p).
      CGAL_precondition(cv1.contains_point(p) && cv2.contains_point(p));

      CGAL_precondition_code(Alg_kernel ker;);

      CGAL_precondition(cv1.is_right_unbounded() ||
                        ker.compare_xy_2_object()(p, cv1.right()) == SMALLER);
      CGAL_precondition(cv2.is_right_unbounded() ||
                        ker.compare_xy_2_object()(p, cv2.right()) == SMALLER);

      // If one of the curves is vertical, it is above the other one.
      if (cv1.is_vertical()) {
        if (cv2.is_vertical()) {
          // Both are vertical:
          DCT_TRACE("Res: " << EQUAL << endl);
          return EQUAL;
        }
        else {
          DCT_TRACE("Res: " << LARGER << endl);
          return LARGER;
        }
      }
      else if (cv2.is_vertical()) {
        DCT_TRACE("Res: " << SMALLER << endl);
        return SMALLER;
      }

      // Compare the two curves immediately to the right of p:
      DCT_TRACE("Res: " << cv1.compare_to_right (cv2, p) << endl);
      return cv1.compare_to_right(cv2, p);
    }
  };

  /*! Get a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(); }

  class Equal_2 {
  public:
    /*! Check if the two x-monotone curves are the same (have the same graph).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are the same; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const {
      DCT_TRACE("Equal_2!!!!1" << endl);
      if (&cv1 == &cv2) return true;
      return cv1.equals(cv2);
    }

    /*! Check if the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    bool operator()(const Point_2& p1, const Point_2& p2) const {
      DCT_TRACE("Equal_2!!!!1" << endl);
      if (&p1 == &p2) return (true);
      Alg_kernel ker;
      return (ker.compare_xy_2_object()(p1, p2) == EQUAL);
    }
  };

  /*! Get an Equal_2 functor object. */
  Equal_2 equal_2_object() const { return Equal_2(); }
  //@}

  /// \name Functor definitions for supporting intersections.
  //@{

  class Make_x_monotone_2 {
    using Self = Arr_conic_traits_2 <Rat_kernel_, Alg_kernel_, Nt_traits_>;
    using Nt_traits = typename Self::Nt_traits;
    using Rational = typename Nt_traits::Rational;
    using Algebraic = typename Nt_traits::Algebraic;
  public:
    /*! Cut the given conic curve (or conic arc) into x-monotone subcurves
     * and insert them to the given output iterator.
     * \param cv The curve.
     * \param oi The output iterator, whose value-type is Object. The returned
     *           objects are all wrappers X_monotone_curve_2 objects.
     * \return The past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator() (const Curve_2& cv, OutputIterator oi) {
      // Increment the serial number of the curve cv, which will serve as its
      // unique identifier.
      unsigned int index = Self::get_index();
      Conic_id conic_id (index);

      // in case that this is a full conic
      if (cv.is_full_conic())
        return _split_full_conic(cv, oi, conic_id);

      // Find the points of vertical tangency to cv and act accordingly.
      using Curve_point_2 = typename Curve_2::Point_2;
      Curve_point_2  vtan_ps[2];
      int n_vtan_ps;

      n_vtan_ps = cv.vertical_tangency_points (vtan_ps);

      if (n_vtan_ps == 0) {
        // In case the given curve is already x-monotone:
        *oi++ = make_object (X_monotone_curve_2 (cv, conic_id));
        return oi;
      }

      if (n_vtan_ps == 1) {
        Alg_kernel ker;

        // Split the arc into two x-monotone sub-curves: one going from the
        // arc source to ps[0], and the other from ps[0] to the target.
        // If one of the arc's ends is unbounded there is a chance that one
        // of the end points is exactly where the v_tan point is. In this case
        // We find another point on the curve to be used as target / source.

        const Curve_point_2* p_point = nullptr;
        Curve_point_2  ps[2];
        if (ker.compare_xy_2_object() (cv.source(), vtan_ps[0]) != EQUAL)
          p_point = &cv.source();
        else {
          CGAL_assertion_msg(cv.is_source_unbounded(),
                             "The points can be equal only"
                             " if we are dealing with unbounded case.");
          // we find another the point like the target, but in the other
          // branch of the curve.
          int n_ps = cv.get_points_at_x(cv.target(), ps);
          (void) n_ps; // disables warning of unused variable.
          CGAL_assertion (n_ps == 2);
          if (ker.compare_xy_2_object() (ps[0], cv.target()) == EQUAL)
            p_point = &ps[1];
          else {
            CGAL_assertion
              (ker.compare_xy_2_object()(ps[1], cv.target()) == EQUAL);
            p_point = &ps[0];
          }
        }
        *oi++ = make_object(X_monotone_curve_2(cv, *p_point, vtan_ps[0],
                                               conic_id,
                                               cv.is_source_unbounded(),
                                               false));

        // same thing for target
        if (ker.compare_xy_2_object()(cv.target(), vtan_ps[0]) != EQUAL)
          p_point = &cv.target();
        else {
          CGAL_assertion_msg(cv.is_target_unbounded(),
                              "The points can be equal only"
                              " if we are dealing with unbounded case.");
          // we find another the point like the source, but on the other
          // branch of the curve.
          int n_ps = cv.get_points_at_x(cv.source(), ps);
          (void) n_ps; // disables warning of unused variable.
          CGAL_assertion(n_ps == 2);
          if (ker.compare_xy_2_object() (ps[0], cv.source()) == EQUAL)
            p_point = &ps[1];
          else {
            CGAL_assertion
              (ker.compare_xy_2_object() (ps[1], cv.source()) == EQUAL);
            p_point = &ps[0];
          }
        }
        *oi++ = make_object(X_monotone_curve_2(cv, vtan_ps[0], *p_point,
                                               conic_id, false,
                                               cv.is_target_unbounded()));
      }
      else {
        CGAL_assertion (n_vtan_ps == 2);

        // Identify the first point we encounter when going from cv's source
        // to its target, and the second point we encounter. Note that the
        // two endpoints must both be below the line connecting the two
        // tangnecy points (or both lies above it).
        int ind_first = 0;
        int ind_second = 1;
        Alg_kernel_ker;
        typename Alg_kernel_::Line_2 line =
          ker.construct_line_2_object()(vtan_ps[0], vtan_ps[1]);
        const Comparison_result start_pos =
          ker.compare_y_at_x_2_object()(cv.source(), line);
        const Comparison_result order_vpts =
          ker.compare_x_2_object()(vtan_ps[0], vtan_ps[1]);

        CGAL_assertion(start_pos != EQUAL &&
                       ker.compare_y_at_x_2_object()(cv.target(),
                                                     line) == start_pos);
        CGAL_assertion (order_vpts != EQUAL);

        if ((cv.orientation() == COUNTERCLOCKWISE &&
             start_pos == order_vpts) ||
            (cv.orientation() == CLOCKWISE &&
             start_pos != order_vpts))
        {
          ind_first = 1;
          ind_second = 0;
        }

        // Split the arc into three x-monotone sub-curves.
        *oi++ = make_object (X_monotone_curve_2(cv, cv.source(),
                                                vtan_ps[ind_first],
                                                conic_id));

        *oi++ = make_object (X_monotone_curve_2(cv, vtan_ps[ind_first],
                                                vtan_ps[ind_second],
                                                conic_id));

        *oi++ = make_object (X_monotone_curve_2(cv, vtan_ps[ind_second],
                                                cv.target(),
                                                conic_id));
      }

      return oi;
    }

  private:
    /*! Cut the given full conic curve into x-monotone subcurves and insert
     * them to the given output iterator. It handle whole kinds of x-monotone
     * curves.
     * \param cv The curve.
     * \param oi The output iterator, whose value-type is Object. The returned
     *           objects are all wrappers X_monotone_curve_2 objects.
     * \param conic_id The conic_id to use for the subcurves (caching perpuses).
     * \return The past-the-end iterator.
     */
    template<typename OutputIterator>
    OutputIterator _split_full_conic(const Curve_2& cv, OutputIterator oi,
                                     const Conic_id &conic_id) {
      using Curve_point_2 = typename Curve_2::Point_2;
      // If this is a full conic, we need to split the curve (and probably
      // calculate end points.
      // For the case that this is a regular line, we just need to calculate
      // end points.
      if (CGAL::sign(cv.r()) == ZERO && CGAL::sign(cv.s()) == ZERO &&
          CGAL::sign(cv.t()) == ZERO) {
        Rational u = cv.u();
        Rational v = cv.v();
        Rational w = cv.w();

        Curve_point_2 start_point, end_point;
        if (CGAL::sign (cv.v()) == ZERO) {
          // in case that the line is vertical:
          Rational x = -w/u;
          start_point = Curve_point_2(x, 0);
          end_point = Curve_point_2(x, 10);
        }
        else {
          start_point = Curve_point_2(0, - w/v);
          end_point = Curve_point_2(1, - (w+u)/v);
        }
        *oi++ = make_object(X_monotone_curve_2(cv, start_point, end_point,
                                               conic_id, true, true));
        return oi;
      }

      // If the curve is not a straight line then the type of the curve is
      // determined by (4rs - t^2):
      // 1. positive - ellipse
      // 2. zero - parabola / parallel lines / coincident lines
      // 3. negative - hyperbola / intersecting lines
      CGAL::Sign type = CGAL::sign(4*cv.r()*cv.s() - cv.t()*cv.t());
      switch (type) {
      case POSITIVE: // ellipse
        _split_full_ellipse(cv, oi, conic_id);
        break;

      case ZERO: // parabola / parallel lines / coincident lines
        _split_full_parabola(cv, oi, conic_id);
        break;

      case NEGATIVE: // hyperbola
        _split_full_hyperbola(cv, oi, conic_id);
        break;
      }

      return oi;
    }

    /*! Cut the given full ellipse into x-monotone subcurves and insert
     * them to the given output iterator.
     * \param cv The ellipse.
     * \param oi The output iterator, whose value-type is Object. The returned
     *           objects are all wrappers X_monotone_curve_2 objects.
     * \param conic_id The conic_id to use for the subcurves (caching perpuses).
     * \return The past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator _split_full_ellipse(const Curve_2& cv, OutputIterator oi,
                                       const Conic_id &conic_id) {
      // Find the points of vertical tangency to cv and act accordingly.
      using Curve_point_2 = typename Curve_2::Point_2;
      Curve_point_2  vtan_ps[2];
      int                        n_vtan_ps;
      n_vtan_ps = cv.vertical_tangency_points (vtan_ps);

      // This is a full ellipse. We make sure that it has exactly 2 vertical
      // tangency points and split it into two arcs in those points.
      CGAL_assertion(n_vtan_ps == 2);

      // In case the curve is a full conic, split it into two x-monotone
      // arcs, one going from ps[0] to ps[1], and the other from ps[1] to
      // ps[0].
      *oi++ = make_object(X_monotone_curve_2(cv, vtan_ps[0], vtan_ps[1],
                                             conic_id));
      *oi++ = make_object(X_monotone_curve_2(cv, vtan_ps[1], vtan_ps[0],
                                             conic_id));

      return oi;
    }

    /*! Cut the given full parabola into x-monotone subcurves and insert
     * them to the given output iterator.
     * \param cv The parabola.
     * \param oi The output iterator, whose value-type is Object. The returned
     *           objects are all wrappers X_monotone_curve_2 objects.
     * \param conic_id The conic_id to use for the subcurves (caching perpuses).
     * \return The past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator _split_full_parabola(const Curve_2& cv, OutputIterator oi,
                                        const Conic_id &conic_id) {
      // Find the points of vertical tangency to cv and act accordingly.
      using Curve_point_2 = typename Curve_2::Point_2;
      Curve_point_2 vtan_ps[2];
      int n_vtan_ps;
      n_vtan_ps = cv.vertical_tangency_points (vtan_ps);

      Nt_traits nt_traits;

      // In this case, there should be no more than 1 vertical tangency point.
      // If there is one we are dealing with a "regular" parabola that has two
      // branches. If there are no vertical tangency points then there are two
      // cases:
      // 1. parabola which has no vertical tangency points.
      // 2. parallel / coincident lines
      CGAL_assertion_msg(n_vtan_ps < 2,
                         "expecting less then 2 vertical tangency points in case of a parabola");

      if (n_vtan_ps == 1) {
        // This is a parabola with vertical tangency point. We try to find 2
        // end points on the parabola left to the vertical tangency point and
        // right to the vertical tangency points. One of them should get 2
        // intersection points with the parabola (the end points of the arc)
        // and the other should get none.
        CGAL_assertion (cv.orientation() != COLLINEAR);

        Curve_point_2 end_points[2];
        const Algebraic &tangency_x = vtan_ps[0].x();

        Rational left =
          nt_traits.rational_in_interval(tangency_x - 10, tangency_x);
        Curve_point_2 left_p (nt_traits.convert(left), 0);
        int n_end_points = cv.get_points_at_x(left_p, end_points);
        CGAL_assertion (n_end_points == 0 || n_end_points == 2);
        if (n_end_points == 0) {
          Rational right =
            nt_traits.rational_in_interval(tangency_x, tangency_x + 10);
          Curve_point_2 right_p (nt_traits.convert(right), 0);
          int n_end_points = cv.get_points_at_x(right_p, end_points);
          (void) n_end_points; // disables warning of unused variable.
          CGAL_assertion (n_end_points == 2);
        }

        // create x_monotone curves, but check that the orientation fits the orientation
        // of the curve.
        Curve_point_2 *source_p = &end_points[0];
        Curve_point_2 *target_p = &end_points[1];
        if (cv.orientation() !=
            CGAL::orientation(end_points[0], vtan_ps[0], end_points[1]))
          std::swap(source_p, target_p);

        *oi++ = make_object(X_monotone_curve_2
                            (cv, *source_p, vtan_ps[0], conic_id,
                              true, false));
        *oi++ = make_object(X_monotone_curve_2
                            (cv, vtan_ps[0], *target_p, conic_id,
                             false, true));
      }
      else if (n_vtan_ps == 0 && cv.orientation() != COLLINEAR) {
        // This a parabola with no vertical tengancy points, so we can get
        // the ends points easily.
        Curve_point_2 end_points[2];
        Curve_point_2 center_point[2];
        int n_points = cv.get_points_at_x(Curve_point_2(0, 0), &end_points[0]);
        CGAL_assertion (n_points == 1);
        n_points = cv.get_points_at_x(Curve_point_2(10, 0), &end_points[1]);
        CGAL_assertion (n_points == 1);

        n_points = cv.get_points_at_x(Curve_point_2(5, 0), center_point);
        CGAL_assertion (n_points == 1);

        if (cv.orientation() ==
            CGAL::orientation(end_points[0], center_point[0], end_points[1]))
        {
          *oi++ = make_object(X_monotone_curve_2(cv,
                                                 end_points[0], end_points[1],
                                                 conic_id, true, true));
        }
        else {
          *oi++ = make_object(X_monotone_curve_2(cv,
                                                 end_points[1], end_points[0],
                                                 conic_id, true, true));
        }
      }
      else {
        // if this is a degenerate parabola, there are 2 cases: parallel lines
        // and coincident lines. We just intersect the parabola with vertical
        // lines to get the end points. If the vertical line doesn't intersect
        // the parabola, we are dealing with a vertical line and then we
        // intersect the parabola with to horizontal lines.
        CGAL_assertion(cv.orientation() == COLLINEAR);
        Curve_point_2 left_points[2];
        Curve_point_2 right_points[2];
        int n_right_points = 0;
        int n_left_points = cv.get_points_at_x(Curve_point_2(0, 0), left_points);
        if (n_left_points != 0) {
          // no vertical lines.
          n_right_points = cv.get_points_at_x(Curve_point_2(10, 0), right_points);

        }
        else {
          n_left_points = cv.get_points_at_y(Curve_point_2(0, 0), left_points);
          n_right_points = cv.get_points_at_y(Curve_point_2(0, 10), right_points);
        }

        CGAL_assertion (n_left_points == n_right_points);
        CGAL_assertion (n_left_points > 0);
        *oi++ = make_object(X_monotone_curve_2(left_points[0], right_points[0],
                                               true, true));
        if (n_left_points > 1) {
          *oi++ = make_object(X_monotone_curve_2(left_points[1], right_points[1],
                                                 true, true));
        }
      }

      return oi;
    }

    /*! Cut the given full hyperbola into x-monotone subcurves and insert
     * them to the given output iterator.
     * \param cv The hyperbola.
     * \param oi The output iterator, whose value-type is Object. The returned
     *           objects are all wrappers X_monotone_curve_2 objects.
     * \param conic_id The conic_id to use for the subcurves (caching perpuses).
     * \return The past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator _split_full_hyperbola(const Curve_2& cv, OutputIterator oi,
                                         const Conic_id &conic_id) {
      // Find the points of vertical tangency to cv and act accordingly.
      using Curve_point_2 = typename Curve_2::Point_2;
      Curve_point_2  vtan_ps[2];
      int n_vtan_ps;
      n_vtan_ps = cv.vertical_tangency_points (vtan_ps);
      Nt_traits nt_traits;

      // There are two main cases:
      // 1. A regular hyperbola.
      // 2. Degenerate hyperbola - Which is actually 2 intersecting lines.

      if (cv.orientation() != COLLINEAR) {
        // regular hyperbola
        // 3 cases:
        // 1. The hyperbola has 2 vertical tangency points - We split the
        //    hyperbola into 4 rays.
        // 2. The hyperbola has no vertical tangency points - We split the
        //    hyperbola into 2 x-monotone branches.
        // 3. The hyperbola has vertical asymptot - We need to split the
        //    hyperbola into 2 branches.
        if (n_vtan_ps == 2) {
          CGAL_assertion (CGAL::compare_x(vtan_ps[0], vtan_ps[1]) == SMALLER);

          // we have left branch and right branch.
          Rational left = nt_traits.rational_in_interval(vtan_ps[0].x() - 10,
                                                         vtan_ps[0].x());
          Rational right = nt_traits.rational_in_interval(vtan_ps[1].x(),
                                                          vtan_ps[1].x() + 10);

          Curve_point_2 left_p (nt_traits.convert(left), 0);
          Curve_point_2 right_p (nt_traits.convert(right), 0);

          Curve_point_2 left_points[2];
          Curve_point_2 right_points[2];

          int n_left_points = cv.get_points_at_x(left_p, left_points);
          (void) n_left_points; // disables warning of unused variable.
          CGAL_assertion (n_left_points == 2);
          int n_right_points = cv.get_points_at_x(right_p, right_points);
          (void) n_right_points; // disables warning of unused variable.
          CGAL_assertion (n_right_points == 2);

          Curve_point_2* left_source = &left_points[0];
          Curve_point_2* left_target = &left_points[1];
          Curve_point_2* right_source = &right_points[1];
          Curve_point_2* right_target = &right_points[0];

          if (cv.orientation() !=
              CGAL::orientation(*left_source, vtan_ps[0], *left_target))
            std::swap(left_source, left_target);

          if (cv.orientation() !=
              CGAL::orientation(*right_source, vtan_ps[1], *right_target))
            std::swap(right_source, right_target);

          *oi++ = make_object(X_monotone_curve_2(cv, *left_source, vtan_ps[0],
                                                 conic_id, true, false));
          *oi++ = make_object(X_monotone_curve_2(cv, vtan_ps[0], *left_target,
                                                 conic_id, false, true));

          *oi++ = make_object(X_monotone_curve_2(cv, *right_source, vtan_ps[1],
                                                 conic_id, true, false));
          *oi++ = make_object(X_monotone_curve_2(cv, vtan_ps[1], *right_target,
                                                 conic_id, false, true));
        }
        else {
          // if we don't have vertical asymptot then we need to split
          // the arc into 2 branches.
          if (cv.has_vertical_asymptote() == false) {
            Curve_point_2 left_p (-10, 0);
            Curve_point_2 center_p (0, 0);
            Curve_point_2 right_p (10, 0);

            Curve_point_2 left_points[2];
            Curve_point_2 center_points[2];
            Curve_point_2 right_points[2];

            int n_left_points = cv.get_points_at_x(left_p, left_points);
            (void) n_left_points; // disables warning of unused variable.
            CGAL_assertion (n_left_points == 2);
            int n_center_points = cv.get_points_at_x(center_p, center_points);
            (void) n_center_points; // disables warning of unused variable.
            CGAL_assertion (n_center_points == 2);
            int n_right_points = cv.get_points_at_x(right_p, right_points);
            (void) n_right_points; // disables warning of unused variable.
            CGAL_assertion (n_right_points == 2);

            Curve_point_2 *bottom_source = &left_points[0];
            Curve_point_2 *bottom_center = &center_points[0];
            Curve_point_2 *bottom_target = &right_points[0];
            Curve_point_2 *top_source = &left_points[1];
            Curve_point_2 *top_center = &center_points[1];
            Curve_point_2 *top_target = &right_points[1];

            // orientation should be correct
            if (cv.orientation() !=
                CGAL::orientation(*bottom_source, *bottom_center,
                                  *bottom_target))
              std::swap(bottom_source, bottom_target);

            if (cv.orientation() !=
                CGAL::orientation(*top_source, *top_center, *top_target))
              std::swap(top_source, top_target);

            *oi++ = make_object(X_monotone_curve_2(cv, *bottom_source,
                                                   *bottom_target,
                                                   conic_id, true, true));
            *oi++ = make_object(X_monotone_curve_2(cv, *top_source,
                                                   *top_target,
                                                   conic_id, true, true));
          }
          else {
            // vertical asymptot.
            typename Rat_kernel::Line_2 vert_asympt =
              cv.get_vertical_asymptote();
            CGAL_assertion (vert_asympt.is_vertical() == true);

            Rational x_asympt = -vert_asympt.c() / vert_asympt.a();
            Curve_point_2 points[2];

            // the right branch
            int n_points = cv.get_points_at_x(Curve_point_2
                                              (nt_traits.convert(x_asympt + 10),
                                               0), points);
            CGAL_assertion (n_points == 1);
            Curve_point_2 right_1 = points[0];

            n_points = cv.get_points_at_x(Curve_point_2(nt_traits.convert
                                                        (x_asympt + 20), 0),
                                          points);
            CGAL_assertion (n_points == 1);
            Curve_point_2 right_2 = points[0];

            n_points = cv.get_points_at_x(Curve_point_2(nt_traits.convert
                                                        (x_asympt + 30), 0),
                                          points);
            CGAL_assertion (n_points == 1);
            Curve_point_2 right_3 = points[0];

            // left branch
            n_points = cv.get_points_at_x(Curve_point_2(nt_traits.convert
                                                        (x_asympt - 10), 0),
                                          points);
            CGAL_assertion (n_points == 1);
            Curve_point_2 left_1 = points[0];

            n_points = cv.get_points_at_x(Curve_point_2(nt_traits.convert
                                                        (x_asympt - 20), 0),
                                          points);
            CGAL_assertion (n_points == 1);
            Curve_point_2 left_2 = points[0];

            // checking the orientation
            if (cv.orientation() == CGAL::orientation(right_1, right_2, right_3))
            {
              *oi++ = make_object(X_monotone_curve_2(cv, right_1, right_2,
                                                     conic_id, true, true));
              *oi++ = make_object(X_monotone_curve_2(cv, left_1, left_2,
                                                     conic_id, true, true));
            }
            else {
              *oi++ = make_object(X_monotone_curve_2(cv, right_2, right_1,
                                                     conic_id, true, true));
              *oi++ = make_object(X_monotone_curve_2(cv, left_2, left_1,
                                                     conic_id, true, true));
            }
          }
        }
      }
      else {
        // Degenerate hyperbola is very similar to degenerate parabola.
        // Degenerate hyperbola could be only intersecting lines. We diffrentiating
        // between a degenerate hyperbola with a vertical line (vertical asymptot)
        // and a degenerate hyperbola with no veritcal asymptot.
        CGAL_assertion(cv.orientation() == COLLINEAR);

        Curve_point_2 left_1, right_1, left_2, right_2;
        if (cv.has_vertical_asymptote() == true) {
          Rat_line_2 asympt = cv.get_vertical_asymptote();
          Rational x_asympt = -asympt.c() / asympt.a();
          left_1 = Curve_point_2(x_asympt, 0);
          right_1 = Curve_point_2(x_asympt, 1);

          Curve_point_2 points[2];
          int n_points = cv.get_points_at_x(Curve_point_2
                                            (nt_traits.convert(x_asympt - 10),
                                             0), points);
          CGAL_assertion (n_points == 1);
          left_2 = points[0];

          n_points = cv.get_points_at_x(Curve_point_2
                                        (nt_traits.convert(x_asympt + 10),
                                         0), points);
          CGAL_assertion (n_points == 1);
          right_2 = points[0];
        }
        else {
          // We act as in the degenerate parabola case.
          Curve_point_2 points[2];
          int n_points = cv.get_points_at_x(Curve_point_2(0, 0), points);
          if (n_points < 2)
            n_points = cv.get_points_at_x(Curve_point_2(-10, 0), points);
          CGAL_assertion (n_points == 2);
          left_1 = points[0]; left_2 = points[1];

          n_points = cv.get_points_at_x(Curve_point_2(10, 0), points);
          if (n_points < 2)
            n_points = cv.get_points_at_x(Curve_point_2(20, 0), points);
          CGAL_assertion (n_points == 2);
          right_1 = points[1]; right_2 = points[0];
        }
        *oi++ = make_object (X_monotone_curve_2 (left_1, right_1, true, true));
        *oi++ = make_object (X_monotone_curve_2 (left_2, right_2, true, true));
      }
      return oi;
    }
  };

  /*! Get a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object ()
  { return Make_x_monotone_2(); }

  class Split_2 {
  public:
    /*! Split a given x-monotone curve at a given point into two sub-curves.
     * \param cv The curve to split
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve (p is its right endpoint).
     * \param c2 Output: The right resulting subcurve (p is its left endpoint).
     * \pre p lies on cv but is not one of its end-points.
     */
    void operator()(const X_monotone_curve_2& cv, const Point_2& p,
                    X_monotone_curve_2& c1, X_monotone_curve_2& c2) const {
      DCT_TRACE("Split_2:" << endl);
      DCT_TRACE("Curve: " << cv << endl);
      DCT_TRACE("Point: " << p << endl);
      cv.split (p, c1, c2);
      DCT_TRACE("Result1: " << c1 << endl);
      DCT_TRACE("Result2: " << c2 << endl);
    }
  };

  /*! Get a Split_2 functor object. */
  Split_2 split_2_object() const { return Split_2(); }

  class Intersect_2 {
  private:
    Intersection_map& m_inter_map;       // The map of intersection points.

  public:
    /*! Constructor. */
    Intersect_2(Intersection_map& map) : m_inter_map(map) {}

    /*! Find the intersections of the two given curves and insert them to the
     * given output iterator. As two segments may itersect only once, only a
     * single will be contained in the iterator.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template <typename OutputIterator>
      OutputIterator operator()(const X_monotone_curve_2& cv1,
                                const X_monotone_curve_2& cv2,
                                OutputIterator oi) {
      DCT_TRACE("Intersect_2:" << endl);
      DCT_TRACE("Curve 1: " << cv1 << endl);
      DCT_TRACE("Curve 2: " << cv2 << endl);
      return (cv1.intersect (cv2, m_inter_map, oi));
    }
  };

  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object()
  {
    // TODO: if mutable - change to const function
    return (Intersect_2 (inter_map));
  }

  class Are_mergeable_2 {
  public:
    /*! Check whether it is possible to merge two given x-monotone curves.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are mergeable - if they are supported
     *         by the same line and share a common endpoint; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const {
      DCT_TRACE("Are_mergeable_2!!!!1" << endl);
      return (cv1.can_merge_with (cv2));
    }
  };

  /*! Get an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object () const { return Are_mergeable_2(); }

  class Merge_2 {
  public:
    /*! Merge two given x-monotone curves into a single curve (segment).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable, that is they are supported by the
     *      same conic curve and share a common endpoint.
     */
    void operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2,
                    X_monotone_curve_2& c) const {
      DCT_TRACE("Merge_2!!!!1" << endl);
      c = cv1;
      c.merge (cv2);
    }
  };

  /*! Get a Merge_2 functor object. */
  Merge_2 merge_2_object() const { return Merge_2(); }

  //@}

  /// \name Functor definitions for the landmarks point-location strategy.
  //@{
  using Approximate_number_type = double;

  class Approximate_2 {
  public:
    /*! Return an approximation of a point coordinate.
     * \param p The exact point.
     * \param i The coordinate index (either 0 or 1).
     * \pre i is either 0 or 1.
     * \return An approximation of p's x-coordinate (if i == 0), or an
     *         approximation of p's y-coordinate (if i == 1).
     */
    Approximate_number_type operator()(const Point_2& p, int i) const {
      CGAL_precondition (i == 0 || i == 1);
      if (i == 0) return (CGAL::to_double(p.x()));
      else return (CGAL::to_double(p.y()));
    }
  };

  /*! Get an Approximate_2 functor object. */
  Approximate_2 approximate_2_object() const { return Approximate_2(); }

  class Construct_x_monotone_curve_2 {
  public:
    /*! Return an x-monotone curve connecting the two given endpoints.
     * \param p The first point.
     * \param q The second point.
     * \pre p and q must not be the same.
     * \return A segment connecting p and q.
     */
    X_monotone_curve_2 operator()(const Point_2& p, const Point_2& q) const
    { return (X_monotone_curve_2(p, q)); }
  };

  /*! Get a Construct_x_monotone_curve_2 functor object. */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object () const
  { return Construct_x_monotone_curve_2(); }
  //@}

  /// \name Functor definitions for the Boolean set-operation traits.
  //@{

  class Compare_endpoints_xy_2 {
  public:
    /*! Compare the endpoints of an $x$-monotone curve lexicographically.
     * (assuming the curve has a designated source and target points).
     * \param cv The curve.
     * \return SMALLER if the curve is directed right;
     *         LARGER if the curve is directed left.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv) {
      if (cv.is_directed_right()) return (SMALLER);
      else return (LARGER);
    }
  };

  /*! Get a Compare_endpoints_xy_2 functor object. */
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(); }

  class Construct_opposite_2 {
  public:
    /*! Construct an opposite x-monotone (with swapped source and target).
     * \param cv The curve.
     * \return The opposite curve.
     */
    X_monotone_curve_2 operator() (const X_monotone_curve_2& cv)
    { return (cv.flip()); }
  };

  /*! Get a Construct_opposite_2 functor object. */
  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(); }
  //@}
};

} //namespace CGAL

#endif
