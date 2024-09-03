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
// $URL:
// $Id:
//
// Author(s)     : Ophir Setter          <ophirset@post.tau.ac.il>

#ifndef CGAL_ARR_CIRCLE_LINEAR_TRAITS_2_H
#define CGAL_ARR_CIRCLE_LINEAR_TRAITS_2_H

/*! \file
 * The header file for the Arr_circle_linear_traits_2<Kenrel> class.
 */

#include <CGAL/tags.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_geometry_traits/Circle_linear_2.h>

#include <fstream>

namespace CGAL {

/*! \class
 * A traits class for maintaining an arrangement of circles.
 */
template <typename Kernel_, bool Filter_ = true>
class Arr_circle_linear_traits_2 {
public:
  using Kernel = Kernel_;
  using NT = typename Kernel::FT;
  using Point_2 = _One_root_point_2<NT, Filter_>;
  using CoordNT = typename Point_2::CoordNT;
  using Curve_2 = _Circle_linear_2<Kernel, Filter_>;
  using X_monotone_curve_2 = _X_monotone_circle_linear_2<Kernel, Filter_>;
  using Multiplicity = unsigned int;
  using Self = Arr_circle_linear_traits_2<Kernel, Filter_>;

  // Category tags:
  using Has_left_category = Tag_true;
  using Has_merge_category = Tag_true;
  using Has_do_intersect_category = Tag_false;

  using Left_side_category = Arr_open_side_tag;
  using Bottom_side_category = Arr_open_side_tag;
  using Top_side_category = Arr_open_side_tag;
  using Right_side_category = Arr_open_side_tag;

protected:
  // Type definition for the intersection points mapping.
  using Intersection_map = typename X_monotone_curve_2::Intersection_map;

  mutable Intersection_map  inter_map;   // Mapping pairs of curve IDs to their
                                 // intersection points.
  bool m_use_cache;

public:
  /*! Default constructor. */
  Arr_circle_linear_traits_2(bool use_intersection_caching = false) :
    m_use_cache(use_intersection_caching)
  {}

  /*! Get the next curve index. */
  static unsigned int get_index() {
    static unsigned int index = 0;
    return ++index;
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
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const {
      if (p1.identical(p2)) return (EQUAL);
      return CGAL::compare(p1.x(), p2.x());
    }
  };

  /*! Get a Compare_x_2 functor object. */
  Compare_x_2 compare_x_2_object() const { return Compare_x_2(); }

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
      if (p1.identical(p2)) return (EQUAL);
      Comparison_result res = CGAL::compare(p1.x(), p2.x());
      if (res != EQUAL) return (res);
      return (CGAL::compare(p1.y(), p2.y()));
    }
  };

  /*! Get a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object() const { return Compare_xy_2(); }

  class Construct_min_vertex_2 {
  public:
    /*! Get the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The left endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2 & cv) const {
      CGAL_precondition(cv.has_left());
      return cv.left();
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(); }

  class Construct_max_vertex_2 {
  public:
    /*! Get the right endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The right endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2& cv) const {
      CGAL_precondition(cv.has_right());
      return cv.right();
    }
  };

  /*! Get a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(); }

  class Is_vertical_2 {
  public:
    /*! Check whether the given x-monotone curve is a vertical segment.
     * \param cv The curve.
     * \return (true) if the curve is a vertical segment; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv) const
    { return (cv.is_vertical()); }
  };

  /*! Get an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object() const
  { return Is_vertical_2(); }

  class Compare_y_at_x_2 {
  public:
    /*! Return the location of the given point with respect to the input curve.
     * \param cv The curve.
     * \param p The point.
     * \pre p is in the x-range of cv.
     * \return SMALLER if y(p) < cv(x(p)), i.e. the point is below the curve;
     *         LARGER if y(p) > cv(x(p)), i.e. the point is above the curve;
     *         EQUAL if p lies on the curve.
     */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& cv) const {
      CGAL_precondition(cv.is_in_x_range(p));
      return cv.point_position(p);
    }
  };

  /*! Get a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(); }

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
      // Make sure that p lies on both curves, and that both are defined to its
      // right (so their right endpoint is lexicographically larger than p).
      CGAL_precondition((cv1.point_position(p) == EQUAL) &&
                        (cv2.point_position(p) == EQUAL));

      // Compare the two curves immediately to the right of p:
      return cv1.compare_to_right(cv2, p);
    }
  };

  /*! Get a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(); }

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
      // Make sure that p lies on both curves, and that both are defined to its
      // left (so their left endpoint is lexicographically smaller than p).

      CGAL_precondition(cv1.point_position(p) == EQUAL &&
                        cv2.point_position(p) == EQUAL);

      if ((CGAL::compare(cv1.left().x(),cv1.right().x()) == EQUAL) &&
          (CGAL::compare(cv2.left().x(),cv2.right().x()) == EQUAL))
      { //both cv1 and cv2 are vertical
         CGAL_precondition(!(cv1.left()).equals(p) && !(cv2.left()).equals(p));
      }
      else if ((CGAL::compare(cv1.left().x(),cv1.right().x()) != EQUAL) &&
               (CGAL::compare(cv2.left().x(),cv2.right().x()) == EQUAL))
      { //only cv1 is vertical
        CGAL_precondition(!(cv1.left()).equals(p));
      }
      else if ((CGAL::compare(cv1.left().x(),cv1.right().x()) == EQUAL) &&
               (CGAL::compare(cv2.left().x(),cv2.right().x()) != EQUAL))
      { //only cv2 is vertical
        CGAL_precondition(!(cv2.left()).equals(p));
      }
      else { //both cv1 and cv2 are non vertical
        CGAL_precondition(CGAL::compare(cv1.left().x(),p.x()) == SMALLER &&
                          CGAL::compare(cv2.left().x(),p.x()) == SMALLER);
      }
      // Compare the two curves immediately to the left of p:
      return (cv1.compare_to_left(cv2, p));
    }
  };

  /*! Get a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(); }

  class Equal_2 {
  public:
    /*! Check if the two x-monotone curves are the same (have the same graph).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are the same; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const {
      if (&cv1 == &cv2) return true;
      return cv1.equals(cv2);
    }

    /*!
     * Check if the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    bool operator()(const Point_2& p1, const Point_2& p2) const
    { return (p1.equals(p2)); }
  };

  /*! Get an Equal_2 functor object. */
  Equal_2 equal_2_object() const { return Equal_2(); }
  //@}

  /// \name Functor definitions to handle boundaries
  //@{

  /*! A function object that determines whether a curve end is bounded.
   */
  class Is_bounded_2 {
  public:
    /*! Is the end of an x-monotone curve bounded?
     * \param xcv The x-monotone curve.
     * \param ce The end of xcv identifier.
     * \return true is the curve end is bounded, and false otherwise
     */
    bool operator()(const X_monotone_curve_2 & xcv, Arr_curve_end ce) const {
      if (xcv.is_circular()) return true;
      return (ce == ARR_MIN_END) ? (xcv.has_left()) : (xcv.has_right());
    }
  };

  /*! Obtain a Is_bounded_2 function object. */
  Is_bounded_2 is_bounded_2_object() const { return Is_bounded_2(); }

  /*! A function object that obtains the parameter space of a geometric
   * entity along the x-axis
   */
  class Parameter_space_in_x_2 {
  public:
    /*! Obtains the parameter space at the end of a line along the x-axis.
     * \param xcv the xcurve.
     * \param ce the line end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the line xcv.
     *   ARR_LEFT_BOUNDARY  - the line approaches the identification arc from
     *                        the right at the line left end.
     *   ARR_INTERIOR       - the line does not approache the identification arc.
     *   ARR_RIGHT_BOUNDARY - the line approaches the identification arc from
     *                        the left at the line right end.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & xcv,
                                   Arr_curve_end ce) const {
      if (xcv.is_circular()) return ARR_INTERIOR;
      if (sign(xcv.b()) == ZERO) return ARR_INTERIOR;
      if (ce == ARR_MIN_END)
        return xcv.has_left() ? ARR_INTERIOR : ARR_LEFT_BOUNDARY;
      else return xcv.has_right() ? ARR_INTERIOR : ARR_RIGHT_BOUNDARY;
    }

    /*! Obtains the parameter space at a point along the x-axis.
     * \param p the point.
     * \return the parameter space at p.
     */
    Arr_parameter_space operator()(const Point_2 ) const { return ARR_INTERIOR; }
  };

  /*! Obtain a Parameter_space_in_x_2 function object */
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const
  { return Parameter_space_in_x_2(); }

  /*! A function object that obtains the parameter space of a geometric
   * entity along the y-axis
   */
  class Parameter_space_in_y_2 {
  public:
    /*! Obtains the parameter space at the end of a line along the y-axis .
     * Note that if the line end coincides with a pole, then unless the line
     * coincides with the identification arc, the line end is considered to
     * be approaching the boundary, but not on the boundary.
     * If the line coincides with the identification arc, it is assumed to
     * be smaller than any other object.
     * \param xcv the xcurve.
     * \param ce the line end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the line xcv.
     *   ARR_BOTTOM_BOUNDARY  - the line approaches the south pole at the line
     *                          left end.
     *   ARR_INTERIOR         - the line does not approache a contraction point.
     *   ARR_TOP_BOUNDARY     - the line approaches the north pole at the line
     *                          right end.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & xcv,
                                   Arr_curve_end ce) const {
      if (xcv.is_circular()) return ARR_INTERIOR;
      if (sign(xcv.a()) == ZERO) return ARR_INTERIOR;
      if (ce == ARR_MIN_END)
        return xcv.has_left() ? ARR_INTERIOR : ARR_BOTTOM_BOUNDARY;
      else return xcv.has_right() ? ARR_INTERIOR : ARR_TOP_BOUNDARY;
    }

    /*! Obtains the parameter space at a point along the y-axis.
     * \param p the point.
     * \return the parameter space at p.
     */
    Arr_parameter_space operator()(const Point_2 ) const { return ARR_INTERIOR; }
  };

  /*! Obtain a Parameter_space_in_y_2 function object */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(); }

  /*! A function object that compares the x-limits of arc ends on the
   * boundary of the parameter space.
   */
  class Compare_x_on_boundary_2 {
  protected:
    typedef Arr_circle_linear_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_on_boundary_2(const Traits& traits) : m_traits(traits) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_circle_linear_traits_2<Kernel>;

  public:
    /*! Compare the x-limit of a vertical line at a point with the x-limit of
     * a line end near the boundary at y = +/- oo.
     * \param p the point direction.
     * \param xcv the line, the endpoint of which is compared.
     * \param ce the line-end indicator -
     *            ARR_MIN_END - the minimal end of xc or
     *            ARR_MAX_END - the maximal end of xc.
     * \return the comparison result:
     *         SMALLER - x(p) < x(xc, ce);
     *         EQUAL   - x(p) = x(xc, ce);
     *         LARGER  - x(p) > x(xc, ce).
     * \pre p lies in the interior of the parameter space.
     * \pre the ce end of the line xcv lies on a boundary.
     */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& xcv,
                                 Arr_curve_end) const {
      CGAL_precondition(xcv.is_linear());
      CGAL_precondition(xcv.is_vertical());

      Kernel kernel;
      typename Kernel::Line_2 line = xcv.supporting_line();
      typename Kernel::FT a = kernel.compute_a_2_object()(line);
      typename Kernel::FT c = kernel.compute_c_2_object()(line);
      return CGAL::compare(p.x(), -c/a);
      //return (kernel.compare_x_at_y_2_object()(p, xcv.supporting_line()));
    }

    /*! Compare the x-limits of 2 arcs ends near the boundary of the
     * parameter space at y = +/- oo.
     * \param xcv1 the first arc.
     * \param ce1 the first arc end indicator -
     *            ARR_MIN_END - the minimal end of xcv1 or
     *            ARR_MAX_END - the maximal end of xcv1.
     * \param xcv2 the second arc.
     * \param ce2 the second arc end indicator -
     *            ARR_MIN_END - the minimal end of xcv2 or
     *            ARR_MAX_END - the maximal end of xcv2.
     * \return the second comparison result:
     *         SMALLER - x(xcv1, ce1) < x(xcv2, ce2);
     *         EQUAL   - x(xcv1, ce1) = x(xcv2, ce2);
     *         LARGER  - x(xcv1, ce1) > x(xcv2, ce2).
     * \pre the ce1 end of the line xcv1 lies on a boundary.
     * \pre the ce2 end of the line xcv2 lies on a boundary.
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 Arr_curve_end /* ce1 */,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end /*! ce2 */) const {
      CGAL_precondition(xcv1.is_linear());
      CGAL_precondition(xcv2.is_linear());
      CGAL_precondition(xcv1.is_vertical());
      CGAL_precondition(xcv2.is_vertical());

      Kernel kernel;
      typename Kernel::Point_2 p = kernel.construct_point_2_object()(ORIGIN);
      return (kernel.compare_x_at_y_2_object()(p,
                                               xcv1.supporting_line(),
                                               xcv2.supporting_line()));
    }
  };

  /*! Obtain a Compare_x_on_boundary_2 function object */
  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const
  { return Compare_x_on_boundary_2(*this); }

  /*! A function object that compares the x-coordinates of arc ends near the
   * boundary of the parameter space
   */
  class Compare_x_near_boundary_2 {
  public:
    /*! Compare the x-coordinates of 2 arcs ends near the boundary of the
     * parameter space at y = +/- oo.
     * \param xcv1 the first arc.
     * \param ce1 the first arc end indicator -
     *            ARR_MIN_END - the minimal end of xcv1 or
     *            ARR_MAX_END - the maximal end of xcv1.
     * \param xcv2 the second arc.
     * \param ce2 the second arc end indicator -
     *            ARR_MIN_END - the minimal end of xcv2 or
     *            ARR_MAX_END - the maximal end of xcv2.
     * \return the second comparison result:
     *         SMALLER - x(xcv1, ce1) < x(xcv2, ce2);
     *         EQUAL   - x(xcv1, ce1) = x(xcv2, ce2);
     *         LARGER  - x(xcv1, ce1) > x(xcv2, ce2).
     * \pre the ce end of the line xcv1 lies on a boundary, implying
     *      that xcv1 is vertical.
     * \pre the ce end of the line xcv2 lies on a boundary, implying
     *      that xcv2 is vertical.
     * \pre the $x$-coordinates of xcv1 and xcv2 at their ce ends are
     *      equal, implying that the curves overlap!
     */
    Comparison_result
    operator()(const X_monotone_curve_2& CGAL_precondition_code(xcv1),
               const X_monotone_curve_2& CGAL_precondition_code(xcv2),
               Arr_curve_end /*! ce2 */) const {
      CGAL_precondition(xcv1.is_linear());
      CGAL_precondition(xcv2.is_linear());
      CGAL_precondition(xcv1.is_vertical());
      CGAL_precondition(xcv2.is_vertical());
      return EQUAL;
    }
  };

  /*! Obtain a Compare_x_near_boundary_2 function object */
  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const
  { return Compare_x_near_boundary_2(); }

  /*! A function object that compares the y-coordinates of arc ends near the
   * boundary of the parameter space.
   */
  class Compare_y_near_boundary_2 {
  public:
    /*! Compare the y-coordinates of 2 lines at their ends near the boundary
     * of the parameter space at x = +/- oo.
     * \param xcv1 the first arc.
     * \param xcv2 the second arc.
     * \param ce the line end indicator.
     * \return the second comparison result.
     * \pre the ce ends of the lines xcv1 and xcv2 lie either on the left
     * boundary or on the right boundary of the parameter space.
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end ce) const {
      // Make sure both curves are defined at x = -oo (or at x = +oo).
      CGAL_precondition ((ce == ARR_MIN_END &&
                          ! xcv1.has_left() && ! xcv2.has_left()) ||
                         (ce == ARR_MAX_END &&
                          ! xcv1.has_right() && ! xcv2.has_right()));

      // Compare the slopes of the two supporting lines.
      Kernel kernel;
      const Comparison_result res_slopes =
        kernel.compare_slope_2_object()(xcv1.supporting_line(),
                                        xcv2.supporting_line());

      if (res_slopes == EQUAL) {
        // In case the two supporting line are parallel, compare their
        // relative position at x = 0, which is the same as their position
        // at infinity.
        typename Kernel::Point_2 p = kernel.construct_point_2_object()(ORIGIN);
        return (kernel.compare_y_at_x_2_object()(p,
                                                 xcv1.supporting_line(),
                                                 xcv2.supporting_line()));
      }

      if (ce == ARR_MIN_END)
        // Flip the slope result if we compare at x = -oo:
        return (res_slopes == LARGER) ? SMALLER : LARGER;

      // If we compare at x = +oo, the slope result is what we need:
      return res_slopes;
    }
  };

  /*! Obtain a Compare_y_near_boundary_2 function object */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const
  { return Compare_y_near_boundary_2(); }

  //@}

  /// \name Functor definitions for supporting intersections.
  //@{

  class Make_x_monotone_2 {
  private:
    using Self = Arr_circle_linear_traits_2<Kernel_, Filter_>;

    bool m_use_cache;

  public:
    Make_x_monotone_2(bool use_cache = false) : m_use_cache(use_cache) {}

    /*! Cut the given conic curve (ocv.is_in_x_range (p)r conic arc) into
     * x-monotone subcurves and insert them to the given output iterator.
     * \param cv The curve.
     * \param oi The output iterator, whose value-type is Object. The returned
     *           objects are all wrcv.is_in_x_range (p)appers
     *           X_monotone_curve_2 objects.
     * \return The past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const {
      using Base = typename Curve_2::Base;

      // Increment the serial number of the curve cv, which will serve as its
      // unique identifier.
      unsigned int index = 0;
      if (m_use_cache) index = Self::get_index();

      if (cv.orientation() == COLLINEAR) {
        if (cv.has_source() && cv.has_target()) {
          // The curve is a line segment.
          *oi++ = make_object(X_monotone_curve_2(cv.supporting_line(),
                                                 cv.source(), cv.target(),
                                                 index));
        }
        else if (cv.has_source() && !cv.has_target()) {
          // The curve is a ray.
          *oi++ = make_object(X_monotone_curve_2(cv.supporting_line(),
                                                 cv.source(), index));
        }
        else {
          // The curve is a line.
          *oi++ = make_object(X_monotone_curve_2(cv.supporting_line(), index));
        }
        return oi;
      }

      // Check the case of a degenrate circle (a point).
      const typename Kernel::Circle_2&  circ = cv.supporting_circle();
      CGAL::Sign   sign_rad = CGAL::sign(circ.squared_radius());
      CGAL_precondition(sign_rad != NEGATIVE);

      if (sign_rad == ZERO) {
        // Create an isolated point.
        *oi++ = make_object(Point_2(circ.center().x(), circ.center().y()));
        return oi;
      }

      // The curve is circular: compute the to vertical tangency points
      // of the supporting circle.
      Point_2 vpts[2];
      unsigned int n_vpts = cv.vertical_tangency_points(vpts);

      if (cv.is_full()) {
        CGAL_assertion(n_vpts == 2);

        // Subdivide the circle into two arcs (an upper and a lower half).
        *oi++ = make_object(X_monotone_curve_2(circ, vpts[0], vpts[1],
                                               cv.orientation(), index));
        *oi++ = make_object(X_monotone_curve_2(circ, vpts[1], vpts[0],
                                               cv.orientation(), index));
      }
      else {
        // Act according to the number of vertical tangency points.
        if (n_vpts == 2) {
          // Subdivide the circular arc into three x-monotone arcs.
          *oi++ = make_object(X_monotone_curve_2(circ, cv.source(), vpts[0],
                                                 cv.orientation(), index));
          *oi++ = make_object(X_monotone_curve_2(circ, vpts[0], vpts[1],
                                                 cv.orientation(),
                                                 index));
          *oi++ = make_object(X_monotone_curve_2(circ, vpts[1], cv.target(),
                                                 cv.orientation(), index));
        }
        else if (n_vpts == 1) {
          // Subdivide the circular arc into two x-monotone arcs.
          *oi++ = make_object(X_monotone_curve_2(circ, cv.source(), vpts[0],
                                                 cv.orientation(), index));
          *oi++ = make_object(X_monotone_curve_2(circ, vpts[0], cv.target(),
                                                 cv.orientation(), index));
        }
        else {
          CGAL_assertion(n_vpts == 0);

          // The arc is already x-monotone:
          *oi++ = make_object(X_monotone_curve_2(circ, cv.source(), cv.target(),
                                                 cv.orientation(), index));
        }
      }

      return oi;
    }
  };

  /*! Get a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(m_use_cache); }

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
      CGAL_precondition(cv.point_position(p)==EQUAL &&
                        (! cv.has_source() || ! p.equals(cv.source())) &&
                        (! cv.has_target() || ! p.equals(cv.target())));

      cv.split(p, c1, c2);
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
                              OutputIterator oi) const
    { return (cv1.intersect(cv2, oi, &m_inter_map)); }
  };

  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object() const { return (Intersect_2(inter_map)); }

  class Are_mergeable_2 {
  public:
    /*! Check whether it is possible to merge two given x-monotone curves.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are mergeable - if they are supported
     *         by the same line and share a common endpoint; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const
    { return (cv1.can_merge_with(cv2)); }
  };

  /*! Get an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object() const { return Are_mergeable_2(); }

  class Merge_2 {
  public:
    /*! Merge two given x-monotone curves into a single curve.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable, that is they are supported by the
     *      same conic curve and share a common endpoint.
     */
    void operator()(const X_monotone_curve_2& cv1, const X_monotone_curve_2& cv2,
                    X_monotone_curve_2& c) const {
      c = cv1;
      c.merge(cv2);
    }
  };

  /*! Get a Merge_2 functor object. */
  Merge_2 merge_2_object() const { return Merge_2(); }

  class Compare_endpoints_xy_2 {
  public:
    /*! compare lexicogrphic the endpoints of a x-monotone curve.
     * \param cv the curve
     * \return SMALLER if the curve is directed right, else return SMALLER
     */
    Comparison_result operator()(const X_monotone_curve_2& cv) const {
      if (cv.is_directed_right()) return SMALLER;
      return LARGER;
    }
  };

  /*! Get a Compare_endpoints_xy_2 functor object. */
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(); }

  class Construct_opposite_2 {
  public:
    /*! construct an opposite x-monotone curve.
     * \param cv the curve
     * \return an opposite x-monotone curve.
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2& cv) const
    { return cv.construct_opposite(); }
  };

  /*! Get a Construct_opposite_2 functor object. */
  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(); }

};

} //namespace CGAL

#endif
