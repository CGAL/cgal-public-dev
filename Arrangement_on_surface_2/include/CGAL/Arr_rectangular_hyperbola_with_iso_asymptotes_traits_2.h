// Copyright (c) 2011  Tel-Aviv University (Israel).
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
// $Id: $
// 
//
// Author(s): Asaf Porat          <asafpor@post.tau.ac.il>
//            Efi Fogel           <efif@post.tau.ac.il>

#ifndef CGAL_ARR_RECTANGULAR_HYPERBOLA_WITH_ISO_ASYMPTOTES_TRAITS_2_H
#define CGAL_ARR_RECTANGULAR_HYPERBOLA_WITH_ISO_ASYMPTOTES_TRAITS_2_H

/*! \file
 * The traits-class for handling rectangular hyperbolas with vertical and
 * horizontal asymptotes.
 */

#include <fstream>

#include <CGAL/tags.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_geometry_traits/Sqrt_extension_point_2.h>
#include <CGAL/Arr_geometry_traits/Rectangular_hyperbola_with_iso_asymptotes_2.h>

namespace CGAL {

/*! \class
 * A traits class for maintaining an arrangement of rectangular hyperbolas
 * with vertical and horiznotal asysmptotes.
 */
template <typename Kernel_, typename Filter_ = true>
class Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2 {
public:
  typedef Kernel_                                     Kernel;
  typedef Filter_                                     Filter;
  
  typedef typename Kernel::FT                         NT;
  typedef typename Kernel::Point_2                    Rational_point_2;
  typedef typename Kernel::Segment_2                  Rational_segment_2;
  typedef Sqrt_extension_point_2<NT, Filter>          Point_2;
  typedef typename Point_2::Coord_NT                  Coord_NT;
  typedef X_monotone_rectangular_hyperbola_with_iso_asymptotes_2<Kernel, Filter>
                                                      X_monotone_curve_2;
  typedef Rectangular_hyperbola_with_iso_asymptotes_2<Kernel, Filter>
                                                      Curve_2;

  // Category tags:
  typedef Tag_true                                    Has_left_category;
  typedef Tag_true                                    Has_merge_category;
  typedef Tag_false                                   Has_do_intersect_category;

  typedef Arr_open_side_tag                           Arr_left_side_category;
  typedef Arr_open_side_tag                           Arr_bottom_side_category;
  typedef Arr_open_side_tag                           Arr_top_side_category;
  typedef Arr_open_side_tag                           Arr_right_side_category;
  
  typedef typename Kernel::Line_2                     Line_2;
  typedef typename Kernel::Segment_2                  Segment_2;

public:

  /*!
   * Default constructor.
   */
  Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2() {}

  /// \name Basic functor definitions.
  //@{

  /*! A functor that compares the x-coordinates of two points */
  class Compare_x_2 {
  public:
    /*!
     * Compare the x-coordinates of two points.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2);
     *         SMALLER if x(p1) < x(p2);
     *         EQUAL if x(p1) = x(p2).
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      if (p1.identical(p2)) return EQUAL;
      return CGAL::compare(p1.x(), p2.x());
    }
  };

  /*! Obtain a Compare_x_2 functor. */
  Compare_x_2 compare_x_2_object() const { return Compare_x_2(); }

  /*! A functor that compares the x-coordinates of two points */
  class Compare_xy_2 {
  public:
    /*!
     * Compare two points lexigoraphically: by x, then by y.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2), or if x(p1) = x(p2) and y(p1) > y(p2);
     *         SMALLER if x(p1) < x(p2), or if x(p1) = x(p2) and y(p1) < y(p2);
     *         EQUAL if the two points are equal.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      if (p1.identical(p2)) return EQUAL;
      Comparison_result res = CGAL::compare(p1.x(), p2.x());
      if (res != EQUAL) return res;
      return CGAL::compare(p1.y(), p2.y());
    }
  };

  /*! Obtain a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object() const { return Compare_xy_2(); }

  /*! A functor that obtains the left endpoint of a segment or a ray. */
  class Construct_min_vertex_2 {
  public:
    /*!
     * Get the left endpoint of the x-monotone curve.
     * \param xcv The curve.
     * \pre The left end of xcv is a valid (bounded) point.
     * \return The left endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2& xcv) const
    {
      CGAL_precondition(xcv.has_right());
      return xcv.right();
    }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(); }

  /*! A functor that obtains the right endpoint of a segment or a ray. */
  class Construct_max_vertex_2 {
  public:
    /*!
     * Get the right endpoint of the x-monotone curve.
     * \param xcv The curve.
     * \pre The right end of xcv is a valid (bounded) point.
     * \return The right endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2& xcv) const
    {
      CGAL_precondition(xcv.has_left());
      return xcv.left();
    }
  };

  /*! Obtain a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(); }

  /*! A functor that checks whether a given linear curve is vertical. */
  class Is_vertical_2 {
  public:
    /*!
     * Checks whether the given x-monotone curve is a vertical segment.
     * \param xcv The curve.
     * \return (true) if the curve is a vertical segment; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& xcv) const
    { return xcv.is_vertical(); }
  };

  /*! Obtain an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object() const { return Is_vertical_2(); }

  /*! A functor that compares the y-coordinates of a point and a line at
   * the point x-coordinate
   */
  class Compare_y_at_x_2 {
  protected:
    typedef Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2<Kernel>
      Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_y_at_x_2(const Traits* traits) : m_traits(traits) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2<Kernel>;
    
  public:
    /*!
     * Return the location of the given point with respect to the input curve.
     * \param xcv The curve.
     * \param p The point.
     * \pre p is in the x-range of xcv.
     * \return SMALLER if y(p) < xcv(x(p)), i.e. the point is below the curve;
     *         LARGER if y(p) > xcv(x(p)), i.e. the point is above the curve;
     *         EQUAL if p lies on the curve.
     */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& xcv) const
    {
      // TODO
      return EQUAL;
    }
  };

  /*! Obtain a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(this); }

  /*! A functor that compares compares the y-coordinates of two linear
   * curves immediately to the left of their intersection point.
   */
  class Compare_y_at_x_left_2 {
  public:
    /*!
     * Compare the y value of two x-monotone curves immediately to the left
     * of their intersection point.
     * \param xcv1 The first curve.
     * \param xcv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its left.
     * \return The relative position of xcv1 with respect to xcv2 immdiately to
     *         the left of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 const Point_2& CGAL_precondition_code(p)) const
    {
      // TODO
      return EQUAL;
    }
  };

  /*! Obtain a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(); }

  /*! A functor that compares compares the y-coordinates of two linear
   * curves immediately to the right of their intersection point.
   */
  class Compare_y_at_x_right_2 {
  public:
    /*!
     * Compare the y value of two x-monotone curves immediately to the right
     * of their intersection point.
     * \param xcv1 The first curve.
     * \param xcv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its right.
     * \return The relative position of xcv1 with respect to xcv2 immdiately to
     *         the right of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 const Point_2& CGAL_precondition_code(p)) const
    {
      // TODO
      return EQUAL;
    }
  };

  /*! Obtain a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(); }

  /*! A functor that checks whether two points and two linear curves are
   * identical.
   */
  class Equal_2 {
  public:
    /*!
     * Check whether the two x-monotone curves are the same (have the same
     * graph).
     * \param xcv1 The first curve.
     * \param xcv2 The second curve.
     * \return (true) if the two curves are the same; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& xcv1,
                    const X_monotone_curve_2& xcv2) const
    { return (xcv1 == xcv2); }

    /*!
     * Check whether the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    bool operator()(const Point_2& p1, const Point_2& p2) const
    { return (p1 == p2); }
  };

  /*! Obtain an Equal_2 functor object. */
  Equal_2 equal_2_object() const { return Equal_2(); }
  //@}

  /// \name Functor definitions to handle boundaries
  //@{

  /*! A function object that obtains the parameter space of a geometric
   * entity along the x-axis
   */
  class Parameter_space_in_x_2 {
  public:
    /*! Obtains the parameter space at the end of a line along the x-axis.
     * \param xcv the line
     * \param ce the line end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the line xcv.
     *   ARR_LEFT_BOUNDARY  - the curve end is unbounded at the left of the
     *                        parameter space.
     *   ARR_INTERIOR       - the curve end is bounded.
     *   ARR_RIGHT_BOUNDARY - the curve end is unbounded at the right of the
     *                        parameter space.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & xcv,
                                   Arr_curve_end ce) const
    {
      return (ce == ARR_MIN_END) ?
        (cv.has_left_x() ? ARR_INTERIOR : ARR_LEFT_BOUNDARY) :
        (cv.has_right_x() ? ARR_INTERIOR : ARR_RIGHT_BOUNDARY);
    }
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
     * \param xcv the line
     * \param ce the line end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the line xcv.
     *   ARR_BOTTOM_BOUNDARY - the curve end is unbounded at the bottom of
     *                         the parameter space.     
     *   ARR_INTERIOR        - the curve end is bounded.
     *   ARR_TOP_BOUNDARY    - the curve end is unbounded at the top of
     *                         the parameter space.     
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & xcv,
                                   Arr_curve_end ce) const
    {
      return (ce == ARR_MIN_END) ?
        (xcv.has_left_y() ? ARR_INTERIOR :
         ((cv.b()*cv.c() < xcv.d()) ? ARR_BOTTOM_BOUNDARY : ARR_TOP_BOUNDARY)) :
        (xcv.has_right_y() ? ARR_INTERIOR :
         ((cv.b()*cv.c() < xcv.d()) ? ARR_BOTTOM_BOUNDARY : ARR_TOP_BOUNDARY)) :
    }
  };

  /*! Obtain a Parameter_space_in_y_2 function object */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(); }

  /*! A function object that compares the x-limits of arc ends on the
   * boundary of the parameter space
   */
  class Compare_x_at_limit_2 {
   public:
    /*! Compare the x-coordinate of a point with the x-coordinate of
     * the vertical asymptote of a hyperbola.
     * a line end on the boundary at y = +/- oo.
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
     * \pre the ce end of the curve xcv lies on a boundary.
     */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2&  xcv, 
                                 Arr_curve_end ce)
    {
      CGAL_precondition(Parameter_space_in_x_2()(xcv,ce) == ARR_INTERIOR);
      CGAL_precondition(Parameter_space_in_y_2()(xcv,ce) != ARR_INTERIOR);
      return CGAL::compare(p.x(),
                           (ce == ARR_MIN_END) ? xcv.left_x() : xcv.right_x());
    }

    /*! Compares the curve end of  xcv1 that is defined by ce1 
     * with the curve end of xcv2 that is defined by ce2
     * at their limits in x. 
     * Returns SMALLER, EQUAL, or LARGER accordingly.
     */
    Comparison_result operator()(const X_monotone_curve_2&  xcv1, 
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2&  xcv2, 
                                 Arr_curve_end ce2)
    {
      CGAL_precondition(Parameter_space_in_x_2()(xcv1,ce1) == ARR_INTERIOR);
      CGAL_precondition(Parameter_space_in_y_2()(xcv1,ce1) != ARR_INTERIOR);
      CGAL_precondition(Parameter_space_in_x_2()(xcv2,ce2) == ARR_INTERIOR);
      CGAL_precondition(Parameter_space_in_y_2()(xcv2,ce2) != ARR_INTERIOR);

      return CGAL::compare((ce1 == ARR_MIN_END) ? xcv1.left_x() : xcv1.right_x(),
                           (ce2 == ARR_MIN_END) ? xcv2.left_x() : xcv2.right_x());
    }

  };

  /*! Obtain a Compare_x_at_limit_2 function object */
  Compare_x_at_limit_2 compare_x_at_limit_2_object() const
  { return Compare_x_at_limit_2(); }

  /*! A function object that compares the x-coordinates of arc ends near the
   * boundary of the parameter space
   */
  class Compare_x_near_limit_2 {
  public:
    /*! Compare the x-coordinates of 2 arcs ends near the boundary of the
     * parameter space at y = +/- oo.
     * \param xcv1 the first arc.
     * \param xcv2 the second arc.
     * \param ce the arc end indicator -
     *           ARR_MIN_END - the minimal end of xcv2 or
     *           ARR_MAX_END - the maximal end of xcv2.
     * \return the second comparison result:
     *         SMALLER - x(xcv1, ce1) < x(xcv2, ce2);
     *         EQUAL   - x(xcv1, ce1) = x(xcv2, ce2);
     *         LARGER  - x(xcv1, ce1) > x(xcv2, ce2).
     * \pre the ce1 end of the line xcv1 lies on a boundary.
     * \pre the ce2 end of the line xcv2 lies on a boundary.
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end ce) const
    {
      return SMALLER;
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
                                 Arr_curve_end ce) const
    {
      // // Make sure both curves are defined at x = -oo (or at x = +oo).
      // CGAL_precondition (! xcv1.is_degenerate());
      // CGAL_precondition (! xcv2.is_degenerate());
      // CGAL_precondition ((ce == ARR_MIN_END &&
      //                     xcv1.left_infinite_in_x() == ARR_LEFT_BOUNDARY &&
      //                     xcv2.left_infinite_in_x() == ARR_LEFT_BOUNDARY) ||
      //                    (ce == ARR_MAX_END &&
      //                     xcv1.right_infinite_in_x() == ARR_RIGHT_BOUNDARY &&
      //                     xcv2.right_infinite_in_x() == ARR_RIGHT_BOUNDARY));

      // // Compare the slopes of the two supporting lines.
      // Kernel                    kernel;
      // const Comparison_result   res_slopes =
      //   kernel.compare_slope_2_object() (xcv1.supp_line(), xcv2.supp_line());

      // if (res_slopes == EQUAL) {
      //   // In case the two supporting line are parallel, compare their
      //   // relative position at x = 0, which is the same as their position
      //   // at infinity.
      //   typename Kernel::Point_2 p = kernel.construct_point_2_object() (ORIGIN);
      //   return (kernel.compare_y_at_x_2_object() (p,
      //                                             xcv1.supp_line(),
      //                                             xcv2.supp_line()));
      // }

      // if (ce == ARR_MIN_END)
      //   // Flip the slope result if we compare at x = -oo:
      //   return ((res_slopes == LARGER) ? SMALLER : LARGER);

      // // If we compare at x = +oo, the slope result is what we need:
      // return (res_slopes);
       //TODO // ASAFP
    }
  };

  /*! Obtain a Compare_y_near_boundary_2 function object */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const
  { return Compare_y_near_boundary_2(); }
  
  //@}
  
  /// \name Functor definitions for supporting intersections.
  //@{

  class Make_x_monotone_2 {
  public:
    /*!
     * Cut the given curve into x-monotone subcurves and insert them into the
     * given output iterator. As segments are always x_monotone, only one
     * object will be contained in the iterator.
     * \param cv The curve.
     * \param oi The output iterator, whose value-type is Object. The output
     *           object is a wrapper of an X_monotone_curve_2 which is
     *           essentially the same as the input curve.
     * \return The past-the-end iterator.
     */
    template<typename OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const
    {
      // Wrap the curve with an object.
      *oi++ = make_object(cv);

      return oi;
    }
  };

  /*! Obtain a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object () const
  { return Make_x_monotone_2(); }

  class Split_2 {
  public:
    /*!
     * Split a given x-monotone curve at a given point into two sub-curves.
     * \param xcv The curve to split
     * \param p The split point.
     * \param xcv1 Output: The left resulting subcurve (p is its right endpoint).
     * \param xcv2 Output: The right resulting subcurve (p is its left endpoint).
     * \pre p lies on xcv but is not one of its end-points.
     */
    void operator()(const X_monotone_curve_2& cv, const Point_2& p,
                    X_monotone_curve_2& xcv1, X_monotone_curve_2& xcv2) const
    {
      // TODO
      return;
    }
  };

  /*! Obtain a Split_2 functor object. */
  Split_2 split_2_object() const { return Split_2(); }

  class Intersect_2 {
  public:
    /*!
     * Find the intersections of the two given curves and insert them into the
     * given output iterator. As two segments may itersect only once, only a
     * single intersection will be contained in the iterator.
     * \param xcv1 The first curve.
     * \param xcv2 The second curve.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template<typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& xcv1,
                              const X_monotone_curve_2& xcv2,
                              OutputIterator oi) const
    {
      return oi;
    }
  };

  /*! Obtain an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () const { return Intersect_2(); }

  class Are_mergeable_2 {
  public:
    /*!
     * Check whether it is possible to merge two given x-monotone curves.
     * \param xcv1 The first curve.
     * \param xcv2 The second curve.
     * \return (true) if the two curves are mergeable - if they are supported
     *         by the same line and share a common endpoint; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& xcv1,
                    const X_monotone_curve_2& xcv2) const
    {
       // TODO
      return false;
    }
  };

  /*! Obtain an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object () const { return Are_mergeable_2(); }

  class Merge_2 {
  public:
    /*!
     * Merge two given x-monotone curves into a single curve (segment).
     * \param xcv1 The first curve.
     * \param xcv2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable, that is they are supported by the
     *      same line and share a common endpoint.
     */
    void operator()(const X_monotone_curve_2& xcv1,
                    const X_monotone_curve_2& xcv2,
                    X_monotone_curve_2& c) const
    {
      // TODO
      return;
    }
  };

  /*! Obtain a Merge_2 functor object. */
  Merge_2 merge_2_object () const { return Merge_2(); }
  //@}

  /// \name Functor definitions for the landmarks point-location strategy.
  //@{
  typedef double Approximate_number_type;

  class Approximate_2 {
  public:

    /*!
     * Return an approximation of a point coordinate.
     * \param p The exact point.
     * \param i The coordinate index (either 0 or 1).
     * \pre i is either 0 or 1.
     * \return An approximation of p's x-coordinate (if i == 0), or an 
     *         approximation of p's y-coordinate (if i == 1).
     */
    Approximate_number_type operator()(const Point_2& p, int i) const
    {
      CGAL_precondition(i == 0 || i == 1);
      return (i == 0) ? CGAL::to_double(p.x()) : CGAL::to_double(p.y());
    }
  };

  /*! Obtain an Approximate_2 functor object. */
  Approximate_2 approximate_2_object () const { return Approximate_2(); }

  // No Construct_x_monotone_curve_2!
  //@}

  /// \name Functor definitions for the Boolean set-operation traits.
  //@{
  class Compare_endpoints_xy_2 {
  public:
    /*!
     * Compare the endpoints of an $x$-monotone curve lexicographically.
     * (assuming the curve has a designated source and target points).
     * \param xcv The curve.
     * \return SMALLER if the curve is directed right;
     *         LARGER if the curve is directed left.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv)
    { return (xcv.is_directed_right()) ? SMALLER : LARGER; }
  };

  /*! Obtain a Compare_endpoints_xy_2 functor object. */
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(); }

  class Construct_opposite_2 {
  public:
    /*!
     * Construct an opposite x-monotone (with swapped source and target).
     * \param xcv The curve.
     * \return The opposite curve.
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2& xcv)
    {
      X_monotone_curve_2(xcv.a(), xcv.b(), xcv.c(), xcv.d(),
                         xcv.target(), xcv.source(), 
                         xcv.has_target(), xcv.has_source(), 
                         !xcv.is_directed_right(),
                         xcv.has_right_x(), xcv.has_left_x());
    }
  };

  /*! Obtain a Construct_opposite_2 functor object. */
  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(); }
  //@}

  /// \name Constructors.
  //@{

  // A constructor object of x-monotone curves
  class Construct_x_monotone_curve_2 {
  protected:
    typedef Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2<Kernel>
      Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Construct_x_monotone_curve_2(const Traits* traits) : m_traits(traits) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2<Kernel>;

  public:
    /*!
     * Constructor from all data fields.
     * \param a The a coefficient.
     * \param b The a coefficient.
     * \param c The a coefficient.
     * \param d The a coefficient.
     * \param sourse The source point.
     * \param target The target point.
     * \param has_source Indicates whether the curve has a valid source point.
     * \param has_source Indicates whether the curve has a valid target point.
     * \param is_directed_right Indicates whether the curve is directed right.
     * \param has_left_x Indicates whether the curve has a vertical asymptote
     *        on the left.
     * \param has_right_x Indicates whether the curve has a vertical asymptote
     *        on the right.
     * \pre The two points must not be the same.
     * \pre The source is on the hyperbola.
     * \pre The target is on the hyperbola.
     * \pre The curve is continueous, which implies one of the following:
     * (1) If has_source and has_target, then
     *     the closed interval bounded by the x-coordinates of source and
     *     target does not contain -c/a (the x-coordinate of the vertical
     *     asymptotes).
     * (2) If has_source and !has_target
     *     (a) If is_directed_right
     *         (i)  If has_left_x, then source.x < -c/a = target.x.
     *         (ii) Otherwise (!has_left_x), -c/a < source.x.
     *     (b) Otherwise (!is_directed_right)
     *         (i)  If has_right_x, then source.x > -c/a = target.x
     *         (ii) Otherwise (!has_right_x), -c/a > source.x.
     * (3) Similar for !has_source and has_target
     * (4) If !has_source and !has_target
     *     (a) has_left_x != has_right_x
     *     (b) If is_directed_right,
     *         (i)   If has_left_x, then -c/a = target.x.
     *         (ii)  If has_right_x, then -c/a = source.x.
     *     (c) Otherwise (!is_directed_right),
     *         (i)   If has_right_x, then -c/a = target.x.
     *         (ii)  If has_left_x, then -c/a = source.x.
     */
    X_monotone_curve_2 operator()(const FT& a, const FT& b,
                                  const FT& c, const FT& d,
                                  const Point_2& source, const Point_2& target,
                                  bool has_source, bool has_target,
                                  bool is_directed_right,
                                  bool has_left_x = true,
                                  bool has_right_x = true)
    {
      X_monotone_curve_2 xcv =
        X_monotone_curve_2(a, b, c, d, source, target,
                           has_source, has_target, is_directed_right,
                           has_left_x, has_right_x, true);
      return xcv;
    }

    /*!
     * Constructor for a curve bounded at the source and at the target.
     * \param a The a coefficient.
     * \param b The a coefficient.
     * \param c The a coefficient.
     * \param d The a coefficient.
     * \param sourse The source point.
     * \pre a != 0.
     * \pre The two points must not be the same.
     * \pre The source is on the underlying hyperbola or the curve source
     *      has a vertical asymptote at the x-coordinate of the source.
     * \pre The target is on the underlying hyperbola or the curve target
     *      has a vertical asymptote at the x-coordinate of the source.
     * \pre The curve is continueous. That is, the open interval bounded by
     *      the x-coordinates of source and target does not contain -c/a (the
     *      x-coordinate of the vertical asymptotes).
     */
    X_monotone_curve_2 operator()(const FT& a, const FT& b,
                                  const FT& c, const FT& d,
                                  const Point_2& source, const Point_2& target)
    {
      CGAL_assertion(!is_zero(a));
      Comparison_result res = compare_xy_2_object()(source, target);
      CGAL_assertion(res != EQUAL);
      bool is_directed_right = (res == SMALLER);
      
      NT x_singular = -c/a;
      bool has_source = x_singular != source.x();
      bool has_target = x_singular != target.x();
      X_monotone_curve_2 xcv =
        X_monotone_curve_2(a, b, c, d, source, target, has_source, has_target,
                           is_directed_right, true, true, true);
      return xcv;
    }

    /*!
     * Constructor for a curve bounded at one endpoint.
     * (a) If is_directed_right
     *     (i)  If source.x < -c/a, then has_left_x <= true
     *     (ii) Otherwise (source.x > -c/a), has_left_x <= false
     * (b) Otherwise (!is_directed_right)
     *     (i)  If source.x > -c/a, then has_right_x <= true
     *     (ii) Otherwise (source.x > -c/a), has_right_x <= false
     * \param a The a coefficient.
     * \param b The a coefficient.
     * \param c The a coefficient.
     * \param d The a coefficient.
     * \param sourse The source point.
     * \param is_directed_right Indicates whether the curve is directed right.
     * \pre The source is on the hyperbola.
     * The curve must be continueous, which implies: 
     */
    X_monotone_curve_2 operator()(const FT& a, const FT& b,
                                  const FT& c, const FT& d,
                                  const Point_2& source, bool is_directed_right)
    {
      NT x;
      NT x_singular = -c/a;
      bool has_left_x, has_right_x;
      if (is_directed_right) {
        has_left_x = (source.x() < -c/a);
        has_right_x = true;
      } else {
        has_left_x = true;
        has_right_x = (source.x() > -c/a);
      }
      Point_2 target = Point_2(x_singular, 0);
      X_monotone_curve_2 xcv =
        X_monotone_curve_2(a, b, c, d, source, target, true, false,
                           is_directed_right, has_left_x, has_right_x, true);
      return xcv;
    }
  };
  
  /*! Obtain a Construct_x_monotone_curve_2 functor object. */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const
  { return Construct_x_monotone_curve_2(this); }

  // A constructor object of curves
  class Construct_curve_2 {
    /*!
     * Constructor an unbounded curve.
     * \param a The a coefficient.
     * \param b The a coefficient.
     * \param c The a coefficient.
     * \param d The a coefficient.
     * \param has_left_x Indicates whether the curve has a vertical asymptote
     *        on the left.
     * \param has_right_x Indicates whether the curve has a vertical asymptote
     *        on the right.
     */
    X_monotone_curve_2 operator()(const FT& a, const FT& b,
                                  const FT& c, const FT& d,
                                  bool has_left_x = false,
                                  bool has_right_x = false)
    {
      Point_2 source = Point_2();
      Point_2 target = Point_2();
      bool is_continuous = has_left_x || has_right_x;
      X_monotone_curve_2 xcv =
        X_monotone_curve_2(a, b, c, d, source, target, false, false, true,
                           has_left_x, has_right_x, is_continuous);
      return xcv;
    }

    /*!
     * Constructor for a curve bounded at the source and at the target.
     * \param a The a coefficient.
     * \param b The a coefficient.
     * \param c The a coefficient.
     * \param d The a coefficient.
     * \param sourse The source point.
     * \pre The two points must not be the same.
     * \pre The source is on the hyperbola.
     * \pre The target is on the hyperbola.
     */
    X_monotone_curve_2 operator()(const FT& a, const FT& b,
                                  const FT& c, const FT& d,
                                  const Point_2& source, const Point_2& target)
    {
      Comparison_result res = compare_xy_2_object()(source, target);
      CGAL_assertion(res != EQUAL);
      bool is_directed_right = (res == SMALLER);
      NT x_singular = -c/a;
      bool is_continuous =
        (is_directed_right) ? ((source.x() < -c/a) && (-c/a < target.x())) :
        ((target.x() < -c/a) && (-c/a < source.x()))
      X_monotone_curve_2 xcv =
        X_monotone_curve_2(a, b, c, d, source, target, true, true,
                           is_directed_right, true, true, is_continuous);
      return xcv;
    }

    /*!
     * Constructor for a curve bounded at one endpoint.
     * \param a The a coefficient.
     * \param b The a coefficient.
     * \param c The a coefficient.
     * \param d The a coefficient.
     * \param sourse The source point.
     * \param is_directed_right Indicates whether the curve is directed right.
     * \pre The source is on the hyperbola.
     */
    X_monotone_curve_2 operator()(const FT& a, const FT& b,
                                  const FT& c, const FT& d,
                                  const Point_2& source,
                                  bool is_directed_right,
                                  bool is_bounded = false)
    {
      NT x_singular = -c/a;
      Point_2 target = Point_2(x_singular, 0);
      bool has_left_x, has_right_x, is_continuous;
      if (is_directed_right) {
        has_right_x = true;
        if (source.x() < x_singular) {
          is_continuous = is_bounded;
          has_left_x = is_bounded
        } else {
          is_continuous = true;
          has_left_x = false;
        }
      } else {
        has_left_x = true;
        if (source.x() > x_singular) {
          is_continuous = is_bounded;
          has_right_x = is_bounded
        } else {
          is_continuous = true;
          has_right_x = false;
        }
      }
      X_monotone_curve_2 xcv =
        X_monotone_curve_2(a, b, c, d, source, target, true, false,
                           is_directed_right, has_left_x, has_right_x,
                           is_continuous);
      return xcv;
    }
  };
  
  /*! Obtain a Construct_curve_2 functor object. */
  Construct_curve_2 construct_curve_2_object() const
  { return Construct_curve_2(); }
  //@}
};

} //namespace CGAL

#endif
