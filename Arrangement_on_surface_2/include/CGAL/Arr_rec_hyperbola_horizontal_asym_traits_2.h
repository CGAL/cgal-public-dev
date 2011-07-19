// Copyright (c) 2006  Tel-Aviv University (Israel).
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
// $URL: svn+ssh://asafpor@scm.gforge.inria.fr/svn/cgal/branches/features/Lines_through_segments_3-asafpor/Arrangement_on_surface_2/include/CGAL/Arr_rec_hyperbola_horizontal_asym_traits_2.h $
// $Id: Arr_rec_hyperbola_horizontal_asym_traits_2.h 56667 2010-06-09 07:37:13Z sloriot $
// 
//
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>

#ifndef CGAL_ARR_REC_HYPERBOLA_HORIZONTAL_ASYM_TRAITS_2_H
#define CGAL_ARR_REC_HYPERBOLA_HORIZONTAL_ASYM_TRAITS_2_H

/*! \file
 * The traits-class for handling rectangular hyperbolas with vertical and horizontal asymptotes.
 * in the arrangement package.
 */

#include <CGAL/tags.h>
#include <CGAL/intersections.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_geometry_traits/Segment_assertions.h>
#include <fstream>

namespace CGAL {

template <class Kernel_> class Arr_hyperbola_object_2;

/*! \class
 * A traits class for maintaining an arrangement of linear objects (lines,
 * rays and segments), aoviding cascading of computations as much as possible.
 */
template <class Kernel_>
class Arr_rec_hyperbola_horizontal_asym_traits_2 : public Kernel_ {
   friend class Arr_rec_hyperbola_horizontal_asym_object_2<Kernel_>;
   

public:

  typedef Kernel_                         Kernel;
  typedef typename Kernel::FT             FT;

  typedef typename Algebraic_structure_traits<FT>::Is_exact 
  Has_exact_division;

  // Category tags:
  typedef Tag_true                        Has_left_category;
  typedef Tag_true                        Has_merge_category;
  typedef Tag_false                       Has_do_intersect_category;

  typedef Arr_open_side_tag               Arr_left_side_category;
  typedef Arr_open_side_tag               Arr_bottom_side_category;
  typedef Arr_open_side_tag               Arr_top_side_category;
  typedef Arr_open_side_tag               Arr_right_side_category;
  
  typedef typename Kernel::Line_2         Line_2;
  typedef typename Kernel::Segment_2      Segment_2;

  typedef CGAL::Segment_assertions<Arr_rec_hyperbola_horizontal_asym_traits_2<Kernel> >
  Segment_assertions;
  
public:

  // Traits objects
  typedef typename Kernel::Point_2                 Point_2;
  typedef Arr_hyperbola_object_2<Kernel>           X_monotone_hyperbola_2;
  typedef Arr_hyperbola_object_2<Kernel>           Hyperbola_2;
  typedef unsigned int                             Multiplicity; // TODO: Ask Efi.

public:

  /*!
   * Default constructor.
   */
  Arr_rec_hyperbola_horizontal_asym_traits_2 ()
  {}

  /// \name Basic functor definitions.
  //@{

  /*! A functor that compares the x-coordinates of two points */
  class Compare_x_2 {
  protected:
    typedef Arr_rec_hyperbola_horizontal_asym_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits * m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_2(const Traits * traits) : m_traits(traits) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_rec_hyperbola_horizontal_asym_traits_2<Kernel>;
    
  public:
    /*!
     * Compare the x-coordinates of two points.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2);
     *         SMALLER if x(p1) < x(p2);
     *         EQUAL if x(p1) = x(p2).
     */
    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      const Kernel * kernel = m_traits;
      return (kernel->compare_x_2_object()(p1, p2));
    }
  };

  /*! Obtain a Compare_x_2 functor. */
  Compare_x_2 compare_x_2_object () const
  {
    return Compare_x_2(this);
  }

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
    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      Kernel    kernel;
      return (kernel.compare_xy_2_object()(p1, p2));
    }
  };

  /*! Obtain a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object () const
  {
    return Compare_xy_2();
  }

  /*! A functor that obtains the left endpoint of a segment or a ray. */
  class Construct_min_vertex_2
  {
  public:
    /*!
     * Get the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \pre The left end of cv is a valid (bounded) point.
     * \return The left endpoint.
     */
    const Point_2& operator() (const X_monotone_curve_2& cv) const
    {
       //TODO: ASAFP
    }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object () const
  {
    return Construct_min_vertex_2();
  }

  /*! A functor that obtains the right endpoint of a segment or a ray. */
  class Construct_max_vertex_2
  {
  public:
    /*!
     * Get the right endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \pre The right end of cv is a valid (bounded) point.
     * \return The right endpoint.
     */
    const Point_2& operator() (const X_monotone_curve_2& cv) const
     {
        //TODO // ASAFP
     }
  };

  /*! Obtain a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object () const
  {
    return Construct_max_vertex_2();
  }

  /*! A functor that checks whether a given linear curve is vertical. */
  class Is_vertical_2
  {
  public:
    /*!
     * Check whether the given x-monotone curve is a vertical segment.
     * \param cv The curve.
     * \return (true) if the curve is a vertical segment; (false) otherwise.
     */
    bool operator() (const X_monotone_curve_2& cv) const
     {
        //TOOD // ASAFP
    }
  };

  /*! Obtain an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object () const
  {
    return Is_vertical_2();
  }

  /*! A functor that compares the y-coordinates of a point and a line at
   * the point x-coordinate
   */
  class Compare_y_at_x_2 {
  protected:
    typedef Arr_rec_hyperbola_horizontal_asym_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits * m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_y_at_x_2(const Traits * traits) : m_traits(traits) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_rec_hyperbola_horizontal_asym_traits_2<Kernel>;
    
  public:
    /*!
     * Return the location of the given point with respect to the input curve.
     * \param cv The curve.
     * \param p The point.
     * \pre p is in the x-range of cv.
     * \return SMALLER if y(p) < cv(x(p)), i.e. the point is below the curve;
     *         LARGER if y(p) > cv(x(p)), i.e. the point is above the curve;
     *         EQUAL if p lies on the curve.
     */
    Comparison_result operator() (const Point_2& p,
                                  const X_monotone_curve_2& cv) const
    {
      // TODO: //ASAFP

      // const Kernel * kernel = m_traits;
      // if (! cv.is_vertical())
      //   // Compare p with the segment's supporting line.
      //   return (kernel->compare_y_at_x_2_object()(p, cv.supp_line()));

      // // Compare with the vertical segment's end-points.
      // typename Kernel::Compare_y_2  compare_y = kernel->compare_y_2_object();
      // const Comparison_result res1 =
      //   cv.has_left() ? compare_y (p, cv.left()) : LARGER;
      // const Comparison_result res2 = 
      //   cv.has_right() ? compare_y (p, cv.right()) : SMALLER;
      return EQUAL;
//      return (res1 == res2) ? res1 : EQUAL;
    }
  };

  /*! Obtain a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object () const
  {
    return Compare_y_at_x_2(this);
  }

  /*! A functor that compares compares the y-coordinates of two linear
   * curves immediately to the left of their intersection point.
   */
  class Compare_y_at_x_left_2
  {
  public:
    /*!
     * Compare the y value of two x-monotone curves immediately to the left
     * of their intersection point.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its left.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the left of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2,
                                  const Point_2& CGAL_precondition_code(p)) const
    {
      // CGAL_precondition (! cv1.is_degenerate());
      // CGAL_precondition (! cv2.is_degenerate());

      // Kernel                        kernel;

      // // Make sure that p lies on both curves, and that both are defined to its
      // // left (so their left endpoint is lexicographically smaller than p).
      // CGAL_precondition_code (
      //   typename Kernel::Compare_xy_2 compare_xy = kernel.compare_xy_2_object();
      // );

      // CGAL_precondition 
      //   (Segment_assertions::_assert_is_point_on (p, cv1, 
      //                                             Has_exact_division()) &&
      //    Segment_assertions::_assert_is_point_on (p, cv2,
      //                                             Has_exact_division()));

      // CGAL_precondition ((! cv1.has_left() ||
      //                     compare_xy(cv1.left(), p) == SMALLER) &&
      //                    (! cv2.has_left() ||
      //                     compare_xy(cv2.left(), p) == SMALLER));

      // // Compare the slopes of the two segments to determine thir relative
      // // position immediately to the left of q.
      // // Notice we use the supporting lines in order to compare the slopes,
      // // and that we swap the order of the curves in order to obtain the
      // // correct result to the left of p.
      // return (kernel.compare_slope_2_object()(cv2.supp_line(), cv1.supp_line()));
       //TODO //asaFP
    }
  };

  /*! Obtain a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object () const
  {
    return Compare_y_at_x_left_2();
  }

  /*! A functor that compares compares the y-coordinates of two linear
   * curves immediately to the right of their intersection point.
   */
  class Compare_y_at_x_right_2
  {
  public:
    /*!
     * Compare the y value of two x-monotone curves immediately to the right
     * of their intersection point.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its right.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the right of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2,
                                  const Point_2& CGAL_precondition_code(p)) const
    {
      // CGAL_precondition (! cv1.is_degenerate());
      // CGAL_precondition (! cv2.is_degenerate());

      // Kernel                        kernel;

      // // Make sure that p lies on both curves, and that both are defined to its
      // // right (so their right endpoint is lexicographically larger than p).
      // CGAL_precondition_code (
      //   typename Kernel::Compare_xy_2 compare_xy = kernel.compare_xy_2_object();
      // );

      // CGAL_precondition
      //   (Segment_assertions::_assert_is_point_on (p, cv1, 
      //                                             Has_exact_division()) &&
      //    Segment_assertions::_assert_is_point_on (p, cv2,
      //                                             Has_exact_division()));

      // CGAL_precondition ((! cv1.has_right() ||
      //                     compare_xy(cv1.right(), p) == LARGER) &&
      //                    (! cv2.has_right() ||
      //                     compare_xy(cv2.right(), p) == LARGER));

      // // Compare the slopes of the two segments to determine thir relative
      // // position immediately to the left of q.
      // // Notice we use the supporting lines in order to compare the slopes.
      // return (kernel.compare_slope_2_object()(cv1.supp_line(),
      //                                         cv2.supp_line()));
       //TODO: ASAFP
    }
  };

  /*! Obtain a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object () const
  {
    return Compare_y_at_x_right_2();
  }

  /*! A functor that checks whether two points and two linear curves are
   * identical.
   */
  class Equal_2
  {
  public:
    /*!
     * Check whether the two x-monotone curves are the same (have the same
     * graph).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are the same; (false) otherwise.
     */
    bool operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2) const
    {
      // CGAL_precondition (! cv1.is_degenerate());
      // CGAL_precondition (! cv2.is_degenerate());

      // Kernel                    kernel;
      // typename Kernel::Equal_2  equal = kernel.equal_2_object();

      // // Check that the two supporting lines are the same.
      // if (! equal (cv1.supp_line(), cv2.supp_line()) &&
      //     ! equal (cv1.supp_line(), 
      //              kernel.construct_opposite_line_2_object()(cv2.supp_line())))
      // {
      //   return (false);
      // }

      // // Check that either the two left endpoints are at infinity, or they
      // // are bounded and equal.
      // if ((cv1.has_left() != cv2.has_left()) ||
      //     (cv1.has_left() && ! equal (cv1.left(), cv2.left())))
      // {
      //   return (false);
      // }

      // // Check that either the two right endpoints are at infinity, or they
      // // are bounded and equal.
      // return ((cv1.has_right() == cv2.has_right()) &&
      //         (! cv1.has_right() || equal (cv1.right(), cv2.right())));
       //TODO: ASAFP
    }

    /*!
     * Check whether the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    bool operator() (const Point_2& p1, const Point_2& p2) const
    {
      Kernel    kernel;
      return (kernel.equal_2_object()(p1, p2));
    }
  };

  /*! Obtain an Equal_2 functor object. */
  Equal_2 equal_2_object () const
  {
    return Equal_2();
  }
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
     *   ARR_LEFT_BOUNDARY  - the line approaches the identification arc from
     *                        the right at the line left end.
     *   ARR_INTERIOR       - the line does not approache the identification arc.
     *   ARR_RIGHT_BOUNDARY - the line approaches the identification arc from
     *                        the left at the line right end.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & xcv,
                                   Arr_curve_end ce) const
    {
      // CGAL_precondition (! xcv.is_degenerate());
      // return (ce == ARR_MIN_END) ?
      //   xcv.left_infinite_in_x() : xcv.right_infinite_in_x();
       //TODO ASAFP
    }

    /*! Obtains the parameter space at a point along the x-axis.
     * \param p the point.
     * \return the parameter space at p.
     */
    Arr_parameter_space operator()(const Point_2 ) const
    {
      return ARR_INTERIOR;
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
     *   ARR_BOTTOM_BOUNDARY  - the line approaches the south pole at the line
     *                          left end.
     *   ARR_INTERIOR         - the line does not approache a contraction point.
     *   ARR_TOP_BOUNDARY     - the line approaches the north pole at the line
     *                          right end.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & xcv,
                                   Arr_curve_end ce) const
    {
      // CGAL_precondition (! xcv.is_degenerate());

      // return (ce == ARR_MIN_END) ?
      //   xcv.left_infinite_in_y() : xcv.right_infinite_in_y();
       //TODO //ASAFP
    }

    /*! Obtains the parameter space at a point along the y-axis.
     * \param p the point.
     * \return the parameter space at p.
     */
    Arr_parameter_space operator()(const Point_2 ) const
    {
      return ARR_INTERIOR;
    }
  };

  /*! Obtain a Parameter_space_in_y_2 function object */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(); }

  /*! A function object that compares the x-coordinates of arc ends near the
   * boundary of the parameter space
   */
  class Compare_x_near_boundary_2 {
  public:
    /*! Compare the x-coordinate of a point with the x-coordinate of
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
    Comparison_result operator()(const Point_2 & p,
                                 const X_monotone_curve_2 & xcv,
                                 Arr_curve_end ) const
    {
      // CGAL_precondition (! xcv.is_degenerate());
      // CGAL_precondition (xcv.is_vertical());

      // Kernel                    kernel;
      // return (kernel.compare_x_at_y_2_object() (p, xcv.supp_line()));
       //TODO: //ASAFP
    }

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
     * \pre the ce1 end of the line xcv1 lies on a boundary.
     * \pre the ce2 end of the line xcv2 lies on a boundary.
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 Arr_curve_end /* ce1 */,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end /*! ce2 */) const
    {
      // CGAL_precondition (! xcv1.is_degenerate());
      // CGAL_precondition (! xcv2.is_degenerate());
      // CGAL_precondition (xcv1.is_vertical());
      // CGAL_precondition (xcv2.is_vertical());

      // Kernel                    kernel;
      // typename Kernel::Point_2 p = kernel.construct_point_2_object() (ORIGIN);
      // return (kernel.compare_x_at_y_2_object() (p,
      //                                           xcv1.supp_line(),
      //                                           xcv2.supp_line()));
       //TODO // ASAFP
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

  class Make_x_monotone_2
  {
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
    template<class OutputIterator>
    OutputIterator operator() (const Curve_2& cv, OutputIterator oi) const
    {
      // Wrap the curve with an object.
      *oi = make_object (cv);
      ++oi;

      return (oi);
    }
  };

  /*! Obtain a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object () const
  {
    return Make_x_monotone_2();
  }

  class Split_2
  {
  public:
    /*!
     * Split a given x-monotone curve at a given point into two sub-curves.
     * \param cv The curve to split
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve (p is its right endpoint).
     * \param c2 Output: The right resulting subcurve (p is its left endpoint).
     * \pre p lies on cv but is not one of its end-points.
     */
    void operator() (const X_monotone_curve_2& cv, const Point_2& p,
                     X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
      // CGAL_precondition (! cv.is_degenerate());

      // // Make sure that p lies on the interior of the curve.
      // CGAL_precondition_code (
      //   Kernel                        kernel;
      //   typename Kernel::Compare_xy_2 compare_xy = kernel.compare_xy_2_object();
      // );

      // CGAL_precondition
      //   (Segment_assertions::_assert_is_point_on (p, cv,
      //                                             Has_exact_division()) &&
      //    (! cv.has_left() || compare_xy(cv.left(), p) == SMALLER) &&
      //    (! cv.has_right() || compare_xy(cv.right(), p) == LARGER));

      // // Perform the split.
      // c1 = cv;
      // c1.set_right (p);

      // c2 = cv;
      // c2.set_left (p);
//TODO // ASAFP
      return;
    }
  };

  /*! Obtain a Split_2 functor object. */
  Split_2 split_2_object () const
  {
    return Split_2();
  }

  class Intersect_2
  {
  public:
    /*!
     * Find the intersections of the two given curves and insert them into the
     * given output iterator. As two segments may itersect only once, only a
     * single intersection will be contained in the iterator.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template<class OutputIterator>
    OutputIterator operator() (const X_monotone_curve_2& cv1,
                               const X_monotone_curve_2& cv2,
                               OutputIterator oi) const
    {
      // CGAL_precondition (! cv1.is_degenerate());
      // CGAL_precondition (! cv2.is_degenerate());

      // // Intersect the two supporting lines.
      // Kernel       kernel;
      // CGAL::Object obj = kernel.intersect_2_object()(cv1.supp_line(),
      //                                                cv2.supp_line());

      // if (obj.is_empty())
      // {
      //   // The supporting line are parallel lines and do not intersect:
      //   return (oi);
      // }

      // // Check whether we have a single intersection point.
      // const Point_2  *ip = object_cast<Point_2> (&obj);
      
      // if (ip != NULL)
      // {
      //   // Check whether the intersection point ip lies on both segments.
      //   const bool    ip_on_cv1 = cv1.is_vertical() ? cv1.is_in_y_range(*ip) :
      //                                                 cv1.is_in_x_range(*ip);

      //   if (ip_on_cv1)
      //   {
      //     const bool  ip_on_cv2 = cv2.is_vertical() ? cv2.is_in_y_range(*ip) :
      //                                                 cv2.is_in_x_range(*ip);

      //     if (ip_on_cv2)
      //     {
      //       // Create a pair representing the point with its multiplicity,
      //       // which is always 1 for line segments.
      //       std::pair<Point_2, unsigned int>   ip_mult (*ip, 1);
      //       *oi = make_object (ip_mult);
      //       oi++;
      //     }
      //   }
      //   return (oi);
      // }

      // // In this case, the two supporting lines overlap.
      // // We start with the entire cv1 curve as the overlapping subcurve,
      // // then clip it to form the true overlapping curve.
      // typename Kernel::Compare_xy_2  compare_xy = kernel.compare_xy_2_object();
      // X_monotone_curve_2             ovlp = cv1;

      // if (cv2.has_left())
      // {
      //   // If the left endpoint of cv2 is to the right of cv1's left endpoint,
      //   // clip the overlapping subcurve.
      //   if (! cv1.has_left())
      //   {
      //     ovlp.set_left (cv2.left(), false);
      //   }
      //   else
      //   {
      //     if (compare_xy (cv1.left(), cv2.left()) == SMALLER)
      //       ovlp.set_left (cv2.left(), false);
      //   }
      // }

      // if (cv2.has_right())
      // {
      //   // If the right endpoint of cv2 is to the left of cv1's right endpoint,
      //   // clip the overlapping subcurve.
      //   if (! cv1.has_right())
      //   {
      //     ovlp.set_right (cv2.right(), false);
      //   }
      //   else
      //   {
      //     if (compare_xy (cv1.right(), cv2.right()) == LARGER)
      //       ovlp.set_right (cv2.right(), false);
      //   }
      // }

      // // Examine the resulting subcurve.
      // Comparison_result        res = SMALLER;

      // if (ovlp.has_left() && ovlp.has_right())
      //   res = compare_xy (ovlp.left(), ovlp.right());

      // if (res == SMALLER)
      // {
      //   // We have discovered a true overlapping subcurve:
      //   *oi = make_object (ovlp);
      //   oi++;
      // }
      // else if (res == EQUAL)
      // {
      //   // The two objects have the same supporting line, but they just share
      //   // a common endpoint. Thus we have an intersection point, but we leave
      //   // the multiplicity of this point undefined.
      //   std::pair<Point_2, unsigned int>   ip_mult (ovlp.left(), 0);
      //   *oi = make_object (ip_mult);
      //   oi++;
      // }
       //TODO //ASAFP
      return (oi);
    }
  };

  /*! Obtain an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () const
  {
    return Intersect_2();
  }

  class Are_mergeable_2
  {
  public:
    /*!
     * Check whether it is possible to merge two given x-monotone curves.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are mergeable - if they are supported
     *         by the same line and share a common endpoint; (false) otherwise.
     */
    bool operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2) const
    {
      // CGAL_precondition (! cv1.is_degenerate());
      // CGAL_precondition (! cv2.is_degenerate());

      // Kernel                    kernel;
      // typename Kernel::Equal_2  equal = kernel.equal_2_object();

      // // Check whether the two curves have the same supporting line.
      // if (! equal (cv1.supp_line(), cv2.supp_line()) && 
      //     ! equal (cv1.supp_line(), 
      //              kernel.construct_opposite_line_2_object()(cv2.supp_line())))
      //   return (false);

      // // Check whether the left endpoint of one curve is the right endpoint of the
      // // other.
      // return ((cv1.has_right() && cv2.has_left() &&
      //          equal (cv1.right(), cv2.left())) ||
      //         (cv2.has_right() && cv1.has_left() &&
      //          equal (cv2.right(), cv1.left())));
       //TODO: ASAFP
    }
  };

  /*! Obtain an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object () const
  {
    return Are_mergeable_2();
  }

  class Merge_2
  {
  public:
    /*!
     * Merge two given x-monotone curves into a single curve (segment).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable, that is they are supported by the
     *      same line and share a common endpoint.
     */
    void operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2,
                     X_monotone_curve_2& c) const
    {
      // CGAL_precondition (! cv1.is_degenerate());
      // CGAL_precondition (! cv2.is_degenerate());

      // Kernel                    kernel;
      // typename Kernel::Equal_2  equal = kernel.equal_2_object();

      // CGAL_precondition
      //   (equal (cv1.supp_line(), 
      //           cv2.supp_line()) ||
      //    equal (cv1.supp_line(),
      //           kernel.construct_opposite_line_2_object()(cv2.supp_line())));

      // // Check which curve extends to the right of the other.
      // if (cv1.has_right() && cv2.has_left() &&
      //     equal (cv1.right(), cv2.left()))
      // {
      //   // cv2 extends cv1 to the right.
      //   c = cv1;

      //   if (cv2.has_right())
      //     c.set_right (cv2.right());
      //   else
      //     c.set_right();      // Unbounded endpoint. 
      // }
      // else
      // {
      //   CGAL_precondition (cv2.has_right() && cv1.has_left() &&
      //                      equal (cv2.right(), cv1.left()));

      //   // cv1 extends cv2 to the right.
      //   c = cv2;

      //   if (cv1.has_right())
      //     c.set_right (cv1.right());
      //   else
      //     c.set_right();      // Unbounded endpoint.
      // }
      //TODO: ASAFP
      return;
    }
  };

  /*! Obtain a Merge_2 functor object. */
  Merge_2 merge_2_object () const
  {
    return Merge_2();
  }
  //@}

  /// \name Functor definitions for the landmarks point-location strategy.
  //@{
  typedef double                          Approximate_number_type;

  class Approximate_2
  {
  public:

    /*!
     * Return an approximation of a point coordinate.
     * \param p The exact point.
     * \param i The coordinate index (either 0 or 1).
     * \pre i is either 0 or 1.
     * \return An approximation of p's x-coordinate (if i == 0), or an 
     *         approximation of p's y-coordinate (if i == 1).
     */
    Approximate_number_type operator() (const Point_2& p,
                                        int i) const
    {
      CGAL_precondition (i == 0 || i == 1);

      if (i == 0)
        return (CGAL::to_double(p.x()));
      else
        return (CGAL::to_double(p.y()));
    }
  };

  /*! Obtain an Approximate_2 functor object. */
  Approximate_2 approximate_2_object () const
  {
    return Approximate_2();
  }

  class Construct_x_monotone_curve_2
  {
  public:

    /*!
     * Return an x-monotone curve connecting the two given endpoints.
     * \param p The first point.
     * \param q The second point.
     * \pre p and q must not be the same.
     * \return A segment connecting p and q.
     */
    X_monotone_curve_2 operator() (const Point_2& p,
                                   const Point_2& q) const
    {
      Kernel     kernel;
      Segment_2  seg = kernel.construct_segment_2_object() (p, q);

      return (X_monotone_curve_2 (seg));
    }
  };

  /*! Obtain a Construct_x_monotone_curve_2 functor object. */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object () const
  {
    return Construct_x_monotone_curve_2();
  }
  //@}

};

/*!
 * \class A representation of a segment, as used by the Arr_segment_traits_2
 * traits-class.
 */
template <class Kernel_>
class Arr_hyperbola_object_2 
{

public:

  typedef Kernel_                                           Kernel;
  typedef typename Kernel::FT                               FT;

  typedef typename Kernel::Point_2                          Point_2;

private:
    Point_2   m_ps;               // The source point (if exists).
    Point_2   m_pt;               // The target point (if exists).
    bool      m_has_source;       // Is the source point valid
    bool      m_has_target;       // Is the target point valid
    bool      m_is_segment;       // Is segment.

   /* y = (a*x + b)/(cx + d) */
   FT m_a, m_b, m_c, m_d;
   
   
public:

  /*!
   * Default constructor.
   */
  Arr_hyperbola_object_2 ()
  {
     has_source = false;
     has_target = false;
  }
    
  /*!
   * Constructor from two points.
   * \param s The source point.
   * \param t The target point.
   * \pre The two points must not be the same.
   */
  Arr_hyperbola_object_2(const FT& a,
                         const FT& b,
                         const FT& c,
                         const FT& d,
                         const Point_2& ps,
                         const Point_2& pt)
   {
      m_ps = ps;
      m_pt = pt;
      m_has_source = true;
      m_has_target = true;
      m_a = a;
      m_b = b;
      m_c = c;
      m_d = d;

      m_is_segment = ((c == 0) && m_has_target && m_has_source);
   }


  /*!
   * Check whether the object is actually a segment.
   */
  bool is_segment () const
  {
     return m_is_segment;
  }


  // /*!
  //  * Check whether the object is actually a line.
  //  */
  // bool is_line () const
  // {
  //   return (! this->is_degen && ! this->has_source && ! this->has_target);
  // }

  // /*!
  //  * Cast to a line.
  //  * \pre The linear object is really a line.
  //  */
  // Line_2 line () const
  // {
  //   CGAL_precondition (is_line());
  //   return (this->l);
  // }


  /*!
   * Get the source point.
   * \pre The object is a point, a segment or a ray.
   */
  const Point_2& source() const
  {
    CGAL_precondition (has_source);

    return (this->ps);
  }

  /*!
   * Get the target point.
   * \pre The object is a point or a segment.
   */
  const Point_2& target() const
  {
    CGAL_precondition (has_target);

    return (this->pt);
  }

};

/*!
 * Exporter for the segment class used by the traits-class.
 */
template <class Kernel, class OutputStream>
OutputStream& operator<< (OutputStream& os,
                          const Arr_hyperbola_object_2<Kernel>& lobj)
{
  // Print a letter identifying the object type, then the object itself.
  if (lobj.is_segment())
    os << " S " << lobj.segment();
  else if (lobj.is_ray())
    os << " R " << lobj.ray();
  else
    os << " L " << lobj.line();

  return (os);
}

/*!
 * Importer for the segment class used by the traits-class.
 */
template <class Kernel, class InputStream>
InputStream& operator>> (InputStream& is, Arr_hyperbola_object_2<Kernel>& lobj)
{
  // Read the object type.
  char        c;

  do
  {
    is >> c;
  } while ((c != 'S' && c != 's') &&
           (c != 'R' && c != 'r') &&
           (c != 'L' && c != 'l'));

  // Read the object accordingly.
  if (c == 'S' || c == 's')
  {
    typename Kernel::Segment_2  seg;
    is >> seg;
    lobj = seg;
  }
  else if (c == 'R' || c == 'r')
  {
    typename Kernel::Ray_2      ray;
    is >> ray;
    lobj = ray;
  }
  else
  {
    typename Kernel::Line_2     line;
    is >> line;
    lobj = line;
  }

  return (is);
}

} //namespace CGAL

#endif
