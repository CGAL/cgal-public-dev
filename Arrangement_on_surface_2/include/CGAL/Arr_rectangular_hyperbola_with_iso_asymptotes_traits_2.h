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
#include <CGAL/Arr_geometry_traits/Rectangular_hyperbola_with_iso_asymptotes_2.h>

namespace CGAL {

/*! \class
 * A traits class for maintaining an arrangement of rectangular hyperbolas
 * with vertical and horiznotal asysmptotes.
 */
template <typename Kernel_>
class Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2 {
public:
  typedef Kernel_                                     Kernel;  
  typedef typename Kernel::FT                         NT;
 
  typedef X_monotone_rectangular_hyperbola_with_iso_asymptotes_2<Kernel>
                                                      X_monotone_curve_2;
  typedef Rectangular_hyperbola_with_iso_asymptotes_2<Kernel>
                                                      Curve_2;
  typedef typename Curve_2::Root_of_2_kernel          Root_of_2_kernel; 
  typedef typename Curve_2::Root_of_2_point           Point_2; 
  
  // Category tags:
  typedef Tag_true                                    Has_left_category;
  typedef Tag_true                                    Has_merge_category;
  typedef Tag_false                                   Has_do_intersect_category;

  typedef Arr_open_side_tag                           Left_side_category;
  typedef Arr_open_side_tag                           Bottom_side_category;
  typedef Arr_open_side_tag                           Top_side_category;
  typedef Arr_open_side_tag                           Right_side_category;
  
  typedef typename Kernel::Point_2                    Rational_point_2;
  typedef typename Kernel::Line_2                     Line_2;
  typedef typename Kernel::Ray_2                      Ray_2;
  typedef typename Kernel::Segment_2                  Segment_2;

private:
  typedef Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2<Kernel>
                                                      Self;
  
  mutable Kernel* m_kernel;
  bool            m_own_kernel;
  
public:
  /*! Default constructor.
   */
  Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2() :
    m_own_kernel(false)
  {
    m_kernel = new Kernel;
  }

  /*! Constructor from a kernel.
   */
  Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2(Kernel* kernel) :
    m_kernel(kernel),
    m_own_kernel(true)
  {}

  /*! Copy constructor.
   */
  Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2(const Self& other) :
    m_own_kernel(other.own_kernel())
  {
    m_kernel = (own_kernel) ?
      new Kernel(*other.kernel()) : other.kernel();
  }

  /*! Destructor
   */
  ~Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2()
  { if (m_own_kernel) delete(m_kernel); }
  
  /*! Obtains the kernel
   * \return the kernel.
   */
  Kernel* kernel() const {return m_kernel;}

  /*! Obtains the flag that indicates whether the kernel is owned by the traits.
   * \return the flag that indicates whether the kernel is owned by the traits.
   */
  bool own_kernel() const {return m_own_kernel;}
  
  /*! \class
   * A functor that checks whether a point is on a curve
   */
  class Is_on_2 {
    /*! Checks whether the point is on the vertical or horizontal line.
     * \param p The point.
     * \param is_vertical Indicates whether the line is vertical or horizontal.
     * \param d The d coefficient.
     * \return If vertical, true, if p.x = -d; false, otherwise.
     *         else, true, if p.y = -d; false, otherwise..
     */
    bool operator()(const Point_2& p, bool is_vertical, const NT& d) const
    { return (is_vertical) ? (p.x == -d) : (p.y == -d); }

    /*! Checks whether the point is on the hyperbola.
     * \param p The point.
     * \param b The b coefficien.
     * \param c The c coefficient.
     * \param d The d coefficient.
     * \return true, if p.x * p.y + b * p.x + c * p.y + d = 0; false, otherwise.
     */
    bool operator()(const Point_2& p,
                    const NT& b, const NT& c, const NT& d) const
    {
      assert(false);
      return true;
    }
    
    /*! Checks whether the point is on the curve.
     * \param p The point.
     * \param is_linear Indicates whether the curve is linear.
     * \param b The b coefficien.
     * \param c The c coefficient.
     * \param d The d coefficient.
     * \return If (is_linear),
     *           true, if b * p.x + c * p.y + d = 0; false, otherwise.
     *         else,
     *           true, if p.x * p.y + b * p.x + c * p.y + d = 0;
     *           false, otherwise.
     * \pre If is_linear, either (b == 0 && c != 0) or (b != 0 && c == 0).
     */
    bool operator()(const Point_2& p,
                    bool is_linear, const NT& b, const NT& c, const NT& d) const
    {
      if (is_linear) {
        CGAL_assertion(((b == 0) && (c != 0)) || ((b != 0) && (c == 0)));
        return (c == 0) ? operator()(true, d/b) : operator()(false, d/c)
      }
      return operator()(p, b, c, d);
    }
  };

  /*! Obtains a Is_on_2 functor. */
  Is_on_2 is_on_2_object() const { return Is_on_2(); }
  
  /// \name Basic functor definitions.
  //@{
  
  /*! \class
   * A functor that compares the x-coordinates of two points
   */
  typedef typename Root_of_2_kernel::Compare_x_2 Compare_x_2; 

  /*! Obtains a Compare_x_2 functor. */
  Compare_x_2 compare_x_2_object() const { return Compare_x_2(); }

  /*! \class
   * A functor that compares the xy-coordinates of two points
   */
  typedef typename Root_of_2_kernel::Compare_xy_2 Compare_xy_2; 
  

  /*! Obtains a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object() const { return Compare_xy_2(); }

  /*! \class
   * A functor that obtains the left endpoint of a segment or a ray.
   */
  class Construct_min_vertex_2 {
  public:
    /*! Obtains the left endpoint of the x-monotone curve.
     * \param xc The curve.
     * \pre The left end of xc is a valid (bounded) point.
     * \return The left endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2& xc) const
    {
      CGAL_precondition(xc.has_left());
      return xc.left();
    }
  };

  /*! Obtains a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(); }

  /*! \class
   * A functor that obtains the right endpoint of a segment or a ray.
   */
  class Construct_max_vertex_2 {
  public:
    /*! Obtains the right endpoint of the x-monotone curve.
     * \param xc The curve.
     * \pre The right end of xc is a valid (bounded) point.
     * \return The right endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2& xc) const
    {
      CGAL_precondition(xc.has_right());
      return xc.right();
    }
  };

  /*! Obtains a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(); }

  /*! \class
   * A functor that checks whether a given linear curve is vertical.
   */
  class Is_vertical_2 {
  public:
    /*! Checks whether the given x-monotone curve is a vertical segment.
     * \param xc The curve.
     * \return (true) if the curve is a vertical segment; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& xc) const
    { return xc.is_vertical(); }
  };

  /*! Obtains an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object() const { return Is_vertical_2(); }

  /*! A functor that compares the y-coordinates of a point and a line at
   * the point x-coordinate.
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
     * \param xc The curve.
     * \param p The point.
     * \pre p is in the x-range of xc.
     * \return SMALLER if y(p) < xc(x(p)), i.e. the point is below the curve;
     *         LARGER if y(p) > xc(x(p)), i.e. the point is above the curve;
     *         EQUAL if p lies on the curve.
     */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& xc) const
    {
      // TODO
      assert(false); 
      return EQUAL;
    }
  };

  /*! Obtains a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(this); }

  /*! \class
   * A functor that compares compares the y-coordinates of two linear
   * curves immediately to the left of their intersection point.
   */
  class Compare_y_at_x_left_2 {
  public:
    /*! Compares the y value of two x-monotone curves immediately to the left
     * of their intersection point.
     * \param xc1 The first curve.
     * \param xc2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its left.
     * \return The relative position of xc1 with respect to xc2 immdiately to
     *         the left of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator()(const X_monotone_curve_2& xc1,
                                 const X_monotone_curve_2& xc2,
                                 const Point_2& CGAL_precondition_code(p)) const
    {
      // TODO
      assert(false);
      return EQUAL;
    }
  };

  /*! Obtains a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(); }

  /*! \class
   * A functor that compares the y-coordinates of two linear
   * curves immediately to the right of their intersection point.
   */
  class Compare_y_at_x_right_2 {
  public:
    /*! Compares the y value of two x-monotone curves immediately to the right
     * of their intersection point.
     * \param xc1 The first curve.
     * \param xc2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its right.
     * \return The relative position of xc1 with respect to xc2 immdiately to
     *         the right of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator()(const X_monotone_curve_2& xc1,
                                 const X_monotone_curve_2& xc2,
                                 const Point_2& CGAL_precondition_code(p)) const
    {
      // TODO
      assert(false);
      return EQUAL;
    }
  };

  /*! Obtains a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(); }

  /*! \class
   * A functor that checks whether two points and two linear curves are
   * identical.
   */
  class Equal_2 {
  public:
    /*! Checks whether the two x-monotone curves are the same (have the same
     * graph).
     * \param xc1 The first curve.
     * \param xc2 The second curve.
     * \return (true) if the two curves are the same; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& xc1,
                    const X_monotone_curve_2& xc2) const
    { return (xc1 == xc2); }

    /*! Checks whether the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    bool operator()(const Point_2& p1, const Point_2& p2) const
    { return (p1 == p2); }
  };

  /*! Obtains an Equal_2 functor object. */
  Equal_2 equal_2_object() const { return Equal_2(); }
  //@}

  /// \name Functor definitions to handle boundaries
  //@{

  /*! \class
   * A function object that obtains the parameter space of a geometric
   * entity along the x-axis
   */
  class Parameter_space_in_x_2 {
  public:
    /*! Obtains the parameter space at the end of a line along the x-axis.
     * \param xc the line
     * \param ce the line end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the line xc.
     *   ARR_LEFT_BOUNDARY  - the curve end is unbounded at the left of the
     *                        parameter space.
     *   ARR_INTERIOR       - the curve end is bounded.
     *   ARR_RIGHT_BOUNDARY - the curve end is unbounded at the right of the
     *                        parameter space.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & xc,
                                   Arr_curve_end ce) const
    {
      return (ce == ARR_MIN_END) ?
        (xc.has_left_x() ? ARR_INTERIOR : ARR_LEFT_BOUNDARY) :
        (xc.has_right_x() ? ARR_INTERIOR : ARR_RIGHT_BOUNDARY);
    }
  };

  /*! Obtains a Parameter_space_in_x_2 function object */
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const
  { return Parameter_space_in_x_2(); }
  
  /*! \class
   * A function object that obtains the parameter space of a geometric
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
     * \param xc the line
     * \param ce the line end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the line xc.
     *   ARR_BOTTOM_BOUNDARY - the curve end is unbounded at the bottom of
     *                         the parameter space.     
     *   ARR_INTERIOR        - the curve end is bounded.
     *   ARR_TOP_BOUNDARY    - the curve end is unbounded at the top of
     *                         the parameter space.     
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & xc,
                                   Arr_curve_end ce) const
    {
      return (ce == ARR_MIN_END) ?
        (xc.has_left_y() ? ARR_INTERIOR :
         ((xc.b()*xc.c() < xc.d()) ? ARR_BOTTOM_BOUNDARY : ARR_TOP_BOUNDARY)) :
        (xc.has_right_y() ? ARR_INTERIOR :
         ((xc.b()*xc.c() < xc.d()) ? ARR_BOTTOM_BOUNDARY : ARR_TOP_BOUNDARY)) ;
    }
  };
  
  /*! Obtains a Parameter_space_in_y_2 function object */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(); }

  /*! \class
   * A function object that compares the x-limits of arc ends on the
   * boundary of the parameter space
   */
  class Compare_x_at_limit_2 {
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
    Compare_x_at_limit_2(const Traits* traits) : m_traits(traits) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2<Kernel>;

  public:
    /*! Compare the x-coordinate of a point with the x-coordinate of
     * the vertical asymptote of a hyperbola or a vertical line.
     * a line end on the boundary at y = +/- oo.
     * \param p the point direction.
     * \param xc the line, the endpoint of which is compared.
     * \param ce the line-end indicator -
     *            ARR_MIN_END - the minimal end of xc or
     *            ARR_MAX_END - the maximal end of xc.
     * \return the comparison result:
     *         SMALLER - x(p) < x(xc, ce);
     *         EQUAL   - x(p) = x(xc, ce);
     *         LARGER  - x(p) > x(xc, ce).     
     * \pre p lies in the interior of the parameter space.
     * \pre the ce end of the curve xc lies on the bottom or top boundary.
     */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2&  xc, 
                                 Arr_curve_end ce)
    {
      CGAL_precondition(m_traits->parameter_space_in_x_2_object()(xc, ce) ==
                        ARR_INTERIOR);
      CGAL_precondition(m_traits->parameter_space_in_y_2_object()(xc, ce) !=
                        ARR_INTERIOR);
      return CGAL::compare(p.x(),
                           (ce == ARR_MIN_END) ? xc.left_x() : xc.right_x());
    }

    /*! Compare the x-coordinates of the vertical asymptote of a hyperbolas or
     * vertical lines.
     * \param xc1 the first arc.
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
     * \pre the ce1 end of the curve xc1 lies on the bottom or top boundary.
     * \pre the ce2 end of the curve xc2 lies on the bottom or top boundary.
     */
    Comparison_result operator()(const X_monotone_curve_2&  xc1, 
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2&  xc2, 
                                 Arr_curve_end ce2)
    {
      CGAL_precondition_code
      (
        Parameter_space_in_x_2 psx = m_traits->parameter_space_in_x_2_object();
        Parameter_space_in_y_2 psy = m_traits->parameter_space_in_y_2_object();
      )
      CGAL_precondition(psx(xc1, ce1) == ARR_INTERIOR);
      CGAL_precondition(psy(xc1, ce1) != ARR_INTERIOR);
      CGAL_precondition(psx(xc2, ce2) == ARR_INTERIOR);
      CGAL_precondition(psy(xc2, ce2) != ARR_INTERIOR);

      return CGAL::compare((ce1 == ARR_MIN_END) ? xc1.left_x() : xc1.right_x(),
                           (ce2 == ARR_MIN_END) ? xc2.left_x() : xc2.right_x());
    }

  };

  /*! Obtains a Compare_x_at_limit_2 function object */
  Compare_x_at_limit_2 compare_x_at_limit_2_object() const
  { return Compare_x_at_limit_2(this); }

  /*! \class
   * A function object that compares the x-coordinates of arc ends near the
   * boundary of the parameter space
   */
  class Compare_x_near_limit_2 {
  public:
    /*! Compares the x-coordinates of 2 arcs ends near the boundary of the
     * parameter space at y = +/- oo.
     * \param xc1 the first arc.
     * \param xc2 the second arc.
     * \param ce the arc end indicator -
     *           ARR_MIN_END - the minimal end of xc2 or
     *           ARR_MAX_END - the maximal end of xc2.
     * \return the second comparison result:
     *         SMALLER - x(xc1, ce) < x(xc2, ce);
     *         EQUAL   - x(xc1, ce) = x(xc2, ce);
     *         LARGER  - x(xc1, ce) > x(xc2, ce).
     * \pre the ce end of the line xc1 lies on the bottom or top boundary.
     * \pre the ce end of the line xc2 lies on the bottom or top boundary.
     * \pre the $x$-coordinates of xc1 and xc2 at their limit at their ce
     *      ends are equal.
     */
    Comparison_result operator()(const X_monotone_curve_2 & xc1,
                                 const X_monotone_curve_2 & xc2,
                                 Arr_curve_end ce) const
    {
      CGAL_precondition_code
      (
        Parameter_space_in_x_2 psx = m_traits->parameter_space_in_x_2_object();
        Parameter_space_in_y_2 psy = m_traits->parameter_space_in_y_2_object();
      )
      CGAL_precondition(psx(xc1, ce) == ARR_INTERIOR);
      CGAL_precondition(psy(xc1, ce) != ARR_INTERIOR);
      CGAL_precondition(psx(xc2, ce) == ARR_INTERIOR);
      CGAL_precondition(psy(xc2, ce) != ARR_INTERIOR);
        
      if (xc1.is_vertical() && xc2.is_vertical()) return EQUAL;
      if (xc2.is_vertical()) return (ce == ARR_MAX_END) ? SMALLER : LARGER;
      if (xc1.is_vertical()) return (ce == ARR_MAX_END) ? LARGER : SMALLER;

      // Both curves are hyperbola.
      // TODO
      assert(false); 
      return SMALLER;
    }
  };

  /*! Obtains a Compare_x_near_limit_2 function object */
  Compare_x_near_limit_2 compare_x_near_limit_2_object() const
  { return Compare_x_near_limit_2(); }
    

  /*! \class
   * A function object that compares the y-coordinates of arc ends near the
   * boundary of the parameter space.
   */
  class Compare_y_near_boundary_2 {
  public:
    /*! Compares the y-coordinates of 2 lines at their ends near the boundary
     * of the parameter space at x = +/- oo.
     * \param xc1 the first arc.
     * \param xc2 the second arc.
     * \param ce the line end indicator.
     * \return the second comparison result:
     *         SMALLER - y(xc1, ce) < y(xc2, ce);
     *         EQUAL   - y(xc1, ce) = y(xc2, ce);
     *         LARGER  - y(xc1, ce) > y(xc2, ce).
     * \pre the ce ends of the lines xc1 and xc2 lie either on the left
     * boundary or on the right boundary of the parameter space.
     */
    Comparison_result operator()(const X_monotone_curve_2 & xc1,
                                 const X_monotone_curve_2 & xc2,
                                 Arr_curve_end ce) const
    {
      CGAL_precondition_code
      (
        Parameter_space_in_x_2 psx = m_traits->parameter_space_in_x_2_object();
        Parameter_space_in_y_2 psy = m_traits->parameter_space_in_y_2_object();
      )
      CGAL_precondition(psx(xc1, ce) != ARR_INTERIOR);
      CGAL_precondition(psy(xc1, ce) == ARR_INTERIOR);
      CGAL_precondition(psx(xc2, ce) != ARR_INTERIOR);
      CGAL_precondition(psy(xc2, ce) == ARR_INTERIOR);

      Comparison_result res =
        CGAL::compare((ce == ARR_MIN_END) ? xc1.left_x() : xc1.right_x(),
                      (ce == ARR_MIN_END) ? xc2.left_x() : xc2.right_x());
      if (res != EQUAL) return res;
      if (xc1.is_horizontal() && xc2.is_horizontal()) return EQUAL;

      // At most one of the curves is horizontal.
      // TODO
      assert(false); 
      return EQUAL;
    }
  };

  /*! Obtains a Compare_y_near_boundary_2 function object */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const
  { return Compare_y_near_boundary_2(); }
  
  //@}
  
  /// \name Functor definitions for supporting intersections.
  //@{

  /*! \class
   */
  class Make_x_monotone_2 {
  public:
    /*! Decomposes the given curve into x-monotone subcurves and insert them
     * into the given output iterator.
     * \param cv The curve.
     * \param oi The output iterator, whose value-type is Object. The output
     *           object is a wrapper of an X_monotone_curve_2 which is
     *           essentially the same as the input curve.
     * \return The past-the-end iterator.
     */
    template<typename OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const
    {
      if (cv.is_continuous()) return *oi++ = make_object(cv);

      const NT& b = cv.b();
      const NT& c = cv.c();
      const NT& d = cv.d();
      bool is_directed_right = cv.is_directed_right();
      Point_2 p(-c, 0);
      *oi++ = make_object(X_monotone_curve_2(true, b, c, d, cv.left(), p,
                                             cv.has_left_x(), cv.has_left_y(),
                                             true, false, is_directed_right));
      *oi++ = make_object(X_monotone_curve_2(true, b, c, d, p, cv.right(), 
                                             true, false,
                                             cv.has_right_x(), cv.has_right_y(),
                                             is_directed_right));

      return oi;
    }
  };

  /*! Obtains a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object () const
  { return Make_x_monotone_2(); }

  /*! \class
   */
  class Split_2 {
  public:
    /*! Splits a given x-monotone curve at a given point into two sub-curves.
     * \param xc The curve to split
     * \param p The split point.
     * \param xc1 Output: The left resulting subcurve (p is its right endpoint).
     * \param xc2 Output: The right resulting subcurve (p is its left endpoint).
     * \pre p lies on xc but is not one of its end-points.
     */
    void operator()(const X_monotone_curve_2& xc, const Point_2& p,
                    X_monotone_curve_2& xc1, X_monotone_curve_2& xc2) const
    {
      bool is_linear = xc.is_linear();
      const NT& b = xc.b();
      const NT& c = xc.c();
      const NT& d = xc.d();

      CGAL_precondition(m_traits->is_on_object_2()(p, is_linear, b, c, d));
      
      xc1 = X_monotone_curve_2(is_linear, b, c, d, xc.left(), p,
                               xc.has_left_x(), xc.has_left_y(), true, true,
                               xc.is_directed_right());
      xc2 = X_monotone_curve_2(is_linear, b, c, d, p, xc.right(), 
                               true, true, xc.has_right_x(), xc.has_right_y(), 
                               xc.is_directed_right());
    }
  };

  /*! Obtains a Split_2 functor object. */
  Split_2 split_2_object() const { return Split_2(); }

  /*! \class
   */
  class Intersect_2 {
  public:
    /*!
     * Find the intersections of the two given curves and insert them into the
     * given output iterator. As two segments may itersect only once, only a
     * single intersection will be contained in the iterator.
     * \param xc1 The first curve.
     * \param xc2 The second curve.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template<typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& xc1,
                              const X_monotone_curve_2& xc2,
                              OutputIterator oi) const
    {
      //TODO
      assert(false);
      return oi;
    }
  };

  /*! Obtains an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () const { return Intersect_2(); }

  /*! \class
   */
  class Are_mergeable_2 {
  public:
    /*!
     * Check whether it is possible to merge two given x-monotone curves.
     * \param xc1 The first curve.
     * \param xc2 The second curve.
     * \return (true) if the two curves are mergeable - if they are supported
     *         by the same line and share a common endpoint; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& xc1,
                    const X_monotone_curve_2& xc2) const
    {
      return ((xc1.is_linear() == xc2.is_linear()) &&
              (xc1.b() == xc2.b())
              (xc1.c() == xc2.c())
              (xc1.d() == xc2.d())
              (xc1.is_directed_right() == xc2.is_directed_right()) &&
              ((xc1.has_right_x() && xc1.has_right_y() &&
                xc2.has_left_x() && xc2.has_left_y() &&
                (xc1.right() == xc2.left())) ||
               (xc2.has_right_x() && xc2.has_right_y() &&
                xc1.has_left_x() && xc1.has_left_y() &&
                (xc2.right() == xc1.left()))));
    }
  };

  /*! Obtains an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object () const { return Are_mergeable_2(); }

  /*! \class
   */
  class Merge_2 {
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
    Merge_2(const Traits* traits) : m_traits(traits) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2<Kernel>;
    
  public:
    /*!
     * Merge two given x-monotone curves into a single curve (segment).
     * \param xc1 The first curve.
     * \param xc2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable, that is they are supported by the
     *      same line and share a common endpoint.
     */
    void operator()(const X_monotone_curve_2& xc1,
                    const X_monotone_curve_2& xc2,
                    X_monotone_curve_2& xc) const
    {
      CGAL_assertion(m_traits->are_mergeable_2_object()(xc1, xc2));

      if (xc1.has_right_x() && xc1.has_right_y() &&
          xc2.has_left_x() && xc2.has_left_y() &&
          (xc1.right() == xc2.left()))
      {
        xc = X_monotone_curve_2(xc1.is_linear(), xc1.b(), xc1.c(), xc1.d(),
                                xc1.left(), xc2.right(),
                                xc1.has_left_x(), xc1.has_left_y(),
                                xc2.has_right_x(), xc2.has_right_y(),
                                xc1.is_directed_right());
        return;
      }
      CGAL_assertion(xc2.has_right_x() && xc2.has_right_y() &&
                     xc1.has_left_x() && xc1.has_left_y() &&
                     (xc2.right() == xc1.left()));
      xc = X_monotone_curve_2(xc2.is_linear(), xc2.b(), xc2.c(), xc2.d(),
                              xc2.left(), xc1.right(),
                              xc2.has_left_x(), xc2.has_left_y(),
                              xc1.has_right_x(), xc1.has_right_y(),
                              xc2.is_directed_right());       
    }
  };

  /*! Obtains a Merge_2 functor object. */
  Merge_2 merge_2_object () const { return Merge_2(this); }
  //@}

  /// \name Functor definitions for the landmarks point-location strategy.
  //@{
  typedef double Approximate_number_type;

  /*! \class
   */
  class Approximate_2 {
  public:
    /*! Obtains an approximation of a point coordinate.
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

  /*! Obtains an Approximate_2 functor object. */
  Approximate_2 approximate_2_object () const { return Approximate_2(); }

  // No Construct_x_monotone_curve_2!
  //@}

  /// \name Functor definitions for the Boolean set-operation traits.
  //@{

  /*! \class
   */
  class Compare_endpoints_xy_2 {
  public:
    /*!
     * Compare the endpoints of an $x$-monotone curve lexicographically.
     * (assuming the curve has a designated source and target points).
     * \param xc The curve.
     * \return SMALLER if the curve is directed right;
     *         LARGER if the curve is directed left.
     */
    Comparison_result operator()(const X_monotone_curve_2& xc)
    { return (xc.is_directed_right()) ? SMALLER : LARGER; }
  };

  /*! Obtains a Compare_endpoints_xy_2 functor object. */
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(); }

  /*! \class
   */
  class Construct_opposite_2 {
  public:
    /*! Constructs an opposite x-monotone.
     * \param xc The curve.
     * \return The opposite curve.
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2& xc)
    {
      return X_monotone_curve_2(xc.is_linear(), xc.b(), xc.c(), xc.d(),
                                xc.right(), xc.left(), 
                                xc.has_right_x(), xc.has_right_y()
                                xc.has_left_x(), xc.has_left_y()
                                !xc.is_directed_right());
    }
  };

  /*! Obtains a Construct_opposite_2 functor object. */
  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(); }
  //@}

  /// \name Constructors.
  //@{

  /*! \class
   * A constructor object of curves.
   */
  template <typename CurveType_>
  class Construct_base_curve_2 {
  private:
    typedef CurveType_                                  Curve_type;    
    
  protected:
    typedef Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2<Kernel>
      Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Construct_base_curve_2(const Traits* traits) : m_traits(traits) {}
    
  public:
    /*! Constructor of either a vertical or a horizontal line.
     * \param is_vertical Indicates whether the curve is vertical or horizontal.
     * \param d The d coefficient.
     * \pre if the curve is a line, it must be either vertical or horizontal.
     */
    Curve_type operator()(bool is_vertical, const NT& d)
    {
      Point_2 p(-d, -d);
      Curve_type xc = (is_vertical) ?
        Curve_type(true, 1, 0, d, p, p, true, false, true, false, true) :
        Curve_type(true, 0, 1, d, p, p, false, true, false, true, true);
      return xc;
    }
    
    /*! Constructor of a vertical or a horizontal segment.
     * \param is_vertical Indicates whether the curve is vertical or horizontal.
     * \param d The d coefficient.
     * \param source The source point.
     * \param target The target point.
     * \pre The two points must not be the same.
     * \pre source is on the underlying line.
     * \pre target is on the underlying line.
     */
    Curve_type operator()(bool is_vertical, const NT& d,
                          const Point_2& source, const Point_2& target)
    {
      Comparison_result res = CGAL::compare_xy(source, target);
      CGAL_assertion(res != EQUAL);
      bool is_directed_right = (res == SMALLER);

      bool is_vertical = (c == 0);
      CGAL_assertion(m_traits->is_on_object_2()(left, is_vertical, d));
      CGAL_assertion(m_traits->is_on_object_2()(right, is_vertical, d));
        
      Curve_type xc = (is_vertical) ?
        ((is_directed_right) ?
         Curve_type(true, 1, 0, d, source, target, true, true, true, true,
                    is_directed_right) :
         Curve_type(true, 1, 0, d, target, source, true, true, true, true,
                    is_directed_right)) :
        ((is_directed_right) ?
         Curve_type(true, 0, 1, d, source, target, true, true, true, true,
                    is_directed_right) :
         Curve_type(true, 0, 1, d, target, source, true, true, true, true,
                    is_directed_right));
      return xc;
    }

    /*! Constructor of a curve bounded at source and target.
     * \param is_linear Indicates whether the curve is linear.
     * \param b The b coefficient.
     * \param c The c coefficient.
     * \param d The d coefficient.
     * \param source The source point.
     * \param target The target point.
     * \pre The two points must not be the same.
     * \pre source is on the underlying hyperbola or the curve has a
     *      vertical asymptote at the x-coordinate of source.
     * \pre target is on the underlying hyperbola or the curve has a
     *      vertical asymptote at the x-coordinate of right.
     * \pre The curve is continueous. That is, the open interval bounded by
     *      the x-coordinates of source and target does not contain -c (the
     *      x-coordinate of the vertical asymptotes).
     */
    Curve_type operator()(bool is_linear, const NT& b, const NT& c, const NT& d,
                          const Point_2& source, const Point_2& target)
    {
      if (is_linear) {
        CGAL_assertion(((b == 0) && (c != 0)) || ((b != 0) && (c == 0)));
        return (c == 0) ?
          operator()(true, d/b, source, target) :
          operator()(false, d/c, source, target) :
      }
      return operator()(b, c, d, source, target);
    }

    /*! Constructor of a curve bounded at source and target.
     * \param a The a coefficient.
     * \param b The b coefficient.
     * \param c The c coefficient.
     * \param d The d coefficient.
     * \param source The source point.
     * \param target The target point.
     * \pre The two points must not be the same.
     * \pre source is on the underlying hyperbola or the curve has a
     *      vertical asymptote at the x-coordinate of source.
     * \pre target is on the underlying hyperbola or the curve has a
     *      vertical asymptote at the x-coordinate of right.
     * \pre The curve is continueous. That is, the open interval bounded by
     *      the x-coordinates of source and target does not contain -c (the
     *      x-coordinate of the vertical asymptotes).
     */
    Curve_type operator()(bool is_linear,
                          const NT& a, const NT& b, const NT& c, const NT& d,
                          const Point_2& source, const Point_2& target)
    { 
      if (a = 0) {
        CGAL_assertion(((b == 0) && (c != 0)) || ((b != 0) && (c == 0)));
        return (c == 0) ?
          operator()(true, d/b, source, target) :
          operator()(false, d/c, source, target) :
      }
      return operator()(b/a, c/a, d/a, source, target);
    }
    
    /*! Constructor of a vertical or horizontal ray.
     * \param is_vertical Indicates whether the ray is vertical or horizontal.
     * \param d The d coefficient.
     * \param source The source point.
     * \param is_directed_right Indicates whether the curve is directed right
     *        (top is right).
     * \pre source is on the underlying line.
     */
    Curve_type operator()(bool is_vertical, const NT& d,
                          const Point_2& source, bool is_directed_right)
    {
      Curve_type xc = (is_vertical) ?
        Curve_type(true, 1, 0, d, source, source,
                   is_directed_right, is_directed_right,
                   !is_directed_right, !is_directed_right, is_directed_right) :
        Curve_type(true, 0, 1, d, source, source,
                   is_directed_right, is_directed_right,
                   !is_directed_right, !is_directed_right, is_directed_right);        
      return xc;
    }

    /*! Constructor of a hyperbola bounded at one endpoint.
     * (a) If is_directed_right
     *     (i)  If source.x < -c, then has_left_x <- true
     *     (ii) Otherwise (source.x > -c), has_left_x <- false
     * (b) Otherwise (!is_directed_right)
     *     (i)  If source.x > -c, then has_right_x <- true
     *     (ii) Otherwise (source.x > -c), has_right_x <- false
     * \param is_linear Indicates whether the curve is linear.
     * \param b The b coefficient.
     * \param c The c coefficient.
     * \param d The d coefficient.
     * \param source The source point.
     * \param is_directed_right Indicates whether the curve is directed right.
     * \pre source is on the underlying hyperbola or the curve has a
     *      vertical asymptote at the x-coordinate of source.
     */
    Curve_type operator()(bool is_linear, const NT& b, const NT& c, const NT& d,
                          const Point_2& source, bool is_directed_right)
    {
      if (is_linear) {
        CGAL_assertion(((b == 0) && (c != 0)) || ((b != 0) && (c == 0)));
        return (c == 0) ?
          operator()(true, d/b, source, is_directed_right) :
          operator()(false, d/c, source, is_directed_right) :
      }
      return operator()(b, c, d, source, is_directed_right);
    }

    /*! Constructor of a hyperbola bounded at one endpoint.
     * (a) If is_directed_right
     *     (i)  If source.x < -c, then has_left_x <- true
     *     (ii) Otherwise (source.x > -c), has_left_x <- false
     * (b) Otherwise (!is_directed_right)
     *     (i)  If source.x > -c, then has_right_x <- true
     *     (ii) Otherwise (source.x > -c), has_right_x <- false
     * \param a The a coefficient.
     * \param b The b coefficient.
     * \param c The c coefficient.
     * \param d The d coefficient.
     * \param source The source point.
     * \param is_directed_right Indicates whether the curve is directed right.
     * \pre source is on the underlying hyperbola or the curve has a
     *      vertical asymptote at the x-coordinate of source.
     */
    Curve_type operator()(const NT& a, const NT& b, const NT& c, const NT& d,
                          const Point_2& source, bool is_directed_right)
    { 
      if (a = 0) {
        CGAL_assertion(((b == 0) && (c != 0)) || ((b != 0) && (c == 0)));
        return (c == 0) ?
          operator()(true, d/b, source, is_directed_right) :
          operator()(false, d/c, source, is_directed_right) :
      }
      return operator()(b/a, c/a, d/a, source, is_directed_right);
    }
    
    /*! Constructor of a curve from a line
     * \param line The line.
     * \pre the line must be either vertical or horizontal
     */
    Curve_type operator()(const Line_2& line)
    {
      CGAL_assertion(is_vertical(line) || is_horizontal(line));

      Kernel* kernel = m_traits->kernel();
      Is_vertical_2 is_vertical = kernel->is_vertical_2_object();
      CGAL_assertion(is_vertical(line) ||
                     kernel->is_horizontal_2_object()(line));
      return (is_vertical) ?
        operator()(true, line.c()/line.a()) :
        operator()(false, line.c()/line.b());
    }

    /*! Constructor of a curve from a ray
     * \param ray The ray.
     * \pre the ray must be either vertical or horizontal
     */
    Curve_type operator()(const Ray_2& ray)
    {
      Kernel* kernel = m_traits->kernel();
      typename Kernel::Construct_point_on_2
        construct_vertex = kernel.construct_point_on_2_object();
      Rational_point_2 source = construct_vertex(ray, 0); // The source point.
      Rational_point_2 target = construct_vertex(ray, 1); // Some point on ray.
      Comparison_result res = kernel.compare_xy_2_object()(source, target);
      CGAL_assertion(res != EQUAL);
      bool is_directed_right = (res == SMALLER);
      Point_2 ps(source);
      Point_2 pt(target);
      Curve_type xc;
      xc = (is_directed_right) ?
        Curve_type(true, b, c, d, ps, pt, true, true, false, false,
                   is_directed_right) :
        Curve_type(true, b, c, d, pt, ps, false, false, true, true, 
                   is_directed_right);      
      return xc;
    }

    /*! Constructor of a curve from a segment
     * \param segment The segment.
     * \pre the segment must be either vertical or horizontal
     */
    Curve_type operator()(const Segment_2& segment)
    {
      Kernel* kernel = m_traits->kernel();
      CGAL_assertion(is_vertical(segment) ||
                     kernel->is_horizontal_2_object()(segment));
      
      typename Kernel::Construct_point_on_2
        construct_vertex = kernel.construct_point_on_2_object();
      Rational_point_2 source = construct_vertex(segment, 0); // source point.
      Rational_point_2 target = construct_vertex(segment, 1); // target point.
      Comparison_result res = kernel.compare_xy_2_object()(source, target);
      CGAL_assertion(res != EQUAL);
      bool is_directed_right = (res == SMALLER);
      Point_2 ps(source);
      Point_2 pt(target);
      return (is_directed_right) ?
        ((is_vertical) ?
         operator()(true, line.c()/line.a(), ps, pt) :
         operator()(false, line.c()/line.b(), ps, pt)) :
        ((is_vertical) ?
         operator()(true, line.c()/line.a(), pt, ps) :
         operator()(false, line.c()/line.b(), pt, ps));
    }
  };

  /*! \class
   * A constructor object of x-monotone curves.
   */
  class Construct_x_monotone_curve_2 :
    public Construct_base_curve_2<X_monotone_curve_2>
  {
  private:
    typedef Construct_base_curve_2<X_monotone_curve_2>  Base;

  protected:
    /*! Constructor
     * \param traits the traits (in case it has state)
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Construct_x_monotone_curve_2(const Traits* traits) : Base(traits) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2<Kernel>;

  public:
    /*! Constructor from all data fields.
     * \param is_linear Indicates whether the curve is linear.
     * \param b The b coefficient.
     * \param c The c coefficient.
     * \param d The d coefficient.
     * \param left The left point.
     * \param right The right point.
     * \param is_directed_right Indicates whether the curve is directed right.
     * \param has_left_x Indicates whether the left endpoint of the curve has a
     *        valid x-coordinate stored as the x-coordinate of m_left.
     * \param has_left_y Indicates whether the left endpoint of the curve has a
     *        valid y-coordinate stored as the y-coordinate of m_left.
     * \param has_right_x Indicates whether the right endpoint of the curve has
     *        a valid x-coordinate stored as the x-coordinate of m_right.
     * \param has_right_y Indicates whether the right endpoint of the curve has
     *        a valid y-coordinate stored as the y-coordinate of m_right.
     * \param is_directed_right Indicates whether the curve is directed right.
     * \param is_continuous Indicates the curve continuous.
     * \pre The two points must not be the same.
     * \pre If has_left_x && has_left_y, m_left is on the underlying hyperbola.
     *      If has_left_x && !has_left_y, the left end of the underlying
     *      hyperbola has a vertical asymptote at the x-coordinate of m_left.
     *      If !has_left_x && has_left_y, the left end of the underlying
     *      hyperbola has a horizontal asymptote at the y-coordinate of m_left.
     * \pre If has_right_x && has_right_y, m_right is on the underlying
     *      hyperbola.
     *      If has_right_x && !has_right_y, the right end of the underlying
     *      hyperbola has a vertical asymptote at the x-coordinate of m_right.
     *      If !has_right_x && has_right_y, the right end of the underlying
     *      hyperbola has a horizontal asymptote at the y-coordinate of m_right.
     * \pre The curve is continueous.
     */
    X_monotone_curve_2 operator()(bool is_linear,
                                  const NT& b, const NT& c, const NT& d,
                                  const Point_2& left, const Point_2& right,
                                  bool has_left_x, bool has_left_y,
                                  bool has_right_x, bool has_right_y,
                                  bool is_directed_right)
    {    
      // Validity check:
      CGAL_assertion(!has_left_x || !has_left_y || 
                     m_traits->is_on_object_2()(left, is_linear, b, c, d));
      CGAL_assertion(!has_left_x || has_left_y || (left.x() == -c));
      CGAL_assertion(has_left_x || !has_left_y || (left.y() == -b));

      CGAL_assertion(!has_right_x || !has_right_y || 
                     m_traits->is_on_object_2()(right, is_linear, b, c, d));
      CGAL_assertion(!has_right_x || has_right_y || (right.x() == -c));
      CGAL_assertion(has_right_x || !has_right_y || (right.y() == -b));

      // Continuity check:
      CGAL_assertion(!has_left_x  || (-c <= left.x()));
      CGAL_assertion(!has_right_x || (right.x() <= -c));
            
      X_monotone_curve_2 xc =
        X_monotone_curve_2(is_linear, b, c, d, left, right,
                           has_left_x, has_left_y, has_right_x, has_right_y,
                           is_directed_right);
      return xc;
    }
    
    /*! Constructor of a hyperbolic arc.
     * \param b The b coefficient.
     * \param c The c coefficient.
     * \param d The d coefficient.
     * \param source The source point.
     * \param target The target point.
     * \pre The two points must not be the same.
     * \pre source is on the underlying hyperbola or the curve has a
     *      vertical asymptote at the x-coordinate of source.
     * \pre target is on the underlying hyperbola or the curve has a
     *      vertical asymptote at the x-coordinate of right.
     * \pre The curve is continueous. That is, the open interval bounded by
     *      the x-coordinates of source and target does not contain -c (the
     *      x-coordinate of the vertical asymptotes).
     */
    X_monotone_curve_2 operator()(const NT& b, const NT& c, const NT& d,
                                  const Point_2& source, const Point_2& target)
    {
      Comparison_result res = CGAL::compare_xy(source, target);
      CGAL_assertion(res != EQUAL);
      bool is_directed_right = (res == SMALLER);
      
      NT asymptote_x = -c;
      const Point_2& left = (is_directed_right) ? source : target;
      const Point_2& right = (is_directed_right) ? target : left;
      bool has_left_y = asymptote_x != left.x();
      bool has_right_y = asymptote_x != right.x();
      CGAL_assertion((asymptote_x <= left.x()) || (right.x() <= asymptote_x));

      CGAL_assertion(!has_left_y || 
                     m_traits->is_on_object_2()(left, is_linear, b, c, d));
      CGAL_assertion(!has_right_y || 
                     m_traits->is_on_object_2()(right, is_linear, b, c, d));
        
      X_monotone_curve_2 xc =
        X_monotone_curve_2(is_linear, b, c, d, left, right,
                           true, has_left_y, true, has_right_y,
                           is_directed_right);
      return xc;
    }

    /*! Constructor of a hyperbola bounded at one endpoint.
     * (a) If is_directed_right
     *     (i)  If source.x < -c, then has_left_x <- true
     *     (ii) Otherwise (source.x > -c), has_left_x <- false
     * (b) Otherwise (!is_directed_right)
     *     (i)  If source.x > -c, then has_right_x <- true
     *     (ii) Otherwise (source.x > -c), has_right_x <- false
     * \param b The b coefficient.
     * \param c The c coefficient.
     * \param d The d coefficient.
     * \param source The source point.
     * \param is_directed_right Indicates whether the curve is directed right.
     * \pre source is on the underlying hyperbola or the curve has a
     *      vertical asymptote at the x-coordinate of source.
     */
    X_monotone_curve_2 operator()(const NT& b, const NT& c, const NT& d,
                                  const Point_2& source, bool is_directed_right)
    {
      NT asymptote_x = -c;
      NT asymptote_y = -b;
      Point_2 target = Point_2(asymptote_x, asymptote_y);
      bool has_left_x, has_left_y, has_right_x, has_right_y;
      if (is_directed_right) {
        has_left_x = true;
        if (asymptote_x < source.x()) {
          has_left_y = true;
          has_right_x = false;
          has_right_y = true;
        } else if (source.x() == asymptote_x) {
          has_left_y = false;
          has_right_x = false;
          has_right_y = true;
        } else {
          has_left_y = true;
          has_right_x = true;
          has_right_y = false;
        }
      } else {
        has_right_x = true;
        if (asymptote_x > source.x()) {
          has_right_y = true;
          has_left_x = false;
          has_left_y = true;
        } else if (source.x() == asymptote_x) {
          has_right_y = false;
          has_left_x = false;
          has_left_y = true;
        } else {
          has_right_y = true;
          has_left_x = true;
          has_left_y = false;
        }
      }
      
      const Point_2& left = (is_directed_right) ? source : target;
      const Point_2& right = (is_directed_right) ? target : left;
      
      X_monotone_curve_2 xc =
        X_monotone_curve_2(false, b, c, d, left, right,
                           has_left_x, has_left_y, has_right_x, has_right_y,
                           is_directed_right);
      return xc;
    }
  };
  
  /*! Obtains a Construct_x_monotone_curve_2 functor object. */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const
  { return Construct_x_monotone_curve_2(this); }

  /*! \class
   * A constructor object of general curves.
   */
  class Construct_curve_2 :
    public Construct_base_curve_2<Curve_2>
  {
  private:
    typedef Construct_base_curve_2<Curve_2>     Base;

  protected:
    /*! Constructor
     * \param traits the traits (in case it has state)
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Construct_curve_2(const Traits* traits) : Base(traits) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2<Kernel>;

  public:
    /*! Constructor from all data fields.
     * \param is_linear Indicates whether the curve is linear.
     * \param b The b coefficient.
     * \param c The c coefficient.
     * \param d The d coefficient.
     * \param left The left point.
     * \param right The right point.
     * \param is_directed_right Indicates whether the curve is directed right.
     * \param has_left_x Indicates whether the left endpoint of the curve has a
     *        valid x-coordinate stored as the x-coordinate of m_left.
     * \param has_left_y Indicates whether the left endpoint of the curve has a
     *        valid y-coordinate stored as the y-coordinate of m_left.
     * \param has_right_x Indicates whether the right endpoint of the curve has
     *        a valid x-coordinate stored as the x-coordinate of m_right.
     * \param has_right_y Indicates whether the right endpoint of the curve has
     *        a valid y-coordinate stored as the y-coordinate of m_right.
     * \param is_directed_right Indicates whether the curve is directed right.
     * \param is_continuous Indicates the curve continuous.
     * \pre The two points must not be the same.
     * \pre If has_left_x && has_left_y, m_left is on the underlying hyperbola.
     *      If has_left_x && !has_left_y, the left end of the underlying
     *      hyperbola has a vertical asymptote at the x-coordinate of m_left.
     *      If !has_left_x && has_left_y, the left end of the underlying
     *      hyperbola has a horizontal asymptote at the y-coordinate of m_left.
     * \pre If has_right_x && has_right_y, m_right is on the underlying
     *      hyperbola.
     *      If has_right_x && !has_right_y, the right end of the underlying
     *      hyperbola has a vertical asymptote at the x-coordinate of m_right.
     *      If !has_right_x && has_right_y, the right end of the underlying
     *      hyperbola has a horizontal asymptote at the y-coordinate of m_right.
     */
    Curve_2 operator()(bool is_linear, const NT& b, const NT& c, const NT& d,
                       const Point_2& left, const Point_2& right,
                       bool has_left_x, bool has_left_y,
                       bool has_right_x, bool has_right_y,
                       bool is_directed_right)
    {    
      // Validity check:
      CGAL_assertion(!has_left_x || !has_left_y || 
                     m_traits->is_on_object_2()(left, is_linear, b, c, d));
      CGAL_assertion(!has_left_x || has_left_y || (left.x() == -c));
      CGAL_assertion(has_left_x || !has_left_y || (left.y() == -b));

      CGAL_assertion(!has_right_x || !has_right_y || 
                     m_traits->is_on_object_2()(right, is_linear, b, c, d));
      CGAL_assertion(!has_right_x || has_right_y || (right.x() == -c));
      CGAL_assertion(has_right_x || !has_right_y || (right.y() == -b));

      // Continuity check:
      bool is_continuous = ((!has_left_x  || (-c <= left.x())) && 
                            (!has_right_x || (right.x() <= -c)));
            
      Curve_2 xc =
        Curve_2(is_linear, b, c, d, left, right,
                has_left_x, has_left_y, has_right_x, has_right_y,
                is_directed_right, is_continuous);
      return xc;
    }

    /*! Constructor of a hyperbola.
     * \param b The b coefficient.
     * \param c The c coefficient.
     * \param d The d coefficient.
     */
    Curve_2 operator()(const NT& b, const NT& c, const NT& d)
    {
      Point_2 p(-c, -b);
      Curve_2 xc =
        Curve_2(false, b, c, d, p, p, false, true, false, true, true, false);
      return xc;
    }
    
    /*! Constructor of a hyperbolic arc.
     * \param b The b coefficient.
     * \param c The c coefficient.
     * \param d The d coefficient.
     * \param source The source point.
     * \param target The target point.
     * \pre The two points must not be the same.
     * \pre source is on the underlying hyperbola or the curve has a
     *      vertical asymptote at the x-coordinate of source.
     * \pre target is on the underlying hyperbola or the curve has a
     *      vertical asymptote at the x-coordinate of right.
     */
    Curve_2 operator()(const NT& b, const NT& c, const NT& d,
                       const Point_2& source, const Point_2& target)
    {
      Comparison_result res = CGAL::compare_xy(source, target);
      CGAL_assertion(res != EQUAL);
      bool is_directed_right = (res == SMALLER);
      
      NT asymptote_x = -c;
      const Point_2& left = (is_directed_right) ? source : target;
      const Point_2& right = (is_directed_right) ? target : left;
      bool has_left_y = asymptote_x != left.x();
      bool has_right_y = asymptote_x != right.x();

      CGAL_assertion(!has_left_y || 
                     m_traits->is_on_object_2()(left, is_linear, b, c, d));
      CGAL_assertion(!has_right_y || 
                     m_traits->is_on_object_2()(right, is_linear, b, c, d));
      
      bool is_continuous =
        (asymptote_x <= left.x()) || (right.x() <= asymptote_x);
      Curve_2 xc =
        Curve_2(is_linear, b, c, d, left, right,
                true, has_left_y, true, has_right_y,
                is_directed_right, is_continuous);
      return xc;
    }
    
    /*! Constructor of a hyperbola bounded at one endpoint.
     * \param b The b coefficient.
     * \param c The c coefficient.
     * \param d The d coefficient.
     * \param sourse The left point.
     * \param is_directed_right Indicates whether the curve is directed right.
     * \pre source is on the underlying hyperbola or the curve has a
     *      vertical asymptote at the x-coordinate of source.
     */
    Curve_2 operator()(const NT& b, const NT& c, const NT& d,
                       const Point_2& source, bool is_directed_right)
    {
      NT asymptote_x = -c;
      NT asymptote_y = -b;
      Point_2 target = Point_2(asymptote_x, asymptote_y);
      bool has_left_x, has_left_y, has_right_x, has_right_y;
      bool is_continuous;
      if (is_directed_right) {
        has_left_x = true;
        has_right_x = false;
        has_right_y = true;
        if (asymptote_x < source.x()) {
          has_left_y = true;
          is_continuous = true;
        } else if (source.x() == asymptote_x) {
          has_left_y = false;
          is_continuous = true;
        } else {
          has_left_y = true;
          is_continuous = false;
        }
      } else {
        has_right_x = true;
        has_left_x = false;
        has_left_y = false;
        if (asymptote_x > source.x()) {
          has_right_y = true;
          is_continuous = true;
        } else if (source.x() == asymptote_x) {
          has_right_y = false;
          is_continuous = true;
        } else {
          has_right_y = true;
          is_continuous = false;
        }
      }

      const Point_2& left = (is_directed_right) ? source : target;
      const Point_2& right = (is_directed_right) ? target : left;
      
      Curve_2 xc =
        Curve_2(false, b, c, d, left, right,
                has_left_x, has_left_y, has_right_x, has_right_y,
                is_directed_right, is_continuous);
      return xc;
    }
  };
  
  /*! Obtains a Construct_curve_2 functor object. */
  Construct_curve_2 construct_curve_2_object() const
  { return Construct_curve_2(this); }
  //@}
};

} //namespace CGAL

#endif
