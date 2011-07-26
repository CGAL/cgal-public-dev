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
// $URL$
// $Id$
// 
//
// Author(s): Asaf Porat          <asafpor@post.tau.ac.il>
//            Efi Fogel           <efif@post.tau.ac.il>

#ifndef CGAL_RECTANGULAR_HYPERBOLA_WITH_ISO_ASYMPTOTES_2_H
#define CGAL_RECTANGULAR_HYPERBOLA_WITH_ISO_ASYMPTOTES_2_H

/*! \file
 * A header file for the class templates that represents a general and an
 * arbitrary rectangular hyperbola with vertical and horizontal asymptotes.
 */

#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Arr_geometry_traits/Sqrt_extension_point_2.h>

namespace CGAL {

// Forward declaration:
template <typename NumberType_, bool Filter_>
class Rectangular_hyperbola_with_iso_asymptotes_2;

/*! \class
 * A representation of a rectangular hyperbola with vertical and horizontal
 * asymptotes.
 */
template <typename Kernel_, bool Filter_>
class Rectangular_hyperbola_with_iso_asymptotes_2_rep {
public:
  typedef Kernel_                                      Kernel;
  typedef Filter_                                      Filter;

  typedef typename Kernel::NT                          NT;
  typedef Sqrt_extension_point_2<NT, Filter>           Point_2;

  friend class Rectangular_hyperbola_with_iso_asymptotes_2<Kernel, Filter>
  
private:
  typedef Rectangular_hyperbola_with_iso_asymptotes_2_rep<Kernel, Filter>
                                                        Self;
  
  bool a;                        // Indicates whether the curve is linear.
  NT m_b, m_c, m_d;              // axy + bx + cy + d = 0
  
  Point_2   m_left;              // The left point (if exists).
  Point_2   m_right;             // The right point (if exists).
  bool      m_has_left_x;        // The left endpoint has a valid x-coordinate. 
  bool      m_has_left_y;        // The left endpoint has a valid y-coordinate. 
  bool      m_has_right_x;       // The right endpoint has a valid x-coordinate. 
  bool      m_has_right_y;       // The right endpoint has a valid y-coordinate. 
  bool      m_is_directed_right; // Is the curve directed right?
  bool      m_is_continuous;     // Is the curve continuous?

public:
  /*! Default constructor.
   */
  Rectangular_hyperbola_with_iso_asymptotes_2_rep() :
    has_left_x(false),
    has_left_y(false),
    has_right_x(false),
    has_right_y(false)
  {}

  /*! Constructor from all data fields.
   * \param a The a coefficient (either 0 or 1).
   * \param b The a coefficient.
   * \param c The a coefficient.
   * \param d The a coefficient.
   * \param left The left point.
   * \param right The right point.
   * \param is_directed_right Indicates whether the curve is directed right.
   * \param has_left_x Indicates whether the left endpoint of the curve has a
   *        valid x-coordinate stored as the x-coordinate of m_left.
   * \param has_left_y Indicates whether the left endpoint of the curve has a
   *        valid y-coordinate stored as the y-coordinate of m_left.
   * \param has_right_x Indicates whether the right endpoint of the curve has a
   *        valid x-coordinate stored as the x-coordinate of m_right.
   * \param has_right_y Indicates whether the right endpoint of the curve has a
   *        valid y-coordinate stored as the y-coordinate of m_right.
   * \param is_directed_right Indicates whether the curve is directed right.
   * \param is_continuous Indicates the curve continuous.
   * \pre The two points must not be the same.
   * \pre If has_left_x && has_left_y, m_left is on the underlying hyperbola.
   *      If has_left_x && !has_left_y, the left end of the underlying
   *      hyperbola has a vertical asymptote at the x-coordinate of m_left.
   *      If !has_left_x && has_left_y, the left end of the underlying
   *      hyperbola has a horizontal asymptote at the y-coordinate of m_left.
   * \pre If has_right_x && has_right_y, m_right is on the underlying hyperbola.
   *      If has_right_x && !has_right_y, the right end of the underlying
   *      hyperbola has a vertical asymptote at the x-coordinate of m_right.
   *      If !has_right_x && has_right_y, the right end of the underlying
   *      hyperbola has a horizontal asymptote at the y-coordinate of m_right.
   */
  Rectangular_hyperbola_with_iso_asymptotes_2(bool a,
                                              const NT& b,
                                              const NT& c,
                                              const NT& d,
                                              const Point_2& left,
                                              const Point_2& rihgt,
                                              bool has_left_x, bool has_left_y,
                                              bool has_right_x, bool has_rigth_y,
                                              bool is_directed_right,
                                              bool is_continuous) :
    m_a(a),
    m_b(b),
    m_c(c),
    m_d(d),
    m_left(left),
    m_right(right),
    m_has_left_x(has_left_x),
    m_has_left_y(has_left_y),
    m_has_right_x(has_right_x),
    m_has_right_y(has_right_y),
    m_is_directed_right(is_directed_right),
    m_is_continuous(is_continuous)
  {}
};

/*! \class
 * A representation of a rectangular hyperbola with vertical and horizontal
 * asymptotes used by the traits class
 * Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2.
 */
template <typename Kernel_, bool Filter_>
class Rectangular_hyperbola_with_iso_asymptotes_2 :
  public Handle_for<Rectangular_hyperbola_with_iso_asymptotes_2_rep<Kernel_,
                                                                    Filter_> >
{
public:
  typedef Kernel_                                       Kernel;
  typedef Filter_                                       Filter;

private:  
  typedef Rectangular_hyperbola_with_iso_asymptotes_2<Kernel, Filter>
                                                        Self;
  typedef Rectangular_hyperbola_with_iso_asymptotes_2_rep<NT, Filter>
                                                        Curve_rep;
  typedef Handle_for<Curve_rep>
                                                        Curve_handle;

public:
  typedef typename Kernel::NT                           NT;
  typedef Sqrt_extension_point_2<NT, Filter>            Point_2;
  typedef typename Curve_rep::Coord_NT                  Coord_NT;

  /*! Default constructor.
   */
  Rectangular_hyperbola_with_iso_asymptotes_2() : Curve_handle(Curve_rep()) {}

  /*! Copy constructor.
   */
  Rectangular_hyperbola_with_iso_asymptotes_2(const Self& cv) :
    Curve_handle(cv) {}
  
  /*! Constructor from all data fields.
   * \param a The a coefficient (either 0 or 1).
   * \param b The a coefficient.
   * \param c The a coefficient.
   * \param d The a coefficient.
   * \param left The left point.
   * \param right The right point.
   * \param has_left_x Indicates whether the left endpoint of the curve has a
   *        valid x-coordinate.
   * \param has_left_y Indicates whether the left endpoint of the curve has a
   *        valid y-coordinate.
   * \param has_right_x Indicates whether the right endpoint of the curve has a
   *        valid x-coordinate.
   * \param has_right_y Indicates whether the right endpoint of the curve has a
   *        valid y-coordinate.
   * \param is_continuous Indicates the curve continuous.
   * \pre The two points must not be the same.
   * \pre The two points must not be the same.
   * \pre If has_left_x && has_left_y, m_left is on the underlying hyperbola.
   *      If has_left_x && !has_left_y, the left end of the underlying
   *      hyperbola has a vertical asymptote at the x-coordinate 'm_left'.
   *      If !has_left_x && has_left_y, the left end of the underlying
   *      hyperbola has a horizontal asymptote at the y-coordinate 'm_left'.
   * \pre If has_right_x && has_right_y, m_right is on the underlying hyperbola.
   *      If has_right_x && !has_right_y, the right end of the underlying
   *      hyperbola has a vertical asymptote at the x-coordinate 'm_right'.
   *      If !has_right_x && has_right_y, the right end of the underlying
   *      hyperbola has a horizontal asymptote at the y-coordinate 'm_right'.
   */
  Rectangular_hyperbola_with_iso_asymptotes_2(bool a,
                                              const NT& b,
                                              const NT& c,
                                              const NT& d,
                                              const Point_2& left,
                                              const Point_2& right,
                                              bool has_left_x, bool has_left_y,
                                              bool has_right_x, bool has_right_y,
                                              bool is_directed_right,
                                              bool is_continuous) :
    Curve_handle(a, b, c, d, left, right,
                 has_left_x, has_left_y, has_right_x, has_right_y,
                 is_directed_right, is_continuous)
  {}

  /*! Assignment operator.
   */
  Self& operator=(const Self& cv)
  {
    if (this == &cv || this->identical(cv)) return *this;
    Curve_handle::operator=(cv);
    return *this;
  }

private:
  /*! Obtains the representation.
   * \return The representation.
   */
  inline const Curve_rep& rep() const { return (*(this->ptr())); }

  inline Curve_rep& rep() { return (*(this->ptr())); }

public:
  /*! Checks for equality.
   * \return If cv1 and cv2 are identical curves, true; otherwise, false.
   */
  bool operator=(const Self& cv1, const Self& cv2) const
  {
    if (cv1.identical(cv2)) return true;
    return ((cv1.a() == cv2.a()) &&
            (cv1.b() == cv2.b())
            (cv1.c() == cv2.c())
            (cv1.d() == cv2.d())
            (cv1.left() == cv2.left())
            (cv1.right() == cv2.right())
            (cv1.has_left_x() == cv2.has_left_x())
            (cv1.has_left_y() == cv2.has_left_y())
            (cv1.has_right_x() == cv2.has_right_x())
            (cv1.has_right_y() == cv2.has_right_y())
            (cv1.is_directed_right() == cv2.is_directed_right())
            (cv1.is_continuous() == cv2.is_continuous()));
  }

  /*! Obtains the 'a' coefficient of the curve.
   * \return If the curve is a linear object, false; Otherwise, true.
   */
  bool a() const { return rep().m_a; }

  /*! Obtains the 'b' coefficient of the hyperbola.
   * \return The 'b' coefficient.
   */
  const NT& b() const { return rep().m_b; }

  /*! Obtains the 'c' coefficient of the hyperbola.
   * \return The 'c' coefficient.
   */
  const NT& c() const { return rep().m_c; }

  /*! Obtains the 'd' coefficient of the hyperbola.
   * \return The 'd' coefficient.
   */
  const NT& d() const { return rep().m_d; }

  /*! Obtains the left point.
   * \return The left point.
   * \pre The curve has a valid left end point.
   */
  const Point_2& left() const
  {
    CGAL_precondition(has_left());
    return rep().m_left;
  }

  /*! Obtains the right point.
   * \return The right point.
   */
  const Point_2& right() const
  {
    CGAL_precondition(has_right());
    return rep().m_right;
  }

  /*! Indicates whether the curve is directed right.
   * \return If the curve is directed right, true; otherwise false.
   */
  bool is_directed_right() const { return rep().m_is_directed_right; }

  /*! Indicates whether the left point has a valid x-coordinate.
   * \return If the left point has a valid x-coordinate, true; otherwise, false.
   */
  bool has_left_x() const { return rep().m_has_left_x; }

  /*! Indicates whether the left point has a valid y-coordinate.
   * \return If the left point has a valid y-coordinate, true; otherwise, false.
   */
  bool has_left_y() const { return rep().m_has_left_y; }

  /*! Indicates whether the right point has a valid x-coordinate.
   * \return If the right point a valid x-coordinate, true; otherwise, false.
   */
  bool has_right_x() const { return rep().m_has_right_x; }

  /*! Indicates whether the right point has a valid y-coordinate.
   * \return If the right point a valid y-coordinate, true; otherwise, false.
   */
  bool has_right_y() const { return rep().m_has_right_y; }
  
  /*! Indicates whether the curve is continuous.
   * \return If the curve is continuous, true; otherwise false.
   */
  bool is_continuous() const { return rep().m_is_continuous; }

  /*! Indicates whether the left point is a valid left endpoint of the curve.
   * \return If the left point is a valid left endpoint of the curve, true;
   *         otherwise, false.
   */
  bool has_left() const { return has_left_x() && has_left_y(); }

  /*! Indicates whether the right point is a valid right endpoint of the curve.
   * \return If the right point is a valid right endpoint of the curve, true;
   *         otherwise, false.
   */
  bool has_right() const { return has_right_x() && has_right_y(); }

  /*! Indicates whether x-coordinate of the left point is the x-coordinate of
   * an asymptote at the left end of the curve.
   * \return If the x-coordinate of the left point is the x-coordinate of
   *         an asymptote at the left end of the curve, true; otherwise, false.
   */
  bool has_left_asymptote() const { return has_left_x() && !has_left_y(); }

  /*! Indicates whether x-coordinate of the right point is the x-coordinate of
   * an asymptote at the right end of the curve.
   * \return If the x-coordinate of the right point is the x-coordinate of
   *         an asymptote at the right end of the curve, true; otherwise, false.
   */
  bool has_right_asymptote() const { return has_right_x() && !has_right_y(); }

  /*! Checks whether the curve is linear.
   * \return If the curve is linear, true; otherwise, false.
   */
  bool is_linear() const { return !a(); }
  
  /*! Checks whether the object is a segment.
   * \return If the object has two valid endpoints, true; otherwise, false.
   */
  bool is_segment() const
  { return (is_linear() && has_left() && has_right()); }

  /*! Checks whether the object is a ray.
   * \return If the object has exactly one valid endpoint, true;
   * otherwise, false.
   */
  bool is_ray() const
  { return (is_linear() && has_right() != has_left())); }

  /*! Check whether the object is a line.
   * \return If the object has no valid endpoint, true; otherwise, false.
   */
  bool is_line() const
  { return (is_linear() && !has_right() && !has_left()); }

  /*! Checks whether the curve is vertical.
   * \return If the curve is vertical, true; otherwise, false.
   */
  bool is_vertical() const { return (is_linear() && is_zero(rep().m_c)); }

  /*! Checks whether the curve is horizontal.
   * \return If the curve is horizontal, true; otherwise, false.
   */
  bool is_vertical() const { return (is_linear() && is_zero(rep().m_b)); }
  
  /*! Obtains the x-coordinate of the left endpoint.
   * \return The x-coordinate of the left endpoint.
   * \pre The left endpoint of the curve has a valid x-coordinate.
   */
  const NT& left_x() const
  {
    CGAL_precondition(has_left_x());
    return rep().m_left.x();
  }

  /*! Obtains the y-coordinate of the left endpoint.
   * \return The y-coordinate of the left endpoint.
   * \pre The left endpoint of the curve has a valid y-coordinate.
   */
  const NT& left_y() const
  {
    CGAL_precondition(has_left_y());
    return rep().m_left.y();
  }

  /*! Obtains the x-coordinate of the right endpoint.
   * \return The x-coordinate of the right endpoint.
   * \pre The right endpoint of the curve has a valid x-coordinate.
   */
  const NT& right_x() const
  {
    CGAL_precondition(has_right_x());
    return rep().m_right.x();
  }

  /*! Obtains the y-coordinate of the right endpoint.
   * \return The y-coordinate of the right endpoint.
   * \pre The right endpoint of the curve has a valid y-coordinate.
   */
  const NT& right_y() const
  {
    CGAL_precondition(has_right_y());
    return rep().m_right.y();
  }
};

/*!
 * Exporter for the segment class used by the traits-class.
 */
template <typename Kernel, bool Filter, class OutputStream>
OutputStream&
operator<<(OutputStream& os,
           const Rectangular_hyperbola_with_iso_asymptotes_2<Kernel, Filter>& cv)
{
  os << cv.a() << " " << cv.b() << " " << cv.c() << " " << cv.d() << " "
     << cv.left << " " << cv.right << " "
     << cv.has_left_x() << " " << cv.has_left_y() << " "
     << cv.has_right_x() << " " << cv.has_right_y() << " "
     << cv.is_directed_right << " " << cv.is_continuous();
  return os;
}

/*!
 * Importer for the segment class used by the traits-class.
 */
template <typename Kernel, bool Filter, typename InputStream>
InputStream&
operator>>(InputStream& is,
           Rectangular_hyperbola_with_iso_asymptotes_2<Kernel, Filter>& cv)
{
  typedef typename Kernel::NT                          NT;
  typedef Sqrt_extension_point_2<NT, Filter>           Point_2;

  bool a;
  NT b, c, d;
  Point_2 left, right;
  bool has_left_x, has_left_y, has_right_x, has_right_y;
  bool is_directed_right, is_continuous;
  is >> a >> b >> c >> d >> left >> right >>
     >> has_left_x >> has_left_y >> has_right_x >> has_right_y
     >> is_directed_right >> is_continuous;
  cv = Curve_2(a, b, c, d, left, right,
               has_left_x, has_left_y, has_right_x, has_right_y,
               is_directed_right, is_continuous);
  return is;
}

}

#endif
