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
#include <CGAL/Root_of_traits.h>
#include <CGAL/Cartesian.h>

namespace CGAL {

// Forward declaration:
template <typename Kernel_>
class Rectangular_hyperbola_with_iso_asymptotes_2;

/*! \class
 * A representation of a rectangular hyperbola with vertical and horizontal
 * asymptotes.
 */
template <typename Kernel_>
class Rectangular_hyperbola_with_iso_asymptotes_2_rep {
public:
  typedef Kernel_                                       Kernel;
  typedef typename Kernel::FT                           NT;
  typedef typename CGAL::Root_of_traits<NT>::Root_of_2  Root_of_2; 
  typedef Cartesian<Root_of_2>                          Root_of_2_kernel; 
  typedef typename Root_of_2_kernel::Point_2            Root_of_2_point; 
  
  friend class Rectangular_hyperbola_with_iso_asymptotes_2<Kernel>;
  
private:
  typedef Rectangular_hyperbola_with_iso_asymptotes_2_rep<Kernel> Self;
  
  bool m_is_linear;             // Indicates whether the curve is linear.
  NT m_b, m_c, m_d;             // xy + bx + cy + d = 0
  
  Root_of_2_point m_left;       // The left point (if exists).
  Root_of_2_point m_right;      // The right point (if exists).
  bool m_has_left_x;            // The left endpoint has a valid x-coordinate. 
  bool m_has_left_y;            // The left endpoint has a valid y-coordinate. 
  bool m_has_right_x;           // The right endpoint has a valid x-coordinate. 
  bool m_has_right_y;           // The right endpoint has a valid y-coordinate. 
  bool m_is_directed_right;     // Is the curve directed right?
  bool m_is_continuous;         // Is the curve continuous?

public:
  /*! Default constructor.
   */
  Rectangular_hyperbola_with_iso_asymptotes_2_rep() :
    m_has_left_x(false),
    m_has_left_y(false),
    m_has_right_x(false),
    m_has_right_y(false)
  {}

  /*! Copy constructor (isn't really used).
   */
  Rectangular_hyperbola_with_iso_asymptotes_2_rep(const Self& other) :
    m_is_linear(other.m_is_linear),
    m_b(other.m_b),
    m_c(other.m_c),
    m_d(other.m_d),
    m_left(other.m_left),
    m_right(other.m_right),
    m_has_left_x(other.m_has_left_x),
    m_has_left_y(other.m_has_left_y),
    m_has_right_x(other.m_has_right_x),
    m_has_right_y(other.m_has_right_y),
    m_is_directed_right(other.m_is_directed_right),
    m_is_continuous(other.m_is_continuous)
  {}
 
  /*! Constructor from all data fields.
   * \param is_linear Indicates whether the curve is linear.
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
  Rectangular_hyperbola_with_iso_asymptotes_2_rep(bool is_linear,
                                                  const NT& b,
                                                  const NT& c,
                                                  const NT& d,
                                                  const Root_of_2_point& left,
                                                  const Root_of_2_point& right,
                                                  bool has_left_x,
                                                  bool has_left_y,
                                                  bool has_right_x,
                                                  bool has_right_y,
                                                  bool is_directed_right,
                                                  bool is_continuous) :
    m_is_linear(is_linear),
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
template <typename Kernel_>
class Rectangular_hyperbola_with_iso_asymptotes_2 :
    public Handle_for<Rectangular_hyperbola_with_iso_asymptotes_2_rep<Kernel_> >
{
public:
  typedef Kernel_                                             Kernel;
  typedef typename Kernel::FT                                 NT; 
  
private:
  typedef Rectangular_hyperbola_with_iso_asymptotes_2<Kernel> Self;

protected:
  typedef Rectangular_hyperbola_with_iso_asymptotes_2_rep<Kernel>
                                                              Curve_rep;
  typedef Handle_for<Curve_rep>                               Curve_handle;
  
public:
  typedef typename Curve_rep::Root_of_2                       Root_of_2; 
  typedef typename Curve_rep::Root_of_2_kernel                Root_of_2_kernel;
  typedef typename Curve_rep::Root_of_2_point                 Root_of_2_point;
  
  /*! Default constructor.
   */
  Rectangular_hyperbola_with_iso_asymptotes_2() : Curve_handle(Curve_rep()) {}
  
  /*! Copy constructor.
   */
  Rectangular_hyperbola_with_iso_asymptotes_2(const Self& cv) :
    Curve_handle(cv) {}
  
  /*! Constructor from all data fields.
   * \param is_linear Indicates whether the curve is linear.
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
  Rectangular_hyperbola_with_iso_asymptotes_2(bool is_linear,
                                              const NT& b,
                                              const NT& c,
                                              const NT& d,
                                              const Root_of_2_point& left,
                                              const Root_of_2_point& right,
                                              bool has_left_x,
                                              bool has_left_y,
                                              bool has_right_x,
                                              bool has_right_y,
                                              bool is_directed_right,
                                              bool is_continuous = true) :
    Curve_handle(Curve_rep(is_linear, b, c, d, left, right,
                           has_left_x, has_left_y, has_right_x, has_right_y,
                           is_directed_right, is_continuous))
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
  bool operator==(const Self& cv2) const
  {
    if (this->identical(cv2)) return true;
    return ((this->is_linear() == cv2.is_linear()) &&
            (this->b() == cv2.b())
            (this->c() == cv2.c())
            (this->d() == cv2.d())
            (this->left() == cv2.left())
            (this->right() == cv2.right())
            (this->has_left_x() == cv2.has_left_x())
            (this->has_left_y() == cv2.has_left_y())
            (this->has_right_x() == cv2.has_right_x())
            (this->has_right_y() == cv2.has_right_y())
            (this->is_directed_right() == cv2.is_directed_right())
            (this->is_continuous() == cv2.is_continuous()));
  }

  /*! Checks whether the curve is linear.
   * \return If the curve is linear, true; otherwise, false.
   */
  bool is_linear() const { return rep().m_is_linear; }

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
  const Root_of_2_point& left() const
  { return rep().m_left; }

  /*! Obtains the right point.
   * \return The right point.
   */
  const Root_of_2_point& right() const
  { return rep().m_right; }

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
  { return (is_linear() && has_right() != has_left()); }

  /*! Check whether the object is a line.
   * \return If the object has no valid endpoint, true; otherwise, false.
   */
  bool is_line() const
  { return (is_linear() && !has_right() && !has_left()); }

  /*! Checks whether the curve is vertical.
   * \return If the curve is vertical, true; otherwise, false.
   */
  bool is_vertical()   const { return (is_linear() && is_zero(rep().m_c)); }

  /*! Checks whether the curve is horizontal.
   * \return If the curve is horizontal, true; otherwise, false.
   */
  bool is_horizontal() const { return (is_linear() && is_zero(rep().m_b)); }
  
  /*! Obtains the x-coordinate of the left endpoint.
   * \return The x-coordinate of the left endpoint.
   * \pre The left endpoint of the curve has a valid x-coordinate.
   */
  const Root_of_2& left_x() const
  {
    CGAL_precondition(has_left_x());
    return rep().m_left.x();
  }

  /*! Obtains the y-coordinate of the left endpoint.
   * \return The y-coordinate of the left endpoint.
   * \pre The left endpoint of the curve has a valid y-coordinate.
   */
  const Root_of_2& left_y() const
  {
    CGAL_precondition(has_left_y());
    return rep().m_left.y();
  }

  /*! Obtains the x-coordinate of the right endpoint.
   * \return The x-coordinate of the right endpoint.
   * \pre The right endpoint of the curve has a valid x-coordinate.
   */
  const Root_of_2& right_x() const
  {
    CGAL_precondition(has_right_x());
    return rep().m_right.x();
  }

  /*! Obtains the y-coordinate of the right endpoint.
   * \return The y-coordinate of the right endpoint.
   * \pre The right endpoint of the curve has a valid y-coordinate.
   */
  const Root_of_2& right_y() const
  {
    CGAL_precondition(has_right_y());
    return rep().m_right.y();
  }
};

/*! \class
 * A representation of an x-monotone rectangular hyperbola with vertical
 * and horizontal asymptotes used by the traits class
 * Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2.
 */
template <typename Kernel_>
class X_monotone_rectangular_hyperbola_with_iso_asymptotes_2 :
    public Rectangular_hyperbola_with_iso_asymptotes_2<Kernel_>
{
private:
  typedef X_monotone_rectangular_hyperbola_with_iso_asymptotes_2<Kernel_> Self;

public:
  typedef typename Curve_rep::Root_of_2                       Root_of_2; 
  typedef typename Curve_rep::Root_of_2_kernel                Root_of_2_kernel;
  typedef typename Curve_rep::Root_of_2_point                 Root_of_2_point;
  
  /*! Default constructor.
   */
  X_monotone_rectangular_hyperbola_with_iso_asymptotes_2() :
    Rectangular_hyperbola_with_iso_asymptotes_2() {}
  
  /*! Copy constructor.
   */
  X_monotone_rectangular_hyperbola_with_iso_asymptotes_2(const Self& cv) :
    Rectangular_hyperbola_with_iso_asymptotes_2(cv) {}
  
  /*! Constructor from all data fields.
   * \param is_linear Indicates whether the curve is linear.
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
  X_monotone_rectangular_hyperbola_with_iso_asymptotes_2
    (bool is_linear,
     const NT& b,
     const NT& c,
     const NT& d,
     const Root_of_2_point& left,
     const Root_of_2_point& right,
     bool has_left_x,
     bool has_left_y,
     bool has_right_x,
     bool has_right_y,
     bool is_directed_right) :
    Rectangular_hyperbola_with_iso_asymptotes_2(is_linear, b, c, d, left, right,
                                                has_left_x, has_left_y,
                                                has_right_x, has_right_y,
                                                is_directed_right, true)
  {}
};

/*!
 * Exporter for the segment class used by the traits-class.
 */
template <typename Kernel, class OutputStream>
OutputStream&
operator<<(OutputStream& os,
           const Rectangular_hyperbola_with_iso_asymptotes_2<Kernel>& cv)
{
  os << cv.is_linear() << " "
     << cv.b() << " " << cv.c() << " " << cv.d() << " "
     << cv.left << " " << cv.right << " "
     << cv.has_left_x() << " " << cv.has_left_y() << " "
     << cv.has_right_x() << " " << cv.has_right_y() << " "
     << cv.is_directed_right << " " << cv.is_continuous();
  return os;
}

/*!
 * Importer for the segment class used by the traits-class.
 */
template <typename Kernel, typename InputStream>
InputStream& 
operator>>(InputStream& is,
           Rectangular_hyperbola_with_iso_asymptotes_2<Kernel>& cv)
{
  typedef typename Kernel::NT                                 NT;
  typedef Rectangular_hyperbola_with_iso_asymptotes_2<Kernel> Curve_2;
  typedef typename Curve_2::Root_of_2_point Root_of_2_point;  

  bool is_linear;
  NT b, c, d;
  Root_of_2_point left, right;
  bool has_left_x, has_left_y, has_right_x, has_right_y;
  bool is_directed_right, is_continuous;
  is >> is_linear >> b >> c >> d >> left >> right
     >> has_left_x >> has_left_y >> has_right_x >> has_right_y
     >> is_directed_right >> is_continuous;
  cv = Curve_2(is_linear, b, c, d, left, right,
               has_left_x, has_left_y, has_right_x, has_right_y,
               is_directed_right, is_continuous);
  return is;
}

}

#endif
