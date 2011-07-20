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
 * A header file for the class templates that represent general and arbitrary
 * rectangular hyperbola with vertical and horizontal asymptotes.
 */

#include <CGAL/Arr_geometry_traits/Sqrt_extension_point_2.h>

namespace CGAL {

/*!
 * \class A representation of a rectangular hyperbola with vertical and
 * horizontal asymptotes used by the traits class
 * Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2.
 */
template <typename Kernel_, bool Filter_>
class Rectangular_hyperbola_with_iso_asymptotes_2 {
public:
  typedef Kernel_                                       Kernel;
  typedef typename Kernel::FT                           NT;
  typedef Sqrt_extension_point_2<NT, Filter_>           Point_2;

private:
    Point_2   m_ps;               // The source point (if exists).
    Point_2   m_pt;               // The target point (if exists).
    bool      m_has_source;       // Is the source point valid
    bool      m_has_target;       // Is the target point valid

   /* y = (a*x + b)/(cx + d) */
   FT m_a, m_b, m_c, m_d;
public:

  /*!
   * Default constructor.
   */
  Rectangular_hyperbola_with_iso_asymptotes_2() :
    has_source(false),
    has_target(false)
  {}
    
  /*!
   * Constructor from two points.
   * \param s The source point.
   * \param t The target point.
   * \pre The two points must not be the same.
   */
  Rectangular_hyperbola_with_iso_asymptotes_2(const FT& a,
                                              const FT& b,
                                              const FT& c,
                                              const FT& d,
                                              const Point_2& ps,
                                              const Point_2& pt) :
    m_ps(ps),
    m_pt(pt),
    m_has_source(true),
    m_has_target(true),
    m_a(a),
    m_b(b),
    m_c(c),
    m_d(d)
  {
    m_is_segment = ((c == 0) && m_has_target && m_has_source);
  }


  /*!
   * Checks whether the object is a segment.
   */
  bool is_segment() const
  { return ((c == 0) && m_has_target && m_has_source); }


  /*!
   * Checks whether the object is a ray.
   */
  bool is_ray() const
  { return ((c == 0) && !m_has_target && !m_has_source); }

  /*!
   * Check whether the object is a line.
   */
  bool is_line() const
  {
    return ((c == 0) &&
            ((m_has_target && !m_has_source) ||
             (!m_has_target && m_has_source)));
  }

  /*!
   * Cast to a line.
   * \pre The linear object is really a line.
   */
  Line_2 line () const
  {
    CGAL_precondition(is_line());
    // TODO: return ...
  }


  /*!
   * Get the source point.
   * \pre The object is a point, a segment or a ray.
   */
  const Point_2& source() const
  {
    CGAL_precondition(has_source);
    return (this->ps);
  }

  /*!
   * Get the target point.
   * \pre The object is a point or a segment.
   */
  const Point_2& target() const
  {
    CGAL_precondition(has_target);
    return (this->pt);
  }
};

/*!
 * Exporter for the segment class used by the traits-class.
 */
template <typename Kernel, class OutputStream>
OutputStream& operator<<(OutputStream& os,
                         const
                         Rectangular_hyperbola_with_iso_asymptotes_2<Kernel>&
                         lobj)
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
template <typename Kernel, typename InputStream>
InputStream& operator>>(InputStream& is,
                        Rectangular_hyperbola_with_iso_asymptotes_2<Kernel>&
                        lobj)
{
  // Read the object type.
  char c;
  do {
    is >> c;
  } while ((c != 'S' && c != 's') &&
           (c != 'R' && c != 'r') &&
           (c != 'L' && c != 'l'));

  // Read the object accordingly.
  if (c == 'S' || c == 's') {
    typename Kernel::Segment_2  seg;
    is >> seg;
    lobj = seg;
  }
  else if (c == 'R' || c == 'r') {
    typename Kernel::Ray_2      ray;
    is >> ray;
    lobj = ray;
  }
  else {
    typename Kernel::Line_2     line;
    is >> line;
    lobj = line;
  }

  return is;
}

}

#endif
