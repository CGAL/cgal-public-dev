// Copyright (c) 2010,2011 Tel-Aviv University (Israel).
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
// Author(s): Efi Fogel          <efif@post.tau.ac.il>

#ifndef CGAL_SQRT_ECTENSION_POINT_2_H
#define CGAL_SQRT_ECTENSION_POINT_2_H

/*! \file
 * This header file contains a type that represent a cartesian point in
 * the plane where the x- and y-coordinates of the point are square-root
 * numbers and are stored using the type Sqrt_extension type.
 */

#include <ostream>

#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Sqrt_extension.h>

namespace CGAL {

// Forward declaration:
template <typename NumberType_, bool Filter_> class Sqrt_extension_point_2;

/*! \class
 * A representation of a point whose coordinates are square-root numbers.
 */
template <typename NumberType_, bool Filter_>
class Sqrt_extension_point_2_rep {
public:
  typedef NumberType_                                            NT;
  typedef Filter_                                                Filter;

  typedef Sqrt_extension_point_2_rep<NT, Filter_>                Self;
  typedef Sqrt_extension<NT, NT, Tag_true, Boolean_tag<Filter> > Coord_NT;

  friend class Sqrt_extension_point_2<NumberType_, Filter_>;
  
private:
  Coord_NT m_x;            // The coordinates.
  Coord_NT m_y;

public:
  /*! Default constructor. */
  Sqrt_extension_point_2_rep() : m_x(0), m_y(0) {}

  /*! Constructor of a point with one-root coefficients. 
     This constructor of a point can also be used with rational coefficients
     thanks to convertor of CoordNT. */
  Sqrt_extension_point_2_rep(const Coord_NT& x, const Coord_NT& y) :
    m_x(x),
    m_y(y)
  {}
};

/*! \class
 * A handle for a point whose coordinates are one-root numbers.
 */
template <typename NumberType_, bool Filter_>
class Sqrt_extension_point_2 :
  public Handle_for<Sqrt_extension_point_2_rep<NumberType_, Filter_> >
{
public:
  typedef NumberType_                                   NT;
  typedef Filter_                                       Filter;
  
private:
  typedef Sqrt_extension_point_2<NT, Filter>            Self;
  typedef Sqrt_extension_point_2_rep<NT, Filter>        Point_rep;
  typedef Handle_for<Point_rep>                         Point_handle;

public:
  typedef typename Point_rep::Coord_NT                  Coord_NT;

  /*! Default constructor. */
  Sqrt_extension_point_2() : Point_handle(Point_rep()) {}

  /*! Copy constructor. */
  Sqrt_extension_point_2(const Self& p) : Point_handle(p) {}

  /*! Constructor of a point with one-root coefficients. 
     This constructor of a point can also be used with rational coefficients
     thanks to convertor of CoordNT. */
  Sqrt_extension_point_2(const Coord_NT& x, const Coord_NT& y) :
    Point_handle(Point_rep(x, y))
  {}

  /*! Obtains the x-coordinate. */
  const Coord_NT& x() const { return (this->ptr()->m_x); }

  /*! Obtains the y-coordinate. */
  const Coord_NT& y() const { return (this->ptr()->m_y); }

  /*! Checks for equality. */
  bool equal(const Self& p) const
  {
    if (this->identical(p)) return true;
    return (CGAL::compare(this->ptr()->m_x, p.ptr()->m_x) == EQUAL &&
            CGAL::compare(this->ptr()->m_y, p.ptr()->m_y) == EQUAL);
  }

  /*! Checks for equality. */
  bool operator=(const Self& p, const Self& q) const
  {
    if (p.identical(q)) return true;
    return (CGAL::compare(p.ptr()->m_x, q.ptr()->m_x) == EQUAL &&
            CGAL::compare(p.ptr()->m_y, q.ptr()->m_y) == EQUAL);
  }
  
  bool operator!=(const Self& p, const Self& q) const { return !(p == q); }

  /*! Sets the point coordinates. */
  void set(const NT& x, const NT& y)
  {
    this->copy_on_write();
    this->ptr()->_x = Coord_NT(x);
    this->ptr()->_y = Coord_NT(y);
  }

  /*! Sets the point coordinates. */
  void set(const Coord_NT& x, const Coord_NT& y)
  {
    this->copy_on_write();
    this->ptr()->m_x = x;
    this->ptr()->m_y = y;
  }
};

/*!
 * An exporter for Sqrt_extension_point_2.
 */
template <typename NT, bool Filter>
std::ostream& operator<<(std::ostream& os,
                         const Sqrt_extension_point_2<NT, Filter>& p)
{
  return (os << CGAL::to_double(p.x()) << ' ' << CGAL::to_double(p.y()));
}

}

#nedif
