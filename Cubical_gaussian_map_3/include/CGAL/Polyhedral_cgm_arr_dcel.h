// Copyright (c) 2005, 2006  Tel-Aviv University (Israel).
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
// Author(s)     : Efi Fogel          <efif@post.tau.ac.il>

#ifndef CGAL_POLYHEDRAL_CGM_ARR_DCEL_H
#define CGAL_POLYHEDRAL_CGM_ARR_DCEL_H

#include <CGAL/basic.h>

#include "CGAL/Cgm_arr_dcel.h"

namespace CGAL {

/*! Extend the arrangement vertex */
template <class Point_2>
class Polyhedral_cgm_arr_vertex : public Cgm_arr_vertex<Point_2> {
public:
  /*! Constructor */
  Polyhedral_cgm_arr_vertex() {}
};

/*! Extend the arrangement halfedge */
template <class X_monotone_curve_2>
class Polyhedral_cgm_arr_halfedge : public Cgm_arr_halfedge<X_monotone_curve_2>
{
private:
  /*! A mask of the ids of the original arrangements that contributed the
   * halfedge while performing the minkowski sum.
   * \todo This should be made optional. It is relevant only if the polytope
   * is the result of a Minkowski sum operation, and it is needed only by the 
   * drawing routines.
   */
  unsigned int m_arr_mask;
  
public:
  /*! Constructor */
  Polyhedral_cgm_arr_halfedge() : m_arr_mask(0x0) {}

  /*! Add a arrangement to the mask of the original arrangements in the
   * minkowski sum.
   * \param arr_id the id of the added arrangement
   */
  void add_arr(unsigned int id) { m_arr_mask |= 0x1 << id; }

  /*! Return true iff a given arrangement contributed this halfedge
   * while performing the minkowski sum
   */
  bool is_arr(unsigned int id) const { return m_arr_mask & (0x1 << id); }

  /*! Obtain the mask of the ids of the original arrangements that contributed
   * the halfedge while performing the minkowski sum
   */
  unsigned int get_arr_mask() const { return m_arr_mask; }
};

/*! Extend the arrangement face */
template <class Point_3>
class Polyhedral_cgm_arr_face : public Cgm_arr_face {
private:
  /*! The original point of the polyhedron */
  Point_3 m_point;

  /*! Indicates that the point has been set already */
  bool m_is_set;

public:
  /*! Constructor */
  Polyhedral_cgm_arr_face() : m_is_set(false) { }
    
  /*! Set the 3D point of the original polyhedron */
  void set_point(const Point_3 & point)
  {
    m_point = point;
    m_is_set = true;
  }

  /*! Obtain the 3D point of the original polyhedron */
  const Point_3 & point() const { return m_point; }

  /*! \brief returns true iff the point has been set already */
  bool get_is_set() const { return m_is_set; }

  /*! \brief resets the flag  */
  void set_is_set(bool flag) { m_is_set = flag; }
};

/*! A new dcel builder with CGM features */
template <class Traits>
class Polyhedral_cgm_arr_dcel :
  public CGAL::Arr_dcel_base<Polyhedral_cgm_arr_vertex<typename Traits::Point_2>,
                             Polyhedral_cgm_arr_halfedge<typename Traits::X_monotone_curve_2>,
                             Polyhedral_cgm_arr_face<typename Traits::Point_3> >
{
public:
  /*! Constructor */
  Polyhedral_cgm_arr_dcel() {}
};

} //namespace CGAL

#endif
