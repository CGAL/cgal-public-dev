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
// $URL$
// $Id$
// 
//
// Author(s)     : Kapelushnik Lior <liorkape@post.tau.ac.il>

/*! \file
 * spherical arrangements of none intersecting arcs of great circles on a sphere
 */

#ifndef CGAL_SPHERE_TRAITS_H
#define CGAL_SPHERE_TRAITS_H

#include <CGAL/Sphere_arc.h>

#include <ostream>

namespace CGAL {

/*
  traits class for a sphere,
  holds curves that represents arcs of great circles on a sphere

  Kernel_ - the kernel with number types for spherical arcs endpoints directions
 */
template <class Kernel_>
class Sphere_traits : public Kernel_ {
public:
  // spherical traits types decleration
  typedef Kernel_                Kernel;
  typedef typename Kernel::Direction_3    Direction_3;
  typedef typename Kernel::Vector_3      Vector_3;

  typedef typename Kernel::Point_3      Point_3;
  // the curve is a spherical arc on a great circle
  typedef Sphere_arc<Kernel_>          X_monotone_curve_2;
  typedef X_monotone_curve_2           Curve_2;

private:
};

} //namespace CGAL

#endif
