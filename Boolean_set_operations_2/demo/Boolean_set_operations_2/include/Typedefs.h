// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Apurva Bhatt <response2apurva@gmail.com>
//             Efi Fogel <efifogel@gmain.com>

#ifndef CGAL_TYPEDEFS_H
#define CGAL_TYPEDEFS_H

#include <iostream>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/General_polygon_set_2.h>
//#include <CGAL/Gps_traits_2.h>
//#include <CGAL/Arr_segment_traits_2.h>
//#include <CGAL/Gps_segment_traits_2.h>
//#include <CGAL/Qt/PolygonWithHolesGraphicsItem.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include "QT5/Polygon_set_2.h"
#include "QT5/BezierCurves.h"
//#include <CGAL/Gps_circle_segment_traits_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel        Kernel;

typedef CGAL::Gps_traits_2<CGAL::Arr_segment_traits_2<Kernel> > Linear_traits;
typedef Linear_traits::General_polygon_with_holes_2   Linear_polygon_with_holes;
typedef Linear_traits::Polygon_2                      Linear_polygon;
typedef CGAL::General_polygon_set_2<Linear_traits>    Linear_polygon_set;
typedef Linear_traits::Curve_2                        Linear_curve;

typedef std::vector<Linear_polygon_with_holes>  Linear_region_source_container;
typedef std::vector<Linear_curve>               Linear_boundary_source ;
typedef std::vector<Linear_boundary_source>     Linear_region_source ;

// Circlular polygons

#ifdef CGAL_USE_GMP

  typedef CGAL::Gmpq                     Base_nt;

#else

  typedef CGAL::Quotient<CGAL::MP_Float> Base_nt;

#endif

typedef CGAL::Lazy_exact_nt<Base_nt> Coord_type;


struct Gps_circular_kernel : public CGAL::Cartesian<Coord_type> {};


//typedef Linear_kernel                             Circular_Linear_kernel;//check out for future may generate a bug
/*
typedef CGAL::Simple_cartesian<double>            Linear_kernel ;
typedef CGAL::Polygon_2<Linear_kernel>            Linear_polygon;
typedef CGAL::Polygon_with_holes_2<Linear_kernel> Linear_polygon_with_holes;

typedef Linear_kernel::Point_2 Linear_point ;
*/

typedef Kernel::Point_2                    Circular_Linear_point;
typedef CGAL::Polygon_2<Kernel>            Circular_Linear_polygon;
typedef CGAL::Polygon_with_holes_2<Kernel> Circular_Linear_polygon_with_holes;

typedef CGAL::Gps_circle_segment_traits_2<Kernel>  Circular_traits;
typedef Circular_traits::Curve_2                   Circular_curve;
typedef Circular_traits::X_monotone_curve_2        Circular_X_monotone_curve;
typedef Circular_traits::Point_2                   Circular_point;
typedef Circular_traits::Polygon_2                 Circular_polygon;

//check out the change
//typedef CGAL::General_polygon_with_holes_2<Circular_polygon>   Circular_polygon_with_holes;
typedef Circular_traits::General_polygon_with_holes_2
  Circular_polygon_with_holes;

typedef CGAL::General_polygon_set_2<Circular_traits>    Circular_polygon_set;

typedef std::vector<Circular_polygon_with_holes>
  Circular_region_source_container ;


// Bezier Curves typedefs

#ifdef CGAL_USE_CORE

typedef CGAL::CORE_algebraic_number_traits            Bezier_nt_traits;
typedef Bezier_nt_traits::Rational                    Bezier_rational;
typedef Bezier_nt_traits::Algebraic                   Bezier_algebraic;

struct Bezier_rat_kernel  : public CGAL::Cartesian<Bezier_rational>  {};
struct Bezier_alg_kernel  : public CGAL::Cartesian<Bezier_algebraic> {};

struct Bezier_traits : public CGAL::Arr_Bezier_curve_traits_2<Bezier_rat_kernel, Bezier_alg_kernel, Bezier_nt_traits> {};
  
typedef Bezier_rat_kernel::Point_2                      Bezier_rat_point;
typedef Bezier_traits::Curve_2                          Bezier_curve;
typedef Bezier_traits::X_monotone_curve_2               Bezier_X_monotone_curve;
typedef Bezier_traits::Point_2                          Bezier_point;
typedef CGAL::Gps_traits_2<Bezier_traits>               Bezier_gps_traits;
typedef Bezier_gps_traits::General_polygon_2            Bezier_polygon;
typedef std::vector<Bezier_polygon>                     Bezier_polygon_vector ;
typedef Bezier_gps_traits::General_polygon_with_holes_2 Bezier_polygon_with_holes;
typedef CGAL::General_polygon_set_2<Bezier_gps_traits>  Bezier_polygon_set ;

typedef CGAL::Qt::Bezier_set_graphics_item<Bezier_polygon_set> Bezier_GI;

typedef std::vector<Bezier_curve>                Bezier_boundary_source ;
typedef std::vector<Bezier_boundary_source>      Bezier_region_source ;
typedef std::vector<Bezier_region_source>        Bezier_region_source_container ;

#endif


#endif // CGAL_TYPEDEFS_H

/*
typedef struct Iterator_and_polygons
{
public:
    typedef Linear_traits::Curve_const_iterator        Curve_const_iterator;
    typedef CGAL::General_polygon_set_2<Linear_traits> Polygons;
} Linear_polygon_set;
*/


/*
#ifdef CGAL_USE_GMP

  typedef CGAL::Gmpq                     Base_nt;

#else

  typedef CGAL::Quotient<CGAL::MP_Float> Base_nt;

#endif
*/
//typedef Kernel::FT Base_nt;
//typedef CGAL::Lazy_exact_nt<Base_nt> Coord_type;

//Linear polygons

//typedef CGAL::Simple_cartesian<Coord_type>                       Linear_kernel;


/*
typedef std::vector<Linear_kernel::Point_2>                  Linear_Container;
typedef CGAL::Arr_segment_traits_2<Linear_kernel>            ArrSegmentTraits;
struct Gps_linear_kernel : public CGAL::Cartesian<Coord_type> {};

struct Gps_linear_kernel : public ,
                           public ,
                           public {};
*/
/*
typedef CGAL::Gps_segment_traits_2<CGAL::Exact_predicates_exact_constructions_kernel,
                                   Linear_Container,
                                   ArrSegmentTraits>         Linear_traits;
*/
