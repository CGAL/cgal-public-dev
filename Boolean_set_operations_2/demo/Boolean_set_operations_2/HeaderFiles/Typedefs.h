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
// Author(s) : Saar Katz <kats.saar@gmail.com>

#ifndef CGAL_TYPEDEFS_H
#define CGAL_TYPEDEFS_H

#include <CGAL/General_polygon_set_2.h>
#include "CGAL/Exact_predicates_exact_constructions_kernel.h"
#include "CGAL/Qt/PolygonWithHolesGraphicsItem.h"

typedef CGAL::Exact_predicates_exact_constructions_kernel         Kernel;
typedef Kernel::Point_2                                           Point;
typedef CGAL::Polygon_2<Kernel>                                   Polygon;
typedef CGAL::Polygon_with_holes_2<Kernel>                        PolygonWithHoles;
typedef CGAL::Qt::PolygonWithHolesGraphicsItem<PolygonWithHoles>  PolygonGraphicsItem;

typedef	std::vector<Point>  Container;

//typedef CGAL::General_polygon_set_2<Kernel> PolygonSet;

#endif // CGAL_TYPEDEFS_H