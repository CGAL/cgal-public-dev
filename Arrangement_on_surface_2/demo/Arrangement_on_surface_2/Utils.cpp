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
// $URL: $
// $Id: $
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#include "Utils.h"
#include "ArrangementTypes.h"

CGAL::Bbox_2
Construct_bbox_2_for_Bezier_point::
operator()( Construct_bbox_2_for_Bezier_point::ArrPointType& pt )
{
  if ( pt.is_exact() )
  {
    std::pair<double, double> x_interval = CGAL::to_interval( pt.x() );
    std::pair<double, double> y_interval = CGAL::to_interval( pt.y() );
    CGAL::Bbox_2 res( x_interval.first, y_interval.first,
      x_interval.second, y_interval.second );
    return res;
  }
  else
  {
    CORE::BigRat x_min, x_max, y_min, y_max;
    pt.get_bbox( x_min, x_max, y_min, y_max );
    CGAL::Bbox_2 res( CGAL::to_double(x_min),
      CGAL::to_double(y_min),
      CGAL::to_double(x_max),
      CGAL::to_double(y_max) );
    return res;
  }
}
