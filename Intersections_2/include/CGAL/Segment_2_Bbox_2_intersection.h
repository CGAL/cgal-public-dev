// Copyright (c) 2013  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Andreas Fabri


#ifndef CGAL_SEGMENT_2_BBOX_2_INTERSECTION_H
#define CGAL_SEGMENT_2_BBOX_2_INTERSECTION_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Iso_rectangle_2.h>

namespace CGAL {


template <class R>
inline bool do_intersect(const Segment_2<R> &s,
                         const Bbox_2 &box)
{
  typename K::Iso_rectangle_2 ir(box);
  return do_intersect(s, ir);
}


template <class R>
inline bool do_intersect(const Bbox_2 &box,
                         const Segment_2<R> &s)
{
  typename K::Iso_rectangle_2 ir(box);
  return do_intersect(s, ir);
}


} //namespace CGAL



#endif
