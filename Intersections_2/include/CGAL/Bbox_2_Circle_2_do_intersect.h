// Copyright (c) 2008  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
// All rights reserved.
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
// Author(s)     : Camille Wormser, Jane Tournois, Pierre Alliez


#ifndef CGAL_INTERNAL_INTERSECTIONS_2_BBOX_2_CIRCLE_2_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_2_BBOX_2_CIRCLE_2_DO_INTERSECT_H

#include <CGAL/Circle_2.h>
#include <CGAL/Bbox_2.h>

#include <CGAL/number_utils.h>


namespace CGAL {

namespace internal {

    template <class K>
    bool do_intersect(const typename K::Circle_2& sphere,
        const CGAL::Bbox_2& bbox,
        const K&)
    {
        typedef typename K::FT FT;
        typedef typename K::Point_2 Point;
        FT d = FT(0);
        FT distance = FT(0);
		Point center = sphere.center();

		if(center.x() < (FT)bbox.xmin())
		{
			d = (FT)bbox.xmin() - center.x();
			distance += d * d;
		}
		else if(center.x() > (FT)bbox.xmax())
		{
			d = center.x() - (FT)bbox.xmax();
			distance += d * d;
		}

		if(center.y() < (FT)bbox.ymin())
		{
			d = (FT)bbox.ymin() - center.y();
			distance += d * d;
		}
		else if(center.y() > (FT)bbox.ymax())
		{
			d = center.y() - (FT)bbox.ymax();
			distance += d * d;
		}

	
		return distance <= sphere.squared_radius();
    }

    template <class K>
    bool do_intersect(const CGAL::Bbox_2& bbox,
                      const typename K::Circle_2& sphere,
                      const K&)
    { return do_intersect(sphere, bbox, K()); }


} // namespace internal

template<typename K>
bool do_intersect(const CGAL::Bbox_2& a,
                  const Circle_2<K>& b) {
  return K().do_intersect_2_object()(a, b);
}

template<typename K>
bool do_intersect(const Circle_2<K>& a,
                  const CGAL::Bbox_2& b) {
  return K().do_intersect_2_object()(a, b);
}


} //namespace CGAL

#endif  // CGAL_INTERNAL_INTERSECTIONS_2_BBOX_2_CIRCLE_2_DO_INTERSECT_H
