#ifndef CGAL_BBOX_3_Iso_cuboid_3_INTERSECTION_H
#define CGAL_BBOX_3_Iso_cuboid_3_INTERSECTION_H


#include <CGAL/Bbox_3.h>
#include <CGAL/Iso_cuboid_3.h>

namespace CGAL {


template <class R>
inline bool do_intersect(const Iso_cuboid_3<R> &s,
                         const Bbox_3 &box)
{
  typename R::Iso_cuboid_3 ir(box);
  return do_intersect(s, ir);
}


template <class R>
inline bool do_intersect(const Bbox_3 &box,
                         const Iso_cuboid_3<R> &s)
{
  typename R::Iso_cuboid_3 ir(box);
  return do_intersect(s, ir);
}


} //namespace CGAL

#endif