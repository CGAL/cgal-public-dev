#ifndef CGAL_ISO_RECTANGLE_2_BBOX_2_INTERSECTION_H_
#define CGAL_ISO_RECTANGLE_2_BBOX_2_INTERSECTION_H_

#include <CGAL/Bbox_2.h>
#include <CGAL/Iso_rectangle_2.h>

namespace CGAL {


template <class R>
inline bool do_intersect(const Iso_rectangle_2<R> &s,
                         const Bbox_2 &box)
{
  typename R::Iso_rectangle_2 ir(box);
  return do_intersect(s, ir);
}


template <class R>
inline bool do_intersect(const Bbox_2 &box,
                         const Iso_rectangle_2<R> &s)
{
  typename R::Iso_rectangle_2 ir(box);
  return do_intersect(s, ir);
}


} //namespace CGAL


#endif /*ISO_RECTANGLE_2_BBOX_2_INTERSECTION_H_*/
