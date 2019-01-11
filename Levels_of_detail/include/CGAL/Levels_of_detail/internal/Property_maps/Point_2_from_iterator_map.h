#ifndef CGAL_LEVEL_OF_DETAIL_POINT_2_FROM_ITERATOR_H
#define CGAL_LEVEL_OF_DETAIL_POINT_2_FROM_ITERATOR_H

// CGAL includes.
#include <CGAL/property_map.h>

namespace CGAL {

namespace Level_of_detail {	
         
template <typename Iterator, typename Point_2, typename PointMap>
struct Point_2_from_iterator_map
{
  typedef Iterator key_type;
  typedef Point_2 value_type;
  typedef const value_type& reference;
  using category   = boost::lvalue_property_map_tag;

  PointMap point_map;

  Point_2_from_iterator_map (PointMap point_map)
    : point_map (point_map)
  { }


  friend reference get (const Point_2_from_iterator_map& pmap, const key_type& k)
  {
    const typename PointMap::value_type& point_3 = get (pmap.point_map, *k);

    // Hack to satisfy Spatial_searching classes which require lvalue
    // property map. Here, as a CGAL::Point_2 is basically a
    // CGAL::Point_3 if we forget the third coordinate (both are based
    // on a cpp11::array<FT, 3/2>), we can just reinterpret the
    // reference and it works well.
    return reinterpret_cast<const Point_2&>(point_3);
  }
    
};

    
} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_POINT_2_FROM_ITERATOR_H
