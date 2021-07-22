#ifndef CGAL_INTERSECTIONS_BENCHMARK_BOXED_QUERY_H
#define CGAL_INTERSECTIONS_BENCHMARK_BOXED_QUERY_H

#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Ray_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/intersections.h>

using CGAL::Bbox_3;
using CGAL::Ray_3;
using CGAL::Segment_3;
using CGAL::Iso_cuboid_3;

template<typename Q>
class Boxed_query {
public:

  const Q &_query;
  Bbox_3 _box;

  Boxed_query(const Q &query) : _query(query), _box(query.bbox()) {}

  template<typename K>
  Boxed_query(const Q &query, const Iso_cuboid_3<K> &boundary) :
          _query(query), _box(bbox(query, boundary)) {}

  template<typename T, typename K>
  Bbox_3 bbox(const T &query, const Iso_cuboid_3<K> &boundary) {
    return boost::apply_visitor([](auto &&boolean_intersection_result)->Bbox_3 {
      return boolean_intersection_result.bbox();
    }, *CGAL::intersection(query, boundary));
  }
};

//template<typename Q>
//Boxed_query<Q> create_boxed_query(const Q &query) {
//
//}

namespace CGAL {
  template<typename Q>
  bool
  inline
  do_intersect(const Boxed_query<Q> &q,
               const CGAL::Bbox_3 &bbox) {
    return CGAL::do_overlap(q._box, bbox) && do_intersect(q._query, bbox);
  }

} //namespace CGAL

#endif //CGAL_INTERSECTIONS_BENCHMARK_BOXED_QUERY_H
