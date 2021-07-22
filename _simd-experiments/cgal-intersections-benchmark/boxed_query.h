#ifndef CGAL_INTERSECTIONS_BENCHMARK_BOXED_QUERY_H
#define CGAL_INTERSECTIONS_BENCHMARK_BOXED_QUERY_H

#include <CGAL/Bbox_3.h>

using CGAL::Bbox_3;
using CGAL::Ray_3;

template<typename Q>
class Boxed_query {
public:

  const Q &_query;
  Bbox_3 _box;

  Boxed_query(const Q &query) : _query(query), _box(query.bbox()) {}

  Boxed_query(const Q &query, const Bbox_3 &boundary) {
    // TODO
  }
};

namespace CGAL {
  template<typename Q>
  bool
  inline
  do_intersect(const Boxed_query<Q>& q,
               const CGAL::Bbox_3& bbox)
  {
    if (CGAL::do_overlap(q._box, bbox)) return do_intersect(q._query, bbox);
    return false;
  }

} //namespace CGAL

#endif //CGAL_INTERSECTIONS_BENCHMARK_BOXED_QUERY_H
