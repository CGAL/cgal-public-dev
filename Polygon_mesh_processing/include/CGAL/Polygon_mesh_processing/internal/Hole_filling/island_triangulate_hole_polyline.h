#ifndef ISLAND_TRIANGULATE_HOLE_POLYLINE_H
#define ISLAND_TRIANGULATE_HOLE_POLYLINE_H

namespace CGAL {
namespace internal {




template <
  typename PointRange1,
  typename PointRange2,
  typename Tracer,
  typename WeightCalculator,
  typename Kernel
>
void triangulate_hole_polyline_islands(const PointRange1& points_b,
                                       const PointRange2& points_h,
                                       Tracer& tracer,
                                       const WeightCalculator& WC,
                                       bool use_delaunay_triangulation,
                                       const Kernel&)
{
  typedef Kernel        K;
  typedef typename K::Point_3    Point_3;

  std::vector<Point_3> B(boost::begin(points_b), boost::end(points_b)); // points on the b.
  std::vector<Point_3> H(boost::begin(points_h), boost::end(points_h)); // points on islands


  std::cout<<"ok."<<std::endl;


  //return w;
}




} // namespace internal
} // namsepace CGAL





#endif // ISLAND_TRIANGULATE_HOLE_POLYLINE_H
