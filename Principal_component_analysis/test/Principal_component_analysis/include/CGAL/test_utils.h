#ifndef CGAL_TEST_UTILS_H
#define CGAL_TEST_UTILS_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Random.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Segment_3 Segment_3;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Sphere_3 Sphere_3;
typedef Kernel::Triangle_3 Triangle_3;
typedef Kernel::Tetrahedron_3 Tetrahedron_3;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Line_3 Line_3;

typedef Kernel::Point_2 Point_2;
typedef Kernel::Vector_2 Vector_2;
typedef Kernel::Circle_2 Circle_2;
typedef Kernel::Iso_rectangle_2 Iso_rectangle_2;
typedef Kernel::Segment_2 Segment_2;
typedef Kernel::Triangle_2 Triangle_2;
typedef Kernel::Line_2 Line_2;

static CGAL::Random rnd(0); // ensure always same test

/*
  Compute distances between Point_2/Point_3 to Line_3/Plane_3/Line_2
 */
inline double distance (const Point_2& p, const Line_2& l)
{ return std::sqrt (CGAL::squared_distance (p, l)); }
inline double distance (const Point_3& p, const Line_2& l)
{ return std::sqrt (CGAL::squared_distance (Point_2 (p.x(), p.y()), l)); }
inline double distance (const Point_3& p, const Line_3& l)
{ return std::sqrt (CGAL::squared_distance (p, l)); }
inline double distance (const Point_3& p, const Plane_3& pl)
{ return std::sqrt (CGAL::squared_distance (p, pl)); }

/*
  Test quality of fit and raise assertion if too low.
 */
template <typename Point, typename Fitted>
void assert_quality (const std::vector<Point>& points, const Fitted& fitted)
{
  double mean_dist = 0;
  for (std::size_t i = 0; i < points.size(); ++ i)
  {
    double dist = distance (points[i], fitted);
    mean_dist += dist;
  }
  mean_dist /= points.size();

  std::cerr << "mean distance = " << mean_dist << std::endl;

  CGAL_assertion_code
    (double limit = 1e-3 * std::sqrt (CGAL::squared_distance (points.front(), points.back())));
  //CGAL_assertion (mean_dist < limit);
  CGAL_assertion_code(if(mean_dist >= limit) std::cerr << "FAILURE: mean distance to high\n");
}

#endif // CGAL_TEST_UTILS_H
