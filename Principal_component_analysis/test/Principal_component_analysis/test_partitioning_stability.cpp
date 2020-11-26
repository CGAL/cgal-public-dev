#include <queue>

#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>
#include <CGAL/Simple_cartesian.h>
#include "include/CGAL/test_utils.h"

inline void generate_random_primitive(Point_2 &target)
{
  target = Point_2 (rnd.get_double(), rnd.get_double());
}

inline void generate_random_primitive(Point_3 &target)
{
  target = Point_3 (rnd.get_double(), rnd.get_double(), rnd.get_double());
}

inline void generate_random_primitive(const Point_2 &min,
                                      const Point_2 &max,
                                      Point_2 &target)
{
  double x = rnd.get_double(min.x(), max.x());
  double y = rnd.get_double(min.y(), max.y());
  target = Point_2(x,y);
}

inline void generate_random_primitive(const Point_3 &min,
                                      const Point_3 &max,
                                      Point_3 &target)
{
  double x = rnd.get_double(min.x(), max.x());
  double y = rnd.get_double(min.y(), max.y());
  double z = rnd.get_double(min.z(), max.z());
  target = Point_3(x,y,z);
}

inline void generate_random_primitive(Segment_2 &target)
{
  Point_2 p0, p1;
  generate_random_primitive(p0);
  generate_random_primitive(p1);
  target = Segment_2 (p0,p1);
}

inline void generate_random_primitive(Segment_3 &target)
{
  Point_3 p0, p1;
  generate_random_primitive(p0);
  generate_random_primitive(p1);
  target = Segment_3 (p0,p1);
}

inline void generate_random_primitive(Triangle_2 &target)
{
  // random equilateral triangle
  Point_2 center;
  generate_random_primitive(center);
  double length = rnd.get_double(0, std::abs(rnd.get_double()));
  double height = length * std::sqrt(3.) / 2.;
  target = Triangle_2 (center+ Vector_2 (0, 2. * height / 3.),
                       center+ Vector_2 (-length/2., -height/3.),
                       center+ Vector_2 (length/2., -height/3.));
}

inline void generate_random_primitive(Triangle_3 &target)
{
  // random equilateral triangle
  Point_3 center;
  generate_random_primitive(center);
  double length = rnd.get_double(0, std::abs(rnd.get_double()));
  double height = length * std::sqrt(3.) / 2.;
  target = Triangle_3 (center+ Vector_3 (0, 2. * height / 3., 0),
                       center+ Vector_3 (-length/2., -height/3., -length/2.),
                       center+ Vector_3 (length/2., -height/3., length/2.));
}

inline void generate_random_primitive(Tetrahedron_3 &target)
{
  // random equilateral tetrahedron
  Point_3 center;
  generate_random_primitive(center);
  double length = rnd.get_double(0, std::abs(rnd.get_double()));
  double base_height = length * std::sqrt(3.) / 2.;
  double height =  length * std::sqrt(2./3.);
  target = Tetrahedron_3 (center+ Vector_3 (0, 2. * base_height / 3., -height/4.),
                          center+ Vector_3 (-length/2., -base_height/3., -height/4.),
                          center+ Vector_3 (length/2., -base_height/3., -height/4.),
                          center+ Vector_3 (0, 0, 3.*height/4.));
}

inline void generate_random_primitive(Iso_rectangle_2 &target)
{
  Point_2 center;
  generate_random_primitive(center);
  double width = rnd.get_double(0, std::abs(rnd.get_double()));
  double height = rnd.get_double(0, std::abs(rnd.get_double()));
  target = Iso_rectangle_2 (center, center+Vector_2(width,height));
}

inline void generate_random_primitive(Iso_cuboid_3 &target)
{
  Point_3 center;
  generate_random_primitive(center);
  double width = rnd.get_double(0, std::abs(rnd.get_double()));
  double height = rnd.get_double(0, std::abs(rnd.get_double()));
  double depth = rnd.get_double(0, std::abs(rnd.get_double()));
  target = Iso_cuboid_3 (center, center+Vector_3(width,height,depth));
}

template <typename OutputIterator>
void random_split_of_primitive(const Segment_2 &source,
                               OutputIterator target)
{
  CGAL::Random_points_on_segment_2<Point_2> g(source[0], source[1], rnd);
  std::vector<Point_2> random;
  random.reserve(1);
  std::copy_n(g, 1, std::back_inserter(random));
  *target++ = Segment_2(source[0], random[0]);
  *target++ = Segment_2(random[0], source[1]);
}

template <typename OutputIterator>
void random_split_of_primitive(const Segment_3 &source,
                               OutputIterator target)
{
  CGAL::Random_points_on_segment_3<Point_3> g(source[0], source[1], rnd);
  std::vector<Point_3> random;
  random.reserve(1);
  std::copy_n(g, 1, std::back_inserter(random));
  *target++ = Segment_3(source[0], random[0]);
  *target++ = Segment_3(random[0], source[1]);
}

template <typename OutputIterator>
void random_split_of_primitive(const Triangle_2 &source,
                               OutputIterator target)
{
  CGAL::Random_points_in_triangle_2<Point_2> g(source, rnd);
  std::vector<Point_2> random;
  random.reserve(1);
  std::copy_n(g, 1, std::back_inserter(random));
  for(int i = 0; i < 3; i++)
    *target++ = Triangle_2(source[i], source[(i+1)%3], random[0]);
}

template <typename OutputIterator>
void random_split_of_primitive(const Triangle_3 &source,
                               OutputIterator target)
{
  CGAL::Random_points_in_triangle_3<Point_3> g(source, rnd);
  std::vector<Point_3> random;
  random.reserve(1);
  std::copy_n(g, 1, std::back_inserter(random));
  for(int i = 0; i < 3; i++)
    *target++ = Triangle_3(source[i], source[(i+1)%3], random[0]);
}

template <typename OutputIterator>
void random_split_of_primitive(const Tetrahedron_3 &source,
                               OutputIterator target)
{
  CGAL::Random_points_in_tetrahedron_3<Point_3> g(source, rnd);
  std::vector<Point_3> random;
  random.reserve(1);
  std::copy_n(g, 1, std::back_inserter(random));
  for(int i = 0; i < 4; i++)
    *target++ = Tetrahedron_3(source[i], source[(i+1)%4],
                              source[(i+2)%4], random[0]);
}

template <typename OutputIterator>
void random_split_of_primitive(const Iso_rectangle_2 &source,
                               OutputIterator target)
{
  Point_2 random;
  generate_random_primitive(source.min(), source.max(), random);
  for(int i = 0; i < 4; i++)
    *target++ = Iso_rectangle_2(source[i], random);
}


template <typename OutputIterator>
void random_split_of_primitive(const Iso_cuboid_3 &source,
                               OutputIterator target)
{
  Point_3 random;
  generate_random_primitive(source.min(), source.max(), random);
  for(int i = 0; i < 8; i++)
    *target++ = Iso_cuboid_3(source[i], random);
}

template <typename Primitive, typename OutputIterator>
void generate_random_points_on_primitive_bounded_by_bbox (
        Primitive primitive,
        const Iso_rectangle_2 &bbox,
        OutputIterator target)
{
  for(int i = 0; i < 100; i++)
  {
    Point_2 random;
    Point_2 bbox_min(bbox.xmin(), bbox.ymin());
    Point_2 bbox_max(bbox.xmax(), bbox.ymax());
    generate_random_primitive(bbox_min, bbox_max, random);
    Point_2 proj = primitive.projection(random);
    *target++ = proj;
  }
}

template <typename Primitive, typename OutputIterator>
void generate_random_points_on_primitive_bounded_by_bbox (
        Primitive primitive,
        const Iso_cuboid_3 &bbox,
        OutputIterator target)
{
  for(int i = 0; i < 100; i++)
  {
    Point_3 random;
    Point_3 bbox_min(bbox.xmin(), bbox.ymin(), bbox.zmin());
    Point_3 bbox_max(bbox.xmax(), bbox.ymax(), bbox.zmax());
    generate_random_primitive(bbox_min, bbox_max, random);
    Point_3 proj = primitive.projection(random);
    *target++ = proj;
  }
}

template <typename Object, typename Fitted, int dim>
void test_partitioning_stability(const Object &initial_primitive,
                                 const std::vector<Object> &partitioning,
                                 const CGAL::Dimension_tag<2>& /* ambient dimension */)
{
  // Compute the initial Fitted
  std::vector<Object> initial_primitive_v = {initial_primitive};
  Fitted initial_fitted;
  linear_least_squares_fitting_2(initial_primitive_v.begin(),
                                 initial_primitive_v.end(),
                                 initial_fitted, CGAL::Dimension_tag<dim>());

  // Compute the partition Fitted
  Fitted partitioning_fitted;
  linear_least_squares_fitting_2(partitioning.begin(),
                                 partitioning.end(),
                                 partitioning_fitted, CGAL::Dimension_tag<dim>());

  // Check if the initial Fitted and the partition Fitted are almost the same
  std::vector<Point_2> random_points_on_partitioning_fitted;
  generate_random_points_on_primitive_bounded_by_bbox(partitioning_fitted,
                                                      bbox_2(partitioning.begin(), partitioning.end()),
                                                      std::back_inserter(random_points_on_partitioning_fitted));
  assert_quality(random_points_on_partitioning_fitted, initial_fitted);
}

template <typename Object, typename Fitted, int dim>
void test_partitioning_stability(const Object &initial_primitive,
                                 const std::vector<Object> &partitioning,
                                 const CGAL::Dimension_tag<3>& /* ambient dimension */)
{
  // Compute the initial Fitted
  std::vector<Object> initial_primitive_v = {initial_primitive};
  Fitted initial_fitted;
  linear_least_squares_fitting_3(initial_primitive_v.begin(),
                                 initial_primitive_v.end(),
                                 initial_fitted, CGAL::Dimension_tag<dim>());

  // Compute the partition Fitted
  Fitted partitioning_fitted;
  linear_least_squares_fitting_3(partitioning.begin(),
                                 partitioning.end(),
                                 partitioning_fitted, CGAL::Dimension_tag<dim>());

  // Check if the initial Fitted and the partition Fitted are almost the same
  std::vector<Point_3> random_points_on_partitioning_fitted;
  generate_random_points_on_primitive_bounded_by_bbox(partitioning_fitted,
                                                      bbox_3(partitioning.begin(), partitioning.end()),
                                                      std::back_inserter(random_points_on_partitioning_fitted));
  assert_quality(random_points_on_partitioning_fitted, initial_fitted);
}

template <typename Object, typename Fitted, int dim>
void test_partitioning_stability()
{
  // Generating an initial random Object and compute its Fitted
  Object initial_primitive;
  generate_random_primitive(initial_primitive);

  // Partitioning of the initial Object
  std::queue<Object> queue;
  queue.push(initial_primitive);
  int iteration = 0;
  while (iteration < 10)
  {
    Object primitive = queue.front();
    queue.pop();
    std::vector<Object> split;
    random_split_of_primitive(primitive, std::back_inserter(split));
    for (Object s : split)
      queue.push(s);
    iteration++;
  }
  std::vector<Object> partitioning;
  while (!queue.empty())
  {
    partitioning.push_back(queue.front());
    queue.pop();
  }

  // Compute both Fitted and check if they are almost identical
  test_partitioning_stability<Object, Fitted, dim>(initial_primitive, partitioning,
                              typename CGAL::Ambient_dimension<Object, Kernel>::type());

}

int main()
{

  std::cerr << "Partitioning stability test with seed " << rnd.get_seed() << std::endl;

  std::cerr << std::endl << "=== 2D ===" << std::endl << std::endl;

  std::cerr << "[Testing line fitting on Segment_2 objects]" << std::endl;
  test_partitioning_stability<Segment_2, Line_2, 1> ();

  std::cerr << "[Testing line fitting on Triangle_2 objects]" << std::endl;
  test_partitioning_stability<Triangle_2, Line_2, 2> ();

  std::cerr << "[Testing line fitting on Iso_rectangle_2 objects]" << std::endl;
  test_partitioning_stability<Iso_rectangle_2, Line_2, 2> ();

  std::cerr << std::endl << "=== 3D ===" << std::endl << std::endl;

  std::cerr << "[Testing line fitting on Segment_3 objects]" << std::endl;
  test_partitioning_stability<Segment_3, Line_3, 1> ();

  std::cerr << "[Testing plane fitting on Segment_3 objects]" << std::endl;
  test_partitioning_stability<Segment_3, Plane_3, 1> ();

  std::cerr << "[Testing line fitting on Triangle_3 objects]" << std::endl;
  test_partitioning_stability<Triangle_3, Line_3, 2> ();

  std::cerr << "[Testing plane fitting on Triangle_3 objects]" << std::endl;
  test_partitioning_stability<Triangle_3, Plane_3, 2> ();

  std::cerr << "[Testing line fitting on Tetrahedron_3 objects]" << std::endl;
  test_partitioning_stability<Tetrahedron_3, Line_3, 3> ();

  std::cerr << "[Testing plane fitting on Tetrahedron_3 objects]" << std::endl;
  test_partitioning_stability<Tetrahedron_3, Plane_3, 3> ();

  std::cerr << "[Testing line fitting on Iso_cuboid_3 objects]" << std::endl;
  test_partitioning_stability<Iso_cuboid_3, Line_3, 3> ();

  std::cerr << "[Testing plane fitting on Iso_cuboid_3 objects]" << std::endl;
  test_partitioning_stability<Iso_cuboid_3, Plane_3, 3> ();

  return EXIT_SUCCESS;
}
