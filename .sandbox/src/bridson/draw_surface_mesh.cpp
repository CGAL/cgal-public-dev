// // Author(s) : Camille Wormser, Pierre Alliez

#include <iostream>
// #include <list>

// #include <CGAL/Simple_cartesian.h>
// #include <CGAL/AABB_tree.h>
// #include <CGAL/AABB_traits.h>
// #include <CGAL/AABB_triangle_primitive.h>

// typedef CGAL::Simple_cartesian<double> K;

// typedef K::FT FT;
// typedef K::Ray_3 Ray;
// typedef K::Line_3 Line;
// typedef K::Point_3 Point;
// typedef K::Triangle_3 Triangle;

// typedef std::list<Triangle>::iterator Iterator;
// typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
// typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
// typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

// int main()
// {
//     Point a(1.0, 0.0, 0.0);
//     Point b(0.0, 1.0, 0.0);
//     Point c(0.0, 0.0, 1.0);
//     Point d(0.0, 0.0, 0.0);

//     std::list<Triangle> triangles;
//     triangles.push_back(Triangle(a,b,c));
//     triangles.push_back(Triangle(a,b,d));
//     triangles.push_back(Triangle(a,d,c));

//     // constructs AABB tree
//     Tree tree(triangles.begin(),triangles.end());

//     // counts #intersections
//     Ray ray_query(a,b);
//     std::cout << tree.number_of_intersected_primitives(ray_query)
//         << " intersections(s) with ray query" << std::endl;

//     // compute closest point and squared distance
//     Point point_query(2.0, 2.0, 2.0);
//     Point closest_point = tree.closest_point(point_query);
//     std::cerr << "closest point is: " << closest_point << std::endl;
//     FT sqd = tree.squared_distance(point_query);
//     std::cout << "squared distance: " << sqd << std::endl;

//     return EXIT_SUCCESS;
// }

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_surface_mesh.h>
#include <bilinear_patch_3.h>
#include <Ray_3_Bilinear_patch_3_do_intersect.h>
#include <fstream>
typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef Kernel::Tetrahedron_3                                Tetrahedron;
typedef CGAL::Surface_mesh<Point>                            Mesh;
typedef CGAL::BilinearPatchC3<Kernel>                        BilinearPatch;   

int main(int argc, char* argv[])
{
  Point a = Point(0, 0, 0);
  Point b = Point(1, 0, 0);
  Point c = Point(0, 1, 0);
  Point d = Point(0, 0, 1);

  BilinearPatch bp = BilinearPatch(a, b, c, d);
  Tetrahedron t = bp.tetrahedron();

  Point outside_point = Point(1, 1, 1);
  Point inside_point = Point(0.25, 0.25, 0.25);
  std::cout << "Orientation is positive: " << t.orientation() << std::endl;
  std::cout << "Inside point is inside: " << t.has_on_bounded_side(inside_point) << std::endl;
  std::cout << "Outside point is outside: " << t.has_on_unbounded_side(outside_point) << std::endl;
  std::cout << "Should be zero: " << bp.phi(d) << std::endl;
  std::cout << "Should be false: " << bp.is_planar() << std::endl;
  // Mesh sm;
  // if(!CGAL::IO::read_polygon_mesh(filename, sm))
  // {
  //   std::cerr << "Invalid input file." << std::endl;
  //   return EXIT_FAILURE;
  // }
  // CGAL::draw(sm);
  return EXIT_SUCCESS;
}