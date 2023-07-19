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
#include <CGAL/intersections.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_surface_mesh.h>
#include <Simulation_mesh.h>
#include <bilinear_patch_3.h>
#include <Ray_3_Bilinear_patch_3_do_intersect.h>
#include <fstream>
#include <collisions.h>
#include <AABB_triangle_pair_primitive.h>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef Kernel::Vector_3                                     Vector;
typedef Kernel::Tetrahedron_3                                Tetrahedron;
typedef Kernel::Ray_3                                        Ray;
typedef Kernel::Aff_transformation_3                         Transform;                        
typedef CGAL::Surface_mesh<Point>                            Mesh;                    
typedef CGAL::Simulation_mesh<Kernel>                        SMesh;
typedef CGAL::BilinearPatchC3<Kernel>                        BilinearPatch;  

int main(int argc, char* argv[])
{
  // Point a = Point(0, 0, 0);
  // Point b = Point(1, 0, 0);
  // Point c = Point(0, 1, 0);
  // Point d = Point(0, 0, 1);

  // BilinearPatch bp = BilinearPatch(a, b, c, d);
  // Tetrahedron t = bp.tetrahedron();

  // Point outside_point = Point(1, 1, 1);
  // Point inside_point = Point(0.25, 0.25, 0.25);
  // std::cout << "Orientation is positive: " << t.orientation() << std::endl;
  // std::cout << "Inside point is inside: " << t.has_on_bounded_side(inside_point) << std::endl;
  // std::cout << "Outside point is outside: " << t.has_on_unbounded_side(outside_point) << std::endl;
  // std::cout << "Should be zero: " << bp.phi(d) << std::endl;
  // std::cout << "Should be false: " << bp.is_planar() << std::endl;
  // Point a = Point(0, 0, 0);
  // Point b = Point(1, 0, 0);
  // Point c = Point(0, 1, 0);
  // Point d = Point(0, 0, 1);
  // Point outside_point = Point(1, 1, 1);
  // Point inside_point = Point(0.25, 0.25, 0.25);

  // BilinearPatch bp = BilinearPatch(a, b, c, d);
  // Tetrahedron t = bp.tetrahedron();
  // Ray r = Ray(a, outside_point);

  // std::cout << "Orientation is positive: " << t.orientation() << std::endl;
  // std::cout << "Inside point is inside: " << t.has_on_bounded_side(inside_point) << std::endl;
  // std::cout << "Outside point is outside: " << t.has_on_unbounded_side(outside_point) << std::endl;
  // std::cout << "Should be true: " << bp.aux_phi(d).do_overlap(0) << std::endl;
  // std::cout << "Should be true: " << bp.has_on(d) << std::endl;
  // std::cout << "Should be false: " << bp.has_on(outside_point) << std::endl;
  // std::cout << "Should be false: " << bp.is_planar() << std::endl;

  const std::string inner_sphere_filename = CGAL::data_file_path("meshes/sphere_internal.off");
  const std::string outer_sphere_filename = CGAL::data_file_path("meshes/sphere_external.off");

  Mesh inner_sphere_mesh;
  Mesh outer_sphere_mesh;

  if(!CGAL::IO::read_polygon_mesh(inner_sphere_filename, inner_sphere_mesh))
  {
    std::cerr << "Invalid input file for internal sphere." << std::endl;
    return EXIT_FAILURE;
  }

  if(!CGAL::IO::read_polygon_mesh(outer_sphere_filename, outer_sphere_mesh))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  SMesh a(inner_sphere_mesh);
  SMesh b(outer_sphere_mesh);

  std::vector<SMesh> meshes{a, b};

  double timestep{0.1};

  bool check = do_collide<SMesh, Kernel>(meshes, timestep);

  std::cout << "Should be false: " << check << std::endl;

  // auto vcm = a.add_property_map<Mesh::Vertex_index, CGAL::IO::Color>("v:color").first;
  // auto ecm = a.add_property_map<Mesh::Edge_index, CGAL::IO::Color>("e:color").first;
  // auto fcm = a.add_property_map<Mesh::Face_index>("f:color", CGAL::IO::white() /*default*/).first; 

  // for(auto v : vertices(a))
  // {
  //   if(v.idx()%2) 
  //   {
  //     Point& p = a.point(v);
  //     Vector& q = a.velocity(v);

  //     Transform transformation(CGAL::TRANSLATION, q);

  //     // swap(q,p.transform(transformation));

  //     std::cout << "P.x current: " << a.point(v).x() << std::endl;
  //     std::cout << "V.x current: " << a.velocity(v).x() << std::endl;
  //   }
  // }

  // for(auto f : a.faces())
  // {
  //   // or the same again, but directly with a range based loop
  //   for(SMesh::Vertex_index vi : a.vertices_around_face(a.halfedge(f)))
  //   {
  //     std::cout << vi << std::endl;
  //   }
  // }

  // CGAL_USE(fcm);



  // CGAL::draw_color(a);

  // Case 0: planar patch
  // Point a = Point(0, 0, 0);
  // Point b = Point(1, 0, 0);
  // Point c = Point(0, 1, 0);
  // Point d = Point(1, 1, 0);
  // Point outside_point = Point(0.25, 0.25, 1);
  // Point outside_point_2 = Point(-0.25, -0.25, 1);
  // Point inside_point = Point(0.25, 0.25, 0);

  // BilinearPatch bp = BilinearPatch(a, b, c, d); 

  // Ray ray_true = Ray(outside_point, inside_point);
  // Ray ray_false = Ray(outside_point, outside_point_2);

  // bool intersects = CGAL::Intersections::internal::do_intersect_odd_parity(bp, ray_true);
  // std::cout << "Case 0, True: " << intersects << std::endl;

  // intersects = CGAL::Intersections::internal::do_intersect_odd_parity(bp, ray_false);
  // std::cout << "Case 0, False: " << intersects << std::endl;

  // // Case 1a: non-planar, ray origin on the bilinear patch
  // a = Point(0, 0, 0);
  // b = Point(1, 0, 0);
  // c = Point(0, 1, 0);
  // d = Point(0, 0, 1);
  // outside_point = Point(-0.1, -0.1, -0.1);
  // outside_point_2 = Point(1, 1, 1);
  // inside_point = Point(0.333, 0.333, 0.333);

  // bp = BilinearPatch(a, b, c, d); 

  // ray_true = Ray(a, outside_point);

  // intersects = CGAL::Intersections::internal::do_intersect_odd_parity(bp, ray_true);
  // std::cout << "Case 1a, True: " << intersects << std::endl;

  // // Case 1b: non-planar, ray origin inside bounding tetrahedron, but not on patch
  // ray_true = Ray(inside_point, outside_point);
  // ray_false = Ray(inside_point, outside_point_2);

  // intersects = CGAL::Intersections::internal::do_intersect_odd_parity(bp, ray_true);
  // std::cout << "Case 1b, True: " << intersects << std::endl;

  // intersects = CGAL::Intersections::internal::do_intersect_odd_parity(bp, ray_false);
  // std::cout << "Case 1b, False: " << intersects << std::endl;

  return EXIT_SUCCESS;
}
