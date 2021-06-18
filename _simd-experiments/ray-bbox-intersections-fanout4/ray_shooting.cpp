
#define FANOUT_4

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Ray_3 Ray;

typedef CGAL::Surface_mesh<Point> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;


int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "tetra.off";

  Mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  Tree tree(faces(mesh).first, faces(mesh).second, mesh);

  Ray ray(Point(5, 0.2, 0.2), Vector(1,0,0));

  bool intersection = tree.do_intersect(ray);
  if(intersection)
    {
      std::cout << "there is an intersection" << std::endl;
    }

  std::cout << "done" << std::endl;

  return 0;
}
