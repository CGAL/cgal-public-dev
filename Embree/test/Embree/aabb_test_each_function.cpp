#include <iostream>
#include <fstream>
#include <cassert>

#include <CGAL/Embree/AABB_tree.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>


typedef CGAL::Simple_cartesian<float> K;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;
typedef K::Ray_3 Ray;
typedef K::Vector_3 Vector;

typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Surface_mesh<Point> SurfaceMesh;


template<typename T>
class Test{

  typedef typename CGAL::Embree::Triangle_mesh_geometry<T, K> TriangleMesh;
  typedef CGAL::Embree::AABB_tree<TriangleMesh, K> Tree;

public:

  Test(){}
  
  int run(const T& mesh){    
    Tree tree;

    assert(tree.empty());

    tree.insert(mesh);

    typename Tree::Bounding_box bb = tree.bbox();
    std::cout<<"Bounding Box success"<<std::endl;

    size_t tree_size = tree.size();
    if(!tree_size){
      std::cout<<"Tree size is zero"<<std::endl;
      return 1;
    }

    Ray ray(Point(0.1f, 0.1f, -1.0f), Vector(0.0f, 0.0f, 1.0f));

    if(!tree.do_intersect(ray)){
      std::cout<<"No Intersection reported"<<std::endl;
      return 1;
    }

    std::cout<<"Number of Intersected primitives : "<<tree.number_of_intersected_primitives(ray)<<std::endl;

    std::vector<boost::optional<typename Tree::Intersection_and_primitive_id>> intersections;

    tree.all_intersections(ray, std::back_inserter(intersections));
    for (int i=0;i<intersections.size();i++){
      if(intersections[i]){
        Point p = intersections[i]->first;
        std::cout<<"Point of intersection : "<<p<<std::endl;
      }
    }

    boost::optional<typename Tree::Intersection_and_primitive_id> any_intersection = tree.any_intersection(ray);
    boost::optional<typename Tree::Intersection_and_primitive_id> first_intersection = tree.first_intersection(ray);

    std::cout<<"Closest Point to ray origin : "<< tree.closest_point(Point(0.1f, 0.1f, -1.0f))<<std::endl; 

    tree.clear();
    assert(tree.empty());
    return 0; /*Success*/
    }
};


int main(int argc, char const *argv[])
{   
  std::ifstream input("./data/sphere.off");
  Polyhedron polyhedron;
  input >> polyhedron;

  std::ifstream input2("./data/sphere.off");
  SurfaceMesh surface_mesh;
  input2 >> surface_mesh;
  
  Test<Polyhedron> test_polyhedron;
  Test<SurfaceMesh> test_sm;

  std::cout<<"Test on CGAL::Polyhedron_3."<<std::endl;
  test_polyhedron.run(polyhedron);

  std::cout<<std::endl<<std::endl;
  
  std::cout<<"Test on CGAL::SurfaceMesh."<<std::endl;
  test_sm.run(surface_mesh);

  return 0;
}
