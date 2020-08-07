#include <iostream>
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

  typedef CGAL::Embree::Triangle_mesh_geometry<T, K> TriangleMesh;
  typedef CGAL::Embree::AABB_tree<TriangleMesh, K> Tree;

public:

  int run(std::ifstream input){    
    T Mesh;
    input >> Mesh;
    Tree tree;

    assert(tree.empty());

    tree.insert(surfaceMesh);

    Tree::Bounding_box bb = tree.bbox();
    if (!bb){
      std::cout<<"Bounding Box is empty"<<std::endl;
      return 1;
    }

    size_t tree_size = tree.size();
    if(!tree_size){
      std::cout<<"Tree size is zero"<<std::endl;
      return 1;
    }

    if(!tree.do_intersect()){
      std::cout<<"No Intersection reported"<<std::endl;
      return 1;
    }

    std::cout<<tree.number_of_intersected_primitives()<<std::endl;

    std::vector<boost::optional<Tree::Intersection_and_primitive_id>> intersections;
    
    Point rayOrigin(0.1f, 0.1f, -1.0f);
    Vector rayDirection(0.0f, 0.0f, 1.0f); 
    Ray ray(rayOrigin, rayDirection);

    tree.all_intersections(ray, std::back_inserter(intersections));
    for (int i=0;i<intersections.size();i++){
      if(intersections[i]){
        Point p = intersections[i]->first;
        std::cout<<"Primtive ID : "<<intersections[i].second<<" Point of intersection : "<<p<<std::endl;
      }
    }

    boost::optional<Primitive_id> any_intersection = tree.any_intersection(ray);
    boost::optional<Primitive_id> first_intersection = tree.first_intersection(ray);

    std::cout<<"Closest Point to ray origin : "<< tree.closest_point(rayOrigin)<<std::endl; 

    tree.clear();
    assert(tree.empty());
    return 0; /*Success*/
    }
};


int main(int argc, char const *argv[])
{   
  std::ifstream input("./data/sphere.off");

  Test test;
  test.run(input);

  return 0;
}
