#include <iostream>
#include <fstream>

#include <CGAL/Embree/AABB_tree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Simple_cartesian<float> K;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;
typedef K::Ray_3 Ray;
typedef K::Vector_3 Vector;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef CGAL::Embree::Triangle_mesh_geometry<Mesh, K> TriangleMesh;
typedef CGAL::Embree::AABB_tree<TriangleMesh, K> Tree;

int main(int argc, char const *argv[])
{
    bool help = false;
    bool offFile = false;

    for (int i = 1; i < argc; i++) {
        if ( strcmp( "-o", argv[i]) == 0)
            offFile = true;
    }

    const char* filename =  argv[1];
    std::ifstream input(filename);

    Mesh surfaceMesh;
    if(offFile)
        input >> surfaceMesh;
    else    /*Support for both Ply and off files. */
        CGAL::read_ply(input, surfaceMesh);

    /*Initialise the tree. The tree constructor automatically creates a single scene and attaches it to the device.*/
    Tree tree;
    tree.insert(surfaceMesh);

    // Configure the ray.
    Point rayOrigin(0.1f, 0.1f, -1.0f);
    Vector rayDirection(0.0f, 0.0f, 1.0f); /*Direction need not be normalized.*/
    Ray ray(rayOrigin, rayDirection);

    std::vector<boost::optional<Tree::Intersection_and_primitive_id>> intersections;

    tree.all_intersections(ray, std::back_inserter(intersections));
        for (int i=0;i<intersections.size();i++){
            if(intersections[i]){
            Point p = intersections[i]->first;
            std::cout<<"Point of intersection : "<<p<<std::endl;
            }
        }        
    return 0;
}
