#include <iostream>
#include <fstream>
#include <sstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <CGAL/Timer.h>

#include "RaysGenerate.h"


typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Ray_3 Ray;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;

int main(int argc, char* argv[])
{   
    bool help = false;
    bool offFile = false;
    for (int i = 1; i < argc; i++) {
        if ( (strcmp( "-h", argv[i]) == 0) || (strcmp( "-help", argv[i]) == 0)) 
            help = true;
        else if ( strcmp( "-o", argv[i]) == 0)
            offFile = true;    
    }
    if(argc == 1 || help){
        std::cerr << "Usage: " << argv[0] << " <infile> <NumberOfRays> <XPoint> <YPoint> <ZPoint> <-o>[if the input file is .off]"<< std::endl;
        return 0;
    }
    else if(argc<5){
        std::cerr << "Too less arguments."<<std::endl;
        return 0;
    }

    const char* filename =  argv[1];
    std::ifstream input(filename);

    std::stringstream ss(argv[2]);
    int _numberOfRays = 0;
    ss >> _numberOfRays ;

    double _xPoint, _yPoint, _zPoint;
    ss = std::stringstream(argv[3]);
    ss >> _xPoint;

    ss = std::stringstream(argv[4]);
    ss >> _yPoint;

    ss = std::stringstream(argv[5]);
    ss >> _zPoint;

    Mesh mesh;
    if(offFile)
        input >> mesh;
    else
        CGAL::read_ply(input, mesh);

    CGAL::Timer time;

    time.start();
    
    Tree tree(faces(mesh).first, faces(mesh).second, mesh);
    tree.build(); /*Constructing the Tree*/
    
    time.stop();
    std::cout << "  Construction time: " << time.time() << std::endl;   
    time.reset();

    Point p(_xPoint, _yPoint, _zPoint); /*POINT FOR SHOOTING RAY QUERIES*/

    int numberOfRays = _numberOfRays; /*NUMBER OF RAY QUERIES*/
    RaysGenerate rg(numberOfRays); 
    time.start();

    for (int i=0; i<numberOfRays;i++){
        Vector v(rg.rayDirections[i]._x, rg.rayDirections[i]._y, rg.rayDirections[i]._z);
        Ray ray(p, v);
        
        Ray_intersection intersection = tree.first_intersection(ray);
        // if(intersection){
        //    if(boost::get<Point>(&(intersection->first))){
        //         const Point* p =  boost::get<Point>(&(intersection->first) );
        //         std::cout <<"Point of intersection : "<<  *p << std::endl;
        //     }
        // }
    }
    
    time.stop();
    std::cout << "  Function() time: " << time.time() << std::endl;   

    return 0;
}