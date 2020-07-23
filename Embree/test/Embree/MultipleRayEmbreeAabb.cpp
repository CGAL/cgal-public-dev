#include <iostream>
#include <fstream>

#include <CGAL/Embree/AABB_tree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Random.h>
#include <CGAL/Timer.h>

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
    bool visual = false;

    for (int i = 1; i < argc; i++) {
        if ( (strcmp( "-h", argv[i]) == 0) || (strcmp( "-help", argv[i]) == 0))
            help = true;
        else if ( strcmp( "-o", argv[i]) == 0)
            offFile = true;
        else if ( strcmp( "-v", argv[i]) == 0)
            visual = true;
    }
    if(argc == 1 || help){
        std::cerr << "Usage: " << argv[0] << " <infile> <NumberOfRays> <XPoint> <YPoint> <ZPoint> <-o>[if the input file is .off]  <-v>[for visualisation]"<< std::endl;
        return 0;
    }
    else if(argc<5){
        std::cerr << "Too less arguments."<<std::endl;
        return 0;
    }

    const char* filename =  argv[1];
    std::ifstream input(filename);

    std::ofstream output;
    if(visual) output.open("MultipleRayEmbreeAabb.xyz");

    std::stringstream ss(argv[2]);
    int _numberOfRays = 0;
    ss >> _numberOfRays ;

    float _xPoint, _yPoint, _zPoint;
    ss = std::stringstream(argv[3]);
    ss >> _xPoint;

    ss = std::stringstream(argv[4]);
    ss >> _yPoint;

    ss = std::stringstream(argv[5]);
    ss >> _zPoint;

    Mesh surfaceMesh;
        if(offFile)
        input >> surfaceMesh;
    else
        CGAL::read_ply(input, surfaceMesh);

    Tree tree;
    tree.insert(surfaceMesh);

    Point rayOrigin(_xPoint, _yPoint, _zPoint);

    CGAL::Random rand;
    CGAL::Timer time;

    time.start();

    for(size_t n=0; n!=_numberOfRays; ++n){

        Vector rayDirection(rand.get_double(-1.0, 1.0), rand.get_double(-1.0, 1.0), rand.get_double(-1.0, 1.0));
        Ray ray(rayOrigin, rayDirection);

        boost::optional<Tree::Intersection_and_primitive_id> intersection = tree.first_intersection(ray);

        if(visual){
            if(intersection){
                Point p = intersection->first;
                output<<p<<std::endl;
            }
        }

    }
    time.stop();
    std::cout << "  Function() time: " << time.time() << std::endl;
    output.close();
    return 0;
}
