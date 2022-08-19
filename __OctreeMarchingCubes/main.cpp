#include <CGAL/Octree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Aff_transformation_3.h>

#include "Octree_mesh_extractor.h"

#include <functional>
#include <algorithm>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT FT;
typedef CGAL::Point_3<Kernel> Point;
typedef CGAL::Vector_3<Kernel> Vector;
typedef CGAL::Point_set_3<Point> Point_set;
typedef Point_set::Point_map Point_map;

typedef Octree_wrapper<Kernel> Octree;
typedef CGAL::Orthtrees::Leaves_traversal Leaves_traversal;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::HalfedgeDS HalfedgeDS;
typedef Polyhedron::Vertex_handle Vertex_handle;

typedef std::function<FT(Vector)> ImplicitFunction;

// Custom refinement predicate, splits cell if mid point of it is "close" to isosurface
struct Split_by_closeness {
    ImplicitFunction func;
    Octree oct;
    Split_by_closeness(ImplicitFunction func, Octree oct) : func(func), oct(oct) {}

    bool operator()(const Octree::Node &n) const {
        CGAL::Bbox_3 b = oct.bbox(oct.root()); // maybe easier to store only this instead of octree
                                               // note: oct.bbox(n) does not give correct result here
        int d = n.depth();

        // calculating midpoint of box of node
        FT x = computeMiddle(b.xmin(), b.xmax(), n.global_coordinates()[0], n.depth());
        FT y = computeMiddle(b.ymin(), b.ymax(), n.global_coordinates()[1], n.depth());
        FT z = computeMiddle(b.zmin(), b.zmax(), n.global_coordinates()[2], n.depth());
        Vector mid{x,y,z}; // note: oct.barycenter(n) does not give correct result here

        return (n.depth() <= 2 || func(Vector(mid.x(), mid.y(), mid.z())) < 0.05) && n.depth() <= 4; // custom predicate, can be different
    }

private:
    FT computeMiddle(FT min, FT max, int c, int d) const {
        return min + (max - min) * (c + 0.5) / (1 << d);
    }
};

int main(int argc, char** argv) {
    ImplicitFunction sphere = [](Vector p) { return sqrt((p.x()-0.1)*(p.x()-0.1) + (p.y()-0.2)*(p.y()-0.2) + p.z()*p.z()) - 0.5; };
    ImplicitFunction cube = [](Vector p) {
        return std::max({std::abs(p.x()), std::abs(p.y()), std::abs(p.z())}) - 0.5;
    };
    ImplicitFunction rotcube = [](Vector p) {
        double cosa = cos(1), sina = sin(1), cosb = cos(0.5), sinb = sin(0.5);
        Kernel::Aff_transformation_3 rot1(
            1.0, 0.0, 0.0,
            0.0, cosa, -sina,
            0.0, sina, cosa);
        Kernel::Aff_transformation_3 rot2(
            cosb, -sinb, 0.0,
            sinb, cosb, 0.0,
            0.0, 0.0, 1.0);
        Vector q = rot2(rot1(p));
        return std::max({std::abs(q.x()), std::abs(q.y()), std::abs(q.z())}) - 0.5;
    };
    ImplicitFunction torus = [](Vector p) {
        return pow(sqrt(p.x()*p.x() + p.y()*p.y()) - 0.5, 2) + p.z()*p.z() - 0.2*0.2;
    };
    ImplicitFunction tanglecube = [](Vector p) {
        double x = 2*p.x(), y = 2*p.y(), z = 2*p.z();
        double x2=x*x, y2=y*y, z2=z*z;
        double x4=x2*x2, y4=y2*y2, z4=z2*z2;
        return x4 - 5*x2 + y4 - 5*y2 + z4 - 5*z2 + 11.8;
    };

    ImplicitFunction func;
    if(argc >= 2 && std::string(argv[1]) == "--function") {
        if(std::string(argv[2]) == "sphere")
            func = sphere;
        else if(std::string(argv[2]) == "cube")
            func = cube;
        else if(std::string(argv[2]) == "rotcube")
            func = rotcube;
        else if(std::string(argv[2]) == "torus")
            func = torus;
        else if(std::string(argv[2]) == "tanglecube")
            func = tanglecube;
        else exit(1);
    }
    else
        func = sphere;

    Octree octree(CGAL::Bbox_3(-1.2,-1.2,-1.2,1.2,1.2,1.2));
    octree.refine(Split_by_closeness(func, octree), [func](const Point& p) { return func(Vector(p[0], p[1], p[2])); });

    CGAL::Octree_mesh_extractor<Kernel> extractor (octree);

    // Traverse octree and process each cell
    // todo: later prepare for parallelization
    for (Octree::Node node : octree.traverse<Leaves_traversal>()) {
        extractor.process_node(node);
    }
    auto vertices = extractor.get_vertices();
    auto faces = extractor.get_faces();
    std::cout << faces.size() << std::endl;

    // writing out resulting points to file
    std::ofstream mesh_out("a.obj");
    int i = 1;
    for(auto p : vertices) {
        mesh_out << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    }
    for (auto f : faces) {
        mesh_out << "f ";
        for(auto p : f)
            mesh_out << (p+1) << " ";
        mesh_out << std::endl;
    }

    return 0;
}