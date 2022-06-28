#include <CGAL/Octree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>

#include <functional>
#include <algorithm>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::FT FT;
typedef CGAL::Point_set_3<Point> Point_set;
typedef Point_set::Point_map Point_map;

typedef CGAL::Octree<Kernel, Point_set, Point_map> Octree;
typedef CGAL::Orthtrees::Preorder_traversal Preorder_traversal;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::HalfedgeDS HalfedgeDS;
typedef Polyhedron::Vertex_handle Vertex_handle;

typedef std::function<double(Point)> ImplicitFunction;

// Custom refinement predicate, splits cell if mid point of it is "close" to isosurface
struct Split_by_closeness {
    ImplicitFunction func;
    Octree oct;
    Split_by_closeness(std::function<double(Point)> func, Octree oct) : func(func), oct(oct) {}
    
    bool operator()(const Octree::Node &n) const {
        auto b = oct.bbox(oct.root()); // maybe easier to store only this instead of octree
                                       // note: oct.bbox(n) does not give correct result here
        int d = n.depth();
        // calculating midpoint of box of node
        Point mid(b.xmin() + (b.xmax() - b.xmin()) * (n.global_coordinates()[0] + 0.5) / (1 << d),
            b.ymin() + (b.ymax() - b.ymin()) * (n.global_coordinates()[1] + 0.5) / (1 << d),
            b.zmin() + (b.zmax() - b.zmin()) * (n.global_coordinates()[2] + 0.5) / (1 << d)); // TODO: into function

        return (n.depth() <= 1 || func(mid) < 0.1) && n.depth() <= 3; // custom predicate, can be different
    }
};

// Extracting the polyhedron from an octree cell
// splitting polyhedron optimally into triangles -> later work
template <class HDS>
class Add_polyhedron : public CGAL::Modifier_base<HDS> {
public:
    Add_polyhedron(Octree::Bbox b, ImplicitFunction f) : box(b), func(f) {}

    void operator()( HDS& hds) {
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        // TODO: extract polyhedron
    }
private:
    Octree::Bbox box;
    ImplicitFunction func;
};

int main()
{
    ImplicitFunction sphere = [](Point p){ return sqrt(p.x()*p.x() + p.y()*p.y() + p.z()*p.z()) - 1; };
    Point_set points;
    // adding some dummy points, so that octree has a nonempty bounding box
    // todo: new constructor for orthtree
    points.insert(Point(-1,0,0));
    points.insert(Point(0,-1,0));
    points.insert(Point(0,0,-1));
    points.insert(Point(1,0,0));
    points.insert(Point(0,1,0));
    points.insert(Point(0,0,1));
    Octree octree(points, points.point_map());
    octree.refine(Split_by_closeness(sphere, octree));

    Polyhedron mesh;

    // Traverse octree and process each cell
    // todo: later prepare for parallelization
    for (Octree::Node node : octree.traverse<Preorder_traversal>()) {
        std::cout << node << std::endl;
        std::cout << octree.bbox(node) << std::endl;

        Add_polyhedron<HalfedgeDS> add(octree.bbox(node), sphere);
        mesh.delegate(add);
    }

    return 0;
}
