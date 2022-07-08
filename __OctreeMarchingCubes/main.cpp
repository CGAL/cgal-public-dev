#include <CGAL/Octree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>

#include "EdgeStore.h"

#include <functional>
#include <algorithm>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Point_3<Kernel> Point;
typedef CGAL::Vector_3<Kernel> Vector;
typedef Kernel::FT FT;
typedef CGAL::Point_set_3<Point> Point_set;
typedef Point_set::Point_map Point_map;

typedef CGAL::Octree<Kernel, Point_set, Point_map> Octree;
typedef CGAL::Orthtrees::Preorder_traversal Preorder_traversal;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::HalfedgeDS HalfedgeDS;
typedef Polyhedron::Vertex_handle Vertex_handle;

typedef std::function<FT(Vector)> ImplicitFunction;

typedef CGAL::EdgeStore<Kernel, Point_set, Point_map> Edges;

// Custom refinement predicate, splits cell if mid point of it is "close" to isosurface
struct Split_by_closeness {
    ImplicitFunction func;
    Octree oct;
    Split_by_closeness(ImplicitFunction func, Octree oct) : func(func), oct(oct) {}

    bool operator()(const Octree::Node &n) const {
        Octree::Bbox b = oct.bbox(oct.root()); // maybe easier to store only this instead of octree
                                               // note: oct.bbox(n) does not give correct result here
        int d = n.depth();

        // calculating midpoint of box of node
        FT x = computeMiddle(b.xmin(), b.xmax(), n.global_coordinates()[0], n.depth());
        FT y = computeMiddle(b.ymin(), b.ymax(), n.global_coordinates()[1], n.depth());
        FT z = computeMiddle(b.zmin(), b.zmax(), n.global_coordinates()[2], n.depth());
        Vector mid{x,y,z}; // note: oct.barycenter(n) does not give correct result here

        return (n.depth() <= 1 || func(Vector(mid.x(), mid.y(), mid.z())) < 0.1) && n.depth() <= 3; // custom predicate, can be different
    }

    private:
        FT computeMiddle(FT min, FT max, int c, int d) const {
            return min + (max - min) * (c + 0.5) / (1 << d);
        }
};

Point vertex_interpolation(const Vector& r1, const Vector& r2, const FT d1, const FT d2) {
    FT mu = -1.f;
    FT iso_value = 0.f;

    // don't divide by 0
    if (abs(d2 - d1) < 0.000001) {
        mu = 0.5f;  // if both points have the same value, assume isolevel is in the middle
    } else {
        mu = (iso_value - d1) / (d2 - d1);
    }

    if (mu < 0.f || mu > 1.f) {
        printf("ERROR: isolevel is not between points\n");  // TODO
    }

    const Vector res = r2 * mu + r1 * (1 - mu);
    return Point(res.x(), res.y(), res.z());
}

// appends the vertices extracted from a given cell to the vector passed as reference
void processNode(const Octree& octree, const Edges& edges, const Octree::Node& node, ImplicitFunction f, std::vector<Point>& points) {

    Octree::Bbox b = octree.bbox(node);
    //std::vector<std::pair<Vector,FT>> corners(8);

    // pre-calculating corner values
    /*for (int i = 0; i < 8; ++i) {
        Vector p (i&4 ? b.xmax() : b.xmin(), i&2 ? b.ymax() : b.ymin(), i&1 ? b.zmax() : b.zmin());
        corners[i] = std::make_pair(p, f(p));
    }*/

    for(auto edge : edges.get(node)){
        std::pair<Point,Point> seg = octree.segment(*edge);
        Vector p1 (seg.first.x(), seg.first.y(), seg.first.z());
        Vector p2 (seg.second.x(), seg.second.y(), seg.second.z());

        if(f(p1) * f(p2) < 0) {
            Octree::Edge minEdge = *(edge->findMinimalEdge());
            std::pair<Point,Point> minSeg = octree.segment(minEdge);
            Vector p1 (minSeg.first.x(), minSeg.first.y(), minSeg.first.z());
            Vector p2 (minSeg.second.x(), minSeg.second.y(), minSeg.second.z());
            points.push_back(vertex_interpolation(p1, p2, f(p1), f(p2)));
        }
    }

    // iterating through edges
    /*for (int i = 0; i < 8; ++i) {
        for (int k = 0; k <= 2; ++k) {
            int j = i | (1 << k);
            if (i != j && corners[i].second * corners[j].second < 0) {
                points.push_back(vertex_interpolation(corners[i].first, corners[j].first, corners[i].second, corners[j].second));
                // todo: traverse edge tree and find proper location
            }
        }
    }*/
    // todo: order points to get a closed polygon
}

int main() {
    ImplicitFunction sphere = [](Vector p){ return sqrt(p.x()*p.x() + p.y()*p.y() + p.z()*p.z()) - 0.5; };
    Point_set points; // This is only here because Orthtree constructor requires it
    Octree octree(Octree::Bbox(-1,-1,-1,1,1,1), points);
    octree.refine(Split_by_closeness(sphere, octree));

    Edges edges(octree, sphere);

    // Add edges to container
    for (Octree::Node node : octree.traverse<Preorder_traversal>()) {
        if(node == octree.root())
            edges.addRoot(node);
        else
            edges.addNode(node);
    }

    std::vector<Point> ps;
    // Traverse octree and process each cell
    // todo: later prepare for parallelization
    for (Octree::Node node : octree.traverse<Preorder_traversal>()) {
        //std::cout << node << std::endl;
        //std::cout << octree.bbox(node) << std::endl;

        processNode(octree, edges, node, sphere, ps);
    }
    //std::cout << octree.depth() << std::endl;
    //std::cout << ps.size() << std::endl;

    // writing out resulting points to file
    // expected result: some points will not match, where cells have different size (seems like that happens)
    std::ofstream mesh_out("a.obj");
    for (auto p : ps) {
        mesh_out << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    }

    return 0;
}