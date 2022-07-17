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

        return (n.depth() <= 1 || func(Vector(mid.x(), mid.y(), mid.z())) < 0.05) && n.depth() <= 3; // custom predicate, can be different
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
        throw("ERROR: isolevel is not between points\n");  // TODO
    }

    const Vector res = r2 * mu + r1 * (1 - mu);
    return Point(res.x(), res.y(), res.z());
}

// returns the two faces along the segment between the two corners
std::pair<int,int> cornersToFaces(std::pair<int,int> corners) {
    int p0 = corners.first, p1 = corners.second;
    int x0 = (p0&4) >> 2, y0 = (p0&2) >> 1, z0 = p0&1;
    int x1 = (p1&4) >> 2, y1 = (p1&2) >> 1, z1 = p1&1;
    std::vector<int> faces;
    // exactly 2 will match
    if(x0 == x1)
        faces.push_back(x0*3);
    if(y0 == y1)
        faces.push_back(y0*3+1);
    if(z0 == z1)
        faces.push_back(z0*3+2);
    assert(faces.size() == 2);
    return std::make_pair(faces[0], faces[1]);
}

// appends the vertices extracted from a given cell to the vector passed as reference
void processNode(const Octree& octree,
    const Edges& edges,
    const Octree::Node& node,
    ImplicitFunction f,
    std::vector<std::vector<Point>>& polyhedron) {

    Octree::Bbox b = octree.bbox(node);

    std::vector<std::pair<Point, std::pair<int,int>>> unorderedPolygonWithFaces;

    for(auto edge : edges.get(node)) {
        std::pair<Point,Point> seg = octree.segment(*edge);
        Vector p1 (seg.first.x(), seg.first.y(), seg.first.z());
        Vector p2 (seg.second.x(), seg.second.y(), seg.second.z());

        if(f(p1) * f(p2) < 0) {
            Octree::Edge minEdge = *(edge->findMinimalEdge());
            std::pair<Point,Point> minSeg = octree.segment(minEdge);
            Vector p1 (minSeg.first.x(), minSeg.first.y(), minSeg.first.z());
            Vector p2 (minSeg.second.x(), minSeg.second.y(), minSeg.second.z());
            unorderedPolygonWithFaces.push_back(std::make_pair(vertex_interpolation(p1, p2, f(p1), f(p2)), 
                cornersToFaces(edge->corners(node))));
        }
    }
    if(unorderedPolygonWithFaces.size() == 0) return;

    std::vector<Point> polygon;

    int ind = 0;
    std::vector<int> visited {ind};
    for(int i = 0; i < unorderedPolygonWithFaces.size(); ++i) {
        polygon.push_back(unorderedPolygonWithFaces[ind].first);
        std::pair<int,int> faces = unorderedPolygonWithFaces[ind].second;
        for(int j = 0; j < unorderedPolygonWithFaces.size(); ++j) {
            auto& el = unorderedPolygonWithFaces[j];
            if(std::find(visited.begin(), visited.end(), j) == visited.end() && 
                (el.second.first == faces.first || el.second.second == faces.first
                || el.second.first == faces.second || el.second.second == faces.second)) {
                ind = j;
                visited.push_back(ind);
                break;
            }
        }
    }

    if(polygon.size() > 0){
        polyhedron.push_back(polygon);
    }
}

int main() {
    ImplicitFunction sphere = [](Vector p){ return sqrt(p.x()*p.x() + p.y()*p.y() + p.z()*p.z()) - 0.5; };
    Point_set points; // This is only here because Orthtree constructor requires it
    Octree octree(Octree::Bbox(-1.1,-1.1,-1.1,1.1,1.1,1.1), points);
    octree.refine(Split_by_closeness(sphere, octree));

    Edges edges(octree, sphere);

    // Add edges to container
    for (Octree::Node node : octree.traverse<Preorder_traversal>()) {
        if(node == octree.root())
            edges.addRoot(node);
        else
            edges.addNode(node);
    }

    std::vector<std::vector<Point>> faces;
    // Traverse octree and process each cell
    // todo: later prepare for parallelization
    for (Octree::Node node : octree.traverse<Preorder_traversal>()) {
        //std::cout << node << std::endl;
        //std::cout << octree.bbox(node) << std::endl;
        if(node.is_leaf())
            processNode(octree, edges, node, sphere, faces);
    }
    //std::cout << octree.depth() << std::endl;
    std::cout << faces.size() << std::endl;

    // writing out resulting points to file
    std::ofstream mesh_out("a.obj");
    int i = 1;
    for (auto f : faces) {
        for(auto p : f)
            mesh_out << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
        mesh_out << "f ";
        for(auto p : f)
            mesh_out << i++ << " ";
        mesh_out << std::endl;
    }

    return 0;
}