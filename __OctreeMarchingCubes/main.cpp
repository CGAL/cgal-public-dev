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

typedef CGAL::Orthtree_traits_3<Kernel>::Adjacency Adjacency;

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

        return (n.depth() <= 1 || func(Vector(mid.x(), mid.y(), mid.z())) < 0.05) && n.depth() <= 4; // custom predicate, can be different
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

std::array<Adjacency, 6> sides {Adjacency::LEFT, Adjacency::DOWN, Adjacency::BACK, Adjacency::RIGHT, Adjacency::UP, Adjacency::FRONT};

bool is_inside(Octree::Bbox b, Point p) {
    return b.xmin() <= p.x() && p.x() <= b.xmax()
        && b.ymin() <= p.y() && p.y() <= b.ymax()
        && b.zmin() <= p.z() && p.z() <= b.zmax();
}

std::vector<Point> processFace(const Octree& octree,
    const Edges& edges,
    const Octree::Node& node,
    ImplicitFunction f,
    Adjacency adj,
    Point start) {

    Octree::Node neighbour = node.adjacent_node(adj);

    if (neighbour.is_leaf())
        return std::vector<Point>{};

    std::vector<std::vector<Point>> segments;

    unsigned mask = (!(adj & ~1) ? 1 : adj & ~1);
    std::queue<Octree::Node> toDivide, leafNeighbours;
    toDivide.push(neighbour);
    while (!toDivide.empty()) {
        Octree::Node n = toDivide.front();
        for (unsigned i = 0; i < 8; ++i) {
            if ( !(adj & 1) != !(mask & i) ) { // if the node's child with index i is along the given side of the cell
                if (n[i].is_leaf())
                    leafNeighbours.push(n[i]);
                else
                    toDivide.push(n[i]);
            }
        }
        toDivide.pop();
    }

    int ind = 0, startv = -1, startp = -1;
    while (!leafNeighbours.empty()) {
        Octree::Node child = leafNeighbours.front();
        std::array<Octree::Edge*, 12> child_edges = edges.get(child);
        std::vector<Point> segment;
        for(int i = 0; i < 12; ++i) {
            auto edge = child_edges[i]->findMinimalEdge();
            auto endpoints = octree.segment(*edge);
            auto corners = child_edges[i]->corners(child);
            // if both corners are on the shared face with `node`
            bool c1in = !(corners.first & (4 / mask)) != !(adj & 1);
            bool c2in = !(corners.second & (4 / mask)) != !(adj & 1);
            if (c1in && c2in) {
                auto vals = edge->values();
                if (vals.first * vals.second <= 0) {
                    auto point = vertex_interpolation(
                        Vector(endpoints.first.x(), endpoints.first.y(), endpoints.first.z()),
                        Vector(endpoints.second.x(), endpoints.second.y(), endpoints.second.z()),
                        vals.first, vals.second);
                    segment.push_back(point);
                    if ((point - start).squared_length() < 1e-6 * (endpoints.first - endpoints.second).squared_length()) {
                        startv = ind;
                        startp = segment.size() - 1;
                    }
                }
            }
        }
        if(segment.size() > 0) {
            segments.push_back(segment);
            ++ind;
        }
        leafNeighbours.pop();
    }

    std::vector<Point> polyline;

    int v = startv, p = !startp; // what if it is just an internal polyline?
    for (int i = 0; i < segments.size(); ++i) {
        polyline.push_back(segments[v][p]);
        for (int j = 0; j < segments.size(); ++j) if (v != j) {
            if (segments[j][0] == polyline.back()) { v = j; p = 1; break; }
            if (segments[j][1] == polyline.back()) { v = j; p = 0; break; }
        }
    }
    if(polyline.size() > 0)
        polyline.pop_back();

    return polyline;

}

// returns the vertices extracted from a given cell, ordered as a polyline
std::vector<Point> processNode(const Octree& octree,
    const Edges& edges,
    const Octree::Node& node,
    ImplicitFunction f) {

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
            auto vals = minEdge.values();
            unorderedPolygonWithFaces.push_back(std::make_pair(vertex_interpolation(p1, p2, vals.first, vals.second),
                cornersToFaces(edge->corners(node))));
        }
    }
    if(unorderedPolygonWithFaces.size() == 0) return std::vector<Point> {};

    std::vector<Point> polygon;

    int ind = 0;
    std::vector<int> visited {ind};
    for(int i = 0; i < unorderedPolygonWithFaces.size(); ++i) {
        polygon.push_back(unorderedPolygonWithFaces[ind].first);
        std::pair<int,int> faces = unorderedPolygonWithFaces[ind].second;
        for(int j = 0; j < unorderedPolygonWithFaces.size(); ++j) {
            auto& el = unorderedPolygonWithFaces[j];
            bool face1 = el.second.first == faces.first || el.second.first == faces.second;
            bool face2 = el.second.second == faces.first || el.second.second == faces.second;
            if((std::find(visited.begin(), visited.end(), j) == visited.end()
                || j == 0 && i == unorderedPolygonWithFaces.size() - 1) && (face1 || face2)) {
                int face = face1 ? el.second.first : el.second.second;
                std::vector<Point> internal_points = processFace(octree, edges, node, f, sides[face], polygon.back());
                for (auto it : internal_points) {
                    polygon.push_back(it);
                }
                ind = j;
                visited.push_back(ind);
                break;
            }
        }
    }

    return polygon;
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
        if(node.is_leaf()) {
            std::vector<Point> polygon = processNode(octree, edges, node, sphere);
            if(polygon.size() > 0)
                faces.push_back(polygon);
        }
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