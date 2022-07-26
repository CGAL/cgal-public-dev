#include <CGAL/Octree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include<CGAL/Aff_transformation_3.h>

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
typedef CGAL::Octree_edge<Kernel> Edge;

typedef CGAL::Orthtree_traits_3<Kernel>::Adjacency Adjacency;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::HalfedgeDS HalfedgeDS;
typedef Polyhedron::Vertex_handle Vertex_handle;

typedef std::function<FT(Vector)> ImplicitFunction;

typedef CGAL::Edge_store<Kernel, Point_set, Point_map> Edges;

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

// returns the two faces along the segment between the two corners
// this function is really messy, but seems working; maybe it would be better to use tables instead
// would be better to use the same ordering as Orthtree_traits_3<>::Adjacency
std::pair<int,int> cornersToFaces(std::pair<int,int> corners) {
    int p0 = corners.first, p1 = corners.second;
    int x0 = (p0&4) >> 2, y0 = (p0&2) >> 1, z0 = p0&1;
    int x1 = (p1&4) >> 2, y1 = (p1&2) >> 1, z1 = p1&1;
    std::vector<int> faces {-1, -1};
    // exactly 2 will match
    if(x0 == x1) {
        if(z0 < z1 || y0 > y1)
            faces[0] = x0*3;
        else
            faces[1] = x0*3;
    }
    if(y0 == y1) {
        if(x0 < x1 || z0 > z1)
            faces[0] = y0*3+1;
        else
            faces[1] = y0*3+1;
    }
    if(z0 == z1) {
        if(y0 < y1 || x0 > x1)
            faces[0] = z0*3+2;
        else
            faces[1] = z0*3+2;

    }
    if(p0 == 0 || p0 == 7 || p1 == 0 || p1 == 7) {
        auto tmp = faces[0];
        faces[0] = faces[1];
        faces[1] = tmp;
    }
    assert(faces[0] != -1 && faces[1] != -1);
    return std::make_pair(faces[0], faces[1]);
}

std::array<Adjacency, 6> sides {Adjacency::LEFT, Adjacency::DOWN, Adjacency::BACK, Adjacency::RIGHT, Adjacency::UP, Adjacency::FRONT};

bool is_inside(Octree::Bbox b, Point p) {
    return b.xmin() <= p.x() && p.x() <= b.xmax()
        && b.ymin() <= p.y() && p.y() <= b.ymax()
        && b.zmin() <= p.z() && p.z() <= b.zmax();
}

std::map<Octree::Node, std::vector<size_t>> pointsInCells;
std::map<Edge, size_t> edge_points_in_mesh;

size_t addPoint(Octree::Node n, const Edge& e, std::vector<Point>& points) {
    auto i = edge_points_in_mesh.find(e);
    if (i != edge_points_in_mesh.end())
        return i->second;

    Point p = e.extract_isovertex();
    for (auto it : sides) {
        Octree::Node neighbour = n.adjacent_node(it);
        for (auto ind : pointsInCells[neighbour]) {
            if ((points[ind] - p).squared_length() < 1e-6) {
                pointsInCells[n].push_back(ind);
                return ind;
            }
        }
    }
    size_t ind = points.size();
    pointsInCells[n].push_back(ind);
    edge_points_in_mesh[e] = ind;
    points.push_back(p);
    return ind;
}

std::vector<size_t> processFace(const Octree& octree,
    const Edges& edges,
    const Octree::Node& node,
    ImplicitFunction f,
    Adjacency adj,
    Point start,
    std::vector<Point>& points) {

    Octree::Node neighbour = node.adjacent_node(adj);

    if (neighbour.is_null() || neighbour.is_leaf())
        return std::vector<size_t>{};

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

    std::vector<std::vector<size_t>> segments;

    int ind = 0, startv = -1, startp = -1;
    while (!leafNeighbours.empty()) {
        Octree::Node child = leafNeighbours.front();
        std::array<Edge*, 12> child_edges = edges.get(child);
        std::vector<size_t> segment;
        for(int i = 0; i < 12; ++i) {
            auto edge = child_edges[i]->find_minimal_edge();
            auto endpoints = edge->segment();
            auto corners = child_edges[i]->corners(child);
            // if both corners are on the shared face with `node`
            bool c1in = !(corners.first & (4 / mask)) != !(adj & 1);
            bool c2in = !(corners.second & (4 / mask)) != !(adj & 1);
            if (c1in && c2in) {
                auto vals = edge->values();
                if (vals.first * vals.second <= 0) {
                    auto point = edge->extract_isovertex();
                    segment.push_back(addPoint(node, *edge, points));
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

    std::vector<size_t> polyline;

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
std::vector<size_t> processNode(const Octree& octree,
    const Edges& edges,
    const Octree::Node& node,
    ImplicitFunction f,
    std::vector<Point>& points) {

    Octree::Bbox b = octree.bbox(node);

    std::vector<std::pair<size_t, std::pair<int,int>>> unorderedPolygonWithFaces;

    for(auto edge : edges.get(node)) {
        std::pair<Point,Point> seg = edge->segment();
        Vector p1 (seg.first.x(), seg.first.y(), seg.first.z());
        Vector p2 (seg.second.x(), seg.second.y(), seg.second.z());

        if(edge->values().first * edge->values().second < 0) {
            Edge minEdge = *(edge->find_minimal_edge());
            std::pair<Point,Point> minSeg = minEdge.segment();
            Vector p1 (minSeg.first.x(), minSeg.first.y(), minSeg.first.z());
            Vector p2 (minSeg.second.x(), minSeg.second.y(), minSeg.second.z());
            auto vals = minEdge.values();
            auto corners = edge->corners(node);
            if (vals.second < 0) { auto tmp = corners.first; corners.first = corners.second; corners.second = tmp;  }
            unorderedPolygonWithFaces.push_back(std::make_pair(
                addPoint(node, minEdge, points),
                cornersToFaces(corners)));
        }
    }
    if(unorderedPolygonWithFaces.size() == 0) return std::vector<size_t> {};

    std::vector<size_t> polygon;

    int ind = 0;
    std::vector<int> visited {ind};
    for(int i = 0; i < unorderedPolygonWithFaces.size(); ++i) {
        polygon.push_back(unorderedPolygonWithFaces[ind].first);
        std::pair<int,int> faces = unorderedPolygonWithFaces[ind].second;
        for(int j = 0; j < unorderedPolygonWithFaces.size(); ++j) {
            auto& el = unorderedPolygonWithFaces[j];
            bool face1 = el.second.first == faces.first || (ind != 0 && el.second.first == faces.second);
            bool face2 = el.second.second == faces.first || (ind != 0 && el.second.second == faces.second);
            if((std::find(visited.begin(), visited.end(), j) == visited.end()
                || j == 0 && i == unorderedPolygonWithFaces.size() - 1) && (face1 || face2)) {
                int face = face1 ? el.second.first : el.second.second;
                std::vector<size_t> internal_points = processFace(octree, edges, node, f, sides[face], points[polygon.back()], points);
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
    ImplicitFunction sphere = [](Vector p) { return sqrt(p.x()*p.x() + p.y()*p.y() + p.z()*p.z()) - 0.5; };
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

    ImplicitFunction& func = sphere;
    Point_set points; // This is only here because Orthtree constructor requires it
    Octree octree(Octree::Bbox(-1.1,-1.1,-1.1,1.1,1.1,1.1), points);
    octree.refine(Split_by_closeness(func, octree));

    Edges edges(octree, func);

    // Add edges to container
    for (Octree::Node node : octree.traverse<Preorder_traversal>()) {
        if(node == octree.root())
            edges.add_root(node);
        else
            edges.add_node(node);
    }

    std::vector<Point> vertices;
    std::vector<std::vector<size_t>> faces;
    // Traverse octree and process each cell
    // todo: later prepare for parallelization
    for (Octree::Node node : octree.traverse<Preorder_traversal>()) {
        //std::cout << node << std::endl;
        //std::cout << octree.bbox(node) << std::endl;
        if(node.is_leaf()) {
            std::vector<size_t> polygon = processNode(octree, edges, node, func, vertices);
            if(polygon.size() > 0)
                faces.push_back(polygon);
        }
    }
    //std::cout << octree.depth() << std::endl;
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