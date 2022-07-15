#ifndef CGAL_EDGESTORE_H
#define CGAL_EDGESTORE_H

#include <CGAL/Orthtree/Edge.h>
#include <CGAL/Orthtree/Node.h>
#include <CGAL/Octree.h>

#include <map>

namespace CGAL {

template<typename Traits_, typename PointRange_,
         typename PointMap_ = Identity_property_map<typename Traits_::Point_d> >
class EdgeStore {
    typedef typename Octree<Traits_, PointRange_, PointMap_>::Node Node;
    typedef typename Octree<Traits_, PointRange_, PointMap_>::Edge Edge;
    typedef typename Octree<Traits_, PointRange_, PointMap_>::Bbox Bbox;
    typedef typename Octree<Traits_, PointRange_, PointMap_>::Point Point;

    typedef CGAL::Vector_3<Traits_> Vector;

    typedef typename Traits_::FT FT;

    std::map<Node, std::array<Edge*, 12>> edgeLists;

    const Octree<Traits_, PointRange_, PointMap_>& tree;

    std::function<FT(Vector)> func;

public:

    EdgeStore(const Octree<Traits_, PointRange_, PointMap_>& tree, std::function<FT(Vector)> func)
    : tree(tree), func(func) {}

    const std::array<Edge*, 12>& get(const Node& n) const { return edgeLists.at(n); }

    void addRoot(const Node& root) {
        int ind = 0;
        Bbox b = tree.bbox(root);
        std::array<Edge*, 12> edges;
        for (int i = 0; i < 8; ++i) {
            for (int k = 0; k <= 2; ++k) {
                int j = i | (1 << k);
                if (i != j){
                    Vector p1 (i&4 ? b.xmax() : b.xmin(), i&2 ? b.ymax() : b.ymin(), i&1 ? b.zmax() : b.zmin());
                    Vector p2 (j&4 ? b.xmax() : b.xmin(), j&2 ? b.ymax() : b.ymin(), j&1 ? b.zmax() : b.zmin());
                    edges[ind++] = new Edge(root, i, j, func(p1), func(p2));
                }
            }
        }
        edgeLists[root] = edges;
    }

    void addNode(const Node& node) {
        Node parent = node.parent();
        auto coords = node.local_coordinates();
        //size_t corner = coords[0]*4 + coords[1]*2 + coords[2]*1;
        size_t corner = coords[0]*1 + coords[1]*2 + coords[2]*4; // this does not feel right, seems contrary to local_coordinates doc

        std::array<Edge*, 12> parentEdges = edgeLists[parent];

        int ind = 0;
        Bbox b = tree.bbox(node);
        std::array<Edge*, 12> edges;
        for (int i = 0; i < 8; ++i) {
            for (int k = 0; k <= 2; ++k) {
                int j = i | (1 << k);
                if (i != j){
                    if (i == corner) {
                        if (parentEdges[ind]->getChild1() == nullptr) {
                            std::pair<Point,Point> seg = tree.segment(*(parentEdges[ind]));
                            Vector p1 (seg.first.x(), seg.first.y(), seg.first.z());
                            Vector p2 (seg.second.x(), seg.second.y(), seg.second.z());
                            Vector p ((p1 + p2) / 2);
                            parentEdges[ind]->divide(node, parent[j], func(p));
                        }
                        edges[ind++] = parentEdges[ind]->getChild1();
                    }
                    else if (j == corner) {
                        if (parentEdges[ind]->getChild2() == nullptr) {
                            std::pair<Point,Point> seg = tree.segment(*(parentEdges[ind]));
                            Vector p1 (seg.first.x(), seg.first.y(), seg.first.z());
                            Vector p2 (seg.second.x(), seg.second.y(), seg.second.z());
                            Vector p ((p1 + p2) / 2);
                            parentEdges[ind]->divide(parent[i], node, func(p));
                        }
                        edges[ind++] = parentEdges[ind]->getChild2();
                    }
                    else {
                        Vector p1 (i&4 ? b.xmax() : b.xmin(), i&2 ? b.ymax() : b.ymin(), i&1 ? b.zmax() : b.zmin());
                        Vector p2 (j&4 ? b.xmax() : b.xmin(), j&2 ? b.ymax() : b.ymin(), j&1 ? b.zmax() : b.zmin());
                        edges[ind++] = new Edge(node, i, j, func(p1), func(p2));
                    }
                }
            }
        }
        edgeLists[node] = edges;
    }

};

}

#endif