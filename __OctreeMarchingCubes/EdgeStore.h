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

    typedef CGAL::Vector_3<Traits_> Vector;

    typedef typename Traits_::FT FT;

    std::map<typename Node::Global_coordinates, std::array<Edge*, 12>> edgeLists;

    const Octree<Traits_, PointRange_, PointMap_>& tree;

    std::function<FT(Vector)> func;

public:

    EdgeStore(const Octree<Traits_, PointRange_, PointMap_>& tree, std::function<FT(Vector)> func)
    : tree(tree), func(func) {}

    const std::array<Edge*, 12>& get(const Node& n) const { return edgeLists.at(n.global_coordinates()); }

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
        edgeLists[root.global_coordinates()] = edges;
    }

    void addNode(const Node& node) {
        Node parent = node.parent();
        auto coords = node.local_coordinates();
        size_t corner = coords[0]*4 + coords[1]*2 + coords[2]*1;

        std::array<Edge*, 12> parentEdges = edgeLists[parent.global_coordinates()];

        int ind = 0;
        Bbox b = tree.bbox(node);
        std::array<Edge*, 12> edges;
        for (int i = 0; i < 8; ++i) {
            for (int k = 0; k <= 2; ++k) {
                int j = i | (1 << k);
                if (i != j){
                    if (i == corner) {
                        if (parentEdges[ind]->getChild1() == nullptr) {
                            Vector p (j&4 ? b.xmax() : b.xmin(), j&2 ? b.ymax() : b.ymin(), j&1 ? b.zmax() : b.zmin());
                            parentEdges[ind]->divide(node, parent[j], func(p));
                        }
                        edges[ind++] = parentEdges[ind]->getChild1();
                    }
                    else if (j == corner) {
                        if (parentEdges[ind]->getChild2() == nullptr) {
                            Vector p (i&4 ? b.xmax() : b.xmin(), i&2 ? b.ymax() : b.ymin(), i&1 ? b.zmax() : b.zmin());
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
        edgeLists[node.global_coordinates()] = edges;
    }

};

}

#endif