// Copyright (c) 2022  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ágoston Sipos

#ifndef CGAL_EDGESTORE_H
#define CGAL_EDGESTORE_H

#include <CGAL/Orthtree/Node.h>
#include <CGAL/Octree.h>

#include "Edge.h"

#include <map>

namespace CGAL {

template<typename Traits_, typename PointRange_,
         typename PointMap_ = Identity_property_map<typename Traits_::Point_3> >
class Edge_store {
public:

    typedef typename Octree<Traits_, PointRange_, PointMap_>::Node Node;
    typedef Octree_edge<Traits_> Edge;
    typedef typename Octree<Traits_, PointRange_, PointMap_>::Bbox Bbox;
    typedef typename Octree<Traits_, PointRange_, PointMap_>::Point Point;

    typedef CGAL::Orthtree_traits_3<Traits_> Traits3;

    typedef CGAL::Vector_3<Traits_> Vector;

    typedef typename Traits_::FT FT;

    Edge_store(const Octree<Traits_, PointRange_, PointMap_>& tree, std::function<FT(Vector)> func)
    : tree(tree), func(func) {}

    const std::array<Edge*, 12>& get(const Node& n) const { return edge_lists.at(n); }

    void add_root(const Node& root) {
        int ind = 0;
        Bbox b = tree.bbox(root);
        std::array<Edge*, 12> edges;
        for (unsigned i = 0; i < 8; ++i) {
            for (unsigned k = 0; k <= 2; ++k) {
                unsigned j = i | (1 << k);
                if (i != j){
                    Vector p1 (i&4 ? b.xmax() : b.xmin(), i&2 ? b.ymax() : b.ymin(), i&1 ? b.zmax() : b.zmin());
                    Vector p2 (j&4 ? b.xmax() : b.xmin(), j&2 ? b.ymax() : b.ymin(), j&1 ? b.zmax() : b.zmin());
                    edges[ind++] = new Edge(tree, {(i&4)>>2, (i&2)>>1, i&1, (j&4)>>2, (j&2)>>1, j&1}, 0, func(p1), func(p2));
                }
            }
        }
        edge_lists[root] = edges;
    }

    void add_node(const Node& node) {
        Node parent = node.parent();
        auto coords1 = node.global_coordinates();
        auto coords = node.local_coordinates();
        size_t corner = coords[0]*4 + coords[1]*2 + coords[2]*1;
        size_t corner2 = coords[0]*1 + coords[1]*2 + coords[2]*4; // node child indexing uses this sceme for some reason...

        std::array<Edge*, 12> parent_edges = edge_lists[parent];

        int ind = 0;
        Bbox b = tree.bbox(node);
        std::array<Edge*, 12> edges;
        for (unsigned i = 0; i < 8; ++i) {
            for (unsigned k = 0; k <= 2; ++k) {
                unsigned j = i | (1 << k);
                if (i != j){
                    if (i == corner) { // edge is the left child of another edge
                        if (parent_edges[ind]->get_child1() == nullptr) {
                            std::pair<Point,Point> seg = parent_edges[ind]->segment();
                            Vector p1 (seg.first.x(), seg.first.y(), seg.first.z());
                            Vector p2 (seg.second.x(), seg.second.y(), seg.second.z());
                            Vector p ((p1 + p2) / 2);
                            parent_edges[ind]->divide(func(p));
                        }
                        edges[ind++] = parent_edges[ind]->get_child1();
                    }
                    else if (j == corner) { // edge is the right child of another edge
                        if (parent_edges[ind]->get_child2() == nullptr) {
                            std::pair<Point,Point> seg = parent_edges[ind]->segment();
                            Vector p1 (seg.first.x(), seg.first.y(), seg.first.z());
                            Vector p2 (seg.second.x(), seg.second.y(), seg.second.z());
                            Vector p ((p1 + p2) / 2);
                            parent_edges[ind]->divide(func(p));
                        }
                        edges[ind++] = parent_edges[ind]->get_child2();
                    }
                    else { // edge is not the child of any other edge
                        Edge* sameEdge = nullptr;
                        auto c = node.global_coordinates();
                        std::array<std::uint32_t, 6> coords {c[0] + ((i&4)>>2), c[1] + ((i&2)>>1), c[2] + (i&1)
                                                , c[0] + ((j&4)>>2), c[1] + ((j&2)>>1), c[2] + (j&1)};
                        std::array<typename Traits3::Adjacency, 6> adjacencies
                            {Traits3::LEFT,Traits3::RIGHT,Traits3::DOWN,Traits3::UP,Traits3::BACK,Traits3::FRONT};
                        for(auto a : adjacencies) {
                            Node neighbour = node.adjacent_node(a);
                            auto edges = edge_lists.find(neighbour);
                            if(edges != edge_lists.end())
                                for(auto edge : edges->second) {
                                    if(edge->global_coordinates() == coords)
                                        sameEdge = edge;
                                }
                        }
                        Vector p1 (i&4 ? b.xmax() : b.xmin(), i&2 ? b.ymax() : b.ymin(), i&1 ? b.zmax() : b.zmin());
                        Vector p2 (j&4 ? b.xmax() : b.xmin(), j&2 ? b.ymax() : b.ymin(), j&1 ? b.zmax() : b.zmin());
                        edges[ind++] = sameEdge == nullptr
                            ? new Edge(tree, coords, node.depth(), func(p1), func(p2))
                            : sameEdge;
                    }
                }
            }
        }
        edge_lists[node] = edges;
    }

private:

    std::map<Node, std::array<Edge*, 12>> edge_lists;

    const Octree<Traits_, PointRange_, PointMap_>& tree;

    std::function<FT(Vector)> func;

};

}

#endif