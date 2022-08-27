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

#include "Octree_edge.h"

#include <map>

template <typename Node>
struct Node_hash
{
    std::size_t operator()(const Node& n) const {
        return n.hash_value();
    }
};


namespace CGAL {

// Stores all the edge-trees for an octree
template<typename Traits_>
class Edge_store {
public:

    typedef Octree_wrapper<Traits_> Octree;
    typedef typename Octree::Node Node;
    typedef Octree_edge<Traits_> Edge;
    typedef typename CGAL::Bbox_3 Bbox;

    typedef CGAL::Orthtree_traits_3<Traits_> Traits3;

    typedef CGAL::Point_3<Traits_> Point;
    typedef CGAL::Vector_3<Traits_> Vector;

    typedef typename Traits_::FT FT;

    Edge_store(const Octree& tree)
    : tree(tree) {}

    ~Edge_store() {
        for (auto it : roots) {
            delete it;
        }
    }

    const std::array<Edge*, 12>& get(const Node& n) const { return edge_lists.at(n); }

    // adds edges of the root node to the edge trees and registers them to the node
    void add_root(const Node& root) {
        int ind = 0;
        Bbox b = tree.bbox(root);
        std::array<Edge*, 12> edges;
        for (unsigned i = 0; i < 8; ++i) {
            for (unsigned k = 0; k <= 2; ++k) {
                unsigned j = i | (1 << k);
                if (i != j){
                    std::bitset<3> c1(i);
                    bool tmp = c1[0]; c1[0] = c1[2]; c1[2] = tmp; // coordinate order has to be reversed
                    std::bitset<3> c2(j);
                    tmp = c2[0]; c2[0] = c2[2]; c2[2] = tmp;
                    FT v1 = tree.value(tree.vertex_uniform_coordinates(root, c1));
                    FT v2 = tree.value(tree.vertex_uniform_coordinates(root, c2));
                    Edge* edge = new Edge(tree,
                        {(i&4)>>2, (i&2)>>1, i&1, (j&4)>>2, (j&2)>>1, j&1}, 0,
                        v1, v2);
                    edges[ind++] = edge;
                    roots.push_back(edge);
                }
            }
        }
        edge_lists[root] = edges;
    }

    // adds edges of a non-root node to the edge trees and registers them to the node
    void add_node(const Node& node) {
        Node parent = node.parent();
        auto coords1 = node.global_coordinates();
        auto coords = node.local_coordinates();
        size_t corner = coords[0]*4 + coords[1]*2 + coords[2]*1;

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
                            std::bitset<3> c(j);
                            bool tmp = c[0]; c[0] = c[2]; c[2] = tmp;
                            FT v = tree.value(tree.vertex_uniform_coordinates(node, c));
                            parent_edges[ind]->divide(v);
                        }
                        edges[ind++] = parent_edges[ind]->get_child1();
                    }
                    else if (j == corner) { // edge is the right child of another edge
                        if (parent_edges[ind]->get_child2() == nullptr) {
                            std::bitset<3> c(i);
                            bool tmp = c[0]; c[0] = c[2]; c[2] = tmp;
                            FT v = tree.value(tree.vertex_uniform_coordinates(node, c));
                            parent_edges[ind]->divide(v);
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
                        // check if edge is already registered for a neighbouring cell
                        for(auto a : adjacencies) {
                            Node neighbour = node.adjacent_node(a);
                            auto edges = edge_lists.find(neighbour);
                            if(edges != edge_lists.end())
                                for(auto edge : edges->second) {
                                    if(edge->global_coordinates() == coords)
                                        sameEdge = edge;
                                }
                        }
                        std::bitset<3> c1(i);
                        bool tmp = c1[0]; c1[0] = c1[2]; c1[2] = tmp;
                        std::bitset<3> c2(j);
                        tmp = c2[0]; c2[0] = c2[2]; c2[2] = tmp;
                        FT v1 = tree.value(tree.vertex_uniform_coordinates(node, c1));
                        FT v2 = tree.value(tree.vertex_uniform_coordinates(node, c2));
                        edges[ind++] = sameEdge == nullptr
                            ? (roots.push_back(new Edge(tree, coords, node.depth(), v1, v2)), roots.back())
                            : sameEdge;
                    }
                }
            }
        }
        edge_lists[node] = edges;
    }

private:

    std::unordered_map<Node, std::array<Edge*, 12>, Node_hash<Node>> edge_lists;

    const Octree& tree;

    std::vector<Edge*> roots; // for destructing

};

}

#endif