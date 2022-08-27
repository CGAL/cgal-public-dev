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

#ifndef CGAL_ORTHTREE_EDGE_H
#define CGAL_ORTHTREE_EDGE_H

#include <CGAL/Orthtree.h>
#include <CGAL/exceptions.h>

#include "../Isosurfacing_3/Octree_wrapper.h"

#include <optional>

namespace CGAL {

/*!

  \brief represents an edge of an Octree.
 */
template <typename Traits_>
class Octree_edge
{
    typedef typename Traits_::FT FT;
    typedef CGAL::Vector_3<Traits_> Vector_3;
    typedef CGAL::Point_3<Traits_> Point_3;
    typedef CGAL::Bbox_3 Bbox;
    typedef Octree_wrapper<Traits_> Tree;
    typedef typename Tree::Node Node;

public:
    // constructs an edge without a parent (a.k.a. a root edge)
    Octree_edge(const Tree& tree, const std::array<std::uint32_t, 6>& coords, int depth, FT v1, FT v2)
    : tree(tree)
    , Global_coordinates(coords), depth_(depth)
    , value1(v1), value2(v2)
    , parent(nullptr), child1(nullptr), child2(nullptr) {}

    // constructs an edge with a specified parent
    Octree_edge(const Tree& tree, const std::array<std::uint32_t, 6>& coords, int depth, Octree_edge* parent, FT v1, FT v2)
    : tree(tree)
    , Global_coordinates(coords), depth_(depth)
    , value1(v1), value2(v2)
    , parent(parent), child1(nullptr), child2(nullptr) {}

    ~Octree_edge() {
        delete child1;
        delete child2;
    }

    // creates two children for the edge and inserts the given value to the midpoint
    void divide(FT midValue) {
        std::array<std::uint32_t, 3> midCoords { (Global_coordinates[0] + Global_coordinates[3]),
                    (Global_coordinates[1] + Global_coordinates[4]),
                    (Global_coordinates[2] + Global_coordinates[5]) };
        child1 = new Octree_edge(
            tree,
            {2*Global_coordinates[0], 2*Global_coordinates[1], 2*Global_coordinates[2], midCoords[0], midCoords[1], midCoords[2]}
            , depth_ + 1, this, value1, midValue);
        child2 = new Octree_edge(
            tree,
            {midCoords[0], midCoords[1], midCoords[2], 2*Global_coordinates[3], 2*Global_coordinates[4], 2*Global_coordinates[5]}
            , depth_ + 1, this, midValue, value2);
    }

    // [Kazhdan et al, 2007] 3.2 "Defining Consistent Isovertices"
    const Octree_edge* find_minimal_edge(FT isovalue) const {
        if(child1 == nullptr) return this;
        else if((child1->value1 - isovalue) * (child1->value2 - isovalue) <= 0) return child1->find_minimal_edge(isovalue);
        else return child2->find_minimal_edge(isovalue);
    }

    // [Kazhdan et al, 2007] 3.2 "Closing the Isopolylines"
    const Octree_edge* twin_edge(FT isovalue) const {
        if (parent == nullptr) return nullptr;
        if ((parent->value1 - isovalue) * (parent->value2 - isovalue) <= 0) return parent->twin_edge(isovalue);
        else if (parent->child1 == this) return parent->child2->find_minimal_edge(isovalue);
        else return parent->child1->find_minimal_edge(isovalue);
    }

    Octree_edge* get_child1() const { return child1; }
    Octree_edge* get_child2() const { return child2; }

    // the two corner points of the given node that are on this edge
    // assumes the edge is an edge of that node, otherwise gives nonsense results
    std::pair<int,int> corners(Node n) const {
        int corner1 = (Global_coordinates[0] > n.global_coordinates()[0])*4
            + (Global_coordinates[1] > n.global_coordinates()[1])*2
            + (Global_coordinates[2] > n.global_coordinates()[2])*1;
        int corner2 = (Global_coordinates[3] > n.global_coordinates()[0])*4
            + (Global_coordinates[4] > n.global_coordinates()[1])*2
            + (Global_coordinates[5] > n.global_coordinates()[2])*1;
        return std::make_pair(corner1, corner2);
    }
    std::pair<FT,FT> values() const { return std::make_pair(value1, value2); }
    std::array<std::uint32_t, 6> global_coordinates() const { return Global_coordinates; }
    int depth() const { return depth_; }
    bool is_root() const { return parent == nullptr; }

    friend std::ostream& operator<<(std::ostream& out, const Octree_edge& edge) {
        auto c = edge.global_coordinates();
        out << "{ ( " << c[0] << " " << c[1] << " " << c[2] << " ) ( " << c[3] << " " << c[4] << " " << c[5] << " ) " << edge.depth() << " }" << std::endl;
        return out;
    }

    std::pair<Point_3, Point_3> segment() const {
        Bbox box = tree.bbox(tree.root());
        std::array<FT, 3> p1;
        std::array<FT, 3> p2;
        for (int i = 0; i < 3; i++) {
            FT length = (box.max(i) - box.min(i)) / (1 << depth());
            p1[i] = box.min(i) + (global_coordinates()[i] * length);
            p2[i] = box.min(i) + (global_coordinates()[i + 3] * length);
        }
        return std::make_pair(Point_3(p1[0], p1[1], p1[2]), Point_3(p2[0], p2[1], p2[2]));
    }

    Point_3 extract_isovertex(FT isovalue) const {
        if (isovertex) { // this apparently never runs
            return *isovertex;
        }
        else {
            Point_3 p1, p2;
            Vector_3 r1, r2;
            FT d1, d2;
            std::tie(p1, p2) = segment();
            r1 = Vector_3(p1.x(), p1.y(), p1.z());
            r2 = Vector_3(p2.x(), p2.y(), p2.z());
            std::tie(d1, d2) = values();
            FT mu = -1.f;

            // don't divide by 0
            if (abs(d2 - d1) < std::numeric_limits<FT>::epsilon()) {
                mu = 0.5f;  // if both points have the same value, assume isolevel is in the middle
            } else {
                mu = (isovalue - d1) / (d2 - d1);
            }

            if (mu < 0.f || mu > 1.f) {
                throw Error_exception("Octree_marching_cubes"
                    , "isolevel is not between points"
                    , "Octree_edge.h", 137);
            }

            Vector_3 res = r2 * mu + r1 * (1 - mu);
            isovertex = std::make_optional(Point_3(res.x(), res.y(), res.z()));
            return *isovertex;
        }
    }

private:
    const Tree &tree;

    std::array<std::uint32_t, 6> Global_coordinates; // 4 would be enough
    int depth_;
    FT value1;
    FT value2;

    mutable std::optional<Point_3> isovertex; // maybe cache index instead? would need setting from outside

    Octree_edge* parent;
    Octree_edge* child1;
    Octree_edge* child2;
};

// needed because it is a key in an std::map, maybe a hash map would be better
template <typename Traits_>
bool operator<(const Octree_edge<Traits_>& lhs, const Octree_edge<Traits_>& rhs) {
    if (lhs.depth() < rhs.depth()) return true;
    else if (lhs.depth() > rhs.depth()) return false;

    return lhs.global_coordinates() < rhs.global_coordinates();
}

}

#endif // CGAL_ORTHTREE_EDGE_H