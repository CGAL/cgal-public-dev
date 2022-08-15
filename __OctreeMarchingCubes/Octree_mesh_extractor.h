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

#ifndef CGAL_OCTREEMESHEXTRACTOR_H
#define CGAL_OCTREEMESHEXTRACTOR_H

#include <CGAL/Octree.h>

#include "Edge_store.h"

namespace CGAL {

template <typename GeomTraits>
class Octree_mesh_extractor {
public:
    typedef CGAL::Point_3<GeomTraits> Point;
    typedef CGAL::Vector_3<GeomTraits> Vector;
    typedef typename GeomTraits::FT FT;
    typedef CGAL::Point_set_3<Point> Point_set;
    typedef typename Point_set::Point_map Point_map;

    typedef CGAL::Octree<GeomTraits, Point_set, Point_map> Octree;
    typedef CGAL::Orthtrees::Preorder_traversal Preorder_traversal;
    typedef typename CGAL::Octree<GeomTraits, Point_set, Point_map>::Bbox Bbox;
    typedef typename CGAL::Octree<GeomTraits, Point_set, Point_map>::Node Node;
    typedef CGAL::Octree_edge<GeomTraits> Edge;

    typedef typename CGAL::Orthtree_traits_3<GeomTraits>::Adjacency Adjacency;

    typedef std::function<FT(Vector)> Implicit_function;

    typedef CGAL::Edge_store<GeomTraits, Point_set, Point_map> Edge_store;

    Octree_mesh_extractor(const Octree& tree, Implicit_function func)
    : octree(tree), func(func), edges(tree, func) {
        // Add edges to container
        for (Node node : octree.template traverse<Preorder_traversal>()) { // ????
            if(node == octree.root())
                edges.add_root(node);
            else
                edges.add_node(node);
        }
    }

    // inserts point to mesh, if not there already
    // and returns its index
    size_t add_point_to_mesh(const Edge* e) {
        auto i = edge_points_in_mesh.find(*e);
        if (i != edge_points_in_mesh.end())
            return i->second;

        Point p = e->extract_isovertex();
        size_t ind = vertices.size();
        edge_points_in_mesh[*e] = ind;
        vertices.push_back(p);
        return ind;
    }

    // this function is to find closed polygons in the inside of
    // cell faces; currently it is not called in any tests
    std::vector<const Edge*> find_closed_polyline(
        std::vector<std::vector<const Edge*>>& segments
    ) {
        std::vector<const Edge*> polyline;

        const Edge* start = segments[0][0];
        const Edge* next = segments[0][1];

        polyline.push_back(start);
        std::cerr << *(polyline.back());
        int v = 0, p = 1;
        std::set<int> vs;
        vs.insert(v);
        while (next != start) {
            polyline.push_back(next);
            std::cerr << *(polyline.back());
            bool found = false;
            for (int j = 0; j < segments.size(); ++j) if (v != j) {
                if (segments[j][0] == polyline.back()) { found = true; v = j; p = 1; break; }
                else if (segments[j][1] == polyline.back()) { found = true; v = j; p = 0; break; }
            }
            if (found) {
                next = segments[v][p];
                std::cerr << "next (segment found):" << *next;
            }
            if (!found) {
                next = next->twin_edge();
                std::cerr << "next (twin edge):" << *next;
                v = -1;
                for (int j = 0; j < segments.size(); ++j) {
                    if (segments[j][0] == next || segments[j][1] == next) { v = j; }
                }
                if (v == -1)
                throw Error_exception("Octree_marching_cubes"
                    , "twin edge not in any segment"
                    , "Octree_mesh_extractor.h", 97);
            }
            vs.insert(v);
        }

        std::vector<std::vector<const Edge*>> new_segments;
        for (int i = 0; i < segments.size(); ++i)
            if (vs.find(i) == vs.end())
                new_segments.push_back(segments[i]);
        segments = std::move(new_segments);

        return polyline;
    }

    // returns the internal points of the finer polyline
    // from the other side of the cell face
    std::vector<size_t> process_face(
        const Node& node,
        Adjacency adj,
        size_t start,
        std::vector<std::vector<size_t>>& additional_polygons) {

        Node neighbour = node.adjacent_node(adj);

        if (neighbour.is_null() || neighbour.is_leaf())
            return std::vector<size_t>{};

        unsigned mask = (!(adj & ~1) ? 1 : adj & ~1);
        std::queue<Node> to_divide, leaf_neighbours;
        to_divide.push(neighbour);
        while (!to_divide.empty()) {
            Node n = to_divide.front();
            for (unsigned i = 0; i < 8; ++i) {
                if ( !(adj & 1) != !(mask & i) ) { // if the node's child with index i is along the given side of the cell
                    if (n[i].is_leaf())
                        leaf_neighbours.push(n[i]);
                    else
                        to_divide.push(n[i]);
                }
            }
            to_divide.pop();
        }

        std::vector<std::vector<const Edge*>> segments;

        int ind = 0, startv = -1, startp = -1;
        while (!leaf_neighbours.empty()) {
            Node child = leaf_neighbours.front();
            std::array<Edge*, 12> child_edges = edges.get(child);
            std::vector<const Edge*> segment;
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
                        segment.push_back(edge);
                        if (add_point_to_mesh(segment.back()) == start) {
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
            leaf_neighbours.pop();
        }

        std::vector<const Edge*> polyline_edges;

        int v = startv, p = !startp;
        std::set<int> vs;
        for (int i = 0; i < segments.size(); ++i) {
            polyline_edges.push_back(segments[v][p]);
            vs.insert(v);
            bool found = false;
            // try to find a segment continuing the polyline
            for (int j = 0; j < segments.size(); ++j) if (v != j) {
                if (segments[j][0] == polyline_edges.back()) { found = true; v = j; p = 1; break; }
                else if (segments[j][1] == polyline_edges.back()) { found = true; v = j; p = 0; break; }
            }
            if (!found) {
                // if none, try connecting to a twin edge
                auto next = polyline_edges.back()->twin_edge();
                if (next == nullptr) break; // no twin edge: polyline to be ended
                v = -1;
                for (int j = 0; j < segments.size(); ++j) {
                    if (segments[j][0] == next || segments[j][1] == next) { v = j; }
                }
                if (v == -1) break; // twin edge exists, but not in this cell: polyline is ended
            }
        }
        std::vector<size_t> polyline(polyline_edges.size());
        std::transform(polyline_edges.cbegin(), polyline_edges.cend(), polyline.begin(),
            [this] (const Edge* e) { return add_point_to_mesh(e); });

        // unsure, if this part is needed, currently it does not run with any of the tests
        if(polyline.size() < segments.size()) {
            std::vector<std::vector<const Edge*>> remaining_segments;
            for (int i = 0; i < segments.size(); ++i) {
                if (vs.find(i) == vs.end()) remaining_segments.push_back(segments[i]);
            }
            std::cout << "Not all points were used." << std::endl
            << "Remaining segments: " << remaining_segments.size() << std::endl;

            std::cerr << "Edges already used: " << std::endl;
            for(auto it : polyline_edges)
                std::cerr << *it << std::endl;

            std::cerr << "Node: " << node << std::endl;
            while (remaining_segments.size() > 0) {
                auto poly = find_closed_polyline(remaining_segments);
                std::vector<size_t> polygon(poly.size());
                std::transform(poly.cbegin(), poly.cend(), polygon.begin(),
                    [this] (const Edge* e) { return add_point_to_mesh(e); });
                additional_polygons.push_back(polygon);
            }
        }
        if(polyline.size() > 0)
            polyline.pop_back();

        return polyline;

    }

    // returns the vertices extracted from a given cell, ordered as a polyline
    void process_node(
        const Node& node) {

        Bbox b = octree.bbox(node);

        std::vector<std::pair<size_t, std::pair<int,int>>> unordered_polygon_with_faces;

        for(auto edge : edges.get(node)) {
            std::pair<Point,Point> seg = edge->segment();
            Vector p1 (seg.first.x(), seg.first.y(), seg.first.z());
            Vector p2 (seg.second.x(), seg.second.y(), seg.second.z());

            if(edge->values().first * edge->values().second < 0) {
                auto minEdge = edge->find_minimal_edge();
                std::pair<Point,Point> minSeg = minEdge->segment();
                Vector p1 (minSeg.first.x(), minSeg.first.y(), minSeg.first.z());
                Vector p2 (minSeg.second.x(), minSeg.second.y(), minSeg.second.z());
                auto vals = minEdge->values();
                auto corners = edge->corners(node);
                if (vals.second < 0) { auto tmp = corners.first; corners.first = corners.second; corners.second = tmp;  }
                unordered_polygon_with_faces.push_back(std::make_pair(
                    add_point_to_mesh(minEdge),
                    corners_to_faces(corners)));
            }
        }
        if(unordered_polygon_with_faces.size() == 0) return;

        std::vector<size_t> polygon;
        std::vector<std::vector<size_t>> additional_polygons;

        int ind = 0;
        std::vector<int> visited {ind};
        for(int i = 0; i < unordered_polygon_with_faces.size(); ++i) {
            polygon.push_back(unordered_polygon_with_faces[ind].first);
            std::pair<int,int> faces = unordered_polygon_with_faces[ind].second;
            for(int j = 0; j < unordered_polygon_with_faces.size(); ++j) {
                auto& el = unordered_polygon_with_faces[j];
                bool face1 = el.second.first == faces.first || (ind != 0 && el.second.first == faces.second); // at the first vertex we want to start in the good direction
                bool face2 = el.second.second == faces.first || (ind != 0 && el.second.second == faces.second);
                if((std::find(visited.begin(), visited.end(), j) == visited.end()
                    || j == 0 && i == unordered_polygon_with_faces.size() - 1) // we want to check the face between the last and first vertex
                    && (face1 || face2)) {
                    int face = face1 ? el.second.first : el.second.second;
                    std::vector<size_t> internal_points = process_face(
                        node, Adjacency(face), polygon.back(),
                        additional_polygons);
                    for (auto it : internal_points) {
                        polygon.push_back(it);
                    }
                    ind = j;
                    visited.push_back(ind);
                    break;
                }
            }
        }
        if(polygon.size() > 0)
            faces.push_back(polygon);
    }

    std::vector<Point> get_vertices() { return vertices; }
    std::vector<std::vector<size_t>> get_faces() { return faces; }

private:

    const Octree& octree;
    Implicit_function func;
    Edge_store edges;

    std::map<Edge, size_t> edge_points_in_mesh;

    std::vector<Point> vertices;
    std::vector<std::vector<size_t>> faces;

    const int corners_to_faces_table[8][3][2] = {
        {{2,0}, {0,4}, {4,2}},
        {{0,2}, {5,0}, {2,5}},
        {{0,3}, {4,0}, {3,4}},
        {{3,0}, {0,5}, {5,3}},
        {{1,2}, {4,1}, {2,4}},
        {{2,1}, {1,5}, {5,2}},
        {{3,1}, {1,4}, {4,3}},
        {{1,3}, {5,1}, {3,5}},
    };

    // returns the two faces along the segment between the two corners
    // such that from the direction of the edge, the left face is first
    std::pair<int,int> corners_to_faces(std::pair<int,int> corners) {
        int p0 = corners.first, p1 = corners.second;
        int r0 = std::abs(p0 - p1) / 2; // 1,2,4 -> 0,1,2
        return std::make_pair(corners_to_faces_table[p0][r0][0], corners_to_faces_table[p0][r0][1]);
    }

};

}

#endif