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

// Builds up the edge trees for the octree, then
// extracts the polygons node-by-node
template <typename GeomTraits>
class Octree_mesh_extractor {
public:
    typedef CGAL::Point_3<GeomTraits> Point;
    typedef CGAL::Vector_3<GeomTraits> Vector;
    typedef typename GeomTraits::FT FT;

    typedef Octree_wrapper<GeomTraits> Octree;
    typedef CGAL::Orthtrees::Preorder_traversal Preorder_traversal;
    typedef typename CGAL::Bbox_3 Bbox;
    typedef typename Octree::Node Node;
    typedef CGAL::Octree_edge<GeomTraits> Edge;

    typedef typename CGAL::Orthtree_traits_3<GeomTraits>::Adjacency Adjacency;

    typedef std::function<FT(Vector)> Implicit_function;

    typedef CGAL::Edge_store<GeomTraits> Edge_store;

    Octree_mesh_extractor(const Octree& tree, FT isovalue)
    : octree(tree), edges(tree), isovalue(isovalue) {
        // Add edges to container
        for (Node node : octree.template traverse<Preorder_traversal>()) {
            if(node == octree.root())
                edges.add_root(node);
            else
                edges.add_node(node);
        }
    }

    // inserts point to mesh, if not there already
    // and returns its index
    size_t add_point_to_mesh(const Edge* e) {
        auto i = edge_points_in_mesh.find(e);
        if (i != edge_points_in_mesh.end())
            return i->second;

        Point p = e->extract_isovertex(isovalue);
        size_t ind = vertices.size();
        edge_points_in_mesh[e] = ind;
        vertices.push_back(p);
        return ind;
    }

    // this function is to find closed polygons in the inside of cell faces
    // it returns a closed loop of isovertices and removes those segments from the parameter
    std::vector<const Edge*> find_closed_polyline(std::vector<std::vector<const Edge*>>& segments) {
        std::vector<const Edge*> polyline;

        const Edge* start = segments[0][0];
        const Edge* next = segments[0][1];

        polyline.push_back(start);
        int v = 0, p = 1;
        std::set<int> vs; // segments flagged for removal
        vs.insert(v);
        while (next != start) {
            polyline.push_back(next);
            bool found = false;
            for (int j = 0; j < segments.size(); ++j) if (v != j) {
                if (segments[j][0] == polyline.back()) { found = true; v = j; p = 1; break; }
                else if (segments[j][1] == polyline.back()) { found = true; v = j; p = 0; break; }
            }
            if (found) {
                next = segments[v][p];
            }
            if (!found) { // if no stored segment continues the line, we move to twin edge
                next = next->twin_edge(isovalue);
                v = -1;
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

    // polygons from find_closed_polyline() are not guaranteed to have a good ordering
    // however, the outside is always towards the bigger cell so we can reorient based on its center point
    void reorient_polygon(std::vector<size_t>& polygon, const Point& center) {
        size_t n = polygon.size();
        int x = 0, y = 0;
        for (size_t i = 0; i < n; ++i) {
            Point a = vertices[polygon[i]];
            Point b = vertices[polygon[(i+1)%n]];
            Point c = vertices[polygon[(i+2)%n]];
            FT f = CGAL::cross_product(b - a, c - b) * (center - b);
            if (f > 0) ++x;
            else ++y;
        }
        if (x < y) {
            for (size_t i = 0; i < n / 2; ++i) {
                std::swap(polygon[i], polygon[n-i-1]);
            }
        }
    }

    // returns the internal points of the finer polyline
    // from the other side of the cell face
    std::vector<size_t> process_face(
        const Node& node, // the big cell which might have subdivided neighbours
        Adjacency adj, // which face of it we are working on
        std::optional<size_t> start, // if there is an edge on the bigger face, its start...
        std::optional<size_t> end, // ...and end index
        std::vector<std::vector<size_t>>& additional_polygons // output parameter to store any additional polygons in the interior
    ) {

        Node neighbour = node.adjacent_node(adj);

        if (neighbour.is_null() || neighbour.is_leaf()) // nothing to do here
            return std::vector<size_t>{};

        Point center = octree.barycenter(node);

        // collecting all the cells along the face
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

        // computing all edge intersections on the face
        std::vector<std::vector<const Edge*>> segments;
        int ind = 0, startv = -1, startp = -1;
        int endv = -1, endp = -1;
        while (!leaf_neighbours.empty()) {
            Node child = leaf_neighbours.front();
            std::array<Edge*, 12> child_edges = edges.get(child);
            std::vector<const Edge*> segment;
            for(int i = 0; i < 12; ++i) {
                auto corners = child_edges[i]->corners(child);
                // if both corners are on the shared face with `node`; 4 among the 12 will be
                bool c1in = !(corners.first & (4 / mask)) != !(adj & 1);
                bool c2in = !(corners.second & (4 / mask)) != !(adj & 1);
                if (c1in && c2in) {
                    auto vals = child_edges[i]->values();
                    if ((vals.first - isovalue) * (vals.second - isovalue) < 0) {
                        auto edge = child_edges[i]->find_minimal_edge(isovalue);
                        segment.push_back(edge);
                        if (start && add_point_to_mesh(segment.back()) == *start) {
                            startv = ind;
                            startp = segment.size() - 1;
                        }
                        if (end && add_point_to_mesh(segment.back()) == *end) {
                            endv = ind;
                            endp = segment.size() - 1;
                        }
                    }
                }
            }
            if(segment.size() > 0) {
                if(segment.size() != 2) { // this happens when there are 4 intersections on a square, currently no such test
                    assert(segment.size() == 4);
                    // we create two edges from the 4 points based on smallest edge length sum
                    FT min = std::numeric_limits<FT>::max();
                    int a;
                    for(int i = 0; i < 3; ++i) {
                        Point p1 = segment[0]->extract_isovertex(isovalue);
                        Point p2 = segment[i%3+1]->extract_isovertex(isovalue);
                        Point q1 = segment[(i+1)%3+1]->extract_isovertex(isovalue);
                        Point q2 = segment[(i+2)%3+1]->extract_isovertex(isovalue);
                        FT d =
                            sqrt((p1 - p2).squared_length())
                            + sqrt((q1 - q2).squared_length());
                        if (d < min) { min = d; a = i; }
                    }
                    segments.push_back(std::vector<const Edge*>{segment[0], segment[a%3+1]});
                    segments.push_back(std::vector<const Edge*>{segment[(a+1)%3+1], segment[(a+2)%3+1]});
                    ++ind;
                }
                else
                    segments.push_back(segment);
                ++ind;
            }
            leaf_neighbours.pop();
        }

        std::vector<const Edge*> polyline_edges;
        std::vector<size_t> polyline;
        std::set<int> vs;
        // if there is a coarser line, find its equivalent finer one
        if(start) {
            int v = startv, p = !startp;

            const Edge* first = segments[startv][startp];
            const Edge* next = segments[v][p];
            const Edge* last = segments[endv][endp];

            std::set<int> vs;
            vs.insert(v);
            while (next != last) {
                polyline_edges.push_back(next);
                bool found = false;
                for (int j = 0; j < segments.size(); ++j) if (v != j) {
                    if (segments[j][0] == polyline_edges.back()) { found = true; v = j; p = 1; break; }
                    else if (segments[j][1] == polyline_edges.back()) { found = true; v = j; p = 0; break; }
                }
                if (found) { // found next segment which is connected
                    next = segments[v][p];
                }
                if (!found) { // no connected segment; moving to twin edge
                    next = next->twin_edge(isovalue);
                    v = -1;
                }
                vs.insert(v);
            }
            polyline_edges.push_back(last);

            polyline.resize(polyline_edges.size());
            std::transform(polyline_edges.cbegin(), polyline_edges.cend(), polyline.begin(),
                [this] (const Edge* e) { return add_point_to_mesh(e); });
        }
        // finding additional segment loops among the remaining intersection points
        if(polyline.size() < segments.size()) {
            std::vector<std::vector<const Edge*>> remaining_segments;
            for (int i = 0; i < segments.size(); ++i) {
                if (vs.find(i) == vs.end()) remaining_segments.push_back(segments[i]);
            }

            while (remaining_segments.size() > 0) {
                auto poly = find_closed_polyline(remaining_segments);
                std::vector<size_t> polygon(poly.size());
                std::transform(poly.cbegin(), poly.cend(), polygon.begin(),
                    [this] (const Edge* e) { return add_point_to_mesh(e); });
                reorient_polygon(polygon, center);
                additional_polygons.push_back(polygon);
            }
        }
        if(polyline.size() > 0)
            polyline.pop_back();

        return polyline;

    }

    // processes a leaf cell of the octree
    // finds the polygon(s) in its inside and any additional polygons on its faces
    void process_node(const Node& node) {

        Bbox b = octree.bbox(node);

        std::vector<std::pair<size_t, std::pair<int,int>>> unordered_polygon_with_faces;

        // go through the edges and collect all intersections
        for(auto edge : edges.get(node)) {
            std::pair<Point,Point> seg = edge->segment();

            if((edge->values().first - isovalue) * (edge->values().second - isovalue) < 0) {
                auto minEdge = edge->find_minimal_edge(isovalue);
                std::pair<Point,Point> minSeg = minEdge->segment();
                auto vals = minEdge->values();
                auto corners = edge->corners(node);
                if (vals.second - isovalue < 0) { auto tmp = corners.first; corners.first = corners.second; corners.second = tmp; }
                unordered_polygon_with_faces.push_back(std::make_pair(
                    add_point_to_mesh(minEdge),
                    corners_to_faces(corners)));
            }
        }

        // the isosurface inside the cell can be composed of multiple disconnected part
        std::vector<std::vector<size_t>> polygons(1);
        std::vector<std::vector<size_t>> additional_polygons;
        std::vector<bool> processed_faces (6, false);

        size_t loop = 0, ind = 0;
        int prevFace = -1, lastFace = -1;
        std::vector<size_t> visited {ind};
        for (size_t i = 0; i < unordered_polygon_with_faces.size(); ++i) {
            polygons[loop].push_back(unordered_polygon_with_faces[ind].first);
            visited.push_back(ind);
            std::pair<int,int> faces = unordered_polygon_with_faces[ind].second;
            if(polygons[loop].size() == 1) lastFace = faces.second;
            bool foundnext = false;
            std::vector<int> sameFace; // collect points also on the next face
            int nextFace = -1;
            for (size_t j = 0; j < unordered_polygon_with_faces.size(); ++j) {
                auto& el = unordered_polygon_with_faces[j];
                bool face1 = el.second.first == faces.first
                    || ((polygons[loop].size() != 1) && el.second.first == faces.second); // at the first vertex we want to start in the good direction
                bool face2 = el.second.second == faces.first
                    || ((polygons[loop].size() != 1) && el.second.second == faces.second);
                if (ind != j
                    && ((face1 && el.second.first != prevFace)
                        || (face2 && el.second.second != prevFace))) {
                    sameFace.push_back(j);
                    nextFace = (face1 ? el.second.first : el.second.second);
                }
            }
            if (nextFace != -1) prevFace = nextFace;
            if (!sameFace.empty()) { // there is a point on the next face
                foundnext = true;
                if(sameFace.size() == 1) { // only one point: choose that
                    ind = sameFace[0];
                }
                else { // 3 other intersection points: choose one based on edge lengths
                    assert(sameFace.size() == 3);
                    FT min = std::numeric_limits<FT>::max();
                    int next;
                    for(int i = 0; i < 3; ++i) {
                        Point curr = vertices[polygons[loop].back()];
                        Point connecting = vertices[unordered_polygon_with_faces[sameFace[i]].first];
                        Point other1 = vertices[unordered_polygon_with_faces[sameFace[(i+1)%3]].first];
                        Point other2 = vertices[unordered_polygon_with_faces[sameFace[(i+2)%3]].first];
                        FT d =
                            sqrt((curr - connecting).squared_length())
                            + sqrt((other1 - other2).squared_length());
                        if (d < min) { min = d; next = sameFace[i]; }
                    }
                    ind = next;
                }
                if (std::find(visited.begin(), visited.end(), ind) != visited.end()) {
                    foundnext = false; // next point has already been visited: loop is closed
                }
                else{
                    // check if other side of the face has a finer polyline
                    std::vector<size_t> internal_points = process_face(
                        node, Adjacency(nextFace), polygons[loop].back(), unordered_polygon_with_faces[ind].first,
                        additional_polygons);
                    processed_faces[nextFace] = true;
                    for (auto it : internal_points) {
                        polygons[loop].push_back(it);
                    }
                    visited.push_back(ind);
                }
            }
            if (!foundnext) { // loop closed
                // check last face (between last and first point)
                std::vector<size_t> internal_points = process_face(
                    node, Adjacency(lastFace), polygons[loop].back(), polygons[loop][0],
                    additional_polygons);
                processed_faces[lastFace] = true;
                for (auto it : internal_points) {
                    polygons[loop].push_back(it);
                }
                // if there are still points find a starting point for a next loop
                if (i < unordered_polygon_with_faces.size() - 1) {
                    ++loop;
                    polygons.resize(loop+1);
                    prevFace = -1;
                    for(size_t j = 0; j < unordered_polygon_with_faces.size(); ++j)
                        if (std::find(visited.begin(), visited.end(), j) == visited.end()) {
                            ind = j;
                            visited.push_back(ind);
                            break;
                        }
                }
            }
        }

        // check for additional polygons on the yet unchecked faces
        for(int i = 0; i < 6; ++i) {
            if(!processed_faces[Adjacency(i)]) {
                process_face(node, Adjacency(i), std::nullopt, std::nullopt, additional_polygons);
            }
        }

        for (auto it : additional_polygons) {
            faces.push_back(it);
        }

        for (auto it : polygons) {
            if (it.size() > 0)
                faces.push_back(it);
        }
    }

    void operator()(const typename Octree::Voxel_handle& vox) {
        process_node(octree.get_node(vox));
    }

    std::vector<Point> get_vertices() const { return vertices; }
    std::vector<std::vector<size_t>> get_faces() const { return faces; }

private:

    const Octree& octree;
    Edge_store edges;
    FT isovalue;

    std::unordered_map<const Edge*, size_t> edge_points_in_mesh;

    std::vector<Point> vertices;
    std::vector<std::vector<size_t>> faces;

    static constexpr int corners_to_faces_table[8][3][2] = {
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
    static std::pair<int,int> corners_to_faces(std::pair<int,int> corners) {
        int p0 = corners.first, p1 = corners.second;
        int r0 = std::abs(p0 - p1) / 2; // 1,2,4 -> 0,1,2
        return std::make_pair(corners_to_faces_table[p0][r0][0], corners_to_faces_table[p0][r0][1]);
    }

};

}

#endif