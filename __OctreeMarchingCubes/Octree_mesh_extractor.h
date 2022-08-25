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

    typedef Octree_wrapper<GeomTraits> Octree;
    typedef CGAL::Orthtrees::Preorder_traversal Preorder_traversal;
    typedef typename CGAL::Bbox_3 Bbox;
    typedef typename Octree::Node Node;
    typedef CGAL::Octree_edge<GeomTraits> Edge;

    typedef typename CGAL::Orthtree_traits_3<GeomTraits>::Adjacency Adjacency;

    typedef std::function<FT(Vector)> Implicit_function;

    typedef CGAL::Edge_store<GeomTraits, Point_set, Point_map> Edge_store;

    Octree_mesh_extractor(const Octree& tree, FT isovalue)
    : octree(tree), edges(tree), isovalue(isovalue) {
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

        Point p = e->extract_isovertex(isovalue);
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
        //std::cerr << *(polyline.back());
        int v = 0, p = 1;
        std::set<int> vs;
        vs.insert(v);
        while (next != start) {
            polyline.push_back(next);
            //std::cerr << *(polyline.back());
            bool found = false;
            for (int j = 0; j < segments.size(); ++j) if (v != j) {
                if (segments[j][0] == polyline.back()) { found = true; v = j; p = 1; break; }
                else if (segments[j][1] == polyline.back()) { found = true; v = j; p = 0; break; }
            }
            if (found) {
                next = segments[v][p];
            }
            if (!found) {
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
        const Node& node,
        Adjacency adj,
        std::optional<size_t> start,
        std::vector<std::vector<size_t>>& additional_polygons) {

        //std::cerr << "processing face #" << adj << " of node " << node << " (startpoint given: " << start.has_value() << ")" << std::endl;

        Node neighbour = node.adjacent_node(adj);

        if (neighbour.is_null() || neighbour.is_leaf())
            return std::vector<size_t>{};

        //std::cerr << "face is subdivided" << std::endl;
        Point center = octree.barycenter(node);

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
                auto corners = child_edges[i]->corners(child);
                // if both corners are on the shared face with `node`
                bool c1in = !(corners.first & (4 / mask)) != !(adj & 1);
                bool c2in = !(corners.second & (4 / mask)) != !(adj & 1);
                if (c1in && c2in) {
                    auto edge = child_edges[i]->find_minimal_edge(isovalue);
                    auto vals = edge->values();
                    if ((vals.first - isovalue) * (vals.second - isovalue) < 0) {
                        segment.push_back(edge);
                        if (start && add_point_to_mesh(segment.back()) == *start) {
                            startv = ind;
                            startp = segment.size() - 1;
                        }
                    }
                }
            }
            if(segment.size() > 0) {
                if(segment.size() != 2) std::cerr << segment.size() << "\n"; // this happens when there are 4 intersections on a square, currently no such test
                segments.push_back(segment);
                ++ind;
            }
            leaf_neighbours.pop();
        }

        std::vector<const Edge*> polyline_edges;
        std::vector<size_t> polyline;
        std::set<int> vs;

        if(start) {
            int v = startv, p = !startp;
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
                    auto next = polyline_edges.back()->twin_edge(isovalue);
                    if (next == nullptr) break; // no twin edge: polyline to be ended
                    v = -1;
                    for (int j = 0; j < segments.size(); ++j) {
                        if (segments[j][0] == next || segments[j][1] == next) { v = j; }
                    }
                    if (v == -1) break; // twin edge exists, but not in this cell: polyline is ended
                }
            }
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
            /*std::cout << "Not all points were used." << std::endl
            << "Remaining segments: " << remaining_segments.size() << std::endl;*/

            //std::cerr << "Edges already used: " << std::endl;
            /*for(auto it : polyline_edges)
                std::cerr << *it << std::endl;*/

            //std::cerr << "Node: " << node << std::endl;
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

    // returns the vertices extracted from a given cell, ordered as a polyline
    void process_node(const Node& node) {

        Bbox b = octree.bbox(node);

        std::vector<std::pair<size_t, std::pair<int,int>>> unordered_polygon_with_faces;

        for(auto edge : edges.get(node)) {
            std::pair<Point,Point> seg = edge->segment();
            Vector p1 (seg.first.x(), seg.first.y(), seg.first.z());
            Vector p2 (seg.second.x(), seg.second.y(), seg.second.z());

            if((edge->values().first - isovalue) * (edge->values().second - isovalue) < 0) {
                auto minEdge = edge->find_minimal_edge(isovalue);
                std::pair<Point,Point> minSeg = minEdge->segment();
                Vector p1 (minSeg.first.x(), minSeg.first.y(), minSeg.first.z());
                Vector p2 (minSeg.second.x(), minSeg.second.y(), minSeg.second.z());
                auto vals = minEdge->values();
                auto corners = edge->corners(node);
                if (vals.second - isovalue < 0) { auto tmp = corners.first; corners.first = corners.second; corners.second = tmp; }
                unordered_polygon_with_faces.push_back(std::make_pair(
                    add_point_to_mesh(minEdge),
                    corners_to_faces(corners)));
            }
        }

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
            std::vector<int> sameFace;
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
            if (!sameFace.empty()) {
                foundnext = true;
                if(sameFace.size() == 1) {
                    ind = sameFace[0];
                }
                else {
                    assert(sameFace.size() == 3);
                    double min = std::numeric_limits<double>::max();
                    int next;
                    for(int i = 0; i < 3; ++i) {
                        Point curr = vertices[polygons[loop].back()];
                        Point connecting = vertices[unordered_polygon_with_faces[sameFace[i]].first];
                        Point other1 = vertices[unordered_polygon_with_faces[sameFace[(i+1)%3]].first];
                        Point other2 = vertices[unordered_polygon_with_faces[sameFace[(i+2)%3]].first];
                        double d =
                            sqrt((curr - connecting).squared_length())
                            + sqrt((other1 - other2).squared_length());
                        if (d < min) { min = d; next = sameFace[i]; }
                    }
                    ind = next;
                }
                if (std::find(visited.begin(), visited.end(), ind) != visited.end()) {
                    foundnext = false;
                }
                else{
                    std::vector<size_t> internal_points = process_face(
                        node, Adjacency(nextFace), polygons[loop].back(),
                        additional_polygons);
                    processed_faces[nextFace] = true;
                    for (auto it : internal_points) {
                        polygons[loop].push_back(it);
                    }
                    visited.push_back(ind);
                }
            }
            if (!foundnext) {
                std::vector<size_t> internal_points = process_face(
                    node, Adjacency(lastFace), polygons[loop].back(),
                    additional_polygons);
                processed_faces[lastFace] = true;
                for (auto it : internal_points) {
                    polygons[loop].push_back(it);
                }
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

        for(int i = 0; i < 6; ++i) {
            if(!processed_faces[Adjacency(i)]) {
                process_face(node, Adjacency(i), std::nullopt, additional_polygons);
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

    std::vector<Point> get_vertices() { return vertices; }
    std::vector<std::vector<size_t>> get_faces() { return faces; }

private:

    const Octree& octree;
    Edge_store edges;
    FT isovalue;

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