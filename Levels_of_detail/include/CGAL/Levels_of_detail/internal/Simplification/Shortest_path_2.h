// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is a part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY, AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri
//

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_SHORTEST_PATH_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_SHORTEST_PATH_2_H

// STL includes.
#include <map>
#include <utility>
#include <iostream>
#include <vector>
#include <string>
#include <list>
#include <limits>
#include <set>
#include <algorithm>
#include <iterator>

// CGAL includes.
#include <CGAL/barycenter.h>
#include <CGAL/property_map.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>

// Spatial search.
#include <CGAL/Levels_of_detail/internal/Spatial_search/Sphere_neighbor_query.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Shortest_path_2 {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Segment_2 = typename Traits::Segment_2;

    using Vertex_t = int;
    using Weight_t = FT;

    using Pair = std::pair<Point_2, std::size_t>;
    using Point_map = CGAL::First_of_pair_property_map<Pair>;

    using Sphere_neighbor_query =
    internal::Sphere_neighbor_query<Traits, std::vector<Pair>, Point_map>;

    const Weight_t max_weight = std::numeric_limits<FT>::infinity();

    struct Neighbor {

      Vertex_t target;
      Weight_t weight;

      Neighbor(
        Vertex_t arg_target, Weight_t arg_weight) :
      target(arg_target), weight(arg_weight) { }
    };

    using Adjacency_list_t = std::vector< std::vector<Neighbor> >;

    struct GSegment {

      std::size_t idx = std::size_t(-1);
      Segment_2 segment;
      std::vector<Point_2> points;

      const FT distance() const {
        return internal::distance(segment.source(), segment.target());
      }
    };

    struct GNode {

      std::size_t idx = std::size_t(-1);
      std::size_t seg_idx = std::size_t(-1);
      Point_2 point;
      std::vector<Neighbor> neighbors;
      bool used = false;
    };

    struct GComponent {
      std::size_t idx = std::size_t(-1);
      std::vector<Point_2> contour;
      std::vector<GNode> gnodes;
    };

    Shortest_path_2(
      const FT noise_level,
      const FT min_length) :
    m_noise_level(noise_level),
    m_min_length(min_length)
    { }

    void find(
      const std::vector< std::vector<Point_2> >& regions,
      const std::vector<Segment_2>& segments,
      std::vector< std::vector<Segment_2> >& contours) const {

      std::vector<GSegment> gsegments;
      create_graph_segments(regions, segments, gsegments);

      std::vector<GNode> gnodes;
      std::map<std::size_t, std::pair<std::size_t, std::size_t> > gmap;
      create_graph_nodes_gsegments(gsegments, gmap, gnodes);

      std::vector<GComponent> gcomponents;
      std::vector< std::vector<Point_2> > out_contours;
      create_connected_components_gnodes(
        gmap, gnodes, gcomponents, out_contours);

      std::vector<Segment_2> out_segments;
      for (const auto& out_contour : out_contours) {
        if (out_contour.size() == 0) continue;

        out_segments.clear();
        for (std::size_t i = 0; i < out_contour.size() - 1; ++i) {
          const std::size_t ip = i + 1;
          const auto& s = out_contour[i];
          const auto& t = out_contour[ip];
          out_segments.push_back(Segment_2(s, t));
        }
        contours.push_back(out_segments);
      }
    }

  private:
    const FT m_noise_level;
    const FT m_min_length;

    void create_graph_segments(
      const std::vector< std::vector<Point_2> >& regions,
      const std::vector<Segment_2>& segments,
      std::vector<GSegment>& gsegments) const {

      CGAL_assertion(segments.size() == regions.size());

      gsegments.clear();
      gsegments.reserve(segments.size());

      GSegment gsegment;
      for (std::size_t i = 0; i < segments.size(); ++i) {
        const auto& segment = segments[i];

        gsegment.segment = segment;
        gsegment.idx     = i;
        gsegment.points  = regions[i];

        gsegments.push_back(gsegment);
      }

      std::sort(gsegments.begin(), gsegments.end(),
      [](const GSegment& a, const GSegment& b) -> bool {
        return a.distance() > b.distance();
      });
    }

    void create_graph_nodes_gsegments(
      const std::vector<GSegment>& gsegments,
      std::map<std::size_t, std::pair<std::size_t, std::size_t> >& gmap,
      std::vector<GNode>& gnodes) const {

      gnodes.clear();
      gnodes.reserve(gsegments.size() * 2);
      gmap.clear();

      GNode gnode1, gnode2; std::size_t count = 0;
      for (std::size_t i = 0; i < gsegments.size(); ++i) {
        const auto& segment = gsegments[i].segment;

        const auto& s = segment.source();
        const auto& t = segment.target();

        gnode1.seg_idx = i; gnode1.point = s;
        gnode1.idx = count; ++count;
        gnodes.push_back(gnode1);

        gnode2.seg_idx = i; gnode2.point = t;
        gnode2.idx = count; ++count;
        gnodes.push_back(gnode2);

        gmap[i] = std::make_pair(gnode1.idx, gnode2.idx);
      }
    }

    void create_connected_components_gnodes(
      const std::map<std::size_t, std::pair<std::size_t, std::size_t> >& gmap,
      std::vector<GNode>& gnodes,
      std::vector<GComponent>& gcomponents,
      std::vector< std::vector<Point_2> >& contours) const {

      contours.clear();
      const FT radius = m_noise_level / FT(2);
      create_graph_neighbors_tree(
        radius, false, -1, -1, gnodes);

      std::vector<GNode> sorted = gnodes;
      std::sort(sorted.begin(), sorted.end(),
      [](const GNode& a, const GNode& b) -> bool {
        return a.neighbors.size() < b.neighbors.size();
      });

      std::size_t count = 0;
      std::vector<GNode> component;
      for (auto& s : sorted) {
        const std::size_t idx = s.idx;
        auto& gnode = gnodes[idx];
        if (gnode.used) continue;

        component.clear();
        component.push_back(gnode);
        gnode.used = true;

        const std::size_t seg_idx = gnode.seg_idx;
        const auto& pair = gmap.at(seg_idx);
        add_neighbors(gmap, gnodes[pair.first], gnodes, component);
        add_neighbors(gmap, gnodes[pair.second], gnodes, component);

        for (std::size_t i = 0; i < component.size(); ++i) {
          auto& gn = component[i];
          gn.idx = i; gn.neighbors.clear();
        }

        std::list<Vertex_t> path;
        FT path_length = FT(-1);

        std::vector<int> clean;
        std::list<Vertex_t> tmp;
        create_graph_neighbors_all(false, -1, -1, component);

        for (std::size_t j = 0; j < component.size(); ++j) {
          for (std::size_t i = 0; i < component.size(); ++i) {
            apply_dijkstra(j, i, component, tmp);

            clean.clear();
            for (const int idx : tmp) clean.push_back(idx);
            if (clean.size() == 0) continue;

            FT length = FT(0);
            for (std::size_t i = 0; i < clean.size() - 1; ++i) {
              const std::size_t ip = (i + 1) % clean.size();
              length += internal::distance(
                component[clean[i]].point, component[clean[ip]].point);
            }

            if (length > path_length) {
              path = tmp; path_length = length;
            }
          }
        }

        if (path_length > m_min_length) {

          std::vector<Point_2> contour;
          get_contour(path, component, contour);
          contours.push_back(contour);
        }
        ++count;
      }
    }

    void add_neighbors(
      const std::map<std::size_t, std::pair<std::size_t, std::size_t> >& gmap,
      GNode& gnode,
      std::vector<GNode>& gnodes,
      std::vector<GNode>& component) const {

      if (gnode.used) return;
      gnode.used = true;
      component.push_back(gnode);
      const auto& neighbors = gnode.neighbors;

      for (const auto& neighbor : neighbors) {
        auto& gnode = gnodes[neighbor.target];
        if (gnode.used) continue;

        component.push_back(gnode);
        gnode.used = true;

        const std::size_t seg_idx = gnode.seg_idx;
        const auto& pair = gmap.at(seg_idx);
        add_neighbors(gmap, gnodes[pair.first], gnodes, component);
        add_neighbors(gmap, gnodes[pair.second], gnodes, component);
      }
    }

    void create_graph_neighbors_tree(
      const FT radius,
      const bool connected,
      const int start, const int end,
      std::vector<GNode>& gnodes) const {

      Point_map point_map;
      std::vector<Pair> pairs;
      pairs.reserve(gnodes.size());
      for (const auto& gnode : gnodes)
        pairs.push_back(std::make_pair(gnode.point, gnode.idx));

      std::vector<std::size_t> neighbors;
      Sphere_neighbor_query snq(pairs, radius, point_map);

      for (auto& gnode_i : gnodes) {
        const std::size_t i = gnode_i.idx;
        const std::size_t idx_i = gnode_i.seg_idx;
        auto& neighbors_i = gnode_i.neighbors;

        const auto& query = gnode_i.point;
        snq(query, neighbors);

        for (const std::size_t idx : neighbors) {
          const auto& gnode_j = gnodes[idx];
          const std::size_t j = gnode_j.idx;
          if (i == j) continue;

          if (connected) {
            if (i == start && j == end) {
              neighbors_i.push_back(Neighbor(j, max_weight));
              continue;
            }
          }

          const std::size_t idx_j = gnode_j.seg_idx;
          if (idx_i == idx_j) {
            neighbors_i.push_back(Neighbor(j, 0.0));
            continue;
          }

          const FT distance = internal::distance(gnode_i.point, gnode_j.point);
          neighbors_i.push_back(Neighbor(j, distance));
        }
      }
    }

    void create_graph_neighbors_all(
      const bool connected,
      const int start, const int end,
      std::vector<GNode>& gnodes) const {

      for (auto& gnode_i : gnodes) {
        const std::size_t i = gnode_i.idx;
        const std::size_t idx_i = gnode_i.seg_idx;
        auto& neighbors_i = gnode_i.neighbors;

        for (const auto& gnode_j : gnodes) {
          const std::size_t j = gnode_j.idx;
          if (i == j) continue;

          if (connected) {
            if (i == start && j == end) {
              neighbors_i.push_back(Neighbor(j, max_weight));
              continue;
            }
          }

          const std::size_t idx_j = gnode_j.seg_idx;
          if (idx_i == idx_j) {
            neighbors_i.push_back(Neighbor(j, 0.0));
            continue;
          }

          const FT distance = internal::distance(gnode_i.point, gnode_j.point);
          neighbors_i.push_back(Neighbor(j, distance));
        }
      }
    }

    void apply_dijkstra(
      const int start, const int end,
      const std::vector<GNode>& gnodes,
      std::list<Vertex_t>& path) const {

      Adjacency_list_t adjacency_list;
      adjacency_list.reserve(gnodes.size());
      for (const auto& gnode : gnodes)
        adjacency_list.push_back(gnode.neighbors);

      std::vector<Weight_t> min_distance;
      std::vector<Vertex_t> previous;

      dijkstra_compute_paths(
        start, adjacency_list, min_distance, previous);
      path = dijkstra_get_shortest_path_to(end, previous);

      /*
      std::cout << "Shortest path found : ";
      std::copy(
        path.begin(), path.end(),
        std::ostream_iterator<Vertex_t>(std::cout, " "));
      std::cout << std::endl; */
    }

    void dijkstra_compute_paths(
      Vertex_t source,
      const Adjacency_list_t& adjacency_list,
      std::vector<Weight_t>& min_distance,
      std::vector<Vertex_t>& previous) const {

      int n = adjacency_list.size();

      min_distance.clear();
      min_distance.resize(n, max_weight);
      min_distance[source] = 0;

      previous.clear();
      previous.resize(n, -1);

      std::set< std::pair<Weight_t, Vertex_t> > vertex_queue;
      vertex_queue.insert(
        std::make_pair(min_distance[source], source));

      while (!vertex_queue.empty()) {

        Weight_t dist = vertex_queue.begin()->first;
        Vertex_t u = vertex_queue.begin()->second;
        vertex_queue.erase(vertex_queue.begin());

        const std::vector<Neighbor>& neighbors = adjacency_list[u];
        for (
          typename std::vector<Neighbor>::const_iterator neighbor_iter = neighbors.begin();
          neighbor_iter != neighbors.end();
          neighbor_iter++) {

          Vertex_t v = neighbor_iter->target;
          Weight_t weight = neighbor_iter->weight;
          Weight_t distance_through_u = dist + weight;

          if (distance_through_u < min_distance[v]) {

            vertex_queue.erase(std::make_pair(min_distance[v], v));
            min_distance[v] = distance_through_u;
            previous[v] = u;
            vertex_queue.insert(std::make_pair(min_distance[v], v));
          }
        }
      }
    }

    std::list<Vertex_t> dijkstra_get_shortest_path_to(
      Vertex_t vertex,
      const std::vector<Vertex_t>& previous) const {

      std::list<Vertex_t> path;
      for ( ; vertex != -1; vertex = previous[vertex])
        path.push_front(vertex);
      return path;
    }

    void get_contour(
      const std::list<Vertex_t>& path,
      const std::vector<GNode>& gnodes,
      std::vector<Point_2>& contour) const {

      contour.clear();
      contour.reserve(path.size());
      for (const int idx : path)
        contour.push_back(gnodes[idx].point);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_SHORTEST_PATH_2_H
