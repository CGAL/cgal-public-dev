// All rights reserved.
// Copyright (c) 2020 SARL GeometryFactory (France).
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
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri

#ifndef CGAL_URBAN_AREA_PROCESSING_INTERNAL_SHORTEST_PATH_CONTOURING_2_H
#define CGAL_URBAN_AREA_PROCESSING_INTERNAL_SHORTEST_PATH_CONTOURING_2_H

// #include <CGAL/license/Urban_area_processing.h>

// STL includes.
#include <vector>
#include <utility>
#include <algorithm>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/squared_distance_2.h>

// Internal includes.
#include <CGAL/Urban_area_processing/internal/Dijkstra.h>
#include <CGAL/Urban_area_processing/internal/Sphere_neighbor_query.h>

// TODO:
// Try to optimize this class. Do we need Dijkstra at all?
// This class should also handle closed contours. For the moment, the last segment
// is missing.

namespace CGAL {
namespace Urban_area_processing {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename SegmentMap>
  class Shortest_path_contouring_2 {
  
  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;
    using Size_pair = std::pair<std::size_t, std::size_t>;
    using Points_2 = std::vector<Point_2>;
    using Indices = std::vector<std::size_t>;

    using Dijkstra = Dijkstra<Traits>;
    using Node = typename Dijkstra::Node;
    using Nodes = typename Dijkstra::Nodes;
    using Neighbor = typename Dijkstra::Neighbor;
    using Node_map = typename Dijkstra::Node_map;
    using Vertex_t = typename Dijkstra::Vertex_t;

    struct Segment {

      Segment_2 segment;
      std::size_t index = std::size_t(-1);
      const FT squared_length() const {
        return segment.squared_length();
      }
    };

    using Sphere_neighbor_query =
      internal::Sphere_neighbor_query<Traits, Nodes, Node_map>;

    Shortest_path_contouring_2(
      const Input_range& input_range,
      const Segment_map segment_map,
      const FT scale,
      const FT min_length_2,
      const bool verbose = true) :
    m_input_range(input_range),
    m_segment_map(segment_map),
    m_scale(scale),
    m_min_length_2(min_length_2),
    m_verbose(verbose) { 

      CGAL_precondition(input_range.size() > 0);
    }

    template<typename OutputIterator>
    void merge(OutputIterator contours) {

      std::vector<Segment> wraps;
      wrap_segments(wraps);
      if (m_verbose)
        std::cout << "- segments are wrapped: " << wraps.size() << std::endl;

      std::vector<Node> nodes;
      std::map<std::size_t, Size_pair> nmap; 
      create_nodes(wraps, nmap, nodes);
      if (m_verbose)
        std::cout << "- nodes are created: " << nodes.size() << std::endl;

      std::vector<Points_2> components;
      create_connected_components(
        nmap, nodes, components);
      if (m_verbose)
        std::cout << "- components are created: " << components.size() << std::endl;

      output_contours(components, contours);
      if (m_verbose)
        std::cout << "- contours are created" << std::endl;
    }

  private:
    const Input_range& m_input_range;
    const Segment_map m_segment_map;
    const FT m_scale;
    const FT m_min_length_2;
    const bool m_verbose;

    void wrap_segments(
      std::vector<Segment>& wraps) const {

      wraps.clear();
      wraps.reserve(m_input_range.size());
      
      Segment wrap;
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const auto& segment = get(
          m_segment_map, *(m_input_range.begin() + i));

        wrap.index = i;
        wrap.segment = segment;
        wraps.push_back(wrap);
      }

      std::sort(wraps.begin(), wraps.end(), 
      [](const Segment& a, const Segment& b) -> bool { 
        return a.squared_length() > b.squared_length();
      });
    }

    void create_nodes(
      const std::vector<Segment>& wraps,
      std::map<std::size_t, Size_pair>& nmap,
      std::vector<Node>& nodes) const {
      
      nmap.clear(); nodes.clear();
      nodes.reserve(wraps.size() * 2);
      
      Node node1, node2; std::size_t count = 0;
      for (const auto& wrap : wraps) {
        const auto& segment = wrap.segment;
        
        const auto& source = segment.source();
        const auto& target = segment.target();

        node1.point = source;
        node1.segment_index = wrap.index;
        node1.index = count; ++count;
        nodes.push_back(node1);

        node2.point = target;
        node2.segment_index = wrap.index;
        node2.index = count; ++count;
        nodes.push_back(node2);

        nmap[wrap.index] = 
          std::make_pair(node1.index, node2.index);
      }
    }

    void create_connected_components(
      const std::map<std::size_t, Size_pair>& nmap,
      std::vector<Node>& nodes,
      std::vector<Points_2>& components) const {

      components.clear();
      create_close_neighbors(nodes);

      auto sorted = nodes;
      std::sort(sorted.begin(), sorted.end(), 
      [](const Node& a, const Node& b) -> bool { 
        return a.neighbors.size() < b.neighbors.size();
      });

      std::size_t count = 0;
      std::vector<Node> component;

      std::vector<Point_2> contour;
      for (const auto& item : sorted) {
        auto& node = nodes[item.index];
        if (node.used) continue;

        component.clear();
        component.push_back(node);
        node.used = true;

        const auto& pair = nmap.at(node.segment_index);
        add_neighbors(nmap, nodes[pair.first], nodes, component);
        add_neighbors(nmap, nodes[pair.second], nodes, component);

        for (std::size_t i = 0; i < component.size(); ++i) {
          auto& cnode = component[i];
          cnode.index = i; 
          cnode.neighbors.clear();
        }
        
        std::list<Vertex_t> path;
        FT path_length = FT(-1);

        std::vector<int> result;
        std::list<Vertex_t> temp;
        create_all_neighbors(component);

        Dijkstra dijkstra(component);
        for (std::size_t i = 0; i < component.size(); ++i) {
          for (std::size_t j = 0; j < component.size(); ++j) {      
            dijkstra.apply(i, j, temp);
            
            result.clear();
            for (const int index : temp) 
              result.push_back(index);
            if (result.size() == 0) continue;

            FT length = FT(0);
            for (std::size_t k = 0; k < result.size() - 1; ++k) {
              const std::size_t kp = (k + 1) % result.size();
              length += internal::distance(
                component[result[k]].point, component[result[kp]].point);
            }

            if (length > path_length) {
              path = temp; path_length = length;
            }
          }
        }

        if (path_length > m_min_length_2) {
          contour.clear();
          get_contour(path, component, contour);
          components.push_back(contour);
        }
        ++count;
      }
    }

    void create_close_neighbors(
      std::vector<Node>& nodes) const {

      // Create a neighbor query.
      Node_map node_map;
      const FT radius = m_scale;
      Sphere_neighbor_query neighbor_query(
        nodes, radius, node_map);

      // Find neighbors.
      Indices indices;
      for (auto& nodei : nodes) {
        const std::size_t i = nodei.index;
        const std::size_t si = nodei.segment_index;
        auto& neighbors = nodei.neighbors;

        const auto& query = nodei.point;
        neighbor_query(query, indices);

        for (const std::size_t index : indices) {
          const auto& nodej = nodes[index];
          const std::size_t j = nodej.index;
          if (i == j) continue;

          const std::size_t sj = nodej.segment_index;
          if (si == sj) {
            const FT weight = FT(0);
            neighbors.push_back(Neighbor(j, weight)); continue;
          }

          const FT weight = CGAL::squared_distance(
            nodei.point, nodej.point);
          neighbors.push_back(Neighbor(j, weight));
        }
      }
    }

    void create_all_neighbors(
      std::vector<Node>& nodes) const {

      for (auto& nodei : nodes) {
        const std::size_t i = nodei.index;
        const std::size_t si = nodei.segment_index;
        auto& neighbors = nodei.neighbors;

        for (const auto& nodej : nodes) {
          const std::size_t j = nodej.index;
          if (i == j) continue;

          const std::size_t sj = nodej.segment_index;
          if (si == sj) {
            const FT weight = FT(0);
            neighbors.push_back(Neighbor(j, weight)); 
            continue;
          }

          const FT weight = CGAL::squared_distance(
            nodei.point, nodej.point);
          neighbors.push_back(Neighbor(j, weight));
        }
      }
    }

    void add_neighbors(
      const std::map<std::size_t, Size_pair>& nmap,
      Node& nodei,
      std::vector<Node>& nodes,
      std::vector<Node>& component) const {

      if (nodei.used) return;
      nodei.used = true;
      component.push_back(nodei);
      const auto& neighbors = nodei.neighbors;

      for (const auto& neighbor : neighbors) {
        auto& nodej = nodes[neighbor.target];
        if (nodej.used) continue;
        
        component.push_back(nodej);
        nodej.used = true;

        const auto& pair = nmap.at(nodej.segment_index);
        add_neighbors(nmap, nodes[pair.first], nodes, component);
        add_neighbors(nmap, nodes[pair.second], nodes, component);
      }
    }

    void get_contour(
      const std::list<Vertex_t>& path,
      const std::vector<Node>& nodes,
      std::vector<Point_2>& contour) const {
      
      contour.clear();
      contour.reserve(path.size());
      for (const int index : path)
        contour.push_back(nodes[index].point);
    }

    template<typename OutputIterator>
    void output_contours(
      const std::vector<Points_2>& components,
      OutputIterator contours) const {

      std::vector<Segment_2> contour;
      for (const auto& component : components) {
        if (component.size() == 0) continue;
        
        contour.clear();
        for (std::size_t i = 0; i < component.size() - 1; ++i) {
          const std::size_t ip = i + 1;
          const auto& source = component[i];
          const auto& target = component[ip];
          contour.push_back(Segment_2(source, target));
        }
        *(++contours) = contour;
      }
    }
  };

} // internal
} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_INTERNAL_SHORTEST_PATH_CONTOURING_2_H
