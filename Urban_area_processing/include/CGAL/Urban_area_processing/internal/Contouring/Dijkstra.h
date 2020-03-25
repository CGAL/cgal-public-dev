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

#ifndef CGAL_URBAN_AREA_PROCESSING_INTERNAL_DIJKSTRA_H
#define CGAL_URBAN_AREA_PROCESSING_INTERNAL_DIJKSTRA_H

// #include <CGAL/license/Urban_area_processing.h>

// STL includes.
#include <set>
#include <list>
#include <vector>
#include <utility>
#include <algorithm>

// CGAL includes.
#include <CGAL/assertions.h>

// TODO:
// Try to use the boost version here.

namespace CGAL {
namespace Urban_area_processing {
namespace internal {

  template<typename GeomTraits>
  class Dijkstra {
  
  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;

    using Vertex_t = int;
    using Weight_t = FT;

    using Pair_t = std::pair<Weight_t, Vertex_t>;

    struct Neighbor {
      
      Vertex_t target;
      Weight_t weight;
      
      Neighbor(
        Vertex_t arg_target, 
        Weight_t arg_weight) : 
      target(arg_target), 
      weight(arg_weight) 
      { }
    };

    using Neighbors = std::vector<Neighbor>;
    using Adjacency_list_t = std::vector<Neighbors>;

    struct Node {

      std::size_t index = std::size_t(-1);
      std::size_t segment_index = std::size_t(-1);
      Point_2 point;
      Neighbors neighbors;
      bool used = false;
    };

    using Nodes = std::vector<Node>;

    struct Node_map {
      using key_type = Node;
      using value_type = Point_2;
      using reference = const value_type&;
      using category = boost::lvalue_property_map_tag;

      value_type& operator[](key_type& key) const { 
        return key.point; 
      }

      friend reference get(
        const Node_map&, const key_type& key) {
        return key.point;
      }
    };

    Dijkstra(
      const std::vector<Node>& nodes) :
    m_nodes(nodes),
    m_max_weight(internal::max_value<FT>()) 
    { }

    void apply(
      const int start, 
      const int end,
      std::list<Vertex_t>& path) const {

      Adjacency_list_t adjacency_list;
      adjacency_list.reserve(m_nodes.size());
      for (const auto& node : m_nodes)
        adjacency_list.push_back(node.neighbors);

      std::vector<Weight_t> min_distance;
      std::vector<Vertex_t> previous;

      dijkstra_compute_paths(
        start, adjacency_list, min_distance, previous);
      path = dijkstra_get_shortest_path_to(end, previous);
    }

  private:
    const std::vector<Node>& m_nodes;
    const Weight_t m_max_weight;

    void dijkstra_compute_paths(
      Vertex_t source,
      const Adjacency_list_t& adjacency_list,
      std::vector<Weight_t>& min_distance,
      std::vector<Vertex_t>& previous) const {
        
      const int n = adjacency_list.size();

      min_distance.clear();
      min_distance.resize(n, m_max_weight);
      min_distance[source] = 0;

      previous.clear();
      previous.resize(n, -1);

      std::set<Pair_t> vertex_queue;
      vertex_queue.insert(
        std::make_pair(min_distance[source], source));
    
      while (!vertex_queue.empty()) {
        
        const Weight_t dist = vertex_queue.begin()->first;
        const Vertex_t u = vertex_queue.begin()->second;
        vertex_queue.erase(vertex_queue.begin());
  
        const auto& neighbors = adjacency_list[u];
        for (
          auto neighbor_iter = neighbors.begin();
          neighbor_iter != neighbors.end();
          neighbor_iter++) {
          
          const Vertex_t v = neighbor_iter->target;
          const Weight_t weight = neighbor_iter->weight;
          const Weight_t distance_through_u = dist + weight;
        
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
  };

} // internal
} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_INTERNAL_DIJKSTRA_H
