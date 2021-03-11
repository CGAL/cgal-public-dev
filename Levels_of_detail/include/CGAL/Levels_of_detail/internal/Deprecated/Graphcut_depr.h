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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_GRAPHCUT_DEPR_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_GRAPHCUT_DEPR_H

// CGAL includes.
#include <CGAL/assertions.h>

/// \cond SKIP_IN_MANUAL
namespace Maxflow {
	#include <CGAL/internal/auxiliary/graph.h>
}
/// \endcond

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename Partition_d>
  class Graphcut_depr {

  public:
    using Traits = GeomTraits;
    using Partition = Partition_d;

    using FT = typename Traits::FT;

    using Face = typename Partition::Face;
    using Edge = typename Partition::Edge;

		using Graph = Maxflow::Graph;
		using Node_id = typename Graph::node_id;

		struct Graph_wrapper {

			Node_id* pNodes;
			Graph* graph;

			Graph_wrapper(const unsigned int num_nodes) {
				// +1 gives an extra infinite node
				pNodes = new Node_id[num_nodes + 1];
				graph = new Graph();
			}
			~Graph_wrapper() {
				delete graph;
				delete[] pNodes;
			}
			void add_node(const std::size_t idx) {
				pNodes[idx] = graph->add_node();
			}
			void add_node_weights(
				const std::size_t idx, const FT cost_in, const FT cost_out) {
				graph->add_tweights(
					pNodes[idx], CGAL::to_double(cost_in), CGAL::to_double(cost_out));
			}
			void add_edge(
				const std::size_t idxi, const std::size_t idxj, const FT cost_value) {
				graph->add_edge(pNodes[idxi], pNodes[idxj],
        	CGAL::to_double(cost_value),
        	CGAL::to_double(cost_value));
			}
			bool is_source(const std::size_t idx) const {
				return graph->what_segment(pNodes[idx]) == Graph::SOURCE;
			}
			void maxflow() {
				graph->maxflow();
			}
		};

    Graphcut_depr(const FT graphcut_beta) :
    m_beta(graphcut_beta)
    { }

    void apply(Partition& partition) const {

      if (partition.empty()) return;
      auto& faces = partition.faces;
      auto& edges = partition.edges;

			compute_weights(faces);
			compute_weights(edges);

			const unsigned int num_nodes =
			static_cast<unsigned int>(faces.size());
			Graph_wrapper graph(num_nodes);
			const unsigned int inf_node_index = num_nodes;

			set_graph_nodes(faces, graph);
			set_graph_edges(inf_node_index, edges, graph);

			graph.maxflow();
			set_solution(graph, faces);
    }

  private:
    const FT m_beta;

		template<typename Object>
		void compute_weights(
			std::vector<Object>& objects) const {

			FT sum = FT(0);
			for (auto& object : objects) {
				object.compute_weight();
				sum += object.weight;
			}
			CGAL_assertion(sum > FT(0));
			for (auto& object : objects)
				object.weight /= sum;
		}

    void set_graph_nodes(
      const std::vector<Face>& faces,
      Graph_wrapper& graph) const {

			for (std::size_t i = 0; i < faces.size(); ++i) {
				const auto& face = faces[i];

				const FT in = face.inside;
				const FT out = face.outside;

				CGAL_assertion(in >= FT(0) && in <= FT(1));
				CGAL_assertion(out >= FT(0) && out <= FT(1));
				CGAL_assertion((in + out) == FT(1));

				const FT node_weight = face.weight;
				CGAL_precondition(node_weight >= FT(0));

				const FT cost_in = get_graph_node_cost(in, node_weight);
				const FT cost_out = get_graph_node_cost(out, node_weight);

				graph.add_node(i);
				graph.add_node_weights(i, cost_in, cost_out);
			}
			set_infinite_node(faces.size(), graph);
		}

		FT get_graph_node_cost(
      const FT node_value,
      const FT node_weight) const {

      return node_weight * node_value;
		}

		void set_infinite_node(
      const std::size_t inf_node_index,
      Graph_wrapper& graph) const {

			const FT cost_in = FT(0);
			const FT cost_out = internal::max_value<FT>();
			graph.add_node(inf_node_index);
			graph.add_node_weights(inf_node_index, cost_in, cost_out);
		}

		void set_graph_edges(
      const int inf_node_index,
      const std::vector<Edge>& edges,
      Graph_wrapper& graph) const {

			for (const auto& edge : edges) {
				const auto& neighbors = edge.neighbors;
				int idx1 = neighbors.first;
				int idx2 = neighbors.second;

				// Boundary edges.
				if (idx1 < 0 && idx2 >= 0)
					idx1 = inf_node_index;
				if (idx2 < 0 && idx1 >= 0)
					idx2 = inf_node_index;

				// Internal edges.
				const FT edge_weight = edge.weight;
				CGAL_assertion(edge_weight >= FT(0));
				add_graph_edge(idx1, idx2, edge_weight, graph);
			}
		}

		void add_graph_edge(
			const int i,
			const int j,
			const FT edge_weight,
			Graph_wrapper& graph) const {

			CGAL_assertion(i >= 0 && j >= 0);
			const FT cost_value = get_graph_edge_cost(edge_weight);
			graph.add_edge(std::size_t(i), std::size_t(j), cost_value);
		}

		FT get_graph_edge_cost(const FT edge_weight) const {
			return m_beta * edge_weight;
		}

		void set_solution(
      const Graph_wrapper& graph,
      std::vector<Face>& faces) const {

			for (std::size_t i = 0; i < faces.size(); ++i) {
				auto& face = faces[i];
				if (graph.is_source(i)) {
          face.visibility = Visibility_label::INSIDE;
          continue;
				}
				face.visibility = Visibility_label::OUTSIDE;
			}
		}
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_GRAPHCUT_DEPR_H
