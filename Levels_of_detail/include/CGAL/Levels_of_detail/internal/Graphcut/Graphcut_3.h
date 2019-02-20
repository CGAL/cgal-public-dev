#ifndef CGAL_LEVELS_OF_DETAIL_GRAPHCUT_3_H
#define CGAL_LEVELS_OF_DETAIL_GRAPHCUT_3_H

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>

namespace Maxflow {
	#include <CGAL/internal/auxiliary/graph.h>
}

// Internal includes.
#include <CGAL/Levels_of_detail/enumerations.h>
#include <CGAL/Levels_of_detail/internal/structures.h>
#include <CGAL/Levels_of_detail/internal/utilities.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Graphcut_3 {
			
  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;

    using Polyhedron_facet_3 = Polyhedron_facet_3<Traits>;
    using Graphcut_face_3 = Graphcut_face_3<Traits>;

		using Graph = Maxflow::Graph;
		using Node_id = typename Graph::node_id;

    Graphcut_3(const FT beta) :
    m_beta(beta) 
    { }

    void apply(
      const std::vector<Graphcut_face_3>& graphcut_faces,
      std::vector<Polyhedron_facet_3>& polyhedrons) const {

			const unsigned int num_nodes = polyhedrons.size();

			Node_id* pNodes = new Node_id[num_nodes + 1]; // +1 gives an extra infinite node
			Graph* graph = new Graph();

			const int inf_node_index = num_nodes;

			set_graph_nodes(polyhedrons, pNodes, graph);
			set_graph_edges(inf_node_index, graphcut_faces, pNodes, graph);

			graph->maxflow();
			set_solution(pNodes, graph, polyhedrons);

			delete graph;
			delete[] pNodes;
    }
    
  private:
    const FT m_beta;

		void set_graph_nodes(
      const std::vector<Polyhedron_facet_3>& polyhedrons, 
      Node_id* pNodes, 
      Graph* graph) const {

			for (std::size_t i = 0; i < polyhedrons.size(); ++i) {
				const auto& polyhedron = polyhedrons[i];

				const FT in = polyhedron.in;
				const FT out = polyhedron.out;

				CGAL_precondition(in >= FT(0) && in <= FT(1));
				CGAL_precondition(out >= FT(0) && out <= FT(1));

				const FT node_weight = polyhedron.weight;
				CGAL_precondition(node_weight >= FT(0));

				const FT cost_in = get_graph_node_cost(in, node_weight);
				const FT cost_out = get_graph_node_cost(out, node_weight);

				pNodes[i] = graph->add_node();
				graph->add_tweights(pNodes[i], CGAL::to_double(cost_in), CGAL::to_double(cost_out));
			}
			set_infinite_node(polyhedrons.size(), pNodes, graph);
		}

		FT get_graph_node_cost(
      const FT node_value, 
      const FT node_weight) const {
			
      return node_weight * node_value;
		} 

		void set_infinite_node(
      const std::size_t node_index, 
      Node_id* pNodes, 
      Graph* graph) const {

			const FT cost_in = FT(0);
			const FT cost_out = internal::max_value<FT>();

			pNodes[node_index] = graph->add_node();
			graph->add_tweights(
        pNodes[node_index], 
        CGAL::to_double(cost_in), 
        CGAL::to_double(cost_out));
		}


		void set_graph_edges(
      const int inf_node_index, 
      const std::vector<Graphcut_face_3>& gc_faces, 
      const Node_id* pNodes, 
      Graph* graph) const {

			for (std::size_t i = 0; i < gc_faces.size(); ++i) {
					
				const auto& gc_face = gc_faces[i];
				const auto& neighbors = gc_face.neighbors;

				const auto& neigh_1 = neighbors.first;
				const auto& neigh_2 = neighbors.second;

				int polyhedron_index_1 = neigh_1.first;
				int polyhedron_index_2 = neigh_2.first;

				// Boundary facets.
				if (polyhedron_index_1 < 0 && polyhedron_index_2 >= 0)
					polyhedron_index_1 = inf_node_index;

				if (polyhedron_index_2 < 0 && polyhedron_index_1 >= 0)
					polyhedron_index_2 = inf_node_index;

				// Internal facets.
				const FT edge_weight = gc_face.weight;
				CGAL_precondition(edge_weight >= FT(0));

				add_graph_edge(
          pNodes, 
          polyhedron_index_1, polyhedron_index_2, 
          edge_weight, 
          graph);
			}
		}

    void add_graph_edge(
			const Node_id* pNodes, 
			const int i,
			const int j,
			const FT edge_weight, 
			Graph* graph) const {

			const FT cost_value = get_graph_edge_cost(edge_weight);
			graph->add_edge(
        pNodes[i], pNodes[j], 
        CGAL::to_double(cost_value), 
        CGAL::to_double(cost_value));
		}

		FT get_graph_edge_cost(const FT edge_weight) const {
			return m_beta * edge_weight;
		}

		void set_solution(
      const Node_id* pNodes, 
      Graph* graph, 
      std::vector<Polyhedron_facet_3>& polyhedrons) const {

			for (std::size_t i = 0; i < polyhedrons.size(); ++i) {
				auto& polyhedron = polyhedrons[i];

				if (graph->what_segment(pNodes[i]) == Graph::SOURCE) { 
          
          polyhedron.visibility = Visibility_label::INSIDE;
          continue; 
        }				
        if (graph->what_segment(pNodes[i]) == Graph::SINK) { 
          
          polyhedron.visibility = Visibility_label::OUTSIDE;
          continue; 
        }
				polyhedron.visibility = Visibility_label::OUTSIDE;
			}
		}

  }; // Graphcut_3

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_GRAPHCUT_3_H
