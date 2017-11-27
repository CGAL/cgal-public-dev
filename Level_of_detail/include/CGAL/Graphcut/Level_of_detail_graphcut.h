#ifndef CGAL_LEVEL_OF_DETAIL_GRAPHCUT_H
#define CGAL_LEVEL_OF_DETAIL_GRAPHCUT_H

// STL includes.
#include <vector>
#include <iostream>
#include <cassert>
#include <string>
#include <map>

// CGAL includes.
#include <CGAL/number_utils.h>
#include <CGAL/utils_classes.h>
#include <CGAL/utils.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/IO/Color.h>

namespace Maxflow {
	#include <CGAL/internal/auxiliary/graph.h>
}

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class CDTInput>
		class Level_of_detail_graphcut {

		public:
			typedef KernelTraits Kernel;
			typedef CDTInput     CDT;

			typedef typename Kernel::FT 	   FT;
			typedef typename Kernel::Point_2   Point_2;
			typedef typename Kernel::Segment_2 Edge;
			typedef typename Kernel::Vector_2  Vector_2;

			typedef Maxflow::Graph 			Graph;
			typedef typename Graph::node_id Node_id;

			typedef typename CDT::Face_handle 			Face_handle;
			typedef typename CDT::Finite_faces_iterator Face_iterator;
			typedef typename CDT::Finite_edges_iterator Edge_iterator;

			using Nodes_map = std::map<Face_handle, int>;
			using Log 		= CGAL::LOD::Mylog;

			enum class Node_value_adapter { NO_ADAPTER, TANH };

			Level_of_detail_graphcut() : 
			m_alpha(FT(1)), 
			m_beta(FT(100000)), 
			m_gamma(FT(1000)), 
			m_save_info(true),
			m_node_value_adapter(Node_value_adapter::NO_ADAPTER) { 

				set_coherence_dictionary();
			}

			void set_alpha_parameter(const FT alpha) {
				m_alpha = alpha;
			}

			void set_beta_parameter(const FT beta) {
				m_beta = beta;
			}

			void set_gamma_parameter(const FT gamma) {
				m_gamma = gamma;
			}

			void save_info(const bool new_state) {
				m_save_info = new_state;
			}

			// Source - inside, sink - outside
			void max_flow(CDT &cdt)  {

				Nodes_map nodes_map;
				// const int num_nodes = ;

				unsigned int nb_faces = cdt.number_of_faces();
				Node_id *pNodes = new Node_id[nb_faces];
				Graph *graph = new Graph();

				set_graph_nodes(cdt, pNodes, nodes_map, graph);
				set_graph_edges(cdt, pNodes, nodes_map, graph);

				graph->maxflow();

				extract_solution(pNodes, graph, cdt);

				// Remove later.
				Log log;
				log.save_visibility_eps(cdt, "tmp/after_cut");

				// cleanup
				delete graph;
			}

		private:
			enum class Edge_coherence { STRUCTURE_COHERENT, FREE_FORM_COHERENT, INCOHERENT };

			FT m_alpha;
			FT m_beta;
			FT m_gamma;

			std::map< std::pair<Structured_label, Structured_label>, Edge_coherence> m_coherence_dict;
			bool m_save_info;

			const Node_value_adapter m_node_value_adapter;

			void set_coherence_dictionary() {

				m_coherence_dict.clear();

				m_coherence_dict[std::make_pair(Structured_label::LINEAR, Structured_label::LINEAR)] = Edge_coherence::STRUCTURE_COHERENT;
				m_coherence_dict[std::make_pair(Structured_label::LINEAR, Structured_label::CORNER)] = Edge_coherence::STRUCTURE_COHERENT;
				m_coherence_dict[std::make_pair(Structured_label::CORNER, Structured_label::LINEAR)] = Edge_coherence::STRUCTURE_COHERENT;

				m_coherence_dict[std::make_pair(Structured_label::CLUTTER, Structured_label::CLUTTER)] = Edge_coherence::FREE_FORM_COHERENT;
				m_coherence_dict[std::make_pair(Structured_label::CLUTTER, Structured_label::LINEAR)]  = Edge_coherence::FREE_FORM_COHERENT;
				m_coherence_dict[std::make_pair(Structured_label::LINEAR , Structured_label::CLUTTER)] = Edge_coherence::FREE_FORM_COHERENT;

				m_coherence_dict[std::make_pair(Structured_label::CORNER , Structured_label::CORNER)]  = Edge_coherence::INCOHERENT;
				m_coherence_dict[std::make_pair(Structured_label::CLUTTER, Structured_label::CORNER)]  = Edge_coherence::INCOHERENT;
				m_coherence_dict[std::make_pair(Structured_label::CORNER , Structured_label::CLUTTER)] = Edge_coherence::INCOHERENT;
			}

			// Favour bigger cost value!
			void set_graph_nodes(const CDT &cdt, Node_id *pNodes, Nodes_map &nodes_map, Graph *graph) const {

				Log log;

				int index = 0;
				for (Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit, ++index) {
					log.add_index(index);

					pNodes[index]   = graph->add_node();
					nodes_map[fit] = index;

					const FT in  = fit->info().in;
					const FT out = FT(1) - in;

					assert(in >= FT(0) && in <= FT(1));

					const FT cost_in  = get_graph_node_cost(in);  // here we add bigger value for the desirable outcome
					const FT cost_out = get_graph_node_cost(out); // e.g. if we want in (or SOURCE) then the cost_in with bigger value is favourable

					graph->add_tweights(pNodes[index], cost_in, cost_out);

					if (m_save_info) {
						log.out << "in: " << cost_in << "; out: " << cost_out << std::endl;
						log.skip_line();
					}
				}
				if (m_save_info) log.save("tmp/graph_nodes");
			}

			FT get_graph_node_cost(const FT node_value) const {
				return get_graph_node_weight() * get_final_node_value(node_value);
			} 

			FT get_final_node_value(const FT node_value) const {

				switch (m_node_value_adapter) {

					case Node_value_adapter::NO_ADAPTER:
						return node_value;
						break;

					case Node_value_adapter::TANH:
						return tanh_adapter(node_value);
						break;

					default:
						assert(!"Wrong node value adapter!");
						return FT(0.0);
				}
			}

			FT tanh_adapter(const FT value) const {
				return FT(1) / FT(2) + tanh(value * FT(8) - FT(4)) / FT(2);
			}

			FT get_graph_node_weight() const {
				return m_beta;
			}

			// Favour lower cost value!
			void set_graph_edges(const CDT &cdt, const Node_id nodes[], const Nodes_map &nodes_map, Graph *graph) const {

				Log log;

				int index = 0;
				for (Edge_iterator eit = cdt.finite_edges_begin(); eit != cdt.finite_edges_end(); ++eit, ++index) {

					// Skip all boundary edges.				
					if (is_boundary_edge(cdt, eit)) continue;


					// Output current edge.
					if (m_save_info) log.add_index(index);
					Edge edge = cdt.segment(eit);
					if (m_save_info) log.out << "Edge: " << edge << " with: " << std::endl;


					// Compute edge weight.
					const FT edge_weight = compute_edge_weight(edge);
					if (m_save_info) log.out << "weight = " << edge_weight << std::endl;


					// Find edge coherence.
					const Edge_coherence edge_coherence = find_edge_coherence(cdt, eit);

					std::string str_coherence = "default"; 
					switch(edge_coherence) {
						case Edge_coherence::STRUCTURE_COHERENT:
							str_coherence = "strucutre coherent";
							break;

						case Edge_coherence::FREE_FORM_COHERENT:
							str_coherence = "free form coherent";
							break;

						case Edge_coherence::INCOHERENT:
							str_coherence = "incoherent";
							break;

						default:
							assert(!"Wrong edge coherence label!");
							break;
					}
					if (m_save_info) log.out << "coherence = " << str_coherence << std::endl;


					// Compute edge quality.
					const FT edge_quality = compute_edge_quality(edge_coherence, cdt, eit);
					if (m_save_info) log.out << "quality penalizer = " << edge_quality << std::endl;


					// Set graph edge.
					add_graph_edge(eit, nodes, nodes_map, edge_weight, edge_quality, graph, log);


					if (m_save_info) log.skip_line();
 				}
 				if (m_save_info) log.save("tmp/graph_edges");

 				/* // Test data.
				graph->add_edge(nodes[0] , nodes[7] , 0.0, 0.0); // for min cut we add big weights to the incorrect cells
				graph->add_edge(nodes[1] , nodes[0] , 0.0, 0.0);
				graph->add_edge(nodes[2] , nodes[1] , 0.0, 0.0);
				graph->add_edge(nodes[6] , nodes[1] , 0.0, 0.0);
				graph->add_edge(nodes[6] , nodes[5] , 0.0, 0.0);
				graph->add_edge(nodes[3] , nodes[6] , 0.0, 0.0);
				graph->add_edge(nodes[2] , nodes[11], 0.0, 0.0);
				graph->add_edge(nodes[11], nodes[13], 0.0, 0.0);
				graph->add_edge(nodes[7] , nodes[4] , 0.0, 0.0);
				graph->add_edge(nodes[11], nodes[4] , 0.0, 0.0);
				graph->add_edge(nodes[4] , nodes[14], 0.0, 0.0);
				graph->add_edge(nodes[13], nodes[12], 0.0, 0.0);
				graph->add_edge(nodes[12], nodes[15], 0.0, 0.0);
				graph->add_edge(nodes[14], nodes[15], 0.0, 0.0);
				graph->add_edge(nodes[15], nodes[8] , 0.0, 0.0);
				graph->add_edge(nodes[8] , nodes[9] , 0.0, 0.0);
				graph->add_edge(nodes[8] , nodes[10], 0.0, 0.0); */
			}

			void add_graph_edge(
				const Edge_iterator &edge_handle, 
				const Node_id nodes[], 
				const Nodes_map &nodes_map, 
				const FT edge_weight, 
				const FT edge_quality, Graph *graph, Log &log) const {

				const int vertex_index = edge_handle->second;

				const Face_handle face_1 = edge_handle->first;
				const Face_handle face_2 = face_1->neighbor(vertex_index);

				const FT cost_value = get_graph_edge_cost(edge_weight, edge_quality);

				const int index_1 = nodes_map.at(face_1);
				const int index_2 = nodes_map.at(face_2);

				graph->add_edge(nodes[index_1], nodes[index_2], cost_value, cost_value);

				if (m_save_info) log.out << "Final edge: " << cost_value << " < === > " << cost_value << std::endl;
			}

			FT get_graph_edge_cost(const FT edge_weight, const FT edge_quality) const {
				return edge_weight * edge_quality;
			}

			bool is_boundary_edge(const CDT &cdt, const Edge_iterator &edge_handle) const {

				const int vertex_index = edge_handle->second;

				const Face_handle face_1 = edge_handle->first;
				const Face_handle face_2 = face_1->neighbor(vertex_index);

				if (cdt.is_infinite(face_1) || cdt.is_infinite(face_2)) return true;
				return false;
			}

			FT compute_edge_weight(const Edge &edge) const {
				return CGAL::sqrt(edge.squared_length());
			}

			Edge_coherence find_edge_coherence(const CDT &cdt, const Edge_iterator &edge_handle) const {
				
				const Face_handle face = edge_handle->first;
				const int vertex_index = edge_handle->second;

				const Structured_label label_source = face->vertex(cdt.ccw(vertex_index))->info().label;
				const Structured_label label_target = face->vertex( cdt.cw(vertex_index))->info().label;

				Edge_coherence edge_coherence = find_coherence_in_dictionary(label_source, label_target);
				if (edge_coherence == Edge_coherence::STRUCTURE_COHERENT && !cdt.is_constrained(*edge_handle)) edge_coherence = Edge_coherence::INCOHERENT;

				return edge_coherence;
			}

			Edge_coherence find_coherence_in_dictionary(Structured_label label_source, Structured_label label_target) const {

				if (!is_valid_label(label_source)) set_invalid_label(label_source); 
				if (!is_valid_label(label_target)) set_invalid_label(label_target);

				return m_coherence_dict.at(std::make_pair(label_source, label_target));
			}

			bool is_valid_label(const Structured_label label) const {

				switch (label) {

					case Structured_label::LINEAR:
						return true;

					case Structured_label::CORNER:
						return true;

					case Structured_label::CLUTTER:
						return true;

					default:
						assert(!"Wrong vertex label!");
						return false;
				}
				return false;
			}

			void set_invalid_label(Structured_label &label) const {

				label = Structured_label::CLUTTER; // maybe I can be smarter than CLUTTER here and below!
			}

			FT compute_edge_quality(const Edge_coherence edge_coherence, const CDT &cdt, const Edge_iterator &edge_handle) const {

				switch(edge_coherence) {

					case Edge_coherence::STRUCTURE_COHERENT:
						return FT(0);

					case Edge_coherence::FREE_FORM_COHERENT:
						return compute_free_form_quality(cdt, edge_handle);

					case Edge_coherence::INCOHERENT:
						return m_gamma;

					default:
						assert(!"Wrong edge coherence!");
						return FT(0.0);
				}
			}

			// Is it ok if cos_in and/or cos_out are negative?
			FT compute_free_form_quality(const CDT &cdt, const Edge_iterator &edge_handle) const {

				const auto mirror_edge = cdt.mirror_edge(*edge_handle);

				const int vertex_index_1 = edge_handle->second;
				const int vertex_index_2 = mirror_edge.second;

				const Face_handle face_1 = edge_handle->first;
				const Face_handle face_2 = face_1->neighbor(vertex_index_1);

				assert(face_2->neighbor(vertex_index_2) == face_1);
				assert(!cdt.is_infinite(face_1) && !cdt.is_infinite(face_2));

				const Edge edge = cdt.segment(edge_handle);

				/* // test data
				const Edge edge = Edge(Point_2(2.72, 5.10), Point_2(7.02, 4.54));
				const Point_2 c = Point_2(5.30, 8.14);
				const Point_2 d = Point_2(5.14, 2.26); */

				const Point_2 &a = edge.source();
				const Point_2 &b = edge.target();
				const Point_2 &c = cdt.triangle(face_1).vertex(vertex_index_1);
				const Point_2 &d = cdt.triangle(face_2).vertex(vertex_index_2);

				const FT cos_in  = compute_cos_value(a, b, c, edge, false);
				const FT cos_out = compute_cos_value(a, d, b, edge, true);

				/*
				std::cout << a << std::endl;
				std::cout << b << std::endl;
				std::cout << c << std::endl;
				std::cout << d << std::endl;

				std::cout << cos_in << " " << cos_out << std::endl; */

				const FT eps = FT(1) / FT(1000000000);

				assert(CGAL::abs(cos_in)  >= FT(0) && (FT(1) - CGAL::abs(cos_in))  > -eps);
				assert(CGAL::abs(cos_out) >= FT(0) && (FT(1) - CGAL::abs(cos_out)) > -eps);

				const FT min_cos = CGAL::min(cos_in, cos_out);
				const FT result = m_alpha - min_cos;

				assert(result >= FT(0));

				return result;
			}

			FT compute_cos_value(const Point_2 &a, const Point_2 &b, const Point_2 &c, const Edge &edge, const bool rotate) const {

				const Vector_2 edge_vector  = Vector_2(edge.source(), edge.target());
				Vector_2 edge_normal        = edge_vector.perpendicular(rotate ? CGAL::CLOCKWISE : CGAL::COUNTERCLOCKWISE);
				edge_normal 		       /= CGAL::sqrt(edge_normal * edge_normal);

				const Point_2 circumcentre      = CGAL::circumcenter(a, b, c);
				Vector_2 circle_tangent_normal  = Vector_2(a, circumcentre);
				circle_tangent_normal          /= CGAL::sqrt(circle_tangent_normal * circle_tangent_normal);

				return rotate ? (edge_normal * circle_tangent_normal) : (circle_tangent_normal * edge_normal);
			}

			void extract_solution(const Node_id *pNodes, Graph *graph, CDT &cdt) const {

				int index = 0;
				for (Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit, ++index) {

					if (graph->what_segment(pNodes[index]) == Graph::SOURCE) fit->info().in = FT(1);
					if (graph->what_segment(pNodes[index]) == Graph::SINK)   fit->info().in = FT(0);

					fit->info().in_color = get_color(fit->info().in);
				}
			}

			CGAL::Color get_color(const FT visibility) const {

				const FT half = FT(1) / FT(2);

				if (visibility > half)      return CGAL::Color(51, 255, 51); // INSIDE
				else if (visibility < half) return CGAL::Color(255, 51, 51); // OUTSIDE
									  
				assert(!"Cannot be here!");
				return CGAL::Color(255, 204, 0); // UNKNOWN
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_GRAPHCUT_H