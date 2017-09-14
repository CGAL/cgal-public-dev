#ifndef CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_0_H
#define CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_0_H

// STL includes.
#include <vector>
#include <iostream>

// CGAL includes.
#include <CGAL/number_utils.h>
#include <CGAL/utils_classes.h>
#include <CGAL/utils.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class CDTInput, class Visibility, class StructuredLabelsInput>
		class Level_of_detail_reconstruction_0 {

		public:
			enum class Edge_label { ST, FF, IN, DEL };

			typedef KernelTraits Kernel;
			typedef CDTInput     CDT;
			typedef Visibility   Visibility_result;

			typedef StructuredLabelsInput Structured_labels;
			typedef typename Kernel::FT FT;

			typedef typename Kernel::Segment_2 Segment;
			typedef std::vector<Segment> Lod_0_result;

			using Log = CGAL::LOD::Mylog;

			Level_of_detail_reconstruction_0() : m_alpha(FT(1)), m_beta(FT(100000)), m_gamma(FT(1000)) { }

			void set_alpha_parameter(const FT alpha) {
				m_alpha = alpha;
			}

			void set_beta_parameter(const FT beta) {
				m_beta = beta;
			}

			void set_gamma_parameter(const FT gamma) {
				m_gamma = gamma;
			}

			void reconstruct(const CDT &cdt, const Visibility_result &visibility, const Structured_labels &str_labels, Lod_0_result &) {

				Log log;
				log.out << "Lod 0 reconstruction...\n\n" << std::endl;


				// (1) Create P_out and P_in predictions.
				create_predictions(visibility);

				log.out << "(1) Predicitions are created with parameter beta = " << m_beta << ". Results are saved in tmp/predictions" << std::endl;


				// (2) Compute weight for each edge in cdt.
				auto number_of_traversed_edges = compute_edge_weights(cdt);

				log.out << "(2) Edge weights are computed. Results are saved in tmp/edge_weights. Number of traversed edges: " << number_of_traversed_edges << std::endl;


				// (3) Compute circumsphere weight for each edge in cdt.
				number_of_traversed_edges = compute_circumsphere_weights(cdt);

				log.out << "(3) Circumsphere weights are computed. Results are saved in tmp/circumsphere_weights. Number of traversed edges: " << number_of_traversed_edges << std::endl;


				// (4) Compute edge coherence.
				number_of_traversed_edges = compute_edge_coherence(cdt, str_labels);

				log.out << "(4) All edges are labeled. Results are saved in tmp/coherence. Number of traversed edges: " << number_of_traversed_edges << std::endl;


				// (5) Compute surface quality function.
				compute_surface_quality_function();

				log.out << "(5) Surface quality function is computed. Results are saved in tmp/surface_quality." << std::endl;


				// (6) Compute graph cut cost function.
				const FT cost_value = compute_graph_cut_cost_function();

				log.out << "(6) Graph cut cost function is computed with cost value: " << cost_value << std::endl;


				// (7) Apply graph cut.
				apply_graph_cut(cost_value);

				log.out << "(7) Graph cut is applied." << std::endl;


				log.out << "\n\n...FINISHED" << std::endl;
				log.save("lod_0");
			}

			void max_flow(const CDT &, const Visibility_result &) {

				std::cout << "\nMAX FLOW\n" << std::endl;

			}

		private:
			FT m_alpha;
			FT m_beta;
			FT m_gamma;

			std::vector<FT> m_p_in;
			std::vector<FT> m_p_out;
			
			std::vector<FT> m_a; // edge length weights
			std::vector<FT> m_c; // edge circumsphere weights

			std::vector<Edge_label> m_coherence;
			std::vector<FT> m_q; // surface quality Q

			void create_predictions(const Visibility_result &visibility) {

				clear_predicitions();
				assert(m_p_in.empty() && m_p_out.empty());

				for (typename Visibility_result::const_iterator it = visibility.begin(); it != visibility.end(); ++it) {

					const int visibility_label = static_cast<int>((*it).second);
					switch(visibility_label) {

						case 0: // inside that is IN
							add_new_inside_predicition();
							break;

						case 1: // outside that is OUT
							add_new_outside_predicition();
							break;

						case 2: // no predicition that is UNKNOWN
							// What to do with cells that have no predicition?
							break;

						default: // any other stuff though should not be here ever
							break;
					}
				}

				// Remove later.
				Log log;

				log.out << "IN predicitions: " << std::endl;
				for (size_t i = 0; i < m_p_in.size(); ++i) log.out << "p_in = "  << m_p_in[i] << std::endl;
				log.out << std::endl;

				log.out << "OUT predicitions: " << std::endl;
				for (size_t i = 0; i < m_p_out.size(); ++i) log.out << "p_out = " << m_p_out[i] << std::endl;
				log.out << std::endl;				

				log.save("tmp/predictions");
			}

			void clear_predicitions() {
				m_p_in.clear();
				m_p_out.clear();
			}

			// Maybe later better to have a count of IN and OUT predicitions so here I could omit the push_back() then.
			void add_new_inside_predicition() {

				const FT predicition = compute_prediction();
				m_p_in.push_back(predicition);
			}

			void add_new_outside_predicition() {

				const FT predicition = compute_prediction();
				m_p_out.push_back(predicition);
			}

			FT compute_prediction() const {
				return m_beta * FT(1);
			}

			int compute_edge_weights(const CDT &cdt) {

				Log log;

				int number_of_traversed_edges = 0;

				const auto first = cdt.finite_edges_begin();
				const auto last  = cdt.finite_edges_end();

				const int number_of_edges = static_cast<int>(std::distance(first, last));
				clear_edge_weights(number_of_edges);

				assert(!m_a.empty());
				for (typename CDT::Finite_edges_iterator it = first; it != last; ++it, ++number_of_traversed_edges) {

					const auto segment = cdt.segment(it);
					m_a[number_of_traversed_edges] = CGAL::sqrt(segment.squared_length());

					log.out << number_of_traversed_edges << ": " << m_a[number_of_traversed_edges] << std::endl;
				}

				log.save("tmp/edge_weights");

				assert(number_of_traversed_edges == number_of_edges);
				return number_of_traversed_edges;
			}

			void clear_edge_weights(const int number_of_edges) {

				assert(number_of_edges >= 0);

				m_a.clear();
				m_a.resize(number_of_edges, FT(1));
			}

			// Compute them only for edges that are FF coherent.
			int compute_circumsphere_weights(const CDT &cdt) {

				Log log;

				int number_of_traversed_edges = 0;

				const auto first = cdt.finite_edges_begin();
				const auto last  = cdt.finite_edges_end();

				const int number_of_edges = static_cast<int>(std::distance(first, last));
				clear_circumsphere_weights(number_of_edges);

				assert(!m_c.empty());
				for (typename CDT::Finite_edges_iterator it = first; it != last; ++it, ++number_of_traversed_edges) {

					const auto segment = cdt.segment(it);
					m_c[number_of_traversed_edges] = compute_circumsphere_weight(segment);

					log.out << number_of_traversed_edges << ": " << m_c[number_of_traversed_edges] << std::endl;
				}

				log.save("tmp/circumsphere_weights");

				assert(number_of_traversed_edges == number_of_edges);
				return number_of_traversed_edges;
			}

			void clear_circumsphere_weights(const int number_of_edges) {

				assert(number_of_edges >= 0);

				m_c.clear();
				m_c.resize(number_of_edges, FT(1));
			}

			FT compute_circumsphere_weight(const Segment &segment) const {

				const FT cos_in  = compute_circumsphere_cosine(segment);
				const FT cos_out = compute_circumsphere_cosine(segment); 

				return m_alpha - CGAL::min(cos_in, cos_out);
			}

			FT compute_circumsphere_cosine(const Segment &) const {
				return FT(0);
			}

			int compute_edge_coherence(const CDT &cdt, const Structured_labels &str_labels) {

				Log log;

				int number_of_traversed_edges = 0;

				const auto first = cdt.finite_edges_begin();
				const auto last  = cdt.finite_edges_end();

				const int number_of_edges = static_cast<int>(std::distance(first, last));
				clear_edge_coherence(number_of_edges);

				assert(!m_coherence.empty());
				for (typename CDT::Finite_edges_iterator it = first; it != last; ++it, ++number_of_traversed_edges) {

					const auto segment = cdt.segment(it);
					m_coherence[number_of_traversed_edges] = discover_edge_label(segment, str_labels);

					std::string label = "default";
					switch(m_coherence[number_of_traversed_edges]) {
						
						case Edge_label::ST:
							label = "structure";
							break;

						case Edge_label::FF:
							label = "free form";
							break;

						case Edge_label::IN:
							label = "incoherent";
							break;

						case Edge_label::DEL:
							label = "to be deleted";
							break;

					}
					log.out << number_of_traversed_edges << ": " << label << std::endl;
				}

				log.save("tmp/coherence");

				assert(number_of_traversed_edges == number_of_edges);
				return number_of_traversed_edges;
			}

			void clear_edge_coherence(const int number_of_edges) {

				assert(number_of_edges >= 0);

				m_coherence.clear();
				m_coherence.resize(number_of_edges);
			}

			Edge_label discover_edge_label(const Segment &, const Structured_labels &) const {
				return Edge_label::IN;

				// cl li co

				// st					ff						in
				// li--li				cl--cl					cl--co
				// li--co				cl--li					co--co

				// cl - clutter
				// li - linear
				// co - corner
			}

			void compute_surface_quality_function() {

				Log log;

				assert(!m_a.empty());
				assert(m_a.size() == m_c.size() && m_c.size() == m_coherence.size());

				clear_surface_quality();
				assert(m_q.size() == m_coherence.size());

				for (size_t i = 0; i < m_coherence.size(); ++i) {
					switch(m_coherence[i]) {

						case Edge_label::ST:
							m_q[i] = FT(0);
							break;

						case Edge_label::FF:
							m_q[i] = m_c[i];
							break;

						case Edge_label::IN:
							m_q[i] = m_gamma;
							break;

						case Edge_label::DEL:
							// Should I add all bbox segments to incoherent to suppress them from the solution?
							// Or better to remove them from m_q and so m_q will not have the same size as m_a?
							m_q[i] = m_gamma; 
							break;
					}
					log.out << i << ": " << m_q[i] << std::endl;
				}
				log.save("tmp/surface_quality");
			}

			void clear_surface_quality() {

				m_q.clear();
				m_q.resize(m_a.size());
			}

			FT compute_graph_cut_cost_function() {

				FT cost_value = FT(0);

				assert(!m_a.empty());
				assert(m_a.size() == m_q.size());
				assert(!m_p_in.empty() && !m_p_out.empty());

				for (size_t i = 0; i < m_a.size();     ++i) cost_value += m_a[i] * m_q[i];
				for (size_t i = 0; i < m_p_in.size();  ++i) cost_value += m_p_in[i];
				for (size_t i = 0; i < m_p_out.size(); ++i) cost_value += m_p_out[i];

				return cost_value;
			}

			void apply_graph_cut(const FT) {

			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_0_H