#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_FOR_ORDINATES_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_FOR_ORDINATES_H

#if defined(WIN32) || defined(_WIN32) 
#define PS "\\"
#else 
#define PS "/" 
#endif 

// STL includes.
#include <map>
#include <set>
#include <cmath>
#include <vector>
#include <memory>
#include <utility>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/number_utils.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// New CGAL includes.
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_debugger.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_parameters.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_max_difference.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_max_difference_test.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_regular_segment_property_map.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_regular_segment.h>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
        class Level_of_detail_segment_regularizer_for_ordinates {

        public:
            typedef KernelTraits Kernel;

            using FT      = typename Kernel::FT;
            using Point   = typename Kernel::Point_2;
            using Segment = typename Kernel::Segment_2;
            using Line    = typename Kernel::Line_2;
            using Vector  = typename Kernel::Vector_2;

            using Regular_segment  = CGAL::LOD::Level_of_detail_segment_regularizer_regular_segment<Kernel>;
            using Regular_segments = std::vector<Regular_segment *>;

            using RegularMap   = CGAL::LOD::Level_of_detail_segment_regularizer_regular_segment_property_map<Regular_segment, Segment>;
            using RegularRange = Regular_segments;
            
            using Debugger   = CGAL::LOD::Level_of_detail_segment_regularizer_debugger;
            using Parameters = CGAL::LOD::Level_of_detail_segment_regularizer_parameters<Kernel>;

            using Differences         = std::vector<FT>;
            using Max_difference      = CGAL::LOD::Level_of_detail_segment_regularizer_max_difference<Parameters>;
            using Max_difference_test = CGAL::LOD::Level_of_detail_segment_regularizer_max_difference_test<Parameters, Regular_segments>;

            using Neighbours_graph_data = CGAL::LOD::Level_of_detail_segment_regularizer_neighbours_graph_data<Kernel>;
            
            using QP_problem_data = CGAL::LOD::Level_of_detail_segment_regularizer_qp_problem_data<Kernel, Neighbours_graph_data>;
            using QP_problem      = CGAL::LOD::Level_of_detail_segment_regularizer_ooqp_problem<Kernel, QP_problem_data>;

            using Tree = CGAL::LOD::Level_of_detail_segment_regularizer_tree<Kernel, QP_problem_data>;
            
            using Parallel_segments           = typename Tree::Parallel_segments;
            using Parallel_segments_iterator  = typename Parallel_segments::iterator;
            using Parallel_segments_tree_node = typename Tree::Parallel_segments_tree_node;

            using Internal_collinear_segments_tree_node = typename Parallel_segments_tree_node::Collinear_segments_tree_node;
            using Internal_collinear_segments           = typename Parallel_segments_tree_node::Collinear_segments;
            
            using Internal_parallel_segments          = typename Parallel_segments_tree_node::Parallel_segments;
            using Internal_parallel_segments_iterator = typename Internal_parallel_segments::const_iterator;

            using Point_with_index  = std::pair<Point, size_t>;
            using Points            = std::vector<Point_with_index>;
            using Point_segment_map = std::map<size_t, Regular_segment *>;

            using VB = CGAL::Triangulation_vertex_base_with_info_2<size_t, Kernel>;
            using DS = CGAL::Triangulation_data_structure_2<VB>;
            using DT = CGAL::Delaunay_triangulation_2<Kernel, DS>;

            using Edge_iterator = typename DT::Finite_edges_iterator;
            using Edge          = typename DT::Edge;

            using Considered_potential  = std::pair<Regular_segment *, Regular_segment *>;
            using Considered_potentials = std::set<Considered_potential>;

            using Segments_to_groups = typename Tree::Segments_to_groups;
            using Groups_to_segments = typename Tree::Groups_to_segments;
            
            using List_element    = typename Tree::List_element;
            using Nodes_to_groups = std::map<Parallel_segments_tree_node *, List_element>;

            using Mus_matrix     = typename QP_problem_data::Mus_matrix;
            using Targets_matrix = typename QP_problem_data::Targets_matrix;

            using Ordinates = std::map<int, FT>;

            Level_of_detail_segment_regularizer_for_ordinates(Regular_segments &segments, Tree *tree_pointer, const Parameters &parameters) :
            m_debug(false), m_silent(false), m_use_test_difference(true), 
            m_input_segments(segments), m_tree_pointer(tree_pointer), m_parameters(parameters) { }

            void regularize() {
                if (m_input_segments.size() == 0) return;
                
                // Set max differences between translated segments.
                set_max_differences();

                // Build a graph of all spatially close segments.
                build_graph_of_neighbours();

                // Prepare all necessary data for the sparse QP solver.
                create_input_data_for_qp_solver();

                // Solve the QP problem.
                solve_qp_problem();

                // Translate all segments based on optimized differences.
                translate_segments();

                // Print debug information if the corresponding flag is on.
                print_debug_information();
            }

            void make_silent(const bool new_state) {
                m_silent = new_state;
            }

        private:
            const bool m_debug;
            bool       m_silent;
            bool       m_use_test_difference;
            
            Regular_segments &m_input_segments;
            const Parameters &m_parameters;
            
            Differences m_max_differences;
            Differences m_final_differences;

            Debugger m_debugger;
            Tree *m_tree_pointer;

            Neighbours_graph_data m_neighbours_graph_data;
            QP_problem_data       m_qp_problem_data;

            void set_max_differences() {
                const size_t num_input_segments = m_input_segments.size();
                assert(num_input_segments > 0);

                Max_difference max_difference(m_parameters);
                Max_difference_test max_difference_test(m_parameters, m_input_segments);

                m_max_differences.clear();
                m_max_differences.resize(num_input_segments);

                if (m_use_test_difference) for (size_t i = 0; i < num_input_segments; ++i) m_max_differences[i] = max_difference_test.get(i);
                else for (size_t i = 0; i < num_input_segments; ++i) m_max_differences[i] = max_difference.get();
            }

            void build_graph_of_neighbours() {
                assert(m_input_segments.size() > 0);
                assert(m_max_differences.size() == m_input_segments.size());

                m_neighbours_graph_data.clear();

                Parallel_segments &parallel_segments = m_tree_pointer->get_parallel_segments();
                Parallel_segments_iterator it_c      = parallel_segments.begin();

                while (it_c != parallel_segments.end()) {
                    Parallel_segments_tree_node &tree_node = it_c->second;
                    
                    tree_node.delete_collinear_segments();
                    if (tree_node.get_parallel_segments().size() == 0) {

                        it_c = parallel_segments.erase(it_c);
                        continue;
                    }

                    // Transform coordinates to a rotated frame.
                    const FT angle = it_c->first;
                    transform_coordinates(angle, tree_node);

                    build_local_graph_of_neighbours(tree_node);
                    ++it_c;
                }
            }

            void transform_coordinates(const FT angle, Parallel_segments_tree_node &tree_node) {

                const Internal_parallel_segments &parallel_segments = tree_node.get_parallel_segments();
                const Point frame_origin = parallel_segments.front()->get_barycentre();

                // We are going to compute the coordinates of the segments in the frame defined by
                // the centre of first segment and the cluster.
                for (Internal_parallel_segments_iterator it_s = parallel_segments.begin(); it_s != parallel_segments.end(); ++it_s) {

                    Regular_segment *segment = (*it_s);
                    const Point &barycentre = segment->get_barycentre();

                    const FT cos_val = static_cast<FT>(cos(CGAL_PI * CGAL::to_double(angle) / 180.0));
                    const FT sin_val = static_cast<FT>(sin(CGAL_PI * CGAL::to_double(angle) / 180.0));

                    const FT x = (barycentre.x() - frame_origin.x()) * cos_val + (barycentre.y() - frame_origin.y()) * sin_val;
                    const FT y = (barycentre.y() - frame_origin.y()) * cos_val - (barycentre.x() - frame_origin.x()) * sin_val;

                    segment->set_reference_coordinates(Point(x, y));
                }
            }

            void build_local_graph_of_neighbours(const Parallel_segments_tree_node &tree_node) {

                Points            points;
                Point_segment_map points_to_segments;

                size_t j = 0;
                const Internal_parallel_segments &parallel_segments = tree_node.get_parallel_segments();
                
                for (Internal_parallel_segments_iterator it_s = parallel_segments.begin(); it_s != parallel_segments.end(); ++it_s) {
                    Regular_segment* segment = (*it_s);
                    
                    const Point &barycentre = segment->get_barycentre();
                    points.push_back(std::make_pair(Point(barycentre.x(), barycentre.y()), j));
                    points_to_segments[j] = segment;
                    ++j;
                }

                DT dt;
	            dt.insert(points.begin(), points.end());

                Considered_potentials considered_potentials;
                for (Edge_iterator it_e = dt.finite_edges_begin(); it_e != dt.finite_edges_end(); ++it_e) {

                    const Edge &edge = *it_e;

                    const size_t e_i = edge.first->vertex((edge.second + 1) % 3)->info();
                    const size_t e_j = edge.first->vertex((edge.second + 2) % 3)->info();

                    Regular_segment* s_i = points_to_segments[e_i];
                    Regular_segment* s_j = points_to_segments[e_j];

                    const size_t i = s_i->get_index();
                    const size_t j = s_j->get_index();

                    if (i == j) continue;

                    Considered_potential p_ij = (i < j ? std::make_pair(s_i, s_j) : std::make_pair(s_j, s_i));
                    if (considered_potentials.find(p_ij) != considered_potentials.end()) continue;
                    considered_potentials.insert(p_ij);

                    const FT mu_ij = m_parameters.get_lambda();
                    const FT y_ij  = s_i->get_reference_coordinates().y() - s_j->get_reference_coordinates().y();
                    
                    if (CGAL::abs(y_ij) < m_max_differences[i] + m_max_differences[j]) {

                        m_neighbours_graph_data.get_mus().push_back(Triplet<FT>(i, j, mu_ij));
                        m_neighbours_graph_data.get_targets().push_back(Triplet<FT>(i, j, y_ij));
                    }
                }
            }

            void create_input_data_for_qp_solver() {

                assert(m_neighbours_graph_data.filled());
                const size_t num_input_segments = m_input_segments.size();
                
                m_qp_problem_data.set_from(m_neighbours_graph_data, num_input_segments);
            }

            void solve_qp_problem() {
                assert(m_qp_problem_data.filled());

                QP_problem qp_problem(m_max_differences, m_qp_problem_data, m_parameters, m_input_segments);
                qp_problem.solve(m_final_differences);
            }

            void translate_segments() {

                assert(m_input_segments.size() > 0);
                assert(m_final_differences.size() >= m_input_segments.size()); 
                assert(m_qp_problem_data.filled());

                build_regularization_tree();
                Parallel_segments &parallel_segments = m_tree_pointer->get_parallel_segments();

                for (Parallel_segments_iterator it_c = parallel_segments.begin(); it_c != parallel_segments.end(); ++it_c) {
                    Parallel_segments_tree_node &tree_node = it_c->second;

                    translate(tree_node);
                    tree_node.delete_parallel_segments();
                }
            }

            void build_regularization_tree() {
                
                const int n = static_cast<int>(m_input_segments.size());

                Segments_to_groups segments_to_groups(n, -1);
                Groups_to_segments groups_to_segments;
                Nodes_to_groups nodes_to_groups;
                
                int g = 0;
                int p = 0;

                for (int k = 0; k < m_qp_problem_data.get_targets_matrix().outerSize(); ++k) {
                    for (typename Targets_matrix::InnerIterator it_tar(m_qp_problem_data.get_targets_matrix(), k); it_tar; ++it_tar) {

                        int i = it_tar.row(), j = it_tar.col();

                        const FT eps = FT(1) / FT(1000000);
                        if (CGAL::abs(m_final_differences[n + p]) < eps) {

                            // Then segments i and j belong to the same group of parallel segments.
                            // For the moment, these groups are materialized by integers.
                            if (segments_to_groups[i] == -1 && segments_to_groups[j] == -1) {
                                
                                // Segments i and j are not assigned to any group of parallel segments
                                // So we create one with them.
                                segments_to_groups[i] = segments_to_groups[j] = g;

                                groups_to_segments[g].push_back(i);
                                groups_to_segments[g].push_back(j);

                                nodes_to_groups[m_input_segments[i]->parallel_node].push_back(g);
                                g++;

                            } else if (segments_to_groups[i] == -1 && segments_to_groups[j] != -1) {

                                // Assigns segment i to the group of the segment j.
                                const int g_j = segments_to_groups[j];

                                segments_to_groups[i] = g_j;
                                groups_to_segments[g_j].push_back(i);

                            } else if (segments_to_groups[i] != -1 && segments_to_groups[j] == -1) {

                                // Assigns segment j to the group of the segment i.
                                const int g_i = segments_to_groups[i];

                                segments_to_groups[j] = g_i;
                                groups_to_segments[g_i].push_back(j);

                            } else {
                                
                                const int g_i = segments_to_groups[i];
                                const int g_j = segments_to_groups[j];

                                if (g_i != g_j) {
                                    
                                    // Segments i and j have been assigned to different groups, but in fact
                                    // they belong to the same group. That's why we merge them.
                                    for (typename List_element::iterator it_l = groups_to_segments[g_j].begin(); it_l != groups_to_segments[g_j].end(); ++it_l) {
                                        
                                        segments_to_groups[*it_l] = g_i;
                                        groups_to_segments[g_i].push_back(*it_l);
                                    }

                                    groups_to_segments[g_j].clear();

                                    // Delete entry g_j from 'nodes_to_groups'.
                                    typename List_element::iterator it_n = nodes_to_groups[m_input_segments[i]->parallel_node].begin();
                                    while (it_n != nodes_to_groups[m_input_segments[i]->parallel_node].end()) {

                                        if ((*it_n) == g_j) {

                                            nodes_to_groups[m_input_segments[i]->parallel_node].erase(it_n);
                                            break;
                                        }
                                        ++it_n;
                                    }
                                }
                            }
                        } ++p;
                    }
                }

                // We prepare the construction of the regularization tree.
                const FT y_eps = 1;
                Ordinates ordinates;

                for (size_t i = 0; i < segments_to_groups.size(); ++i) {
                    const int g_i = segments_to_groups[i];

                    Parallel_segments_tree_node *node_i = m_input_segments[i]->parallel_node;
                    if (g_i != -1) {

                        if (ordinates.find(g_i) == ordinates.end()) {
                            const FT y = m_input_segments[i]->get_reference_coordinates().y() + m_final_differences[i];

                            // Check if this ordinate seems to be associated to another group of segments.
                            int g_j = -1;
                            for (typename Ordinates::iterator it_m = ordinates.begin(); it_m != ordinates.end(); ++it_m) {
                                if (CGAL::abs(it_m->second - y) < y_eps) {

                                    // We found a value close to it_m, but does it correspond to the same group of parallel segments?
                                    typename List_element::iterator it_n = nodes_to_groups[node_i].begin();
                                    while (it_n != nodes_to_groups[node_i].end()) {
                                        
                                        if ((*it_n) == it_m->first) break;
                                        ++it_n;
                                    }
                                    if (it_n != nodes_to_groups[node_i].end()) g_j = it_m->first;
                                }
                            }

                            if (g_j == -1) ordinates[g_i] = y;
                            else {
                                
                                // Merge groups.
                                for (typename List_element::iterator it = groups_to_segments[g_i].begin(); it != groups_to_segments[g_i].end(); ++it) {
                                    
                                    segments_to_groups[*it] = g_j;
                                    groups_to_segments[g_j].push_back(*it);
                                }
                                groups_to_segments[g_i].clear();
                            }
                        }
                    }
                }

                // Try to assign segments whose orientation has not been optimized thanks to the regularization process, to an existing group.
                for (size_t i = 0; i < segments_to_groups.size(); ++i) {
                    int g_i = segments_to_groups[i];

                    Parallel_segments_tree_node* node_i = m_input_segments[i]->parallel_node;
                    if (g_i == -1) {

                        const FT y = m_input_segments[i]->get_reference_coordinates().y();

                        int g_j = -1;
                        for (typename Ordinates::iterator it_m = ordinates.begin(); it_m != ordinates.end(); ++it_m) {
                            const FT y_j = it_m->second;

                            if (CGAL::abs(y_j - y) < y_eps) {
                                
                                // We found a value close to it_m, but does it correspond to the same group of parallel segments?
                                typename List_element::iterator it_n = nodes_to_groups[node_i].begin();
                                while (it_n != nodes_to_groups[node_i].end()) {

                                    if ((*it_n) == it_m->first) break;
                                    ++it_n;
                                }
                                if (it_n != nodes_to_groups[node_i].end()) g_j = it_m->first;
                            }
                            if (g_j != -1) break;
                        }

                        if (g_j == -1) {
                            g_i = ordinates.rbegin()->first + 1;

                            ordinates[g_i] = y;
                            nodes_to_groups[node_i].push_back(g_i);

                        } else g_i = g_j;

                        segments_to_groups[i] = g_i;
                        groups_to_segments[g_i].push_back(i);
                    }
                }

                // Finally build the regularization tree.
                for (typename Nodes_to_groups::iterator it_m = nodes_to_groups.begin(); it_m != nodes_to_groups.end(); ++it_m) {
                    Parallel_segments_tree_node* node = it_m->first;

                    List_element &groups = it_m->second;
                    for (typename List_element::iterator it_n = groups.begin(); it_n != groups.end(); ++it_n) {

                        const FT y = ordinates[*it_n];
                        node->create_collinear_node(y);
                    }
                }

                // Assign segments.
                for (size_t i = 0; i < segments_to_groups.size(); ++i) {

                    const int g_i = segments_to_groups[i];
                    m_input_segments[i]->parallel_node->assign_to_collinear_node(ordinates[g_i], m_input_segments[i]);
                }
            }

            void translate(Parallel_segments_tree_node &tree_node) {

                typename Internal_collinear_segments::iterator it_m = tree_node.get_collinear_segments().begin();
                while (it_m != tree_node.get_collinear_segments().end()) {

                    const FT dt = it_m->first;
                    Internal_collinear_segments_tree_node* tree_subnode = it_m->second;

                    if (tree_subnode->get_collinear_segments().empty()) {

                        it_m = tree_node.get_collinear_segments().erase(it_m);
                        continue;
                    }

                    // Get the longest segment.
                    FT l_max = -FT(1000000000000);
                    Regular_segment* s_longest = NULL;

                    using Local_collinear_segments = std::list<Regular_segment *>;
                    for (typename Local_collinear_segments::iterator it_s = tree_subnode->get_collinear_segments().begin(); it_s != tree_subnode->get_collinear_segments().end(); ++it_s) {
                        if ((*it_s)->get_length() > l_max) {

                            l_max = (*it_s)->get_length();
                            s_longest = (*it_s);
                        }
                    }

                    // Translate the longest segment and get the line equation.
                    s_longest->set_difference(dt - s_longest->get_reference_coordinates().y());
                    
                    const FT a = s_longest->get_a();
                    const FT b = s_longest->get_b();
                    const FT c = s_longest->get_c();

                    const Vector &direction = s_longest->get_direction();

                    // Translate the other segments, so that they rest upon the line ax + by + c = 0.
                    for (typename Local_collinear_segments::iterator it_s = tree_subnode->get_collinear_segments().begin(); it_s != tree_subnode->get_collinear_segments().end(); ++it_s) {
                        if ((*it_s) != s_longest)
                            (*it_s)->set_difference(dt - (*it_s)->get_reference_coordinates().y(), a, b, c, direction);
                    }
                    ++it_m;
                }
            }

            void print_debug_information() {

                if (!m_debug) return;
                m_debugger.print_values(m_max_differences, "differences threshold in meters");
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_FOR_ORDINATES_H