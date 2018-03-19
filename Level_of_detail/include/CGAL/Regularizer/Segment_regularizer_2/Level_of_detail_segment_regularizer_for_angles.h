#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_FOR_ANGLES_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_FOR_ANGLES_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\"
#else 
#define PSR "/" 
#endif 

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <memory>
#include <utility>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/property_map.h>

// New CGAL includes.
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_tree.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_debugger.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_parameters.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_regular_segment.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_max_orientation.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_max_orientation_test.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_neighbours_graph_data.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_regular_segment_property_map.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_delaunay_neighbours_graph_builder.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_qp_problem_data.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_ooqp_problem.h>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
        class Level_of_detail_segment_regularizer_for_angles {

        public:
            typedef KernelTraits Kernel;

            using FT      = typename Kernel::FT;
            using Point   = typename Kernel::Point_2;
            using Segment = typename Kernel::Segment_2;
            using Line    = typename Kernel::Line_2;

            using Regular_segment  = CGAL::LOD::Level_of_detail_segment_regularizer_regular_segment<Kernel>;
            using Regular_segments = std::vector<Regular_segment *>;

            using RegularMap   = CGAL::LOD::Level_of_detail_segment_regularizer_regular_segment_property_map<Regular_segment, Segment>;
            using RegularRange = Regular_segments;
            
            using Debugger   = CGAL::LOD::Level_of_detail_segment_regularizer_debugger;
            using Parameters = CGAL::LOD::Level_of_detail_segment_regularizer_parameters<Kernel>;

            using Orientations         = std::vector<FT>;
            using Max_orientation      = CGAL::LOD::Level_of_detail_segment_regularizer_max_orientation<Parameters>;
            using Max_orientation_test = CGAL::LOD::Level_of_detail_segment_regularizer_max_orientation_test<Parameters, Regular_segments>;

            using Neighbours_graph_data    = CGAL::LOD::Level_of_detail_segment_regularizer_neighbours_graph_data<Kernel>;
            using Neighbours_graph_builder = CGAL::LOD::Level_of_detail_segment_regularizer_delaunay_neighbours_graph_builder<Kernel, Neighbours_graph_data>;

            using QP_problem_data = CGAL::LOD::Level_of_detail_segment_regularizer_qp_problem_data<Kernel, Neighbours_graph_data>;
            using QP_problem      = CGAL::LOD::Level_of_detail_segment_regularizer_ooqp_problem<Kernel, QP_problem_data>;

            using Tree = CGAL::LOD::Level_of_detail_segment_regularizer_tree<Kernel, QP_problem_data>;

            Level_of_detail_segment_regularizer_for_angles(Regular_segments &segments, const Parameters &parameters) :
            m_debug(false), m_silent(false), m_use_test_orientation(true), m_input_segments(segments), m_parameters(parameters) { }

            ~Level_of_detail_segment_regularizer_for_angles() {
                delete m_tree_pointer;
            }

            void regularize() {
                if (m_input_segments.size() == 0) return;
                
                // Set max orientations for all segments.
                set_max_orientations();

                // Build a graph of all spatially close segments.
                build_graph_of_neighbours();
                if (!m_neighbours_graph_data.filled()) return;

                // Prepare all necessary data for the sparse QP solver.
                create_input_data_for_qp_solver();

                // Solve the QP problem.
                solve_qp_problem();
                
                // Rotate all segments based on optimized orientations.
                reorient_segments();

                // Print debug information if the corresponding flag is on.
                print_debug_information();
            }

            Tree *get_tree_pointer() {
                return m_tree_pointer;
            }

            void make_silent(const bool new_state) {
                m_silent = new_state;
            }

        private:
            const bool m_debug;
            bool       m_silent;
            bool       m_use_test_orientation;
            
            Regular_segments &m_input_segments;
            const Parameters &m_parameters;
            
            Orientations m_max_orientations;
            Orientations m_final_orientations;
            
            Neighbours_graph_data m_neighbours_graph_data;
            QP_problem_data       m_qp_problem_data;

            Debugger m_debugger;
            Tree *m_tree_pointer;

            void set_max_orientations() {
                const size_t num_input_segments = m_input_segments.size();
                assert(num_input_segments > 0);

                Max_orientation max_orientation(m_parameters);
                Max_orientation_test max_orientation_test(m_parameters, m_input_segments);

                m_max_orientations.clear();
                m_max_orientations.resize(num_input_segments);

                if (m_use_test_orientation) for (size_t i = 0; i < num_input_segments; ++i) m_max_orientations[i] = max_orientation_test.get(i);
                else for (size_t i = 0; i < num_input_segments; ++i) m_max_orientations[i] = max_orientation.get();
            }

            void build_graph_of_neighbours() {
                assert(m_input_segments.size() > 0);
                assert(m_max_orientations.size() == m_input_segments.size());

                Neighbours_graph_builder neighbours_graph_builder(m_input_segments, m_max_orientations, m_parameters);

                neighbours_graph_builder.make_silent(m_silent);
                neighbours_graph_builder.build_graph_data(m_neighbours_graph_data);
            }

            void create_input_data_for_qp_solver() {
                
                assert(m_neighbours_graph_data.filled());
                const size_t num_input_segments = m_input_segments.size();
                
                m_qp_problem_data.set_from(m_neighbours_graph_data, num_input_segments);
            }

            void solve_qp_problem() {
                assert(m_qp_problem_data.filled());

                QP_problem qp_problem(m_max_orientations, m_qp_problem_data, m_parameters, m_input_segments);
                qp_problem.solve(m_final_orientations);
            }

            void reorient_segments() {

                assert(m_input_segments.size() > 0);
                assert(m_final_orientations.size() >= m_input_segments.size()); 
                assert(m_qp_problem_data.filled());

                m_tree_pointer = new Tree(m_input_segments, m_final_orientations, m_qp_problem_data, m_parameters);
                m_tree_pointer->apply_new_orientations();
            }

            void print_debug_information() {

                if (!m_debug) return;
                m_debugger.print_values(m_max_orientations, "orientations threshold in degrees");
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_FOR_ANGLES_H