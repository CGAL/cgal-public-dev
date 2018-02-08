#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_2_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_2_H

#if defined(WIN32) || defined(_WIN32) 
#define PS "\\"
#else 
#define PS "/" 
#endif 

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <utility>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/property_map.h>

// New CGAL includes.
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_tree.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_debugger.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_parameters.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_max_orientation.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_regular_segment.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_neighbours_graph_data.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_regular_segment_property_map.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_delaunay_neighbours_graph_builder.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_qp_problem_data.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_ooqp_problem.h>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
        class Level_of_detail_segment_regularizer_2 {

        public:
            typedef KernelTraits Kernel;

            using FT      = typename Kernel::FT;
            using Point   = typename Kernel::Point_2;
            using Segment = typename Kernel::Segment_2;
            using Line    = typename Kernel::Line_2;

            using Regular_segment  = CGAL::LOD::Level_of_detail_segment_regularizer_regular_segment<Kernel>;
            using Regular_segments = std::vector<Regular_segment>;

            using RegularMap   = CGAL::LOD::Level_of_detail_segment_regularizer_regular_segment_property_map<Regular_segment, Segment>;
            using RegularRange = Regular_segments;
            
            using Logger     = CGAL::LOD::Level_of_detail_segment_regularizer_debugger;
            using Debugger   = CGAL::LOD::Level_of_detail_segment_regularizer_debugger;
            using Parameters = CGAL::LOD::Level_of_detail_segment_regularizer_parameters<Kernel>;

            using Max_orientation = CGAL::LOD::Level_of_detail_segment_regularizer_max_orientation<Parameters>;
            using Orientations    = std::vector<FT>;

            using Neighbours_graph_data    = CGAL::LOD::Level_of_detail_segment_regularizer_neighbours_graph_data<Kernel>;
            using Neighbours_graph_builder = CGAL::LOD::Level_of_detail_segment_regularizer_delaunay_neighbours_graph_builder<Kernel, Neighbours_graph_data>;

            using QP_problem_data = CGAL::LOD::Level_of_detail_segment_regularizer_qp_problem_data<Kernel, Neighbours_graph_data>;
            using QP_problem      = CGAL::LOD::Level_of_detail_segment_regularizer_ooqp_problem<Kernel, QP_problem_data>;

            using Tree = CGAL::LOD::Level_of_detail_segment_regularizer_tree<Kernel, QP_problem_data>;

            Level_of_detail_segment_regularizer_2() : m_debug(false), m_silent(false), m_parameters(), m_max_orientation(m_parameters) { }

            template<typename SegmentRange, typename SegmentMap>
            void regularize(SegmentRange &input_segments, SegmentMap segment_map) {
                if (input_segments.size() == 0) return;

                // Copy input segments in the internal storage.
                copy_input_segments(input_segments, segment_map);
                
                // Set max orientations for all segments.
                set_max_orientations();

                // Build a graph of all spatially close segments.
                build_graph_of_neighbours();

                // Prepare all necessary data for the sparse QP solver.
                create_input_data_for_qp_solver();

                // Solve the QP problem.
                solve_qp_problem();
                
                // Rotate all segments based on optimized orientations.
                reorient_segments();

                // Update orientations of input segments.
                update_input_segments(input_segments, segment_map);

                // Print debugging information if the corresponding flag is on.
                print_debug_information();

                // Export final regularized segments.
                if (!m_silent) 
                    m_logger.export_segments<SegmentRange, SegmentMap, Kernel>(input_segments, segment_map, "regularized_segments_segment_regularizer");
            }

            void make_silent(const bool new_state) {
                m_silent = new_state;
            }

            // LOD function, can be removed.
            void get_lines_from_segments(const std::vector<Segment> &segments, std::vector<Line> &lines) const {
                
                lines.clear();
                lines.resize(segments.size());

                for (size_t i = 0; i < segments.size(); ++i) {

                    const Point &source = segments[i].source();
                    const Point &target = segments[i].target();

                    lines[i] = Line(source, target);
                }
            }

        private:
            const bool m_debug;
            bool       m_silent;
            
            Logger   m_logger;
            Debugger m_debugger;
            
            Regular_segments m_input_segments, m_final_segments;
            Orientations     m_max_orientations, m_final_orientations;
            
            Parameters            m_parameters;
            Neighbours_graph_data m_neighbours_graph_data;
            QP_problem_data       m_qp_problem_data;
            Max_orientation       m_max_orientation;
            
            template<typename SegmentRange, typename SegmentMap>
            void copy_input_segments(const SegmentRange &input_segments, const SegmentMap &segment_map) {
                
                using Segment_iterator = typename SegmentRange::const_iterator;
                assert(input_segments.size() > 0);

                m_input_segments.clear();
                m_input_segments.resize(input_segments.size());

                size_t i = 0;
                for (Segment_iterator sit = input_segments.begin(); sit != input_segments.end(); ++sit, ++i) {
                    
                    const Segment &segment = get(segment_map, *sit);
                    m_input_segments[i] = Regular_segment(segment);
                }

                assert(i == m_input_segments.size());
            }

            void set_max_orientations() {
                const size_t num_input_segments = m_input_segments.size();
                assert(num_input_segments > 0);

                m_max_orientations.clear();
                m_max_orientations.resize(num_input_segments);

                for (size_t i = 0; i < num_input_segments; ++i) m_max_orientations[i] = m_max_orientation.get();
            }

            void build_graph_of_neighbours() {
                assert(m_input_segments.size() > 0);
                assert(m_max_orientations.size() == m_input_segments.size());

                Neighbours_graph_builder neighbours_graph_builder(m_input_segments, m_max_orientations, m_parameters);
                neighbours_graph_builder.build_graph_data(m_neighbours_graph_data);
            }

            void create_input_data_for_qp_solver() {
                
                assert(m_neighbours_graph_data.filled());
                const size_t num_input_segments = m_input_segments.size();
                
                m_qp_problem_data.set_from(m_neighbours_graph_data, num_input_segments);
            }

            void solve_qp_problem() {
                assert(m_qp_problem_data.filled());

                QP_problem qp_problem(m_max_orientations, m_qp_problem_data, m_parameters);
                qp_problem.solve(m_final_orientations);
            }

            void reorient_segments() {

                assert(m_input_segments.size() > 0);
                assert(m_final_orientations.size() >= m_input_segments.size()); 
                assert(m_qp_problem_data.filled());

                m_final_segments.clear();
                m_final_segments = m_input_segments;

                Tree tree(m_final_segments, m_final_orientations, m_qp_problem_data, m_parameters);
                tree.apply_new_orientations();
            }

            template<typename SegmentRange, typename SegmentMap>
            void update_input_segments(SegmentRange &input_segments, const SegmentMap &segment_map) {

                using Segment_iterator = typename SegmentRange::iterator;
                assert(m_final_segments.size() == input_segments.size());

                size_t i = 0;
                
                for (Segment_iterator sit = input_segments.begin(); sit != input_segments.end(); ++sit, ++i)
                    put(segment_map, *sit, m_final_segments[i].get());

                assert(i == m_final_segments.size());
            }
            
            void print_debug_information() {

                if (!m_debug) return;
                m_debugger.print_values(m_max_orientations, "orientations threshold in degrees");

                RegularMap regular_map;
                m_debugger.print_segments<RegularRange, RegularMap, Kernel>(m_input_segments, regular_map, "segments_before_regularization");
                m_debugger.print_segments<RegularRange, RegularMap, Kernel>(m_final_segments, regular_map, "segments_after_regularization");
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_2_H