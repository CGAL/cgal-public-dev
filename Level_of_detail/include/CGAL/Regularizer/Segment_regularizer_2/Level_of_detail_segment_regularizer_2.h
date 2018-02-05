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
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_regular_segment.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_max_initial_orientation.h>

#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_neighbours_graph_data.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_delaunay_neighbours_graph_builder.h>

#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_qp_solver_data.h>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
        class Level_of_detail_segment_regularizer_2 {

        public:
            typedef KernelTraits Kernel;

            using FT      = typename Kernel::FT;
            using Point   = typename Kernel::Point_2;
            using Segment = typename Kernel::Segment_2;

            using Regular_segment  = CGAL::LOD::Level_of_detail_segment_regularizer_regular_segment<Kernel>;
            using Regular_segments = std::vector<Regular_segment>;
            
            using Debugger   = CGAL::LOD::Level_of_detail_segment_regularizer_debugger;
            using Parameters = CGAL::LOD::Level_of_detail_segment_regularizer_parameters<Kernel>;
            using Tree       = CGAL::LOD::Level_of_detail_segment_regularizer_tree<Kernel>;

            using Initial_orientation = CGAL::LOD::Level_of_detail_segment_regularizer_max_initial_orientation<Parameters>;
            using Orientations        = std::vector<FT>;

            using Neighbours_graph_data    = CGAL::LOD::Level_of_detail_segment_regularizer_neighbours_graph_data<Kernel>;
            using Neighbours_graph_builder = CGAL::LOD::Level_of_detail_segment_regularizer_delaunay_neighbours_graph_builder<Kernel>;

            using QP_solver_data = CGAL::LOD::Level_of_detail_segment_regularizer_qp_solver_data<Kernel>;

            Level_of_detail_segment_regularizer_2() : m_debug(true), m_parameters(), m_initial_orientation(m_parameters) { }

            template<typename SegmentRange, typename SegmentMap>
            void regularize(SegmentRange &input_segments, SegmentMap segment_map) {
                if (input_segments.size() == 0) return;

                // Copy input segments in the internal storage.
                copy_input_segments(input_segments, segment_map);
                
                // Set initial orientations for all segments. These orientations will later be modified.
                set_initial_orientations();

                // Build a graph of all spatially close segments.
                build_graph_of_neighbours();

                // Prepare all necessary data for the sparse QP solver.
                create_input_data_for_qp_solver();
                
                // Print debugging information if the corresponding flag is on.
                print_debug_information();
            }

        private:
            bool m_debug;
            Debugger m_debugger;
            
            Regular_segments m_segments;
            Orientations     m_x;
            
            Parameters            m_parameters;
            Neighbours_graph_data m_neighbours_graph_data;
            QP_solver_data        m_qp_solver_data;

            Tree                m_tree;
            Initial_orientation m_initial_orientation;
            
            template<typename SegmentRange, typename SegmentMap>
            void copy_input_segments(const SegmentRange &input_segments, const SegmentMap &segment_map) {
                
                using Segment_iterator = typename SegmentRange::const_iterator;
                assert(input_segments.size() > 0);

                m_segments.clear();
                m_segments.resize(input_segments.size());

                size_t i = 0;
                for (Segment_iterator sit = input_segments.begin(); sit != input_segments.end(); ++sit, ++i) {
                    
                    const Segment &segment = get(segment_map, *sit);
                    m_segments[i] = Regular_segment(segment);
                }

                assert(i == m_segments.size());
            }

            void set_initial_orientations() {
                const size_t num_segments = m_segments.size();
                assert(num_segments > 0);

                m_x.clear();
                m_x.resize(num_segments);

                for (size_t i = 0; i < num_segments; ++i) m_x[i] = m_initial_orientation.get();
            }

            void build_graph_of_neighbours() {
                assert(m_segments.size() > 0 && m_x.size() == m_segments.size());

                Neighbours_graph_builder neighbours_graph_builder(m_segments, m_x, m_parameters);
                neighbours_graph_builder.build_graph(m_neighbours_graph_data);
            }

            void create_input_data_for_qp_solver() {
                assert(m_neighbours_graph_data.filled());

            }

            void print_debug_information() {
                if (!m_debug) return;

                m_debugger.print_segments(m_segments, "segments_before_regularization");
                m_debugger.print_values(m_x, "initial orientations");
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_2_H