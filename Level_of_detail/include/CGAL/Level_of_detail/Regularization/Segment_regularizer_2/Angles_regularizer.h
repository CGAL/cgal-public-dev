#ifndef CGAL_LEVEL_OF_DETAIL_ANGLES_REGULARIZER_H
#define CGAL_LEVEL_OF_DETAIL_ANGLES_REGULARIZER_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <memory>
#include <utility>

// CGAL includes.
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Segment_regularizer_tree.h>
#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Segment_regularizer_parameters.h>

#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Max_orientation_local.h>
#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Max_orientation_global.h>

#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Regular_segment.h>
#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Regular_segment_property_map.h>

#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Neighbours_graph_data.h>
#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Delaunay_neighbours_graph_builder.h>

#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/OOQP_problem.h>
#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/OOQP_problem_data.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

        template<class InputKernel>
        class Angles_regularizer {

        public:
            using Kernel = InputKernel;

            using FT      = typename Kernel::FT;
            using Point   = typename Kernel::Point_2;
            using Segment = typename Kernel::Segment_2;
            using Line    = typename Kernel::Line_2;

            using Regular_segment  = LOD::Regular_segment<Kernel>;
            using Regular_segments = std::vector<Regular_segment *>;

            using RegularMap   = LOD::Regular_segment_property_map<Regular_segment, Segment>;
            using RegularRange = Regular_segments;
            
            using Parameters   = LOD::Segment_regularizer_parameters<FT>;
            using Orientations = std::vector<FT>;

            using Max_orientation_local  = LOD::Max_orientation_local<Kernel, Parameters, Regular_segments>;
            using Max_orientation_global = LOD::Max_orientation_global<Parameters>;
            
            using Neighbours_graph_data    = LOD::Neighbours_graph_data<Kernel>;
            using Neighbours_graph_builder = LOD::Delaunay_neighbours_graph_builder<Kernel, Neighbours_graph_data>;

            using QP_problem_data = LOD::OOQP_problem_data<Kernel, Neighbours_graph_data>;
            using QP_problem      = LOD::OOQP_problem<Kernel, QP_problem_data>;

            using Tree = LOD::Segment_regularizer_tree<Kernel, QP_problem_data>;

            Angles_regularizer(Regular_segments &segments, const Parameters &parameters) :
            m_input_segments(segments), 
            m_parameters(parameters) 
            { }

            ~Angles_regularizer() {
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
            }

            Tree *get_tree_pointer() {
                return m_tree_pointer;
            }

        private:
            Regular_segments &m_input_segments;
            const Parameters &m_parameters;
            
            Orientations m_max_orientations;
            Orientations m_final_orientations;
            
            Neighbours_graph_data m_neighbours_graph_data;
            QP_problem_data       m_qp_problem_data;

            Tree *m_tree_pointer;

            void set_max_orientations() {

                const size_t num_input_segments = m_input_segments.size();
                CGAL_precondition(num_input_segments > 0);

                const Max_orientation_local   max_orientation_local(m_parameters, m_input_segments);
                const Max_orientation_global max_orientation_global(m_parameters);
                
                m_max_orientations.clear();
                m_max_orientations.resize(num_input_segments);

                if (m_parameters.use_local_orientation()) for (size_t i = 0; i < num_input_segments; ++i) m_max_orientations[i] = max_orientation_local.get(i);
                else for (size_t i = 0; i < num_input_segments; ++i) m_max_orientations[i] = max_orientation_global.get();
            }

            void build_graph_of_neighbours() {
                
                CGAL_precondition(m_input_segments.size() > 0);
                CGAL_precondition(m_max_orientations.size() == m_input_segments.size());

                Neighbours_graph_builder neighbours_graph_builder(m_input_segments, m_max_orientations, m_parameters);
                neighbours_graph_builder.build_graph_data(m_neighbours_graph_data);
            }

            void create_input_data_for_qp_solver() {
                
                CGAL_precondition(m_neighbours_graph_data.filled());
                const size_t num_input_segments = m_input_segments.size();
                
                m_qp_problem_data.set_from(m_neighbours_graph_data, num_input_segments);
            }

            void solve_qp_problem() {
                CGAL_precondition(m_qp_problem_data.filled());

                const QP_problem qp_problem(m_max_orientations, m_qp_problem_data, m_parameters, m_input_segments);
                qp_problem.solve(m_final_orientations);
            }

            void reorient_segments() {

                CGAL_precondition(m_input_segments.size() > 0);
                CGAL_precondition(m_final_orientations.size() >= m_input_segments.size()); 
                CGAL_precondition(m_qp_problem_data.filled());

                m_tree_pointer = new Tree(m_input_segments, m_final_orientations, m_qp_problem_data, m_parameters);
                m_tree_pointer->apply_new_orientations();
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_ANGLES_REGULARIZER_H