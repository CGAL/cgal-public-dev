#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_H

// STL includes.
#include <map>
#include <list>
#include <cassert>

// New CGAL includes.
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_parameters.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_regular_segment.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_tree_parallel_segments_node.h>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits, class QPProblemData>
		class Level_of_detail_segment_regularizer_tree {

        public:
            typedef KernelTraits  Kernel;
            typedef QPProblemData QP_problem_data;

            using FT = typename Kernel::FT;

            using Regular_segment             = Level_of_detail_segment_regularizer_regular_segment<Kernel>;
            using Regular_segments            = std::vector<Regular_segment>;
            using Parallel_segments_tree_node = Level_of_detail_segment_regularizer_tree_parallel_segments_node<Kernel>;

            using Parallel_segments = std::map<FT, Parallel_segments_tree_node>;
            using Other_segments    = std::list<Regular_segment>;

            using Orientations = std::vector<FT>;
            using Parameters   = CGAL::LOD::Level_of_detail_segment_regularizer_parameters<Kernel>;

            Level_of_detail_segment_regularizer_tree(const Orientations &orientations, const QP_problem_data &qp_data, const Parameters &parameters) 
            : m_orientations(orientations), m_qp_data(qp_data), m_parameters(parameters) { 

                build_tree();
            }

            void apply_new_orientations(Regular_segments &segments) {
                assert(segments.size() > 0);

                // to be added!
            }

        private:
            Parallel_segments m_parallel_segments;
            Other_segments    m_other_segments;

            const Orientations    &m_orientations;
            const QP_problem_data &m_qp_data;
            const Parameters      &m_parameters;

            void build_tree() {

                // to be added!
            }
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_H