#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_PARALLEL_SEGMENTS_NODE_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_PARALLEL_SEGMENTS_NODE_H

// STL includes.
#include <map>
#include <list>
#include <memory>
#include <cassert>

// New CGAL includes.
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_regular_segment.h>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
		class Level_of_detail_segment_regularizer_tree_parallel_segments_node {

        public:
            typedef KernelTraits Kernel;
            using FT = typename Kernel::FT;

            using Regular_segment         = Level_of_detail_segment_regularizer_regular_segment<Kernel>;
            using Regular_segment_pointer = std::shared_ptr<Regular_segment>;
            using Parallel_segments       = std::list<Regular_segment_pointer>;

            Level_of_detail_segment_regularizer_tree_parallel_segments_node() { }

            inline void add(const Regular_segment &segment) {
                m_parallel_segments.push_back(std::make_shared<Regular_segment>(segment));
            }

        private:
            Parallel_segments m_parallel_segments;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_PARALLEL_SEGMENTS_NODE_H