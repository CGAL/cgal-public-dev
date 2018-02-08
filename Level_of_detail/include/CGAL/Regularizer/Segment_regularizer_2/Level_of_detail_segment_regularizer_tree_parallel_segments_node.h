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

            using Regular_segment            = Level_of_detail_segment_regularizer_regular_segment<Kernel>;
            using Parallel_segments          = std::list<Regular_segment *>;
            using Parallel_segments_iterator = typename Parallel_segments::const_iterator;

            Level_of_detail_segment_regularizer_tree_parallel_segments_node() { 
                allocate_memory();
            }

            ~Level_of_detail_segment_regularizer_tree_parallel_segments_node() { 
                deallocate_memory();
            }

            inline const Parallel_segments &get_parallel_segments() const {
                return m_parallel_segments;
            }

            inline void add(Regular_segment *segment_pointer) {
                m_parallel_segments.push_back(segment_pointer);
            }

        private:
            Parallel_segments m_parallel_segments;

            inline void allocate_memory() {
                m_parallel_segments = Parallel_segments();
            }

            void deallocate_memory() {
                m_parallel_segments.clear();
            }
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_PARALLEL_SEGMENTS_NODE_H