#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_PARALLEL_SEGMENTS_NODE_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_PARALLEL_SEGMENTS_NODE_H

// STL includes.
#include <map>
#include <list>
#include <cassert>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
		class Level_of_detail_segment_regularizer_tree_parallel_segments_node {

        public:
            typedef KernelTraits Kernel;
            using FT = typename Kernel::FT;

            Level_of_detail_segment_regularizer_tree_parallel_segments_node() { }

        private:

		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_PARALLEL_SEGMENTS_NODE_H