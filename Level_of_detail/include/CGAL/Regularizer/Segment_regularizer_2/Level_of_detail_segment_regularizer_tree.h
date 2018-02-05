#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_H

// STL includes.
#include <map>
#include <list>
#include <cassert>

// New CGAL includes.
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_regular_segment.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_tree_parallel_segments_node.h>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
		class Level_of_detail_segment_regularizer_tree {

        public:
            typedef KernelTraits Kernel;
            using FT = typename Kernel::FT;

            using Regular_segment             = Level_of_detail_segment_regularizer_regular_segment<Kernel>;
            using Parallel_segments_tree_node = Level_of_detail_segment_regularizer_tree_parallel_segments_node<Kernel>;

            using Parallel_segments = std::map<FT, Parallel_segments_tree_node>;
            using Other_segments    = std::list<Regular_segment>;

            Level_of_detail_segment_regularizer_tree() { }

        private:
            Parallel_segments m_parallel_segments;
            Other_segments    m_other_segments;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_H