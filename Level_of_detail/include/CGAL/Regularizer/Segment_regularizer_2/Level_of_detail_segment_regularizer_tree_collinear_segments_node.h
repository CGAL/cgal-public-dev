#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_COLLINEAR_SEGMENTS_NODE_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_COLLINEAR_SEGMENTS_NODE_H

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
		class Level_of_detail_segment_regularizer_tree_collinear_segments_node {

        public:
            typedef KernelTraits Kernel;
            using FT = typename Kernel::FT;

            using Regular_segment    = Level_of_detail_segment_regularizer_regular_segment<Kernel>;
            using Collinear_segments = std::list<Regular_segment *>;

            Level_of_detail_segment_regularizer_tree_collinear_segments_node() { 
                allocate_memory();
            }

            ~Level_of_detail_segment_regularizer_tree_collinear_segments_node() { 
                deallocate_memory();
            }

            inline const Collinear_segments &get_collinear_segments() const {
                return m_collinear_segments;
            }

            inline Collinear_segments &get_collinear_segments() {
                return m_collinear_segments;
            }

            inline void add(Regular_segment *segment_pointer) {
                m_collinear_segments.push_back(segment_pointer);
            }

        private:
            Collinear_segments m_collinear_segments;

            inline void allocate_memory() {
                m_collinear_segments = Collinear_segments();
            }

            inline void deallocate_memory() {
                m_collinear_segments.clear();
            }
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_COLLINEAR_SEGMENTS_NODE_H