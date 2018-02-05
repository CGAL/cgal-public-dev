#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_REGULAR_SEGMENT_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_REGULAR_SEGMENT_H

// STL includes.
#include <map>
#include <list>
#include <cassert>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
		class Level_of_detail_segment_regularizer_regular_segment {

        public:
            typedef KernelTraits Kernel;
            
            using FT      = typename Kernel::FT;
            using Segment = typename Kernel::Segment_2;

            Level_of_detail_segment_regularizer_regular_segment() : m_is_set(false) { }
            Level_of_detail_segment_regularizer_regular_segment(const Segment &segment) : m_segment(segment), m_is_set(true) { }

            const Segment &get() const {
                assert(m_is_set);
                return m_segment;
            }

        private:
            Segment m_segment;
            bool m_is_set;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_REGULAR_SEGMENT_H