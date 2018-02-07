#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_REGULAR_SEGMENT_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_REGULAR_SEGMENT_H

// STL includes.
#include <map>
#include <list>
#include <cmath>
#include <cassert>

// CGAL includes.
#include <CGAL/number_utils.h>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
		class Level_of_detail_segment_regularizer_regular_segment {

        public:
            typedef KernelTraits Kernel;
            
            using FT      = typename Kernel::FT;
            using Segment = typename Kernel::Segment_2;
            using Vector  = typename Kernel::Vector_2;

            Level_of_detail_segment_regularizer_regular_segment() : 
            m_orientation(-FT(1)), m_is_set(false) { }

            Level_of_detail_segment_regularizer_regular_segment(const Segment &segment) : 
            m_segment(segment), m_orientation(-FT(1)), m_is_set(true) { 

                compute_initial_orientation();
            }

            Segment &get() {
                assert(m_is_set);
                return m_segment;
            }

            const Segment &get() const {
                assert(m_is_set);
                return m_segment;
            }

            FT get_orientation() const {
                assert(m_is_set);
                return m_orientation;
            }

        private:
            Segment m_segment;
            FT      m_orientation;
            bool    m_is_set;

            void compute_initial_orientation() {

                Vector direction = m_segment.to_vector();
                if (direction.y() < FT(0) || (direction.y() == FT(0) && direction.x() < FT(0))) direction = -direction;

                const FT atan = static_cast<FT>(std::atan2(CGAL::to_double(direction.y()), CGAL::to_double(direction.x())));
                m_orientation = atan * FT(180) / static_cast<FT>(CGAL_PI);
            }
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_REGULAR_SEGMENT_H