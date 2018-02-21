#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_MAX_DIFFERENCE_TEST_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_MAX_DIFFERENCE_TEST_H

// STL includes.
#include <vector>
#include <cassert>

// CGAL includes.
#include <CGAL/number_utils.h>

// New CGAL includes.
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_length_threshold_estimator.h>

namespace CGAL {

	namespace LOD {

        template<class InputParameters, class InputSegments>
		class Level_of_detail_segment_regularizer_max_difference_test {

        public:
            typedef InputParameters Parameters;
            typedef InputSegments   Segments;

            typedef typename Parameters::Kernel Kernel;
            using FT = typename Kernel::FT;

            typedef CGAL::LOD::Level_of_detail_segment_regularizer_length_threshold_estimator<Kernel, Segments> Length_threshold_estimator;

            Level_of_detail_segment_regularizer_max_difference_test(const Parameters &parameters, const Segments &segments) : 
            m_parameters(parameters), m_segments(segments), m_small_fixed_difference(FT(1) / FT(10)) { }

            FT get(const size_t segment_index) const {

                const FT value = m_parameters.get_max_difference_in_meters();
                assert(value > FT(0));

                Length_threshold_estimator length_threshold_estimator(m_segments);
                if (length_threshold_estimator.is_too_long_segment(segment_index)) return m_small_fixed_difference;
                
                return value;
            }

        private:
            const Parameters &m_parameters;
            const Segments   &m_segments;
        
            const FT m_small_fixed_difference;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_MAX_DIFFERENCE_TEST_H