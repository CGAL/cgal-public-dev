#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_MAX_ORIENTATION_TEST_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_MAX_ORIENTATION_TEST_H

namespace CGAL {

	namespace LOD {

        template<class InputParameters, class InputSegments>
		class Level_of_detail_segment_regularizer_max_orientation_test {

        public:
            typedef InputParameters Parameters;
            typedef InputSegments   Segments;

            typedef typename Parameters::Kernel Kernel;
            using FT = typename Kernel::FT;

            Level_of_detail_segment_regularizer_max_orientation_test(const Parameters &parameters, const Segments &segments) : 
            m_parameters(parameters), m_segments(segments), m_length_threshold(FT(3)) { }

            FT get(const size_t segment_index) const {
                
                const FT value = m_parameters.get_max_angle_in_degrees();
                assert(value > FT(0));

                const FT decreased_value = value / FT(4);
                if (is_too_long_segment(segment_index)) return decreased_value;
                
                const FT increased_value = value * FT(4);
                return increased_value;
            }

        private:
            const Parameters &m_parameters;
            const Segments   &m_segments;
            const FT m_length_threshold;

            bool is_too_long_segment(const size_t segment_index) const {
                if (m_segments[segment_index].get().squared_length() < m_length_threshold * m_length_threshold) return false;
                return true;
            }
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_MAX_ORIENTATION_TEST_H