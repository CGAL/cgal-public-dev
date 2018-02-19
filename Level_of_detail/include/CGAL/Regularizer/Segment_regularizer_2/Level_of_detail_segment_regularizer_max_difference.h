#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_MAX_DIFFERENCE_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_MAX_DIFFERENCE_H

namespace CGAL {

	namespace LOD {

        template<class InputParameters>
		class Level_of_detail_segment_regularizer_max_difference {

        public:
            typedef InputParameters Parameters;
            typedef typename Parameters::Kernel Kernel;

            using FT = typename Kernel::FT;

            Level_of_detail_segment_regularizer_max_difference(const Parameters &parameters) : m_parameters(parameters) { }

            FT get() const {
                
                const FT value = m_parameters.get_max_difference_in_meters();
                assert(value > FT(0));

                return value;
            }

        private:
            const Parameters &m_parameters;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_MAX_DIFFERENCE_H