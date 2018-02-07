#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_MAX_ORIENTATION_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_MAX_ORIENTATION_H

namespace CGAL {

	namespace LOD {

        template<class InputParameters>
		class Level_of_detail_segment_regularizer_max_orientation {

        public:
            typedef InputParameters Parameters;
            typedef typename Parameters::Kernel Kernel;

            using FT = typename Kernel::FT;

            Level_of_detail_segment_regularizer_max_orientation(const Parameters &parameters) : m_parameters(parameters) { }

            inline FT get() const {
                return m_parameters.get_max_angle_in_degrees();
            }

        private:
            const Parameters &m_parameters;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_MAX_ORIENTATION_H