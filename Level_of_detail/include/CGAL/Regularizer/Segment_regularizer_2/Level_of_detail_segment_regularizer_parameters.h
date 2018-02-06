#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_PARAMETERS_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_PARAMETERS_H

// STL includes.
#include <cassert>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
		class Level_of_detail_segment_regularizer_parameters {

        public:
            typedef KernelTraits Kernel;
            using FT = typename Kernel::FT;

            Level_of_detail_segment_regularizer_parameters() : 
            m_theta_max_deg(5), 
            m_lambda(FT(4) / FT(5)),
            m_epsilon(FT(1) / FT(4))
            { }

            // Setters.
            void set_max_angle_in_degrees(const FT new_value) {
                assert(new_value >= FT(0));
                m_theta_max_deg = new_value;
            }

            void set_lambda(const FT new_value) {
                assert(new_value >= FT(0) && new_value <= FT(1));
                m_lambda = new_value;
            }

            void set_epsilon(const FT new_value) {
                assert(new_value >= FT(0));
                m_epsilon = new_value;
            }

            // Getters.
            inline FT get_max_angle_in_degrees() const {
                return m_theta_max_deg;
            }

            inline FT get_lambda() const {
                return m_lambda;
            }

            inline FT get_epsilon() const {
                return m_epsilon;
            }

        private:
            FT m_theta_max_deg;
            FT m_lambda;
            FT m_epsilon;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_PARAMETERS_H