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
            m_theta_max_deg(5),              // 45 for the final version of the code
            m_d_max_meters(1),               // 5 for the final version of the code 
            m_lambda(FT(4) / FT(5)),
            m_epsilon(FT(1) / FT(4)),
            m_num_intervals_per_segment(10),
            m_optimize_parallelizm(true),
            m_optimize_orthogonality(true)
            { }


            // Setters.
            void set_max_angle_in_degrees(const FT new_value) {
                assert(new_value > FT(0));
                m_theta_max_deg = new_value;
            }

            void set_max_difference_in_meters(const FT new_value) {
                assert(new_value > FT(0));
                m_d_max_meters = new_value;
            }

            void set_lambda(const FT new_value) {
                assert(new_value >= FT(0) && new_value <= FT(1));
                m_lambda = new_value;
            }

            void set_epsilon(const FT new_value) {
                assert(new_value >= FT(0));
                m_epsilon = new_value;
            }

            void set_number_of_intervals_per_segment(const size_t new_value) {
                assert(new_value >= 0);
                m_num_intervals_per_segment = new_value;
            }


            // Getters.
            inline FT get_max_angle_in_degrees() const {
                return m_theta_max_deg;
            }

            inline FT get_max_difference_in_meters() const {
                return m_d_max_meters;
            }

            inline FT get_lambda() const {
                return m_lambda;
            }

            inline FT get_epsilon() const {
                return m_epsilon;
            }

            inline size_t get_number_of_intervals_per_segment() const {
                return m_num_intervals_per_segment;
            }


            // Flags.
            inline bool optimize_parallelizm() const {
                return m_optimize_parallelizm;
            }

            inline bool optimize_orthogonality() const {
                return m_optimize_orthogonality;
            }

        private:

            // Important.
            FT m_theta_max_deg;
            FT m_d_max_meters;


            // Others.
            FT m_lambda;
            FT m_epsilon;

            size_t m_num_intervals_per_segment;


            // Flags.
            bool m_optimize_parallelizm;
            bool m_optimize_orthogonality;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_PARAMETERS_H