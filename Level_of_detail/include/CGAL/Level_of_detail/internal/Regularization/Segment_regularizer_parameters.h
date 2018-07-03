#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_PARAMETERS_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_PARAMETERS_H

namespace CGAL {

	namespace Level_of_detail {

        template<typename NumberType>
		class Segment_regularizer_parameters {

        public:
			using FT = NumberType;

            Segment_regularizer_parameters() : 
            m_theta_max_deg(FT(45)),   
            m_d_max_meters(FT(1) / FT(2)),

			// these are weak parameters
            m_lambda(FT(4) / FT(5)),        		   
            m_epsilon(FT(1) / FT(4)),       		   
			m_small_fixed_orientation(FT(1) / FT(10)), 
			m_small_fixed_difference(FT(1) / FT(10)),  
			m_tolerance(FT(1) / FT(1000000)),		   
			m_ooqp_problem_weight(FT(100000)),		   
			
			// these are options
            m_optimize_parallelizm(true),   		   
            m_optimize_orthogonality(true), 		   
            m_optimize_angles(true),        		   
            m_optimize_ordinates(true),     		   
			m_use_local_orientation(true),			   
			m_use_local_difference(true),			   
            m_verbose(false)                		   
            { }


            // Important.
            inline FT& max_angle_in_degrees() {
				return m_theta_max_deg;
			}

			inline const FT& max_angle_in_degrees() const {
				return m_theta_max_deg;
			}

            inline FT& max_difference_in_meters() {
				return m_d_max_meters;
			}

			inline const FT& max_difference_in_meters() const {
				return m_d_max_meters;
			}


            // Less important.
            inline FT& lambda() {
				return m_lambda;
			}

			inline const FT& lambda() const {
				return m_lambda;
			}

            inline FT& epsilon() {
				return m_epsilon;
			}

			inline const FT& epsilon() const {
				return m_epsilon;
			}

			inline FT& small_fixed_orientation() {
				return m_small_fixed_orientation;
			}

			inline const FT& small_fixed_orientation() const {
				return m_small_fixed_orientation;
			}

			inline FT& small_fixed_difference() {
				return m_small_fixed_difference;
			}

			inline const FT& small_fixed_difference() const {
				return m_small_fixed_difference;
			}

			inline FT& tolerance() {
				return m_tolerance;
			}

			inline const FT& tolerance() const {
				return m_tolerance;
			}

			inline FT& ooqp_problem_weight() {
				return m_ooqp_problem_weight;
			}

			inline const FT& ooqp_problem_weight() const {
				return m_ooqp_problem_weight;
			}


            // Flags.
			inline bool& optimize_parallelizm() {
				return m_optimize_parallelizm;
			}

			inline const bool& optimize_parallelizm() const {
				return m_optimize_parallelizm;
			}

			inline bool& optimize_orthogonality() {
				return m_optimize_orthogonality;
			}

			inline const bool& optimize_orthogonality() const {
				return m_optimize_orthogonality;
			}

			inline bool& optimize_angles() {
				return m_optimize_angles;
			}

			inline const bool& optimize_angles() const {
				return m_optimize_angles;
			}

			inline bool& optimize_ordinates() {
				return m_optimize_ordinates;
			}

			inline const bool& optimize_ordinates() const {
				return m_optimize_ordinates;
			}

			inline bool& use_local_orientation() {
				return m_use_local_orientation;
			}

			inline const bool& use_local_orientation() const {
				return m_use_local_orientation;
			}

			inline bool& use_local_difference() {
				return m_use_local_difference;
			}

			inline const bool& use_local_difference() const {
				return m_use_local_difference;
			}

			inline bool& verbose() {
				return m_verbose;
			}

			inline const bool& verbose() const {
				return m_verbose;
			}

        private:

            // Important.
            FT m_theta_max_deg;
            FT m_d_max_meters;

            // Less important.
            FT m_lambda;
            FT m_epsilon;
			FT m_small_fixed_orientation;
			FT m_small_fixed_difference;
			FT m_tolerance;
			FT m_ooqp_problem_weight;

            // Flags.
            bool m_optimize_parallelizm;
            bool m_optimize_orthogonality;
            bool m_optimize_angles;
            bool m_optimize_ordinates;
			bool m_use_local_orientation;
			bool m_use_local_difference;
            bool m_verbose;
		};

	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_PARAMETERS_H