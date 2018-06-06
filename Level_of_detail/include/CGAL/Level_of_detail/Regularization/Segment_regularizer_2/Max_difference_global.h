#ifndef CGAL_LEVEL_OF_DETAIL_MAX_DIFFERENCE_GLOBAL_H
#define CGAL_LEVEL_OF_DETAIL_MAX_DIFFERENCE_GLOBAL_H

namespace CGAL {

	namespace Level_of_detail {

        template<class InputParameters>
		class Max_difference_global {

        public:
            using Parameters = InputParameters;
            using FT         = typename Parameters::FT;

            Max_difference_global(const Parameters &parameters) : 
            m_parameters(parameters) 
            { }

            FT get() const {
                
                const FT value = m_parameters.max_difference_in_meters();
                CGAL_precondition(value > FT(0));

                return value;
            }

        private:
            const Parameters &m_parameters;
		};

	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_MAX_DIFFERENCE_GLOBAL_H