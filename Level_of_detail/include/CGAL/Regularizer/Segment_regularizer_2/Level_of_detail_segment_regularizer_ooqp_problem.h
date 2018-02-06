#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_OOQP_PROBLEM_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_OOQP_PROBLEM_H

// STL includes.
#include <map>
#include <list>
#include <cassert>

// New CGAL includes.
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_parameters.h>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits, class QPProblemData>
		class Level_of_detail_segment_regularizer_ooqp_problem {

        public:
            typedef KernelTraits  Kernel;
            typedef QPProblemData QP_problem_data;

            using FT           = typename Kernel::FT;
            using Orientations = std::vector<FT>;

            using Parameters = CGAL::LOD::Level_of_detail_segment_regularizer_parameters<Kernel>;

            Level_of_detail_segment_regularizer_ooqp_problem(const Orientations &x_initial, const QP_problem_data &qp_data, const Parameters &parameters) 
            : m_x_initial(x_initial), m_qp_data(qp_data), m_parameters(parameters) { }

            void solve(Orientations &x) {
                assert(x.size() == 0);

                // to be added!
            }

        private:
            const Orientations    &m_x_initial;
            const QP_problem_data &m_qp_data;
            const Parameters      &m_parameters;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_OOQP_PROBLEM_H