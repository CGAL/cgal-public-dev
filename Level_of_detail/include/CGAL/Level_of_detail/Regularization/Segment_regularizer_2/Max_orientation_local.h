#ifndef CGAL_LEVEL_OF_DETAIL_MAX_ORIENTATION_LOCAL_H
#define CGAL_LEVEL_OF_DETAIL_MAX_ORIENTATION_LOCAL_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/number_utils.h>

// LOD includes.
#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Regular_segments_info_estimator.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

        template<class InputParameters, class InputSegments>
		class Max_orientation_local {

        public:
            using Parameters = InputParameters;
            using Segments   = InputSegments;

            using Kernel = typename Parameters::Kernel;
            using FT     = typename Kernel::FT;

            using Segments_info_estimator = LOD::Regular_segments_info_estimator<Kernel, Segments>;

            Max_orientation_local(const Parameters &parameters, const Segments &segments) : 
            m_parameters(parameters), 
            m_segments(segments), 
            m_small_fixed_orientation(parameters.small_fixed_orientation()) 
            { }

            FT get(const size_t segment_index) const {
                
                const FT value = m_parameters.max_angle_in_degrees();
                CGAL_precondition(value > FT(0));

                const Segments_info_estimator segments_info_estimator(m_segments);
                if (segments_info_estimator.is_too_long_segment(segment_index)) return m_small_fixed_orientation;

                return value;
            }

        private:
            const Parameters &m_parameters;
            const Segments   &m_segments;
            
            const FT m_small_fixed_orientation;
		};

	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_MAX_ORIENTATION_LOCAL_H