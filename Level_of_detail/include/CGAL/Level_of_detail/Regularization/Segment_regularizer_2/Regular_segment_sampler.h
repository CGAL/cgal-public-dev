#ifndef CGAL_LEVEL_OF_DETAIL_REGULAR_SEGMENT_SAMPLER_H
#define CGAL_LEVEL_OF_DETAIL_REGULAR_SEGMENT_SAMPLER_H

// STL includes.
#include <map>
#include <list>

// CGAL includes.
#include <CGAL/number_utils.h>

// LOD includes.
#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Regular_segment.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

        template<class InputKernel>
		class Regular_segment_sampler {

        public:
            using Kernel = InputKernel;
            
            using FT     = typename Kernel::FT;
            using Point  = typename Kernel::Point_2;
            using Vector = typename Kernel::Vector_2;

            using Regular_segment  = LOD::Regular_segment<Kernel>;
            using Regular_segments = std::vector<Regular_segment *>;

            Regular_segment_sampler(const Regular_segments &segments) : 
            m_segments(segments) 
            { }

            template<class Points, class Point_segment_map>
            void sample(Points &points, Point_segment_map &points_to_segments) const {

                points.clear();
                points_to_segments.clear();

                size_t j = 0;
                for (size_t i = 0; i < m_segments.size(); ++i) {
			        const Regular_segment *segment = m_segments[i];

                    const Point &barycentre = segment->get_barycentre();
                    points.push_back(std::make_pair(Point(barycentre.x(), barycentre.y()), j));

			        points_to_segments[j] = i;
			        ++j;
                }
            }

        private:
            const Regular_segments &m_segments;

            void normalize(Vector &vector) const {
                const FT length = static_cast<FT>(CGAL::sqrt(CGAL::to_double(vector * vector)));

				CGAL_precondition(length > FT(0));
				vector /= length;
            }
		};

	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_REGULAR_SEGMENT_SAMPLER_H