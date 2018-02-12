#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_SEGMENT_SAMPLER_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_SEGMENT_SAMPLER_H

// STL includes.
#include <map>
#include <list>
#include <cassert>

// CGAL includes.
#include <CGAL/number_utils.h>

// New CGAL includes.
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_debugger.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_regular_segment.h>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
		class Level_of_detail_segment_regularizer_segment_sampler {

        public:
            typedef KernelTraits Kernel;
            
            using FT     = typename Kernel::FT;
            using Point  = typename Kernel::Point_2;
            using Vector = typename Kernel::Vector_2;

            using Regular_segment  = CGAL::LOD::Level_of_detail_segment_regularizer_regular_segment<Kernel>;
            using Regular_segments = std::vector<Regular_segment>;

            using Debugger = CGAL::LOD::Level_of_detail_segment_regularizer_debugger;

            Level_of_detail_segment_regularizer_segment_sampler(const Regular_segments &segments) 
            : m_debug(false), m_segments(segments) { }

            template<class Points, class Point_segment_map>
            void sample(Points &points, Point_segment_map &points_to_segments, const size_t /* num_intervals_per_segment */) {

                points.clear();
                points_to_segments.clear();

                size_t j = 0;
                for (size_t i = 0; i < m_segments.size(); ++i) {
			        const Regular_segment &segment = m_segments[i];

                    /*
                    const Point &source = segment.get().source();
			        const Point &target = segment.get().target();

			        Vector direction = segment.get_direction();
                    normalize(direction);

                    const FT segment_length = segment.get_length();
                    */

                    // FIX THIS!
                    // What is wrong with this part?
                    // If I sample segments, the connectivity graph looks absolutely wrong!

                    // <---------------- new part
                    /*

                    const FT num_steps = static_cast<FT>(num_intervals_per_segment) / FT(2);
                    const FT h         = segment_length / (FT(2) * num_steps);

			        for (size_t k = 0; k < num_steps; ++k) {

                        points.push_back(std::make_pair(Point(source.x() + k * direction.x() * h, source.y() + k * direction.y() * h), j));
				        points_to_segments[j] = i;
				        ++j;

                        points.push_back(std::make_pair(Point(target.x() - k * direction.x() * h, target.y() - k * direction.y() * h), j));
				        points_to_segments[j] = i;
				        ++j;
                    }

                    */
                    // ---------------- old part
                    /*

                    const FT ds = static_cast<FT>(num_intervals_per_segment);
                    const size_t num_subdivs = static_cast<size_t>(std::floor(CGAL::to_double(segment_length / (FT(2) * ds))));

                    for (size_t k = 0; k < num_subdivs; ++k) {

                        points.push_back(std::make_pair(Point(source.x() + k * direction.x(), source.y() + k * direction.y()), j));
                        points_to_segments[j] = i;
                        ++j;

                        points.push_back(std::make_pair(Point(target.x() - k * direction.x(), target.y() - k * direction.y()), j));
                        points_to_segments[j] = i;
                        ++j;
                    }
                    
                    */
                    // ---------------->

                    const Point &barycentre = segment.get_barycentre();

                    points.push_back(std::make_pair(Point(barycentre.x(), barycentre.y()), j));
			        points_to_segments[j] = i;
			        ++j;
                }

                print_debug_information(points);
            }

        private:
            const bool m_debug;
            Debugger   m_debugger;

            const Regular_segments &m_segments;

            void normalize(Vector &vector) const {
                const FT length = static_cast<FT>(CGAL::sqrt(CGAL::to_double(vector * vector)));

				assert(length > FT(0));
				vector /= length;
            }

            template<class Points>
            void print_debug_information(const Points &points) {
                if (!m_debug) return;

                m_debugger.print_sampled_segments(points, "regularizer_sampled_segments");
            }
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_SEGMENT_SAMPLER_H