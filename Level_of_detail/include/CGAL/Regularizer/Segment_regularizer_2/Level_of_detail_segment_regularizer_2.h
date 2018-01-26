#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_2_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_2_H

#if defined(WIN32) || defined(_WIN32) 
#define PS "\\"
#else 
#define PS "/" 
#endif 

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <utility>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/property_map.h>

// New CGAL includes.
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_log.h>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
        class Level_of_detail_segment_regularizer_2 {

        public:
            typedef KernelTraits Kernel;

            using FT      = typename Kernel::FT;
            using Point   = typename Kernel::Point_2;
            using Segment = typename Kernel::Segment_2;
			
            using Segments = std::vector<Segment>;
            using Log = CGAL::LOD::Level_of_detail_segment_regularizer_log;

            Level_of_detail_segment_regularizer_2() : m_debug(true) { }

            template<typename SegmentRange, typename SegmentMap>
            void regularize(SegmentRange &input_segments, SegmentMap segment_map) {

                copy_input_segments(input_segments, segment_map);
                if (m_debug) print_segments("segments_before_regularization");

                
            }

        private:
            bool m_debug;
            Segments m_segments;
            
            template<typename SegmentRange, typename SegmentMap>
            void copy_input_segments(SegmentRange &input_segments, SegmentMap &segment_map) {
                using Segment_iterator = typename SegmentRange::const_iterator;

                m_segments.clear();
                m_segments.resize(input_segments.size());

                size_t i = 0;
                for (Segment_iterator sit = input_segments.begin(); sit != input_segments.end(); ++sit, ++i) {
                    
                    const Segment &segment = get(segment_map, *sit);
                    m_segments[i] = segment;
                }

                assert(i == m_segments.size());
            }

            void print_segments(const std::string &name) const {
                assert(!m_segments.empty()); 
                Log log; log.export_segments(m_segments, name);
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_2_H