#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_LENGTH_THRESHOLD_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_LENGTH_THRESHOLD_ESTIMATOR_H

// STL includes.
#include <vector>
#include <cassert>

// CGAL includes.
#include <CGAL/number_utils.h>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits, class InputSegments>
		class Level_of_detail_segment_regularizer_length_threshold_estimator {

        public:
            typedef KernelTraits  Kernel;
            typedef InputSegments Segments;
            
            using FT = typename Kernel::FT;

            Level_of_detail_segment_regularizer_length_threshold_estimator(const Segments &segments) : 
            m_segments(segments), m_length_threshold(-FT(1)) { 

                compute_length_threshold();
            }

            inline FT get_length_threshold() const {
                assert(m_length_threshold > FT(0));
                return m_length_threshold;
            }

            bool is_too_long_segment(const size_t segment_index) const {
                assert(m_length_threshold > FT(0));

                if (m_segments[segment_index]->get().squared_length() < m_length_threshold * m_length_threshold) return false;
                return true;
            }

        private:
            const Segments &m_segments;
            FT m_length_threshold;

            void compute_length_threshold() {
                
                std::vector<FT> segment_lengths;
                compute_segment_lengths(segment_lengths);

                const FT mean = compute_mean(segment_lengths);
                const FT stde = compute_standard_deviation(segment_lengths, mean);

                const FT estimated_length_threshold = estimate_length_threshold(mean, stde);
                m_length_threshold = estimated_length_threshold;
            }

            void compute_segment_lengths(std::vector<FT> &segment_lengths) const {
                assert(m_segments.size() > 0);

                segment_lengths.clear();
                segment_lengths.resize(m_segments.size());

                for (size_t i = 0; i < m_segments.size(); ++i)
                    segment_lengths[i] = compute_segment_length(i);
            }

            FT compute_segment_length(const size_t segment_index) const {
                return static_cast<FT>(CGAL::sqrt(CGAL::to_double(m_segments[segment_index]->get().squared_length())));
            } 

            FT compute_mean(const std::vector<FT> &values) const {
                assert(values.size() > 0);

                FT mean = FT(0);
                for (size_t i = 0; i < values.size(); ++i)
                    mean += values[i];
                mean /= static_cast<FT>(values.size());

                return mean;
            }

            FT compute_standard_deviation(const std::vector<FT> &values, const FT mean) const {
                assert(values.size() > 0);

                FT sum = FT(0);
                for (size_t i = 0; i < values.size(); ++i)
                    sum += (values[i] - mean) * (values[i] - mean);
                sum /= static_cast<FT>(values.size());

                return static_cast<FT>(CGAL::sqrt(CGAL::to_double(sum)));
            }

            FT estimate_length_threshold(const FT mean, const FT stde) const {
                return mean + stde;
            }
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_LENGTH_THRESHOLD_ESTIMATOR_H