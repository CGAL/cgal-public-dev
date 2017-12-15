#ifndef CGAL_LEVEL_OF_DETAIL_QUALITY_H
#define CGAL_LEVEL_OF_DETAIL_QUALITY_H

// STL includes.
#include <cassert>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits>
		class Level_of_detail_quality {

		public:
			typedef KernelTraits Kernel;
			typedef typename Kernel::FT FT;

			Level_of_detail_quality(const FT initial_complexity, const FT initial_distortion) : m_initial_complexity(initial_complexity), m_initial_distortion(initial_distortion),
			m_complexity(-FT(1)), m_distortion(-FT(1)), m_total_quality(-FT(1)) { 

				assert(initial_complexity >= FT(0));
				assert(initial_distortion >= FT(0));

				compute_final_metrics();
			}

			FT get_complexity() const {
				assert(m_complexity >= FT(0));
				return m_complexity;
			}

			FT get_distortion() const {
				assert(m_distortion >= FT(0));
				return m_distortion;
			}

			FT get_total_quality() const {
				assert(m_total_quality >= FT(0));
				return m_total_quality;
			}

		private:
			const FT m_initial_complexity;
			const FT m_initial_distortion;

			FT m_complexity;
			FT m_distortion;
			FT m_total_quality;

			void compute_final_metrics() {

				const FT sum = m_initial_complexity + m_initial_distortion;
				m_complexity = m_initial_complexity / sum;
				m_distortion = m_initial_distortion / sum;

				m_total_quality = m_complexity * m_distortion;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_QUALITY_H