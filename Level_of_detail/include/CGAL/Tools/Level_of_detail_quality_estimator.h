#ifndef CGAL_LEVEL_OF_DETAIL_QUALITY_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_QUALITY_ESTIMATOR_H

// STL includes.
#include <map>
#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>

namespace CGAL {

	namespace LOD {

		template<class LodQuality>
		class Level_of_detail_quality_estimator{

		public:
			typedef LodQuality Lod_quality;
			using Params = char**;

			Level_of_detail_quality_estimator(const int num_params, const Params params) 
			: m_lod_quality(num_params, params) { }

			void run_quality_test() {
				m_lod_quality.compute_xy_data();

				
			}

		private:
			Lod_quality m_lod_quality;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_QUALITY_ESTIMATOR_H