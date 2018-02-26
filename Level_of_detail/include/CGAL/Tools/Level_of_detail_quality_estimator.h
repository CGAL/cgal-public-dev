#ifndef CGAL_LEVEL_OF_DETAIL_QUALITY_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_QUALITY_ESTIMATOR_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\"
#else 
#define PSR "/" 
#endif

// STL includes.
#include <map>
#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/IO/Color.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class LodQuality>
		class Level_of_detail_quality_estimator {

		public:
			typedef LodQuality Lod_quality;
			
			using Kernel = typename Lod_quality::Kernel;
			using FT 	 = typename Kernel::FT;
			
			using Params = char**;
			using Data = std::vector<FT>;

			using Log = CGAL::LOD::Mylog;

			Level_of_detail_quality_estimator(const int num_params, const Params params) 
			: m_lod_quality(num_params, params), m_debug(false) { }

			void run_quality_test() {
				m_lod_quality.compute_data();

				// X data.
				const Data &x_data = m_lod_quality.retreive_x_data();
				
				// Y data: complexity.
				save_complexity_y_data(x_data);

				// Y data: distortion.
				save_distortion_y_data(x_data);

				// Y data: coverage.
				save_coverage_y_data(x_data);
			}

		private:
			Lod_quality m_lod_quality;
			const bool m_debug;

			void save_complexity_y_data(const Data &x_data) {
				
				const Data &cmp_roofs_y = m_lod_quality.retreive_cmp_y_data(Quality_data_type::CMP_ROOFS);
				const Data &cmp_walls_y = m_lod_quality.retreive_cmp_y_data(Quality_data_type::CMP_WALLS);
				const Data &cmp_total_y = m_lod_quality.retreive_cmp_y_data(Quality_data_type::CMP);

				std::vector<Data> y_data(3);
				y_data[0] = cmp_roofs_y;
				y_data[1] = cmp_walls_y;
				y_data[2] = cmp_total_y;

				Log saver;
				saver.save_quality_data_final("cmp_data", x_data, y_data);
			}

			void save_distortion_y_data(const Data &x_data) {

				save_distortion_roofs_y_data(x_data);
				save_distortion_walls_y_data(x_data);
			}

			void save_distortion_roofs_y_data(const Data &x_data) {
				
				const std::vector<Data> &y_data = m_lod_quality.retreive_dst_y_data(Quality_data_type::DST_ROOFS);

				Log saver;
				saver.save_quality_data_intermediate("dst_roofs_data", x_data, y_data);
			}

			void save_distortion_walls_y_data(const Data &x_data) {
				
				const std::vector<Data> &y_data = m_lod_quality.retreive_dst_y_data(Quality_data_type::DST_WALLS);

				Log saver;
				saver.save_quality_data_intermediate("dst_walls_data", x_data, y_data);
			}

			void save_coverage_y_data(const Data &x_data) {
				
				const Data &cov_roofs_y = m_lod_quality.retreive_cov_y_data(Quality_data_type::COV_ROOFS);
				const Data &cov_walls_y = m_lod_quality.retreive_cov_y_data(Quality_data_type::COV_WALLS);
				const Data &cov_total_y = m_lod_quality.retreive_cov_y_data(Quality_data_type::COV);

				std::vector<Data> y_data(3);
				y_data[0] = cov_roofs_y;
				y_data[1] = cov_walls_y;
				y_data[2] = cov_total_y;

				Log saver;
				saver.save_quality_data_final("cov_data", x_data, y_data);
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_QUALITY_ESTIMATOR_H