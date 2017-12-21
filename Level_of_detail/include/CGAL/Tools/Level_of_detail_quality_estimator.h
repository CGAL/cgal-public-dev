#ifndef CGAL_LEVEL_OF_DETAIL_QUALITY_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_QUALITY_ESTIMATOR_H

#if defined(WIN32) || defined(_WIN32) 
#define PS "\\"
#else 
#define PS "/" 
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
		class Level_of_detail_quality_estimator{

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
				m_lod_quality.compute_xy_data(false);

				const Data &x = m_lod_quality.retreive_x_data();
				const Data &min_y = m_lod_quality.retreive_y_data(Distortion_fitting_type::MIN);
				const Data &avg_y = m_lod_quality.retreive_y_data(Distortion_fitting_type::AVG);
				const Data &max_y = m_lod_quality.retreive_y_data(Distortion_fitting_type::MAX);

				std::vector<Data> data(3);
				data[0] = min_y;
				data[1] = avg_y;
				data[2] = max_y;

				std::vector<Color> colors(3);
				colors[0] = Color(255, 153, 153);
				colors[1] = Color(255, 80, 80);
				colors[2] = Color(128, 0, 0);

				if (m_debug) {
					Log plotter;
					plotter.plot_2d("quality", x, data, colors);
				}

				Log saver;
				saver.save_quality_data("data", x, data);
			}

		private:
			Lod_quality m_lod_quality;
			const bool m_debug;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_QUALITY_ESTIMATOR_H