#ifndef CGAL_LEVEL_OF_DETAIL_QUALITY_H
#define CGAL_LEVEL_OF_DETAIL_QUALITY_H

// STL includes.
#include <map>
#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>

// New CGAL includes.
#include <CGAL/Tools/Level_of_detail_parameters.h>
#include <CGAL/Base/Level_of_detail_base.h>
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class LodTraits>
		class Level_of_detail_quality {

		public:
			using Kernel = typename LodTraits::Kernel;
			using FT 	 = typename Kernel::FT;

			using Lod_parameters = CGAL::LOD::Level_of_detail_parameters<FT>;
			using Lod_base 		 = CGAL::LOD::Level_of_detail_base<LodTraits>;
			using Params 		 = char**;

			using Data = std::vector<FT>;
			using Lod_complexity = typename Lod_base::Lod_complexity;
			using Lod_distortion = typename Lod_base::Lod_distortion;

			typedef typename Lod_parameters::Input_parameters Parameters;
			
			Level_of_detail_quality(const int num_params, const Params params) : m_num_runs(0), m_debug(false) { 

				Lod_parameters lod_parameters(num_params, params);

				m_lod_base.set_optimal_configuration();
				m_lod_base.set_user_defined_parameters(lod_parameters);

				m_lod_base.create_lods();
				update_number_of_runs(lod_parameters);
			}

			void compute_xy_data(const bool test = false) {
				if (test) compute_test_data();
				else compute_data();
			}

			const Data &retreive_x_data() const {
				assert(!m_x_data.empty());
				return m_x_data;
			}

			const Data &retreive_y_data(const Distortion_fitting_type type) const {
				
				switch (type) {
					case Distortion_fitting_type::MIN:
						assert(!m_min_y_data.empty());
						return m_min_y_data;

					case Distortion_fitting_type::AVG:
						assert(!m_avg_y_data.empty());
						return m_avg_y_data;

					case Distortion_fitting_type::MAX:
						assert(!m_max_y_data.empty());
						return m_max_y_data;

					default:
						assert(!"Wrong fitting type!");
						return m_avg_y_data;
				}
			}

		private:
			size_t m_num_runs;
			Lod_base m_lod_base;

			Data m_x_data;
			Data m_min_y_data, m_avg_y_data, m_max_y_data;

			const bool m_debug;
			
			std::shared_ptr<Lod_complexity> m_lod_complexity;
			std::shared_ptr<Lod_distortion> m_lod_distortion;

			void update_number_of_runs(const Lod_parameters &lod_parameters) {
				const Parameters &parameters = lod_parameters.get();
				add_val_parameter("-num_runs", m_num_runs, parameters);
			}

			template<typename Scalar>
			void add_val_parameter(const std::string &parameter_name, Scalar &variable_value, const Parameters &parameters) {
				
				if (!does_parameter_exist(parameter_name, parameters)) return;
				const std::string parameter_value = parameters.at(parameter_name);

				if (parameter_value != "default")
					variable_value = static_cast<Scalar>(std::stod(parameter_value.c_str()));

				std::cout << parameter_name << " : " << variable_value << std::endl;
			}

			bool does_parameter_exist(const std::string &parameter_name, const Parameters &parameters) {
				
				for (typename Parameters::const_iterator param = parameters.begin(); param != parameters.end(); ++param)
					if ((*param).first == parameter_name) return true;

				return false;
			}

			void compute_test_data() {

				const size_t size = 1000;
				clear_and_resize(size);

				for (size_t i = 0; i < m_x_data.size(); ++i) {
					
					m_x_data[i] = static_cast<FT>(i) / static_cast<FT>(m_x_data.size());
					m_min_y_data[i] = m_x_data[i] / FT(8);
					m_avg_y_data[i] = m_x_data[i] * m_x_data[i];
					m_max_y_data[i] = m_x_data[i] * FT(8) / FT(7);
				}
			}

			void compute_data() {
				
				const size_t size = m_num_runs * 2 + 1;
				clear_and_resize(size);

				const FT init_x = m_lod_base.get_scale();
				compute_initial_x_data(init_x);
				compute_final_data();
			}

			void clear_and_resize(const size_t size) {
				m_x_data.clear();
				m_x_data.resize(size);

				m_min_y_data.clear(); m_avg_y_data.clear(); m_max_y_data.clear();
				m_min_y_data.resize(m_x_data.size()); m_avg_y_data.resize(m_x_data.size()); m_max_y_data.resize(m_x_data.size());
			}

			void compute_initial_x_data(const FT init_x) {
				
				assert(init_x > FT(0));
				assert(!m_x_data.empty());

				const FT h = init_x / (static_cast<FT>(m_num_runs) * FT(2));

				FT count = static_cast<FT>(m_num_runs);
				for (size_t i = 0; i < m_num_runs; ++i, count -= FT(1)) m_x_data[i] = init_x - h * count;
				m_x_data[m_num_runs] = init_x;

				count = FT(1);
				for (size_t i = m_num_runs + 1; i < m_num_runs * 2 + 1; ++i, count += FT(1)) m_x_data[i] = init_x + h * count;

				if (m_debug) for (size_t i = 0; i < m_x_data.size(); ++i) std::cout << m_x_data[i] << std::endl;
			}

			void compute_final_data() {

				Data tmp_data(m_x_data.size());
				for (size_t i = 0; i < m_x_data.size(); ++i)
					tmp_data[i] = compute_complexity_and_distortion(i, m_x_data[i]);
				m_x_data = tmp_data;
			}

			FT compute_complexity_and_distortion(const size_t index, const FT scale) {

				m_lod_base.estimate_parameters(false);
				m_lod_base.set_scale(scale);
				m_lod_base.create_lods();

				m_lod_complexity = m_lod_base.get_lod_complexity_ptr();
				m_lod_distortion = m_lod_base.get_lod_distortion_ptr();

				const FT complexity = m_lod_complexity->get();
				const FT min_distortion = m_lod_distortion->get(Distortion_fitting_type::MIN);
				const FT avg_distortion = m_lod_distortion->get(Distortion_fitting_type::AVG);
				const FT max_distortion = m_lod_distortion->get(Distortion_fitting_type::MAX);

				m_min_y_data[index] = min_distortion;
				m_avg_y_data[index] = avg_distortion;
				m_max_y_data[index] = max_distortion;

				return complexity;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_QUALITY_H