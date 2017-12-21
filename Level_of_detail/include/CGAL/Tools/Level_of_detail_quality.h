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

			typedef typename Lod_parameters::Input_parameters Parameters;
			
			Level_of_detail_quality(const int num_params, const Params params) : m_num_runs(2) { 

				Lod_parameters lod_parameters(num_params, params);

				m_lod_base.set_optimal_configuration();
				m_lod_base.set_user_defined_parameters(lod_parameters);

				m_lod_base.create_lods();
				update_number_of_runs(lod_parameters);
			}

			void compute_xy_data() {
				for (size_t i = 0; i < m_num_runs; ++i) m_lod_base.create_lods();
			}

			std::vector<FT> &retreive_x_data() const {
				assert(!m_x_data.empty());
				return m_x_data;
			}

			std::vector<FT> &retreive_y_data() const {
				assert(!m_y_data.empty());
				return m_y_data;
			}

		private:
			size_t m_num_runs;
			Lod_base m_lod_base;

			std::vector<FT> m_x_data;
			std::vector<FT> m_y_data;

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
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_QUALITY_H