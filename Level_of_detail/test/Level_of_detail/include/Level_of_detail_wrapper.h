#ifndef CGAL_LEVEL_OF_DETAIL_WRAPPER_H
#define CGAL_LEVEL_OF_DETAIL_WRAPPER_H

// STL includes.
#include <cassert>

// CGAL new includes.
#include <CGAL/Level_of_detail.h>

// Local includes.
#include "terminal/Level_of_detail_terminal.h"

namespace CGAL {

	namespace LOD {

		template<class LodTraits>
		class Level_of_detail_wrapper {

		public:
			using Kernel = typename LodTraits::Kernel;
			using FT 	 = typename Kernel::FT;

			using Lod_base       = CGAL::LOD::Level_of_detail_base<LodTraits>;
			using Lod_parameters = CGAL::LOD::Level_of_detail_parameters<FT>;
			using Params         = char**;

			Level_of_detail_wrapper(const int num_params, const Params params) {

				Lod_parameters lod_parameters(num_params, params);

				m_lod_base.set_optimal_configuration();
				m_lod_base.set_user_defined_parameters(lod_parameters);
			}

			void run_lod_pipeline() {
				m_lod_base.create_lods();
			}

		private:
			Lod_base m_lod_base;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_WRAPPER_H