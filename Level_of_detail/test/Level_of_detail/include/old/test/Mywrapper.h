#ifndef CGAL_LEVEL_OF_DETAIL_MYWRAPPER_H
#define CGAL_LEVEL_OF_DETAIL_MYWRAPPER_H

// STL includes.
#include <cassert>

// Local includes.
#include "../Level_of_detail.h"
#include "terminal/Myterminal_parser.h"

namespace CGAL {

	namespace Level_of_detail {

		template<class LodTraits>
		class Mywrapper {

		public:
			using Kernel = typename LodTraits::Kernel;
			using FT 	 = typename Kernel::FT;

			using Lod_base       = CGAL::Level_of_detail::Level_of_detail_base<LodTraits>;
			using Lod_parameters = CGAL::Level_of_detail::Myterminal_parser<FT>;
			using Params         = char**;

			Mywrapper(const int num_params, const Params params) {

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

#endif // CGAL_LEVEL_OF_DETAIL_MYWRAPPER_H