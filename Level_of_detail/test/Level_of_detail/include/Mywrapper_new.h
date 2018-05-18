#ifndef CGAL_LEVEL_OF_DETAIL_MYWRAPPER_H
#define CGAL_LEVEL_OF_DETAIL_MYWRAPPER_H

// STL includes.
#include <string>

// LOD includes.
#include <CGAL/Level_of_detail/Level_of_detail_new.h>
#include <CGAL/Level_of_detail/Level_of_detail_traits_new.h>

// Local includes.
#include "terminal/Myterminal_parser.h"

namespace CGAL {

	namespace Level_of_detail {

		namespace LOD = CGAL::Level_of_detail;

		template<class InputKernel>
		class Mywrapper {

		public:
			using Kernel = InputKernel;
			
			using Parameters = char**;
			using FT 		 = typename Kernel::FT;
			
			using Terminal_parser = LOD::Myterminal_parser<FT>;

			using LOD_traits = LOD::Level_of_detail_traits<Kernel>;
			using LOD_base   = LOD::Level_of_detail_reconstruction<LOD_traits>;

			Mywrapper(const int num_parameters, const Parameters parameters, const std::string &logs_path)
			: m_terminal_parser(num_parameters, parameters, logs_path) { }

			void run_lod_pipeline() {
				LOD_base lod_base;
				
				lod_base.build_lod0();
				lod_base.build_lod1();
			}

		private:
			Terminal_parser m_terminal_parser;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_MYWRAPPER_H