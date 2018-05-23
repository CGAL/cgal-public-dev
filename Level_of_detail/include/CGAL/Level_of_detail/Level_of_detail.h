#ifndef CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_H
#define CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_H

// STL includes.
#include <iostream>

// LOD includes.
#include <CGAL/Level_of_detail/Level_of_detail_parameters.h>
#include <CGAL/Level_of_detail/Level_of_detail_enumerations.h>
#include <CGAL/Level_of_detail/Level_of_detail_data_structure.h>

namespace CGAL {

	namespace Level_of_detail {

		namespace LOD = CGAL::Level_of_detail;

		template<class InputKernel, class InputRange, class InputNormalMap, class InputLabelMap>
		class Level_of_detail_reconstruction {

		public:
			using Kernel      = InputKernel;
			using Input_range = InputRange;
			using Normal_map  = InputNormalMap;
			using Label_map   = InputLabelMap;

			using FT = typename Kernel::FT;

			using Parameters 	 = LOD::Level_of_detail_parameters<FT>;
			using Data_structure = LOD::Level_of_detail_data_structure<Kernel, Input_range, Normal_map, Label_map>;

			Level_of_detail_reconstruction(const Input_range &input_range, const Normal_map &normal_map, const Label_map &label_map, const Parameters &parameters) :
			m_data_structure(input_range, normal_map, label_map),
			m_parameters(parameters)
			{ }

			void build() {

				if (m_parameters.verbose()) std::cout << std::endl << "... building LOD data ..." << std::endl << std::endl;
			}

			void get_lod0() {

				if (m_parameters.verbose()) std::cout << "* constructing LOD0 ..." << std::endl;
			}

			void get_lod1() {

				if (m_parameters.verbose()) std::cout << "* constructing LOD1 ..." << std::endl;
			}

			void get_lod2() {

				if (m_parameters.verbose()) std::cout << "* constructing LOD2 ..." << std::endl;
			}

		private:
			const Data_structure m_data_structure;
			const Parameters &m_parameters;
		};
	
	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_H