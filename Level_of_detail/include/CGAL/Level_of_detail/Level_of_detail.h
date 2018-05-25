#ifndef CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_H
#define CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_H

// STL includes.
#include <iostream>

// LOD includes.
#include <CGAL/Level_of_detail/Level_of_detail_include.h>

namespace CGAL {

	namespace Level_of_detail {

		namespace LOD = CGAL::Level_of_detail;

		template<class InputKernel, class InputRange, class InputPointMap, class InputNormalMap, class InputLabelMap>
		class Level_of_detail_reconstruction {

		public:
			using Kernel      = InputKernel;
			using Input_range = InputRange;
			using Point_map   = InputPointMap;
			using Normal_map  = InputNormalMap;
			using Label_map   = InputLabelMap;

			using FT = typename Kernel::FT;

			using Parameters 	 = LOD::Level_of_detail_parameters<FT>;
			using Data_structure = LOD::Level_of_detail_data_structure<Kernel, Input_range, Point_map, Normal_map, Label_map>;

			// Default template classes.

			// Selection.
			using Default_ground_selection_strategy 		   = LOD::Ground_selection_strategy;
			using Default_building_boundary_selection_strategy = LOD::Building_boundary_selection_strategy;
			using Default_building_interior_selection_strategy = LOD::Building_interior_selection_strategy;

			using Default_ground_selector 			 = LOD::Selector<Default_ground_selection_strategy>;
			using Default_building_boundary_selector = LOD::Selector<Default_building_boundary_selection_strategy>;
			using Default_building_interior_selector = LOD::Selector<Default_building_interior_selection_strategy>;

			// Add here other default template classes.

			Level_of_detail_reconstruction(const Input_range &input_range, const Point_map &point_map, const Normal_map &normal_map, const Label_map &label_map, const Parameters &parameters) :
			m_data_structure(input_range, point_map, normal_map, label_map),
			m_parameters(parameters) { 
				
				CGAL_assertion(input_range.size() != 0);
			}

			//////////////////////////////
			// Functions to be documented!

			void build() {

				if (m_parameters.verbose()) std::cout << std::endl << "... building LOD data ..." << std::endl << std::endl;

				// * Split semantic labels into separate groups.
				split_semantic_labels();
			}

			void get_lod0() {

				if (m_parameters.verbose()) std::cout << "* constructing LOD0 ..." << std::endl;
			}

			void get_lod1() {

				if (m_parameters.verbose()) std::cout << "* constructing LOD1 ..." << std::endl;
			}

			template<
			class Ground_selector = Default_ground_selector, 
			class Building_boundary_selector = Default_building_boundary_selector, 
			class Building_interior_selector = Default_building_interior_selector
			>
			void split_semantic_labels() {

				if (m_parameters.verbose()) std::cout << "* splitting semantic labels" << std::endl;

				// For the moment we find only ground and building points and save their indices.

				Ground_selector 		   ground_selector;
				Building_boundary_selector building_boundary_selector;
				Building_interior_selector building_interior_selector;

				ground_selector.select_elements(m_data_structure.input(), m_data_structure.label_map(), std::back_inserter(m_data_structure.ground_indices()));
				building_boundary_selector.select_elements(m_data_structure.input(), m_data_structure.label_map(), std::back_inserter(m_data_structure.building_boundary_indices()));
				building_interior_selector.select_elements(m_data_structure.input(), m_data_structure.label_map(), std::back_inserter(m_data_structure.building_interior_indices()));
			}

			//////////////////////////////////
			// Functions to be not documented!

			inline const Data_structure& get_internal_data_structure() const {
				return m_data_structure;
			}

			// to be added!

		private:
			Data_structure m_data_structure;
			const Parameters &m_parameters;
		};
	
	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_H