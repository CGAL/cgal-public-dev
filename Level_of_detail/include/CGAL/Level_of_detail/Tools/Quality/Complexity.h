#ifndef CGAL_LEVEL_OF_DETAIL_COMPLEXITY_H
#define CGAL_LEVEL_OF_DETAIL_COMPLEXITY_H

// STL includes.
#include <map>
#include <string>
#include <vector>

// CGAL includes.
/*
#include <CGAL/number_utils.h> */

// LOD includes.
/*
#include <Region_growing/Faces_based_3/Faces_based_region_growing_3.h> */

namespace CGAL {

	namespace Level_of_detail {

		namespace LOD = CGAL::Level_of_detail;

		template<class InputKernel, class InputContainer, class InputReconstruction>
		class Complexity {

		public:
			using Kernel 		 = InputKernel;
			using Container 	 = InputContainer;
			using Reconstruction = InputReconstruction;

			using FT = typename Kernel::FT;

			/*
			using Mesh 		   = typename Reconstruction::Mesh;
			using Facet_handle = typename Mesh::Facet_const_handle;

			using Index   = int;
			using Indices = std::vector<Index>;
			using Faces   = std::vector< std::vector<Facet_handle> >;

			using Building_region  = std::vector<Facet_handle>;
			using Building_regions = std::vector<Building_region>;
			using Regions 		   = std::vector<Building_regions>;

			typedef LOD::Level_of_detail_building_interior<Kernel, Container> Roofs_points_selection_strategy;
			typedef LOD::Level_of_detail_building_boundary<Kernel, Container> Walls_points_selection_strategy;
			typedef LOD::Level_of_detail_selector<Kernel, Roofs_points_selection_strategy> Roofs_points_selector;
			typedef LOD::Level_of_detail_selector<Kernel, Walls_points_selection_strategy> Walls_points_selector;
			typedef LOD::Level_of_detail_planar_region_growing<Kernel, Mesh, Faces> Planar_region_growing; */

			Complexity(const Container &input, const Reconstruction &reconstruction) : 
			m_input(input),
			m_reconstruction(reconstruction),
			m_roofs_complexity(-FT(1)),
			m_walls_complexity(-FT(1)),
			m_complexity(-FT(1))
			{ }

			/*
			void estimate() {

				const FT num_roofs = get_number_of_roofs();
				const FT num_walls = get_number_of_walls();

				m_roofs_complexity = num_roofs;
				m_walls_complexity = num_walls;
				m_complexity 	   = m_roofs_complexity + m_walls_complexity;
			}

			FT get_for_roofs() const {
				
				CGAL_precondition(m_roofs_complexity >= FT(0));
				return m_roofs_complexity;
			}

			FT get_for_walls() const {
				
				CGAL_precondition(m_walls_complexity >= FT(0));
				return m_walls_complexity;
			}

			FT get() const {

				CGAL_precondition(m_complexity >= FT(0));
				return m_complexity;
			} */

		private:
			const Container 	 &m_input;
			const Reconstruction &m_reconstruction;
			
			FT m_roofs_complexity;
			FT m_walls_complexity;
			FT m_complexity;

			/*
			FT get_number_of_roofs() const {

				const int num_roofs = m_lods.get_number_of_roofs();
				CGAL_precondition(num_roofs >= 0);

				return static_cast<FT>(num_roofs);
			}

			FT get_number_of_walls() const {

				Faces walls_faces;
				m_lods.get_walls_faces(walls_faces);

				Planar_region_growing planar_region_growing(walls_faces);

				Regions regions;
				planar_region_growing.find_regions(regions);

				const int num_walls = planar_region_growing.get_number_of_regions();
				CGAL_precondition(num_walls >= 0);

				return static_cast<FT>(num_walls);
			}

			FT get_roofs_normalization_factor() const {
				
				Indices roofs_point_indices;
				Roofs_points_selector selector;
				selector.select_elements(m_input, std::back_inserter(roofs_point_indices));

				const int roofs_num_points = roofs_point_indices.size();
				return static_cast<FT>(roofs_num_points);
			}

			FT get_walls_normalization_factor() const {			
				
				Indices walls_point_indices;
				Walls_points_selector selector;
				selector.select_elements(m_input, std::back_inserter(walls_point_indices));

				const int walls_num_points = walls_point_indices.size();
				return static_cast<FT>(walls_num_points);
			} */
		};

	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_COMPLEXITY_H