#ifndef CGAL_LEVEL_OF_DETAIL_COMPLEXITY_H
#define CGAL_LEVEL_OF_DETAIL_COMPLEXITY_H

#if defined(WIN32) || defined(_WIN32) 
#define PS "\\"
#else 
#define PS "/"
#endif

// STL includes.
#include <map>
#include <string>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/number_utils.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Selector/Level_of_detail_selector.h>
#include <CGAL/Selector/Level_of_detail_selection_strategy.h>
#include <CGAL/Region_growing/Level_of_detail_planar_region_growing.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class LodReconstruction>
		class Level_of_detail_complexity {

		public:
			typedef KernelTraits   	  Kernel;
			typedef InputContainer    Container;
			typedef LodReconstruction LODS;

			typedef typename Kernel::FT FT;

			typedef typename LodReconstruction::Mesh Mesh;
			using Facet_handle = typename Mesh::Facet_const_handle;

			using Index   = int;
			using Indices = std::vector<Index>;
			using Log 	  = CGAL::LOD::Mylog;
			using Faces   = std::vector< std::vector<Facet_handle> >;

			using Building_region  = std::vector<Facet_handle>;
			using Building_regions = std::vector<Building_region>;
			using Regions 		   = std::vector<Building_regions>;

			typedef CGAL::LOD::Level_of_detail_building_interior<Kernel, Container> Roofs_points_selection_strategy;
			typedef CGAL::LOD::Level_of_detail_building_boundary<Kernel, Container> Walls_points_selection_strategy;
			
			typedef CGAL::LOD::Level_of_detail_selector<Kernel, Roofs_points_selection_strategy> Roofs_points_selector;
			typedef CGAL::LOD::Level_of_detail_selector<Kernel, Walls_points_selection_strategy> Walls_points_selector;

			typedef CGAL::LOD::Level_of_detail_planar_region_growing<Kernel, Mesh, Faces> Planar_region_growing;

			Level_of_detail_complexity(const Container &input, const LODS &lods) : 
			m_input(input), m_lods(lods), 
			m_roofs_complexity(-FT(1)), m_walls_complexity(-FT(1)), m_complexity(-FT(1)),
			m_debug(false) { }

			void estimate() {

				const FT num_roofs = get_number_of_roofs();
				const FT num_walls = get_number_of_walls();

				m_roofs_complexity = num_roofs;
				m_walls_complexity = num_walls;
				m_complexity 	   = m_roofs_complexity + m_walls_complexity;

				if (m_debug) {
					get_roofs_normalization_factor();
					get_walls_normalization_factor();
				}
			}

			FT get_for_roofs() const {
				assert(m_roofs_complexity >= FT(0));
				return m_roofs_complexity;
			}

			FT get_for_walls() const {
				assert(m_walls_complexity >= FT(0));
				return m_walls_complexity;
			}

			FT get() const {
				assert(m_complexity >= FT(0));
				return m_complexity;
			}

		private:
			const Container &m_input;
			const LODS 		&m_lods;
			
			FT 	 m_roofs_complexity, m_walls_complexity, m_complexity;
			bool m_debug;

			FT get_number_of_roofs() const {

				const int num_roofs = m_lods.get_number_of_roofs();
				assert(num_roofs >= 0);

				return static_cast<FT>(num_roofs);
			}

			FT get_number_of_walls() const {

				Faces walls_faces;
				m_lods.get_walls_faces(walls_faces);

				Planar_region_growing planar_region_growing(walls_faces);

				Regions regions;
				planar_region_growing.find_regions(regions);

				const int num_walls = planar_region_growing.get_number_of_regions();
				assert(num_walls >= 0);

				if (m_debug) {
					Log log;
					log.save_mesh_as_ply<Mesh, Regions>(regions, "tmp" + std::string(PS) + "detected_walls");
				}

				return static_cast<FT>(num_walls);
			}

			FT get_roofs_normalization_factor() const {
				
				Indices roofs_point_indices;
				Roofs_points_selector selector;
				selector.select_elements(m_input, std::back_inserter(roofs_point_indices));
	
				if (m_debug) {
					Log roofs_points_saver;
					roofs_points_saver.export_points_using_indices(m_input, roofs_point_indices, "tmp" + std::string(PS) + "roofs_points_for_complexity");
				}

				const int roofs_num_points = roofs_point_indices.size();
				if (m_debug) std::cout << std::endl << "num roofs points " << roofs_num_points;

				return static_cast<FT>(roofs_num_points);
			}

			FT get_walls_normalization_factor() const {			
				
				Indices walls_point_indices;
				Walls_points_selector selector;
				selector.select_elements(m_input, std::back_inserter(walls_point_indices));
	
				if (m_debug) {
					Log walls_points_saver;
					walls_points_saver.export_points_using_indices(m_input, walls_point_indices, "tmp" + std::string(PS) + "walls_points_for_complexity");
				}

				const int walls_num_points = walls_point_indices.size();
				if (m_debug) std::cout << ", num walls points " << walls_num_points << std::endl << std::endl;

				return static_cast<FT>(walls_num_points);
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_COMPLEXITY_H