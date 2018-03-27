#ifndef CGAL_LEVEL_OF_DETAIL_LOD2_H
#define CGAL_LEVEL_OF_DETAIL_LOD2_H

// STL includes.
#include <map>
#include <string>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/Random.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/number_utils.h>
#include <CGAL/utils_classes.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>


// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Utils/Level_of_detail_utils.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputHDS, class InputCDT, class InputBuildings, class FacetHandle>
		class LOD2_builder : public Modifier_base<InputHDS> {
		
		public:
			typedef KernelTraits 	 Kernel;
    		typedef InputHDS 		 HDS;
    		typedef InputCDT 		 CDT;
    		typedef InputBuildings 	 Buildings;
    		typedef FacetHandle 	 Color_facet_handle;
    		typedef FacetHandle 	 Facet_handle;

    		typedef typename HDS::Vertex   Vertex;
    		typedef typename Vertex::Point Point;

    		using FT = typename Kernel::FT;

			using Color 	   = CGAL::Color;
			using Facet_colors = std::map<Color_facet_handle, Color>;
			using Builder 	   = CGAL::Polyhedron_incremental_builder_3<HDS>;
    		
    		using Log = CGAL::LOD::Mylog;
    		using Building_iterator = typename Buildings::iterator;
    		using Ground = std::vector<Point>;
        };

        template<class KernelTraits, class InputCDT, class InputBuildings, class OutputMesh>
		class Level_of_detail_lod2 {

		public:
			typedef KernelTraits   Kernel;
			typedef InputCDT 	   CDT;
			typedef InputBuildings Buildings;
			typedef OutputMesh 	   Mesh;

            using HDS = typename Mesh::HalfedgeDS;

            using Mesh_facet_handle = typename Mesh::Facet_const_handle;
            using Mesh_builder 		= LOD2_builder<Kernel, HDS, CDT, Buildings, Mesh_facet_handle>;
			using Mesh_facet_colors = typename Mesh_builder::Facet_colors;

			using Point  = typename Mesh_builder::Point;
			using Ground = typename Mesh_builder::Ground;

            Level_of_detail_lod2(const CDT &cdt, const Buildings &buildings, const Ground &ground) :
            m_cdt(cdt), m_buildings(buildings), m_ground(ground) { }

            void reconstruct(Mesh &mesh, Mesh_facet_colors &mesh_facet_colors) const {
                std::cout << std::endl << "to be added" << std::endl;
            }

        private:
            const CDT       &m_cdt;
            const Buildings &m_buildings;
            const Ground    &m_ground;
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_LOD2_H