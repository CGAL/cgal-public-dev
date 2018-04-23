#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ENVELOPE_PLANE_ASSOCIATER_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ENVELOPE_PLANE_ASSOCIATER_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputBuilding, class InputDiagram>
		class Level_of_detail_building_envelope_plane_associater {
            
        public:
            typedef KernelTraits  Kernel;
            typedef InputBuilding Building;
            typedef InputDiagram  Diagram;

            using Roof              = typename Building::Roof;
            using Associated_planes = typename Roof::Associated_planes;

            using Face_iterator    = typename Diagram::Face_const_iterator;
            using Surface_iterator = typename Diagram::Surface_const_iterator;
            
            Level_of_detail_building_envelope_plane_associater() { }

            void find_associated_planes(const Face_iterator &fit, Associated_planes &associated_planes) const {
                
                associated_planes.clear();
                for (Surface_iterator sit = fit->surfaces_begin(); sit != fit->surfaces_end(); ++sit) {
                    
                    const size_t index = sit->data().index;
                    associated_planes.push_back(index);
                }
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ENVELOPE_PLANE_ASSOCIATER_H