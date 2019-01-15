#ifndef CGAL_LEVELS_OF_DETAIL_H
#define CGAL_LEVELS_OF_DETAIL_H

// LOD includes.
#include <CGAL/Levels_of_detail/internal/Data_structure.h>
#include <CGAL/Levels_of_detail/property_maps.h>

namespace CGAL {

	namespace Levels_of_detail {

    template<typename GeometricTraits,
             typename InputRange,
             typename PointMap,
             typename SemanticMap,
             typename VisibilityMap = Visibility_from_semantic_map<SemanticMap>,
             typename Verbose = CGAL::Tag_false>
		class Levels_of_detail {

		public:

      using Data_structure = internal::Data_structure<
      GeometricTraits, 
      InputRange, 
      PointMap, 
      SemanticMap, 
      VisibilityMap>;

      Levels_of_detail(
        const InputRange &input_range,
        PointMap point_map,
        SemanticMap semantic_map,
        VisibilityMap visibility_map = VisibilityMap()) : 
      m_data_structure(input_range, point_map, semantic_map, visibility_map) {

      }

      void compute_planar_ground() {

      }

    private:
      Data_structure m_data_structure;

    }; // end class

	} // Levels_of_detail

} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_H
