#ifndef CGAL_LEVELS_OF_DETAIL_PROPERTY_MAPS_H
#define CGAL_LEVELS_OF_DETAIL_PROPERTY_MAPS_H

// CGAL includes.
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Levels_of_detail/enumerations.h>

namespace CGAL {

  namespace Levels_of_detail {

    template<class SemanticMap>
    struct Visibility_from_semantic_map {

    public:
    
      using key_type = typename boost::property_traits<SemanticMap>::key_type;
      using value_type = double;
      using reference = value_type;
      using category = boost::readable_property_map_tag;

      SemanticMap semantic_map;

      Visibility_from_semantic_map() { }

      Visibility_from_semantic_map(SemanticMap semantic_map) : 
      semantic_map(semantic_map) 
      { }

      friend value_type get(
        const Visibility_from_semantic_map &visibility_map, 
        const key_type &key) {

        const Semantic_label label = get(visibility_map.semantic_map, key);
        
        if (label == Semantic_label::BUILDING_INTERIOR)
          return 1.0;
        if (label == Semantic_label::BUILDING_BOUNDARY)
          return 0.5;

        return 0.0; // ground, unassigned, vegetation
      }
      
    }; // Visibility_from_semantic_map

  } // Levels_of_detail

} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_PROPERTY_MAPS_H
