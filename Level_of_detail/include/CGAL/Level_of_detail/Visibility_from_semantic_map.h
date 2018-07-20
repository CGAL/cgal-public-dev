#ifndef CGAL_LEVEL_OF_DETAIL_VISIBILITY_FROM_SEMANTIC_MAP_H
#define CGAL_LEVEL_OF_DETAIL_VISIBILITY_FROM_SEMANTIC_MAP_H

#include <CGAL/Level_of_detail/Enumerations.h>

namespace CGAL {

namespace Level_of_detail {
  
template <typename SemanticMap>
struct Visibility_from_semantic_map
{
  typedef typename boost::property_traits<SemanticMap>::key_type key_type;
  typedef double value_type;
  typedef double reference;
  typedef boost::readable_property_map_tag category;

  SemanticMap semantic_map;

  Visibility_from_semantic_map () { }
  Visibility_from_semantic_map (SemanticMap semantic_map) : semantic_map (semantic_map) { }

  friend value_type get (const Visibility_from_semantic_map& map, const key_type& key)
  {
    Semantic_label l = get (map.semantic_map, key);
    
    if (l == Semantic_label::BUILDING_INTERIOR)
      return 1.;
    if (l == Semantic_label::BUILDING_BOUNDARY)
      return 0.5;

    return 0.; // ground, unassigned, vegetation
  }
};

} // namespace CGAL

} // namespace Level_of_detail


#endif // CGAL_LEVEL_OF_DETAIL_VISIBILITY_FROM_SEMANTIC_MAP_H

