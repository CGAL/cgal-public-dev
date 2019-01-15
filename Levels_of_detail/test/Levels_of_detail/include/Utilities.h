#ifndef CGAL_LOD_UTILITIES_H
#define CGAL_LOD_UTILITIES_H

// STL includes.
#include <map>
#include <memory>
#include <string>

// CGAL includes.
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Levels_of_detail/enumerations.h>

namespace CGAL {

  namespace Levels_of_detail {

    template<class Point_set>
    struct Semantic_map_from_labels {
      
    public:

      using key_type = typename Point_set::Index;
      using value_type = Semantic_label;
      using reference = value_type;
      using category = boost::readable_property_map_tag;

      using Map_l2sl = std::map<int, Semantic_label>;
      using Map_l2sl_ptr = std::shared_ptr<Map_l2sl>;

      using Label_map = typename Point_set:: template Property_map<int>;

      Point_set *m_point_set;
      Label_map label_map;
      Map_l2sl_ptr map_l2sl;
      
      Semantic_map_from_labels() { }

      Semantic_map_from_labels(Point_set *point_set) : 
      m_point_set(point_set), 
      map_l2sl(new Map_l2sl()) {

        label_map = m_point_set->template property_map<int>("label").first;
      }

      friend value_type get(const Semantic_map_from_labels &semantic_map, const key_type &key) {

        const int label = semantic_map.label_map[key];
        const auto found = semantic_map.map_l2sl->find(label);

        if (found == semantic_map.map_l2sl->end())
          return Semantic_label::UNASSIGNED;
        return found->second;
      }

    }; // Semantic_map_from_labels

  } // Levels_of_detail

} // CGAL

#endif // CGAL_LOD_UTILITIES_H
