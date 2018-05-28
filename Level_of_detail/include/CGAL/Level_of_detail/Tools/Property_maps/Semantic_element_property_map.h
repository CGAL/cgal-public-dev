#ifndef CGAL_LEVEL_OF_DETAIL_SEMANTIC_ELEMENT_PROPERTY_MAP_H
#define CGAL_LEVEL_OF_DETAIL_SEMANTIC_ELEMENT_PROPERTY_MAP_H

// LOD includes.
#include <CGAL/Level_of_detail/Enumerations.h>

namespace CGAL {

	namespace Level_of_detail {

		template<typename KeyType, class LabelMap>
		class Semantic_element_property_map {
			
        public:
            using Label_map = LabelMap;

            Semantic_element_property_map(const LabelMap &label_map) : 
            m_label_map(label_map) 
            { }

            inline const Label_map& label_map() const {
                return m_label_map;
            }

		private:
            using ValueType = Semantic_label;

            using key_type   = KeyType;
            using value_type = ValueType;
            using reference  = const ValueType&;
            
            using Self = Semantic_element_property_map<KeyType, LabelMap>;

            friend value_type get(const Self &self, const key_type &key) {                
                const Label_map &label_map = self.label_map();

                const auto label = get(label_map, key);
                const size_t label_value = static_cast<size_t>(label);

                switch (label_value) {
                    case 0: // ground
                        return Semantic_label::GROUND;

                    case 1: // facade
                        return Semantic_label::BUILDING_BOUNDARY;

                    case 2: // roof
                        return Semantic_label::BUILDING_INTERIOR;

                    default: // everything else
                        return Semantic_label::UNASSIGNED;
                }
            }

            friend void put(const Self &, key_type &, const value_type &) { }
            const Label_map &m_label_map;
		};

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_SEMANTIC_ELEMENT_PROPERTY_MAP_H