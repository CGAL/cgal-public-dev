#ifndef CGAL_LEVEL_OF_DETAIL_VISIBILITY_FROM_CLASSIFICATION_PROPERTY_MAP_2_H
#define CGAL_LEVEL_OF_DETAIL_VISIBILITY_FROM_CLASSIFICATION_PROPERTY_MAP_2_H

namespace CGAL {

	namespace Level_of_detail {

		template<typename KeyType, class TreeWrapper, class LabelMap>
		class Visibility_from_classification_property_map_2 {
			
        public:
            using Tree      = TreeWrapper;
            using Label_map = LabelMap;

            Visibility_from_classification_property_map_2(const Tree &tree, const LabelMap &label_map) : 
            m_tree(tree),
            m_label_map(label_map) { }

            inline const Label_map& label_map() const {
                return m_label_map;
            }

            using key_type = KeyType;
            using Self     = Visibility_from_classification_property_map_2<KeyType, TreeWrapper, LabelMap>;
            using Facet    = key_type;

            friend void put(const Self &self, const key_type &facet) {
                self.compute(facet);
            }

        private:
            const Tree      &m_tree;
            const Label_map &m_label_map;

            void compute(const Facet &facet) {

            }
		};

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_VISIBILITY_FROM_CLASSIFICATION_PROPERTY_MAP_2_H