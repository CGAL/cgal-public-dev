#ifndef CGAL_LEVEL_OF_DETAIL_VISIBILITY_FROM_CLASSIFICATION_PROPERTY_MAP_2_H
#define CGAL_LEVEL_OF_DETAIL_VISIBILITY_FROM_CLASSIFICATION_PROPERTY_MAP_2_H

namespace CGAL {

	namespace Level_of_detail {

		template<typename KeyType, class InputKernel, class InputRange, class LabelMap>
		class Visibility_from_classification_property_map_2 {
			
        public:
            using Kernel      = InputKernel;
            using Input_range = InputRange;
            using Label_map   = LabelMap;

            Visibility_from_classification_property_map_2(const Input_range &input_range, const Label_map &label_map) :
            m_input_range(input_range),
            m_label_map(label_map) { }

            using key_type = KeyType;
            using Self     = Visibility_from_classification_property_map_2<key_type, Kernel, Input_range, Label_map>;
            using Facet    = key_type;

            friend void put(const Self &self, key_type &facet) {
                self.compute_shepard_based_label(facet);
            }

        private:
            const Input_range &m_input_range;
            const Label_map   &m_label_map;

            void compute_shepard_based_label(Facet &facet) const {

            }
		};

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_VISIBILITY_FROM_CLASSIFICATION_PROPERTY_MAP_2_H