#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_BOUNDARY_SELECTION_STRATEGY_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_BOUNDARY_SELECTION_STRATEGY_H

namespace CGAL {

	namespace Level_of_detail {

		class Building_boundary_selection_strategy {

		public:
			template<class Element_identifier, class Label_map>
			bool satisfies_condition(const Element_identifier &element_identifier, const Label_map &label_map) const {

				using Label = int;
				const Label facade = 1;

				const Label element_label = static_cast<Label>(get(label_map, element_identifier));
				if (element_label == facade) return true;
				
				return false;
			}
		};

	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_BOUNDARY_SELECTION_STRATEGY_H