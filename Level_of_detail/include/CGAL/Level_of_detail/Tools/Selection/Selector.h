#ifndef CGAL_LEVEL_OF_DETAIL_SELECTOR_H
#define CGAL_LEVEL_OF_DETAIL_SELECTOR_H

#include <iostream>

namespace CGAL {

	namespace Level_of_detail {

		template<class SelectionStrategy>
		class Selector {

		public:
			using Selection_strategy = SelectionStrategy;

			template<class Input_range, class Label_map, class Output_iterator>
			void select_elements(const Input_range &input_range, const Label_map &label_map, Output_iterator output) {

				using Const_range_iterator = typename Input_range::const_iterator;
				for (Const_range_iterator cr_it = input_range.begin(); cr_it != input_range.end(); ++cr_it) {
					
					if (m_selection_strategy.satisfies_condition(*cr_it, label_map)) {
						*output = cr_it; ++output;
					}
				}
			}

		private:
			Selection_strategy m_selection_strategy;
		};
	
	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_SELECTOR_H