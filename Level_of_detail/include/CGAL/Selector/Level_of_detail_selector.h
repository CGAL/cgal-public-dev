#ifndef CGAL_LEVEL_OF_DETAIL_SELECTOR_H
#define CGAL_LEVEL_OF_DETAIL_SELECTOR_H

// STL includes.
#include <vector>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class SelectionStrategy>
		class Level_of_detail_selector {

		public:
			typedef KernelTraits      Traits;
			typedef SelectionStrategy Strategy;

			template<class InputContainer, class OutputIterator>
			void select_elements(const InputContainer &input, OutputIterator output) {

				using Const_iterator = typename InputContainer::const_iterator;
				m_strategy.set_input(input);

				for (Const_iterator it = input.begin(); it != input.end(); ++it) {

					const int elementIndex = static_cast<int>(*it);
					if (m_strategy.satisfies_condition(elementIndex)) {
						
						 *output = *it;
						++output;
					}
				}
			}

		private:
			Strategy m_strategy;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SELECTOR_H