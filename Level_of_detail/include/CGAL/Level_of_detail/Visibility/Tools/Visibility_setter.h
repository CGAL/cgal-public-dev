#ifndef CGAL_LEVEL_OF_DETAIL_VISIBILITY_SETTER_H
#define CGAL_LEVEL_OF_DETAIL_VISIBILITY_SETTER_H

namespace CGAL {

	namespace Level_of_detail {

		class Visibility_setter {

		public:
			template<class Visibility_map, class Facets_range>
			void set_labels(const Visibility_map &visibility_map, Facets_range &facets_range) const {
				using Facets_iterator = typename Facets_range::iterator;

				for (Facets_iterator facet = facets_range.begin(); facet != facets_range.end(); ++facet)
					put(visibility_map, *facet);
			}
		};
	
	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_VISIBILITY_SETTER_H