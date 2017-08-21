#ifndef CGAL_LEVEL_OF_DETAIL_PREPROCESSOR_H
#define CGAL_LEVEL_OF_DETAIL_PREPROCESSOR_H

// CGAL includes.
#include <CGAL/Point_set_3.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits>
		class Level_of_detail_preprocessor {

		public:
			typedef KernelTraits Traits;

			template<class Container, class Planes>
			auto get_planes(const Container &input, Planes &planes) const {

				auto number_of_planes= -1;

				using Index          = int;
				using Const_iterator = typename Container::const_iterator;
				using Index_map      = typename Container:: template Property_map<Index>;

				Index_map indices;
				boost::tie(indices,  boost::tuples::ignore) = input. template property_map<Index>("index");
				
				for (Const_iterator it = input.begin(); it != input.end(); ++it)
					if (indices[*it] >= 0) 
						planes[indices[*it]].push_back(*it);

				number_of_planes = planes.size();

				return number_of_planes;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_PREPROCESSOR_H