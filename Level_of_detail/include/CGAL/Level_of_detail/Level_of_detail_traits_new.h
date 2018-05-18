#ifndef CGAL_LEVEL_OF_DETAIL_TRAITS_H
#define CGAL_LEVEL_OF_DETAIL_TRAITS_H

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel>
		struct Level_of_detail_traits {

			using Kernel = InputKernel;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_TRAITS_H