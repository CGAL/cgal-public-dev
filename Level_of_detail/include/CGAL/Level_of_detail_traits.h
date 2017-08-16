#ifndef CGAL_LEVEL_OF_DETAIL_TRAITS_H
#define CGAL_LEVEL_OF_DETAIL_TRAITS_H

#include <CGAL/Loader/Level_of_detail_loader.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class Container>
		struct Level_of_detail_traits {

			typedef KernelTraits Kernel;
			typedef Level_of_detail_loader<Kernel, Container> Loader;

		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_TRAITS_H