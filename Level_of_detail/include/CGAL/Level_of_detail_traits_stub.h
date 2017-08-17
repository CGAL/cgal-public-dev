#ifndef CGAL_LEVEL_OF_DETAIL_TRAITS_STUB_H
#define CGAL_LEVEL_OF_DETAIL_TRAITS_STUB_H

// New CGAL includes.
#include <CGAL/Loader/Level_of_detail_loader_stub.h>
#include <CGAL/Preprocessor/Level_of_detail_preprocessor.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class OutputContainer>
		struct Level_of_detail_traits_stub {

			typedef KernelTraits 								   Kernel;
			typedef OutputContainer 							   Container;
			typedef Level_of_detail_loader_stub<Kernel, Container> Loader;
			typedef Level_of_detail_preprocessor<Kernel>           Preprocessor;

		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_TRAITS_STUB_H