#ifndef CGAL_LEVEL_OF_DETAIL_LOADER_H
#define CGAL_LEVEL_OF_DETAIL_LOADER_H

// STL includes.
#include <string>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class OutputContainer>
		class Level_of_detail_loader {
		
		public:
			typedef KernelTraits    Traits;
			typedef OutputContainer Container;

			Level_of_detail_loader(Traits traits = Traits()) : m_traits(traits) { }

			virtual void get_data(const std::string &, Container &) const {

				// To be implemented later!
			}

		private:
			Traits m_traits;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_LOADER_H