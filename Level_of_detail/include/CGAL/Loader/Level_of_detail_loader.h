#ifndef CGAL_LEVEL_OF_DETAIL_LOADER_H
#define CGAL_LEVEL_OF_DETAIL_LOADER_H

// STL includes.
#include <string>

namespace CGAL {

	namespace LOD {

		template<class Traits, class OutputContainer>
		class Level_of_detail_loader {
		
		public:
			typedef OutputContainer Container;

			Level_of_detail_loader(Traits traits = Traits()) : m_traits(traits) { }

			virtual void load(const std::string &, Container &) const {

				// To be implemented later!
			}

		private:
			Traits m_traits;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_LOADER_H