#ifndef CGAL_LEVEL_OF_DETAIL_BASE_H
#define CGAL_LEVEL_OF_DETAIL_BASE_H

// STL includes.
#include <string>
#include <iostream>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>

namespace CGAL {

	namespace LOD {

		// Facade class that accumulates all necessary objects and operations
		// related to the level of detail (LOD) reconstruction.
		template<class LodTraits>
		class Level_of_detail_base
		{

		public:
			typedef LodTraits 				   Traits;
			typedef typename Traits::Kernel    Kernel;
			typedef typename Traits::Loader    Loader;
			
			typedef typename Loader::Container Container;
			using Log = CGAL::Mylog;

			Level_of_detail_base(Traits traits = Traits()) : m_traits(traits) { }

			template<class OutputIterator>
			void create_lod_0(const std::string &filePath, OutputIterator) const {

				Container point_set;

				m_loader.load(filePath, point_set);

				Log log;
				log.out << point_set.number_of_points() << std::endl;

				log.save("num_points");
			}
			
		private:
			Traits m_traits;
			Loader m_loader;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_BASE_H