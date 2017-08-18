#ifndef CGAL_LEVEL_OF_DETAIL_BASE_H
#define CGAL_LEVEL_OF_DETAIL_BASE_H

// STL includes.
#include <string>
#include <iostream>
#include <map>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>

namespace CGAL {

	namespace LOD {

		// Facade class that accumulates all necessary objects and operations
		// related to the level of detail (LOD) reconstruction.
		template<class LodTraits>
		class Level_of_detail_base {

		public:
			typedef LodTraits 				      Traits;
			typedef typename Traits::Kernel       Kernel;
			typedef typename Traits::Container    Container;
			typedef typename Traits::Loader       Loader;
			typedef typename Traits::Preprocessor Preprocessor;
			typedef typename Traits::Selector     Selector;     
			
			// Custom types, should be changed later (or removed).
			using Log = CGAL::LOD::Mylog;
			using Planes = std::map<int, std::vector<int> >;

			Level_of_detail_base(Traits traits = Traits()) : m_traits(traits) { } // Do I need to create an instance of these traits here?

			template<class OutputIterator>
			void create_lod_0(const std::string &filePath, OutputIterator &&) const {

				// (START) Create log.
				Log log; log.out << "START EXECUTION\n\n";

				// (1) Read data.
				Container input;
				m_loader.get_data(filePath, input);

				log.out << "(1) Data loaded. Number of points: " << input.number_of_points() << std::endl;

				// (2) Find a set of planes related to the points. Basically here we emulate RANSAC.
				// For each plane we store indices of all points contained in this plane.
				Planes planes;
				const auto number_of_planes = m_preprocessor.get_planes(input, planes);
				assert(number_of_planes >= 0);

				log.out << "(2) Planes found. Number of planes: " << number_of_planes << std::endl;

				// (3) Split data wrt to the labels that define a building.
				// To be implemented!

				// (END) Save log.
				log.out << "\n\nFINISH EXECUTION";
				log.save("create_lod_0");
			}
			
		private:
			Traits       m_traits;
			Loader       m_loader;
			Preprocessor m_preprocessor;
			Selector     m_selector;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_BASE_H