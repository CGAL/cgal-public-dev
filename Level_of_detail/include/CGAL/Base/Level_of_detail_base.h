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

			typedef typename Traits::BuildingBoundarySelector BuildingBoundarySelector; // Maybe use a factory here? 
			typedef typename Traits::BuildingInteriorSelector BuildingInteriorSelector;
			typedef typename Traits::ClutterSelector 		  ClutterSelector;
			typedef typename Traits::GroundSelector 		  GroundSelector;
			
			// Custom types, should be changed later (or removed).
			using Log     = CGAL::LOD::Mylog;
			using Index   = typename Container::Index;
			using Indices = std::vector<Index>;
			using Planes  = std::map<int, Indices>;

			Level_of_detail_base(Traits traits = Traits()) : m_traits(traits) { } // Do I need to create an instance of these traits here?

			template<class OutputIterator>
			void create_lod_0(const std::string &filePath, OutputIterator &&) {

				// (START) Create log.
				Log log; log.out << "START EXECUTION\n\n\n";


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


				// (3) Split data with respect to 4 different semantic labels.
				Indices building_boundary, building_interior, clutter, ground;

				m_clutter_selector.select_elements(input,           std::back_inserter(clutter));
				m_ground_selector.select_elements( input,           std::back_inserter(ground));
				m_building_boundary_selector.select_elements(input, std::back_inserter(building_boundary));
				m_building_interior_selector.select_elements(input, std::back_inserter(building_interior));

				log.out << "(3) Clutter found. Number of elements: " 			 << clutter.size() << std::endl;
				log.out << "(3) Ground found. Number of elements: " 			 << ground.size() << std::endl;
				log.out << "(3) Building boundaries found. Number of elements: " << building_boundary.size() << std::endl;
				log.out << "(3) Building interior found. Number of elements: "   << building_interior.size() << std::endl;


				// (4) ....


				// (END) Save log.
				log.out << "\n\nFINISH EXECUTION";
				log.save("create_lod_0");
			}
			
		private:
			Traits       m_traits;
			Loader       m_loader;
			Preprocessor m_preprocessor;
			
			BuildingBoundarySelector m_building_boundary_selector;
			BuildingInteriorSelector m_building_interior_selector;
			ClutterSelector          m_clutter_selector;
			GroundSelector           m_ground_selector;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_BASE_H