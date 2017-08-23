#ifndef CGAL_LEVEL_OF_DETAIL_BASE_H
#define CGAL_LEVEL_OF_DETAIL_BASE_H

// STL includes.
#include <string>
#include <iostream>
#include <map>

// CGAL includes.
#include <CGAL/linear_least_squares_fitting_3.h>

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
			
			typedef typename Kernel::Plane_3 Plane;
			typedef typename Kernel::Point_3 Point_3;
			typedef typename Kernel::Point_2 Point_2;

			typedef typename Traits::Building_boundary_selector Building_boundary_selector; // Maybe use a factory here? 
			typedef typename Traits::Building_interior_selector Building_interior_selector;
			typedef typename Traits::Clutter_selector 		    Clutter_selector;
			typedef typename Traits::Ground_selector 		    Ground_selector;

			typedef typename Traits::Vertical_regularizer Vertical_regularizer;
			typedef typename Traits::Ground_projector     Ground_projector;
			typedef typename Traits::Projected            Projected_points;
			typedef typename Traits::Planes               Planes_mapping;
			
			// Custom types, should be changed later (or removed).
			using Log     = CGAL::LOD::Mylog;
			using Index   = typename Container::Index;
			using Indices = std::vector<Index>;

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
				Planes_mapping all_planes;
				auto number_of_planes = m_preprocessor.get_planes(input, all_planes);
				assert(number_of_planes >= 0);

				log.out << "(2) Planes found. Number of planes: " << number_of_planes << std::endl;


				// (3) Split data with respect to 4 different semantic labels.
				Indices building_boundary_idxs, building_interior_idxs, clutter_idxs, ground_idxs;

				m_clutter_selector.select_elements(input, std::back_inserter(clutter_idxs));
				m_ground_selector.select_elements( input, std::back_inserter(ground_idxs));
				m_building_boundary_selector.select_elements(input, std::back_inserter(building_boundary_idxs));
				m_building_interior_selector.select_elements(input, std::back_inserter(building_interior_idxs));

				log.out << "(3) Clutter found. Number of elements: " 			 << clutter_idxs.size() << std::endl;
				log.out << "(.) Ground found. Number of elements: " 			 << ground_idxs.size() << std::endl;
				log.out << "(.) Building boundaries found. Number of elements: " << building_boundary_idxs.size() << std::endl;
				log.out << "(3) Building interior found. Number of elements: "   << building_interior_idxs.size() << std::endl;


				// (4) Create plane from the ground points.
				Plane ground_plane;
				fit_ground_plane(input, ground_idxs, ground_plane);

				log.out << "(4) Ground plane fitted: " << ground_plane << std::endl;


				// (5) Map indices from all detected planes to the ones that are a part of the given facades.
				Planes_mapping building_boundary_planes;
				number_of_planes = m_preprocessor.get_planes(input, building_boundary_idxs, building_boundary_planes);

				log.out << "(5) Planes for building's boundary are found. Number of planes: " << number_of_planes << std::endl;


				// (6) Make all nearly vertical planes in the building's boundary exactly vertical.
				const auto number_of_regularized_planes = m_vertical_regularizer.regularize(building_boundary_planes, input);

				log.out << "(6) Building's nearly vertical planes are regularized. Number of regularized planes: " << number_of_regularized_planes << std::endl;


				// (7) Project all vertical building's boundaries onto the ground plane.
				Projected_points building_boundary_projected;
				const auto number_of_projected_points = m_ground_projector.project(input, building_boundary_planes, ground_plane, building_boundary_projected);

				log.out << "(7) All building's boundary points are projected. Number of projected points: " << number_of_projected_points << std::endl;


				// (END) Save log.
				log.out << "\n\nFINISH EXECUTION";
				log.save("create_lod_0");

				// log.save_ply<Kernel, Container>(input, "regularized", true);
			}
			
		private:
			Traits       m_traits;
			Loader       m_loader;
			Preprocessor m_preprocessor;
			
			Building_boundary_selector m_building_boundary_selector;
			Building_interior_selector m_building_interior_selector;
			Clutter_selector           m_clutter_selector;
			Ground_selector            m_ground_selector;
			Vertical_regularizer	   m_vertical_regularizer;
			Ground_projector 		   m_ground_projector;

			// Not efficient since I need to copy all ground points.
			void fit_ground_plane(const Container &input, const Indices &ground_idxs, Plane &ground_plane) const {

				std::vector<Point_3> tmpGround(ground_idxs.size());
				for (size_t i = 0; i < ground_idxs.size(); ++i) tmpGround[i] = input.point(ground_idxs[i]);

				CGAL::linear_least_squares_fitting_3(tmpGround.begin(), tmpGround.end(), ground_plane, CGAL::Dimension_tag<0>()); 
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_BASE_H