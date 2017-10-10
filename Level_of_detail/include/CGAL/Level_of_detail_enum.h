#ifndef CGAL_LEVEL_OF_DETAIL_ENUM_H
#define CGAL_LEVEL_OF_DETAIL_ENUM_H

// STL includes.
#include <map>
#include <vector>

// CGAL includes.
#include <CGAL/IO/Color.h>

namespace CGAL {

	namespace LOD {

		enum class Structured_label { 
			LINEAR, 
			CORNER, 
			CLUTTER 
		};

		// Visibility labels.
		enum class Visibility_label { 
			IN,  	// inside a triangle
			OUT, 	// outside a triangle
			UNKNOWN // 50/50 inside outside
		};

		// Visibility methods.
		enum class Visibility_method {
			POINT_BASED_CLASSIFICATION,    // this is a repeatable method
			FACE_BASED_BARYCENTRIC,        // does not work
			FACE_BASED_NATURAL_NEIGHBOURS, // every time can give different result since it is randomized
			FACE_BASED_COUNT			   // 
		};

		// Visibility approaches.
		enum class Visibility_approach {
			POINT_BASED, // here we go over all input points
			FACE_BASED   // here we go over all input faces
		};

		// Visibility samplers.
		enum class Visibility_sampler {
			UNIFORM_0, // a bit faster than UNIFORM_1 but the result is similar
			UNIFORM_1
		};

		// Face info class.
		template<typename FT>
		class My_face_info {

		public:
			FT in = FT(1) / FT(2); 				  			 // visibility label (in - inside) or (out - outside)
			CGAL::Color in_color = CGAL::Color(255, 204, 0); // visibility color (in < 1/2 - red), (in = 1/2 - yellow), (in > 1/2 - green)
			
			int bu = -1; 		   						 	   // building's index - (0, 1, 2 etc.) where (-1 means not a building)
			CGAL::Color bu_color = CGAL::Color(169, 169, 169); // building's color - random color is used per building
		};

		// Vertex info class.
		template<typename Label>
		class My_vertex_info {

		public:
			Label label;
			CGAL::Color color = CGAL::Color(0, 0, 0);
		};

		// Building structure.
		template<class FT, class Vertex_handle, class Face_handle>
		struct Building {

		public:
			FT height 		  = FT(0); 				  // height of the building
			CGAL::Color color = CGAL::Color(0, 0, 0); // color of the building
			
			std::vector< std::vector<Vertex_handle> > 		   				  boundaries; // boundary vertices of the building ordered counterclockwise (may store multiple boundaries)
			std::vector< std::map<Vertex_handle, std::vector<Face_handle> > > wedges;     // all faces adjacent to each boundary vertex above - must be unique face handles
			std::vector<Face_handle>   						   				  faces;	  // all faces that belong to this building
		};

		// Type of the roof fitter.
		enum class Roof_fitter_type {
			MIN, // fit data to the minimum height
			AVG, // fit data to the average height
			MAX  // fit data to the maximum height
		};

		// Main test data.
		enum class Main_test_data_type {
			BASIC,   // basic data set from the Loader_stub class.
			COMPLEX, // sketch up generated simple data set with square buildings
			PARIS,   // paris real data set
			P10      // p10 data set from the original LOD paper of Yannick
		};

	} // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_ENUM_H