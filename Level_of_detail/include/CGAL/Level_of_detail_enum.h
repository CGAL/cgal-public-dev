#ifndef CGAL_LEVEL_OF_DETAIL_ENUM_H
#define CGAL_LEVEL_OF_DETAIL_ENUM_H

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
		};

	} // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_ENUM_H