#ifndef CGAL_LEVEL_OF_DETAIL_ENUM_H
#define CGAL_LEVEL_OF_DETAIL_ENUM_H

#if defined(WIN32) || defined(_WIN32) 
#define PS "\\"
#define PN "\r\n"
#else 
#define PS "/" 
#define PN "\n"
#endif

// STL includes.
#include <map>
#include <vector>
#include <unordered_set>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

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
			FACE_BASED_AFFINE,			   // use affine coordinates instead of natural neighbours
			FACE_BASED_COUNT			   // hard to find the proper circle radius
		};

		// Visibility approaches.
		enum class Visibility_approach {
			POINT_BASED, // here we go over all input points
			FACE_BASED   // here we go over all input faces
		};

		// Visibility samplers.
		enum class Visibility_sampler {
			RANDOM_UNIFORM_0,    // a bit faster than UNIFORM_1 but the result is similar, randomized
			RANDOM_UNIFORM_1,    // see above
			UNIFORM_SUBDIVISION, // determenistic sampler based on midpoint subdivision
			BARYCENTRE			 // barycentre of the given triangle
		};

		// Type of the building's boudnary.
		enum class Building_boundary_type {
			ORIENTED,  // clockwise-oriented boundary
			UNORIENTED // not oriented boundary, which is the set of segments
		};

		// Face info class.
		template<typename FT>
		class My_face_info {

		public:
			FT in = FT(1) / FT(2); 				  			 // visibility label (in - inside) or (out - outside); or alternatively label A - everything that is bigger than 0.5, label B < 0.5
			CGAL::Color in_color = CGAL::Color(255, 205, 0); // visibility color (in < 1/2 - red), (in = 1/2 - yellow), (in > 1/2 - green)
			
			int bu = -1; 		   						 	   // building's index - (0, 1, 2 etc.) where (-1 means not a building)
			CGAL::Color bu_color = CGAL::Color(169, 169, 169); // building's color - random color is used per building
		};

		// Vertex info class.
		template<typename Label>
		class My_vertex_info {

		public:
			Label label = Label::CLUTTER;
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

			bool is_oriented = true;    		// flag to check if the computed boundary is oriented or not, see Building_boundary_type above
			std::unordered_set<int> neighbours; // indices of all neighbouring buildings of the given building
		};

		// Type of the roof fitter.
		enum class Roof_fitter_type {
			MIN, // fit data to the minimum height
			AVG, // fit data to the average height
			MAX  // fit data to the maximum height
		};

		// Main test data.
		enum class Main_test_data_type {
			BASIC,           // basic data set from the Loader_stub class.
			COMPLEX,         // sketch up generated simple data set with square buildings
			PARIS,           // half of the paris real data set
			P10,             // p10 data set from the original LOD paper of Yannick
			PARIS_FULL,      // full paris data set
			PARIS_ETH,       // half paris data set classified with ETH random forest
			PARIS_FULL_ETH,  // full paris data set classified with ETH random forest
			RESIDENT_TILE_1, // three different residential tiles
			RESIDENT_TILE_2,
			RESIDENT_TILE_3,
			PARIS_TILE_1, 	 // two different Paris tiles
			PARIS_TILE_2,
			PARIS_BIG        // 1 km - 9 tiles stitched together
		};

		// Nearest neighbour search.
		enum class Neighbour_search_type {
			KNN,    // use k nearest neighbours
			CIRCLE, // use all points from the circle of the given radius
			SQUARE  // use all points from the square of the given radius
		};

		// Fitter type used in thinning.
		enum class Thinning_fitter_type { 
			LINE // fit to the line
		};

		// New point type used in the grid simplify algorithm.
		enum class Grid_new_point_type { 
			CENTROID,   // centroid of the cell
			BARYCENTRE, // barycentre of the given set of samples
			CLOSEST     // point closest to the barycentre above
		};

		// Thinning type.
		enum class Thinning_type {
			NAIVE,	// naive thinning where we project all points onto a line
			COMPLEX // a complex version, where we perform many different optimization steps and preserve features
		};

		// Thinning scale type.
		enum class Thinning_scale_type {
			FIXED,			// fixed manually set scale
			ADAPTIVE_FIXED, // automatically chosen scale
			PROGRESSIVE		// scale changes progressively
		};

		// Structuring corner algorithm.
		enum class Structuring_corner_algorithm {
			NO_CORNERS,         // we do not insert any corners
			GRAPH_BASED,	    // in this algorithm, we build an adjacency graph and insert corners based on this graph
			INTERSECTION_BASED, // in this algorithm, we intersect all segments and insert the best intersections
			NO_T_CORNERS		// we do not insert T-like corners
		};

		// Structuring adjacency threshold method.
		enum class Structuring_adjacency_threshold_method {
			LOCAL, // internal local epsilon is chosen
			GLOBAL // user-defined value is chosen
		};

		// Method to estimate normals in 2D region growing.
		enum class Region_growing_normal_estimation {
			PROJECTED, // project exact normals if they exist
			LOCAL      // estimate normals using PCA
		};

		// Quality data type.
		enum class Quality_data_type { 	
			DST_MIN,   // distortion types
			DST_AVG, 
			DST_MAX,
			DST_ROOFS,
			DST_WALLS,
			CMP_ROOFS, // complexity metric
			CMP_WALLS,
			CMP,
			COV_ROOFS, // coverage metric
			COV_WALLS,
			COV,       
		};

	} // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_ENUM_H