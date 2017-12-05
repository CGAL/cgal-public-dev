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
			CGAL::Color in_color = CGAL::Color(255, 204, 0); // visibility color (in < 1/2 - red), (in = 1/2 - yellow), (in > 1/2 - green)
			
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

		// LOD parameters.
		template<class FT>
		class LOD_parameters {

		private:
			using Params = char**;

		public:
			using Input_parameters = std::map<std::string, std::string>;

			LOD_parameters(const int num_params, const Params params) { 

				// Help.
				show_help(num_params, params);


				// Handle all input parameters.
				Input_parameters input_parameters;
				set_input_parameters(num_params, params, input_parameters);


				// Set here all required parameters.
				std::vector<std::string> required(1);
				required[0] = "-data";


				// Set lod parameters.
				set_lod_parameters(input_parameters, required);


				// Set here all parameters that should not be saved.
				std::vector<std::string> exceptions(1);
				exceptions[0] = "-load_params";


				// Save parameters.
				save_parameters_to_file(input_parameters, exceptions);
			}

			inline Input_parameters& get() {
				return m_lod_parameters;
			}

			inline const Input_parameters& get() const {
				return m_lod_parameters;
			}

		private:
			Input_parameters m_lod_parameters;

			void show_help(const int num_params, const Params params) {
				if (!is_asked_for_help(num_params, params)) return;
	
				print_help();
				exit(EXIT_SUCCESS);
			}

			bool is_asked_for_help(const int num_params, const Params params) {
				for (int i = 0; i < num_params; ++i)
					if (std::strcmp(params[i], "-help") == 0) return true;
				return false;
			}

			void print_help() {

				std::cout << std::endl << "HELP:" << std::endl;

				
				std::cout << std::endl << "EXAMPLE:" << std::endl;
				std::cout << "your terminal name $ ." << std::string(PS) << "lod -data path_to_data" << std::string(PS) << "data_name.ply -other_param_name -other_param_value" << std::endl << std::endl;


				std::cout << std::endl << "REQUIRED PARAMETERS:" << std::endl << std::endl;

				std::cout << 
				"param name: -data" 			   				<< std::endl 	   <<
				"param values: path_to_data" << std::string(PS) << "data_name.ply" << std::endl <<
				"description: path to the file with input data" << std::endl       << std::endl;


				std::cout << std::endl << "OPTIONAL PARAMETERS:" << std::endl << std::endl;

				std::cout << 
				"param name: -help" 	 << std::endl <<
				"description: show help" << std::endl << std::endl;

				std::cout << 
				"param name: -silent" 						   							   << std::endl <<
				"description: supress any intermediate output except for the final result" << std::endl << std::endl;

				std::cout << 
				"param name: -load_params" 					 << std::endl    <<
				"param value: path_to" << std::string(PS)    << "params.lod" << std::endl <<
				"description: load parameters from the file" << std::endl    << std::endl;
			}

			void set_input_parameters(const int num_params, const Params params, Input_parameters &input_parameters) {

				assert(num_params > 0);
				for (int i = 1; i < num_params; ++i) {

					std::string str   = static_cast<std::string>(params[i]);
					auto first_letter = str[0];

					if (first_letter == '-') {
						if (i + 1 < num_params) {

							str = static_cast<std::string>(params[i + 1]);
							first_letter = str[0];

							if (first_letter != '-') input_parameters[params[i]] = params[i + 1];
							else input_parameters[params[i]] = "default";
						} else input_parameters[params[i]] = "default";
					}
				}
			}

			void set_lod_parameters(Input_parameters &input_parameters, const std::vector<std::string> &required) {
				
				if (parameters_should_be_loaded(input_parameters))
					load_parameters_from_file(input_parameters);

				if (are_required_parameters_set(input_parameters, required)) m_lod_parameters = input_parameters;
				else {

					std::cerr << std::endl << "ERROR: Not all required parameters are provided!" << std::endl << std::endl;
					exit(EXIT_FAILURE);
				}
			}

			bool are_required_parameters_set(const Input_parameters &input_parameters, const std::vector<std::string> &required) {

				bool are_all_set = true;
				for (size_t i = 0; i < required.size(); ++i)
					if (!is_required_parameter_set(required[i], input_parameters)) are_all_set = false;
				
				return are_all_set;
			}

			bool is_required_parameter_set(const std::string &parameter_name, const Input_parameters &input_parameters) {
				
				const bool is_set = does_parameter_exist(parameter_name, input_parameters) && !does_parameter_have_default_value(parameter_name, input_parameters);
				if (!is_set) std::cerr << std::endl << parameter_name << " parameter IS REQUIRED!" << std::endl;

				return is_set;
			}

			bool does_parameter_exist(const std::string &parameter_name, const Input_parameters &input_parameters) {
				
				for (Input_parameters::const_iterator param = input_parameters.begin(); param != input_parameters.end(); ++param)
					if ((*param).first == parameter_name) return true;

				return false;
			}

			bool does_parameter_have_default_value(const std::string &parameter_name, const Input_parameters &input_parameters) {

				assert(does_parameter_exist(parameter_name, input_parameters));
				return input_parameters.at(parameter_name) == "default";
			}

			void print_input_parameters(const Input_parameters &input_parameters) {

				std::cout << std::endl << "User defined parameter values (all other values are default): " << std::endl;
				for (Input_parameters::const_iterator param = input_parameters.begin(); param != input_parameters.end(); ++param)
					std::cout << (*param).first << " : " << (*param).second << std::endl;
			}

			bool parameters_should_be_loaded(const Input_parameters &input_parameters) {
				
				if (does_parameter_exist("-load_params", input_parameters)) return true;
				return false;
			}

			void load_parameters_from_file(Input_parameters &input_parameters) {

				const std::string filePath = input_parameters.at("-load_params");
				if (filePath == "default") {
				
					std::cerr << std::endl << "ERROR: Path to the file with parameters is not defined!" << std::endl << std::endl;
					exit(EXIT_FAILURE);
				}

				std::ifstream file(filePath.c_str(), std::ios_base::in);
            	if (!file) {

                	std::cerr << std::endl << "ERROR: Error loading file with LOD parameters!" << std::endl << std::endl;
                	exit(EXIT_FAILURE);
            	}

            	Input_parameters tmp_params;
            	while (!file.eof()) {

            		std::string param_name, param_value;
            		file >> param_name >> param_value;

            		if (param_name == "" || param_value == "") {
            		
            			continue;
            			std::cerr << "ERROR: Empty parameter values!" << std::endl << std::endl;
            			exit(EXIT_FAILURE);
            		}

            		tmp_params[param_name] = param_value;
            	}

            	input_parameters = tmp_params;
            	file.close();
			}

			void save_parameters_to_file(const Input_parameters &input_parameters, const std::vector<std::string> &exceptions) {
				
				if (auto path = get_log_path()) {
					const std::string data_path = static_cast<std::string>(path);

					const std::string full_path = data_path + "logs" + std::string(PS) + "params.lod";
					save_input_parameters(full_path, input_parameters, exceptions);

					std::cout << std::endl << "Parameters are saved in : " << full_path << std::endl << std::endl;
					return;
				}

				std::cerr << std::endl << "ERROR: It is not possible to save parameters, because the LOD_LOG_PATH environment variable is not defined!" << std::endl << std::endl;
				exit(EXIT_FAILURE);
			}

			char* get_log_path() {
				return std::getenv("LOD_LOG_PATH");
			}

			void save_input_parameters(const std::string &filePath, const Input_parameters &input_parameters, const std::vector<std::string> &exceptions) {

				std::ofstream file(filePath.c_str(), std::ios_base::out);

				if (!file) {
					std::cerr << std::endl << "ERROR: Error saving log file with the name " << filePath << std::endl << std::endl;
					exit(EXIT_FAILURE);
				}

				for (Input_parameters::const_iterator param = input_parameters.begin(); param != input_parameters.end(); ++param)
					if (parameter_should_be_saved((*param).first, exceptions))
						file << (*param).first << " " << (*param).second << std::endl;

				file.close();
			}

			bool parameter_should_be_saved(const std::string &parameter_name, const std::vector<std::string> &exceptions) {

				for (size_t i = 0; i < exceptions.size(); ++i)
					if (exceptions[i] == parameter_name) return false;

				return true;
			}
		};

	} // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_ENUM_H