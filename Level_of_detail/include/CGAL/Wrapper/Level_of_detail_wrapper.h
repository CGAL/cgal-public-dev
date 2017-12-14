#ifndef CGAL_LEVEL_OF_DETAIL_WRAPPER_H
#define CGAL_LEVEL_OF_DETAIL_WRAPPER_H

// STL includes.
#include <string>
#include <iostream>
#include <cstdlib>
#include <map>
#include <cassert>

// CGAL new includes.
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Level_of_detail_traits.h>
#include <CGAL/Tools/Level_of_detail_tools.h>
#include <CGAL/Base/Level_of_detail_base.h>

namespace CGAL {

	namespace LOD {

		template<class LodTraits>
		class Level_of_detail_wrapper {

		public:
			using Kernel = typename LodTraits::Kernel;
			using FT 	 = typename Kernel::FT;

			using Lod_base         = CGAL::LOD::Level_of_detail_base<LodTraits>;
			using Input_parameters = std::map<std::string, std::string>;
			using Lod_parameters   = CGAL::LOD::Level_of_detail_parameters<FT>;
			using Params           = char**;

			enum class Wrapper_mode { TEST, NATIVE };

			Level_of_detail_wrapper(const int num_params, const Params params, const Wrapper_mode wrapper_mode = Wrapper_mode::NATIVE) : m_data_path("default") {

				switch (wrapper_mode) {

					case Wrapper_mode::NATIVE:
						use_native_mode(num_params, params);
						break;

					// DO NOT FORGET TO CHANGE THE BASE CLASS TO THE BASE_VER0 CLASS!
					case Wrapper_mode::TEST:
						use_test_mode(num_params, params);
						break;

					default:
						assert(!"Wrong wrapper mode!");
						break;
				}
			}

			void run_lod_pipeline() {
				m_lod_base.create_lods();
			}

		private:
			Lod_base m_lod_base;
			std::string m_data_path;

			void use_native_mode(const int num_params, const Params params) {

				Lod_parameters lod_parameters(num_params, params);

				m_lod_base.set_optimal_configuration();
				m_lod_base.set_user_defined_parameters(lod_parameters);
			}

			void use_test_mode(const int num_params, const Params params) {
				
				need_help(num_params, params);

				Input_parameters input_parameters;
				set_input_parameters(num_params, params, input_parameters);

				set_data_path(input_parameters);
				std::cout << "data path : " << m_data_path << std::endl;

				set_user_defined_parameters(input_parameters);
			}

			void set_input_parameters(const int num_params, const Params params, Input_parameters &input_parameters) {

				if (is_valid_input(num_params)) {
		
					set_default_input_parameters(input_parameters);
					std::cout << std::endl << "User defined parameter values (all other values are default): " << std::endl;

					for (int i = 1; i < num_params; i += 2) input_parameters[params[i]] = params[i + 1];
					return;
				}

				std::cerr << std::endl << "ERROR: Missing parameter values. Check your input!" << std::endl << std::endl;
				exit(EXIT_FAILURE);
			}

			void set_data_path(const Input_parameters &input_parameters) {

				bool data_path_found = try_input_path(input_parameters);
				if (data_path_found) return;

				data_path_found = try_default_settings();
				if (data_path_found) return;

				std::cerr << std::endl << "ERROR: You did not provide any input data!" << std::endl;
				std::cerr << std::endl << "Set LOD_DATA_PATH environment variable to your data folder and use the -data_type option or simply use the -data option with the full path." << std::endl << std::endl;
				exit(EXIT_FAILURE);
			}

			void need_help(const int num_params, const Params params) {
				if (!is_asked_for_help(num_params, params)) return;
	
				print_help();
				exit(EXIT_SUCCESS);
			}

			bool is_asked_for_help(const int num_params, const Params params) {
				return num_params == 2 && std::strcmp(params[1], "-help") == 0;
			}

			bool is_valid_input(const int num_params) {
				return num_params % 2 != 0;
			}

			bool try_input_path(const Input_parameters &input_parameters) {
				if (input_parameters.at("-data") != "default") {

					m_data_path = static_cast<std::string>(input_parameters.at("-data").c_str());
					return true;
				}
				return false;
			}

			bool try_default_settings() {
				
				if (const char* data_path = std::getenv("LOD_DATA_PATH")) {
					m_data_path = static_cast<std::string>(data_path);

					m_lod_base.set_prefix_path(m_data_path);
					m_lod_base.set_default_parameters();

					return true;
				}
				return false;
			}

			void set_default_input_parameters(Input_parameters &input_parameters) {
				
				input_parameters["-data"]      = "default";
				input_parameters["-data_type"] = "default";
			    input_parameters["-silent"]    = "default";

			    input_parameters["-clutter"] = "default";
			    input_parameters["-cell"]    = "default";

			    input_parameters["-rg_eps"] 		= "default";
			    input_parameters["-rg_clust_eps"]   = "default";
			    input_parameters["-rg_norm_thresh"] = "default";
			    input_parameters["-rg_min_points"]  = "default";

				input_parameters["-str_eps"] = "default";
				input_parameters["-str_adj"] = "default";
				input_parameters["-str_all"] = "default";

				input_parameters["-gc_beta"]  = "default";
				input_parameters["-gc_gamma"] = "default";
			}

			void set_user_defined_parameters(const Input_parameters &input_parameters) {

				// Input data. This parameter must be assigned first! It sets default values for all other parameters!
				if (input_parameters.at("-data_type") != "default") {

					const size_t data_type = static_cast<size_t>( std::stoi(input_parameters.at("-data_type").c_str()) );
					std::cout << "data type (5 - PARIS_HALF_TILE, 6 - PARIS_FULL_TILE, 10 - PARIS_TILE_SPARCE, 11 - PARIS_TILE_DENSE, 12 - PARIS_9_TILES) : " << data_type << std::endl;
					m_lod_base.set_data_type(data_type);
				}


				// Cout.
				if (input_parameters.at("-silent") != "default") {

					const bool silent = static_cast<bool>( std::stoi(input_parameters.at("-silent").c_str()) );
					std::cout << "use silent mode (0 - false, 1 - true) : " << silent << std::endl;
					m_lod_base.make_silent(silent);
				}


				// Clutter.
				if (input_parameters.at("-clutter") != "default") {

					const bool clutter = static_cast<bool>( std::stoi(input_parameters.at("-clutter").c_str()) );
					std::cout << "use clutter (0 - false, 1 - true) : " << clutter << std::endl;
					m_lod_base.add_clutter(clutter);
				}

				if (input_parameters.at("-cell") != "default") {

					const FT cell = static_cast<FT>( std::stod(input_parameters.at("-cell").c_str()) );
					std::cout << "clutter cell side length : " << cell << std::endl;
					m_lod_base.set_clutter_cell_side_length(cell);
				}


				// Region growing.
				if (input_parameters.at("-rg_eps") != "default") {

					const FT rg_eps = static_cast<FT>( std::stod(input_parameters.at("-rg_eps").c_str()) );
					std::cout << "region growing eps : " << rg_eps << std::endl;
					m_lod_base.set_region_growing_epsilon(rg_eps);
				}

				if (input_parameters.at("-rg_clust_eps") != "default") {

					const FT rg_clust_eps = static_cast<FT>( std::stod(input_parameters.at("-rg_clust_eps").c_str()) );
					std::cout << "region growing cluster eps : " << rg_clust_eps << std::endl;
					m_lod_base.set_region_growing_cluster_epsilon(rg_clust_eps);
				}

				if (input_parameters.at("-rg_norm_thresh") != "default") {

					const FT rg_norm_thresh = static_cast<FT>( std::stod(input_parameters.at("-rg_norm_thresh").c_str()) );
					std::cout << "region growing normal threshold : " << rg_norm_thresh << std::endl;
					m_lod_base.set_region_growing_normal_threshold(rg_norm_thresh);
				}

				if (input_parameters.at("-rg_min_points") != "default") {

					const size_t rg_min_points = static_cast<size_t>( std::stoi(input_parameters.at("-rg_min_points").c_str()) );
					std::cout << "region growing min points : " << rg_min_points << std::endl;
					m_lod_base.set_region_growing_min_points(rg_min_points);
				}


				// Structuring.
				if (input_parameters.at("-str_eps") != "default") {

					const FT str_eps = static_cast<FT>( std::stod(input_parameters.at("-str_eps").c_str()) );
					std::cout << "structuring eps : " << str_eps << std::endl;
					m_lod_base.set_structuring_epsilon(str_eps);
				}

				if (input_parameters.at("-str_adj") != "default") {

					const FT str_adj = static_cast<FT>( std::stod(input_parameters.at("-str_adj").c_str()) );
					std::cout << "structuring adjacency threshold : " << str_adj << std::endl;
					m_lod_base.set_structuring_adjacency_value(str_adj);
				}

				if (input_parameters.at("-str_all") != "default") {

					const bool str_all = static_cast<bool>( std::stoi(input_parameters.at("-str_all").c_str()) );
					std::cout << "structuring get all points (0 - false, 1 - true) : " << str_all << std::endl;
					m_lod_base.get_all_structuring_points(str_all);
				}


				// Graph cut.
				if (input_parameters.at("-gc_beta") != "default") {

					const FT gc_beta = static_cast<FT>( std::stod(input_parameters.at("-gc_beta").c_str()) );
					std::cout << "graph cut main parameter (beta) : " << gc_beta << std::endl;
					m_lod_base.set_graph_cut_beta(gc_beta);
				}

				if (input_parameters.at("-gc_gamma") != "default") {

					const FT gc_gamma = static_cast<FT>( std::stod(input_parameters.at("-gc_gamma").c_str()) );
					std::cout << "graph cut penalty (gamma) : " << gc_gamma << std::endl;
					m_lod_base.set_graph_cut_gamma(gc_gamma);
				}
			}

			void print_help() {

				std::cout << std::endl << "HELP:" << std::endl;

				
				std::cout << std::endl << "EXAMPLE:" << std::endl;
				std::cout << "your terminal name $ ./lod -data_type 5" << std::endl;
				std::cout << "your terminal name $ ./lod -data path_to_data/data_name.ply -str_eps 5.0 -silent 1" << std::endl << std::endl;


				std::cout << std::endl << "BASIC:" << std::endl;

				std::cout << 
				"param name: -data" 					   << std::endl <<
				"param values: path_to_data/data_name.ply" << std::endl <<
				"description: set data path" 			   << std::endl << std::endl;

				std::cout << 
				"param name: -data_type" 																								    << std::endl <<
				"param values: 5 - PARIS_HALF_TILE, 6 - PARIS_FULL_TILE, 10 - PARIS_TILE_SPARCE, 11 - PARIS_TILE_DENSE, 12 - PARIS_9_TILES" << std::endl <<
				"description: set data type, LOD_DATA_PATH environment variable must be set" 												<< std::endl << std::endl;

				std::cout << 
				"param name: -silent" 						   << std::endl <<
				"param values: 0, 1" 						   << std::endl <<
				"description: save temporary steps 1 or not 0" << std::endl << std::endl;


				std::cout << std::endl << "CLUTTER:" << std::endl;

				std::cout << 
				"param name: -clutter" 						 << std::endl <<
				"param values: 0, 1" 						 << std::endl <<
				"description: add clutter points 1 or not 0" << std::endl << std::endl;

				std::cout << 
				"param name: -cell" 						 		  << std::endl <<
				"param values: > 0.0" 						 		  << std::endl <<
				"description: cell side length used in grid simplify" << std::endl << std::endl;


				std::cout << std::endl << "REGION GROWING:" << std::endl;

				std::cout << 
				"param name: -rg_eps" 						  			   << std::endl <<
				"param values: > 0.0" 						  			   << std::endl <<
				"description: distance from the point to the optimal line" << std::endl << std::endl;

				std::cout << 
				"param name: -rg_clust_eps" 				      << std::endl <<
				"param values: > 0.0" 						      << std::endl <<
				"description: distance among neighbouring points" << std::endl << std::endl;

				std::cout << 
				"param name: -rg_norm_thresh" 				  								  << std::endl <<
				"param values: > 0.0 && < 1.0" 				  								  << std::endl <<
				"description: cosine between the point normal and normal of the optimal line" << std::endl << std::endl;

				std::cout << 
				"param name: -rg_min_points" 						  					<< std::endl <<
				"param values: > 0" 						  							<< std::endl <<
				"description: minimum number of points that can contribute to the line" << std::endl << std::endl;


				std::cout << std::endl << "STRUCTURING:" << std::endl;

				std::cout << 
				"param name: -str_eps" 						  						  << std::endl <<
				"param values: > 0.0" 						  						  << std::endl <<
				"description: distance between adjacent points in the resampled line" << std::endl << std::endl;

				std::cout << 
				"param name: -str_adj" 						  						<< std::endl <<
				"param values: > 0.0" 						  						<< std::endl <<
				"description: max distance between end points of adjacent segments" << std::endl << std::endl;

				std::cout << 
				"param name: -str_all" 						  							  << std::endl <<
				"param values: 0, 1" 						  							  << std::endl <<
				"description: use all resampled points or only end points of the segment" << std::endl << std::endl;


				std::cout << std::endl << "GRAPH CUT:" << std::endl;

				std::cout << 
				"param name: -gc_beta" 						   << std::endl <<
				"param values: > 0.0" 						   << std::endl <<
				"description: main parameter of the graph cut" << std::endl << std::endl;

				std::cout << 
				"param name: -gc_gamma" 		 << std::endl <<
				"param values: > 0.0" 			 << std::endl <<
				"description: graph cut penalty" << std::endl;

				std::cout << std::endl;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_WRAPPER_H