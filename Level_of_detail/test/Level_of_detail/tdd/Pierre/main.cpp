#if defined(WIN32) || defined(_WIN32) 
#define PS "\\"
#define PN "\r\n"
#else 
#define PS "/" 
#define PN "\n"
#endif

// STL includes.
#include <string>
#include <iostream>
#include <cstdlib>

// CGAL includes.
// #include <CGAL/Simple_cartesian.h>
// #include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>

// CGAL new includes.
#include <CGAL/Level_of_detail_traits.h>
#include <CGAL/Base/Level_of_detail_base.h>

// using Kernel     = CGAL::Simple_cartesian<double>;
// using Kernel     = CGAL::Exact_predicates_exact_constructions_kernel;

using Kernel     = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT 		 = Kernel::FT;
using Point      = Kernel::Point_3;
using Container  = CGAL::Point_set_3<Point>;
using LodTraits  = CGAL::LOD::Level_of_detail_traits<Kernel, Container>;
using LodBase    = CGAL::LOD::Level_of_detail_base<LodTraits>;
using Parameters = std::map<std::string, std::string>;


void set_user_defined_parameters(const Parameters &parameters, LodBase &lodBase) {

	// Input data. This parameter must be assigned first! It sets default values for all other parameters!
	if (parameters.at("-data") != "default") {

		const size_t data_type = static_cast<size_t>( std::stoi(parameters.at("-data").c_str()) );
		std::cout << "data type (5 - PARIS_HALF_TILE, 6 - PARIS_FULL_TILE, 10 - PARIS_TILE_SPARCE, 11 - PARIS_TILE_DENSE, 12 - PARIS_9_TILES) : " << data_type << std::endl;
		lodBase.set_data_type(data_type);
	}


	// Cout.
	if (parameters.at("-silent") != "default") {

		const bool silent = static_cast<bool>( std::stoi(parameters.at("-silent").c_str()) );
		std::cout << "use silent mode (0 - false, 1 - true) : " << silent << std::endl;
		lodBase.make_silent(silent);
	}


	// Clutter.
	if (parameters.at("-clutter") != "default") {

		const bool clutter = static_cast<bool>( std::stoi(parameters.at("-clutter").c_str()) );
		std::cout << "use clutter (0 - false, 1 - true) : " << clutter << std::endl;
		lodBase.add_clutter(clutter);
	}

	if (parameters.at("-cell") != "default") {

		const FT cell = static_cast<FT>( std::stod(parameters.at("-cell").c_str()) );
		std::cout << "clutter cell side length : " << cell << std::endl;
		lodBase.set_clutter_cell_side_length(cell);
	}


	// Region growing.
	if (parameters.at("-rg_eps") != "default") {

		const FT rg_eps = static_cast<FT>( std::stod(parameters.at("-rg_eps").c_str()) );
		std::cout << "region growing eps : " << rg_eps << std::endl;
		lodBase.set_region_growing_epsilon(rg_eps);
	}

	if (parameters.at("-rg_clust_eps") != "default") {

		const FT rg_clust_eps = static_cast<FT>( std::stod(parameters.at("-rg_clust_eps").c_str()) );
		std::cout << "region growing cluster eps : " << rg_clust_eps << std::endl;
		lodBase.set_region_growing_cluster_epsilon(rg_clust_eps);
	}

	if (parameters.at("-rg_norm_thresh") != "default") {

		const FT rg_norm_thresh = static_cast<FT>( std::stod(parameters.at("-rg_norm_thresh").c_str()) );
		std::cout << "region growing normal threshold : " << rg_norm_thresh << std::endl;
		lodBase.set_region_growing_normal_threshold(rg_norm_thresh);
	}

	if (parameters.at("-rg_min_points") != "default") {

		const size_t rg_min_points = static_cast<size_t>( std::stoi(parameters.at("-rg_min_points").c_str()) );
		std::cout << "region growing min points : " << rg_min_points << std::endl;
		lodBase.set_region_growing_min_points(rg_min_points);
	}


	// Structuring.
	if (parameters.at("-str_eps") != "default") {

		const FT str_eps = static_cast<FT>( std::stod(parameters.at("-str_eps").c_str()) );
		std::cout << "structuring eps : " << str_eps << std::endl;
		lodBase.set_structuring_epsilon(str_eps);
	}

	if (parameters.at("-str_adj") != "default") {

		const FT str_adj = static_cast<FT>( std::stod(parameters.at("-str_adj").c_str()) );
		std::cout << "structuring adjacency threshold : " << str_adj << std::endl;
		lodBase.set_structuring_adjacency_value(str_adj);
	}

	if (parameters.at("-str_all") != "default") {

		const bool str_all = static_cast<bool>( std::stoi(parameters.at("-str_all").c_str()) );
		std::cout << "structuring get all points (0 - false, 1 - true) : " << str_all << std::endl;
		lodBase.get_all_structuring_points(str_all);
	}


	// Graph cut.
	if (parameters.at("-gc_beta") != "default") {

		const FT gc_beta = static_cast<FT>( std::stod(parameters.at("-gc_beta").c_str()) );
		std::cout << "graph cut main parameter (beta) : " << gc_beta << std::endl;
		lodBase.set_graph_cut_beta(gc_beta);
	}

	if (parameters.at("-gc_gamma") != "default") {

		const FT gc_gamma = static_cast<FT>( std::stod(parameters.at("-gc_gamma").c_str()) );
		std::cout << "graph cut penalty (gamma) : " << gc_gamma << std::endl;
		lodBase.set_graph_cut_gamma(gc_gamma);
	}
}

void print_help() {

	std::cout << std::endl << "HELP:" << std::endl;

	
	std::cout << std::endl << "EXAMPLE:" << std::endl;
	std::cout << "your terminal name $ ./lod -data 5 -str_eps 2.0 -silent 1" << std::endl << std::endl;


	std::cout << std::endl << "BASIC:" << std::endl;

	std::cout << 
	"param name: -data" 																										<< std::endl <<
	"param values: 5 - PARIS_HALF_TILE, 6 - PARIS_FULL_TILE, 10 - PARIS_TILE_SPARCE, 11 - PARIS_TILE_DENSE, 12 - PARIS_9_TILES" << std::endl <<
	"description: set data type" 																								<< std::endl << std::endl;

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

int main(int argc, char** argv) {
   
   	LodBase lodBase;
   	const std::string data_path = std::string(std::getenv("LOD_DATA_PATH"));

   	lodBase.set_prefix_path(data_path);
   	lodBase.set_default_parameters();


    // Standard input.
    if (argc == 1) {
    
    	lodBase.create_lods();
    	return 0;
    }


    // Help.
    if (argc == 2 && strcmp(argv[1], "-help") == 0) {
    
    	print_help();
    	return 0;
    }


    // Input with parameters.
    Parameters parameters;
    parameters["-data"]   = "default";
    parameters["-silent"] = "default";

    parameters["-clutter"] = "default";
    parameters["-cell"]    = "default";

    parameters["-rg_eps"] 		  = "default";
    parameters["-rg_clust_eps"]   = "default";
    parameters["-rg_norm_thresh"] = "default";
    parameters["-rg_min_points"]  = "default";

	parameters["-str_eps"] = "default";
	parameters["-str_adj"] = "default";
	parameters["-str_all"] = "default";

	parameters["-gc_beta"]  = "default";
	parameters["-gc_gamma"] = "default";

	if (argc % 2 != 0) {
		std::cout << "" + std::string(PN) + "User defined parameter values (all other values are default): " << std::endl;

		for (int i = 1; i < argc; i += 2)
			parameters[argv[i]] = argv[i + 1];

	} else {
		
		std::cout << "" + std::string(PN) + "Missing parameter values. Check your input!" << std::endl << std::endl;
		return 1;
	}


	// Set user defined parameters.
	set_user_defined_parameters(parameters, lodBase);


	// Run lod.
	std::cout << std::endl;
	lodBase.create_lods();
}