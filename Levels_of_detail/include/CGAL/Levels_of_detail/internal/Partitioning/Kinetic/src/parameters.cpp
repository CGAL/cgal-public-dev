#include "../include/parameters.h"
#include <iostream>
#include <boost/filesystem.hpp>


namespace Skippy {
	Kinetic_Parameters::Kinetic_Parameters()
	{
		reset();

/*#ifdef KINETIC_PARTITION_DEBUG
		print_schedule = true;
#endif*/
	}


	Kinetic_Parameters::~Kinetic_Parameters()
	{
	}


	void Kinetic_Parameters::reset()
	{
		rand_n = 0;
		rand_p = 4;
		rand_d = 0.1;

		path = "";
		location = "";
		basename = "";

		boxes = 10;
		f_margin = FT(1) / FT(1000);

		fast_schedule = true;
		D = FT(1);
		D_2 = FT(1); // D^2
		D_inf_2 = FT(9801) / FT(10000); // (0.99 * D)^2
		D_sup_2 = FT(10201) / FT(10000);  // (1.01 * D)^2
		D_is_set = false;

		output_facets = false;
		output_polyhedrons = false;
		output_graph = false;

		check = false;
		rt_check = false;
		print_schedule = false;
		print_drawings = false;

		stopping_condition = 0;
		K = 0;

		provide_graph_definition = false;
		discretize_orientations = true;

		cl_flush = -1;
		use_landmarks = false;
		process = NONE;

		animator = false;
		t_animator = FT(0);

		p_dist_min = FT(1) / FT(10000);
		reorient = false;
		reorient_2 = false;
		reorient_3 = false;
		
		use_grid = false;
		grid_x = 1;
		grid_y = 1;
		grid_z = 1;

		perfs = false;
	}


	void Kinetic_Parameters::message() const
	{
		std::cout << "Kinetic-Partition-3D.exe [--rand N | --input PATH]" << std::endl;
	}



	void Kinetic_Parameters::set_options(int argc, char *argv[])
	{
		if (argc == 0) return;

		// Tests if the first parameter is the name of an executable
		int r = 0;
		std::string binary = std::string(argv[0]);
		if (binary.size() > 4 && binary.substr(binary.size() - 4, 4) == ".exe") {
			++r;
		}

		while (r < argc) {
			if (!strcmp(argv[r], "--rn") && r + 1 < argc) {
				rand_n = atoi(argv[r + 1]);
				r += 2;
			} else if (!strcmp(argv[r], "--rp") && r + 1 < argc) {
				rand_p = atoi(argv[r + 1]);
				r += 2;
			} else if (!strcmp(argv[r], "--rd") && r + 1 < argc) {
				rand_d = atof(argv[r + 1]);
				r += 2;
			} else if (!strcmp(argv[r], "--input") && r + 1 < argc) {
				path = argv[r + 1];
				location = boost::filesystem::path(path).parent_path().string();
				basename = boost::filesystem::path(path).stem().string();
				r += 2;
			} else if (!strcmp(argv[r], "--planes") && r + 1 < argc) {
				path_planes = argv[r + 1];
				r += 2;
			} else if (!strcmp(argv[r], "--point-cloud") && r + 1 < argc) {
				path_point_cloud = argv[r + 1];
				r += 2;
			} else if (!strcmp(argv[r], "--box") && r + 1 < argc) {
				boxes = atoi(argv[r + 1]);
				r += 2;
			} else if (!strcmp(argv[r], "--margin") && r + 1 < argc) {
				std::stringstream stream;
				stream << std::string(argv[r + 1]);
				stream >> f_margin;
				r += 2;
			} else if (!strcmp(argv[r], "--tau") && r + 1 < argc) {
				std::stringstream stream;
				stream << std::string(argv[r + 1]);
				stream >> D;
				D_2 = D * D;
				FT D_inf = FT(99) * D / FT(100), D_sup = FT(101) * D / FT(100);
				D_inf_2 = D_inf * D_inf;
				D_sup_2 = D_sup * D_sup;
				D_is_set = true;
				r += 2;
			} else if (!strcmp(argv[r], "--polyhedrons")) {
				output_polyhedrons = true;
				r += 1;
			} else if (!strcmp(argv[r], "--facets")) {
				output_facets = true;
				r += 1;
			} else if (!strcmp(argv[r], "--graph")) {
				output_graph = true;
				r += 1;
			} else if (!strcmp(argv[r], "--check")) {
				check = true;
				r += 1;
			} else if (!strcmp(argv[r], "--rt-check")) {
				rt_check = true;
				r += 1;
			} else if (!strcmp(argv[r], "--print-schedule")) {
				print_schedule = true;
				r += 1;
			} else if (!strcmp(argv[r], "--print-drawings")) {
				print_drawings = true;
				r += 1;
			} else if (!strcmp(argv[r], "--K") && r + 1 < argc) {
				stopping_condition = 0;
				K = atoi(argv[r + 1]) - 1;
				r += 2;
			} else if (!strcmp(argv[r], "--density-box") && r + 5 < argc) {
				stopping_condition = 1;
				std::stringstream stream;
				stream << argv[r + 1] << " " << argv[r + 2] << " " << argv[r + 3] << " " << argv[r + 4] << " " << argv[r + 5];
				stream >> density_box_length >> density_box_width >> density_box_height >> density_pts >> density_deviation;
				r += 6;
			} else if (!strcmp(argv[r], "--density-cone") && r + 4 < argc) {
				stopping_condition = 2;
				std::stringstream stream;
				stream << argv[r + 1] << " " << argv[r + 2] << " " << argv[r + 3] << " " << argv[r + 4];
				stream >> density_cone_base >> density_cone_height >> density_pts >> density_deviation;
				r += 5;
			} else if (!strcmp(argv[r], "--graph-definition")) {
				provide_graph_definition = true;
				r += 1;
			} else if (!strcmp(argv[r], "--no-discretization")) {
				discretize_orientations = false;
				r += 1;
			} else if (!strcmp(argv[r], "--no-fast-schedule")) {
				fast_schedule = false;
				r += 1;
			} else if (!strcmp(argv[r], "--landmarks")) {
				use_landmarks = true;
				r += 1;
			} else if (!strcmp(argv[r], "--cl-flush") && r + 1 < argc) {
				cl_flush = atoi(argv[r + 1]);
				r += 2;
			} else if (!strcmp(argv[r], "--preprocess") && r + 1 < argc) {
				if (!strcmp(argv[r + 1], "RECTANGLES")) {
					process = Preprocess::OPTIMAL_RECTANGLES;
				} else {
					goto unknown_option;
				}
				r += 2;
			} else if (!strcmp(argv[r], "--animator") && r + 1 < argc) {
				animator = true;
				std::stringstream stream;
				stream << argv[r + 1];
				stream >> t_animator;
				r += 2;
			} else if (!strcmp(argv[r], "--p-dist-min") && r + 1 < argc) {
				std::stringstream stream;
				stream << argv[r + 1];
				stream >> p_dist_min;
				r += 2;
			} else if (!strcmp(argv[r], "--reorient-2d")) {
				reorient_2 = true;
				reorient = true;
				r += 1;
			} else if (!strcmp(argv[r], "--reorient-3d")) {
				reorient_3 = true;
				reorient = true;
				r += 1;
			} else if (!strcmp(argv[r], "--grid") && r + 3 < argc) {
				use_grid = true;
				grid_x = atoi(argv[r + 1]);
				grid_y = atoi(argv[r + 2]);
				grid_z = atoi(argv[r + 3]);
				r += 4;
			} else if (!strcmp(argv[r], "--perfs")) {
				perfs = true;
				r += 1;
			} else {

			unknown_option:
				std::cout << "Unknown option : " << argv[r] << std::endl;
				reset();
				message();
				break;
			}
		}
	}



	void Kinetic_Parameters::set_option(const std::string & option)
	{
		char* args[1] = { strdup(option.c_str()) };
		set_options(1, args);

		free(args[0]);
	}



	void Kinetic_Parameters::set_option(const std::string & option, const std::string & argument_1)
	{
		if (option == "--basename") {
			basename = argument_1;
		} else {
			char* args[2] = { strdup(option.c_str()), strdup(argument_1.c_str()) };
			set_options(2, args);

			free(args[0]);
			free(args[1]);
		}
	}



	void Kinetic_Parameters::set_option(const std::string & option, const std::string & argument_1, const std::string & argument_2)
	{
		char* args[3] = { strdup(option.c_str()), strdup(argument_1.c_str()), strdup(argument_2.c_str()) };
		set_options(3, args);

		free(args[0]);
		free(args[1]);
		free(args[2]);
	}


	void Kinetic_Parameters::set_option(const std::string & option, const std::vector<std::string> & arguments)
	{
		size_t n = arguments.size();
		char **args = new char*[n + 1];
		args[0] = strdup(option.c_str());
		for (size_t i = 0 ; i < n ; ++i) {
			args[i + 1] = strdup(arguments[i].c_str());
		}

		set_options(int(n + 1), args);
		for (size_t i = 0 ; i < n ; ++i) {
			free(args[i]);
		}
	}
}