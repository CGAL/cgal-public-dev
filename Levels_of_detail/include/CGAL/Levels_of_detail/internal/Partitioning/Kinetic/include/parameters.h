#pragma once
#include "defs.h"
#include "defs_cgal.h"
#include <string>


namespace Skippy {

	class Kinetic_Parameters
	{
	public:
		Kinetic_Parameters();

		~Kinetic_Parameters();

		void message() const;

		void set_options(int argc, char *argv[]);

		void set_option(const std::string & option);

		void set_option(const std::string & option, const std::string & argument_1);

		void set_option(const std::string & option, const std::string & argument_1, const std::string & argument_2);

		void set_option(const std::string & option, const std::vector<std::string> & arguments);

	protected:
		void reset();

	public:
		int rand_n;
		int rand_p;
		double rand_d;

		std::string path;
		std::string location;
		std::string basename;
		std::string path_planes;
		std::string path_point_cloud;

		uint boxes;
		FT f_margin;

		bool fast_schedule;
		FT D;
		FT D_2;
		FT D_inf_2;
		FT D_sup_2;
		bool D_is_set;

		bool output_facets;
		bool output_polyhedrons;
		bool output_graph;

		bool check;
		bool rt_check;
		bool print_schedule;
		bool print_drawings;

		int stopping_condition;
		int K;

		FT density_box_length;
		FT density_box_width;
		FT density_box_height;

		FT density_cone_base;
		FT density_cone_height;

		int density_pts;
		double density_deviation;

		bool provide_graph_definition;
		bool discretize_orientations;

		int cl_flush;
		bool use_landmarks;

		Preprocess process;

		bool animator;
		FT t_animator;

		FT p_dist_min;
		bool reorient;
		bool reorient_2;
		bool reorient_3;

		bool use_grid;
		int grid_x;
		int grid_y;
		int grid_z;

		bool perfs;
	};

}