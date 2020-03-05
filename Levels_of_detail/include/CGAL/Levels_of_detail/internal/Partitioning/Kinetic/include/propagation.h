#pragma once
#include "defs.h"
#include "parameters.h"
#include "support_plane_objects.h"
#include "partition.h"


namespace Skippy {
	
	class Kinetic_Propagation
	{
	public:
		KINETIC_PARTITION_API Kinetic_Propagation();

		/**
		 * Constructs a Kinetic_Propagation object by generating N random polygons
		 *
		 * @param N Number of polygons to generate.
		 */
		KINETIC_PARTITION_API Kinetic_Propagation(const int N);

		/**
		 * Constructs a Kinetic_Propagation object by parsing the command line arguments.
		 *
		 * @param argc Number of arguments
		 * @param argv Command-line arguments
		 */
		KINETIC_PARTITION_API Kinetic_Propagation(int argc, char *argv[]);

		/**
		 * Constructs a Kinetic_Propagation object by reading the contents of an input file
		 *
		 * @param filename Address of the filename to read.
		 * @param process Specifies a task to perform on the data
		 */
		KINETIC_PARTITION_API Kinetic_Propagation(const std::string & filename, Preprocess process = NONE);

		/**
		 * Constructs a Kinetic_Propagation object by copying a sequence of polygons.
		 * An inexact float type is used to define such polygons.
		 *
		 * @param primitives The list of polygons to process.
		 * @param process Specifies a task to perform on the data
		 */
		KINETIC_PARTITION_API Kinetic_Propagation(const std::vector<std::vector<CGAL_Inexact_Point_3> > & primitives, Preprocess process = NONE);

		/**
		 * Constructs a Kinetic_Propagation object by copying a sequence of polygons.
		 * An exact float type is used to define such polygons.
		 *
		 * @param primitives The list of polygons to process.
		 * @param process Specifies a task to perform on the data
		 */
		KINETIC_PARTITION_API Kinetic_Propagation(const std::vector<std::vector<CGAL_Point_3> > & primitives, Preprocess process = NONE);

		/**
		 * Destroys the object
		 */
		KINETIC_PARTITION_API virtual ~Kinetic_Propagation();



		/**
		 * Sets options of the kinetic framework.
		 *
		 * @param argc Number of command-line arguments.
		 * @param argv Table of arguments.
		 */
		KINETIC_PARTITION_API void set_options(int argc, char *argv[]) const;

		/**
		 * Sets an option of the kinetic framework.
		 *
		 * @param option Selected option.
		 */
		KINETIC_PARTITION_API void set_option(const std::string & option) const;

		/**
		 * Sets an option of the kinetic framework.
		 *
		 * @param option Selected option.
		 * @param argument_1 Parameter for the selected option.
		 */
		KINETIC_PARTITION_API void set_option(const std::string & option, const std::string & argument_1) const;

		/**
		 * Sets an option of the kinetic framework.
		 *
		 * @param option Selected option.
		 * @param argument_1 Parameter 1 for the selected option.
		 * @param argument_2 Parameter 2 for the selected option.
		 */
		KINETIC_PARTITION_API void set_option(const std::string & option, const std::string & argument_1, const std::string & argument_2) const;

		/**
		 * Sets an option of the kinetic framework.
		 *
		 * @param option Selected option.
		 * @param arguments Vector of parameters for the selected option.
		 */
		KINETIC_PARTITION_API void set_option(const std::string & option, const std::vector<std::string> & arguments) const;

		/**
		 * Tests the existence of (a priori) valid data
		 * on which the partition process can be performed.
		 */
		KINETIC_PARTITION_API bool data() const;

		/**
		 * Processes the data.
		 * This function initializes a kinetic data structure, processes the associated queue of events,
		 * and finally creates a partition that can later be used for future tasks.
		 */
		KINETIC_PARTITION_API virtual void run() = 0;

		KINETIC_PARTITION_API static bool simple_mode(int argc, char *argv[]);

		KINETIC_PARTITION_API static bool perfs_mode(int argc, char *argv[]);

	protected:

		// Generates, reads or processes data

		void get_polygons_from_input_file(const std::string & filename);

		void get_polygons_from_txt_file(const std::string & filename);

		void get_polygons_from_ply_file(const std::string & filename);

		void generate_polygons_randomly(const int N, const int P, const double D);

		void generate_polygons_randomly(const int N);

		void init_point_cloud_octree(const std::string & filename);

		// Preprocesses data

		void apply(Preprocess process);

		void make_convex_envelops();

		void make_fit_rectangles();

		// Initialization of the algorithm

		void pre_init_main_bounding_box();

		void init_main_bounding_box(CGAL_Point_3 & pt_min,
			CGAL_Point_3 & pt_max,
			std::vector<CGAL_Point_3> & box_corners,
			std::vector<std::pair<int, int> > & box_edges);

		void init_plane_equations();

		void init_plane_equations_from_file();

		void init_plane_equations_without_discretization();

		void init_plane_equations_with_discretization(const int subdivs);

		void init_colors(const int slices);

		void init_supporting_planes(const CGAL_Point_3 & pt_min,
			const CGAL_Point_3 & pt_max,
			const std::vector<CGAL_Point_3> & box_corners,
			const std::vector<std::pair<int, int> > & box_edges,
			const std::map<int, std::map<int, CGAL_Line_3> > & lines) const;

		void init_queue_of_events() const;

		void discretize_normal(const int subdivs, const CGAL_Vector_3 & N, int & k_longitude, int & k_latitude) const;

		void discretize_normal(const int subdivs, const CGAL_Vector_3 & N, CGAL_Vector_3 & N_disc) const;

		void init_intersection_lines(std::map<int, std::map<int, CGAL_Line_3> > & lines) const;

		void init_schedule() const;

		FT get_approx_diagonal() const;

		void set_queue_parameters() const;


		// Processing the events

		void unstack();

		void export_partition();

		

		void write_polygons(int seq, const double t) const;

		void write_kinetic_frame (int seq, const FT & t) const;

		void print_performances(const bool write) const;

	public:
		std::vector<std::vector<CGAL_Point_3> > polygons;
		
		std::vector<int> polygons_to_planes;
		std::vector<std::list<int> > planes_to_polygons;
		
		std::vector<CGAL_Plane> planes;
		std::vector<CGAL_Color> colors;

		Partition* partition;
		std::string definition;
	};
}