#include "../include/propagation_simple.h"
#include "../include/propagation.h"
#include "../include/universe.h"
#include "../include/parameters.h"
#include "../include/vars.h"
#include "../include/stats.h"
#include "../include/support_plane.h"
#include "../include/event_queue.h"
#include "../include/partition_simple.h"
#include "../include/oriented_bounding_box.h"


namespace Skippy
{
	Kinetic_Propagation_Simple::Kinetic_Propagation_Simple()
		: Kinetic_Propagation()
	{
	
	}



	Kinetic_Propagation_Simple::Kinetic_Propagation_Simple(const int N)
		: Kinetic_Propagation(N)
	{
	
	}



	Kinetic_Propagation_Simple::Kinetic_Propagation_Simple(int argc, char *argv[])
		: Kinetic_Propagation(argc, argv)
	{
	
	}



	Kinetic_Propagation_Simple::Kinetic_Propagation_Simple(const std::string & filename, Preprocess process)
		: Kinetic_Propagation(filename, process)
	{
	
	}



	Kinetic_Propagation_Simple::Kinetic_Propagation_Simple(const std::vector<std::vector<CGAL_Inexact_Point_3> > & primitives, Preprocess process)
		: Kinetic_Propagation(primitives, process)
	{
	
	}



	Kinetic_Propagation_Simple::Kinetic_Propagation_Simple(const std::vector<std::vector<CGAL_Point_3> > & primitives, Preprocess process)
		: Kinetic_Propagation(primitives, process)
	{
	
	}



	Kinetic_Propagation_Simple::~Kinetic_Propagation_Simple()
	{
	
	}



	void Kinetic_Propagation_Simple::run()
	{
		// bool perfs = Universe::params->perfs;

		// std::cout << "Basename : " << Universe::params->basename << std::endl;

		if (polygons.empty()) return;

		// clock_t t_0 = clock();

		try {
			// Initializes a kinetic data structure
			init_unique_kinetic_data_structure();

			// Processes all events
			unstack();

			// Builds the partition
			build_partition();

			// Cleans memory
			// delete_unique_kinetic_data_structure();

		} catch (std::exception & except) {
			std::cout << except.what() << std::endl;
		}

		// KP_Stats::running_time = (clock() - t_0) / double(CLOCKS_PER_SEC);
		// KP_Stats::p_initialization_time = 100 * KP_Stats::initialization_time / KP_Stats::running_time;
		// KP_Stats::p_unstack_time = 100 * KP_Stats::unstack_time / KP_Stats::running_time;
		// KP_Stats::p_partition_time = 100 * KP_Stats::partition_time / KP_Stats::running_time;
		// KP_Stats::p_destruction_time = 100 * KP_Stats::destruction_time / KP_Stats::running_time;

		// std::cout << "**" << std::endl;
		// std::cout << "Planes           : " << planes.size() << std::endl;
		// std::cout << "Vertices (t = 0) : " << KP_Stats::initial_vertices << std::endl;
		// std::cout << "Vertices (peak)  : " << KP_Stats::peak_vertices << std::endl;
		// std::cout << "Events           : " << KP_Stats::unstacked_events << std::endl;

		// std::cout << "**" << std::endl
		// 	<< "In average, for each vertex : " << std::endl
		// 	<< "Events computed : " << double(KP_Stats::schedule_events_lines) / KP_Stats::schedule_events_vertices << std::endl
		// 	<< "Events popped   : " << 1 + double(KP_Stats::life_expectancy_lines) / KP_Stats::life_expectancy_terms << std::endl
		// 	<< "Distance ran    : " << double(KP_Stats::life_expectancy_distance) / KP_Stats::life_expectancy_terms << std::endl;

		// print_performances(perfs);
	}



	void Kinetic_Propagation_Simple::init_unique_kinetic_data_structure()
	{
		clock_t t_0 = clock();

		CGAL_Point_3 pt_min, pt_max;
		std::vector<CGAL_Point_3> box_corners;
		std::vector<std::pair<int, int> > box_edges;
		std::map<int, std::map<int, CGAL_Line_3> > lines;

		try {
			const std::string & filename_point_cloud = Universe::params->path_point_cloud;
			if (!filename_point_cloud.empty() && Universe::params->stopping_condition != 0) {
				init_point_cloud_octree(filename_point_cloud);
			}

			// Initialization, part 1 : geometric data structures

			pre_init_main_bounding_box();

			init_plane_equations();
			init_colors(0);

			/*if (Universe::params->reorient_3) {
				Universe::bounding_box->get_bounding_box(Universe::params->path_point_cloud, planes);
			}*/

			init_main_bounding_box(pt_min, pt_max, box_corners, box_edges);
			init_intersection_lines(lines);

			init_supporting_planes(pt_min, pt_max, box_corners, box_edges, lines);
			init_queue_of_events();
			build_polygons();

			// Initialization, part 2 : optimized structure of queue
			init_schedule();

		} catch (std::exception & except) {
			throw except;
		}

		// std::cout << "** Initialized structure" << std::endl;

		KP_Stats::initial_vertices = 0;
		for (size_t i = 0 ; i < Universe::map_of_planes.size() ; ++i) {
			Support_Plane* SP = Universe::map_of_planes[i];
			KP_Stats::initial_vertices += int(SP->vertices_r.size());
		}

		KP_Stats::initialization_time = (clock() - t_0) / double(CLOCKS_PER_SEC);
	}



	void Kinetic_Propagation_Simple::delete_unique_kinetic_data_structure()
	{
		clock_t t_0 = clock();

		// Deletes the support planes and their objects
		for (size_t i = 0; i < Universe::map_of_planes.size(); i++) {
			delete Universe::map_of_planes[i];
		}
		Universe::map_of_planes.clear();

		// Deletes the queue
		delete Universe::event_queue;

		// Reset counters
		Counters::id_planes = -1;
		Counters::id_objects = std::vector<int>();

		Counters::id_partition_vertex = -1;
		Counters::id_partition_edge = -1;
		Counters::id_partition_facet = -1;
		Counters::id_partition_polyhedron = -1;

		Counters::par_v_local_ids = std::vector<int>();
		Counters::par_e_local_ids = std::vector<int>();

		KP_Stats::destruction_time = (clock() - t_0) / double(CLOCKS_PER_SEC);

		delete Universe::params;
	}



	void Kinetic_Propagation_Simple::build_polygons() const
	{
		// Over a first phase, polygons are simply bound to their support plane
		// Segments corresponding to this polygon are created on the other planes
		for (size_t i = 0; i < polygons.size(); i++) {

			const std::vector<CGAL_Point_3> & P = polygons[i];

			Support_Plane* S = Universe::map_of_planes[polygons_to_planes[i]];
			S->init_polygon(P);
		}
	}



	void Kinetic_Propagation_Simple::build_partition()
	{
		clock_t t_0 = clock();

		partition = new Partition_Simple(planes);

		Partition_Simple* this_partition = dynamic_cast<Partition_Simple*>(partition);
		this_partition->build();

		KP_Stats::final_nb_vertices = Counters::id_partition_vertex;
		KP_Stats::final_nb_edges = Counters::id_partition_edge;
		KP_Stats::final_nb_facets = Counters::id_partition_facet;
		KP_Stats::final_nb_polyhedrons = Counters::id_partition_polyhedron;

		// std::cout << "** Built a partition with " << partition->polyhedrons_size() << " polyhedrons" << std::endl;

		KP_Stats::partition_time = (clock() - t_0) / double(CLOCKS_PER_SEC);

		export_partition();
	}
}