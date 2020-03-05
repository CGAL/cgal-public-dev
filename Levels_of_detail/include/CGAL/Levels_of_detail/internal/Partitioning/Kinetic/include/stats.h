#pragma once
#include <ctime>

namespace Skippy {
	namespace KP_Stats
	{
		extern int insert_temporally_calls;
		extern int insert_temporally_loop_case;
		extern int insert_temporally_loop_comparisons;

		extern clock_t schedule_events_search_time;
		extern clock_t schedule_events_computation_time;

		extern int schedule_events_vertices;
		extern int schedule_events_lines;

		extern int life_expectancy_lines;
		extern double life_expectancy_distance;
		extern int life_expectancy_terms;

		extern double visible_planes_mean;
		extern double visible_planes_sd;

		extern double running_time;
		extern double initialization_time;
		extern double unstack_time;
		extern double partition_time;
		extern double destruction_time;

		extern double p_initialization_time;
		extern double p_unstack_time;
		extern double p_partition_time;
		extern double p_destruction_time;

		extern int initial_vertices;
		extern int peak_vertices;
		extern int unstacked_events;

		extern int final_nb_vertices;
		extern int final_nb_edges;
		extern int final_nb_facets;
		extern int final_nb_polyhedrons;
	}
}