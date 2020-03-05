#include "../include/stats.h"

namespace Skippy 
{
	int KP_Stats::insert_temporally_calls = 0;
	int KP_Stats::insert_temporally_loop_case = 0;
	int KP_Stats::insert_temporally_loop_comparisons = 0;

	clock_t KP_Stats::schedule_events_search_time = 0;
	clock_t KP_Stats::schedule_events_computation_time = 0;

	int KP_Stats::schedule_events_vertices = 0;
	int KP_Stats::schedule_events_lines = 0;

	int KP_Stats::life_expectancy_lines = 0;
	double KP_Stats::life_expectancy_distance = 0;
	int KP_Stats::life_expectancy_terms = 0;


	double KP_Stats::visible_planes_mean = 0;
	double KP_Stats::visible_planes_sd = 0;
	
	double KP_Stats::running_time = 0;
	double KP_Stats::initialization_time = 0;
	double KP_Stats::unstack_time = 0;
	double KP_Stats::partition_time = 0;
	double KP_Stats::destruction_time = 0;

	double KP_Stats::p_initialization_time = 0;
	double KP_Stats::p_unstack_time = 0;
	double KP_Stats::p_partition_time = 0;
	double KP_Stats::p_destruction_time = 0;

	int KP_Stats::initial_vertices = 0;
	int KP_Stats::peak_vertices = 0;
	int KP_Stats::unstacked_events = 0;

	int KP_Stats::final_nb_vertices = 0;
	int KP_Stats::final_nb_edges = 0;
	int KP_Stats::final_nb_facets = 0;
	int KP_Stats::final_nb_polyhedrons = 0;
}