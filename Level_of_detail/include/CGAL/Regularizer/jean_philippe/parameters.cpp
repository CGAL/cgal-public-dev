#include "parameters.h"

Parameters::Parameters()
{
	reset();
}


Parameters::~Parameters()
{

}


void Parameters::reset()
{
	verbose_level = 0;

	path_input_image = "";
	path_input_directory = "";
	path_ground_truth_directory = "";

	path_output_directory = ".";
	suffix = "";
	zoom = 1.0;
	lsd_test_mode = false;

	// LSD default parameters
	reset_lsd();

	// Regularization default parameters (1)
	reset_rega();

	// Regularization parameters (2)
	reset_regp();

	// Propagation parameters
	reset_prop();
}


void Parameters::reset_lsd()
{
	lsd_scale = 0.8;
	lsd_sigma_scale = 0.6;
	lsd_quant = 2.0;
	lsd_angle = 22.5;
	lsd_log_eps = 0.0;
	lsd_density = 0.7;
	lsd_create_additional = false;
	lsd_additional_size = 4.0;
}


void Parameters::reset_rega()
{
	// enable or not the regularization
	rega_regp_enabled = true;
	
	// mean shift based or quadratic regularization => set to 1 for quadratic
	rega_method = 1;

	// if two segments make an angle smaller than 0.25 degrees they are considered as having the same orientation
	rega_epsilon = 0.25;

	rega_ms_sigma = 4;
	rega_ms_epsilon = 0.125;
	rega_ms_distance = 200;
	rega_ms_min_terms = 6;
	rega_ms_smooth_dist = 3;

	// balancing term
	rega_quad_lambda = 0.8;

	// case when pairwise potentials are considered using a naive euclidian graph with an underlying euclidian distance
	rega_quad_distance = 40;

	// 0 for euclidian i think, 1 for delaunay
	rega_quad_graph = 1;

	// relations to optimize
	rega_quad_optimize_para = true;
	rega_quad_optimize_ortho = true;

	rega_quad_distance_considered = false;

	// delaunay graph : discretize points regularly and gets neighboring relationships
	rega_quad_discretize = true;
	rega_quad_discretization_step = 10;

	// different ways to set the maximal angle for a rotation
	// three choices : constant theta max, linear function (approx), or a value that depends on orthogonal shift (in pixels)
	// set angle function to 0
	rega_angle_function = 0;
	rega_angle_const = 5;
	rega_angle_offset = 1.8;
}

// see parameters in rega above
void Parameters::reset_regp()
{
	regp_ms_sigma = 4;
	regp_ms_epsilon = 0.125;
	regp_ms_min_terms = 2;
	regp_ms_distx = 4500;
	regp_ms_disty = 6;
	regp_ms_smooth_dist = 3;

	regp_quad_lambda = 0.8;
	regp_quad_distance = 50;

	regp_trans_function = 0;
	regp_trans_const = 1;
}

void Parameters::reset_prop()
{
	prop_policy = 1; 	// different types of algorithm
	prop_ttl = 1; 		// number of intersections
	prop_distance = 50; // max distance in pixels to propogate
	prop_min_edge = 0.00001;
	prop_range = 50;

	prop_extra_enabled = false; // turn it off because it uses image

	prop_region_length = 10;
	prop_region_width = 1;

	prop_sub_region_length = 10;
	prop_sub_region_width = 1;

	prop_compared = 0;
	prop_ratio = 0.8;
	prop_dot_th = 0.05;

	prop_check_m_enabled = true;

	prop_m_compare_to = 0;
	prop_m_factor = 0.5;
	prop_m_fixed_magnitude = 0.05;

	prop_check_t_enabled = true;

	prop_t_compare_to = 0;
	prop_t_tolerance = 12.5;

	merge_enabled = true; 	  // use it to merge thin cells
	merge_min_thinness = 4.0; // if less than this value the cell is assumed to be thin and it is merged
}
