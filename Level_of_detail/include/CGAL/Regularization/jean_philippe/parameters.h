#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>


class Parameters
{
public:
    Parameters();

    ~Parameters();

    void reset();

	void reset_lsd();

	void reset_rega();

	void reset_regp();

	void reset_prop();

public:
    /** Level 0 : user inputs */

	int verbose_level;

    std::string path_input_image;
	std::string path_input_directory;
	std::string path_ground_truth_directory;

	std::string path_output_directory;
	double zoom;
	std::string suffix;
	bool lsd_test_mode;

    /** Level 1 : LSD parameters */

    double lsd_scale;
    double lsd_sigma_scale;
	double lsd_quant;
    double lsd_angle;
    double lsd_log_eps;
    double lsd_density;

    bool lsd_create_additional;
    double lsd_additional_size;

    /** Level 2 : Regularization parameters (1) */

    bool rega_regp_enabled;

    int rega_method;
	double rega_epsilon;

    double rega_ms_sigma;
    double rega_ms_epsilon;
    int rega_ms_distance;
    int rega_ms_min_terms;
    double rega_ms_smooth_dist;

    double rega_quad_lambda;
    int rega_quad_distance;
	int rega_quad_graph;
	bool rega_quad_optimize_para;
	bool rega_quad_optimize_ortho;
	bool rega_quad_distance_considered;
	bool rega_quad_discretize;
	int rega_quad_discretization_step;

    int rega_angle_function;
    double rega_angle_const;
    double rega_angle_offset;

    /** Level 3 : Regularization parameters (2) */

    double regp_ms_sigma;
    double regp_ms_epsilon;
    int regp_ms_distx;
    double regp_ms_disty;
    int regp_ms_min_terms;
    double regp_ms_smooth_dist;

	double regp_quad_lambda;
	int regp_quad_distance;

	int regp_trans_function;
	double regp_trans_const;

    /** Level 4 : Propagation */

    int prop_policy;
    int prop_ttl;
    int prop_distance;
    double prop_min_edge;
    double prop_range;

	bool prop_extra_enabled;

	double prop_region_length;
	double prop_region_width;

	double prop_sub_region_length;
	double prop_sub_region_width;

	int prop_compared;
	double prop_ratio;
	double prop_dot_th;

	bool prop_check_m_enabled;
	int prop_m_compare_to;
	double prop_m_factor;
	double prop_m_fixed_magnitude;

	bool prop_check_t_enabled;
	int prop_t_compare_to;
	double prop_t_tolerance;

	bool merge_enabled;
	double merge_min_thinness;
};

#endif // PARAMETERS_H
