#ifndef MODEL_H
#define MODEL_H

#include <random>
#include <vector>
#include <opencv2/core/core.hpp>
#include <set>

using std::default_random_engine;
using std::vector;
using std::set;
using cv::Mat;

#include "parameters.h"
#include "matrix.h"
#include "segment_ray.h"
#include "segment_tree.h"
#include "indexed_event.h"
#include "partition.h"
#include "line_item.h"



class Kinetic_Model
{
public:
    Kinetic_Model();

    Kinetic_Model(const Kinetic_Model & m);

    virtual ~Kinetic_Model();


    void set_path_input_image(std::string & _path);

	void switch_to_base_matrix();

	void set_time_string();

	void set_basename();

    void reinit();

	void add_line(list<LineItem *> & container, double x1, double y1, double x2, double y2, uchar r, uchar g, uchar b);

	void clear_line_items(list<LineItem *> & objects);

private:
    void clear();

    void init();

#if NOT_MEASURING_PERFORMANCES
    void set_gradient_maps();
#endif

	void hsv_to_rgb(double & h, double & s, double & v, uchar & r, uchar & g, uchar & b);

public:
    void set_lsd_scale(double z);
    void set_lsd_sigma_scale(double z);
    void set_lsd_quant(double z);
    void set_lsd_angle(double z);
    void set_lsd_log_eps(double z);
    void set_lsd_density(double z);
    void set_lsd_create_additional(bool z);
    void set_lsd_additional_size(double z);

    void set_rega_method(int z);
	void set_rega_epsilon(double z);
    void set_rega_ms_sigma(double z);
    void set_rega_ms_epsilon(double z);
    void set_rega_ms_distance(int z);
    void set_rega_ms_min_terms(int z);
    void set_rega_ms_smooth_dist(double z);

    void set_rega_quad_lambda(double z);
    void set_rega_quad_distance(int z);
	void set_rega_quad_graph(int z);
	void set_rega_quad_optimize_para(bool z);
	void set_rega_quad_optimize_ortho(bool z);
	void set_rega_quad_distance_considered(bool z);
	void set_rega_quad_discretize(bool z);
	void set_rega_quad_discretization_step(int z);

	void set_rega_angle_function(int z);
    void set_rega_angle_const(double z);
    void set_rega_angle_offset(double z);

    void set_regp_ms_sigma(double z);
    void set_regp_ms_epsilon(double z);
    void set_regp_ms_distx(int z);
    void set_regp_ms_disty(double z);
    void set_regp_ms_min_terms(int z);
    void set_regp_ms_smooth_dist(double z);
	void set_regp_quad_lambda(double z);
    void set_regp_quad_distance(int z);
	void set_regp_trans_function(int z);
    void set_regp_trans_const(double z);

    void set_prop_policy(int z);
    void set_prop_ttl(int z);
    void set_prop_distance(int z);
    void set_prop_min_edge(double z);
    void set_prop_range(int z);

	void set_prop_extra_enabled(bool z);

	void set_prop_region_length(double z);
	void set_prop_region_width(double z);
	void set_prop_sub_region_length(double z);
	void set_prop_sub_region_width(double z);
	void set_prop_compared_quantity(int z);
	void set_prop_ratio(double z);
	void set_prop_dot_th(double z);

	void set_prop_check_m_enabled(bool z);
	void set_prop_m_compare_to(int z);
	void set_prop_m_factor(double z);
	void set_prop_m_fixed_magnitude(double z);

	void set_prop_check_t_enabled(bool z);
	void set_prop_t_compare_to(int z);
	void set_prop_t_tolerance(double z);

	void set_prop_merge_enabled(bool z);
	void set_prop_merge_min_thinness(double z);

private:
    template <typename T>
    void reallocate_byte_array(uint size, uint & byte_array_size, T* & byte_array)
    {
        if (size != byte_array_size) {
            if (byte_array != NULL) delete[] byte_array;
            byte_array_size = size;
            byte_array = new T[byte_array_size];
        }
    }

    template <typename T>
    void delete_byte_array(uint & byte_array_size, T* & byte_array)
    {
        if (byte_array_size != 0) {
            delete[] byte_array;
            byte_array = NULL;
            byte_array_size = 0;
        }
    }

public:
    default_random_engine generator;
    Parameters* params;

    /** Level 0 : user inputs */

    std::string basename;
	std::string time_string;

    /** Level 1 : LSD */

    Matrix<uchar> I;
#if NOT_MEASURING_PERFORMANCES
    Matrix<double> I_grad_m;
    Matrix<double> I_grad_t;
	Matrix<uchar> I_grad_m_uchar;
	Matrix<uchar> I_grad_t_uchar;
#endif
    uint I_data_size;
    double* I_data;
	int I_data_rows;
	int I_data_cols;

    vector<Segment *> segments;
    Segment_Regularization_Tree* tree;

    list<LineItem *> L_lsd;

    /** Level 2 : Regularization (1) */

#if NOT_MEASURING_PERFORMANCES
    list<LineItem *> L_rega;
	list<LineItem *> L_ag;
#endif
	int rega_quad_potentials;
    bool applied_regularization_angles;

    /** Level 3 : Regularization (2) */

#if NOT_MEASURING_PERFORMANCES
    list<LineItem *> L_regp;
#endif
    bool applied_regularization_ordinates;

    /** Level 4 : Propagation */

    vector<Ray *> rays;
    IndexedEvent* schedule;
    Partition* graph;
#if NOT_MEASURING_PERFORMANCES
    list<LineItem *> L_prop;
#endif
    /** Level 5 : Splitting cells */

    double elapsed_time_lsd;
    double elapsed_time_regularization;
    double elapsed_time_building_graph;

#if NOT_MEASURING_PERFORMANCES
	void segments_to_svg(std::string & directory, int r, int g, int b);
	void segments_to_svg(std::string & directory, int i_min, int i_max, int j_min, int j_max, int r, int g, int b);

	void rega_to_svg(std::string & directory);
	void rega_to_svg(std::string & directory, int i_min, int i_max, int j_min, int j_max);

	void regp_to_svg(std::string & directory);
	void regp_to_svg(std::string & directory, int i_min, int i_max, int j_min, int j_max);

	void partition_to_svg(std::string & directory);
	void partition_to_svg(std::string & directory, int i_min, int i_max, int j_min, int j_max);

	void harlequin_to_svg(std::string & directory);
#endif
};

#endif // MODEL_H
