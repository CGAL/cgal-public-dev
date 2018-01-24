#include "kinetic_model.h"
#include <opencv2/highgui/highgui.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include "svg.h"
#include "gdal_priv.h"
#include "cpl_conv.h"

using jp::Matrix;

Kinetic_Model::Kinetic_Model()
{
	GDALAllRegister();

    params = new Parameters();

    basename = "";
	time_string = "";

    I = Matrix<uchar>();
#if NOT_MEASURING_PERFORMANCES
    I_grad_m = Matrix<double>();
    I_grad_t = Matrix<double>();
	I_grad_m_uchar = Matrix<uchar>();
    I_grad_t_uchar = Matrix<uchar>();
#endif
    I_data_size = 0;
    I_data = NULL;

    segments = vector<Segment *>();
    tree = NULL;
#if NOT_MEASURING_PERFORMANCES
    L_lsd = list<LineItem *>();
    L_rega = list<LineItem *>();
	L_ag = list<LineItem *>();
    L_regp = list<LineItem *>();
#endif
	rega_quad_potentials = 0;

    rays = vector<Ray *>();
    schedule = NULL;
    graph = NULL;

#if NOT_MEASURING_PERFORMANCES
    L_prop = list<LineItem *>();
#endif
    elapsed_time_lsd = 0;
    elapsed_time_regularization = 0;
    elapsed_time_building_graph = 0;
}


Kinetic_Model::Kinetic_Model(const Kinetic_Model & m)
{
    generator = m.generator;
    params = m.params;

    basename = m.basename;
	time_string = m.time_string;

    I = m.I;
#if NOT_MEASURING_PERFORMANCES
    I_grad_m = m.I_grad_m;
    I_grad_t = m.I_grad_t;
	I_grad_m_uchar = m.I_grad_m_uchar;
	I_grad_t_uchar = m.I_grad_t_uchar;
#endif
    I_data_size = m.I_data_size;
    I_data = m.I_data;

    segments = std::vector<Segment *>(m.segments.begin(), m.segments.end());
    tree = m.tree;

#if NOT_MEASURING_PERFORMANCES
    L_lsd = list<LineItem *>(m.L_lsd.begin(), m.L_lsd.end());
    L_rega = list<LineItem *>(m.L_rega.begin(), m.L_rega.end());
	L_ag = list<LineItem *>(m.L_ag.begin(), m.L_ag.end());
    L_regp = list<LineItem *>(m.L_regp.begin(), m.L_regp.end());
#endif
	rega_quad_potentials = m.rega_quad_potentials;

    rays = vector<Ray *>(m.rays.begin(), m.rays.end());
    schedule = m.schedule;
    graph = m.graph;
#if NOT_MEASURING_PERFORMANCES
    L_prop = list<LineItem *>(m.L_prop.begin(), m.L_prop.end());
#endif
}


Kinetic_Model::~Kinetic_Model()
{
	delete params;

    clear();
}


void Kinetic_Model::add_line(list<LineItem *> &container, double x1, double y1, double x2, double y2, uchar r, uchar g, uchar b)
{
	LineItem* line = new LineItem(x1, y1, x2, y2, r, g, b);
    container.push_back(line);
}


void Kinetic_Model::set_path_input_image(std::string & _path)
{
    params->path_input_image = _path;
}


void Kinetic_Model::set_time_string()
{
	// Get time identifier
    time_t rawtime;
    struct tm * timeinfo;
    char _time_string [18];
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    strftime (_time_string, 18, "%Y-%m-%d-%H%M%S", timeinfo);
	time_string = std::string(_time_string);
}


void Kinetic_Model::set_basename()
{
    boost::filesystem::path path_processed_image (params->path_input_image);
    basename = path_processed_image.stem().string();
}


#if NOT_MEASURING_PERFORMANCES
void Kinetic_Model::set_gradient_maps()
{
    I_grad_m = Matrix<double>(I.rows, I.cols, 1);
    I_grad_t = Matrix<double>(I.rows, I.cols, 1);

    // 1. Smoothes I and computes luminance at the same time
	uint tmprows = I.rows;
	uint tmpcols = I.cols;
    Matrix<uchar> I_smoothed(tmprows, tmpcols, 3);
    int radius = 2;
    Matrix<double> G;
    Matrix<double>::gaussian_matrix(G, radius);

    for (uint i = 0 ; i < I.rows ; i++) {
        for (uint j = 0 ; j < I.cols ; j++) {
            double s_r = 0, s_g = 0, s_b = 0;
            for (int k = -radius ; k <= radius ; k++) {
                for (int l = -radius ; l <= radius ; l++) {
                    double g_kl = G(k + radius, l + radius);
                    uint i_kl = uint(jclamp(0, int(i) + k, int(I.rows) - 1));
                    uint j_kl = uint(jclamp(0, int(j) + l, int(I.cols) - 1));
                    s_r += g_kl * I(i_kl, j_kl, 0);
                    s_g += g_kl * I(i_kl, j_kl, 1);
                    s_b += g_kl * I(i_kl, j_kl, 2);
                }
            }
            I_smoothed(i, j, 0) = uchar(jclamp(0, s_r, 255));
            I_smoothed(i, j, 1) = uchar(jclamp(0, s_g, 255));
            I_smoothed(i, j, 2) = uchar(jclamp(0, s_b, 255));
        }
    }

    // 2. Computes luminance (between 0 and 1)
    Matrix<double> I_lum(I.rows, I.cols, 1);
    for (uint i = 0 ; i < I.rows ; i++) {
        for (uint j = 0 ; j < I.cols ; j++) {
            I_lum(i, j) = (0.2126 * I_smoothed(i, j, 0) + 0.7152 * I_smoothed(i, j, 1) + 0.0722 * I_smoothed(i, j, 2)) / 255.0;
        }
    }
    I_smoothed.release();

    // 3. Computes gradient
    int sobel_x[3][3] = { {1, 0, -1}, {2, 0, -2}, {1, 0, -1}};
    int sobel_y[3][3] = { {1, 2, 1}, {0, 0, 0}, {-1, -2, -1}};
	int mask_x[2][2] = { {-1, 1}, {-1, 1} };
	int mask_y[2][2] = { {-1, -1}, {1, 1} };
    for (uint i = 0 ; i < I.rows ; i++) {
        for (uint j = 0 ; j < I.cols ; j++) {
            double g_x = 0;
            double g_y = 0;
			/*for (int k = 0 ; k <= 1 ; k++) {
				for (int l = 0 ; l <= 1 ; l++) {
					uint i_kl = uint(jclamp(0, int(i) + k, int(I.rows) - 1));
                    uint j_kl = uint(jclamp(0, int(j) + l, int(I.cols) - 1));
                    g_x += mask_x[k][l] * I_lum(i_kl, j_kl);
                    g_y += mask_y[k][l] * I_lum(i_kl, j_kl);
				}
			}
			I_grad_m(i, j) = sqrt((g_x * g_x + g_y * g_y) / 4);*/
            for (int k = -1 ; k <= 1 ; k++) {
                for (int l = -1 ; l <= 1 ; l++) {
                    uint i_kl = uint(jclamp(0, int(i) + k, int(I.rows) - 1));
                    uint j_kl = uint(jclamp(0, int(j) + l, int(I.cols) - 1));
                    g_x += sobel_x[k + 1][l + 1] * I_lum(i_kl, j_kl);
                    g_y += sobel_y[k + 1][l + 1] * I_lum(i_kl, j_kl);
                }
            }
			I_grad_m(i, j) = sqrt((g_x * g_x + g_y * g_y));
            if (g_x == 0 && g_y == 0) {
                I_grad_t(i, j) = -FLT_MAX;
            } else {
				I_grad_t(i, j) = atan2(g_x, -g_y);
            }
        }
    }

	// 4. Displayable map
	I_grad_m_uchar = Matrix<uchar>(I_grad_m.rows, I_grad_m.cols, 3);
	I_grad_t_uchar = Matrix<uchar>(I_grad_t.rows, I_grad_t.cols, 3);
	double h, s = 1, v;
    uchar r, g, b;
    for (uint i = 0 ; i < I_grad_t.rows ; i++) {
        for (uint j = 0 ; j < I_grad_t.cols ; j++) {

            double m_ij = I_grad_m(i, j);
            double t_ij = I_grad_t(i, j);

            // I_m : magnitude in black and white
            //r = uchar(jclamp(0, 4 * 255 * m_ij, 255));
			r = uchar(jclamp(0, 255 * m_ij, 255));
            I_grad_m_uchar (i, j, 0) = I_grad_m_uchar (i, j, 1) = I_grad_m_uchar (i, j, 2) = r;

            // I_t : magnitude and angle together
            if (t_ij != -FLT_MAX) {
                h = 180 * (t_ij / PI + 1);
                //v = jclamp(0, 4 * m_ij, 1);
				v = jclamp(0, m_ij, 1);
                hsv_to_rgb(h, s, v, r, g, b);
                I_grad_t_uchar(i, j, 0) = r;
                I_grad_t_uchar(i, j, 1) = g;
                I_grad_t_uchar(i, j, 2) = b;
            }
        }
    }
}
#endif


void Kinetic_Model::hsv_to_rgb(double & h, double & s, double & v, uchar & r, uchar & g, uchar & b)
{
    int t_i = int(h / 60) % 6;
    double f = h / 60.0 - t_i;
    double l = v * (1 - s);
    double m = v * (1 - f * s);
    double n = v * (1 - (1 - f) * s);
    switch (t_i) {
    case 0: r = uchar(255 * v); g = uchar(255 * n); b = uchar(255 * l); break;
    case 1: r = uchar(255 * m); g = uchar(255 * v); b = uchar(255 * l); break;
    case 2: r = uchar(255 * l); g = uchar(255 * v); b = uchar(255 * n); break;
    case 3: r = uchar(255 * l); g = uchar(255 * m); b = uchar(255 * v); break;
    case 4: r = uchar(255 * n); g = uchar(255 * l); b = uchar(255 * v); break;
    case 5: r = uchar(255 * v); g = uchar(255 * l); b = uchar(255 * m); break;
    }
}


void Kinetic_Model::clear_line_items(list<LineItem *> & objects)
{
    for (list<LineItem *>::iterator it_l = objects.begin() ; it_l != objects.end() ; it_l++) {
        delete (*it_l);
    }
    objects.clear();
}


void Kinetic_Model::set_lsd_scale(double z)
{
    params->lsd_scale = z;
}

void Kinetic_Model::set_lsd_sigma_scale(double z)
{
    params->lsd_sigma_scale = z;
}

void Kinetic_Model::set_lsd_quant(double z)
{
    params->lsd_quant = z;
}

void Kinetic_Model::set_lsd_angle(double z)
{
    params->lsd_angle = z;
}

void Kinetic_Model::set_lsd_log_eps(double z)
{
    params->lsd_log_eps = z;
}

void Kinetic_Model::set_lsd_density(double z)
{
    params->lsd_density = z;
}

void Kinetic_Model::set_lsd_create_additional(bool z)
{
    params->lsd_create_additional = z;
}

void Kinetic_Model::set_lsd_additional_size(double z)
{
    params->lsd_additional_size = z;
}


void Kinetic_Model::set_rega_method(int z)
{
    params->rega_method = z;
}

void Kinetic_Model::set_rega_epsilon(double z)
{
	params->rega_epsilon = z;
}

void Kinetic_Model::set_rega_ms_sigma(double z)
{
    params->rega_ms_sigma = z;
}

void Kinetic_Model::set_rega_ms_epsilon(double z)
{
    params->rega_ms_epsilon = z;
}

void Kinetic_Model::set_rega_ms_distance(int z)
{
    params->rega_ms_distance = z;
}

void Kinetic_Model::set_rega_ms_min_terms(int z)
{
    params->rega_ms_min_terms = z;
}

void Kinetic_Model::set_rega_ms_smooth_dist(double z)
{
    params->rega_ms_smooth_dist = z;
}

void Kinetic_Model::set_rega_quad_lambda(double z)
{
    params->rega_quad_lambda = z;
}

void Kinetic_Model::set_rega_quad_distance(int z)
{
    params->rega_quad_distance = z;
}

void Kinetic_Model::set_rega_quad_graph(int z)
{
	params->rega_quad_graph = z;
}

void Kinetic_Model::set_rega_quad_optimize_para(bool z)
{
	params->rega_quad_optimize_para = z;
}

void Kinetic_Model::set_rega_quad_optimize_ortho(bool z)
{
	params->rega_quad_optimize_ortho = z;
}

void Kinetic_Model::set_rega_quad_distance_considered(bool z)
{
	params->rega_quad_distance_considered = z;
}

void Kinetic_Model::set_rega_quad_discretize(bool z)
{
	params->rega_quad_discretize = z;
}

void Kinetic_Model::set_rega_quad_discretization_step(int z)
{
	params->rega_quad_discretization_step = z;
}


void Kinetic_Model::set_rega_angle_function(int z)
{
    params->rega_angle_function = z;
}

void Kinetic_Model::set_rega_angle_const(double z)
{
    params->rega_angle_const = z;
}

void Kinetic_Model::set_rega_angle_offset(double z)
{
    params->rega_angle_offset = z;
}

void Kinetic_Model::set_regp_ms_sigma(double z)
{
    params->regp_ms_sigma = z;
}

void Kinetic_Model::set_regp_ms_epsilon(double z)
{
    params->regp_ms_epsilon = z;
}

void Kinetic_Model::set_regp_ms_distx(int z)
{
    params->regp_ms_distx = z;
}

void Kinetic_Model::set_regp_ms_disty(double z)
{
    params->regp_ms_disty = z;
}

void Kinetic_Model::set_regp_ms_min_terms(int z)
{
    params->regp_ms_min_terms = z;
}

void Kinetic_Model::set_regp_ms_smooth_dist(double z)
{
    params->regp_ms_smooth_dist = z;
}

void Kinetic_Model::set_regp_quad_lambda(double z)
{
	params->regp_quad_lambda = z;
}

void Kinetic_Model::set_regp_quad_distance(int z)
{
	params->regp_quad_distance = z;
}

void Kinetic_Model::set_regp_trans_function(int z)
{
	params->regp_trans_function = z;
}

void Kinetic_Model::set_regp_trans_const(double z)
{
	params->regp_trans_const = z;
}

void Kinetic_Model::set_prop_policy(int z)
{
    params->prop_policy = z;
}

void Kinetic_Model::set_prop_ttl(int z)
{
    params->prop_ttl = z;
}

void Kinetic_Model::set_prop_distance(int z)
{
    params->prop_distance = z;
}

void Kinetic_Model::set_prop_min_edge(double z)
{
    params->prop_min_edge = z;
}

void Kinetic_Model::set_prop_range(int z)
{
    params->prop_range = double(z);
}

void Kinetic_Model::set_prop_extra_enabled(bool z)
{
	params->prop_extra_enabled = z;
}

void Kinetic_Model::set_prop_region_length(double z)
{
	params->prop_region_length = z;
}

void Kinetic_Model::set_prop_region_width(double z)
{
	params->prop_region_width = z;
}

void Kinetic_Model::set_prop_sub_region_length(double z)
{
	params->prop_sub_region_length = z;
}

void Kinetic_Model::set_prop_sub_region_width(double z)
{
	params->prop_sub_region_width = z;
}

void Kinetic_Model::set_prop_compared_quantity(int z)
{
	params->prop_compared = z;
}

void Kinetic_Model::set_prop_ratio(double z)
{
	params->prop_ratio = z;
}

void Kinetic_Model::set_prop_dot_th(double z)
{
	params->prop_dot_th = z;
}

void Kinetic_Model::set_prop_check_m_enabled(bool z)
{
	params->prop_check_m_enabled = z;
}

void Kinetic_Model::set_prop_m_compare_to(int z)
{
	params->prop_m_compare_to = z;
}

void Kinetic_Model::set_prop_m_factor(double z)
{
	params->prop_m_factor = z;
}

void Kinetic_Model::set_prop_m_fixed_magnitude(double z)
{
	params->prop_m_fixed_magnitude = z;
}

void Kinetic_Model::set_prop_check_t_enabled(bool z)
{
	params->prop_check_t_enabled = z;
}

void Kinetic_Model::set_prop_t_compare_to(int z)
{
	params->prop_t_compare_to = z;
}

void Kinetic_Model::set_prop_t_tolerance(double z)
{
	params->prop_t_tolerance = z;
}

void Kinetic_Model::set_prop_merge_enabled(bool z)
{
	params->merge_enabled = z;
}

void Kinetic_Model::set_prop_merge_min_thinness(double z)
{
	params->merge_min_thinness = z;
}


void Kinetic_Model::reinit()
{
    clear();
    init();
}


void Kinetic_Model::clear()
{
#if NOT_MEASURING_PERFORMANCES
    clear_line_items(L_prop);
#endif
    if (graph != NULL) {
        delete graph;
        graph = NULL;
    }
    Ray::clear_rays(rays);
    //IndexedEvent::clear_schedule(schedule);

#if NOT_MEASURING_PERFORMANCES
    clear_line_items(L_regp);
    clear_line_items(L_rega);
	clear_line_items(L_ag);
    clear_line_items(L_lsd);
#endif
    I.release();
    delete_byte_array(I_data_size, I_data);

    if (tree != NULL) {
        delete tree;
        tree = NULL;
    }
    Segment::clear_segments(segments);
}


void Kinetic_Model::switch_to_base_matrix()
{
	delete_byte_array(I_data_size, I_data);

	GDALDataset *poDataset = (GDALDataset *)GDALOpen(params->path_input_image.c_str(), GA_ReadOnly);
	if (poDataset != NULL) {
		I = Matrix<uchar>(poDataset);
	} else {
		std::cerr << "Error : GDAL couldn't load the specified file" << std::endl;
		exit(-1);
	}
	GDALClose(poDataset);
}


void Kinetic_Model::init()
{
	// Dmitry Fills I_data
	/* 
	GDALDataset* im_dataset = (GDALDataset *)GDALOpen(params->path_input_image.c_str(), GA_ReadOnly);
	if (im_dataset == NULL) {
		throw std::logic_error("Error : a problem occured while loading the image");
	}
	I_data_cols = im_dataset->GetRasterXSize();
	I_data_rows = im_dataset->GetRasterYSize();
	int channels = im_dataset->GetRasterCount();

	reallocate_byte_array<double>(I_data_rows * I_data_cols, I_data_size, I_data);
	if (channels == 1) {
		// Case of a grayscale image
		// We read all the rows of the image, and then we sequentially copy these rows to I_data
		GDALRasterBand* im_grayscale_band = im_dataset->GetRasterBand(1);
		double* buffer = (double*)CPLMalloc(I_data_cols * sizeof(double));
		for (int i = 0 ; i < I_data_rows ; i++) {
			im_grayscale_band->RasterIO(GF_Read, 0, i, I_data_cols, 1, buffer, I_data_cols, 1, GDT_Float64, 0, 0);
			memcpy(&(I_data[i * I_data_cols]), buffer, I_data_cols * sizeof(double));
		}
		CPLFree(buffer);

	} else if (channels == 3) {
		// Case of a RGB image
		// It is almost the same process as before, except that we must compute the luminance
		GDALRasterBand* im_r_band = im_dataset->GetRasterBand(1);
		GDALRasterBand* im_g_band = im_dataset->GetRasterBand(2);
		GDALRasterBand* im_b_band = im_dataset->GetRasterBand(3);
		double* buffer_r = (double*)CPLMalloc(I_data_cols * sizeof(double));
		double* buffer_g = (double*)CPLMalloc(I_data_cols * sizeof(double));
		double* buffer_b = (double*)CPLMalloc(I_data_cols * sizeof(double));
		double* buffer = (double*)CPLMalloc(I_data_cols * sizeof(double));
		for (int i = 0 ; i < I_data_rows ; i++) {
			im_r_band->RasterIO(GF_Read, 0, i, I_data_cols, 1, buffer_r, I_data_cols, 1, GDT_Float64, 0, 0);
			im_g_band->RasterIO(GF_Read, 0, i, I_data_cols, 1, buffer_g, I_data_cols, 1, GDT_Float64, 0, 0);
			im_b_band->RasterIO(GF_Read, 0, i, I_data_cols, 1, buffer_b, I_data_cols, 1, GDT_Float64, 0, 0);
			for (int j = 0 ; j < I_data_cols ; j++) {
				buffer[j] = 0.2126 * buffer_r[j] + 0.7152 * buffer_g[j] + 0.0722 * buffer_b[j];
			}
			memcpy(&(I_data[i * I_data_cols]), buffer, I_data_cols * sizeof(double));
		}
		CPLFree(buffer_r);
		CPLFree(buffer_g);
		CPLFree(buffer_b);
		CPLFree(buffer);

	} else {
		throw std::logic_error("Error : the image should have one or three channels");
	}
	GDALClose(im_dataset);

#if NOT_MEASURING_PERFORMANCES
	GDALDataset *poDataset = (GDALDataset *)GDALOpen(params->path_input_image.c_str(), GA_ReadOnly);
	if (poDataset != NULL) {
		I = Matrix<uchar>(poDataset);
	} else {
		std::cerr << "Error : GDAL couldn't load the specified file" << std::endl;
		exit(-1);
	}
	GDALClose(poDataset);
#endif
	*/

    // Initializes Kinetic_Model by opening the image
    /*Mat I_0 = cv::imread(params->path_input_image, CV_LOAD_IMAGE_GRAYSCALE);
    if (I_0.empty()) {
        throw std::logic_error("Error : a problem occured while loading the image");
    } else if (I_0.channels() != 3) {
        I_0 = Mat();
        throw std::logic_error("Error : the image should have three channels");
    }*/
	
	// Dmitry
	/*
	set_time_string();
    set_basename(); */

	
    // Initializes data needed for running the LSD
    /*I = Matrix<uchar>(I_0.rows, I_0.cols, 3);
    for (uint i = 0 ; i < I.rows ; i++) {
        for (uint j = 0 ; j < I.cols ; j++) {
            for (uint c = 0 ; c < I.channels ; c++) {
                I(i, j, c) = uchar(I_0.at<cv::Vec3b>(i, j)[2 - c]);
            }
        }
    }*/

    /*reallocate_byte_array<double>(I_0.rows * I_0.cols, I_data_size, I_data);
	I_data_rows = I_0.rows;
	I_data_cols = I_0.cols;
	for (int i = 0 ; i < I_0.rows * I_0.cols; i++) {
		I_data[i] = double(I_0.data[i]);
	}
	I_0.release();*/

    /*reallocate_byte_array<double>(I.rows * I.cols, I_data_size, I_data);
	for (uint r = 0 ; r < I.rows ; r++) {
        for (uint c = 0 ; c < I.cols ; c++) {
            //I_data[r * I.cols + c] = 0.2126 * I(r, c, 0) + 0.7152 * I(r, c, 1) + 0.0722 * I(r, c, 2);
        }
    }*/

#if NOT_MEASURING_PERFORMANCES
    set_gradient_maps();
#endif

#if NOT_MEASURING_PERFORMANCES
    L_lsd = list<LineItem *>();
#endif
    segments = vector<Segment *>();
    tree = new Segment_Regularization_Tree();

#if NOT_MEASURING_PERFORMANCES
	L_ag = list<LineItem *>();
    L_rega = list<LineItem *>();
#endif
    applied_regularization_angles = false;

#if NOT_MEASURING_PERFORMANCES
    L_regp = list<LineItem *>();
#endif
    applied_regularization_ordinates = false;

    rays = vector<Ray *>();
    schedule = NULL;
    graph = NULL;
#if NOT_MEASURING_PERFORMANCES
    L_prop = list<LineItem *>();
#endif
    elapsed_time_lsd = 0;
    elapsed_time_regularization = 0;
    elapsed_time_building_graph = 0;
}

#if NOT_MEASURING_PERFORMANCES
void Kinetic_Model::segments_to_svg(std::string & directory, int r, int g, int b)
{
	std::string filename = directory + "\\" + basename + "_segments.svg";
	std::filebuf fb;

	fb.open(filename, std::ios::out);
	std::ostream os(&fb);
	
	SVG::markup_header(os, int(I.rows), int(I.cols));
	SVG::print_line_items(os, L_lsd, r, g, b);
	SVG::markup_footer(os);

	fb.close();
}


void Kinetic_Model::segments_to_svg(std::string & directory, int i_min, int i_max, int j_min, int j_max, int r, int g, int b)
{
	int rows = i_max - i_min + 1;
	int cols = j_max - j_min + 1;
	
	// First part : crop the image and save it
	Matrix<uchar> J(I, rows, cols, i_min, j_min);
	std::string filename_crop = directory + "\\" + basename + "_segments_crop.tiff";
	J.write_uchar(filename_crop);

	// Second part : the SVG file
	std::string filename = directory + "\\" + basename + "_segments.svg";
	std::filebuf fb;

	fb.open(filename, std::ios::out);
	std::ostream os(&fb);

	SVG::markup_header(os, rows, cols);
	SVG::markup_image(os, filename_crop, rows, cols, 1);

	double xa = j_min, xb = j_max + 1;
	double ya = int(I.rows) - (i_max + 1), yb = int(I.rows) - i_min;

	SVG::print_line_items(os, L_lsd, int(I.rows), xa, xb, ya, yb, r, g, b);

	SVG::markup_footer(os);
	fb.close();
}


void Kinetic_Model::rega_to_svg(std::string & directory)
{
	std::string filename = directory + "\\" + basename + "_rega.svg";
	std::filebuf fb;

	fb.open(filename, std::ios::out);
	std::ostream os(&fb);
	
	SVG::markup_header(os, int(I.rows), int(I.cols));
	SVG::print_line_items(os, L_rega);
	SVG::markup_footer(os);

	fb.close();
}


void Kinetic_Model::rega_to_svg(std::string & directory, int i_min, int i_max, int j_min, int j_max)
{
	int rows = i_max - i_min + 1;
	int cols = j_max - j_min + 1;
	
	// First part : crop the image and save it
	Matrix<uchar> J(I, rows, cols, i_min, j_min);
	std::string filename_crop = directory + "\\" + basename + "_rega_crop.tiff";
	J.write_uchar(filename_crop);

	// Second part : the SVG file
	std::string filename = directory + "\\" + basename + "_rega.svg";
	std::filebuf fb;

	fb.open(filename, std::ios::out);
	std::ostream os(&fb);

	SVG::markup_header(os, rows, cols);
	SVG::markup_image(os, filename_crop, rows, cols, 1);

	double xa = j_min, xb = j_max + 1;
	double ya = int(I.rows) - (i_max + 1), yb = int(I.rows) - i_min;

	SVG::print_line_items(os, L_rega, int(I.rows), xa, xb, ya, yb);

	SVG::markup_footer(os);
	fb.close();
}


void Kinetic_Model::regp_to_svg(std::string & directory)
{
	std::string filename = directory + "\\" + basename + "_regp.svg";
	std::filebuf fb;

	fb.open(filename, std::ios::out);
	std::ostream os(&fb);
	
	SVG::markup_header(os, int(I.rows), int(I.cols));
	SVG::print_line_items(os, L_regp);
	SVG::markup_footer(os);

	fb.close();
}


void Kinetic_Model::regp_to_svg(std::string & directory, int i_min, int i_max, int j_min, int j_max)
{
	int rows = i_max - i_min + 1;
	int cols = j_max - j_min + 1;
	
	// First part : crop the image and save it
	Matrix<uchar> J(I, rows, cols, i_min, j_min);
	std::string filename_crop = directory + "\\" + basename + "_regp_crop.tiff";
	J.write_uchar(filename_crop);

	// Second part : the SVG file
	std::string filename = directory + "\\" + basename + "_regp.svg";
	std::filebuf fb;

	fb.open(filename, std::ios::out);
	std::ostream os(&fb);

	SVG::markup_header(os, rows, cols);
	SVG::markup_image(os, filename_crop, rows, cols, 1);

	double xa = j_min, xb = j_max + 1;
	double ya = int(I.rows) - (i_max + 1), yb = int(I.rows) - i_min;

	SVG::print_line_items(os, L_regp, int(I.rows), xa, xb, ya, yb);

	SVG::markup_footer(os);
	fb.close();
}


void Kinetic_Model::partition_to_svg(std::string & directory)
{
	std::string filename = directory + "\\" + basename + "_partition_" + std::to_string(graph->faces.size()) + ".svg";
	std::filebuf fb;

	fb.open(filename, std::ios::out);
	std::ostream os(&fb);
	
	SVG::markup_header(os, int(I.rows), int(I.cols));
	SVG::print_line_items(os, L_prop);
	SVG::markup_footer(os);

	fb.close();
}


void Kinetic_Model::partition_to_svg(std::string & directory, int i_min, int i_max, int j_min, int j_max)
{
	int rows = i_max - i_min + 1;
	int cols = j_max - j_min + 1;
	
	// First part : crop the image and save it
	Matrix<uchar> J(I, rows, cols, i_min, j_min);
	std::string filename_crop = directory + "\\" + basename + "_partition_crop.tiff";
	J.write_uchar(filename_crop);

	// Second part : the SVG file
	std::string filename = directory + "\\" + basename + "_partition.svg";
	std::filebuf fb;

	fb.open(filename, std::ios::out);
	std::ostream os(&fb);

	SVG::markup_header(os, rows, cols);
	SVG::markup_image(os, filename_crop, rows, cols, 1);

	double xa = j_min, xb = j_max + 1;
	double ya = int(I.rows) - (i_max + 1), yb = int(I.rows) - i_min;

	SVG::print_line_items(os, L_prop, int(I.rows), xa, xb, ya, yb);

	SVG::markup_footer(os);
	fb.close();
}


void Kinetic_Model::harlequin_to_svg(std::string & directory)
{
	if (graph == nullptr || !graph->is_valid) return;

	std::string filename = directory + "\\" + basename + "_harlequin_" + std::to_string(graph->faces.size()) + ".svg";
	std::filebuf fb;

	fb.open(filename, std::ios::out);
	std::ostream os(&fb);
	
	SVG::markup_header(os, int(I.rows), int(I.cols));
	SVG::markup_image(os, params->path_input_image, int(I.rows), int(I.cols), 0.6);

	std::default_random_engine generator;
	std::uniform_int_distribution<int> uniform_dist (128, 255);
	
	for (list<Face*>::iterator it_f = graph->faces.begin() ; it_f != graph->faces.end() ; it_f++) {
		int r = uniform_dist(generator);
		int g = uniform_dist(generator);
		int b = uniform_dist(generator);
		SVG::markup_polygon(os, int(I.rows), int(I.cols), (*it_f), r, g, b, 2, 0.4);
	}
	SVG::print_line_items(os, L_prop);
	SVG::markup_footer(os);

	fb.close();
}
#endif