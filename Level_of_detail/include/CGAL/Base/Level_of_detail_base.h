#ifndef CGAL_LEVEL_OF_DETAIL_BASE_H
#define CGAL_LEVEL_OF_DETAIL_BASE_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\"
#define PN "\r\n"
#else 
#define PSR "/" 
#define PN "\n"
#endif

// STL includes.
#include <map>
#include <memory>
#include <string>
#include <iostream>
#include <vector>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/Timer.h>

// New CGAL includes.
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Mylog/Mylog.h>

namespace CGAL {

	namespace LOD {

		// Facade class that accumulates all necessary objects and operations
		// related to the level of detail (LOD) reconstruction.
		template<class LodTraits>
		class Level_of_detail_base {

		public:

			// Main typedefs.
			typedef LodTraits 				Traits;
			typedef typename Traits::Kernel Kernel;

			typedef typename Traits::Container_2D Container_2D;
			typedef typename Traits::Container_3D Container_3D;

			typedef typename Traits::Loader       Loader;
			typedef typename Traits::Preprocessor Preprocessor;

			typedef typename Kernel::FT 	   FT;
			typedef typename Kernel::Point_2   Point_2;
			typedef typename Kernel::Point_3   Point_3;
			typedef typename Kernel::Plane_3   Plane_3;
			typedef typename Kernel::Line_2    Line_2;
			typedef typename Kernel::Segment_2 Segment_2;

			typedef typename Traits::Building_boundary_selector Building_boundary_selector;
			typedef typename Traits::Building_interior_selector Building_interior_selector;
			typedef typename Traits::Clutter_selector 		    Clutter_selector;
			typedef typename Traits::Ground_selector 		    Ground_selector;

			typedef typename Traits::Vertical_regularizer Vertical_regularizer;
			typedef typename Traits::Ground_projector     Ground_projector;
			typedef typename Traits::Projected_points     Projected_points;
			typedef typename Traits::Planes        		  Planes;

			typedef typename Traits::Line_regularizer Line_regularizer;

			typedef typename Traits::Building_splitter Building_splitter;
			typedef typename Traits::Building_outliner Building_outliner;

			typedef typename Traits::Building_min_roof_fitter Building_min_roof_fitter;
			typedef typename Traits::Building_avg_roof_fitter Building_avg_roof_fitter;
			typedef typename Traits::Building_max_roof_fitter Building_max_roof_fitter;

			typedef Planes Boundary_data;

			typedef typename Traits::Structuring_2 	  Structuring_2;
			typedef typename Traits::Visibility_2  	  Visibility_2;
			typedef typename Traits::Region_growing_2 Region_growing_2;
			
			typedef typename Traits::Utils Utils;
			
			typedef typename Traits::CDT        		CDT;
			typedef typename CDT::Vertex_handle 		Vertex_handle;
			typedef typename CDT::Face_handle   		Face_handle;
			typedef typename CDT::Finite_edges_iterator Edge_iterator;
			typedef typename CDT::Finite_faces_iterator Face_iterator;

			typedef typename Traits::Graph_cut Graph_cut;
			typedef typename Traits::Lods 	   Lods;

			typedef typename Traits::Mesh 			   Mesh;
			typedef typename Traits::Mesh_facet_colors Mesh_facet_colors;

			typedef typename Traits::Buildings Buildings;

			typedef typename Lods::Point  Ground_point;
			typedef typename Lods::Ground Ground;

			typedef typename Traits::Grid_simplifier   Grid_simplifier;
			typedef typename Traits::Thinning 		   Thinning;
			typedef typename Traits::Clutter_filtering Clutter_filtering;

			typedef typename Traits::Clutter_processor Clutter_processor;
			typedef Thinning_fitter_type    		   Clutter_fitter_type;
			typedef Grid_new_point_type 			   Clutter_new_point_type;

			typedef typename Traits::Level_of_detail_parameters Parameters_wrapper;
			typedef typename Traits::Parameters 		  		Parameters;
			typedef typename Traits::Parameters_estimator 		Parameters_estimator;

			typedef typename Traits::Lod_complexity Lod_complexity;
			typedef typename Traits::Lod_distortion Lod_distortion;
			typedef typename Traits::Lod_coverage   Lod_coverage;

			typedef typename Traits::Polygonizer 			   Polygonizer;
			typedef typename Traits::Lod_data_structure        Data_structure;
			typedef typename Traits::Polygon_based_visibility  Polygon_based_visibility;
			typedef typename Traits::Inside_buildings_selector Inside_buildings_selector;
			
			typedef typename Traits::Region_growing_3  Region_growing_3;
			typedef typename Traits::Roof_cleaner      Roof_cleaner;
			typedef typename Traits::Partition_input   Partition_input;
			typedef typename Traits::Partition_creator Partition_creator;
			typedef typename Traits::Roof_estimator    Roof_estimator;
			
			typedef typename Traits::LOD2_reconstruction LOD2_reconstruction;


			// Extra typedefs.
			using Plane_iterator = typename Planes::const_iterator;

			using Index   = int;
			using Indices = std::vector<Index>;

			using Structured_points  = std::vector< std::vector<Point_2> >; 			  
			using Structured_labels  = std::vector< std::vector<Structured_label> >;  
			using Structured_anchors = std::vector< std::vector<std::vector<int> > >;
			
			using Log 		  = CGAL::LOD::Mylog;
			using Segment_map = CGAL::Identity_property_map<Segment_2>;

			using Lines    = std::vector<Line_2>;
			using Segments = std::vector<Segment_2>;

			using Label     = typename Traits::Label;
			using Label_map = typename Container_3D:: template Property_map<Label>;

			using Point_index 	  = typename Container_3D::Index;
			using Face_points_map = std::map<Face_handle, std::vector<Point_index> >;

			using Box = typename Lod_distortion::Bounding_box;

			enum class Program_version  { VER0 };
			enum class Pipeline_version { WITH_SHAPE_DETECTION, WITHOUT_SHAPE_DETECTION };


			//////////////
			// Main class with all default parameters.
			Level_of_detail_base() :
			m_prefix_path("default"),
			m_default_path("default"),
			m_preprocessor_scale(-FT(1)) ,
			m_structuring_epsilon(-FT(1)),
			m_structuring_log(false),
			m_structuring_resample(true),
			m_structuring_get_all_points(false),
			m_add_cdt_clutter(true),
			m_add_cdt_bbox(false),
			m_visibility_save_info(false),
			m_visibility_approach(Visibility_approach::POINT_BASED),
			m_visibility_method(Visibility_method::POINT_BASED_CLASSIFICATION),
			m_visibility_num_samples(0),
			m_graph_cut_save_info(false),
			m_graph_cut_alpha(-FT(1)),
			m_graph_cut_beta(-FT(1)),
			m_graph_cut_gamma(-FT(1)),
			m_building_boundaries_save_internal_info(false), 
			m_building_boundaries_max_inner_iters(0),
			m_building_boundaries_max_outer_iters(0),
			m_roof_fitter_type(Roof_fitter_type::MAX), 
			m_clean_projected_points(true),
			m_max_reg_angle(-FT(1)),
			m_regularizer_reject_planes(true), 
			m_use_boundaries(true),
			m_prog_version(Program_version::VER0),
			m_pipeline_version(Pipeline_version::WITHOUT_SHAPE_DETECTION),
			m_visibility_show_progress(false),
			m_visibility_norm_threshold(-FT(1)),
			m_clutter_knn(0),
			m_clutter_cell_length(-FT(1)),
			m_clutter_fitter_type(Clutter_fitter_type::LINE),
			m_clutter_new_point_type(Clutter_new_point_type::BARYCENTRE),
			m_visibility_num_neighbours(0),
			m_visibility_sampler(Visibility_sampler::RANDOM_UNIFORM_0),
			m_visibility_rays_per_side(0),
			m_visibility_small_edge_threshold(FT(0)),
			m_building_boundary_type(Building_boundary_type::UNORIENTED),
			m_visibility_angle_eps(-FT(1)),
			m_thinning_neighbour_search_type(Neighbour_search_type::KNN),
			m_thinning_fuzzy_radius(-FT(1)),
			m_thinning_type(Thinning_type::NAIVE),
			m_region_growing_epsilon_2d(-FT(1)),
			m_region_growing_cluster_epsilon_2d(-FT(1)),
			m_region_growing_normal_threshold_2d(-FT(1)),
			m_region_growing_min_points_2d(0),
			m_with_region_growing(true),
			m_use_grid_simplifier_first(true),
			m_alpha_shape_size(-FT(1)),
			m_use_alpha_shapes(false),
			m_structuring_corner_algorithm(Structuring_corner_algorithm::GRAPH_BASED),
			m_structuring_adjacency_method(Structuring_adjacency_threshold_method::LOCAL),
			m_structuring_adjacency_value(-FT(1)),
			m_structuring_global_everywhere(true),
			m_silent(false),
			m_test_data_type(Main_test_data_type::PARIS_ETH),
			m_region_growing_normal_estimation_method(Region_growing_normal_estimation::PROJECTED),
			m_imp_eps(-FT(1)),
			m_imp_scale(-FT(1)),
			m_estimate_parameters(false),
			m_estimate_quality(false),
			m_complexity(-FT(1)),
			m_distortion(-FT(1)),
			m_coverage(-FT(1)),
			m_clutter_filtering_scale(-FT(1)),
			m_clutter_filtering_mean(-FT(1)),
			m_regularize_lines(false),
			m_polygonize(false),
			m_line_regularizer_optimize_angles(true),
			m_line_regularizer_optimize_ordinates(true),
			m_polygonizer_number_of_intersections(0),
			m_polygonizer_min_face_width(-FT(1)),
			m_building_splitter_use_custom_constraints(true),
			m_building_splitter_constraints_threshold(-FT(1)),
			m_line_regularizer_max_angle_in_degrees(-FT(1)),
			m_line_regularizer_max_difference_in_meters(-FT(1)),
			m_region_growing_epsilon_3d(-FT(1)),
			m_region_growing_cluster_epsilon_3d(-FT(1)),
			m_region_growing_normal_threshold_3d(-FT(1)),
			m_region_growing_min_points_3d(0),
			m_roof_cleaner_scale(-FT(1)),
			m_roof_cleaner_max_percentage(-FT(1))
			{ }


			//////////////////
			// Parameter functions!

			void set_prefix_path(const std::string &new_path) {
				m_prefix_path = new_path;
			}

			void set_data_type(const size_t) {
				assert(!"Should not be used here!");
			}

			void make_silent(const bool) {
				assert(!"Should not be used here!");
			}

			void add_clutter(const bool) {
				assert(!"Should not be used here!");
			}

			void set_clutter_cell_side_length(const FT) {
				assert(!"Should not be used here!");
			}

			void set_region_growing_epsilon(const FT) {
				assert(!"Should not be used here!");
			}

			void set_region_growing_cluster_epsilon(const FT) {
				assert(!"Should not be used here!");
			}

			void set_region_growing_normal_threshold(const FT) {
				assert(!"Should not be used here!");
			}

			void set_region_growing_min_points(const size_t) {
				assert(!"Should not be used here!");
			}

			void set_structuring_epsilon(const FT) {
				assert(!"Should not be used here!");
			}

			void set_structuring_adjacency_value(const FT) {
				assert(!"Should not be used here!");
			}

			void get_all_structuring_points(const bool) {
				assert(!"Should not be used here!");
			}

			void set_graph_cut_beta(const FT) {
				assert(!"Should not be used here!");
			}

			void set_graph_cut_gamma(const FT) {
				assert(!"Should not be used here!");
			}

			void set_default_parameters() {
				
				set_optimal_configuration();
				m_default_path = "/Users/danisimo/Documents/pipeline/data/clean/paris_different/paris_half_tile.ply";
			}


			//////////////////
			// Important public functions!

			FT get_complexity() const {

				assert(m_complexity >= FT(0));
				return m_complexity;
			}

			FT get_distortion() const {

				assert(m_distortion >= FT(0));
				return m_distortion;
			}

			FT get_coverage() const {

				assert(m_coverage >= FT(0));
				return m_coverage;
			}

			std::shared_ptr<Lod_complexity> get_lod_complexity_ptr() const {

				assert(m_lod_complexity != nullptr);
				return m_lod_complexity;
			}

			std::shared_ptr<Lod_distortion> get_lod_distortion_ptr() const {

				assert(m_lod_distortion != nullptr);
				return m_lod_distortion;
			}

			std::shared_ptr<Lod_coverage> get_lod_coverage_ptr() const {

				assert(m_lod_coverage != nullptr);
				return m_lod_coverage;
			}

			FT get_scale() const {
				return m_imp_scale;
			}

			void set_scale(const FT new_scale) {
				assert(new_scale > FT(0));

				m_imp_scale = new_scale;
				set_automatically_defined_options();
			}

			void estimate_parameters(const bool new_state) {
				m_estimate_parameters = new_state;
			}


			//////////////////
			// Main functions!

			void set_optimal_configuration() {

				set_not_important_options();
				set_more_important_options();
				set_the_most_important_options();
			}

			void set_not_important_options() {
				
				m_pipeline_version = Pipeline_version::WITHOUT_SHAPE_DETECTION; // for the moment the 3D shape based version gives worse results, so do not use it

				m_prefix_path = "stub"; // never used
				m_add_cdt_bbox = false; // never used
				
				m_clean_projected_points = false; // not a good algorithm, can be turned off
				m_preprocessor_scale 	 = 2.0;   // used in cleaning above, useless				
				
				m_clutter_new_point_type = Clutter_new_point_type::CLOSEST; // this is the only method that keeps original points untouched
				m_clutter_fitter_type 	 = Clutter_fitter_type::LINE; 		// only LINE, other types are not implemented
				m_clutter_knn 			 = 12; // never used, since we use CIRCLE neighbour search below

				m_thinning_type 	  			 = Thinning_type::NAIVE; 		  // this is the only one that is currently fully implemented
				m_thinning_neighbour_search_type = Neighbour_search_type::CIRCLE; // in practice, this is the best one, no need to choose any other one

				m_regularizer_reject_planes = true; // in general, rejecting gives more plausible result, used only in the version with 3D shape detection
				m_max_reg_angle          	= 10.0; // in general, 10 is enough, used only in the version with 3D shape detection

				m_use_alpha_shapes 	  = true; // this is the only way to get missing walls, so it is necessary
				m_with_region_growing = true; // this is the only way to apply structuring afterwards and get correct buildings, so it is necessary

				m_structuring_resample 		   = true;  // always resample, visually it is better
				m_structuring_log 	   		   = false; // debug info
				m_structuring_get_all_points   = true;  // in general, better to use all resampled points, since it gives more freedom for the graph cut in CDT
				m_structuring_corner_algorithm = Structuring_corner_algorithm::GRAPH_BASED; // this is the only one that should be used

				m_visibility_small_edge_threshold = -1000000.0; // never used, since we do not use this argument anymore (see rayshooting)
				m_visibility_num_neighbours       = 6;  		// never used, since it is used in barycentric visibility
				m_visibility_rays_per_side        = 10; 		// used in ray shooting visibility, 10 is enough for all cases
				m_visibility_norm_threshold       = 1000.0;     // used in the face based visibility to verify normal of the natural neighbours, 1000 is always enough
				m_visibility_show_progress        = true;       // shows the percentage
				m_visibility_save_info 			  = false; 		// debug info
				m_visibility_sampler  			  = Visibility_sampler::UNIFORM_SUBDIVISION; // this is the only one to use, other samplers create randomness

				m_graph_cut_alpha 	  = 1.0;   // soft parameter
				m_graph_cut_save_info = false; // debug info

				m_building_boundaries_max_inner_iters    = 1000;  	// simple stopping criteria if smth goes wrong
				m_building_boundaries_max_outer_iters    = 1000000; // simple stopping criteria if smth goes wrong
				m_building_boundaries_save_internal_info = false;   // debug info
				m_use_boundaries 		 				 = false; 	// when using UNORIENTED below, this value does not make any difference, but false is preferable
				m_building_boundary_type 				 = Building_boundary_type::UNORIENTED; // this is the most robust method, works with any corrupted data

				m_roof_fitter_type 						  = Roof_fitter_type::AVG; 					 // gives the best visual result, the most robust to the corrupted data
				m_region_growing_normal_estimation_method = Region_growing_normal_estimation::LOCAL; // in general, both work well
				
				m_use_grid_simplifier_first = true; // better to use it, to make the code faster, but if removing it we have one parameter less: m_clutter_cell_length

				m_line_regularizer_optimize_angles    = true; // do we use parallelism and orthogonality in the segment regularization
				m_line_regularizer_optimize_ordinates = true; // do we use collinearity in the segment regularization

				m_polygonizer_number_of_intersections 	   = 2;    // how far should we propagate segments
				m_building_splitter_use_custom_constraints = true; // should we use constraints from the initially detected segments when splitting buildings
			}

			void set_more_important_options() {

				m_structuring_adjacency_method  = Structuring_adjacency_threshold_method::LOCAL; // global is, in general, better, if using local, we have one parameter less: m_str_adj_value
				m_structuring_global_everywhere = false; // better to have false, since in this case, I use global adjacency graph and global corner insertion consistently
				m_structuring_adjacency_value   = 5.0;   // closest distance between two segments for adjacency graph, probably can be removed

				m_visibility_num_samples = 3;     // number of subdivision steps when sampling triangles, 1 or 2 is enough
				m_add_cdt_clutter 		 = false; // better to avoid clutter since it will pollute the final CDT

				m_regularize_lines = false; // regularize lines after region growing
				m_polygonize 	   = false; // prolongate all segments and get a set of polygons.

				m_visibility_approach  = Visibility_approach::FACE_BASED; 				   // face based is, in general, a better but slower option
				m_visibility_method    = Visibility_method::FACE_BASED_NATURAL_NEIGHBOURS; // face based is, in general, a better but slower option
				m_visibility_angle_eps = 0.18; // do not use this ad-hoc, but when using the value about 0.15 - 0.20 is enough

				m_line_regularizer_max_angle_in_degrees = FT(25); // max regularization angle
			}

			void set_the_most_important_options() {

				// Important.
				m_imp_eps   = 3.2; // global distance to the optimal line, (meters)
				m_imp_scale = 5.0; // global distance between adjacent points, (meters)

				m_graph_cut_beta  = 100000.0; // controls how many red and green triangle we will get, (magic)
				m_graph_cut_gamma = 10000.0;  // controls if we should keep constraints satisfied or not, it is the penalty, (magic)


				// Less important.
				m_region_growing_normal_threshold_2d = 0.7; // normal deviation between the point normal and the normal of the optimal line, necessary, (cosine)
				m_region_growing_min_points_2d 	     = 10;  // minimum number of points in the line, probably can be removed, but it will create some noise, (points)

				m_region_growing_normal_threshold_3d = m_region_growing_normal_threshold_2d; // 3d region growing parameters
				m_region_growing_min_points_3d 	     = m_region_growing_min_points_2d;		 // 3d region growing parameters

				m_roof_cleaner_max_percentage = FT(80); // percentage of points that we keep when filtering 3D roof regions after region growing


				// Automatically defined.
				set_automatically_defined_options();
			}

			void set_automatically_defined_options() {

				m_alpha_shape_size 	  			    = m_imp_scale; 	   // does not change often, size in meters to get the boundary of the set of points, necessary, (meters)
				m_structuring_epsilon 			    = m_imp_scale; 	   // distance between adjacent points in the resampled line, (meters)
				m_region_growing_epsilon_2d 		= m_imp_eps; 		   // distance to the optimal line, necessary, (meters)
				m_region_growing_cluster_epsilon_2d = 0.58 * m_imp_scale; // distance between adjacent points, necessary, (meters)
				m_clutter_cell_length 			    = 0.26 * m_imp_scale; // used in the grid simplify, probably can be removed, (meters)

				m_thinning_fuzzy_radius   = m_imp_scale; 	 // radius of the region of points that should be thinned
				m_clutter_filtering_scale = m_imp_scale; 	 // radius of the region of points that should be considered for filtering
				m_clutter_filtering_mean  = m_imp_eps / 5.0; // value of the required mean that should be satisfied by the chosen points in the filtering

				m_polygonizer_min_face_width 			  = m_imp_scale / FT(2); 		  // width of the face that is assumed to be thin and hence is merged with other faces
				m_building_splitter_constraints_threshold = m_polygonizer_min_face_width; // min distance between two segments that are assumed to be the same

				m_line_regularizer_max_difference_in_meters = m_imp_scale / FT(4); // max collinearity difference

				m_region_growing_epsilon_3d 		= m_region_growing_epsilon_2d; 		   // 3d region growing parameters
				m_region_growing_cluster_epsilon_3d = m_region_growing_cluster_epsilon_2d; // 3d region growing parameters

				m_roof_cleaner_scale = m_imp_scale; // all roofs with the size <= than this value are removed
			}

			void set_required_parameters() {
				
				add_str_parameter("-data", m_default_path, m_parameters);
			}

			void set_optional_parameters() {

				// Flags.
				add_bool_parameter("-silent"     , m_silent 			, m_parameters);
				add_bool_parameter("-auto_params", m_estimate_parameters, m_parameters);
				add_bool_parameter("-quality"	 , m_estimate_quality   , m_parameters);
				add_bool_parameter("-clutter"	 , m_add_cdt_clutter    , m_parameters);
				add_bool_parameter("-regularize" , m_regularize_lines   , m_parameters);
				add_bool_parameter("-polygonize" , m_polygonize   		, m_parameters);
				
				bool remove_alpha = false;
				add_bool_parameter("-remove_alpha"   , remove_alpha   		, m_parameters);
				if (remove_alpha && m_use_alpha_shapes) m_use_alpha_shapes = !m_use_alpha_shapes;


				// Important.
				add_val_parameter("-eps"  , m_imp_eps  , m_parameters);
				add_val_parameter("-scale", m_imp_scale, m_parameters);

				set_automatically_defined_options();

				add_val_parameter("-gc_beta" , m_graph_cut_beta , m_parameters);
				add_val_parameter("-gc_gamma", m_graph_cut_gamma, m_parameters);


				// Less important.
				add_val_parameter("-rg_nt_2d" , m_region_growing_normal_threshold_2d, m_parameters);
				add_val_parameter("-rg_min_2d", m_region_growing_min_points_2d      , m_parameters);

				add_val_parameter("-rg_nt_3d" , m_region_growing_normal_threshold_3d, m_parameters);
				add_val_parameter("-rg_min_3d", m_region_growing_min_points_3d      , m_parameters);

				add_val_parameter("-roof_max", m_roof_cleaner_max_percentage, m_parameters);


				// Automatically defined.
				add_val_parameter("-alpha"    , m_alpha_shape_size                 , m_parameters);
				add_val_parameter("-str_eps"  , m_structuring_epsilon              , m_parameters);
				add_val_parameter("-rg_eps_2d", m_region_growing_epsilon_2d 	   , m_parameters);
				add_val_parameter("-rg_ce_2d" , m_region_growing_cluster_epsilon_2d, m_parameters);
				add_val_parameter("-cell"     , m_clutter_cell_length              , m_parameters);

				add_val_parameter("-th_scale", m_thinning_fuzzy_radius  , m_parameters);
				add_val_parameter("-cf_scale", m_clutter_filtering_scale, m_parameters);
				add_val_parameter("-cf_mean" , m_clutter_filtering_mean , m_parameters);

				add_val_parameter("-angle", m_line_regularizer_max_angle_in_degrees, m_parameters);

				add_val_parameter("-rg_eps_3d", m_region_growing_epsilon_3d 	   , m_parameters);
				add_val_parameter("-rg_ce_3d" , m_region_growing_cluster_epsilon_3d, m_parameters);

				add_val_parameter("-roof_scale", m_roof_cleaner_scale, m_parameters);
			}

			void set_user_defined_parameters(const Parameters_wrapper &parameters_wrapper) {

				m_parameters = parameters_wrapper.get();
				std::cout << "Parameters: " << std::endl;

				set_required_parameters();
				set_optional_parameters();
			}

			// All versions.
			void create_lods() {

				assert_global_parameters();
				switch (m_prog_version) {

					case Program_version::VER0:
						create_lods_ver0();
						break;

					default:
						assert(!"Wrong program version!");
						break;
				}
			}

		private:

			template<typename Scalar>
			void add_val_parameter(const std::string &parameter_name, Scalar &variable_value, const Parameters &parameters) {
				
				if (!does_parameter_exist(parameter_name, parameters)) return;
				const std::string parameter_value = parameters.at(parameter_name);

				if (parameter_value != "default")
					variable_value = static_cast<Scalar>(std::stod(parameter_value.c_str()));

				std::cout << parameter_name << " : " << variable_value << std::endl;
			}

			void add_str_parameter(const std::string &parameter_name, std::string &variable_value, const Parameters &parameters) {
				
				if (!does_parameter_exist(parameter_name, parameters)) return;
				const std::string parameter_value = parameters.at(parameter_name);

				if (parameter_value != "default") 
					variable_value = parameter_value;

				std::cout << parameter_name << " : " << variable_value << std::endl;
			}

			void add_bool_parameter(const std::string &parameter_name, bool &variable_value, const Parameters &parameters) {
				
				if (!does_parameter_exist(parameter_name, parameters)) return;

				variable_value = true;
				std::cout << parameter_name << " : " << (variable_value ? "true" : "false") << std::endl;
			}

			bool does_parameter_exist(const std::string &parameter_name, const Parameters &parameters) {
				
				for (typename Parameters::const_iterator param = parameters.begin(); param != parameters.end(); ++param)
					if ((*param).first == parameter_name) return true;

				return false;
			}


			// Main pipeline functions!
			// ------------------------

			void starting_execution() {
				std::cout << "starting ..." << std::endl;
			}

			void loading_data(Container_3D &input, const size_t exec_step) {

				// Load input data.
				std::cout << "(" << exec_step << ") loading; ";

				assert(m_default_path != "default");
				m_loader.get_data(m_default_path, input);

				std::cout << "number of points: " << input.number_of_points() << ";" << std::endl;
			}

			void estimating_initial_parameters(const Container_3D &input, const size_t exec_step) {

				// Estimate some of the required parameters.
				assert(!m_parameters.empty());
				std::cout << "(" << exec_step << ") estimating initial parameters;" << std::endl;

				Parameters_estimator parameters_estimator = Parameters_estimator(input, m_parameters);
				parameters_estimator.estimate();

				std::cout << std::endl << "Updated parameters: " << std::endl;
				set_optional_parameters();
				std::cout << std::endl;
			}

			void applying_selection(
				Indices &ground_idxs,
				Indices &building_boundary_idxs,
				Indices &building_interior_idxs,
				const Container_3D &input, const size_t exec_step) {

				// Split data with respect to different semantic labels.
				std::cout << "(" << exec_step << ") selection; ";

				m_ground_selector.select_elements(input, std::back_inserter(ground_idxs));
				m_building_boundary_selector.select_elements(input, std::back_inserter(building_boundary_idxs));
				m_building_interior_selector.select_elements(input, std::back_inserter(building_interior_idxs));

				if (!m_silent) {
					Log log; 
					log.export_points_using_indices(input, ground_idxs, "tmp" + std::string(PSR) + "lod_0_1" + std::string(PSR) + "ground");
					log.export_points_using_indices(input, building_boundary_idxs, "tmp" + std::string(PSR) + "lod_0_1" + std::string(PSR) + "building_boundary");
					log.export_points_using_indices(input, building_interior_idxs, "tmp" + std::string(PSR) + "lod_0_1" + std::string(PSR) + "building_interior");
				}

				std::cout << "ground: " << ground_idxs.size() << "; boundary: " << building_boundary_idxs.size() << "; interior: " << building_interior_idxs.size() << "; " << std::endl;
			}

			void ground_fitting(
				Plane_3 &base_ground_plane, 
				Plane_3 &fitted_ground_plane,
				Box &fitted_ground_box,
				const Indices &ground_idxs, const Container_3D &input, const size_t exec_step) {

				// Fit plane to the ground points.
				std::cout << "(" << exec_step << ") ground plane fitting;" << std::endl;

				base_ground_plane = Plane_3(FT(0), FT(0), FT(1), FT(0));
				m_utils.fit_ground_plane(input, ground_idxs, fitted_ground_plane);

				using Bounding_boxes = std::vector<Box>;
				Bounding_boxes boxes(1);

				m_utils.return_bounding_box(base_ground_plane, input, ground_idxs, fitted_ground_box);
				boxes[0] = fitted_ground_box;

				if (!m_silent) {
					Log log; log.save_bounding_boxes_as_ply<Bounding_boxes, Point_3>(boxes, "tmp" + std::string(PSR) + "lod_0_1" + std::string(PSR) + "base_ground_plane");
				}

				m_utils.return_bounding_box(fitted_ground_plane, input, ground_idxs, fitted_ground_box);
				boxes[0] = fitted_ground_box;

				if (!m_silent) {
					Log log; log.save_bounding_boxes_as_ply<Bounding_boxes, Point_3>(boxes, "tmp" + std::string(PSR) + "lod_0_1" + std::string(PSR) + "fitted_ground_plane");
				}
			}

			void getting_boundary_points(
				Boundary_data &boundary_clutter,
				const Indices &building_boundary_idxs, const Indices &building_interior_idxs, const Container_3D &input, const size_t exec_step) {

				// Map indices from all detected planes to the ones that are a part of the given facades.
				std::cout << "(" << exec_step << ") getting boundaries; ";

				m_preprocessor.use_alpha_shapes(m_use_alpha_shapes);
				m_preprocessor.set_alpha(m_alpha_shape_size);
				m_preprocessor.make_silent(m_silent);

				const bool stub_state = false;
				Boundary_data stub;

				const auto number_of_boundary_points = m_preprocessor.get_boundary_points(input, building_boundary_idxs, building_interior_idxs, stub_state, stub, boundary_clutter);
				std::cout << std::endl;

				if (!m_silent) {

					Log log;
					log.export_clutter_as_xyz("tmp" + std::string(PSR) + "lod_0_1" + std::string(PSR) + "full_boundary_points", boundary_clutter, input);
				}
			}

			void projecting(
				Projected_points &boundary_clutter_projected,
				const Plane_3 &base_ground_plane, const Boundary_data &boundary_clutter, const Container_3D &input, const size_t exec_step) {

				// Project all points in 2D.
				std::cout << "(" << exec_step << ") projecting; ";

				const int number_of_projected_points = m_ground_projector.project(input, boundary_clutter, base_ground_plane, boundary_clutter_projected);
				std::cout << "clutter projected: " << number_of_projected_points << ";" << std::endl;

				Log points_exporter;
				if (!m_silent && !boundary_clutter_projected.empty())
					points_exporter.export_projected_points_as_xyz("tmp" + std::string(PSR) + "lod_0_1" + std::string(PSR) + "projected_clutter", boundary_clutter_projected, m_default_path);
			}

			void applying_grid_simplification(Projected_points &boundary_clutter_projected, const size_t exec_step) {

				// Remove all unnecessary points.
				assert(m_clutter_new_point_type == Grid_new_point_type::CLOSEST);
				std::cout << "(" << exec_step << ") applying grid simplification; ";

				m_grid_simplifier.set_grid_cell_length(m_clutter_cell_length);
				m_grid_simplifier.set_new_point_type(m_clutter_new_point_type);
				m_grid_simplifier.make_silent(m_silent);

				Boundary_data stub;
				const auto number_of_removed_points = m_grid_simplifier.process(stub, boundary_clutter_projected);

				std::cout << "removed points: " << number_of_removed_points << ";" << std::endl;
			}

			void detecting_2d_lines(
				Boundary_data &boundary_clutter   , Projected_points &boundary_clutter_projected,
				Boundary_data &building_boundaries, Projected_points &building_boundaries_projected,
				const Container_3D &input, const size_t exec_step) {

				// Detect lines in 2D using region growing.
				std::cout << "(" << exec_step << ") detecting 2d lines; ";

				m_region_growing_2.set_epsilon(m_region_growing_epsilon_2d);
				m_region_growing_2.set_cluster_epsilon(m_region_growing_cluster_epsilon_2d);
				m_region_growing_2.set_normal_threshold(m_region_growing_normal_threshold_2d);
				m_region_growing_2.set_minimum_shape_points(m_region_growing_min_points_2d);

				m_region_growing_2.make_silent(m_silent);
				m_region_growing_2.set_normal_estimation_method(m_region_growing_normal_estimation_method);

				const auto number_of_detected_lines = m_region_growing_2.detect(boundary_clutter, boundary_clutter_projected, building_boundaries, building_boundaries_projected, input);
				std::cout << "detected lines: " << number_of_detected_lines << ";" << std::endl;
			}

			void line_fitting(Lines &lines, const Boundary_data &building_boundaries, const Projected_points &building_boundaries_projected, const size_t exec_step) {

				// Fit lines to all detected linear boundary points.
				std::cout << "(" << exec_step << ") line fitting; ";
				const auto number_of_fitted_lines = m_utils.fit_lines_to_projected_points(building_boundaries_projected, building_boundaries, lines);
				std::cout << "number of fitted lines: " << number_of_fitted_lines << ";" << std::endl;

				// if (!m_silent) {
				//    Log lines_exporter; lines_exporter.export_lines_as_obj("tmp" + std::string(PSR) + "fitted_lines", lines, m_default_path);
				// }
			}

			void creating_segments(
				Segments &segments,
				const Lines &lines, const Boundary_data &building_boundaries, const Projected_points &building_boundaries_projected, const std::string &name, const size_t exec_step) {

				// Find segments from the given lines.
				std::cout << "(" << exec_step << ") creating segments; ";
				const auto number_of_segments = m_utils.create_segments_from_lines(building_boundaries_projected, building_boundaries, lines, segments);
				std::cout << "number of segments: " << number_of_segments << ";" << std::endl;
				
				if (!m_silent) {
					Log segments_exporter; segments_exporter.export_segments_as_obj("tmp" + std::string(PSR) + "lod_0_1" + std::string(PSR) + name, segments, m_default_path);
				}
			}

			void regularizing_lines(Segments &segments, Lines &lines, const size_t exec_step) {

				// Regularize lines.
				std::cout << "(" << exec_step << ") regularizing lines; " << std::endl;

				// New regularizer.
				Segment_map segment_map;
				m_line_regularizer.make_silent(m_silent);

				m_line_regularizer.optimize_angles(m_line_regularizer_optimize_angles);
				m_line_regularizer.optimize_ordinates(m_line_regularizer_optimize_ordinates);
				
				m_line_regularizer.set_max_angle_in_degrees(m_line_regularizer_max_angle_in_degrees);
				m_line_regularizer.set_max_difference_in_meters(m_line_regularizer_max_difference_in_meters);

				m_line_regularizer.regularize(segments, segment_map);
				m_line_regularizer.get_lines_from_segments(segments, lines);

				if (!m_silent) {
				    Log segments_exporter; segments_exporter.export_segments_as_obj("tmp" + std::string(PSR) + "lod_0_1" + std::string(PSR) + "regularized_segments", segments, m_default_path);
				}
			}

			void polygonizing(Segments &segments, Data_structure &data_structure, const size_t exec_step) {

				// Prolongate segments and get a set of polygons.
				std::cout << "(" << exec_step << ") polygonizing; " << std::endl;

				Polygonizer polygonizer;
				polygonizer.make_silent(m_silent);

				polygonizer.set_number_of_intersections(m_polygonizer_number_of_intersections);
				polygonizer.set_min_face_width(m_polygonizer_min_face_width);

				polygonizer.polygonize(segments, data_structure);
			}

			void updating_constraints(Data_structure &data_structure, const Segments &segments, const size_t exec_step) {

				// Update constraints in the final partition.
				std::cout << "(" << exec_step << ") updating constraints;" << std::endl;
				
				m_utils.update_constraints(data_structure, segments, m_silent);
			}

			void computing_polygon_based_visibility(Data_structure &data_structure, const Container_3D &input, const size_t exec_step) {

				// Compute visibility for each polygon in the given data structure.
				std::cout << "(" << exec_step << ") computing visibility; " << std::endl;				
				Polygon_based_visibility polygon_based_visibility(input, data_structure);

				polygon_based_visibility.make_silent(m_silent);
				polygon_based_visibility.compute();
			}

			void creating_polygon_based_cdt(CDT &cdt, const Data_structure &data_structure, const Projected_points &boundary_clutter_projected, const size_t exec_step) {

				// Compute constrained Delaunay triangulation from the given data structure.
				std::cout << "(" << exec_step << ") creating cdt;" << std::endl;
				
				m_utils.compute_cdt(cdt, data_structure, boundary_clutter_projected, m_add_cdt_clutter, m_silent);
			}

			void applying_2d_structuring(const Lines &lines, const Boundary_data &building_boundaries, const Projected_points &building_boundaries_projected, const size_t exec_step) {

				// Apply 2D structuring algorithm.
				std::cout << "(" << exec_step << ") 2d structuring; ";

				m_structuring = std::make_shared<Structuring_2>(building_boundaries_projected, building_boundaries, lines);
				
				m_structuring->set_epsilon(m_structuring_epsilon);
				m_structuring->save_log(m_structuring_log);
				m_structuring->resample(m_structuring_resample);
				m_structuring->set_corner_algorithm(m_structuring_corner_algorithm);
				m_structuring->set_adjacency_threshold_method(m_structuring_adjacency_method);
				m_structuring->set_adjacency_threshold(m_structuring_adjacency_value);
				m_structuring->set_global_everywhere(m_structuring_global_everywhere);
				m_structuring->make_silent(m_silent);

				const auto number_of_structured_segments = m_structuring->structure_point_set();
				std::cout << "number of structured segments: " << number_of_structured_segments << ";" << std::endl;
			}

			void applying_clutter_thinning(Boundary_data &boundary_clutter, Projected_points &boundary_clutter_projected, const Container_3D &input, const size_t exec_step) {

				// Apply thinning.
				std::cout << "(" << exec_step << ") applying thinning to clutter; ";

				m_thinning.set_thinning_type(m_thinning_type);
				m_thinning.set_neighbour_search_type(m_thinning_neighbour_search_type);
				m_thinning.set_fuzzy_radius(m_thinning_fuzzy_radius);
				m_thinning.make_silent(m_silent);

				const auto number_of_thinned_points = m_thinning.process(boundary_clutter, boundary_clutter_projected, input);
				std::cout << "thinned points: " << number_of_thinned_points << ";" << std::endl;
			}

			void filtering_clutter(Boundary_data &boundary_clutter, Projected_points &boundary_clutter_projected, const size_t exec_step) {

				// Filter clutter.
				std::cout << "(" << exec_step << ") clutter filtering; ";

				m_clutter_filtering.set_scale(m_clutter_filtering_scale);
				m_clutter_filtering.set_mean(m_clutter_filtering_mean);
				m_clutter_filtering.make_silent(m_silent);

				std::cout << "number of input points: " << boundary_clutter_projected.size();
				const auto number_of_clutter_points = m_clutter_filtering.filter(boundary_clutter, boundary_clutter_projected);
				std::cout << "; number of output points: " << number_of_clutter_points << std::endl;
			}

			void creating_cdt(CDT &cdt, const Boundary_data &boundary_clutter, const Projected_points &boundary_clutter_projected, const Segments &segments, const size_t exec_step) {

				// Compute constrained Delaunay triangulation of the structured points.
				std::cout << "(" << exec_step << ") creating cdt;" << std::endl;
				auto number_of_faces = -1;
				
				// Version with polygonization.
				if (m_polygonize) {

					number_of_faces = m_utils.compute_cdt(cdt, segments, m_silent);
					assert(number_of_faces != -1);

					return;
				}

				// Version without polygonization.
				if (m_structuring != nullptr && !m_structuring->is_empty()) {
					
					const Structured_points   &structured_points = m_structuring_get_all_points ?  m_structuring->get_structured_points()  :  m_structuring->get_segment_end_points();
					const Structured_labels   &structured_labels = m_structuring_get_all_points ?  m_structuring->get_structured_labels()  :  m_structuring->get_segment_end_labels();
					const Structured_anchors &structured_anchors = m_structuring_get_all_points ?  m_structuring->get_structured_anchors() :  m_structuring->get_segment_end_anchors();

					if (!m_structuring_global_everywhere) m_structuring_adjacency_value = m_structuring->get_local_adjacency_value();
					number_of_faces = m_utils.compute_cdt(structured_points, structured_labels, structured_anchors, m_structuring_adjacency_value, cdt, 
														  m_add_cdt_clutter, boundary_clutter, boundary_clutter_projected, m_silent);
				} else {

					number_of_faces = m_utils.compute_cdt(Structured_points(), Structured_labels(), Structured_anchors(), m_structuring_adjacency_value, cdt, 
													 	  m_add_cdt_clutter, boundary_clutter, boundary_clutter_projected, m_silent);
				}
				assert(number_of_faces != -1);
			}

			void converting_3d_to_2d(Container_2D &input_2d, Face_points_map &fp_map, const CDT &cdt, const Container_3D &input, const size_t exec_step) {

				// Convert 3D input to 2D input.
				std::cout << "(" << exec_step << ") converting 3d input into 2d input and setting face to points map; ";
				const auto number_of_converted_points = m_utils.get_2d_input_and_face_points_map(cdt, input, input_2d, fp_map, m_silent);
				std::cout << "number of converted points: " << number_of_converted_points << ";" << std::endl;
			}

			void computing_visibility(CDT &cdt, const Container_2D &input_2d, const size_t exec_step) {

				// Compute visibility (0 - outside or 1 - inside) for each triangle in CDT above.
				if (m_visibility.name() == "ray shooting" && m_pipeline_version == Pipeline_version::WITHOUT_SHAPE_DETECTION && !m_with_region_growing) assert(!"Ray shooting requires constrained edges!");
				if (m_visibility.name() == "blend") assert(!"Blend visibility is not worth trying!");

				std::cout << "(" << exec_step << ") computing visibility; " << std::endl;

				m_visibility.save_info(m_visibility_save_info);
				m_visibility.set_approach(m_visibility_approach);
				m_visibility.set_method(m_visibility_method);
				m_visibility.set_number_of_samples(m_visibility_num_samples);
				m_visibility.show_progress(m_visibility_show_progress);
				m_visibility.set_norm_threshold(m_visibility_norm_threshold);
				m_visibility.set_number_of_neighbours(m_visibility_num_neighbours);
				m_visibility.set_sampler_type(m_visibility_sampler);
				m_visibility.set_number_of_rays_per_side(m_visibility_rays_per_side);
				m_visibility.set_small_edge_threshold(m_visibility_small_edge_threshold);
				m_visibility.set_angle_eps(m_visibility_angle_eps);

				const auto number_of_traversed_faces = m_visibility.compute(input_2d, cdt);
				std::cout << "number of traversed faces: " << number_of_traversed_faces << ";" << std::endl << std::endl;
				
				if (!m_silent) {
					Log eps_saver; eps_saver.save_visibility_eps(cdt);
					Log ply_vis_saver; ply_vis_saver.save_cdt_ply(cdt, "tmp" + std::string(PSR) + "visibility", "in");
				}
			}

			void applying_graph_cut(CDT &cdt, const size_t exec_step) {

				// Apply graph cut.
				std::cout << "(" << exec_step << ") applying graph cut;" << std::endl;

				m_graph_cut.save_info(m_graph_cut_save_info);
				m_graph_cut.set_alpha_parameter(m_graph_cut_alpha);
				m_graph_cut.set_beta_parameter(m_graph_cut_beta);
				m_graph_cut.set_gamma_parameter(m_graph_cut_gamma);
				m_graph_cut.make_silent(m_silent);

				m_graph_cut.max_flow(cdt);

				if (!m_silent) {
					Log ply_cdt_in; ply_cdt_in.save_cdt_ply(cdt, "tmp" + std::string(PSR) + "after_cut", "in");
				}
			}

			void splitting_buildings(Buildings &buildings, CDT &cdt, const Segments &segments, const size_t exec_step) {

				// Split all buildings.
				std::cout << "(" << exec_step << ") splitting buildings; ";

				m_building_splitter.make_silent(m_silent);
				m_building_splitter.use_custom_constraints(m_building_splitter_use_custom_constraints);
				m_building_splitter.set_constraints_threshold(m_building_splitter_constraints_threshold);
				
				const auto number_of_buildings = m_building_splitter.split(cdt, buildings, segments);
				std::cout << "number of buildings: " << number_of_buildings << ";" << std::endl;
			}

			void finding_buildings_boundaries(Buildings &buildings, const CDT &cdt, const size_t exec_step) {

				// Find building's walls.
				std::cout << "(" << exec_step << ") finding boundaries;" << std::endl;

				m_building_outliner.save_info(m_building_boundaries_save_internal_info);
				m_building_outliner.set_max_inner_iterations(m_building_boundaries_max_inner_iters);
				m_building_outliner.set_max_outer_iterations(m_building_boundaries_max_outer_iters);
				m_building_outliner.set_boundary_type(m_building_boundary_type);
					
				m_building_outliner.find_boundaries(cdt, buildings);
			}

			void fitting_roofs(Buildings &buildings, const Plane_3 &fitted_ground_plane, const Face_points_map &fp_map, const Container_3D &input, const CDT &cdt, const size_t exec_step) {

				// Fit roofs for all buildings.
				std::cout << "(" << exec_step << ") fitting roofs;" << std::endl;
				switch (m_roof_fitter_type) {
					
					case Roof_fitter_type::MIN: 
						m_building_min_roof_fitter.fit_roof_heights(cdt, input, fp_map, fitted_ground_plane, buildings);
						break;

					case Roof_fitter_type::AVG: 
						m_building_avg_roof_fitter.fit_roof_heights(cdt, input, fp_map, fitted_ground_plane, buildings);
						break;

					case Roof_fitter_type::MAX: 
						m_building_max_roof_fitter.fit_roof_heights(cdt, input, fp_map, fitted_ground_plane, buildings);
						break;			

					default:
						assert(!"Wrong roof fitter type!");
						break;
				}
			}

			void reconstructing_lod0(Ground &ground_bbox, const CDT &cdt, Buildings &buildings, Mesh &mesh_0, Mesh_facet_colors &mesh_facet_colors_0, const Container_3D &input, const size_t exec_step) {

				// LOD0 reconstruction.
				m_lods.use_boundaries(m_use_boundaries);
				m_utils. template compute_ground_bbox<Ground, Ground_point>(input, ground_bbox);

				assert(!ground_bbox.empty());

				std::cout << "(" << exec_step << ") reconstructing lod0;" << std::endl;
				m_lods.reconstruct_lod0(cdt, buildings, ground_bbox, mesh_0, mesh_facet_colors_0);

				Log lod_0_saver;
				lod_0_saver.save_mesh_as_ply(mesh_0, mesh_facet_colors_0, "LOD0");
			}

			void reconstructing_lod1(const CDT &cdt, Buildings &buildings, Mesh &mesh_1, Mesh_facet_colors &mesh_facet_colors_1, const Ground &ground_bbox, const size_t exec_step) {

				// LOD1 reconstruction.
				std::cout << "(" << exec_step << ") reconstructing lod1;" << std::endl;
				m_lods.reconstruct_lod1(cdt, buildings, ground_bbox, mesh_1, mesh_facet_colors_1);

				Log lod_1_saver; 
				lod_1_saver.save_mesh_as_ply(mesh_1, mesh_facet_colors_1, "LOD1");
			}

			void estimating_lod1_quality(const Container_3D &input, const size_t exec_step) {

				// LOD1 quality estimation.
				std::cout << "(" << exec_step << ") estimating quality of lod1;" << std::endl;

				std::cout << std::endl << "quality statistics: " << std::endl;

				estimate_lod1_complexity(input);
				estimate_lod1_distortion(input);
				estimate_lod1_coverage(input);
			}

			void translating_lod0_and_lod1(
			const FT ground_height, 
			Mesh &mesh_0, const Mesh_facet_colors &mesh_facet_colors_0, 
			Mesh &mesh_1, const Mesh_facet_colors &mesh_facet_colors_1, 
			const size_t exec_step) {

				// Translate LOD0 and LOD1 results to the ground plane.
				std::cout << "(" << exec_step << ") translating lod0 and lod1;" << std::endl;

				m_utils.translate_mesh_along_z(mesh_0, ground_height);
				m_utils.translate_mesh_along_z(mesh_1, ground_height);

				Log lod_saver;
				lod_saver.save_mesh_as_ply(mesh_0, mesh_facet_colors_0, "LOD0");
				lod_saver.save_mesh_as_ply(mesh_1, mesh_facet_colors_1, "LOD1");
			}

			void adding_points_inside_buildings(const Container_3D &input, const CDT &cdt, const Indices &indices, Buildings &buildings, const std::string &name, const size_t exec_step) {

				// Select all points that lie inside each building.
				std::cout << "(" << exec_step << ") adding " << name << " points inside buildings;" << std::endl;

				m_inside_buildings_selector = std::make_shared<Inside_buildings_selector>(input, cdt, indices);
				
				m_inside_buildings_selector->make_silent(m_silent);
				m_inside_buildings_selector->add_indices(buildings);
			}

			void applying_3d_region_growing(const Container_3D &input, Buildings &buildings, const size_t exec_step) {
				
				// Run region growing on all points that form roofs of the given buildings.
				std::cout << "(" << exec_step << ") detecting 3D planes on roof points;" << std::endl;

				m_region_growing_3 = std::make_shared<Region_growing_3>(input);
				m_region_growing_3->make_silent(m_silent);

				m_region_growing_3->set_epsilon(m_region_growing_epsilon_3d);
				m_region_growing_3->set_cluster_epsilon(m_region_growing_cluster_epsilon_3d);
				m_region_growing_3->set_normal_threshold(m_region_growing_normal_threshold_3d);
				m_region_growing_3->set_minimum_shape_points(m_region_growing_min_points_3d);

				m_region_growing_3->detect(buildings);
			}

			void roofs_filtering(const Container_3D &input, const FT ground_height, Buildings &buildings, const size_t exec_step) {

				// Remove all regions detected before that do not satisfy the correct criteria.
				std::cout << "(" << exec_step << ") filtering roof regions;" << std::endl;
				m_roof_cleaner = std::make_shared<Roof_cleaner>(input, ground_height);
				
				m_roof_cleaner->set_scale_upper_bound(m_roof_cleaner_scale);
				m_roof_cleaner->set_max_percentage(m_roof_cleaner_max_percentage);
				
				m_roof_cleaner->make_silent(m_silent);
				m_roof_cleaner->clean_shapes(buildings);
			}

			void creating_partition_input(const Container_3D &input, Buildings &buildings, const size_t exec_step) {
				
				// Fit a plane to each found region of the given roof and compute its bounding box.
				std::cout << "(" << exec_step << ") creating partition input;" << std::endl;

				m_partition_input = std::make_shared<Partition_input>(input);
				m_partition_input->set_alpha(m_alpha_shape_size);

				m_partition_input->create(buildings);
				if (!m_silent) {

					Log exporter; 
					exporter.save_building_roofs_without_faces(buildings, "tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "fitted_roof_planes", true);
					exporter.save_partition_input(buildings, "tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "partition_input");
				}
			}

			void applying_partitioning(const FT ground_height, Buildings &buildings, const size_t exec_step) {
				
				// Apply 3D partitioning and get a set of filtered 2D faces.
				std::cout << "(" << exec_step << ") applying 3D partitioning;" << std::endl;

				m_partition_creator = std::make_shared<Partition_creator>(ground_height);
				m_partition_creator->set_min_face_width(m_polygonizer_min_face_width);

				m_partition_creator->create(buildings);
				if (!m_silent) {
					Log exporter; exporter.save_partition_diagram<Buildings, FT, Point_3>(buildings, ground_height, "tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "lifted_partition_diagram", true);
				}
			}

			void estimating_initial_roofs(const FT ground_height, Buildings &buildings, const size_t exec_step) {

				// Estimate initial roofs by lifting up the 2D diagram obtained from the 3D envelope.
				std::cout << "(" << exec_step << ") estimating initial roofs;" << std::endl;

				m_roof_estimator = std::make_shared<Roof_estimator>(ground_height);
				m_roof_estimator->estimate(buildings);
				
				if (!m_silent) {
					Log exporter; exporter.save_building_roofs_without_faces(buildings, "tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "estimated_roofs", true);
				}
			}

			void reconstructing_lod2(const Buildings &buildings, const Ground &ground_bbox, const FT ground_height, Mesh &mesh, Mesh_facet_colors &mesh_facet_colors, const size_t exec_step) {

				// LOD2 reconstruction.
				std::cout << "(" << exec_step << ") reconstructing lod2;" << std::endl;
				
				m_lod2 = std::make_shared<LOD2_reconstruction>(buildings, ground_bbox, ground_height, mesh_facet_colors);
				m_lod2->reconstruct(mesh);

				Log lod2_saver; 
				lod2_saver.save_mesh_as_ply(mesh, mesh_facet_colors, "LOD2");
			}

			void finishing_execution() {
				std::cout << "... finishing" + std::string(PN) + "" << std::endl;
			}

		public:

			// Quality metrics.
			void estimate_lod1_complexity(const Container_3D &input) {
				
				m_lod_complexity = std::make_shared<Lod_complexity>(input, m_lods);
				m_lod_complexity->estimate();
				m_complexity = m_lod_complexity->get();

				std::cout << "complexity = " << m_complexity << " elements." << std::endl;
			}

			void estimate_lod1_distortion(const Container_3D &input) {
				
				m_lod_distortion = std::make_shared<Lod_distortion>(input, m_lods);
				m_lod_distortion->estimate();
				m_distortion = m_lod_distortion->get();

				std::cout << "distortion = " << m_distortion << " meters." << std::endl;
			}

			void estimate_lod1_coverage(const Container_3D &input) {
				
				m_lod_coverage = std::make_shared<Lod_coverage>(input, m_lods);
				m_lod_coverage->estimate(m_imp_eps);
				m_coverage = m_lod_coverage->get();

				std::cout << "coverage = " << m_coverage << " percents." << std::endl << std::endl;
			}


			// Version 0.
			void create_lods_ver0() {

				CGAL::Timer timer;
				
				timer.start();
				run_pipeline_ver0();
				timer.stop();

				std::cout << "Running time: " << timer.time() << " seconds." << std::endl << std::endl << std::endl;
			}


			void creating_lod0_and_lod1(const Container_3D &input, 
			Indices &building_interior_idxs, Box &fitted_ground_box, Ground &ground_bbox, CDT &cdt, Buildings &buildings,
			Mesh &mesh_0, 
			Mesh &mesh_1, 
			Mesh_facet_colors &mesh_facet_colors_0, 
			Mesh_facet_colors &mesh_facet_colors_1) {

				size_t exec_step = 0;
				std::cout << std::endl << "...CREATING LOD0 AND LOD1..." << std::endl << std::endl;
				

				// (01) ----------------------------------
				Indices ground_idxs, building_boundary_idxs;
				applying_selection(ground_idxs, building_boundary_idxs, building_interior_idxs, input, ++exec_step);


				// (02) ----------------------------------
				Plane_3 base_ground_plane, fitted_ground_plane;
				ground_fitting(base_ground_plane, fitted_ground_plane, fitted_ground_box, ground_idxs, input, ++exec_step);


				// (03) ----------------------------------
				Boundary_data boundary_clutter;
				getting_boundary_points(boundary_clutter, building_boundary_idxs, building_interior_idxs, input, ++exec_step);


				// (04) ----------------------------------
				Projected_points boundary_clutter_projected;
				projecting(boundary_clutter_projected, base_ground_plane, boundary_clutter, input, ++exec_step);


				// (05) ----------------------------------
				applying_grid_simplification(boundary_clutter_projected, ++exec_step);


				// (06) ----------------------------------
				Boundary_data building_boundaries; Projected_points building_boundaries_projected;
				detecting_2d_lines(boundary_clutter, boundary_clutter_projected, building_boundaries, building_boundaries_projected, input, ++exec_step);


				// (07) ----------------------------------
				Lines lines;
				line_fitting(lines, building_boundaries, building_boundaries_projected, ++exec_step);


				// (08) ----------------------------------
				Segments segments;
				creating_segments(segments, lines, building_boundaries, building_boundaries_projected, "original_segments", ++exec_step);


				// (09) ----------------------------------
				if (m_regularize_lines) regularizing_lines(segments, lines, ++exec_step);


				// (10) ----------------------------------
				Segments original_segments = segments;

				if (m_polygonize) {
					Data_structure data_structure;

					polygonizing(segments, data_structure, ++exec_step);
					computing_polygon_based_visibility(data_structure, input, ++exec_step);

					if (m_add_cdt_clutter) {
						applying_clutter_thinning(boundary_clutter, boundary_clutter_projected, input, ++exec_step);
						filtering_clutter(boundary_clutter, boundary_clutter_projected, ++exec_step);
					}

					creating_polygon_based_cdt(cdt, data_structure, boundary_clutter_projected, ++exec_step);
				}


				// (11) ----------------------------------
				if (!m_polygonize) 
					applying_2d_structuring(lines, building_boundaries, building_boundaries_projected, ++exec_step);


				if (m_add_cdt_clutter && !m_polygonize) {
					
					// (12) ----------------------------------
					applying_clutter_thinning(boundary_clutter, boundary_clutter_projected, input, ++exec_step);

				
					// (13) ----------------------------------
					filtering_clutter(boundary_clutter, boundary_clutter_projected, ++exec_step);
				}


				// (14) ----------------------------------
				if (!m_polygonize) 
					creating_cdt(cdt, boundary_clutter, boundary_clutter_projected, segments, ++exec_step);


				// (15) ----------------------------------
				Container_2D input_2d; Face_points_map fp_map;
				converting_3d_to_2d(input_2d, fp_map, cdt, input, ++exec_step);


				// (16) ----------------------------------
				if (!m_polygonize) 
					computing_visibility(cdt, input_2d, ++exec_step);


				// (17) ----------------------------------
				if (!m_polygonize) 
					applying_graph_cut(cdt, ++exec_step);


				// From now on we handle each building separately.

				// (18) ----------------------------------				
				splitting_buildings(buildings, cdt, original_segments, ++exec_step);


				// (19) ----------------------------------				
				finding_buildings_boundaries(buildings, cdt, ++exec_step);


				// (20) ----------------------------------
				fitting_roofs(buildings, fitted_ground_plane, fp_map, input, cdt, ++exec_step);


				// (21) ----------------------------------
				reconstructing_lod0(ground_bbox, cdt, buildings, mesh_0, mesh_facet_colors_0, input, ++exec_step);


				// (22) ----------------------------------	
				reconstructing_lod1(cdt, buildings, mesh_1, mesh_facet_colors_1, ground_bbox, ++exec_step);
			}


			void creating_lod2(const Container_3D &input, const Indices &building_interior_idxs, const CDT &cdt, const Box &fitted_ground_box, const Ground &ground_bbox,
			Buildings &buildings, 
			Mesh &mesh_0, 
			Mesh &mesh_1, 
			Mesh &mesh_2, 
			const Mesh_facet_colors &mesh_facet_colors_0, 
			const Mesh_facet_colors &mesh_facet_colors_1,
				  Mesh_facet_colors &mesh_facet_colors_2) {
				
				size_t exec_step = 0;
				std::cout << std::endl << "...CREATING LOD2..." << std::endl << std::endl;


				// (01) ----------------------------------
				const FT ground_height = get_average_target_height(fitted_ground_box);
				translating_lod0_and_lod1(ground_height, mesh_0, mesh_facet_colors_0, mesh_1, mesh_facet_colors_1, ++exec_step);


				// (02) ----------------------------------
				clear_interior_indices(buildings);
				adding_points_inside_buildings(input, cdt, building_interior_idxs, buildings, "interior", ++exec_step);


				// (03) ----------------------------------
				applying_3d_region_growing(input, buildings, ++exec_step);
				clear_interior_indices(buildings);


				// (04) ----------------------------------
				roofs_filtering(input, ground_height, buildings, ++exec_step);

				
				// (05) ----------------------------------
				creating_partition_input(input, buildings, ++exec_step);
				clear_shapes(buildings);


				// (06) ----------------------------------
				applying_partitioning(ground_height, buildings, ++exec_step);


				// (07) ----------------------------------
				estimating_initial_roofs(ground_height, buildings, ++exec_step);


				// (08) ----------------------------------
				reconstructing_lod2(buildings, ground_bbox, ground_height, mesh_2, mesh_facet_colors_2, ++exec_step);
			}


			void run_pipeline_ver0() {

				// --- Start ----------->
				size_t exec_step = 0;
				std::cout << std::endl << "...PREPARING DATA..." << std::endl << std::endl;

				// (--) -------------------------------------				
				starting_execution();

				// (initial) --------------------------------
				Container_3D input;
				loading_data(input, ++exec_step);

				// (extra) ----------------------------------
				if (m_estimate_parameters) estimating_initial_parameters(input, ++exec_step);


				// --- Main part ----------->

				// (lod0 and lod1) --------------------------
				Indices building_interior_idxs;
				Box fitted_ground_box;
				
				CDT cdt;
				Buildings buildings;
				Ground ground_bbox;

				Mesh mesh_0, mesh_1, mesh_2; Mesh_facet_colors mesh_facet_colors_0, mesh_facet_colors_1, mesh_facet_colors_2;
				creating_lod0_and_lod1(input, building_interior_idxs, fitted_ground_box, ground_bbox, cdt, buildings, mesh_0, mesh_1, mesh_facet_colors_0, mesh_facet_colors_1);

				// (lod2) -----------------------------------
				creating_lod2(input, building_interior_idxs, cdt, fitted_ground_box, ground_bbox, buildings, mesh_0, mesh_1, mesh_2, mesh_facet_colors_0, mesh_facet_colors_1, mesh_facet_colors_2);


				// --- End ----------->
				exec_step = 0;
				std::cout << std::endl << "...POSTPROCESSING..." << std::endl << std::endl;
			
				// (extra) ----------------------------------
				if (m_estimate_quality) estimating_lod1_quality(input, ++exec_step);

				// (--) ----------------------------------
				finishing_execution();
			}

		private:

			// Main components.
			Loader       m_loader;
			Preprocessor m_preprocessor;
			
			Building_boundary_selector m_building_boundary_selector;
			Building_interior_selector m_building_interior_selector;
			Clutter_selector           m_clutter_selector;
			Ground_selector            m_ground_selector;

			Vertical_regularizer m_vertical_regularizer;
			Ground_projector 	 m_ground_projector;
			Visibility_2 		 m_visibility;
			Utils 		 		 m_utils;
			Region_growing_2   	 m_region_growing_2;
			
			Graph_cut m_graph_cut;
			Lods m_lods;

			Building_splitter m_building_splitter;
			Building_outliner m_building_outliner;
			
			Building_min_roof_fitter m_building_min_roof_fitter;
			Building_avg_roof_fitter m_building_avg_roof_fitter;
			Building_max_roof_fitter m_building_max_roof_fitter;

			std::shared_ptr<Structuring_2> m_structuring;
			Clutter_processor m_clutter_processor;
			
			Grid_simplifier   m_grid_simplifier;
			Thinning 		  m_thinning;
			Clutter_filtering m_clutter_filtering;

			Line_regularizer m_line_regularizer;

			std::shared_ptr<Inside_buildings_selector> m_inside_buildings_selector;
			std::shared_ptr<Region_growing_3> 		   m_region_growing_3;
			std::shared_ptr<Roof_cleaner> 		       m_roof_cleaner;
			std::shared_ptr<Partition_input> 		   m_partition_input;
			std::shared_ptr<Partition_creator> 		   m_partition_creator;
			std::shared_ptr<Roof_estimator> 		   m_roof_estimator;
			
			std::shared_ptr<LOD2_reconstruction> m_lod2;


			// Global parameters.
			std::string m_prefix_path;
			std::string m_default_path;
			
			FT m_preprocessor_scale;
			FT m_structuring_epsilon;
			
			bool m_structuring_log;
			bool m_structuring_resample;
			bool m_structuring_get_all_points;

			bool m_add_cdt_clutter;
			bool m_add_cdt_bbox;

			bool m_visibility_save_info;
			Visibility_approach m_visibility_approach;
			Visibility_method m_visibility_method;
			size_t m_visibility_num_samples;

			bool m_graph_cut_save_info;
			FT m_graph_cut_alpha;
			FT m_graph_cut_beta;
			FT m_graph_cut_gamma;

			bool m_building_boundaries_save_internal_info;
			size_t m_building_boundaries_max_inner_iters;
			size_t m_building_boundaries_max_outer_iters;

			Roof_fitter_type m_roof_fitter_type;
			bool m_clean_projected_points;

			FT m_max_reg_angle;
			bool m_regularizer_reject_planes;
			bool m_use_boundaries;

			const Program_version  m_prog_version;
			
			Pipeline_version m_pipeline_version;

			bool m_visibility_show_progress;
			FT   m_visibility_norm_threshold;

			size_t 				   m_clutter_knn;
			FT 					   m_clutter_cell_length;
			Clutter_fitter_type    m_clutter_fitter_type;
			Clutter_new_point_type m_clutter_new_point_type;

			size_t m_visibility_num_neighbours;
			Visibility_sampler m_visibility_sampler;
			size_t m_visibility_rays_per_side;
			FT m_visibility_small_edge_threshold;

			Building_boundary_type m_building_boundary_type;
			FT m_visibility_angle_eps;

			Neighbour_search_type m_thinning_neighbour_search_type;
			FT 				      m_thinning_fuzzy_radius;
			Thinning_type 		  m_thinning_type;
			
			FT 	   m_region_growing_epsilon_2d;
			FT     m_region_growing_cluster_epsilon_2d;
			FT 	   m_region_growing_normal_threshold_2d;
			size_t m_region_growing_min_points_2d;

			bool m_with_region_growing;
			bool m_use_grid_simplifier_first;

			FT m_alpha_shape_size;
			bool m_use_alpha_shapes;

			Structuring_corner_algorithm 		   m_structuring_corner_algorithm;
			Structuring_adjacency_threshold_method m_structuring_adjacency_method;
			FT 									   m_structuring_adjacency_value;
			
			bool m_structuring_global_everywhere;
			bool m_silent;

			Main_test_data_type m_test_data_type;
			Region_growing_normal_estimation m_region_growing_normal_estimation_method;

			FT m_imp_eps;
			FT m_imp_scale;

			bool m_estimate_parameters;
			bool m_estimate_quality;

			Parameters m_parameters;

			FT m_complexity;
			FT m_distortion;
			FT m_coverage;

			std::shared_ptr<Lod_complexity> m_lod_complexity;
			std::shared_ptr<Lod_distortion> m_lod_distortion;
			std::shared_ptr<Lod_coverage>   m_lod_coverage;

			FT m_clutter_filtering_scale;
			FT m_clutter_filtering_mean;

			bool m_regularize_lines;
			bool m_polygonize;

			bool m_line_regularizer_optimize_angles;
			bool m_line_regularizer_optimize_ordinates;

			size_t m_polygonizer_number_of_intersections;
			FT 	 m_polygonizer_min_face_width;

			bool m_building_splitter_use_custom_constraints;
			FT   m_building_splitter_constraints_threshold;

			FT m_line_regularizer_max_angle_in_degrees;
			FT m_line_regularizer_max_difference_in_meters;

			FT 	   m_region_growing_epsilon_3d;
			FT     m_region_growing_cluster_epsilon_3d;
			FT 	   m_region_growing_normal_threshold_3d;
			size_t m_region_growing_min_points_3d;

			FT m_roof_cleaner_scale;
			FT m_roof_cleaner_max_percentage;


			// Assert default values of all global parameters.
			void assert_global_parameters() {
			
				assert(m_default_path != "default");

				assert(m_preprocessor_scale  != -FT(1));
				assert(m_structuring_epsilon != -FT(1));

				assert(!m_add_cdt_bbox);
				assert(m_visibility_num_samples != 0);

				assert(!(m_visibility_approach == Visibility_approach::FACE_BASED  && m_visibility_method == Visibility_method::POINT_BASED_CLASSIFICATION));
				assert(!(m_visibility_approach == Visibility_approach::POINT_BASED && m_visibility_method == Visibility_method::FACE_BASED_COUNT));
				assert(!(m_visibility_approach == Visibility_approach::POINT_BASED && m_visibility_method == Visibility_method::FACE_BASED_NATURAL_NEIGHBOURS));
				assert(!(m_visibility_approach == Visibility_approach::POINT_BASED && m_visibility_method == Visibility_method::FACE_BASED_BARYCENTRIC));

				assert(m_graph_cut_alpha != -FT(1));
				assert(m_graph_cut_beta  != -FT(1));
				assert(m_graph_cut_gamma != -FT(1));

				assert(m_building_boundaries_max_inner_iters != 0);
				assert(m_building_boundaries_max_outer_iters != 0);

				assert(m_max_reg_angle != -FT(1));
				assert(m_visibility_norm_threshold != -FT(1));
				
				assert(m_clutter_knn > 1);
				assert(m_clutter_cell_length != -FT(1));

				assert(m_visibility_num_neighbours > 1);
				assert(m_visibility_rays_per_side  > 0);
				assert(m_visibility_angle_eps != -FT(1));

				assert(m_thinning_fuzzy_radius != -FT(1));
				assert(m_structuring_adjacency_value > FT(0));

				if (m_with_region_growing) {
					assert(m_region_growing_epsilon_2d 		    != -FT(1));
					assert(m_region_growing_cluster_epsilon_2d  != -FT(1));
					assert(m_region_growing_normal_threshold_2d != -FT(1));
					assert(m_region_growing_min_points_2d 		!=     0);
				}

				assert(m_region_growing_epsilon_3d 		    != -FT(1));
				assert(m_region_growing_cluster_epsilon_3d  != -FT(1));
				assert(m_region_growing_normal_threshold_3d != -FT(1));
				assert(m_region_growing_min_points_3d 		!=     0);

				if (m_use_alpha_shapes) {
					assert(m_alpha_shape_size > FT(0));
				}

				assert(m_prefix_path != "path_to_the_data_folder" && m_prefix_path != "default");

				assert(m_imp_eps   > FT(0));
				assert(m_imp_scale > FT(0));

				assert(m_clutter_filtering_scale > FT(0));
				assert(m_clutter_filtering_mean  > FT(0));

				assert(m_polygonizer_number_of_intersections > 0);
				assert(m_polygonizer_min_face_width          > FT(0));

				assert(m_building_splitter_constraints_threshold > FT(0));

				assert(m_line_regularizer_max_angle_in_degrees     > FT(0));
				assert(m_line_regularizer_max_difference_in_meters > FT(0));

				assert(m_roof_cleaner_scale > FT(0));
				assert(m_roof_cleaner_max_percentage >= FT(0) && m_roof_cleaner_max_percentage <= FT(100));
			}

			void clear_interior_indices(Buildings &buildings) {
				for (typename Buildings::iterator bit = buildings.begin(); bit != buildings.end(); ++bit)
					bit->second.clear_interior_indices();
			}

			void clear_shapes(Buildings &buildings) {
				for (typename Buildings::iterator bit = buildings.begin(); bit != buildings.end(); ++bit)
					bit->second.clear_shapes();
			}

			FT get_average_target_height(const Box &fitted_ground_box) const {
				FT target_height = FT(0);

				target_height += fitted_ground_box.bbmin.z();
				target_height += fitted_ground_box.bbmax.z();

				target_height /= FT(2);
				return target_height;
			}


			// Some extra functions. Can be removed!
			void export_segments(const Segments &segments) {

				const std::string path = "/Users/danisimo/Documents/pipeline/logs/regularizer-data/segments.data";
				std::ofstream file(path.c_str(), std::ios_base::out);

				if (!file) std::cerr << std::endl << "ERROR: Error saving file with segments!" << std::endl << std::endl;

				file << segments.size() << std::endl;
				for (size_t i = 0; i < segments.size(); ++i)
					file << segments[i].source() << " " << segments[i].target() << std::endl;

				file.close();
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_BASE_H