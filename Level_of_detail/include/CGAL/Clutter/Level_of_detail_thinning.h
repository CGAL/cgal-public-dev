#ifndef CGAL_LEVEL_OF_DETAIL_THINNING_H
#define CGAL_LEVEL_OF_DETAIL_THINNING_H

#if defined(WIN32) || defined(_WIN32) 
#define PS "\\" 
#else 
#define PS "/" 
#endif 

// STL includes.
#include <map>
#include <vector>
#include <cassert>
#include <cmath>
#include <array>
#include <tuple>
#include <stdio.h> 
#include <stdlib.h>
#include <time.h>  

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/number_utils.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/property_map.h>
#include <CGAL/Simple_cartesian.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Clutter/Level_of_detail_grid_simplify.h>
#include <CGAL/Utils/Level_of_detail_utils.h>

// Extra libs.
#include <CGAL/Clutter/kmeans/dkm.hpp>
#include <CGAL/Clutter/kmeans/dkm_utils.hpp>

namespace CGAL {

	namespace LOD {

		// Thinning.
		template<class KernelTraits, class BoundaryData, class ProjectedPoints, class InputContainer>
		class Level_of_detail_thinning {

		public:
			// Fields.
			typedef KernelTraits 	Kernel;
			typedef BoundaryData 	Boundary_data;
			typedef ProjectedPoints Projected_points;
			typedef InputContainer  Container;

			typedef typename Kernel::FT 	  FT;
			typedef typename Kernel::Point_2  Point_2;
			typedef typename Kernel::Line_2   Line_2;
			typedef typename Kernel::Vector_2 Normal;

			typedef std::pair<int, Point_2> 					 				Projected_point;
			typedef typename CGAL::Second_of_pair_property_map<Projected_point> Point_map;

			typedef CGAL::Search_traits_2<Kernel>                       					 Search_traits_2;
			typedef CGAL::Search_traits_adapter<Projected_point, Point_map, Search_traits_2> Search_traits;
			typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> 					     Neighbor_search;
			typedef CGAL::Fuzzy_sphere<Search_traits>                    					 Fuzzy_circle;
			typedef CGAL::Kd_tree<Search_traits>					                         Fuzzy_tree;
			typedef CGAL::Fuzzy_iso_box<Search_traits> 										 Fuzzy_square;

			typedef typename Boundary_data::const_iterator    Boundary_iterator;
			typedef typename Projected_points::const_iterator Point_iterator;
			typedef typename Neighbor_search::iterator  	  Neighbour_iterator;

			typedef Level_of_detail_utils_simple<Kernel> Simple_utils;

			typename Kernel::Compute_squared_distance_2 squared_distance;
			typename Kernel::Compute_scalar_product_3   dot_product;


			// Extra.
			using Log = CGAL::LOD::Mylog;

			using Cell_id  = std::pair<int, int>;
			using Grid_map = std::map<Cell_id, std::vector<int> >;

			using Grid_cell_iterator = Grid_map::const_iterator;
			using Normals 			 = std::map<int, Normal>;
			using Corner 			 = std::pair<Point_2, Normal>;
			using Corners			 = std::vector<Corner>;
			using Array_point 		 = std::array<double, 2>;
			using Kmeans_data 		 = std::vector<Array_point>;
			using Kmeans_map 		 = std::vector<int>;
			using Kmeans_result 	 = std::tuple< std::vector< std::array<double, 2> >, std::vector< uint32_t > >;


			// Constructor.
			Level_of_detail_thinning() : 
			m_k(3), 
			m_fitter_type(Thinning_fitter_type::LINE), 
			m_neighbour_search_type(Neighbour_search_type::CIRCLE),
			m_radius(FT(1) / FT(4)),
			m_thinning_type(Thinning_type::NAIVE),
			m_scale_type(Thinning_scale_type::FIXED),
			m_min_knn(2),
			m_dim_0_max(FT(5) / FT(100)),
			m_dim_2_sparsity_num_cells(FT(5)),
			m_dim_2_sparsity_max(FT(8) / FT(10)),
			m_max_num_clusters(5),
			m_silent(false) { }


			// Input functions.
			void set_number_of_neighbours(const size_t new_value) {

				assert(new_value >= 2);
				m_k = new_value;
			}

			void set_fuzzy_radius(const FT new_value) {

				assert(new_value > FT(0));
				m_radius = new_value;
			}

			void set_thinning_type(const Thinning_type new_type) {
				m_thinning_type = new_type;
			}

			void set_neighbour_search_type(const Neighbour_search_type new_type) {
				m_neighbour_search_type = new_type;
			}

			void set_fitter_type(const Thinning_fitter_type new_fitter_type) {
				m_fitter_type = new_fitter_type;
			}

			void set_dim_0_max(const FT new_value) {
				
				assert(new_value > FT(0));
				m_dim_0_max = new_value;
			}

			void set_dim_2_sparsity_num_cells(const FT new_value) {
				
				assert(new_value > FT(0));
				m_dim_2_sparsity_num_cells = new_value;
			}

			void set_dim_2_sparsity_max(const FT new_value) {
				
				assert(new_value > FT(0));
				m_dim_2_sparsity_max = new_value;
			}

			void set_max_number_of_clusters(const size_t new_value) {

				assert(new_value > 0);
				m_max_num_clusters = new_value;
			}

			void make_silent(const bool silent) {
				m_silent = silent;
			}


			// Main.
			int process(Boundary_data &boundary_clutter, Projected_points &boundary_clutter_projected, const Container &input) const {

				switch (m_scale_type) {

					case Thinning_scale_type::FIXED:
						return fixed_process(boundary_clutter, boundary_clutter_projected, input);

					case Thinning_scale_type::ADAPTIVE_FIXED:
						choose_scale_adaptively();
						return fixed_process(boundary_clutter, boundary_clutter_projected, input);

					case Thinning_scale_type::PROGRESSIVE:
						return progressive_process(boundary_clutter, boundary_clutter_projected, input);

					default:
						assert(!"Wrong thinning scale type!");
						return -1;
				}
			}


		// Fields and extra functions.
		private:
			size_t 		   		  m_k;
			Thinning_fitter_type  m_fitter_type;
			Neighbour_search_type m_neighbour_search_type;
			FT 					  m_radius;
			Thinning_type 		  m_thinning_type;
			Thinning_scale_type   m_scale_type;
			
			const size_t m_min_knn;
			
			FT 	   m_dim_0_max; 		       // total percentage
			FT 	   m_dim_2_sparsity_num_cells; // number of cells
			FT 	   m_dim_2_sparsity_max; 	   // total percentage
			size_t m_max_num_clusters; 		   // max number of clusters

			Simple_utils m_simple_utils;
			bool m_silent;

			void choose_scale_adaptively() const {

				// To be implemented later.
				assert(!"Adaptive method is not implemented!");
			}

			int fixed_process(Boundary_data &boundary_clutter, Projected_points &boundary_clutter_projected, const Container &input) const {

				FT stub;
				switch (m_thinning_type) {

					case Thinning_type::NAIVE:
						return apply_naive_thinning(boundary_clutter, boundary_clutter_projected, stub);

					case Thinning_type::COMPLEX:
						assert(!"This method is not finished!");
						return apply_complex_thinning(boundary_clutter, boundary_clutter_projected, stub, input);

					default:
						assert(!"Wrong thinning type!");
						return -1;
				}
			}

			int progressive_process(Boundary_data &, Projected_points &, const Container &) const {

				assert(!"Progressive process is not implemented!");

				// To be implemented later.
				return -1;
			}

			FT compute_scale_error() const {

				// To be implemented later.
				return -FT(1);
			}


		// Naive thinning.
		private:
			int apply_naive_thinning(Boundary_data &, Projected_points &boundary_clutter_projected, FT &error) const {

				int number_of_processed_points = 0;
				if (boundary_clutter_projected.empty()) {

					std::cerr << "WARNING: Thinning empty input!" << std::endl;
					return number_of_processed_points;
				}
				
				Fuzzy_tree tree;
				create_tree(tree, boundary_clutter_projected);

				Projected_points thinned_points;
				thin_all_points(thinned_points, tree, boundary_clutter_projected);

				assert(boundary_clutter_projected.size() == thinned_points.size());
				boundary_clutter_projected = thinned_points;

				number_of_processed_points = static_cast<int>(thinned_points.size());
				assert(number_of_processed_points >= 0);

				if (!m_silent) {
					Log log; log.export_projected_points_as_xyz("tmp" + std::string(PS) + "thinned_clutter", boundary_clutter_projected, "unused path");
				}

				error = compute_scale_error();
				return number_of_processed_points;
			}

			void create_tree(Fuzzy_tree &tree, const Projected_points &boundary_clutter_projected) const {

				assert(!boundary_clutter_projected.empty());
				tree.insert(boundary_clutter_projected.begin(), boundary_clutter_projected.end());
			}

			void thin_all_points(Projected_points &thinned_points, const Fuzzy_tree &tree, const Projected_points &boundary_clutter_projected) const {

				assert(!boundary_clutter_projected.empty());
				thinned_points.clear();

				for (Point_iterator pit = boundary_clutter_projected.begin(); pit != boundary_clutter_projected.end(); ++pit) {
					
					const Projected_point &projected = *pit;
					handle_projected_point_naive(thinned_points, tree, projected);
				}
			}

			void handle_projected_point_naive(Projected_points &thinned_points, const Fuzzy_tree &tree, const Projected_point &query) const {

				Projected_points neighbours;
				find_nearest_neighbours(neighbours, tree, query);

				switch (m_fitter_type) {

					case Thinning_fitter_type::LINE:
						handle_with_line_naive(thinned_points, query, neighbours);
						break;

					default:
						assert(!"Wrong fitter type!");
						break;
				}
			}

			void find_nearest_neighbours(Projected_points &neighbours, const Fuzzy_tree &tree, const Projected_point &query) const {
				
				switch (m_neighbour_search_type) {

					case Neighbour_search_type::KNN:
						find_nearest_neighbours_with_knn(neighbours, tree, query, m_k);
						break;

					case Neighbour_search_type::CIRCLE:
						find_nearest_neighbours_with_circle(neighbours, tree, query);
						break;

					case Neighbour_search_type::SQUARE:
						assert(!"Square neighbour search does not work!");
						find_nearest_neighbours_with_square(neighbours, tree, query);
						break;

					default:
						assert(!"Wrong neighbour search type!");
						break;
				}

				if (neighbours.size() < m_min_knn)
					find_nearest_neighbours_with_knn(neighbours, tree, query, m_min_knn);

				assert(neighbours.size() >= m_min_knn);
			}

			void find_nearest_neighbours_with_knn(Projected_points &neighbours, const Fuzzy_tree &tree, const Projected_point &query, const size_t num_knn) const {

				assert(num_knn >= m_min_knn);
				
				Neighbor_search search(tree, query.second, num_knn);
				const size_t num_points = static_cast<size_t>(std::distance(search.begin(), search.end()));

				neighbours.clear();
				size_t count = 0;

				for (Neighbour_iterator nit = search.begin(); nit != search.end(); ++nit, ++count) neighbours[nit->first.first] = nit->first.second;
				assert(count == num_points);

				if (count != num_points) {
					std::cerr << "Error: count != num_points, find_nearest_neighbours_with_knn function thinning!" << std::endl;
					exit(EXIT_FAILURE);
				}

				neighbours[query.first] = query.second;
			}

			void find_nearest_neighbours_with_circle(Projected_points &neighbours, const Fuzzy_tree &tree, const Projected_point &query) const {

				neighbours.clear();

				Fuzzy_circle circle;
				compute_circle(query.second, circle);

				tree.search(std::inserter(neighbours, neighbours.end()), circle);
				neighbours[query.first] = query.second;
			}

			void find_nearest_neighbours_with_square(Projected_points &neighbours, const Fuzzy_tree &tree, const Projected_point &query) const {

				neighbours.clear();

				Fuzzy_square square;
				compute_square(query.second, square);

				tree.search(std::inserter(neighbours, neighbours.end()), square);
				neighbours[query.first] = query.second;
			}

			void compute_circle(const Point_2 &centre, Fuzzy_circle &circle) const {

				assert(m_radius > FT(0));
				circle = Fuzzy_circle(centre, m_radius);
			}

			void compute_square(const Point_2 &centre, Fuzzy_square &square) const {

				assert(m_radius > FT(0));

				const FT cx = centre.x();
				const FT cy = centre.y();

				const Point_2 p = Point_2(cx - m_radius, cy - m_radius);
				const Point_2 q = Point_2(cx + m_radius, cy + m_radius);

				square = Fuzzy_square(p, q);
			}

			void handle_with_line_naive(Projected_points &thinned_points, const Projected_point &query, const Projected_points &neighbours) const {

				assert(neighbours.size() >= m_min_knn);
				Line_2 line;

				fit_line_to_points(neighbours, line);
				create_thin_point_with_line(thinned_points, query, line);
			}

			void fit_line_to_points(const Projected_points &neighbours, Line_2 &line) const {

				// You can add here a method to fit line based on the centre of mass of the given points + an average normal.

      			using Local_Kernel = CGAL::Simple_cartesian<double>;
				using Point_2ft    = Local_Kernel::Point_2;
				using Line_2ft     = Local_Kernel::Line_2;

				const size_t num_points = neighbours.size();
				std::vector<Point_2ft> tmp_points(num_points);

				size_t count = 0;
				for (Point_iterator pit = neighbours.begin(); pit != neighbours.end(); ++pit, ++count) {
				
					const Point_2 &p = (*pit).second;

					const double x = CGAL::to_double(p.x());
					const double y = CGAL::to_double(p.y());

					tmp_points[count] = Point_2ft(x, y);
				}

				assert(num_points == count);

				Line_2ft tmp_line;
				CGAL::linear_least_squares_fitting_2(tmp_points.begin(), tmp_points.end(), tmp_line, CGAL::Dimension_tag<0>());
				line = Line_2(static_cast<FT>(tmp_line.a()), static_cast<FT>(tmp_line.b()), static_cast<FT>(tmp_line.c()));
			}

			void create_thin_point_with_line(Projected_points &thinned_points, const Projected_point &query, const Line_2 &line) const {

				thinned_points[query.first] = m_simple_utils.project_onto_line(line, query.second);
			}


		// Complex thinning.
		private:
			struct Debug_data {
				std::vector<size_t> dims;   // detected dimensions
				std::vector<size_t> clusts; // number of clusters
			};

			int apply_complex_thinning(Boundary_data &, Projected_points &boundary_clutter_projected, FT &error, const Container &input) const {
				if (boundary_clutter_projected.empty()) {

					std::cerr << "WARNING: Thinning empty input!" << std::endl;
					return 0;
				}
				
				// Create kd tree.
				Fuzzy_tree tree;
				create_tree(tree, boundary_clutter_projected);

				// Get normals.
				Normals normals;
				get_normals(normals, boundary_clutter_projected, input);

				// Thin all points.
				Corners corners;
				Projected_points thinned_points;

				thin_all_points_with_normals(thinned_points, corners, tree, normals, boundary_clutter_projected);
				boundary_clutter_projected = thinned_points;

				// Save log.
				Log log; 
				log.export_projected_points_as_xyz("tmp" + std::string(PS) + "thinning_complex_result", boundary_clutter_projected, "unused path");

				error = compute_scale_error();
				return boundary_clutter_projected.size();
			}

			void get_normals(Normals &normals, const Projected_points &boundary_clutter_projected, const Container &input) const {

				// Add here new methods to estimate normals later!

				m_simple_utils.estimate_2d_normals_from_3d(normals, boundary_clutter_projected, input);
				assert(normals.size() == boundary_clutter_projected.size());

				Log log; 
				log.export_projected_points_with_normals_as_xyz("tmp" + std::string(PS) + "complex_with_normals", boundary_clutter_projected, normals, "unused path");
			}

			void thin_all_points_with_normals(Projected_points &thinned_points, Corners &corners, 
				const Fuzzy_tree &tree, const Normals &normals, const Projected_points &boundary_clutter_projected) const {

				assert(!boundary_clutter_projected.empty());
				thinned_points.clear();

				Debug_data debug_data;
				for (Point_iterator pit = boundary_clutter_projected.begin(); pit != boundary_clutter_projected.end(); ++pit) {
					
					const Projected_point &projected = *pit;
					handle_projected_point_with_normals(thinned_points, corners, debug_data, tree, normals, projected);
				}

				assert(!debug_data.dims.empty() && debug_data.dims.size() == debug_data.clusts.size());
				assert(debug_data.dims.size() == boundary_clutter_projected.size());

				// Save log data.
				Log log;
				log.save_dimensions_as_ply("tmp" + std::string(PS) + "complex_dimensions", boundary_clutter_projected, debug_data.dims  , "unused_path");
				log.save_num_clusters_as_ply("tmp" + std::string(PS) + "complex_num_clusters", boundary_clutter_projected, debug_data.clusts, "unused_path");
			}

			void handle_projected_point_with_normals(Projected_points &thinned_points, Corners &corners, Debug_data &debug_data,
			const Fuzzy_tree &tree, const Normals &normals, const Projected_point &query) const {

				Projected_points neighbours;
				find_nearest_neighbours(neighbours, tree, query);

				switch (m_fitter_type) {

					case Thinning_fitter_type::LINE:
						handle_with_line_complex(thinned_points, corners, debug_data, query, neighbours, normals);
						break;

					default:
						assert(!"Wrong fitter type!");
						break;
				}
			}

			// The data structure "points" also includes the query point!
			void handle_with_line_complex(Projected_points & thinned_points, Corners &corners, Debug_data &debug_data,
			const Projected_point &query, const Projected_points &points, const Normals &normals) const {

				assert(points.size() >= m_min_knn);
				assert(!normals.empty());

				std::vector<Projected_points> clusters;
				const size_t dimension = detect_local_dimension(query, points, normals, clusters, debug_data);

				switch (dimension) {
					case 0:
						handle_with_line_dimension_0(thinned_points, points, query);
						break;

					case 1:
						handle_with_line_dimension_1(thinned_points, corners, clusters, query);
						break;

					case 2:
						handle_with_line_dimension_2(thinned_points, points, clusters, query);
						break;

					default:
						assert(!"Wrong dimension in thinning!"); // can be commented out, see also below!
						break;
				}
			}

			size_t detect_local_dimension(const Projected_point &query, const Projected_points &points, const Normals &normals, 
			std::vector<Projected_points> &clusters, Debug_data &debug_data) const {
				
				const size_t default_dimension    = 100;
				const size_t default_num_clusters = 0;

				size_t dimension    = default_dimension;
				size_t num_clusters = default_num_clusters;

				if (is_dimension_0(query, points)) {
				
					dimension    = 0;
					num_clusters = 1;
				}

				if (is_dimension_2(points, clusters)) {

					dimension    = 2;
					num_clusters = clusters.size() > 0 ? clusters.size() : 1;
				}

				if (is_dimension_1(points, normals, clusters)) {

					dimension    = 1;
					num_clusters = clusters.size();
				}

				add_debug_data(debug_data, dimension, num_clusters);
				assert(dimension != default_dimension); // can be commented out, see also above!

				return dimension;
			}

			bool is_dimension_0(const Projected_point &query, const Projected_points &points) const {

				assert(!points.empty());

				const FT average_spacing = compute_average_spacing_from(query, points);
				const FT width = FT(2) * m_radius;

				if (average_spacing < m_dim_0_max * width) return true;
				return false;
			}

			bool is_dimension_1(const Projected_points &points, const Normals &normals, std::vector<Projected_points> &clusters) const {

				assert(!points.empty() && !normals.empty());
				clusters.clear();

				if (points.empty() || normals.empty()) {
					std::cerr << "Error: points/normals are empty, is_dimension_1 function thinning!" << std::endl;
					exit(EXIT_FAILURE);
				}

				// to be implemented!
				// implement clustering based on normals here or use minshift instead!
				// dimension_0 and dimension_2 checks can be turned off!

				return false;
			}

			bool is_dimension_2(const Projected_points &points, std::vector<Projected_points> &clusters) const {

				assert(!points.empty());
				clusters.clear();

				if (is_sparse(points)) 	 							 return true;
				if (has_multiple_spatial_clusters(points, clusters)) return true;

				return false;
			}

			bool is_sparse(const Projected_points &neighbours) const {

				Projected_points points = neighbours;

				Boundary_data stub;
				Level_of_detail_grid_simplify<Kernel, Boundary_data, Projected_points> grid_simplify;

				const FT width 			  = FT(2) * m_radius;
				const FT grid_cell_length = width / m_dim_2_sparsity_num_cells; 

				grid_simplify.save_result(false);
				grid_simplify.set_grid_cell_length(grid_cell_length);
				grid_simplify.process(stub, points);

				const FT num_filled_cells = static_cast<FT>(points.size());
				const FT num_all_cells 	  = estimate_number_of_cells();
				
				assert(num_all_cells >= num_filled_cells);
				assert(num_all_cells != FT(0));

				const FT ratio = num_filled_cells / num_all_cells;
				assert(ratio >= FT(0) && ratio <= FT(1));

				if (ratio > m_dim_2_sparsity_max) return true;
				return false;
			}

			FT estimate_number_of_cells() const {

				switch (m_neighbour_search_type) {

					case Neighbour_search_type::CIRCLE:
						return estimate_number_of_cells_for_circle();

					case Neighbour_search_type::SQUARE:
						return estimate_number_of_cells_for_square();

					case Neighbour_search_type::KNN:
						return estimate_number_of_cells_for_min_box();

					default:
						assert(!"Wrong neighbour search type!");
						break;
				}
				return FT(0);
			}

			FT estimate_number_of_cells_for_circle() const {
				return estimate_number_of_cells_for_square();
			}

			FT estimate_number_of_cells_for_square() const {

				const  FT adjust = FT(4);
				return m_dim_2_sparsity_num_cells * m_dim_2_sparsity_num_cells + adjust;
			}

			FT estimate_number_of_cells_for_min_box() const {
				return estimate_number_of_cells_for_square();
			}

			bool has_multiple_spatial_clusters(const Projected_points &points, std::vector<Projected_points> &clusters) const {

				assert(m_max_num_clusters > 0);
				assert(!points.empty() && clusters.empty());
				
				Kmeans_data kmeans_data; Kmeans_map kmeans_map;
				set_kmeans_data(kmeans_data, kmeans_map, points);

				std::vector<Kmeans_result> means(m_max_num_clusters);
				std::vector<FT> 		   errors(m_max_num_clusters);

				for (size_t i = 1; i <= m_max_num_clusters; ++i) {	
					if (points.size() < i) return false;

					run_kmeans(kmeans_data, i, means[i-1]);
					errors[i-1] = compute_kmeans_error(kmeans_data, means[i-1], i);
				}

				const size_t num_clusters_found = get_kmeans_number_of_clusters(errors);
				assert(num_clusters_found > 0 && num_clusters_found <= m_max_num_clusters);

				const size_t cluster_index = num_clusters_found - 1;
				assert(cluster_index >= 0 && cluster_index < m_max_num_clusters);

				auto indices = std::get<1>(means[cluster_index]);
				set_kmeans_clusters(clusters, indices, kmeans_map, points, num_clusters_found);

				return true;
			}

			size_t return_global_index() const {
				return rand() % 1000000 + 1;
			}

			size_t get_kmeans_number_of_clusters(const std::vector<FT> &errors) const {

				// Plot errors.
				// const std::string name = "tmp" + std::string(PS) + "plots" + std::string(PS) + "errors_" + std::to_string(return_global_index());
				// Log log; log.plot_2d(name, errors);

				// Add other methods here!

				return get_kmeans_number_of_clusters_min_error(errors);
			}

			size_t get_kmeans_number_of_clusters_min_error(const std::vector<FT> &errors) const {

				const size_t num_errors = errors.size();

				assert(num_errors > 0);				
				if (num_errors == 1) return num_errors;
				assert(num_errors > 1);

				FT min_error = errors[0];
				size_t best_error_index = 0;

				for (size_t i = 1; i < num_errors; ++i) {
					if (errors[i] < min_error) {

						best_error_index = i;
						min_error = errors[i];
					}
				}

				return best_error_index + 1;
			}

			size_t get_kmeans_number_of_clusters_descend_error(const std::vector<FT> &errors) const {

				const size_t num_errors = errors.size();
				
				assert(num_errors > 0);				
				if (num_errors == 1) return num_errors;

				assert(num_errors > 1);
				FT curr_error = errors[0];

				for (size_t i = 2; i <= num_errors; ++i) {
					const FT next_error = errors[i-1];
					
					if (next_error > curr_error) return i-1;
					else curr_error = next_error;
				}

				return num_errors;
			}

			FT compute_kmeans_error(const Kmeans_data &kmeans_data, const Kmeans_result &means, const size_t num_clusters_expected) const {

				// Add other methods here!

				return compute_kmeans_error_elbow(kmeans_data, means, num_clusters_expected);
			}

			// Here we basically use elbow method.
			// http://www.sthda.com/english/articles/29-cluster-validation-essentials/96-determining-the-optimal-number-of-clusters-3-must-know-methods/
			FT compute_kmeans_error_elbow(const Kmeans_data &kmeans_data, const Kmeans_result &means, const size_t num_clusters_expected) const {

				const FT inertia = static_cast<FT>(dkm::means_inertia(kmeans_data, means, num_clusters_expected));

				const FT tmp = static_cast<FT>(num_clusters_expected);
				const FT error = inertia * tmp;

				assert(error >= FT(0));
				return error;
			}

			void set_kmeans_data(Kmeans_data &kmeans_data, Kmeans_map &kmeans_map, const Projected_points &points) const {

				const size_t num_points = points.size();

				kmeans_data.clear(); 			kmeans_map.clear();
				kmeans_data.resize(num_points); kmeans_map.resize(num_points);

				size_t count = 0;
				for (Point_iterator pit = points.begin(); pit != points.end(); ++pit, ++count) {
					const Point_2 &p = (*pit).second;

					kmeans_data[count] = {{CGAL::to_double(p.x()), CGAL::to_double(p.y())}};
					kmeans_map[count]  = (*pit).first;
				}
				assert(count == num_points);
			}

			void run_kmeans(const Kmeans_data &kmeans_data, const size_t num_clusters_expected, Kmeans_result &means) const {

				assert(num_clusters_expected > 0);
				assert(!kmeans_data.empty());

				const size_t num_points = kmeans_data.size();
				assert(num_points >= num_clusters_expected);

				if (num_points < num_clusters_expected) {
					std::cerr << "Error: num_points < num_clusters_expected, run_kmeans function thinning!" << std::endl;
					exit(EXIT_FAILURE);
				}

				means = dkm::kmeans_lloyd(kmeans_data, num_clusters_expected);
			}

			template<class Indices>
			void set_kmeans_clusters(std::vector<Projected_points> &clusters, 
				const Indices &indices, const Kmeans_map &kmeans_map, const Projected_points &points, const size_t num_clusters_expected) const {

				clusters.clear();
				clusters.resize(num_clusters_expected);

				assert(indices.size() == kmeans_map.size());
				for (size_t i = 0; i < indices.size(); ++i) {

					switch (indices[i]) {

						case 0: 
							assert(clusters.size() > 0);
							clusters[0][kmeans_map[i]] = points.at(kmeans_map[i]); 
							break;
						
						case 1:
							assert(clusters.size() > 1);
							clusters[1][kmeans_map[i]] = points.at(kmeans_map[i]); 
							break;
						
						case 2:
							assert(clusters.size() > 2);
							clusters[2][kmeans_map[i]] = points.at(kmeans_map[i]); 
							break;
						
						case 3:
							assert(clusters.size() > 3);
							clusters[3][kmeans_map[i]] = points.at(kmeans_map[i]); 
							break;

						case 4:
							assert(clusters.size() > 4);
							clusters[4][kmeans_map[i]] = points.at(kmeans_map[i]); 
							break;

						default: 
							assert(!"Wrong cluster index! Number of clusters should be [1, 5].");
							break;
					}
				}

				// Print clusters.
				// const std::string name = "tmp" + std::string(PS) + "plots" + std::string(PS) + "clusters_" + std::to_string(return_global_index()) + "_" + std::to_string(num_clusters_expected);
				// Log log; log.draw_clusters(name, clusters);
			}

			FT compute_average_spacing_from(const Projected_point &query, const Projected_points &points) const {

				assert(!points.empty());

				FT average_spacing = FT(0);
				FT num_dists 	   = FT(0);

				for (Point_iterator pit = points.begin(); pit != points.end(); ++pit) {
					
					const FT dist = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance(query.second, (*pit).second))));
					if (dist != FT(0)) {

						average_spacing += dist;
						num_dists       += FT(1);
					}
				}

				assert(average_spacing >= FT(0) && num_dists != FT(0));
				average_spacing /= num_dists;

				return average_spacing;
			}

			void add_debug_data(Debug_data &debug_data, const size_t dimension, const size_t num_clusters) const {

				debug_data.dims.push_back(dimension);
				debug_data.clusts.push_back(num_clusters);
			}

			void handle_with_line_dimension_0(Projected_points & /* thinned_points */, 
			const Projected_points & /* points */, const Projected_point & /* query */) const {

				// To be implemented later.
				return;
			}

			void handle_with_line_dimension_1(Projected_points &thinned_points, Corners &corners, 
				const std::vector<Projected_points> &clusters, const Projected_point &query) const {

				const size_t num_clusters = clusters.size();
				assert(num_clusters > 0);
				
				int best_cluster_index = -1; FT min_dist = FT(10000000);
				std::vector<Line_2> lines(num_clusters);

				fit_line_to_points(clusters[0], lines[0]);
				best_cluster_index = get_best_cluster_index(query, lines[0], 0, best_cluster_index, min_dist);

				for (size_t i = 1; i < num_clusters; ++i) {
					assert(i < num_clusters);

					fit_line_to_points(clusters[i], lines[i]);
					best_cluster_index = get_best_cluster_index(query, lines[i], i, best_cluster_index, min_dist);
				}

				assert(best_cluster_index >= 0 && best_cluster_index < static_cast<int>(lines.size()));
				create_thin_point_with_line(thinned_points, query, lines[best_cluster_index]);

				intersect_clusters(corners, clusters, lines);
			}

			int get_best_cluster_index(const Projected_point &query, const Line_2 &line, 
				const int cluster_index, const int best_cluster_index, FT &min_dist) const {

				const Point_2 projected = m_simple_utils.project_onto_line(line, query.second);
				const FT dist = squared_distance(projected, query.second);

				if (dist < min_dist) {
					min_dist = dist;

					return cluster_index;
				}
				return best_cluster_index;
			}

			void intersect_clusters(Corners & /* corners */, 
			const std::vector<Projected_points> & /* clusters */, const std::vector<Line_2> & /* lines */) const {

				// To be implemented later.
				return;
			}

			void handle_with_line_dimension_2(Projected_points & /* thinned_points */, 
			const Projected_points & /* points */, const std::vector<Projected_points> & /* clusters */, const Projected_point & /* query */) const {
				
				// To be implemented later.
				return;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_THINNING_H