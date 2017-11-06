#ifndef CGAL_LEVEL_OF_DETAIL_CLUTTER_PROCESSOR_H
#define CGAL_LEVEL_OF_DETAIL_CLUTTER_PROCESSOR_H

// STL includes.
#include <map>
#include <vector>
#include <cassert>
#include <cmath>

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

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		// Grid based simplification.
		template<class KernelTraits, class BoundaryData, class ProjectedPoints>
		class Level_of_detail_grid_simplify {

		public:
			// Fields.
			typedef KernelTraits 	Kernel;
			typedef BoundaryData 	Boundary_data;
			typedef ProjectedPoints Projected_points;

			typedef typename Kernel::FT 	 FT;
			typedef typename Kernel::Point_2 Point_2;
			typedef typename Kernel::Line_2  Line_2;

			typedef std::pair<int, Point_2> 					 				Projected_point;
			typedef typename CGAL::Second_of_pair_property_map<Projected_point> Point_map;

			typedef CGAL::Search_traits_2<Kernel>                       					 Search_traits_2;
			typedef CGAL::Search_traits_adapter<Projected_point, Point_map, Search_traits_2> Search_traits;
			typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> 					     Neighbor_search;
			typedef CGAL::Fuzzy_sphere<Search_traits>                    					 Fuzzy_circle;
			typedef CGAL::Kd_tree<Search_traits>					                         Fuzzy_tree;

			typedef typename Boundary_data::const_iterator    Boundary_iterator;
			typedef typename Projected_points::const_iterator Point_iterator;
			typedef typename Neighbor_search::iterator  	  Neighbour_iterator;

			typename Kernel::Compute_squared_distance_2 squared_distance;


			// Extra.
			using Log = CGAL::LOD::Mylog;

			using Cell_id  = std::pair<int, int>;
			using Grid_map = std::map<Cell_id, std::vector<int> >;

			using Grid_cell_iterator = Grid_map::const_iterator;


			// Constructor.
			Level_of_detail_grid_simplify() : 
			m_grid_cell_length(FT(1) / FT(100)),
			m_new_point_type(Grid_new_point_type::BARYCENTRE) { }


			// Input functions.		
			void set_grid_cell_length(const FT new_cell_length) {

				assert(new_cell_length > FT(0));
				m_grid_cell_length = new_cell_length;
			}

			void set_new_point_type(const Grid_new_point_type new_point_type) {
				m_new_point_type = new_point_type;
			}


			// Main.
			int process(Boundary_data &, Projected_points &boundary_clutter_projected) const {

				int number_of_removed_points = 0;
				if (boundary_clutter_projected.empty()) {
					
					std::cerr << "WARNING: Grid simplify empty input!" << std::endl;
					return number_of_removed_points;
				}

				Grid_map grid_map;
				create_grid_map(grid_map, boundary_clutter_projected);

				Projected_points cleaned_points;
				simplify_clutter_using_grid(cleaned_points, grid_map, boundary_clutter_projected);

				number_of_removed_points = static_cast<int>(boundary_clutter_projected.size() - cleaned_points.size());
				boundary_clutter_projected = cleaned_points;

				assert(number_of_removed_points >= 0);
				// set_new_boundary_data(boundary_clutter, boundary_clutter_projected);

				Log log; 
				log.export_projected_points_as_xyz("tmp/grid_simplify_result", boundary_clutter_projected, "unused path");

				return number_of_removed_points;
			}


		// Fields.
		private:
			FT 			   		m_grid_cell_length;
			Grid_new_point_type m_new_point_type;


		// Grid simplify.
		private:
			void create_grid_map(Grid_map &grid_map, const Projected_points &boundary_clutter_projected) const {

				assert(grid_map.empty());
				assert(!boundary_clutter_projected.empty());

				Cell_id cell_id;
				for (Point_iterator pit = boundary_clutter_projected.begin(); pit != boundary_clutter_projected.end(); ++pit) {
					
					const Projected_point &projected = *pit;
					get_cell_id(cell_id, projected);

					grid_map[cell_id].push_back(projected.first);
				}
			}

			void get_cell_id(Cell_id &cell_id, const Projected_point &projected) const {

				const Point_2 &point = projected.second;

				const int id_x = get_id_value(point.x());
				const int id_y = get_id_value(point.y());

				cell_id = std::make_pair(id_x, id_y);
			}

			int get_id_value(const FT value) const {

				assert(m_grid_cell_length > FT(0));
				const int id = static_cast<int>(value / m_grid_cell_length);

				if (value >= 0) return id;
				return id - 1;
			}

			void simplify_clutter_using_grid(Projected_points &cleaned_points, const Grid_map &grid_map, const Projected_points &boundary_clutter_projected) const {

				assert(cleaned_points.empty());
				assert(!grid_map.empty() && !boundary_clutter_projected.empty());

				int point_index = 0; Point_2 new_point;
				for (Grid_cell_iterator git = grid_map.begin(); git != grid_map.end(); ++git, ++point_index) {
					
					const auto &grid_cell = *git;
					create_new_point(new_point, grid_cell, boundary_clutter_projected);

					cleaned_points[point_index] = new_point;
				}		
			}

			template<class Grid_cell>
			void create_new_point(Point_2 &new_point, const Grid_cell &grid_cell, const Projected_points &boundary_clutter_projected) const {

				switch (m_new_point_type) {

					case Grid_new_point_type::CENTROID:
						create_new_centroid(new_point, grid_cell);
						break;

					case Grid_new_point_type::BARYCENTRE:
						create_new_barycentre(new_point, grid_cell, boundary_clutter_projected);
						break;

					case Grid_new_point_type::CLOSEST:
						create_new_closest(new_point, grid_cell, boundary_clutter_projected);
						break;

					default:
						assert(!"Wrong new point type!");
						break;
				}
			}

			template<class Grid_cell>
			void create_new_centroid(Point_2 &new_point, const Grid_cell &grid_cell) const {

				const Cell_id &cell_id = grid_cell.first;

				const FT half_cell_length = m_grid_cell_length / FT(2);

				const FT x = static_cast<FT>(cell_id.first)  * m_grid_cell_length + half_cell_length;
				const FT y = static_cast<FT>(cell_id.second) * m_grid_cell_length + half_cell_length;

				new_point = Point_2(x, y);
			}

			template<class Grid_cell>
			void create_new_barycentre(Point_2 &new_point, const Grid_cell &grid_cell, const Projected_points &boundary_clutter_projected) const {

				const std::vector<int> &point_idxs = grid_cell.second;
				const size_t num_point_idxs = point_idxs.size();

				assert(num_point_idxs != 0);

				FT x = FT(0), y = FT(0);
				for (size_t i = 0; i < num_point_idxs; ++i) {

					const Point_2 &old_point = boundary_clutter_projected.at(point_idxs[i]);

					x += old_point.x();
					y += old_point.y();
				}

				x /= static_cast<FT>(num_point_idxs);
				y /= static_cast<FT>(num_point_idxs);

				new_point = Point_2(x, y);
			}

			template<class Grid_cell>
			void create_new_closest(Point_2 &new_point, const Grid_cell &grid_cell, const Projected_points &boundary_clutter_projected) const {

				Point_2 barycentre;
				create_new_barycentre(barycentre, grid_cell, boundary_clutter_projected);

				const std::vector<int> &point_idxs = grid_cell.second;
				const size_t num_point_idxs = point_idxs.size();

				assert(num_point_idxs != 0);
				int new_point_index = -1;

				FT min_dist = FT(1000000000);
				for (size_t i = 0; i < num_point_idxs; ++i) {
					
					const Point_2 &old_point = boundary_clutter_projected.at(point_idxs[i]);
					const FT dist = CGAL::sqrt(squared_distance(old_point, barycentre));

					if (dist < min_dist) {
						
						min_dist = dist;
						new_point_index = i;
					}
				}

				assert(new_point_index >= 0 && new_point_index < static_cast<int>(num_point_idxs));
				new_point = boundary_clutter_projected.at(point_idxs[new_point_index]);
			}

			void set_new_boundary_data(Boundary_data &boundary_clutter, const Projected_points &boundary_clutter_projected) const {

				boundary_clutter.clear();

				Boundary_data new_data;
				assert(boundary_clutter.empty() && !boundary_clutter_projected.empty());

				size_t count = 0;
				std::vector<int> idxs(boundary_clutter_projected.size());

				for (Point_iterator pit = boundary_clutter_projected.begin(); pit != boundary_clutter_projected.end(); ++pit, ++count) {	
					const Projected_point &projected = *pit;
					
					assert(count < idxs.size());
					idxs[count] = projected.first;
				}
				
				assert(count == idxs.size());

				new_data[0] = idxs;
				boundary_clutter = new_data;
			}
		};


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
			typedef typename Kernel::Vector_3 Normal;

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


			// Constructor.
			Level_of_detail_thinning() : 
			m_k(3), 
			m_fitter_type(Thinning_fitter_type::LINE), 
			m_neighbour_search_type(Neighbour_search_type::CIRCLE),
			m_radius(FT(1) / FT(4)),
			m_thinning_type(Thinning_type::NAIVE),
			m_scale_type(Thinning_scale_type::FIXED),
			m_min_knn(2),
			m_dim_0_thresh(FT(5) / FT(100)),
			m_dim_2_sparse_thresh_1(FT(5)),
			m_dim_2_sparse_thresh_2(FT(7) / FT(10)) { }


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

			void set_param_1(const FT new_value) {
				m_dim_2_sparse_thresh_1 = new_value;
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
						break;
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
			
			FT m_dim_0_thresh; 		    // percentage
			FT m_dim_2_sparse_thresh_1; // number of cells
			FT m_dim_2_sparse_thresh_2; // total percentage

			void choose_scale_adaptively() const {

				// To be implemented later.
				assert(!"Adaptive method is not implemented!");
			}

			int fixed_process(Boundary_data &boundary_clutter, Projected_points &boundary_clutter_projected, const Container &input) const {

				FT stub;
				switch (m_thinning_type) {

					case Thinning_type::NAIVE:
						return apply_naive_thinning(boundary_clutter, boundary_clutter_projected, stub);
						break;

					case Thinning_type::COMPLEX:
						return apply_complex_thinning(boundary_clutter, boundary_clutter_projected, stub, input);
						break;

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

				Log log; 
				log.export_projected_points_as_xyz("tmp/thinning_naive_result", boundary_clutter_projected, "unused path");

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

			void handle_projected_point_naive(Projected_points &thinned_points, const Fuzzy_tree &tree, const Projected_point &projected) const {

				Projected_points neighbours;
				const Point_2 &query = projected.second;
				find_nearest_neighbours(neighbours, tree, query);

				switch (m_fitter_type) {

					case Thinning_fitter_type::LINE:
						handle_with_line_naive(thinned_points, projected, neighbours);
						break;

					default:
						assert(!"Wrong fitter type!");
						break;
				}
			}

			void find_nearest_neighbours(Projected_points &neighbours, const Fuzzy_tree &tree, const Point_2 &query) const {
				
				switch (m_neighbour_search_type) {

					case Neighbour_search_type::KNN:
						find_nearest_neighbours_with_knn(neighbours, tree, query, m_k);
						break;

					case Neighbour_search_type::CIRCLE:
						find_nearest_neighbours_with_circle(neighbours, tree, query);
						break;

					case Neighbour_search_type::SQUARE:
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

			void find_nearest_neighbours_with_knn(Projected_points &neighbours, const Fuzzy_tree &tree, const Point_2 &query, const size_t num_knn) const {

				assert(num_knn >= m_min_knn);
				
				Neighbor_search search(tree, query, num_knn);
				const size_t num_points = static_cast<size_t>(std::distance(search.begin(), search.end()));

				neighbours.clear();
				size_t count = 0;

				for (Neighbour_iterator nit = search.begin(); nit != search.end(); ++nit, ++count) neighbours[nit->first.first] = nit->first.second;
				assert(count == num_points);
			}

			void find_nearest_neighbours_with_circle(Projected_points &neighbours, const Fuzzy_tree &tree, const Point_2 &query) const {

				neighbours.clear();

				Fuzzy_circle circle;
				compute_circle(query, circle);

				tree.search(std::inserter(neighbours, neighbours.end()), circle);
			}

			void find_nearest_neighbours_with_square(Projected_points &neighbours, const Fuzzy_tree &tree, const Point_2 &query) const {

				neighbours.clear();

				Fuzzy_square square;
				compute_square(query, square);

				tree.search(std::inserter(neighbours, neighbours.end()), square);
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

				const size_t num_points = neighbours.size();
				std::vector<Point_2> points(num_points);

				size_t count = 0;
				for (Point_iterator pit = neighbours.begin(); pit != neighbours.end(); ++pit, ++count) points[count] = (*pit).second;
				assert(num_points == count);

				CGAL::linear_least_squares_fitting_2(points.begin(), points.end(), line, CGAL::Dimension_tag<0>());
			}

			void create_thin_point_with_line(Projected_points &thinned_points, const Projected_point &query, const Line_2 &line) const {

				thinned_points[query.first] = project_onto_line(query.second, line);
			}

			// My custom function to handle precision problems when projecting points.
			Point_2 project_onto_line(const Point_2 &p, const Line_2 &line) const {

				const auto a = line.point(0);
				const auto b = line.point(1);

				const auto projected = a + CGAL::scalar_product(p - a, b - a) / CGAL::scalar_product(b - a, b - a) * (b - a);
				
				if (std::isnan(projected.x()) || std::isnan(projected.y())) return line.projection(p);
				else return projected;
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
				log.export_projected_points_as_xyz("tmp/thinning_complex_result", boundary_clutter_projected, "unused path");

				error = compute_scale_error();
				return boundary_clutter_projected.size();
			}

			void get_normals(Normals &normals, const Projected_points &boundary_clutter_projected, const Container &input) const {

				// Add here new methods to estimate normals later!

				project_normals(normals, boundary_clutter_projected, input);
				assert(normals.size() == boundary_clutter_projected.size());

				Log log; 
				log.export_projected_points_with_normals_as_xyz("tmp/complex_with_normals", boundary_clutter_projected, normals, "unused path");
			}

			void project_normals(Normals &normals, const Projected_points &boundary_clutter_projected, const Container &input) const {

				normals.clear();
				assert(!boundary_clutter_projected.empty() && input.number_of_points() != 0);

				for (Point_iterator pit = boundary_clutter_projected.begin(); pit != boundary_clutter_projected.end(); ++pit) {
					
					const Projected_point &projected = *pit;
					const int point_index = projected.first;
					
					project_normal(normals, point_index, input);
				}
			}

			void project_normal(Normals &normals, const int point_index, const Container &input) const {

				const Normal plane_normal  = Normal(FT(0), FT(0), FT(1));
				const Normal &point_normal = input.normal(point_index);

				const Normal projected_normal = point_normal - dot_product(point_normal, plane_normal) * plane_normal;
				normals[point_index] = projected_normal;
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
				log.save_dimensions_as_ply(    "tmp/complex_dimensions", boundary_clutter_projected, debug_data.dims  , "unused_path");
				log.save_num_clusters_as_ply("tmp/complex_num_clusters", boundary_clutter_projected, debug_data.clusts, "unused_path");
			}

			void handle_projected_point_with_normals(Projected_points &thinned_points, Corners &corners, Debug_data &debug_data,
			const Fuzzy_tree &tree, const Normals &normals, const Projected_point &projected) const {

				Projected_points neighbours;
				const Point_2 &query = projected.second;
				find_nearest_neighbours(neighbours, tree, query);

				switch (m_fitter_type) {

					case Thinning_fitter_type::LINE:
						handle_with_line_complex(thinned_points, corners, debug_data, projected, neighbours, normals);
						break;

					default:
						assert(!"Wrong fitter type!");
						break;
				}
			}

			void handle_with_line_complex(Projected_points & thinned_points, Corners &corners, Debug_data &debug_data,
			const Projected_point &query, const Projected_points &neighbours, const Normals &normals) const {

				assert(neighbours.size() >= m_min_knn);
				assert(!normals.empty());

				std::vector<Projected_points> clusters;
				const size_t dimension = detect_local_dimension(query, neighbours, normals, clusters, debug_data);

				switch (dimension) {
					case 0:
						handle_with_line_dimension_0(thinned_points, neighbours, query);
						break;

					case 1:
						handle_with_line_dimension_1(thinned_points, corners, clusters, query);
						break;

					case 2:
						handle_with_line_dimension_2(thinned_points, neighbours, query);
						break;

					default:
						// assert(!"Wrong dimension in thinning!"); // put back!
						break;
				}
			}

			size_t detect_local_dimension(const Projected_point &query, const Projected_points &neighbours, const Normals &normals, 
			std::vector<Projected_points> &clusters, Debug_data &debug_data) const {
				
				const size_t default_dimension    = 100;
				const size_t default_num_clusters = 0;

				size_t dimension    = default_dimension;
				size_t num_clusters = default_num_clusters;

				if (is_dimension_0(query, neighbours)) {
				
					dimension    = 0;
					num_clusters = 1;
				}

				if (is_dimension_2(query, neighbours, clusters)) {

					dimension    = 2;
					num_clusters = clusters.size() > 0 ? clusters.size() : 1;
				}

				if (is_dimension_1(query, neighbours, normals, clusters)) {

					dimension    = 1;
					num_clusters = clusters.size();
				}

				add_debug_data(debug_data, dimension, num_clusters);
				// assert(dimension != default_dimension); // put back!

				return dimension;
			}

			bool is_dimension_0(const Projected_point &query, const Projected_points &neighbours) const {

				assert(!neighbours.empty());

				const FT average_spacing = compute_average_spacing_from(query, neighbours);
				const FT width = FT(2) * m_radius;

				if (average_spacing < m_dim_0_thresh * width) return true;
				return false;
			}

			bool is_dimension_1(const Projected_point & /* query */, const Projected_points &neighbours, const Normals &normals, 
			std::vector<Projected_points> &clusters) const {

				assert(!neighbours.empty() && !normals.empty());
				clusters.clear();

				// to be implemented!

				return false;
			}

			bool is_dimension_2(const Projected_point &query, const Projected_points &neighbours, std::vector<Projected_points> &clusters) const {

				assert(!neighbours.empty());
				clusters.clear();

				if (is_sparse(query, neighbours)) 	 							return true;
				if (has_multiple_spatial_clusters(query, neighbours, clusters)) return true;

				return false;
			}

			bool is_sparse(const Projected_point &query, const Projected_points &neighbours) const {

				Level_of_detail_grid_simplify<Kernel, Boundary_data, Projected_points> grid_simplify;
				
				Boundary_data stub;
				Projected_points points = neighbours;
				points[query.first]     = query.second;

				const FT width 			  = FT(2) * m_radius;
				const FT grid_cell_length = width / m_dim_2_sparse_thresh_1; 

				grid_simplify.set_grid_cell_length(grid_cell_length);
				grid_simplify.process(stub, points);

				const FT num_filled_cells = static_cast<FT>(points.size());
				const FT num_all_cells 	  = estimate_number_of_cells();
				
				assert(num_all_cells >= num_filled_cells);
				assert(num_all_cells != FT(0));

				const FT ratio = num_filled_cells / num_all_cells;
				if (ratio > m_dim_2_sparse_thresh_2) return true;

				return false;
			}

			FT estimate_number_of_cells() const {

				switch (m_neighbour_search_type) {

					case Neighbour_search_type::CIRCLE:
						return estimate_number_of_cells_for_circle();
						break;

					case Neighbour_search_type::SQUARE:
						return estimate_number_of_cells_for_square();
						break;

					case Neighbour_search_type::KNN:
						return estimate_number_of_cells_for_min_box();
						break;

					default:
						assert(!"Wrong neighbour search type!");
						break;
				}
				return 0;
			}

			FT estimate_number_of_cells_for_circle() const {
				return estimate_number_of_cells_for_square();
			}

			FT estimate_number_of_cells_for_square() const {
				return m_dim_2_sparse_thresh_1 * m_dim_2_sparse_thresh_1 + 4;
			}

			FT estimate_number_of_cells_for_min_box() const {
				return estimate_number_of_cells_for_square();
			}

			bool has_multiple_spatial_clusters(const Projected_point & /* query */, const Projected_points & /* points */, std::vector<Projected_points> & /* clusters */) const {

				return false;
			}

			FT compute_average_spacing_from(const Projected_point &query, const Projected_points &points) const {

				assert(!points.empty());

				FT average_spacing = FT(0);
				FT num_dists 	   = FT(0);

				for (Point_iterator pit = points.begin(); pit != points.end(); ++pit) {
					
					const FT dist = CGAL::sqrt(squared_distance(query.second, (*pit).second));
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

			void handle_with_line_dimension_0(Projected_points & /* thinned_points */, const Projected_points & /* neighbours */, const Projected_point & /* query */) const {

				// To be implemented later.
				return;
			}

			void handle_with_line_dimension_1(Projected_points &thinned_points, Corners &corners, 
				const std::vector<Projected_points> &clusters, const Projected_point &query) const {

				const size_t num_clusters = clusters.size();
				assert(num_clusters != 0);

				Line_2 line;
				for (size_t i = 0; i < num_clusters; ++i) {

					fit_line_to_points(clusters[i], line);
					create_thin_point_with_line(thinned_points, query, line);
				}

				intersect_clusters(corners, clusters);
			}

			void intersect_clusters(Corners & /* corners */, const std::vector<Projected_points> & /* clusters */) const {

				// To be implemented later.
				return;
			}

			void handle_with_line_dimension_2(Projected_points & /* thinned_points */, const Projected_points & /* neighbours */, const Projected_point & /* query */) const {
				
				// To be implemented later.
				return;
			}
		};


		// Main clutter class.
		template<class KernelTraits, class BoundaryData, class ProjectedPoints, class InputContainer>
		class Level_of_detail_clutter_processor {

		public:
			// Fields.
			typedef KernelTraits 	Kernel;
			typedef BoundaryData 	Boundary_data;
			typedef ProjectedPoints Projected_points;
			typedef InputContainer  Container;

			typedef typename Kernel::FT 	 FT;
			typedef typename Kernel::Point_2 Point_2;
			typedef typename Kernel::Line_2  Line_2;

			typedef std::pair<int, Point_2> 					 				Projected_point;
			typedef typename CGAL::Second_of_pair_property_map<Projected_point> Point_map;

			typedef CGAL::Search_traits_2<Kernel>                       					 Search_traits_2;
			typedef CGAL::Search_traits_adapter<Projected_point, Point_map, Search_traits_2> Search_traits;
			typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> 					     Neighbor_search;
			typedef CGAL::Fuzzy_sphere<Search_traits>                    					 Fuzzy_circle;
			typedef CGAL::Kd_tree<Search_traits>					                         Fuzzy_tree;

			typedef typename Boundary_data::const_iterator    Boundary_iterator;
			typedef typename Projected_points::const_iterator Point_iterator;
			typedef typename Neighbor_search::iterator  	  Neighbour_iterator;

			typename Kernel::Compute_squared_distance_2 squared_distance;


			// Extra.
			using Log = CGAL::LOD::Mylog;

			using Cell_id  = std::pair<int, int>;
			using Grid_map = std::map<Cell_id, std::vector<int> >;

			using Grid_cell_iterator = Grid_map::const_iterator;


			// Constructor.
			Level_of_detail_clutter_processor() : 
			m_k(3), 
			m_grid_cell_length(FT(1) / FT(100)),
			m_fitter_type(Thinning_fitter_type::LINE), 
			m_new_point_type(Grid_new_point_type::BARYCENTRE),
			m_neighbour_search_type(Neighbour_search_type::CIRCLE),
			m_radius(FT(1) / FT(4)),
			m_thinning_type(Thinning_type::NAIVE) { }


			// Grid simplify.			
			void set_grid_cell_length(const FT new_cell_length) {

				assert(new_cell_length > FT(0));
				m_grid_cell_length = new_cell_length;
			}

			void set_new_point_type(const Grid_new_point_type new_point_type) {
				m_new_point_type = new_point_type;
			}


			// Thinning.
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


			// Main.
			int process(Boundary_data &boundary_clutter, Projected_points &boundary_clutter_projected, const Container &input) const {

				auto number_of_removed_points = 0;

				apply_thinning(boundary_clutter, boundary_clutter_projected, input);
				number_of_removed_points = apply_grid_simplify(boundary_clutter, boundary_clutter_projected);

				return number_of_removed_points;
			}

			// Required parameters: only m_k - number of natural neighbours; m_fitter_type = LINE in general.
			int apply_thinning(Boundary_data &boundary_clutter, Projected_points &boundary_clutter_projected, const Container &input) const {

				Level_of_detail_thinning<Kernel, Boundary_data, Projected_points, Container> thinning;

				thinning.set_number_of_neighbours(m_k);
				thinning.set_fuzzy_radius(m_radius);
				thinning.set_thinning_type(m_thinning_type);
				thinning.set_neighbour_search_type(m_neighbour_search_type);
				thinning.set_fitter_type(m_fitter_type);

				return thinning.process(boundary_clutter, boundary_clutter_projected, input);
			}

			// This is like a moving least squares local algorithm. Later we can solve a global optimization here.
			// Required parameters: m_grid_cell_length - length of the cell side in the grid; here you can also use m_new_point_type.
			int apply_grid_simplify(Boundary_data &boundary_clutter, Projected_points &boundary_clutter_projected) const {

				Level_of_detail_grid_simplify<Kernel, Boundary_data, Projected_points> grid_simplify;

				grid_simplify.set_grid_cell_length(m_grid_cell_length);
				grid_simplify.set_new_point_type(m_new_point_type);

				return grid_simplify.process(boundary_clutter, boundary_clutter_projected);
			}


		// Fields.
		private:
			size_t 		   		  m_k;
			FT 			   		  m_grid_cell_length;
			Thinning_fitter_type  m_fitter_type;
			Grid_new_point_type   m_new_point_type;
			Neighbour_search_type m_neighbour_search_type;
			FT 					  m_radius;
			Thinning_type 		  m_thinning_type;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_CLUTTER_PROCESSOR_H