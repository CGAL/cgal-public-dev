#ifndef CGAL_LEVEL_OF_DETAIL_REGION_GROWING_2_H
#define CGAL_LEVEL_OF_DETAIL_REGION_GROWING_2_H

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

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Utils/Level_of_detail_utils.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class BoundaryData, class ProjectedPoints, class InputContainer>
		class Level_of_detail_region_growing_2 {

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

			typedef typename Boundary_data::const_iterator    Boundary_iterator;
			typedef typename Projected_points::const_iterator Point_iterator;
			typedef typename Neighbor_search::iterator  	  Neighbour_iterator;

			typedef Level_of_detail_utils_simple<Kernel> Simple_utils;

			typename Kernel::Compute_squared_distance_2 squared_distance;
			typename Kernel::Compute_scalar_product_3   dot_product;


			// Extra.
			using Log 	  		 = CGAL::LOD::Mylog;
			using Normals 		 = std::map<int, Normal>;
			using Sorted_indices = std::vector<int>;


			// Constructor.
			Level_of_detail_region_growing_2() : 
			m_epsilon(-FT(1)),
			m_cluster_epsilon(-FT(1)),
			m_normal_threshold(-FT(1)),
			m_min_points(0),
			m_save_info(false) { }


			// Public functions.
			void set_epsilon(const FT new_value) {

				assert(new_value > FT(0));
				m_epsilon = new_value;
			}

			void set_cluster_epsilon(const FT new_value) {

				assert(new_value > FT(0));
				m_cluster_epsilon = new_value;
			}

			void set_normal_threshold(const FT new_value) {

				assert(new_value > FT(0) && new_value < FT(1));
				m_normal_threshold = new_value;
			}

			void set_minimum_shape_points(const size_t new_value) {

				assert(new_value > 0);
				m_min_points = new_value;
			}

			void save_info(const bool new_state) {
				m_save_info = new_state;
			}


			// Main.
			int detect(
				Boundary_data &   				  , Projected_points &boundary_clutter_projected,
				Boundary_data &building_boundaries, Projected_points &building_boundaries_projected,  
				const Container &input) {

				
				// Some preconditions.
				assert(!boundary_clutter_projected.empty());
				assert(building_boundaries.empty() && building_boundaries_projected.empty());

				assert(input.number_of_points() != 0);
				int number_of_detected_lines = -1;

				Log log;


				// (0) Create a kd tree.
				Fuzzy_tree tree;
				create_tree(tree, boundary_clutter_projected);


				// (1) Compute 2D bounding box.
				Point_2 bbmin, bbmax;
				compute_bounding_box(bbmin, bbmax, log, boundary_clutter_projected);


				// (2) First set all parameters to default if no user values are provided.
				const size_t total_points = boundary_clutter_projected.size();
				set_all_parameters(log, bbmin, bbmax, total_points);


				// (3) Estimate normals.
				Normals normals;
				estimate_normals(normals, boundary_clutter_projected, input);


				// (4) Sort all points by the quality of the locally fitted line.
				Sorted_indices sorted_indices;
				sort_projected_points(sorted_indices, log, boundary_clutter_projected, tree);


				// (5) Enter main loop with all sorted points and start growing local regions.
				grow_regions(building_boundaries, building_boundaries_projected, boundary_clutter_projected, log,
					normals, sorted_indices, tree);


				// Save log.
				if (m_save_info) log.save("tmp/region_growing_log");


				// Return number of detected lines.
				return number_of_detected_lines;
			}

		private:
			FT 	   m_epsilon;
			FT 	   m_cluster_epsilon;
			FT 	   m_normal_threshold;
			size_t m_min_points;
			bool   m_save_info;

			Simple_utils m_simple_utils;

			class Sort_by_linearity {

			public:
				Sort_by_linearity(const Projected_points &points, const Fuzzy_tree &tree, const FT cluster_epsilon) : 
				m_points(points),
				m_tree(tree),
				m_cluster_epsilon(cluster_epsilon), 
				m_default_score(-FT(1)),
				m_scores(m_points.size(), m_default_score) { }

				bool operator() (const int i, const int j) {
			        
					if (m_scores[i] == m_default_score) compute_score(i);
					if (m_scores[j] == m_default_score) compute_score(j);

			        return m_scores[i] > m_scores[j];
			    }

			private:
				const Projected_points &m_points;
				const Fuzzy_tree 	   &m_tree;

				const FT m_cluster_epsilon;
				const FT m_default_score;
				
				std::vector<FT> 	 m_scores;
				Projected_points     m_neighbours;
				std::vector<Point_2> m_local_points;

				void compute_score(const int point_index) {
					compute_nearest_neighbours(point_index);

					const size_t num_neighbours = m_neighbours.size();
					if (num_neighbours < 2) {

						m_scores[point_index] = FT(1000000);
						return;
					}

					map_neighbours_to_local_points();
					m_scores[point_index] = estimate_score();
			    }

			    void compute_nearest_neighbours(const int point_index) {

			    	const Point_2 &centre = m_points.at(point_index);
					const FT radius 	  = m_cluster_epsilon * FT(2);

					Fuzzy_circle circle(centre, radius, FT(0), m_tree.traits());

					m_neighbours.clear();
					m_tree.search(std::inserter(m_neighbours, m_neighbours.end()), circle);
			    }

			    void map_neighbours_to_local_points() {

			    	const size_t num_neighbours = m_neighbours.size();
			    	assert(num_neighbours != 0);

			    	m_local_points.clear();
			    	m_local_points.resize(num_neighbours);

					size_t count = 0;
			    	for (Point_iterator nit = m_neighbours.begin(); nit != m_neighbours.end(); ++nit, ++count)
						m_local_points[count] = (*nit).second;

					assert(count == num_neighbours);
			    }

			    FT estimate_score() {

			    	Line_2 stub;
			    	assert(m_local_points.size() >= 2);

			        return CGAL::linear_least_squares_fitting_2(m_local_points.begin(), m_local_points.end(), stub, CGAL::Dimension_tag<0>());			    	
			    }
			};


			// Functions.
			void create_tree(Fuzzy_tree &tree, const Projected_points &boundary_clutter_projected) {

				assert(!boundary_clutter_projected.empty());
				tree.insert(boundary_clutter_projected.begin(), boundary_clutter_projected.end());
			}

			void compute_bounding_box(Point_2 &bbmin, Point_2 &bbmax, Log &log, const Projected_points &boundary_clutter_projected) {

				m_simple_utils.compute_bounding_box_in_2d(bbmin, bbmax, boundary_clutter_projected);
				if (m_save_info) log.out << "Bounding box: " << bbmin << " -- " << bbmax << std::endl;
			}

			void set_all_parameters(Log &log, const Point_2 &bbmin, const Point_2 &bbmax, const size_t total_points) {

				const FT diagonal = m_simple_utils.compute_2d_bounding_box_diagonal(bbmin, bbmax);
				if (m_save_info) log.out << "Bounding box diagonal: " << diagonal << std::endl;

				estimate_epsilon(diagonal);
				estimate_cluster_epsilon(diagonal);
				estimate_normal_threshold();
				estimate_min_points(total_points);

				if (m_save_info) 
				log.out << "Params: " << m_epsilon 			<< ", " 
									  << m_cluster_epsilon  << ", "
									  << m_normal_threshold << ", "
									  << m_min_points 		<< std::endl;
			}

			void estimate_epsilon(const FT diagonal) {
				if (m_epsilon == -FT(1)) m_epsilon = diagonal / FT(100);
			}

			void estimate_cluster_epsilon(const FT diagonal) {
				if (m_cluster_epsilon == -FT(1)) m_cluster_epsilon = diagonal / FT(100);
			}

			void estimate_normal_threshold() {
				if (m_normal_threshold == -FT(1)) m_normal_threshold = FT(9) / FT(10); // about 25 degrees
			}

			void estimate_min_points(const size_t total_points) {

				if (m_min_points == 0) m_min_points = total_points / 100;
				if (m_min_points < 10) m_min_points = 10;

				assert(m_min_points <= total_points);
			}

			void estimate_normals(Normals &normals, const Projected_points &boundary_clutter_projected, const Container &input) {

				// Add here new methods to estimate normals later!

				m_simple_utils.estimate_normals_from_3d(normals, boundary_clutter_projected, input);
				assert(normals.size() == boundary_clutter_projected.size());

				Log log; 
				log.export_projected_points_with_normals_as_xyz("tmp/estimated_normals", boundary_clutter_projected, normals, "unused path");
			}

			void find_nearest_neighbours(Projected_points &neighbours, const Fuzzy_tree &tree, const Projected_point &query, const Fuzzy_circle &circle) const {

				neighbours.clear();
				tree.search(std::inserter(neighbours, neighbours.end()), circle);
				neighbours[query.first] = query.second;
			}

			void compute_circle(Fuzzy_circle &circle, const Point_2 &centre, const FT radius) const {

				assert(radius > FT(0));
				circle = Fuzzy_circle(centre, radius);
			}

			void sort_projected_points(Sorted_indices &sorted_indices, Log &log, const Projected_points &boundary_clutter_projected, const Fuzzy_tree &tree) {

				assert(!boundary_clutter_projected.empty());
				assert(tree.size() != 0);
				assert(m_cluster_epsilon > FT(0));

				fill_unsorted_indices(sorted_indices, log, boundary_clutter_projected);
				std::sort(sorted_indices.begin(), sorted_indices.end(), Sort_by_linearity(boundary_clutter_projected, tree, m_cluster_epsilon));

				// Save log.
				if (m_save_info) {

					log.out << "Sorted indices: ";
					for (size_t i = 0; i < sorted_indices.size(); ++i) log.out << sorted_indices[i] << " ";
					log.out << std::endl;
				}
			}

			void fill_unsorted_indices(Sorted_indices &sorted_indices, Log &log, const Projected_points &boundary_clutter_projected) {

				sorted_indices.clear();
				sorted_indices.resize(boundary_clutter_projected.size());

				size_t count = 0;
				for (Point_iterator pit = boundary_clutter_projected.begin(); pit != boundary_clutter_projected.end(); ++pit, ++count) {
					
					const Projected_point &projected = (*pit);
					sorted_indices[count] = projected.first;
				}
				assert(count == sorted_indices.size());

				// Save log.
				if (m_save_info) {

					log.out << "Unsorted indices: ";
					for (size_t i = 0; i < sorted_indices.size(); ++i) log.out << sorted_indices[i] << " ";
					log.out << std::endl;
				}
			}

			void grow_regions(Boundary_data &boundaries, Projected_points &boundaries_projected, Projected_points &points, Log &,
				const Normals &normals, const Sorted_indices &sorted_indices, const Fuzzy_tree &tree) {

				assert(boundaries.empty() && boundaries_projected.empty());
				assert(!points.empty());

				assert(normals.size() == points.size() && sorted_indices.size() == points.size());
				assert(tree.size() != 0);

				// int line_index = -1;
				
				// const size_t total_points = points.size();
				// size_t available_points   = total_points;

				// std::vector<int> shape_index(available_points, -1);
				// std::vector<int> index_container;
      			// std::vector<int> index_container_former_ring;
      			// std::set<int>    index_container_current_ring;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_REGION_GROWING_2_H