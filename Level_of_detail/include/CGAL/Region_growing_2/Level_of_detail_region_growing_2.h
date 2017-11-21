#ifndef CGAL_LEVEL_OF_DETAIL_REGION_GROWING_2_H
#define CGAL_LEVEL_OF_DETAIL_REGION_GROWING_2_H

// STL includes.
#include <map>
#include <vector>
#include <cassert>
#include <cmath>
#include <unordered_set>
#include <algorithm>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_2.h>

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
			typedef typename Kernel::Vector_2 Normal;

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
			using Index_iterator = std::vector<int>::const_iterator;
			using Set_iterator   = std::unordered_set<int>::const_iterator;


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
				sort_projected_points(sorted_indices, log, boundary_clutter_projected, tree, input.number_of_points());


				// (5) Enter main loop with all sorted points and start growing local regions.
				grow_regions(building_boundaries, building_boundaries_projected, boundary_clutter_projected, log,
					normals, sorted_indices, tree, input.number_of_points());


				// Save log.
				if (m_save_info) log.save("tmp/region_growing_log");

				log.clear();
				log.save_2d_region_growing("tmp/region_growing", building_boundaries, building_boundaries_projected, boundary_clutter_projected);


				// Return number of detected lines.
				number_of_detected_lines = building_boundaries.size();
				return number_of_detected_lines;
			}

		private:
			FT 	   m_epsilon;
			FT 	   m_cluster_epsilon;
			FT 	   m_normal_threshold;
			size_t m_min_points;
			bool   m_save_info;

			Simple_utils m_simple_utils;

			std::vector<int> 		m_index_container_former_ring;
      		std::unordered_set<int> m_index_container_current_ring;
      		Projected_points 		m_neighbours;
      		std::vector<Point_2> 	m_local_points;


			class Sort_by_linearity {

			public:
				Sort_by_linearity(const Projected_points &points, const Fuzzy_tree &tree, const FT cluster_epsilon, const size_t num_input_points) : 
				m_points(points),
				m_tree(tree),
				m_cluster_epsilon(cluster_epsilon), 
				m_default_score(-FT(1)),
				m_num_input_points(num_input_points),
				m_scores(new std::vector<FT>(m_num_input_points, m_default_score)) { 

					assert(m_tree.size() == m_points.size());
				}

				bool operator() (const int i, const int j) const {

					if ((*m_scores)[i] == m_default_score) compute_score(i);
					if ((*m_scores)[j] == m_default_score) compute_score(j);

					assert((*m_scores)[i] >= FT(0) && (*m_scores)[j] >= FT(0)); 
			        return (*m_scores)[i] > (*m_scores)[j];
			    }

			private:
				const Projected_points &m_points;
				const Fuzzy_tree 	   &m_tree;

				const FT m_cluster_epsilon;
				const FT m_default_score;

				const size_t m_num_input_points;
				mutable boost::shared_ptr< std::vector<FT> > m_scores;

				void compute_score(const int point_index) const {
					
					static Projected_points neighbours;
					compute_nearest_neighbours(point_index, neighbours);

					static std::vector<Point_2> local_points;
					map_neighbours_to_local_points(neighbours, local_points);

					estimate_score(local_points, point_index);
			    }

			    void compute_nearest_neighbours(const int point_index, Projected_points &neighbours) const {

			    	const Point_2 &centre = m_points.at(point_index);
					const FT radius 	  = m_cluster_epsilon * FT(2);

					Fuzzy_circle circle(centre, radius, FT(0), m_tree.traits());

					neighbours.clear();
					m_tree.search(std::inserter(neighbours, neighbours.end()), circle);
			    }

			    void map_neighbours_to_local_points(const Projected_points &neighbours, std::vector<Point_2> &local_points) const {

			    	const size_t num_neighbours = neighbours.size();
			    	assert(num_neighbours != 0);

			    	local_points.clear();
			    	local_points.resize(num_neighbours);

					size_t count = 0;
			    	for (Point_iterator nit = neighbours.begin(); nit != neighbours.end(); ++nit, ++count)
						local_points[count] = (*nit).second;

					assert(count == num_neighbours);
			    }

			    void estimate_score(const std::vector<Point_2> &local_points, const int point_index) const {
					
					Line_2 stub;
					assert(point_index >= 0 && point_index < static_cast<int>(m_num_input_points));

			        (*m_scores)[point_index] = CGAL::linear_least_squares_fitting_2(local_points.begin(), local_points.end(), stub, CGAL::Dimension_tag<0>());
			        assert((*m_scores)[point_index] >= FT(0));
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
				if (m_min_points < 2) m_min_points = 2;

				assert(m_min_points <= total_points);
			}

			void estimate_normals(Normals &normals, const Projected_points &boundary_clutter_projected, const Container &input) {

				// Add here new methods to estimate normals later!

				m_simple_utils.estimate_2d_normals_from_3d(normals, boundary_clutter_projected, input);
				assert(normals.size() == boundary_clutter_projected.size());

				Log log; 
				log.export_projected_points_with_normals_as_xyz("tmp/estimated_normals", boundary_clutter_projected, normals, "unused path");
			}

			void sort_projected_points(Sorted_indices &sorted_indices, Log &log, const Projected_points &boundary_clutter_projected, const Fuzzy_tree &tree, const size_t num_input_points) {

				assert(!boundary_clutter_projected.empty());
				assert(tree.size() != 0);
				assert(m_cluster_epsilon > FT(0));

				fill_unsorted_indices(sorted_indices, log, boundary_clutter_projected);

				assert(sorted_indices.size() == boundary_clutter_projected.size());
				std::stable_sort(sorted_indices.begin(), sorted_indices.end(), Sort_by_linearity(boundary_clutter_projected, tree, m_cluster_epsilon, num_input_points));

				// Save log.
				if (m_save_info) {

					log.out << "\nSorted indices (" << sorted_indices.size() << ") : ";
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

					log.out << "\nUnsorted indices (" << sorted_indices.size() << ") : ";
					for (size_t i = 0; i < sorted_indices.size(); ++i) log.out << sorted_indices[i] << " ";
					log.out << std::endl;
				}
			}

			void grow_regions(Boundary_data &boundaries, Projected_points &boundaries_projected, Projected_points &points, Log &log,
				const Normals &normals, const Sorted_indices &sorted_indices, const Fuzzy_tree &tree, const size_t num_input_points) {


				// Some preconditions.
				const size_t total_points = points.size();

				assert(boundaries.empty() && boundaries_projected.empty());
				assert(total_points != 0);

				assert(normals.size() == total_points && sorted_indices.size() == total_points);
				assert(tree.size() != 0);
				

				// Main structures.
				size_t available_points = total_points;

				std::vector<int> shape_index(num_input_points, -1);
				std::vector<int> index_container;

				Log internal;
				if (m_save_info) internal.out << "\nMain loop\n" <<std::endl;


				// Main loop.
				int class_index = -1;
      			for (size_t i = 0; i < sorted_indices.size(); ++i) {
      				if (m_save_info) internal.out<< "\nSorted index: " << i << std::endl << std::endl;
      				
      				const int point_index = sorted_indices[i];
      				if (shape_index[point_index] >= 0) continue;


      				// Get query point and its normal.
      				const Point_2 &query_pos = points.at(point_index);
      				const Normal &query_norm = normals.at(point_index);


      				// Get optimal line and its normal.
      				Line_2 optimal_line; Normal line_normal;
      				get_line_with_normal(query_pos, query_norm, optimal_line, line_normal);


      				// Initialize containers.
					shape_index[point_index] = ++class_index;

      				index_container.clear();
        			index_container.push_back(point_index);


        			// Propagate region through all neighbouring points.
        			propagate(shape_index, index_container, optimal_line, line_normal, internal, 
        				class_index, points, normals, point_index, tree);


        			// Update input data.
      				if (index_container.size() >= m_min_points) {

      					update_input_data(boundaries, boundaries_projected, class_index, index_container, points);
      					available_points -= index_container.size();

      				} else {

      					--class_index;
          				shape_index[point_index] = -1;

          				for (Index_iterator it = index_container.begin(); it != index_container.end(); ++it)
            				shape_index[*it] = -1;
      				}
      			}


      			// Update all unassigned points.
      			update_unassigned_points(points, shape_index);
      			assert(points.size() == available_points);


      			// Save internal log.
				// if (m_save_info) internal.save("tmp/internal_rg");


				// Save found shape indices.
				if (m_save_info) {

					log.out << "\nShape indices (" << shape_index.size() << ") :";
					for (size_t i = 0; i < shape_index.size(); ++i) log.out << shape_index[i] << " ";
					log.out << std::endl;
				}
			}

			void propagate(std::vector<int> &shape_index, std::vector<int> &index_container, Line_2 &optimal_line, Normal &line_normal, Log &log, 
				const int class_index, const Projected_points &points, const Normals &normals, const int point_index, const Fuzzy_tree &tree) {


				// Some preconditions.
				assert(m_cluster_epsilon > FT(0));
				assert(index_container.size() == 1);
				assert(point_index >= 0 && point_index < static_cast<int>(shape_index.size()));


				// Initialize containers.
				m_index_container_former_ring.clear();
        		m_index_container_former_ring.push_back(point_index);

        		m_index_container_current_ring.clear();


        		// Main loop.
        		if (m_save_info) {
        			log.out << "Propagation:\n" << std::endl;
					log.out << "line: " << optimal_line.point(0) << " -- " << optimal_line.point(1) 
							<< "; norm: " << line_normal << " with neighbours: " << std::endl;
				}

        		bool propagation = true;
        		do {

        			propagation = false;
        			for (Index_iterator fit = m_index_container_former_ring.begin(); fit != m_index_container_former_ring.end(); ++fit) {
        				const int former_index = (*fit);

        				const Point_2 &former_point = points.at(former_index);
        				const FT radius = m_cluster_epsilon;


        				// Compute circle.
        				Fuzzy_circle circle;
        				compute_circle(circle, former_point, radius);


        				// Find nearest neighbours.
        				find_nearest_neighbours(tree, circle);


        				// Handle neighbours.
        				handle_neighbours(propagation, shape_index, log, class_index, optimal_line, line_normal, points, normals);
        			}


        			// Update containers.
        			update_containers(index_container);


        			// Refit the line.
					if (index_container.size() < 2) continue;
					refit_the_line(optimal_line, line_normal, index_container, points);


        		} while (propagation);
			}

			void handle_neighbours(bool &propagation, std::vector<int> &shape_index, Log &log,
				const int class_index, const Line_2 &optimal_line, const Normal &line_normal, const Projected_points &points, const Normals &normals) {

				assert(!m_neighbours.empty());
				for (Point_iterator nit = m_neighbours.begin(); nit != m_neighbours.end(); ++nit) {
					const int neighbour_index = (*nit).first;

					assert(neighbour_index >= 0 && neighbour_index < static_cast<int>(shape_index.size()));
					if (shape_index[neighbour_index] >= 0) continue;


					// Get neighbour's position and normal.
					const Point_2 &neighbour_pos = points.at(neighbour_index);
					const Normal &neighbour_norm = normals.at(neighbour_index);


					// Compute distance to the line and cos.
					const FT eps = FT(1) / FT(1000000);
					assert(CGAL::abs(neighbour_norm.squared_length() - FT(1)) < eps && 
						   CGAL::abs(line_normal.squared_length() - FT(1)) < eps);

					const FT squared_distance = compute_squared_distance_to_the_line(neighbour_pos, optimal_line);
					const FT cos_angle 		  = CGAL::abs(neighbour_norm * line_normal);

					if (m_save_info) log.out << "\npos: "     << neighbour_pos 
											 << ", norm: "    << neighbour_norm 
											 << ", sq dist: " << squared_distance 
											 << ", cos: "     << cos_angle;


					// Check region growing conditions.
					assert(m_epsilon > FT(0));
					assert(m_normal_threshold > FT(0) && m_normal_threshold < FT(1));

					if (squared_distance > m_epsilon * m_epsilon || cos_angle < m_normal_threshold) continue;
					if (m_save_info) log.out << ", passed" << std::endl;


					// If all is good, add point to the shape.
					shape_index[neighbour_index] = class_index;
              		m_index_container_current_ring.insert(neighbour_index);
              		propagation = true;
				}
			}

			void update_containers(std::vector<int> &index_container) {

				m_index_container_former_ring.clear();
          		m_index_container_former_ring.resize(m_index_container_current_ring.size());

          		index_container.reserve(index_container.size() + m_index_container_current_ring.size());

          		size_t count = 0;
          		for (Set_iterator cit = m_index_container_current_ring.begin(); cit != m_index_container_current_ring.end(); ++cit, ++count) {
	            	
	            	m_index_container_former_ring[count] = *cit;
	            	index_container.push_back(*cit);
	          	}

	          	assert(count == m_index_container_current_ring.size());
          		m_index_container_current_ring.clear();
			}

			void refit_the_line(Line_2 &optimal_line, Normal &line_normal, const std::vector<int> &index_container, const Projected_points &points) {

				map_indices_to_local_points(index_container, points);

				Line_2 adjusted_line;
			    assert(m_local_points.size() >= 2);
			    CGAL::linear_least_squares_fitting_2(m_local_points.begin(), m_local_points.end(), adjusted_line, CGAL::Dimension_tag<0>());

			    optimal_line = adjusted_line;
			    get_line_normal(line_normal, optimal_line);
			}

			void map_indices_to_local_points(const std::vector<int> &index_container, const Projected_points &points) {

			    const size_t num_points = index_container.size();
			    assert(num_points >= 2);

			    m_local_points.clear();
			    m_local_points.resize(num_points);

				size_t count = 0;
			    for (Index_iterator it = index_container.begin(); it != index_container.end(); ++it, ++count)
					m_local_points[count] = points.at(*it);

				assert(count == num_points);
			}

			void get_line_with_normal(const Point_2 &p, const Normal &n, Line_2 &line, Normal &line_normal) {
				
				line = (Line_2(p, n)).perpendicular(p);
				line_normal = n / CGAL::sqrt(n * n);
			}

			FT compute_squared_distance_to_the_line(const Point_2 &p, const Line_2 &line) {
				return CGAL::squared_distance(p, line);
			}

			void get_line_normal(Normal &n, const Line_2 &line) {
				
				n = (line.perpendicular(line.point(0))).to_vector();
				n /= CGAL::sqrt(n * n);
			}

			void find_nearest_neighbours(const Fuzzy_tree &tree, const Fuzzy_circle &circle) {

				m_neighbours.clear();
				tree.search(std::inserter(m_neighbours, m_neighbours.end()), circle);
			}

			void compute_circle(Fuzzy_circle &circle, const Point_2 &centre, const FT radius) {

				assert(radius > FT(0));
				circle = Fuzzy_circle(centre, radius);
			}

			void update_input_data(Boundary_data &boundaries, Projected_points &boundaries_projected, 
				const int class_index, const std::vector<int> &index_container, const Projected_points &points) {

				boundaries[class_index] = index_container;
				for (size_t i = 0; i < index_container.size(); ++i) boundaries_projected[index_container[i]] = points.at(index_container[i]);
			}

			void update_unassigned_points(Projected_points &points, const std::vector<int> &shape_index) {
				
				assert(shape_index.size() >= points.size());
				Projected_points updated;

				for (Point_iterator pit = points.begin(); pit != points.end(); ++pit) {
					const int point_index = (*pit).first;

					if (shape_index[point_index] < 0) {

						const Point_2 &point_pos = (*pit).second;
						updated[point_index] = point_pos;
					}
				}
				points = updated;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_REGION_GROWING_2_H