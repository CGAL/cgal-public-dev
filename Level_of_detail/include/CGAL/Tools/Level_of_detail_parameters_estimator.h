#ifndef CGAL_LEVEL_OF_DETAIL_PARAMETERS_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_PARAMETERS_ESTIMATOR_H

#if defined(WIN32) || defined(_WIN32) 
#define PS "\\"
#else 
#define PS "/"
#endif

// STL includes
#include <map>
#include <string>
#include <vector>
#include <cassert>
#include <sstream>
#include <iostream>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/Simple_cartesian.h>

// New CGAL includes.
#include <CGAL/Selector/Level_of_detail_selector.h>
#include <CGAL/Selector/Level_of_detail_selection_strategy.h>
#include <CGAL/Projector/Level_of_detail_projector.h>
#include <CGAL/Utils/Level_of_detail_utils_simple.h>
#include <CGAL/Mylog/Mylog.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class InputParameters>
		class Level_of_detail_parameters_estimator {

		public:
			typedef KernelTraits    Kernel;
			typedef InputContainer  Container;
			typedef InputParameters Parameters;

			typedef typename Kernel::FT 	 FT;
			typedef typename Kernel::Line_2  Line_2;
			typedef typename Kernel::Point_2 Point_2;
			typedef typename Kernel::Plane_3 Plane_3;

			using Index   = int;
			using Indices = std::vector<Index>;

			typedef CGAL::LOD::Level_of_detail_building_boundary_and_interior<Kernel, Container> Boundary_and_interior_strategy;
			typedef CGAL::LOD::Level_of_detail_selector<Kernel, Boundary_and_interior_strategy>  Boundary_and_interior_selector;

			typedef std::pair<int, Point_2> 					 									 Projected_point;
			typedef std::map<int, Point_2>           	   									  		 Projected_points;
			typedef CGAL::LOD::Level_of_detail_simple_projector<Kernel, Container, Projected_points> Projector;

			typedef typename CGAL::Second_of_pair_property_map<Projected_point> Point_map;

			typedef CGAL::Search_traits_2<Kernel>                       					 Search_traits_2;
			typedef CGAL::Search_traits_adapter<Projected_point, Point_map, Search_traits_2> Search_traits;
			typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> 					     Neighbor_search;
			typedef CGAL::Kd_tree<Search_traits>					                         Fuzzy_tree;

			typedef typename Projected_points::const_iterator Point_iterator;
			typedef typename Neighbor_search::iterator  	  Neighbour_iterator;

			typename Kernel::Compute_squared_distance_2 squared_distance;
			using Log = CGAL::LOD::Mylog;

			typedef Level_of_detail_utils_simple<Kernel> Simple_utils;

			enum class Estimation_method { FULL };
			enum class Estimated_parameter { SCALE, EPS };


			Level_of_detail_parameters_estimator(const Container &input, Parameters &parameters) :
			m_input(input), m_parameters(parameters),
			m_estimation_method(Estimation_method::FULL),
			m_num_knn(12), m_save_info(false) { }


			void estimate() {

				initialize_estimation_data();

				estimate_scale();
				estimate_eps();
			}

		private:
			const Container   &m_input;
			Parameters 	      &m_parameters;
			Estimation_method m_estimation_method;
			int 			  m_num_knn;
			bool 			  m_save_info;
			Projected_points  m_points;

			Simple_utils m_simple_utils;


			void initialize_estimation_data() {
				create_points();
			}

			void create_points() {

				Indices indices;
				Boundary_and_interior_selector selector;
				selector.select_elements(m_input, std::back_inserter(indices));

				const Plane_3 ground = Plane_3(FT(0), FT(0), FT(1), FT(0));
				m_points.clear();

				Projector projector;
				projector.project_with_indices(m_input, indices, ground, m_points);

				if (m_save_info) {
					Log exporter;
					exporter.export_projected_points_as_xyz("tmp" + std::string(PS) + "points_for_estimation", m_points, "stub");
				}
			}

			void estimate_scale() {
				switch (m_estimation_method) {
					
					case Estimation_method::FULL:
						estimate_scale_full();
						break;

					default:
						assert(!"Wrong scale estimation method");
						break;
				}
			}

			void estimate_eps() {
				switch (m_estimation_method) {
					
					case Estimation_method::FULL:
						estimate_eps_full();
						break;

					default:
						assert(!"Wrong eps estimation method");
						break;
				}
			}

			void estimate_scale_full() {
				
				FT scale = -FT(1);
				Estimated_parameter est_param = Estimated_parameter::SCALE;

				go_over_all_points_and_estimate_parameter(est_param, scale);
				assert(scale >= FT(0));

				const FT scale_factor = static_cast<FT>(m_num_knn);
				m_parameters["-scale"] = to_string(scale_factor * scale);
			}

			void estimate_eps_full() {
				
				FT eps = -FT(1);
				Estimated_parameter est_param = Estimated_parameter::EPS;

				go_over_all_points_and_estimate_parameter(est_param, eps);
				assert(eps >= FT(0));

				const FT eps_factor = static_cast<FT>(m_num_knn);
				m_parameters["-eps"] = to_string(eps_factor * eps);
			}

			std::string to_string(const FT value) {
				
				std::string result;
	        	std::stringstream ss;
	        
	        	ss << value;
	        	ss >> result;

	        	return result;
			}

			void go_over_all_points_and_estimate_parameter(const Estimated_parameter est_param, FT &param_value) {
				
				Fuzzy_tree tree;
				create_tree(tree);

				param_value = FT(0);
				Projected_points neighbours;

				for (Point_iterator pit = m_points.begin(); pit != m_points.end(); ++pit) {
					const Projected_point &query = *pit;
	
					find_nearest_neighbours(tree, neighbours, query);
					param_value += estimate_parameter(neighbours, query, est_param);
				}
				param_value /= static_cast<FT>(m_points.size());
			}

			void create_tree(Fuzzy_tree &tree) {
				tree.clear();

				assert(!m_points.empty());
				tree.insert(m_points.begin(), m_points.end());
			}

			void find_nearest_neighbours(Fuzzy_tree &tree, Projected_points &neighbours, const Projected_point &query) {
				
				Neighbor_search search(tree, query.second, m_num_knn);
				const size_t num_points = static_cast<size_t>(std::distance(search.begin(), search.end()));

				neighbours.clear();
				size_t count = 0;

				for (Neighbour_iterator nit = search.begin(); nit != search.end(); ++nit, ++count) neighbours[nit->first.first] = nit->first.second;
				assert(count == num_points && num_points > 0);

				if (count != num_points) {

					std::cerr << std::endl << "Error: count != num_points, find_nearest_neighbours function parameter estimator!" << std::endl << std::endl;
					exit(EXIT_FAILURE);
				}

				neighbours[query.first] = query.second;
			}

			FT estimate_parameter(const Projected_points &points, const Projected_point &query, const Estimated_parameter est_param) {
				switch (est_param) {

					case Estimated_parameter::SCALE:
						return estimate_scale_parameter(points, query);

					case Estimated_parameter::EPS:
						return estimate_eps_parameter(points);

					default:
						assert(!"Wrong estimated parameter!");
						return -FT(1);
				}
			}

			FT estimate_scale_parameter(const Projected_points &points, const Projected_point &query) {
				assert(points.size() > 0);

				FT scale_value = FT(0);
				FT count = FT(0);

				for (Point_iterator pit = points.begin(); pit != points.end(); ++pit) { 
					const Projected_point &point = (*pit);

					const FT squared_dist = squared_distance(query.second, point.second);
					scale_value += static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_dist)));

					if (squared_dist != FT(0)) count += FT(1);
				}
				scale_value /= count;

				return scale_value;
			}

			FT estimate_eps_parameter(const Projected_points &points) {
				assert(points.size() > 0);

				Line_2 line;
				fit_line_to_points(points, line);

				FT eps_value = FT(0);
				FT count = FT(0);

				for (Point_iterator pit = points.begin(); pit != points.end(); ++pit) {
					
					const Point_2 &point 		  = (*pit).second;
					const Point_2 projected_point = m_simple_utils.project_onto_line(line, point);

					const FT squared_dist = squared_distance(point, projected_point);
					eps_value += static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_dist)));

					if (squared_dist != FT(0)) count += FT(1);
				}
				eps_value /= count;

				return eps_value;
			}

			void fit_line_to_points(const Projected_points &points, Line_2 &line) {

      			using Local_Kernel = CGAL::Simple_cartesian<double>;
				using Point_2ft    = Local_Kernel::Point_2;
				using Line_2ft     = Local_Kernel::Line_2;

				const size_t num_points = points.size();
				std::vector<Point_2ft> tmp_points(num_points);

				size_t count = 0;
				for (Point_iterator pit = points.begin(); pit != points.end(); ++pit, ++count) {
					const Point_2 &point = (*pit).second;

					const double x = CGAL::to_double(point.x());
					const double y = CGAL::to_double(point.y());

					tmp_points[count] = Point_2ft(x, y);
				}
				assert(num_points == count);

				Line_2ft tmp_line;
				CGAL::linear_least_squares_fitting_2(tmp_points.begin(), tmp_points.end(), tmp_line, CGAL::Dimension_tag<0>());
				line = Line_2(static_cast<FT>(tmp_line.a()), static_cast<FT>(tmp_line.b()), static_cast<FT>(tmp_line.c()));
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_PARAMETERS_ESTIMATOR_H