#ifndef CGAL_LEVEL_OF_DETAIL_PARAMETERS_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_PARAMETERS_ESTIMATOR_H

#if defined(WIN32) || defined(_WIN32) 
#define PS "\\"
#define PN "\r\n"
#else 
#define PS "/" 
#define PN "\n"
#endif

// STL includes
#include <string>
#include <vector>
#include <cassert>
#include <map>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/Fuzzy_sphere.h>

// New CGAL includes.
#include <CGAL/Selector/Level_of_detail_selector.h>
#include <CGAL/Selector/Level_of_detail_selection_strategy.h>
#include <CGAL/Projector/Level_of_detail_projector.h>
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
			typedef typename Kernel::Point_2 Point_2;
			typedef typename Kernel::Plane_3 Plane_3;
			typedef typename Kernel::Line_2  Line_2;

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

			enum class Estimation_method { FULL };
			enum class Estimated_parameter { SCALE, EPS };


			Level_of_detail_parameters_estimator(const Container &input, Parameters &parameters) :
			m_input(input), m_parameters(parameters),
			m_estimation_method(Estimation_method::FULL),
			m_num_knn(12) { }


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
			Projected_points  m_points;


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
				go_over_all_points_and_compute_average(est_param, scale);

				m_parameters["-scale"] = to_string(scale);
			}

			void estimate_eps_full() {
				
				FT eps = -FT(1);
				Estimated_parameter est_param = Estimated_parameter::EPS;
				go_over_all_points_and_compute_average(est_param, eps);

				m_parameters["-eps"] = to_string(eps);
			}

			// FINISH!
			std::string to_string(const FT /* value */) {
				return "5.0";
			}

			void go_over_all_points_and_compute_average(const Estimated_parameter &est_param, FT &param_value) {
				
				Fuzzy_tree tree;
				create_tree(tree);

				param_value = FT(0);
				for (Point_iterator pit = m_points.begin(); pit != m_points.end(); ++pit) {
					const Projected_point &query = *pit;

					Projected_points neighbours;
					find_nearest_neighbours(tree, neighbours, query);

					param_value += estimate_parameter(neighbours, est_param);
				}
				param_value /= static_cast<FT>(m_points.size());
			}

			void create_tree(Fuzzy_tree &tree) {
				
				assert(!m_points.empty());
				tree.insert(m_points.begin(), m_points.end());
			}

			void find_nearest_neighbours(Fuzzy_tree &tree, Projected_points &neighbours, const Projected_point &query) {
				
				Neighbor_search search(tree, query.second, m_num_knn);
				const size_t num_points = static_cast<size_t>(std::distance(search.begin(), search.end()));

				neighbours.clear();
				size_t count = 0;

				for (Neighbour_iterator nit = search.begin(); nit != search.end(); ++nit, ++count) neighbours[nit->first.first] = nit->first.second;
				assert(count == num_points);

				if (count != num_points) {

					std::cerr << std::endl << "Error: count != num_points, find_nearest_neighbours function parameter estimator!" << std::endl << std::endl;
					exit(EXIT_FAILURE);
				}

				neighbours[query.first] = query.second;
			}

			FT estimate_parameter(const Projected_points &points, const Estimated_parameter &est_param) {
				switch (est_param) {

					case Estimated_parameter::SCALE:
						return estimate_scale_value(points);

					case Estimated_parameter::EPS:
						return estimate_eps_value(points);

					default:
						assert(!"Wrong estimated parameter!");
						return -FT(1);
				}
			}

			FT estimate_scale_value(const Projected_points & /* points */) {
				return -FT(1);
			}

			FT estimate_eps_value(const Projected_points & /* points */) {
				return -FT(1);
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_PARAMETERS_ESTIMATOR_H