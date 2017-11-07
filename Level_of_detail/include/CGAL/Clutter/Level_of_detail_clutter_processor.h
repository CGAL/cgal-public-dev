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
#include <CGAL/Clutter/Level_of_detail_grid_simplify.h>
#include <CGAL/Clutter/Level_of_detail_thinning.h>

namespace CGAL {

	namespace LOD {

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