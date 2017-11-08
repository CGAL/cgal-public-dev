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

			typename Kernel::Compute_squared_distance_2 squared_distance;
			typename Kernel::Compute_scalar_product_3   dot_product;


			// Extra.
			using Log 	  = CGAL::LOD::Mylog;
			using Normals = std::map<int, Normal>;


			// Constructor.
			Level_of_detail_region_growing_2() { }


			// Main.
			int process(
				Boundary_data &boundary_clutter   , Projected_points &boundary_clutter_projected,
				Boundary_data &building_boundaries, Projected_points &building_boundaries_projected,  
				const Container &input) const {

				assert(!boundary_clutter.empty() && !boundary_clutter_projected.empty());
				assert(building_boundaries.empty() && building_boundaries_projected.empty());

				assert(input.number_of_points() != 0);
				int number_of_detected_lines = -1;

				// to be implemented!

				return number_of_detected_lines;
			}

		private:
			
			// Functions.

		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_REGION_GROWING_2_H