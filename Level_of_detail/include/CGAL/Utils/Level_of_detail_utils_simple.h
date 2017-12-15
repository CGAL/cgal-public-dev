#ifndef CGAL_LEVEL_OF_DETAIL_UTILS_SIMPLE_H
#define CGAL_LEVEL_OF_DETAIL_UTILS_SIMPLE_H

#if defined(WIN32) || defined(_WIN32) 
#define PS "\\" 
#else 
#define PS "/"
#endif 

// STL includes.
#include <map>
#include <cassert>
#include <vector>
#include <memory>
#include <string>

// CGAL includes.
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/utils.h>
#include <CGAL/Random.h>
#include <CGAL/number_utils.h>
#include <CGAL/constructions_d.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/property_map.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>

// Boost includes.
#include <boost/tuple/tuple.hpp>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits>
		class Level_of_detail_utils_simple {

		public:
			typedef KernelTraits   			  Kernel;
			typedef typename Kernel::FT 	  FT;
			typedef typename Kernel::Point_2  Point_2;
			typedef typename Kernel::Point_3  Point_3;
			typedef typename Kernel::Vector_2 Vector_2;
			typedef typename Kernel::Vector_3 Vector_3;
			typedef typename Kernel::Line_2   Line_2;

			using Projected_point = std::pair<int, Point_2>;
			using Point_map       = CGAL::Second_of_pair_property_map<Projected_point>;

			typedef CGAL::Search_traits_2<Kernel>                       					 Search_traits_2;
			typedef CGAL::Search_traits_adapter<Projected_point, Point_map, Search_traits_2> Search_traits;
			typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> 					     Neighbor_search;
			typedef CGAL::Fuzzy_sphere<Search_traits>                    					 Fuzzy_circle;
			typedef CGAL::Kd_tree<Search_traits>					                         Fuzzy_tree;

			typename Kernel::Compute_squared_distance_2 squared_distance;
			typename Kernel::Compute_scalar_product_3   dot_product;


			//////////////////////////////////
			// Bounding box in 2D.
			template<class Projected_points>
			void compute_bounding_box_in_2d(Point_2 &bbmin, Point_2 &bbmax, const Projected_points &points) const {

				using Point_iterator = typename Projected_points::const_iterator;
				const FT big_value   = FT(1000000000);

				FT minx =  big_value, miny =  big_value;
				FT maxx = -big_value, maxy = -big_value;

				for (Point_iterator pit = points.begin(); pit != points.end(); ++pit) {
					const Point_2 &point = (*pit).second;

					minx = CGAL::min(minx, point.x());
					miny = CGAL::min(miny, point.y());

					maxx = CGAL::max(maxx, point.x());
					maxy = CGAL::max(maxy, point.y());
				}

				bbmin = Point_2(minx, miny);
				bbmax = Point_2(maxx, maxy);
			}

			FT compute_2d_bounding_box_diagonal(const Point_2 &bbmin, const Point_2 &bbmax) const {

				return static_cast<FT>(CGAL::sqrt(CGAL::to_double(
          			(bbmax.x() - bbmin.x()) * (bbmax.x() - bbmin.x()) +
          			(bbmax.y() - bbmin.y()) * (bbmax.y() - bbmin.y())
          		)));
			}


			//////////////////////////////////
			// Normal estimation.

			template<class Normals, class Projected_points, class Container>
			void estimate_2d_normals_from_3d(Normals &normals, const Projected_points &points, const Container &input) const {

				using Point_iterator = typename Projected_points::const_iterator;

				normals.clear();
				assert(!points.empty() && input.number_of_points() != 0);

				// Project normals.
				for (Point_iterator pit = points.begin(); pit != points.end(); ++pit) {
					
					const Projected_point &projected = *pit;
					const int point_index = projected.first;
					
					project_normal(normals, point_index, input);
				}
			}

			template<class Normals, class Projected_points>
			void estimate_2d_normals_using_pca(Normals &normals, const Projected_points &points, const FT circle_radius) const {

				using Point_iterator = typename Projected_points::const_iterator;

				normals.clear();
				assert(!points.empty());

				// Create a tree.
				Fuzzy_tree tree;
				tree.insert(points.begin(), points.end());

				// Estimate normals.
				for (Point_iterator pit = points.begin(); pit != points.end(); ++pit) {
					
					const Projected_point &projected = *pit;
					estimate_local_normal<Normals, Projected_points>(normals, projected, circle_radius, tree);
				}
			}

			template<class Normals, class Projected_points>
			void estimate_local_normal(Normals &normals, const Projected_point &query, const FT circle_radius, const Fuzzy_tree &tree) const {
				
				using Point_iterator = typename Projected_points::const_iterator;

				// Find neighbours.
				Projected_points neighbours;
				Fuzzy_circle circle(query.second, circle_radius);

				tree.search(std::inserter(neighbours, neighbours.end()), circle);
				neighbours[query.first] = query.second;

				// Estimate normal.
				using Local_Kernel = CGAL::Simple_cartesian<double>;
				using Point_2ft    = Local_Kernel::Point_2;
				using Line_2ft     = Local_Kernel::Line_2;
				using Vector_2ft   = Local_Kernel::Vector_2;

				const size_t num_neighbours = neighbours.size();
				std::vector<Point_2ft> tmp_points(num_neighbours);

				size_t ind = 0;
				for (Point_iterator nit = neighbours.begin(); nit != neighbours.end(); ++nit, ++ind) {
					const Point_2 &p = (*nit).second;

					const double x = CGAL::to_double(p.x());
					const double y = CGAL::to_double(p.y());

					tmp_points[ind] = Point_2ft(x, y);
				}
				
				assert(num_neighbours   != 0);
				assert(ind == num_neighbours);

				Line_2ft tmp_line;
				CGAL::linear_least_squares_fitting_2(tmp_points.begin(), tmp_points.end(), tmp_line, CGAL::Dimension_tag<0>());
				
				const Vector_2ft vector = tmp_line.to_vector(); 						// Vector_2ft(tmp_line.point(0), tmp_line.point(1)); - this version works worse!
				Vector_2ft normal       = vector.perpendicular(CGAL::COUNTERCLOCKWISE);
				const auto length       = CGAL::sqrt(normal * normal);

				assert(length != 0.0);
				normal /= length;

				assert(std::isfinite(normal.x()) && std::isfinite(normal.y()));
				normals[query.first] = Vector_2(static_cast<FT>(normal.x()), static_cast<FT>(normal.y()));
			}
			

			//////////////////////////////////
			// Projection.

			template<class Normals, class Container>
			void project_normal(Normals &normals, const int point_index, const Container &input) const {

				typedef Vector_3 Normal;

				const Normal plane_normal  = Normal(FT(0), FT(0), FT(1));
				const Normal &point_normal = input.normal(point_index);

				const Normal projected_normal = point_normal - dot_product(point_normal, plane_normal) * plane_normal;
				assert(projected_normal.z() == FT(0));

				normals[point_index]  = Vector_2(projected_normal.x(), projected_normal.y());
				normals[point_index] /= static_cast<FT>(CGAL::sqrt(CGAL::to_double(normals[point_index] * normals[point_index])));
			}

			// My custom function to handle precision problems when projecting points.
			Point_2 project_onto_line(const Line_2 &line, const Point_2 &p) const {

				const auto a = line.point(0);
				const auto b = line.point(1);

				const auto projected = a + CGAL::scalar_product(p - a, b - a) / CGAL::scalar_product(b - a, b - a) * (b - a);
				
				if (std::isnan(CGAL::to_double(projected.x())) || std::isnan(CGAL::to_double(projected.y()))) return line.projection(p);
				else return projected;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_UTILS_SIMPLE_H