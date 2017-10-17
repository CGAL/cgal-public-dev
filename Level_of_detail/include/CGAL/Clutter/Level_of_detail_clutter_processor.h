#ifndef CGAL_LEVEL_OF_DETAIL_CLUTTER_PROCESSOR_H
#define CGAL_LEVEL_OF_DETAIL_CLUTTER_PROCESSOR_H

// STL includes.
#include <map>
#include <vector>
#include <cassert>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/linear_least_squares_fitting_2.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class BoundaryData, class ProjectedPoints>
		class Level_of_detail_clutter_processor {

		public:
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
			typedef typename Neighbor_search::Tree 					    					 Tree;

			typedef typename Boundary_data::const_iterator    Boundary_iterator;
			typedef typename Projected_points::const_iterator Point_iterator;
			typedef typename Neighbor_search::iterator  	  Neighbour_iterator;

			// Extra.
			using Log = CGAL::LOD::Mylog;

			using Cell_id  = std::pair<int, int>;
			using Grid_map = std::map<Cell_id, std::vector<int> >;

			using Grid_cell_iterator = Grid_map::const_iterator;

			enum class Fitter_type { LINE };
			enum class New_point_type { CENTROID };

			Level_of_detail_clutter_processor() : 
			m_k(3), 
			m_fitter_type(Fitter_type::LINE), 
			m_grid_cell_length(FT(1) / FT(100)),
			m_new_point_type(New_point_type::CENTROID) { }

			void set_number_of_neighbours(const size_t new_value) {
				m_k = new_value;
			}

			void set_fitter_type(const Fitter_type new_fitter_type) {
				m_fitter_type = new_fitter_type;
			}

			void set_grid_cell_length(const FT new_cell_length) {
				m_grid_cell_length = new_cell_length;
			}

			void set_new_point_type(const New_point_type new_point_type) {
				m_new_point_type = new_point_type;
			}

			int process(Boundary_data &boundary_clutter, Projected_points &boundary_clutter_projected) const {

				auto number_of_removed_points = 0;

				apply_thinning(boundary_clutter, boundary_clutter_projected);
				number_of_removed_points = apply_grid_simplify(boundary_clutter, boundary_clutter_projected);

				return number_of_removed_points;
			}

			int apply_thinning(Boundary_data &boundary_clutter, Projected_points &boundary_clutter_projected) const {

				int number_of_processed_points = 0;
				if (boundary_clutter.empty() || boundary_clutter_projected.empty()) return number_of_processed_points;
				
				Tree tree;
				create_tree(tree, boundary_clutter_projected);

				Projected_points thinned_points;
				thin_all_points(thinned_points, tree, boundary_clutter_projected);

				assert(boundary_clutter_projected.size() == thinned_points.size());
				boundary_clutter_projected = thinned_points;

				number_of_processed_points = static_cast<int>(thinned_points.size());
				assert(number_of_processed_points >= 0);

				Log log; 
				log.export_projected_points_as_xyz("tmp/thinning", boundary_clutter_projected, "unused path");

				return number_of_processed_points;
			}

			int apply_grid_simplify(Boundary_data &boundary_clutter, Projected_points &boundary_clutter_projected) const {

				int number_of_removed_points = 0;
				if (boundary_clutter.empty() || boundary_clutter_projected.empty()) return number_of_removed_points;

				Grid_map grid_map;
				create_grid_map(grid_map, boundary_clutter_projected);

				Projected_points cleaned_points;
				simplify_clutter_using_grid(cleaned_points, grid_map, boundary_clutter_projected);

				number_of_removed_points = static_cast<int>(boundary_clutter_projected.size() - cleaned_points.size());
				boundary_clutter_projected = cleaned_points;

				assert(number_of_removed_points >= 0);
				if (number_of_removed_points != 0) set_new_boundary_data(boundary_clutter, boundary_clutter_projected);

				Log log; 
				log.export_projected_points_as_xyz("tmp/grid_simplified", boundary_clutter_projected, "unused path");

				return number_of_removed_points;
			}

		// Fields.
		private:
			size_t 		   m_k;
			Fitter_type    m_fitter_type;
			FT 			   m_grid_cell_length;
			New_point_type m_new_point_type;

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
				return static_cast<int>(value / m_grid_cell_length);
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
			void create_new_point(Point_2 &new_point, const Grid_cell &grid_cell, const Projected_points &) const {

				switch (m_new_point_type) {

					case New_point_type::CENTROID:
						create_new_centroid(new_point, grid_cell);
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

			void set_new_boundary_data(Boundary_data &boundary_clutter, const Projected_points &boundary_clutter_projected) const {

				Boundary_data new_data;
				assert(!boundary_clutter.empty() && !boundary_clutter_projected.empty());

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

		// Thinning.
		private:
			void create_tree(Tree &tree, const Projected_points &boundary_clutter_projected) const {

				assert(!boundary_clutter_projected.empty());
				tree.insert(boundary_clutter_projected.begin(), boundary_clutter_projected.end());
			}

			void thin_all_points(Projected_points &thinned_points, const Tree &tree, const Projected_points &boundary_clutter_projected) const {

				assert(!boundary_clutter_projected.empty());
				thinned_points.clear();

				for (Point_iterator pit = boundary_clutter_projected.begin(); pit != boundary_clutter_projected.end(); ++pit) {
					
					const Projected_point &projected = *pit;
					handle_projected_point(thinned_points, tree, projected);
				}
			}

			void handle_projected_point(Projected_points &thinned_points, const Tree &tree, const Projected_point &projected) const {

				std::vector<Point_2> neighbours;
				const Point_2 &query = projected.second;
				find_nearest_neighbours(neighbours, tree, query);

				switch (m_fitter_type) {

					case Fitter_type::LINE:
						handle_with_line(thinned_points, neighbours, projected);
						break;

					default:
						assert(!"Wrong fitter type!");
						break;
				}
			}

			void find_nearest_neighbours(std::vector<Point_2> &neighbours, const Tree &tree, const Point_2 &query) const {
				
				neighbours.clear();
				Neighbor_search search(tree, query, m_k);

				for (Neighbour_iterator nit = search.begin(); nit != search.end(); ++nit) neighbours.push_back(nit->first.second);
			}

			void handle_with_line(Projected_points &thinned_points, const std::vector<Point_2> &neighbours, const Projected_point &projected) const {

				Line_2 line;

				fit_line_to_points(line, neighbours);
				create_thin_point_with_line(thinned_points, line, projected);
			}

			void fit_line_to_points(Line_2 &line, const std::vector<Point_2> &points) const {

				CGAL::linear_least_squares_fitting_2(points.begin(), points.end(), line, CGAL::Dimension_tag<0>());
			}

			void create_thin_point_with_line(Projected_points &thinned_points, const Line_2 &line, const Projected_point &projected) const {

				thinned_points[projected.first] = project_onto_line(line, projected.second);
			}

			// My custom function to handle precision problems when projecting points.
			Point_2 project_onto_line(const Line_2 &line, const Point_2 &p) const {

				const auto a = line.point(0);
				const auto b = line.point(1);

				const auto projected = a + CGAL::scalar_product(p - a, b - a) / CGAL::scalar_product(b - a, b - a) * (b - a);
				
				if (std::isnan(projected.x()) || std::isnan(projected.y())) return line.projection(p);
				else return projected;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_CLUTTER_PROCESSOR_H