#ifndef CGAL_LEVEL_OF_DETAIL_GRID_BASED_FILTERING_H
#define CGAL_LEVEL_OF_DETAIL_GRID_BASED_FILTERING_H

// STL includes.
#include <map>
#include <utility>

// CGAL includes.
#include <CGAL/number_utils.h>

// LOD includes.
#include <CGAL/Level_of_detail/Tools/Estimations/Barycentre_estimator.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

		template<class InputKernel, class PointIdentifier>
		class Grid_based_filtering {

        public:
            using Kernel           = InputKernel;
            using Point_identifier = PointIdentifier;

            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2; 

            using Cell_id   = std::pair<int, int>;
            using Cell_data = std::list<Point_identifier>;
			using Grid      = std::map<Cell_id, Cell_data>;
            
            using Grid_cell_iterator   = typename Grid::const_iterator;
            using Cell_data_iterator   = typename Cell_data::const_iterator;
            using Barycentre_estimator = LOD::Barycentre_estimator<Kernel>;

            typename Kernel::Compute_squared_distance_2 squared_distance_2;

            Grid_based_filtering(const FT grid_cell_width) : 
            m_grid_cell_width(grid_cell_width),
            m_big_value(FT(100000000000000)) 
            { }

            template<class Elements, class Point_map, class Output>
            void apply(const Elements &elements, const Point_map &point_map, Output &output) const {
                CGAL_precondition(elements.size() > 0);
                
                Grid grid;
                output.clear();

				create_grid(elements, point_map, grid);
				clean_points_using_grid(grid, point_map, output);
            }

        private:
            const FT m_grid_cell_width;
            const FT m_big_value;

            template<class Elements, class Point_map>
            void create_grid(const Elements &elements, const Point_map &point_map, Grid &grid) const {
                using Const_elements_iterator = typename Elements::const_iterator;
                
                Cell_id cell_id;
                grid.clear();

                for (Const_elements_iterator ce_it = elements.begin(); ce_it != elements.end(); ++ce_it) {                    
                    const Point_2 &point = get(point_map, *ce_it);

                    get_cell_id(point, cell_id);
                    grid[cell_id].push_back(*ce_it);
                }
			}

            void get_cell_id(const Point_2 &point, Cell_id &cell_id) const {

				const int id_x = get_id_value(point.x());
				const int id_y = get_id_value(point.y());

				cell_id = std::make_pair(id_x, id_y);
			}

			int get_id_value(const FT value) const {

				CGAL_precondition(m_grid_cell_width > FT(0));
				const int id = static_cast<int>(CGAL::to_double(value / m_grid_cell_width));

				if (value >= 0) return id;
				return id - 1;
			}

            template<class Point_map, class Output>
            void clean_points_using_grid(const Grid &grid, const Point_map &point_map, Output &output) const {
				
                CGAL_precondition(grid.size() != 0);
				for (Grid_cell_iterator gc_it = grid.begin(); gc_it != grid.end(); ++gc_it) {
					
					const Cell_data &cell_data = (*gc_it).second;
					find_cell_representative_point(cell_data, point_map, output);
				}		
			}

            template<class Point_map, class Output>
			void find_cell_representative_point(const Cell_data &cell_data, const Point_map &point_map, Output &output) const {

				Point_2 barycentre;
                Barycentre_estimator barycentre_estimator;
				
                barycentre_estimator.compute_barycentre_2(cell_data, point_map, barycentre);
                find_point(cell_data, point_map, barycentre, output);
			}

            template<class Point_map, class Output>
            void find_point(const Cell_data &cell_data, const Point_map &point_map, const Point_2 &barycentre, Output &output) const {

                FT min_distance = m_big_value; Point_identifier point_id;
                for (Cell_data_iterator cd_it = cell_data.begin(); cd_it != cell_data.end(); ++cd_it) {                    
                    
                    const Point_2 &point = get(point_map, *cd_it);
                    const FT distance    = compute_distance_2(point, barycentre);

					if (distance < min_distance) {
						
						min_distance = distance;
						point_id = *cd_it;
					}
                }
                output.push_back(point_id);
            }

            inline FT compute_distance_2(const Point_2 &p, const Point_2 &q) const {
                return static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance_2(p, q))));
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_GRID_BASED_FILTERING_H