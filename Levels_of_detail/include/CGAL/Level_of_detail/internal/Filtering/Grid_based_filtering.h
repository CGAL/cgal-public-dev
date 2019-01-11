#ifndef CGAL_LEVEL_OF_DETAIL_GRID_BASED_FILTERING_H
#define CGAL_LEVEL_OF_DETAIL_GRID_BASED_FILTERING_H

// STL includes.
#include <map>
#include <utility>

// CGAL includes.
#include <CGAL/number_utils.h>
#include <CGAL/barycenter.h>

// LOD includes.
#include <CGAL/Level_of_detail/internal/utils.h>

namespace CGAL {

namespace Level_of_detail {

template<class GeomTraits>
class Grid_based_filtering {

public:
  using Kernel           = GeomTraits;

  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2; 

  using Cell_id   = std::pair<int, int>;
  using Cell_data = std::vector<std::size_t>;
  using Grid      = std::map<Cell_id, Cell_data>;
            
  using Grid_cell_iterator   = typename Grid::const_iterator;
  using Cell_data_iterator   = typename Cell_data::const_iterator;

  typename Kernel::Compute_squared_distance_2 squared_distance_2;

  Grid_based_filtering(const FT grid_cell_width) : 
    m_grid_cell_width(grid_cell_width),
    m_big_value(FT(100000000000000)) 
  { }

  void apply(const std::vector<Point_2>& range, std::vector<Point_2> &output) const {
    CGAL_precondition(range.size() > 0);
                
    Grid grid;
    output.clear();

    create_grid(range, grid);
    clean_points_using_grid(range, grid, output);
  }

private:
  
  struct Cell_data_to_point_unary_function
  {
    typedef std::size_t argument_type;
    typedef std::pair<Point_2, FT> result_type;
    
    const std::vector<Point_2>& range;

    Cell_data_to_point_unary_function (const std::vector<Point_2>& range) : range (range) { }

    result_type operator() (const argument_type& idx) const
    {
      return std::make_pair (range[idx], FT(1));
    }
    
  };

  const FT m_grid_cell_width;
  const FT m_big_value;

  void create_grid(const std::vector<Point_2>& range, Grid &grid) const {
                
    Cell_id cell_id;
    grid.clear();

    for (std::size_t i = 0; i < range.size(); ++ i) {
      get_cell_id(range[i], cell_id);
      grid[cell_id].push_back(i);
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

  void clean_points_using_grid(const std::vector<Point_2>& range, const Grid &grid, std::vector<Point_2> &output) const {
				
    CGAL_precondition(grid.size() != 0);
    for (Grid_cell_iterator gc_it = grid.begin(); gc_it != grid.end(); ++gc_it) {
					
      const Cell_data &cell_data = (*gc_it).second;
      Point_2 barycentre = CGAL::barycenter
        (boost::make_transform_iterator
         (cell_data.begin(), Cell_data_to_point_unary_function(range)),
         boost::make_transform_iterator
         (cell_data.end(), Cell_data_to_point_unary_function(range)));
      
      FT min_distance = std::numeric_limits<FT>::max();
      std::size_t min_id = 0;
      
      for (Cell_data_iterator cd_it = cell_data.begin(); cd_it != cell_data.end(); ++cd_it) {                    
                    
        const Point_2 &point = range[*cd_it];
        const FT distance    = compute_distance_2(point, barycentre);

        if (distance < min_distance) {
						
          min_distance = distance;
          min_id = *cd_it;
        }
      }
      output.push_back(range[min_id]);
    }		
  }

  inline FT compute_distance_2(const Point_2 &p, const Point_2 &q) const {
    return static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance_2(p, q))));
  }
};

} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_GRID_BASED_FILTERING_H
