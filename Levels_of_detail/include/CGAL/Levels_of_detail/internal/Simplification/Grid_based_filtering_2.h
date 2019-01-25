#ifndef CGAL_LEVELS_OF_DETAIL_GRID_BASED_FILTERING_2_H
#define CGAL_LEVELS_OF_DETAIL_GRID_BASED_FILTERING_2_H

// STL includes.
#include <map>
#include <utility>

// CGAL includes.
#include <CGAL/barycenter.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Grid_based_filtering_2 {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2; 

    using Cell_id = std::pair<long, long>;
    using Cell_data = std::vector<std::size_t>;
    using Grid = std::map<Cell_id, Cell_data>;

    Grid_based_filtering_2(const FT grid_cell_width) : 
    m_grid_cell_width(grid_cell_width) 
    { }

    void apply(std::vector<Point_2>& range) const {
      CGAL_precondition(range.size() > 0);
                  
      Grid grid;
      create_grid(range, grid);
      downsample(grid, range);
    }

  private:
    const FT m_grid_cell_width;
    
    struct Cell_data_to_point {
      
    public:
      typedef std::size_t argument_type;
      typedef std::pair<Point_2, FT> result_type;
      
      const std::vector<Point_2>& m_range;

      Cell_data_to_point(const std::vector<Point_2>& range) : 
      m_range(range) 
      { }

      result_type operator()(const argument_type idx) const {
        return std::make_pair(m_range[idx], FT(1));
      }
    };

    void create_grid(
      const std::vector<Point_2>& range, 
      Grid& grid) const {
                  
      Cell_id cell_id;
      grid.clear();

      for (std::size_t i = 0; i < range.size(); ++i) {
        
        get_cell_id(range[i], cell_id);
        grid[cell_id].push_back(i);
      }
    }

    void get_cell_id(const Point_2& point, Cell_id& cell_id) const {

      const long id_x = get_id_value(point.x());
      const long id_y = get_id_value(point.y());

      cell_id = std::make_pair(id_x, id_y);
    }

    long get_id_value(const FT value) const {

      CGAL_precondition(m_grid_cell_width > FT(0));
      const long id = static_cast<long>(
        CGAL::to_double(value / m_grid_cell_width));

      if (value >= FT(0)) return id;
      return id - 1;
    }

    void downsample(
      const Grid& grid, 
      std::vector<Point_2>& range) const {
          
      std::vector<Point_2> output;
      CGAL_precondition(grid.size() != 0);

      for (auto it = grid.begin(); it != grid.end(); ++it) {        
        
        const Cell_data& cell_data = it->second;
        output.push_back(get_cell_representative(range, cell_data));
      }

      CGAL_postcondition(output.size() > 0);
      range = output;
    }

    const Point_2& get_cell_representative(
      const std::vector<Point_2>& range, 
      const Cell_data& cell_data) const {

      const Point_2 barycenter = CGAL::barycenter(
        boost::make_transform_iterator(
          cell_data.begin(), Cell_data_to_point(range)),
        boost::make_transform_iterator(
          cell_data.end(), Cell_data_to_point(range)));
        
      FT min_distance = std::numeric_limits<FT>::max();
      std::size_t min_id = 0;
        
      for (std::size_t i = 0; i < cell_data.size(); ++i) {            
        const Point_2& point = range[cell_data[i]];

        const FT distance = internal::compute_distance_2(point, barycenter);
        if (distance < min_distance) {
              
          min_distance = distance;
          min_id = cell_data[i];
        }
      }
      return range[min_id];
    }

  }; // Grid_based_filtering_2

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_GRID_BASED_FILTERING_2_H
