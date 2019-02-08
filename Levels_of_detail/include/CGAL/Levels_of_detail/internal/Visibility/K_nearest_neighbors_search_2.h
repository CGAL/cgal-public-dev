#ifndef CGAL_LEVELS_OF_DETAIL_K_NEAREST_NEIGHBORS_SEARCH_2_H
#define CGAL_LEVELS_OF_DETAIL_K_NEAREST_NEIGHBORS_SEARCH_2_H

// STL includes.
#include <vector>
#include <memory>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/assertions.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class K_nearest_neighbors_search_2 {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using Iterator = typename Input_range::const_iterator;
    using Point_2 = typename Traits::Point_2;
    
    using Point_2_from_iterator_map = 
    internal::Point_2_from_iterator_map<Iterator, Point_2, Point_map>;

    using Search_base = 
    CGAL::Search_traits_2<Traits>;

    using Search_traits = 
    CGAL::Search_traits_adapter<Iterator, Point_2_from_iterator_map, Search_base>;

    using Search_tree = 
    CGAL::Kd_tree<Search_traits>;

    using Neighbor_search = 
    CGAL::Orthogonal_k_neighbor_search<Search_traits>;

    using Splitter = 
    typename Neighbor_search::Splitter;

    using Distance = 
    typename Neighbor_search::Distance;

    K_nearest_neighbors_search_2(
      const Input_range& input_range,
      const Point_map point_map, 
      const std::size_t number_of_neighbors) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_number_of_neighbors(number_of_neighbors),
    m_distance(m_point_map) { 

      std::vector<Iterator> input_iterators;
      input_iterators.reserve(m_input_range.size());

      for (auto it = m_input_range.begin(); it != m_input_range.end(); ++it)
        input_iterators.push_back(it);
      
      m_tree = std::make_shared<Search_tree>(
        input_iterators.begin(),
        input_iterators.end(),
        Splitter(),
        Search_traits(m_point_map));
    }

    void get_neighbors(
      const Point_2& query_point, 
      std::vector<Iterator>& neighbors) const {

      Neighbor_search neighbor_search(
        *m_tree, 
        query_point, 
        m_number_of_neighbors, 
        0, 
        true, 
        m_distance);

      const std::size_t num_neighbors = 
      std::distance(neighbor_search.begin(), neighbor_search.end());

      neighbors.clear();
      neighbors.reserve(num_neighbors);

      for (auto it = neighbor_search.begin(); it != neighbor_search.end(); ++it)
        neighbors.push_back(it->first);
    }

    void clear() {
      m_tree->clear();
    }

  private:

    // Fields.
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const std::size_t m_number_of_neighbors;

    Distance m_distance;
    std::shared_ptr<Search_tree> m_tree;

  }; // K_nearest_neighbors_search_2

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_K_NEAREST_NEIGHBORS_SEARCH_2_H
