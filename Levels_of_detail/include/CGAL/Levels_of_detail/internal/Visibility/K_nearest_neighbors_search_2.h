#ifndef CGAL_LEVELS_OF_DETAIL_K_NEAREST_NEIGHBORS_SEARCH_2_H
#define CGAL_LEVELS_OF_DETAIL_K_NEAREST_NEIGHBORS_SEARCH_2_H

// STL includes.
#include <vector>

// Boost includes.
#include <CGAL/boost/iterator/counting_iterator.hpp>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Splitters.h>
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

    using Point_2 = typename Traits::Point_2;
    
    using Index_to_point_map = 
    internal::Point_2_from_index_map<Traits, Input_range, Point_map>;

    using Search_base = 
    CGAL::Search_traits_2<Traits>;

    using Search_traits = 
    CGAL::Search_traits_adapter<std::size_t, Index_to_point_map, Search_base>;

    using Distance = 
    CGAL::Distance_adapter<
      std::size_t, 
      Index_to_point_map, 
      CGAL::Euclidean_distance<Search_base> >;

    using Splitter = 
    CGAL::Sliding_midpoint<Search_traits>;

    using Search_tree = 
    CGAL::Kd_tree<Search_traits, Splitter, CGAL::Tag_true>;

    using Neighbor_search = 
    CGAL::Orthogonal_k_neighbor_search<
      Search_traits, 
      Distance, 
      Splitter, 
      Search_tree>;

    using Tree = 
    typename Neighbor_search::Tree;

    K_nearest_neighbors_search_2(
      const Input_range& input_range,
      const Point_map point_map, 
      const std::size_t number_of_neighbors) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_number_of_neighbors(number_of_neighbors),
    m_index_to_point_map(m_input_range, m_point_map),
    m_distance(m_index_to_point_map),
    m_tree(
      boost::counting_iterator<std::size_t>(0),
      boost::counting_iterator<std::size_t>(m_input_range.size()),
      Splitter(),
      Search_traits(m_index_to_point_map)) { 

      CGAL_precondition(m_input_range.size() > 0);

      m_tree.build();
      CGAL_precondition(m_number_of_neighbors >= 0);
    }

    void get_neighbors(
      const Point_2& query_point, 
      std::vector<std::size_t>& neighbors) const {

      Neighbor_search neighbor_search(
        m_tree, 
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
      m_tree.clear();
    }

  private:

    // Fields.
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const std::size_t m_number_of_neighbors;
    const Index_to_point_map m_index_to_point_map;

    Distance m_distance;
    Tree m_tree;

  }; // K_nearest_neighbors_search_2

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_K_NEAREST_NEIGHBORS_SEARCH_2_H
