#ifndef CGAL_LEVELS_OF_DETAIL_POINTS_3_K_NEAREST_NEIGHBORS_CONNECTIVITY_H
#define CGAL_LEVELS_OF_DETAIL_POINTS_3_K_NEAREST_NEIGHBORS_CONNECTIVITY_H

// STL includes.
#include <vector>

// Boost includes.
#include <CGAL/boost/iterator/counting_iterator.hpp>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Splitters.h>
#include <CGAL/assertions.h>
#include <CGAL/Search_traits_3.h>
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
  class Points_3_k_nearest_neighbors_connectivity {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;

    using Index_to_point_map = 
    internal::Point_3_from_index_map<Traits, Input_range, Point_map>;

    using Search_base = 
    CGAL::Search_traits_3<Traits>;

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

    Points_3_k_nearest_neighbors_connectivity(
      const Input_range& input_range,
      const Point_map point_map, 
      const FT search_size) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_number_of_neighbors(static_cast<std::size_t>(
      CGAL::to_double(search_size))),
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
      const std::size_t query_index, 
      std::vector<std::size_t>& neighbors) const {
      
      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_input_range.size());

      Neighbor_search neighbor_search(
        m_tree, 
        get(m_index_to_point_map, query_index), 
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

  }; // Points_3_k_nearest_neighbors_connectivity

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_POINTS_3_K_NEAREST_NEIGHBORS_CONNECTIVITY_H
