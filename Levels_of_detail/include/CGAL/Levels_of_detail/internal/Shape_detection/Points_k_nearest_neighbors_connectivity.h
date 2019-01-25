#ifndef CGAL_LEVELS_OF_DETAIL_POINTS_K_NEAREST_NEIGHBORS_CONNECTIVITY_H
#define CGAL_LEVELS_OF_DETAIL_POINTS_K_NEAREST_NEIGHBORS_CONNECTIVITY_H

// STL includes.
#include <typeinfo>
#include <type_traits>

// Boost includes.
#include <CGAL/boost/iterator/counting_iterator.hpp>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Splitters.h>
#include <CGAL/assertions.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>

namespace CGAL {
namespace Levels_of_detail {

  template<
  typename GeomTraits, 
  typename InputRange, 
  typename PointMap>
  class Points_k_nearest_neighbors_connectivity {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using Point = typename Point_map::value_type;
    
    using Index_to_point_map = 
    internal::Item_property_map<Input_range, Point_map>;

    using Search_base = typename std::conditional<
      std::is_same<typename Traits::Point_2, Point>::value, 
      CGAL::Search_traits_2<Traits>, 
      CGAL::Search_traits_3<Traits> >::type;

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

    Points_k_nearest_neighbors_connectivity(
      const InputRange& input_range, 
      const std::size_t number_of_neighbors = 12, 
      const PointMap point_map = PointMap()) :
    m_input_range(input_range),
    m_number_of_neighbors(number_of_neighbors),
    m_point_map(point_map),
    m_index_to_point_map(m_input_range, m_point_map),
    m_distance(m_index_to_point_map),
    m_tree(
      boost::counting_iterator<std::size_t>(0),
      boost::counting_iterator<std::size_t>(m_input_range.size()),
      Splitter(),
      Search_traits(m_index_to_point_map)) { 

      m_tree.build();
      CGAL_precondition(number_of_neighbors >= 0);
    }

    template<typename OutputIterator>
    void get_neighbors(
      const std::size_t query_index, 
      OutputIterator neighbors) const {
      
      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_input_range.size());

      Neighbor_search neighbor_search(
        m_tree, 
        get(m_index_to_point_map, query_index), 
        m_number_of_neighbors, 
        0, 
        true, 
        m_distance);
                
      for (auto it = neighbor_search.begin(); it != neighbor_search.end(); ++it)
        *(neighbors++) = it->first;
    }

    void clear() {
      m_tree.clear();
    }

  private:

    // Fields.
    const Input_range& m_input_range;
    
    const std::size_t m_number_of_neighbors;

    const Point_map m_point_map;
    const Index_to_point_map m_index_to_point_map;

    Distance m_distance;
    Tree m_tree;
  };

} // namespace Levels_of_detail
} // namespace CGAL

#endif // CGAL_LEVELS_OF_DETAIL_POINTS_K_NEAREST_NEIGHBORS_CONNECTIVITY_H
