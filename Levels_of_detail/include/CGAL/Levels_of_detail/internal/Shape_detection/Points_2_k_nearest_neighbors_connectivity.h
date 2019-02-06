#ifndef CGAL_LEVELS_OF_DETAIL_POINTS_2_K_NEAREST_NEIGHBORS_CONNECTIVITY_H
#define CGAL_LEVELS_OF_DETAIL_POINTS_2_K_NEAREST_NEIGHBORS_CONNECTIVITY_H

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

  template<typename GeomTraits>
  class Points_2_k_nearest_neighbors_connectivity {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;

    using Index_to_point_map = 
    internal::Index_to_point_map<Point_2>;

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

    Points_2_k_nearest_neighbors_connectivity(
      const std::vector<Point_2>& points, 
      const FT search_size) :
    m_points(points),
    m_number_of_neighbors(static_cast<std::size_t>(
      CGAL::to_double(search_size))),
    m_index_to_point_map(m_points),
    m_distance(m_index_to_point_map),
    m_tree(
      boost::counting_iterator<std::size_t>(0),
      boost::counting_iterator<std::size_t>(m_points.size()),
      Splitter(),
      Search_traits(m_index_to_point_map)) { 

      CGAL_precondition(m_points.size() > 0);

      m_tree.build();
      CGAL_precondition(m_number_of_neighbors >= 0);
    }

    void get_neighbors(
      const std::size_t query_index, 
      std::vector<std::size_t>& neighbors) const {
      
      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_points.size());

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
    const std::vector<Point_2>& m_points;
    const std::size_t m_number_of_neighbors;
    const Index_to_point_map m_index_to_point_map;

    Distance m_distance;
    Tree m_tree;

  }; // Points_2_k_nearest_neighbors_connectivity

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_POINTS_2_K_NEAREST_NEIGHBORS_CONNECTIVITY_H
