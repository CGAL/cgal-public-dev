#ifndef CGAL_LEVELS_OF_DETAIL_POINTS_K_NEAREST_NEIGHBORS_CONNECTIVITY_H
#define CGAL_LEVELS_OF_DETAIL_POINTS_K_NEAREST_NEIGHBORS_CONNECTIVITY_H

// STL includes.
#include <vector>
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

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits, 
  typename PointType>
  class Points_k_nearest_neighbors_connectivity {

  public:
    using Traits = GeomTraits;
    using Point = PointType;

    using FT = typename Traits::FT;
    
    using Points = std::vector<Point>;
    struct Index_to_point_map {
                        
    public: 
      using value_type = Point;
      using reference = const value_type&;
      using key_type = std::size_t;
      using category = boost::lvalue_property_map_tag;

      Index_to_point_map(const Points& points) : 
      m_points(points) { 

        CGAL_precondition(m_points.size() > 0);
      }

      reference operator[](key_type index) const { 
                
        CGAL_precondition(index >= 0);
        CGAL_precondition(index < m_points.size());
        return m_points[index];
      }

      friend inline reference get(
        const Index_to_point_map& index_to_point_map, 
        key_type key) { 
        return index_to_point_map[key];
      }
                
    private:
      const Points& m_points;
    };

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
      const Points& points, 
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
      neighbors.resize(num_neighbors);

      std::size_t i = 0;
      for (auto it = neighbor_search.begin(); 
      it != neighbor_search.end(); ++it, ++i)
        neighbors[i] = it->first;
    }

    void clear() {
      m_tree.clear();
    }

  private:

    // Fields.
    const Points& m_points;
    const std::size_t m_number_of_neighbors;
    const Index_to_point_map m_index_to_point_map;

    Distance m_distance;
    Tree m_tree;

  }; // Points_k_nearest_neighbors_connectivity

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_POINTS_K_NEAREST_NEIGHBORS_CONNECTIVITY_H
