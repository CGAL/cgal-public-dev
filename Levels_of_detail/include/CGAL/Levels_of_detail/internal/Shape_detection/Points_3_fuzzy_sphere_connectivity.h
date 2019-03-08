#ifndef CGAL_LEVELS_OF_DETAIL_POINTS_3_FUZZY_SPHERE_CONNECTIVITY_H
#define CGAL_LEVELS_OF_DETAIL_POINTS_3_FUZZY_SPHERE_CONNECTIVITY_H

// STL includes.
#include <vector>

// Boost includes.
#include <CGAL/boost/iterator/counting_iterator.hpp>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Splitters.h>
#include <CGAL/assertions.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Points_3_fuzzy_sphere_connectivity {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;

    using Index_to_point_map = 
    internal::Point_3_from_index_map<Traits, Input_range, Point_map>;

    using Search_base = 
    CGAL::Search_traits_3<Traits>;
                    
    using Search_traits = 
    CGAL::Search_traits_adapter<std::size_t, Index_to_point_map, Search_base>;
      
    using Splitter = 
    CGAL::Sliding_midpoint<Search_traits>;
      
    using Fuzzy_sphere 
    = CGAL::Fuzzy_sphere<Search_traits>;
      
    using Tree 
    = CGAL::Kd_tree<Search_traits, Splitter, CGAL::Tag_true>;

    Points_3_fuzzy_sphere_connectivity(
      const Input_range& input_range,
      const Point_map point_map, 
      const FT search_size) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_search_radius(search_size),
    m_index_to_point_map(m_input_range, m_point_map),
    m_tree(
      boost::counting_iterator<std::size_t>(0),
      boost::counting_iterator<std::size_t>(m_input_range.size()),
      Splitter(),
      Search_traits(m_index_to_point_map)) { 

      CGAL_precondition(m_input_range.size() > 0);

      m_tree.build();
      CGAL_precondition(m_search_radius >= FT(0));
    }

    void get_neighbors(
      const std::size_t query_index, 
      std::vector<std::size_t>& neighbors) const {
      
      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_input_range.size());
      
      const Fuzzy_sphere sphere(
        query_index, 
        m_search_radius, 
        FT(0), 
        m_tree.traits());

      neighbors.clear();
      m_tree.search(std::back_inserter(neighbors), sphere);
    }

    void clear() {
      m_tree.clear();
    }

  private:

    // Fields.
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const FT m_search_radius;
    const Index_to_point_map m_index_to_point_map;

    Tree m_tree;

  }; // Points_3_fuzzy_sphere_connectivity

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_POINTS_3_FUZZY_SPHERE_CONNECTIVITY_H
