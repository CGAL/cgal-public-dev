#ifndef CGAL_LEVELS_OF_DETAIL_POINTS_FUZZY_SPHERE_CONNECTIVITY_H
#define CGAL_LEVELS_OF_DETAIL_POINTS_FUZZY_SPHERE_CONNECTIVITY_H

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
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits, 
  typename InputRange, 
  typename PointMap>
  class Points_fuzzy_sphere_connectivity {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using Point = typename Point_map::value_type;

    using FT = typename GeomTraits::FT;

    using Index_to_point_map = 
    internal::Item_property_map<Input_range, Point_map>;

    using Search_base = typename std::conditional<
      std::is_same<typename Traits::Point_2, Point>::value, 
      CGAL::Search_traits_2<Traits>, 
      CGAL::Search_traits_3<Traits> >::type;
                    
    using Search_traits = 
    CGAL::Search_traits_adapter<std::size_t, Index_to_point_map, Search_base>;
      
    using Splitter = 
    CGAL::Sliding_midpoint<Search_traits>;
      
    using Fuzzy_sphere 
    = CGAL::Fuzzy_sphere<Search_traits>;
      
    using Tree 
    = CGAL::Kd_tree<Search_traits, Splitter, CGAL::Tag_true>;

    Points_fuzzy_sphere_connectivity(
      const InputRange& input_range, 
      const FT search_radius = FT(1), 
      const PointMap point_map = PointMap()) :
    m_input_range(input_range),
    m_search_radius(search_radius),
    m_point_map(point_map),
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

    template<typename OutputIterator>
    void get_neighbors(
      const std::size_t query_index, 
      OutputIterator neighbors) const {
                
      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_input_range.size());
      
      const Fuzzy_sphere sphere(
        query_index, 
        m_search_radius, 
        FT(0), 
        m_tree.traits());

      m_tree.search(neighbors, sphere);
    }

    void clear() {
      m_tree.clear();
    }

  private:

    // Fields.
    const Input_range& m_input_range;
    
    const FT m_search_radius;

    const Point_map m_point_map;
    const Index_to_point_map m_index_to_point_map;

    Tree m_tree;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_POINTS_FUZZY_SPHERE_CONNECTIVITY_H
