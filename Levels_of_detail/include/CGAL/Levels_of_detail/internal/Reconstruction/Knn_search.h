#ifndef CGAL_LEVELS_OF_DETAIL_KNN_SEARCH_H
#define CGAL_LEVELS_OF_DETAIL_KNN_SEARCH_H

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/property_map.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
typename Iterator, 
typename Point_2, 
typename PointMap>
struct Point_2_from_iterator_map {
  
  typedef Iterator key_type;
  typedef Point_2 value_type;
  typedef const value_type& reference;
  using category = boost::lvalue_property_map_tag;

  PointMap point_map;
  Point_2_from_iterator_map(PointMap point_map) : 
  point_map(point_map)
  { }

  friend reference get(
    const Point_2_from_iterator_map& pmap, 
    const key_type& k) {
    
    const typename PointMap::value_type& point_3 = get(pmap.point_map, *k);
    return reinterpret_cast<const Point_2&>(point_3);
  }
}; // Point_2_from_iterator_map

template<
typename GeomTraits,
typename ValueType,
typename PointMap>
class Knn_search {
  
public:
  using Traits = GeomTraits;

  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Point_3 = typename Traits::Point_3;

  using Neighbors = std::vector<ValueType>;
  using Search_traits_2 = CGAL::Search_traits_2<Traits>;
            
  using Search_traits = CGAL::Search_traits_adapter<ValueType, PointMap, Search_traits_2>;
  using Search_circle = CGAL::Fuzzy_sphere<Search_traits>;
  using Search_tree   = CGAL::Kd_tree<Search_traits>;

  using Knn = CGAL::Orthogonal_k_neighbor_search<Search_traits>;
  using Distance = typename Knn::Distance;
  using Splitter = typename Knn::Splitter;

  template<typename InputRange>
  Knn_search(
    const InputRange& input, 
    PointMap point_map, 
    const std::size_t nb_neighbors) :
  m_tree(NULL),
  m_tree_point_map(point_map),
  m_nb_neighbors(nb_neighbors) { 

    create_tree_2(input);
  }

  ~Knn_search() {
    if (m_tree != NULL)
      delete m_tree;
  }

  void get_neighbors(
    const Point_2& query, 
    Neighbors& neighbors) const {
    
    neighbors.clear();
    Distance distance(m_tree_point_map);
    Knn search(*m_tree, query, m_nb_neighbors, 0, true, distance);
    for (typename Knn::iterator it = search.begin(); it != search.end(); ++ it)
      neighbors.push_back (it->first);
  }

  inline const PointMap& point_map() const {
    return m_tree_point_map;
  }

private:
  Search_tree *m_tree;
  PointMap m_tree_point_map;
  const std::size_t m_nb_neighbors;

  template<typename InputRange>
  void create_tree_2(const InputRange& input) {

    std::vector<typename InputRange::const_iterator> input_iterators;
    input_iterators.reserve(input.size());
    for (typename InputRange::const_iterator it = input.begin();
    it != input.end(); ++it)
      input_iterators.push_back(it);
      
    m_tree = new Search_tree(
      input_iterators.begin(),
      input_iterators.end(),
      Splitter(),
      Search_traits(m_tree_point_map));
  }

}; // Knn_search

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_KNN_SEARCH_H
