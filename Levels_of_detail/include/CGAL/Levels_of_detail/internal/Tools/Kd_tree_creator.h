#ifndef CGAL_LEVELS_OF_DETAIL_KD_TREE_WITH_DATA_CREATOR_H
#define CGAL_LEVELS_OF_DETAIL_KD_TREE_WITH_DATA_CREATOR_H

// STL includes.
#include <vector>
#include <utility>

// Boost includes.
#include <boost/iterator/counting_iterator.hpp>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/property_map.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<class GeometricTraits,
         class ValueType = std::size_t,
         class PointMap = 
         typename Pointer_property_map<typename GeometricTraits::Point_2>::const_type>
class Kd_tree_creator {
  
public:
  using Traits = GeometricTraits;
  using Value_type = ValueType;
  using Point_map = PointMap;

  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Point_3 = typename Traits::Point_3;

  using Neighbors = std::vector<Value_type>;
  
  using Search_traits_2 = CGAL::Search_traits_2<Traits>;
  using Search_traits = CGAL::Search_traits_adapter<Value_type, Point_map, Search_traits_2>;
  using Search_circle = CGAL::Fuzzy_sphere<Search_traits>;
  using Search_tree = CGAL::Kd_tree<Search_traits>;

  using Knn = CGAL::Orthogonal_k_neighbor_search<Search_traits>;
  using Distance = typename Knn::Distance;
  using Splitter = typename Knn::Splitter;
          
  Kd_tree_creator(const std::vector<Point_2> &points, const FT local_search_radius) : 
    m_tree(NULL),
    m_local_search_radius(local_search_radius),
    m_nb_neighbors(1) { 

    m_tree_point_map = make_property_map(points);
    create_tree_2(points);
  }
          
  Kd_tree_creator(const std::vector<Point_2> &points, const int nb_neighbors) : 
    m_tree(NULL),
    m_local_search_radius(-FT(1)),
    m_nb_neighbors(nb_neighbors) { 

    m_tree_point_map = make_property_map(points);
    create_tree_2(points);
  }

  template<typename InputRange>
  Kd_tree_creator(
    const InputRange &input_range, 
    Point_map point_map, 
    const int nb_neighbors = 1) :
  m_tree(NULL),
  m_tree_point_map(point_map),
  m_local_search_radius(-FT(1)),
  m_nb_neighbors(nb_neighbors) { 

    create_tree_2(input_range);
  }

  ~Kd_tree_creator() {
    
    if (m_tree != NULL)
      delete m_tree;
  }

  void set_local_search_radius(const FT local_search_radius) {
                
    CGAL_precondition(local_search_radius > FT(0));
    m_local_search_radius = local_search_radius;
  }

  void search_2(const Point_2 &query, Neighbors &neighbors) const {
    neighbors.clear();

    Search_circle circle(query, m_local_search_radius, 0, m_tree->traits());
    m_tree->search(std::back_inserter(neighbors), circle);
  }

  void search_knn_2(const Point_2 &query, Neighbors &neighbors) const {
    
    neighbors.clear();
    Distance distance(m_tree_point_map);
    Knn search(*m_tree, query, m_nb_neighbors, 0, true, distance);
    
    for (auto it = search.begin(); it != search.end(); ++it)
      neighbors.push_back(it->first);
  }

  inline const Point_map &point_map() const {
    return m_tree_point_map;
  }

private:

  Search_tree *m_tree;
  Point_map m_tree_point_map;

  FT m_local_search_radius;
  int m_nb_neighbors;

  void create_tree_2(const std::vector<Point_2> &points) {

    m_tree = new Search_tree(
      boost::counting_iterator<std::size_t>(0),
      boost::counting_iterator<std::size_t>(points.size()),
      Splitter(),
      Search_traits(m_tree_point_map));
  }

  template<typename Input_range>
  void create_tree_2(const Input_range &input_range) {

    std::vector<typename Input_range::const_iterator> input_iterators;
    input_iterators.reserve(input_range.size());
    
    for (auto it = input_range.begin(); it != input_range.end(); ++it)
      input_iterators.push_back(it);
      
    m_tree = new Search_tree(
      input_iterators.begin(),
      input_iterators.end(),
      Splitter(),
      Search_traits(m_tree_point_map));
  }

}; // Kd_tree_creator

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_KD_TREE_WITH_DATA_CREATOR_H
