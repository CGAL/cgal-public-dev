#ifndef CGAL_LEVEL_OF_DETAIL_KD_TREE_WITH_DATA_CREATOR_H
#define CGAL_LEVEL_OF_DETAIL_KD_TREE_WITH_DATA_CREATOR_H

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/property_map.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

namespace CGAL {

namespace Level_of_detail {

template<class GeomTraits,
         class ValueType = std::size_t,
         class PointMap = typename Pointer_property_map<typename GeomTraits::Point_2>::const_type>
class Kd_tree_with_data_creator {
  
public:
  using Kernel             = GeomTraits;

  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;

  using Neighbors = std::vector<ValueType>;
  
  using Search_traits_2 = CGAL::Search_traits_2<Kernel>;
            
  using Search_traits   = CGAL::Search_traits_adapter<ValueType, PointMap, Search_traits_2>;
  using Search_circle   = CGAL::Fuzzy_sphere<Search_traits>;
  using Search_tree     = CGAL::Kd_tree<Search_traits>;

  // using Splitter = Sliding_midpoint<Search_traits>;
  // using Distance = Distance_adapter<ValueType, PointMap, Euclidean_distance<Search_traits_2> >;
  // using Knn = CGAL::Orthogonal_k_neighbor_search<Search_traits, Distance, Splitter, Search_tree>;

  using Knn = CGAL::Orthogonal_k_neighbor_search<Search_traits>;
  using Distance = typename Knn::Distance;
  using Splitter = typename Knn::Splitter;
          
  Kd_tree_with_data_creator(const std::vector<Point_2>& input, const FT local_search_radius) : 
    m_tree(NULL),
    m_local_search_radius(local_search_radius),
    m_nb_neighbors(1) { 

    m_tree_point_map = make_property_map (input);
    create_tree_2(input);
  }
          
  Kd_tree_with_data_creator(const std::vector<Point_2>& input, const int nb_neighbors) : 
    m_tree(NULL),
    m_local_search_radius(-1.),
    m_nb_neighbors(nb_neighbors) { 

    m_tree_point_map = make_property_map (input);
    create_tree_2(input);
  }

  template <typename InputRange>
  Kd_tree_with_data_creator(const InputRange& input, PointMap point_map, const int nb_neighbors = 1) :
    m_tree(NULL),
    m_tree_point_map (point_map),
    m_local_search_radius(-1.),
    m_nb_neighbors(nb_neighbors) { 

    create_tree_2(input);
  }

  ~Kd_tree_with_data_creator()
  {
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
    Distance distance (m_tree_point_map);
    Knn search (*m_tree, query, m_nb_neighbors, 0, true, distance);
    for (typename Knn::iterator it = search.begin(); it != search.end(); ++ it)
      neighbors.push_back (it->first);
  }

  inline const PointMap& point_map() const {
    return m_tree_point_map;
  }

private:

  Search_tree                 *m_tree;
  PointMap              m_tree_point_map;

  FT m_local_search_radius;
  int m_nb_neighbors;

  void create_tree_2(const std::vector<Point_2>& input) {

    m_tree = new Search_tree(
      boost::counting_iterator<std::size_t>(0),
      boost::counting_iterator<std::size_t>(input.size()),
      Splitter(),
      Search_traits(m_tree_point_map));
  }


  template <typename InputRange>
  void create_tree_2(const InputRange& input) {

    std::vector<typename InputRange::const_iterator> input_iterators;
    input_iterators.reserve (input.size());
    for (typename InputRange::const_iterator it = input.begin();
         it != input.end(); ++ it)
      input_iterators.push_back (it);
      
    m_tree = new Search_tree(
      input_iterators.begin(),
      input_iterators.end(),
      Splitter(),
      Search_traits(m_tree_point_map));
  }
};

} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_KD_TREE_WITH_DATA_CREATOR_H
