#ifndef CGAL_LEVELS_OF_DETAIL_TREES_H
#define CGAL_LEVELS_OF_DETAIL_TREES_H

// STL includes.
#include <vector>
#include <utility>

// Boost includes.
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/iterator/transform_iterator.hpp>

// CGAL includes.
#include <CGAL/barycenter.h>
#include <CGAL/number_utils.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>
#include <CGAL/Levels_of_detail/internal/structures.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
typename GeomTraits, 
typename InputRange, 
typename PointMap>
class Trees {

public:
  using Traits = GeomTraits;
  using Input_range = InputRange;
  using Point_map = PointMap;

  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Point_3 = typename Traits::Point_3;
  using Circle_2 = typename Traits::Circle_2;

  using Iterator = typename Input_range::const_iterator;
  using Iterators = std::vector<Iterator>;
  using Tree = Tree<Traits>;
  
private:
  const std::vector<Iterators>& m_clusters;
  const Point_map m_point_map;

  struct Point_from_iterator_and_pmap {

  public: 
    using argument_type = Iterator;
    using return_type = std::pair<Point_3, FT>;

    Point_map m_point_map;

    Point_from_iterator_and_pmap(Point_map point_map) : 
    m_point_map(point_map) 
    { }

    return_type operator()(const argument_type& arg) const {
      return std::make_pair(get(m_point_map, *arg), FT(1));
    }
  };

public:
  Trees(
    const std::vector<Iterators>& clusters, 
    Point_map point_map) : 
  m_clusters(clusters), 
  m_point_map(point_map)
  { }

  void estimate(const FT min_radius, std::vector<Tree>& trees) {

    trees.clear();

    // Estimate.
    std::vector<Tree> tmp_trees(m_clusters.size());
    for (std::size_t i = 0; i < m_clusters.size(); ++i) {
      const Iterators& cluster = m_clusters[i];

      const Point_3 center = CGAL::barycenter(
        boost::make_transform_iterator(
          cluster.begin(), Point_from_iterator_and_pmap(m_point_map)),
        boost::make_transform_iterator(
          cluster.end(), Point_from_iterator_and_pmap(m_point_map)));

      FT radius = FT(0);
      FT height = -internal::max_value<FT>();
      
      for (std::size_t j = 0; j < cluster.size(); ++j) {
        radius += CGAL::squared_distance(
          center, get(m_point_map, *(cluster[j])));
        height = CGAL::max(
          height, get(m_point_map, *(cluster[j])).z());
      }

      radius = static_cast<FT>(
        CGAL::sqrt(
          CGAL::to_double(
            radius / static_cast<FT>(cluster.size()))));

      tmp_trees[i].center = Point_2(center.x(), center.y());
      tmp_trees[i].radius = radius;
      tmp_trees[i].height = height;
    }

    // Sort.
    std::sort(
      tmp_trees.begin(), tmp_trees.end(),
      [](const Tree& a, const Tree& b) -> bool {
        return a.radius > b.radius;
      });

    for (std::size_t i = 0; i < tmp_trees.size(); ++i) {
      const Tree& tmp_tree_a = tmp_trees[i];

      if (tmp_tree_a.radius < min_radius)
        continue;
      
      Circle_2 circle_a(
        tmp_tree_a.center, tmp_tree_a.radius * tmp_tree_a.radius);

      bool okay = true;
      for (std::size_t j = 0; j < trees.size(); ++j) {
        const Tree& tmp_tree_b = trees[j];

        Circle_2 circle_b(
          tmp_tree_b.center, tmp_tree_b.radius * tmp_tree_b.radius);

        if (CGAL::do_intersect(circle_a, circle_b) ||
          circle_b.has_on_bounded_side(circle_a.center())) {
          
          okay = false;
          break;
        }
      }

      if (okay)
        trees.push_back(tmp_tree_a);
    }
  }

}; // Trees
  
} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_TREES_H
