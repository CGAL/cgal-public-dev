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
#include <CGAL/Levels_of_detail/enumerations.h>
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

    void estimate(
      const FT min_radius, 
      std::vector<Tree>& trees) const {

      std::vector<Tree> tmp_trees;
      estimate_center_and_radius(tmp_trees);
      set_clean_trees(min_radius, tmp_trees, trees);
    }

    void compute_heights(
      const Extrusion_type extrusion_type,
      std::vector<Tree>& trees) const {

      switch (extrusion_type) {

        case Extrusion_type::MIN:
        compute_min_heights(trees);
        return;

        case Extrusion_type::AVERAGE:
        compute_average_heights(trees);
        return;

        case Extrusion_type::MAX:
        compute_max_heights(trees);
        return;

        default:
        compute_max_heights(trees);
        return;
      }
    }

  private:
    void estimate_center_and_radius(std::vector<Tree>& trees) const {

      trees.clear();
      trees.resize(m_clusters.size());

      for (std::size_t i = 0; i < m_clusters.size(); ++i) {
        const Iterators& cluster = m_clusters[i];

        const Point_3 center = CGAL::barycenter(
          boost::make_transform_iterator(
            cluster.begin(), Point_from_iterator_and_pmap(m_point_map)),
          boost::make_transform_iterator(
            cluster.end(), Point_from_iterator_and_pmap(m_point_map)));

        FT radius = FT(0);
        for (std::size_t j = 0; j < cluster.size(); ++j)
          radius += CGAL::squared_distance(
            center, get(m_point_map, *(cluster[j])));
        radius = static_cast<FT>(
          CGAL::sqrt(
            CGAL::to_double(
              radius / static_cast<FT>(cluster.size()))));

        trees[i].center = Point_2(center.x(), center.y());
        trees[i].radius = radius;

        trees[i].cluster_index = i;
      }
    }

    void set_clean_trees(
      const FT min_radius,
      std::vector<Tree> &trees,
      std::vector<Tree>& clean) const {

      clean.clear();

      std::sort(
        trees.begin(), trees.end(),
        [](const Tree& a, const Tree& b) -> bool {
          return a.radius > b.radius;
        });

      for (std::size_t i = 0; i < trees.size(); ++i) {
        const Tree& tree_a = trees[i];

        if (tree_a.radius < min_radius)
          continue;
        
        Circle_2 circle_a(tree_a.center, tree_a.radius * tree_a.radius);

        bool okay = true;
        for (std::size_t j = 0; j < clean.size(); ++j) {
          const Tree& tree_b = clean[j];

          Circle_2 circle_b(tree_b.center, tree_b.radius * tree_b.radius);

          if (CGAL::do_intersect(circle_a, circle_b) ||
            circle_b.has_on_bounded_side(circle_a.center())) {
            
            okay = false;
            break;
          }
        }

        if (okay)
          clean.push_back(tree_a);
      }
    }

    void compute_min_heights(std::vector<Tree>& trees) const {

      for (std::size_t i = 0; i < trees.size(); ++i) {
        const Iterators& cluster = m_clusters[trees[i].cluster_index];

        FT height = internal::max_value<FT>();
        for (std::size_t j = 0; j < cluster.size(); ++j)
          height = CGAL::min(
            height, get(m_point_map, *(cluster[j])).z());
        trees[i].height = height;
      }
    }

    void compute_average_heights(std::vector<Tree>& trees) const {

      for (std::size_t i = 0; i < trees.size(); ++i) {
        const Iterators& cluster = m_clusters[trees[i].cluster_index];

        FT height = FT(0);
        for (std::size_t j = 0; j < cluster.size(); ++j)
          height += get(m_point_map, *(cluster[j])).z();
        trees[i].height = height / static_cast<FT>(cluster.size());
      }
    }

    void compute_max_heights(std::vector<Tree>& trees) const {

      for (std::size_t i = 0; i < trees.size(); ++i) {
        const Iterators& cluster = m_clusters[trees[i].cluster_index];

        FT height = -internal::max_value<FT>();
        for (std::size_t j = 0; j < cluster.size(); ++j)
          height = CGAL::max(
            height, get(m_point_map, *(cluster[j])).z());
        trees[i].height = height;
      }
    }

  }; // Trees
    
} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_TREES_H
