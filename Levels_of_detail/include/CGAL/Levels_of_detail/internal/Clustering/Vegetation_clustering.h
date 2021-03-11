#ifndef CGAL_LEVEL_OF_DETAIL_VEGETATION_CLUSTERING_H
#define CGAL_LEVEL_OF_DETAIL_VEGETATION_CLUSTERING_H

// STL includes.
#include <set>
#include <vector>
#include <utility>
#include <algorithm>

// Boost includes.
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/iterator/transform_iterator.hpp>

// CGAL includes.
#include <CGAL/bounding_box.h>
#include <CGAL/property_map.h>
#include <CGAL/Classification/Image.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Vegetation_clustering {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;

    using Iterator = typename Input_range::const_iterator;
    using Iterators = std::vector<Iterator>;
    using Image = CGAL::Classification::Image<Iterators>;
    using Coordinate = std::pair<std::size_t, std::size_t>;

  private:
    class Persistent_component {

    public:
      Persistent_component(
        const Coordinate& coordinate,
        FT height) :
      m_inliers(1, coordinate),
      m_start_height(height)
      { }

      std::vector<Coordinate>& inliers() {
        return m_inliers;
      }
      void add(const Coordinate& coordinate) {
        m_inliers.push_back(coordinate);
      }
      const FT& start_height() const {
        return m_start_height;
      }
      FT& end_height() {
        return m_end_height;
      }
      FT size() const {
        return m_start_height - m_end_height;
      }

    private:
      std::vector<Coordinate> m_inliers;
      FT m_start_height;
      FT m_end_height;
    };

    using Persistent_component_ptr = boost::shared_ptr<Persistent_component>;
    using Tree_image = CGAL::Classification::Image<Persistent_component_ptr>;

    const Input_range& m_input_range;
    const Point_map m_point_map;

  public:
    Vegetation_clustering(
      const Input_range& input_range,
      const Point_map& point_map) :
    m_input_range(input_range),
    m_point_map(point_map)
    { }

    void detect(
      const FT grid_cell_width_2, const FT min_height,
      std::vector<Iterators>& clusters) {
      clusters.clear();

      // Compute bounding box.
      CGAL::Bbox_3 bbox =
      CGAL::bbox_3(
        boost::make_transform_iterator(
          m_input_range.begin(),
          Property_map_to_unary_function<Point_map>(m_point_map)),
        boost::make_transform_iterator(
          m_input_range.end(),
          Property_map_to_unary_function<Point_map>(m_point_map)));

      // Initialize image.
      const std::size_t width =
      static_cast<std::size_t>(
        CGAL::to_double((bbox.xmax() - bbox.xmin()) / grid_cell_width_2)) + 1;
      const std::size_t height =
      static_cast<std::size_t>(
        CGAL::to_double((bbox.ymax() - bbox.ymin()) / grid_cell_width_2)) + 1;
      Image image(width, height);

      // Fill in this image.
      for (auto it = m_input_range.begin(); it != m_input_range.end(); ++it) {
        const Point_3& p = get(m_point_map, *it);

        const std::size_t x =
        static_cast<std::size_t>(
          CGAL::to_double((p.x() - bbox.xmin()) / grid_cell_width_2));
        const std::size_t y =
        static_cast<std::size_t>(
          CGAL::to_double((p.y() - bbox.ymin()) / grid_cell_width_2));

        image(x, y).push_back(it);
      }

      // Sort image cells by height.
      std::vector<Coordinate> sorted;
      for (std::size_t x = 0; x < image.width(); ++x) {
        for (std::size_t y = 0; y < image.height(); ++y) {
          if (!image(x, y).empty()) {

            std::sort(
              image(x, y).begin(),
              image(x, y).end(),
              [&](const Iterator& a, const Iterator& b) -> bool {
                return (get(m_point_map, *a).z() > get(m_point_map, *b).z());
              });

            sorted.push_back(std::make_pair(x, y));
          }
        }
      }

      std::sort(
        sorted.begin(),
        sorted.end(),
        [&](const Coordinate& a, const Coordinate& b) -> bool {
          return (
            get(m_point_map, *(image(a.first, a.second).front())).z() >
            get(m_point_map, *(image(b.first, b.second).front())).z());
        });

      // Initialize a tree map.
      Tree_image tree_map(width, height);
      for (std::size_t x = 0; x < image.width(); ++x)
        for (std::size_t y = 0; y < image.height(); ++y)
          tree_map(x, y) = Persistent_component_ptr();

      std::vector<Persistent_component_ptr> components;

      // Pick cells one by one and track life time of each component.
      for (std::size_t i = 0; i < sorted.size(); ++i) {
        const std::size_t x = sorted[i].first;
        const std::size_t y = sorted[i].second;

        if (tree_map(x, y) != Persistent_component_ptr()) // already handled
          continue;

        Iterators& cell = image(x, y);
        std::set<Persistent_component_ptr> local_trees;

        int xmin = int(x) - 1; if (xmin < 0) xmin = 0;
        int xmax = int(x) + 1; if (xmax == image.width()) xmax = image.width() - 1;
        int ymin = int(y) - 1; if (ymin < 0) ymin = 0;
        int ymax = int(y) + 1; if (ymax == image.height()) ymax = image.height() - 1;

        for (int xx = xmin; xx <= xmax; ++xx) {
          for (int yy = ymin; yy <= ymax; ++yy) {

            // The same cell.
            if (xx == x && yy == y)
              continue;

            Iterators& neighbor = image(xx, yy);
            if (neighbor.empty())
              continue;

            // Not unused cell.
            if (tree_map(xx, yy) != Persistent_component_ptr())
              local_trees.insert(tree_map(xx, yy));
          }
        }

        if (local_trees.empty()) {

          components.push_back(
            boost::make_shared<Persistent_component>(
              sorted[i], get(m_point_map, *(image(x, y).front())).z()));
          tree_map(x, y) = components.back();

        } else if (local_trees.size() == 1) {

          Persistent_component_ptr local_tree = *(local_trees.begin());
          tree_map(x, y) = local_tree;
          local_tree->add(sorted[i]);

        } else { // merge happens

          const FT height = get(m_point_map, *(image(x, y).front())).z();
          Persistent_component_ptr chosen;
          FT size_max = -internal::max_value<FT>();

          // Keep the highest component.
          for (auto it = local_trees.begin(); it != local_trees.end(); ++it) {
            Persistent_component_ptr neighbor = *it;

            neighbor->end_height() = height;
            if (neighbor->size() > size_max) {

              size_max = neighbor->size();
              chosen = neighbor;
            }
          }

          // Add current cell.
          chosen->add(sorted[i]);
          tree_map(x, y) = chosen;

          // Merge other components.
          for (auto it = local_trees.begin(); it != local_trees.end(); ++it) {
            Persistent_component_ptr neighbor = *it;
            if (neighbor == chosen)
              continue;

            // Components with the size above threshold are trees.
            if (neighbor->size() < min_height) {
              for (std::size_t n = 0; n < neighbor->inliers().size(); ++n) {
                const Coordinate& inlier = neighbor->inliers()[n];

                tree_map(inlier.first, inlier.second) = chosen;
                chosen->add(inlier);
              }
              neighbor->inliers().clear();
            }
          }
        }
      }

      for (std::size_t i = 0; i < components.size(); ++i) {
        if (components[i]->inliers().empty())
          continue;

        clusters.push_back(Iterators());
        for (std::size_t j = 0; j < components[i]->inliers().size(); ++j) {
          const Coordinate& coordinate = components[i]->inliers()[j];

          const Iterators& cell = image(coordinate.first, coordinate.second);
          std::copy(
            cell.begin(), cell.end(),
            std::back_inserter(clusters.back()));
        }
      }
    }

  }; // Vegetation_clustering

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_VEGETATION_CLUSTERING_H
