#ifndef CGAL_LEVELS_OF_DETAIL_TREE_ICON_CREATOR_H
#define CGAL_LEVELS_OF_DETAIL_TREE_ICON_CREATOR_H

// STL includes.
#include <map>
#include <vector>
#include <utility>
#include <algorithm>

// CGAL includes.
#include <CGAL/array.h>
#include <CGAL/number_utils.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>
#include <CGAL/Levels_of_detail/internal/structures.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputCluster,
  typename PointMap>
  class Tree_icon_creator {

  public:
    using Traits = GeomTraits;
    using Cluster = InputCluster;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Plane_3 = typename Traits::Plane_3;
    using Vector_2 = typename Traits::Vector_2;

    using Tree = Tree<Traits>;

    typename Traits::Compute_squared_distance_2 squared_distance_2;

  public:

    Tree_icon_creator(
      const std::vector<Cluster>& clusters,
      const Point_map point_map,
      const Plane_3& ground_plane,
      const FT precision) : 
    m_clusters(clusters),
    m_point_map(point_map),
    m_ground_plane(ground_plane),
    m_precision(precision)
    { }

    void create(std::vector<Tree>& trees) const {

      create_models(trees);
      create_triangles(trees);
    }

  private:
    const std::vector<Cluster>& m_clusters;
    const Point_map m_point_map;
    const Plane_3& m_ground_plane;
    const FT m_precision;

    void create_models(std::vector<Tree>& trees) const {
      
      for (std::size_t i = 0; i < trees.size(); ++i) {
        auto& tree = trees[i];

        const Point_2& center = tree.center;
        const FT radius = tree.radius;
        const FT height = tree.height;

        Cluster cluster = m_clusters[tree.cluster_index];
        using Value_type = typename Cluster::value_type;

        std::sort(cluster.begin(), cluster.end(), 
        [&](const Value_type& a, const Value_type& b) -> bool {
          
          return get(m_point_map, a).z() < get(m_point_map, b).z();
        });

        auto& model = tree.model;

        model.height[0] = 
        get(m_point_map, cluster.front()).z();
        model.height[1] = 
        get(m_point_map, cluster[cluster.size() / 10]).z();
        model.height[2] = 
        get(m_point_map, cluster[cluster.size() / 2]).z();
        model.height[3] = 
        get(m_point_map, cluster[static_cast<std::size_t>(9.0 * cluster.size() / 10.0)]).z();
        model.height[4] = 
        get(m_point_map, cluster.back()).z();

        cpp11::array<FT, 4> width = 
          make_array(FT(0), FT(0), FT(0), FT(0));

        cpp11::array<std::size_t, 4> nb = 
          make_array(std::size_t(0), std::size_t(0), std::size_t(0), std::size_t(0));

        for (std::size_t j = 0; j < cluster.size(); ++j) {
          
          const Point_3& p = get(m_point_map, cluster[j]);
          const Point_2 q = Point_2(p.x(), p.y());

          std::size_t idx = 0;
          if (p.z() < model.height[1])
            idx = 0;
          else if (p.z() < model.height[2])
            idx = 1;
          else if (p.z() < model.height[3])
            idx = 2;
          else
            idx = 3;

          width[idx] += squared_distance_2(q, center);
          nb[idx]++;
        }

        for (std::size_t j = 0; j < width.size(); ++j)
          if (nb[j] != 0)
            width[j] = static_cast<FT>(
              CGAL::sqrt(
                CGAL::to_double(width[j] / nb[j])));

        model.width[0] = (width[0] + width[1]) / FT(2);
        model.width[1] = (width[1] + width[2]) / FT(2);
        model.width[2] = (width[2] + width[3]) / FT(2);
      }
    }

    void create_triangles(std::vector<Tree>& trees) const {

      for (std::size_t i = 0; i < trees.size(); ++i) {
        auto& tree = trees[i];

        auto& vertices = tree.vertices;
        auto& faces = tree.faces;

        vertices.clear();
        faces.clear();

        const auto& footprint = tree.triangles;
    
        std::map<Point_2, FT> ground_points;
        for (std::size_t j = 0; j < footprint.size(); ++j) {
          for (std::size_t k = 0; k < 3; ++k) {
            const Point_2& p = footprint[j][k];
            
            const Point_3 q = internal::position_on_plane_3(p, m_ground_plane);
            ground_points.insert(std::make_pair(p, q.z()));
          }
        }

        std::size_t nb_pts = ground_points.size();

        const Point_2& center = tree.center;
        const FT radius = tree.radius / FT(6);
        const auto& model = tree.model;

        Point_2 p2_max = center;
        const FT h_max = model.height[4];
        Point_3 p_max(p2_max.x(), p2_max.y(), h_max);
        vertices.push_back(p_max);

        for (std::size_t j = 0; j < nb_pts; ++j) {
          
          const FT angle = FT(2) * static_cast<FT>(CGAL_PI) * 
          ( static_cast<FT>(j) / static_cast<FT>(nb_pts) );
          
          const Point_2 p2_trunk = center + Vector_2(
            radius * static_cast<FT>(std::cos(CGAL::to_double(angle))), 
            radius * static_cast<FT>(std::sin(CGAL::to_double(angle)))
            );

          const FT h_trunk = 
          (internal::position_on_plane_3(p2_trunk, m_ground_plane)).z();

          const Point_3 p_trunk(p2_trunk.x(), p2_trunk.y(), h_trunk);
          vertices.push_back(p_trunk);
          
          const Point_2 p2_min = center + Vector_2(
            radius * static_cast<FT>(std::cos(CGAL::to_double(angle))), 
            radius * static_cast<FT>(std::sin(CGAL::to_double(angle)))
            );

          const FT h_min = model.height[0] + FT(2);
          Point_3 p_min(p2_min.x(), p2_min.y(), h_min);
          vertices.push_back(p_min);

          for (std::size_t k = 1; k < 4; ++k) {
            
            const Point_2 p2 = center + Vector_2(
              model.width[k-1] * static_cast<FT>(std::cos(CGAL::to_double(angle))), 
              model.width[k-1] * static_cast<FT>(std::sin(CGAL::to_double(angle)))
              );

            const FT h = model.height[k] + FT(2);
            Point_3 p(p2.x(), p2.y(), h);
            vertices.push_back(p);
          }
        }

        for (std::size_t j = 0; j < nb_pts; ++j) {
          
          const std::size_t l = 1 + 5 * j;
          const std::size_t r = 1 + 5 * ((j + 1) % nb_pts);

          for (std::size_t k = 0; k < 4; ++k) {
            
            faces.push_back(make_array(l + k, r + k, l + 1 + k));
            faces.push_back(make_array(l + 1 + k, r + k, r + 1 + k));
          }
          faces.push_back(make_array(l + 4, r + 4, 0));
        }
      }
    }

  }; // Tree_icon_creator
    
} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_TREE_ICON_CREATOR_H
