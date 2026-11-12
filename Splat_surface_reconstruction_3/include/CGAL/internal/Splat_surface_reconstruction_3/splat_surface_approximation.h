#ifndef CGAL_INTERNAL_SPLAT_SURFACE_RECONSTRUCTION_3_SPLAT_SURFACE_APPROXIMATION_H
#define CGAL_INTERNAL_SPLAT_SURFACE_RECONSTRUCTION_3_SPLAT_SURFACE_APPROXIMATION_H

#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Kernel/global_functions_3.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

#include "box_grid.h"
#include "geometry_utils.h"
#include "surface_approximation.h"

namespace CGAL {
namespace internal {
namespace Splat_surface_reconstruction_3 {

template <typename Kernel_>
class Splat_surface_approximation : public Surface_approximation<Kernel_> {
  public:
    using Kernel   = Kernel_;
    using FT       = typename Kernel::FT;
    using Point_3  = typename Kernel::Point_3;
    using Point_2  = typename Kernel::Point_2;
    using Vector_3 = typename Kernel::Vector_3;
    using Index    = std::size_t;
    using Base = Surface_approximation<Kernel_>;
    using Initial_seed = typename Base::Initial_seed;
    using Feature_segment = typename Base::Feature_segment;

    Splat_surface_approximation(
      const std::vector<Point_3>& points,
      const std::vector<Vector_3>& normals,
      const Box_grid<Kernel_>& grid
    ) : 
      grid_(grid)
    { 
      points_ = points;
      normals_ = normals;
      cells_.resize(grid_.number_of_cells());

      CGAL_assertion(points_.size() == normals_.size());

      // Add each point to the grid boxes
      for (Index i = 0; i < points_.size(); ++i) {
        double len2 = CGAL::to_double(normals_[i].squared_length());
        insert(i, points_[i]);
      }

      filter_points_by_normals();

      compute_block_normals();

      estimate_individual_splat_sizes();
    }



    /**
     * @brief Selects an initial pair of mesh vertices from a large splat.
     *
     * The seed is placed symmetrically around the splat center along one tangent
     * direction. A central region of the input bounding box is preferred to avoid
     * starting the reconstruction near the boundary.
     *
     * @return Initial seed information.
    */
    Initial_seed get_initial_seed() const override {

      const FT box_size = grid_.get_box_size();

      std::cout<<"Searching for initial seed with splat size larger than box size..." << std::endl;
      for (Index i = 0; i < splat_sizes_.size(); ++i) {
        // CHANGED: from requiring splat_sizes_[i] > box_size_
        // to splat_sizes_[i] > box_size_ * 0.95
        if (splat_sizes_[i] > box_size * 0.95 && i < normals_.size() && normals_[i] != CGAL::NULL_VECTOR) {
          // Define a central region: keep points that are at least 10% away from each face.
          const FT margin_x = FT(0.1) * (grid_.get_max_x() - grid_.get_min_x());
          const FT margin_y = FT(0.1) * (grid_.get_max_y() - grid_.get_min_y());
          const FT margin_z = FT(0.1) * (grid_.get_max_z() - grid_.get_min_z());
          const FT mid_min_x = grid_.get_min_x() + margin_x;
          const FT mid_max_x = grid_.get_max_x() - margin_x;
          const FT mid_min_y = grid_.get_min_y() + margin_y;
          const FT mid_max_y = grid_.get_max_y() - margin_y;
          const FT mid_min_z = grid_.get_min_z() + margin_z;
          const FT mid_max_z = grid_.get_max_z() - margin_z;

          auto is_in_middle = [&](const Point_3& p) -> bool {
            return (p.x() >= mid_min_x && p.x() <= mid_max_x &&
                    p.y() >= mid_min_y && p.y() <= mid_max_y &&
                    p.z() >= mid_min_z && p.z() <= mid_max_z);
          };

          // Require the seed point to be somewhere in the middle of the cloud.
          if (!is_in_middle(points_[i])) {
            continue;
          }

          auto frame = compute_local_tangent_frame<Point_3, Vector_3, FT>(normals_[i]);
          const Vector_3& u = frame[0];

          const Point_3 seed1 = points_[i] + FT(0.5) * box_size * u;
          const Point_3 seed2 = points_[i] - FT(0.5) * box_size * u;

          // Make sure the two initial vertices are also not near the boundary.
          if (!is_in_middle(seed1) || !is_in_middle(seed2)) {
            continue;
          }

          return Initial_seed{seed1, seed2, normals_[i]};
        }
      }

      std::cerr << "Warning: no suitable initial seed found with splat size larger than box size.\n";

      // Fallback to a default initial seed.
      auto frame = compute_local_tangent_frame<Point_3, Vector_3, FT>(normals_[0]);
      const Vector_3& u = frame[0];
      const Point_3 seed1 = points_[1000] + FT(0.5) * box_size * u;
      const Point_3 seed2 = points_[1000] - FT(0.5) * box_size * u;
      return Initial_seed{seed1, seed2, normals_[0]};
    }



    /**
     * @brief Returns the averaged normal associated with the cell containing p.
     *
     * @param p Query point.
     *
     * @return Cell normal, or `CGAL::NULL_VECTOR` if the point is outside the grid.
    */
    Vector_3 box_normal(const Point_3& p) const override {
      int ix, iy, iz;

      CGAL_assertion(grid_.to_grid_coordinates(p, ix, iy, iz));

      if (!grid_.to_grid_coordinates(p, ix, iy, iz)) {
        return CGAL::NULL_VECTOR;
      }

      return block_normals_[grid_.flat_index(ix,iy,iz)];
    }

    std::vector<std::pair<Point_3, Vector_3>> intersect_circle(const Point_3& circle_center, const Vector_3& circle_normal, const FT circle_radius) const override{
      std::vector<std::pair<Point_3, Vector_3>> result;

      int cx, cy, cz;

      if (!grid_.to_grid_coordinates(circle_center, cx, cy, cz)) {
        return result;
      }

      // Iterate over all nearby splats
      int search_radius = static_cast<int>(std::ceil(circle_radius / grid_.get_box_size()));
      for (int dx = -search_radius; dx <= search_radius; ++dx) {
        for (int dy = -search_radius; dy <= search_radius; ++dy) {
          for (int dz = -search_radius; dz <= search_radius; ++dz) {
            int nx = cx + dx;
            int ny = cy + dy;
            int nz = cz + dz;

            if (!grid_.valid_coordinates(nx, ny, nz)) {
              continue;
            }

            const std::vector<Index>& box_point_ids = point_ids_in_box(
              static_cast<std::size_t>(nx),
              static_cast<std::size_t>(ny),
              static_cast<std::size_t>(nz)
            );

            for (const Index& splat_idx : box_point_ids) {
              std::vector<Point_3> hits = circle_intersects_splat<Point_3, Vector_3, FT>(
                circle_center, circle_normal, circle_radius, 
                splat_position(splat_idx),
                splat_normal(splat_idx),
                splat_size(splat_idx)
              );

              for (const Point_3& hit : hits) {
                result.emplace_back(hit, splat_normal(splat_idx));
              }
            }
          }
        }
      }

      return result;
    }


    std::vector<Feature_segment> find_feature_segments(FT minimum_angle) override {
      std::vector<Feature_segment> feature_segments;

      // Iterate over all splats and identify feature segments
      for (std::size_t i = 0; i < points_.size(); ++i) {
        const Point_3& center = splat_position(i);
        const Vector_3 normal = splat_normal(i);
        const FT radius = grid_.get_box_size();

        // Attempt intersection with other splats
        std::vector<Index> splats_in_block = nearby_point_ids(center, radius+radius);
        for (Index j : splats_in_block) {
          if (j < i) continue;
          Point_3 intersection_point_1, intersection_point_2;

          // If the splats intersect, add the feature segment
          if (splat_intersects_splat(i, j, intersection_point_1, intersection_point_2, minimum_angle)) {
            // Find the normal of the intersecting splat
            Vector_3 normal = splat_normal(i);
            normal += splat_normal(j);
            normal = normal / CGAL::sqrt(normal.squared_length());

            std::size_t feature_segment_index = feature_segments.size();

            feature_segments.emplace_back(Feature_segment{intersection_point_1, intersection_point_2, normal});
          }
        }
      }

      return feature_segments;
    }


    std::vector<std::pair<Point_3, Vector_3>> find_initial_feature_vertex_candidates(FT minimum_angle) override {
      std::vector<std::pair<Point_3, Vector_3>> initial_candidates;

      // Iterate over all splats to find triple intersections
      for (std::size_t i = 0; i < points_.size(); ++i) {
        const Point_3& point_i = splat_position(i);
        const Vector_3 normal_i = splat_normal(i);
        const FT radius = grid_.get_box_size();

        std::vector<Index> splats_in_block = nearby_point_ids(point_i, radius+radius);
        for (Index j : splats_in_block) {
          if (i == j) continue;

          for (Index k : splats_in_block) {
            if (i == k || j == k) continue;

            const Vector_3& normal_j = splat_normal(j);
            const Point_3& point_j = splat_position(j);

            const Vector_3& normal_k = splat_normal(k);
            const Point_3& point_k = splat_position(k);

            Vector_3 normal = normal_i + normal_j + normal_k;
            normal = normal / CGAL::sqrt(normal.squared_length());


            Point_3 intersection_point;
            if (triple_plane_intersection<Point_3, Vector_3, FT, Kernel>(point_i, normal_i, radius,
              point_j, normal_j, radius,
              point_k, normal_k, radius,
              minimum_angle,
              intersection_point))
            {
              initial_candidates.push_back({intersection_point, normal});
            }
          }
        }
      }

      return initial_candidates;
    }

    void write_to_obj(const char * filename) const {

      typename CGAL::Surface_mesh<Point_3> mesh;

      for (std::size_t i = 0; i < points_.size(); ++i) {
        FT global_spacing = 1.05*grid_.get_box_size();

        FT radius = splat_sizes_.size() == 0 ? global_spacing : splat_sizes_[i];

        add_splat(mesh, points_[i], normals_[i], radius);
      }

      CGAL::IO::write_polygon_mesh(filename, mesh);

      std::cout << "Splat mesh written to " << filename << std::endl;
    }

  private:

    /**
     * @brief Filters points within each grid cell based on their normals.
     */
    void filter_points_by_normals() {
      constexpr double min_normal_sum_length = 0.1;
      constexpr double min_normal_sum_length_squared =
        min_normal_sum_length * min_normal_sum_length;

      for (Cell& cell : cells_) {
        auto& ids = cell.point_ids_within_distance_d;

        if (ids.empty()) {
          continue;
        }

        // Step 1: Sum the associated point normals without normalization.
        Vector_3 normal_sum = CGAL::NULL_VECTOR;

        for (Index pid : ids) {
          normal_sum = normal_sum + splat_normal(pid);
        }

        // Step 2: If the sum is too small, keep all associated points.
        const double length_squared =
          CGAL::to_double(normal_sum.squared_length());

        if (length_squared < min_normal_sum_length_squared) {
          continue;
        }

        // Step 3: Remove points whose normals oppose the summed normal.
        ids.erase(
          std::remove_if(
            ids.begin(),
            ids.end(),
            [&](Index pid) {
              return CGAL::to_double(splat_normal(pid) * normal_sum) < 0.0;
            }),
          ids.end());
      }
    }

    /**
     * @brief Computes box normals using the finite-sampling optimization
     *        described in the Feature-aware Manifold Meshing paper.
     *
     * For each box, selects the sampled unit direction N that maximizes
     * the minimum dot product between N and all associated point normals.
     *
     * Assumes point_ids_within_distance_d contains the distance-associated
     * points after the preliminary normal-filtering step.
     */
    std::vector<Vector_3> compute_block_normals()
    {
      // Number of directions in the finite sampling of S^2.
      constexpr std::size_t sample_count = 1000;

      // Construct a deterministic, approximately uniform sampling of S^2.
      std::vector<Vector_3> sphere_samples;
      sphere_samples.reserve(sample_count);

      const double golden_angle =
        CGAL_PI * (3.0 - std::sqrt(5.0));

      for (std::size_t i = 0; i < sample_count; ++i) {
        const double z =
          1.0 - 2.0 * (static_cast<double>(i) + 0.5)
                      / static_cast<double>(sample_count);

        const double radius = std::sqrt(1.0 - z * z);
        const double theta = golden_angle * static_cast<double>(i);

        sphere_samples.emplace_back(
          radius * std::cos(theta),
          radius * std::sin(theta),
          z);
      }

      block_normals_.assign(cells_.size(), CGAL::NULL_VECTOR);

      for (std::size_t idx = 0; idx < cells_.size(); ++idx) {
        const Cell& cell = cells_[idx];

        // Use the distance-associated point list, not cell membership.
        const auto& point_ids = cell.point_ids_within_distance_d;

        if (point_ids.empty()) {
          continue;
        }

        // Find the sampled direction with the largest minimum dot product.
        double best_score = -std::numeric_limits<double>::infinity();
        Vector_3 best_normal = CGAL::NULL_VECTOR;

        for (const Vector_3& candidate : sphere_samples) {
          double min_dot = std::numeric_limits<double>::infinity();

          for (Index pid : point_ids) {
            const double dot = CGAL::to_double(
              splat_normal(pid) * candidate);

            min_dot = (std::min)(min_dot, dot);
          }

          if (min_dot > best_score) {
            best_score = min_dot;
            best_normal = candidate;
          }
        }

        block_normals_[idx] = best_normal;
      }

      return block_normals_;
    }

    /**
     * @brief Estimates one splat radius for every input point.
     *
     * A local 2D Delaunay triangulation is constructed in the tangent plane of
     * each input point. The largest incident circumcenter distance is used as
     * the local splat radius and is clamped to the global spacing.
     *
     * @return Per-point splat radii.
    */
    std::vector<FT> estimate_individual_splat_sizes() {
      FT global_spacing = grid_.get_box_size();
      splat_sizes_.resize(points_.size(), global_spacing);

      std::cout << "Estimating individual splat sizes with global spacing: " << global_spacing << std::endl;

      std::size_t computation_count = 0;

      for (Index i = 0; i < points_.size(); ++i) {
        int cx, cy, cz;
        if (!grid_.to_grid_coordinates(splat_position(i), cx, cy, cz)) {
          std::cerr << "Warning: point " << i << " is out of grid bounds, skipping splat size estimation.\n";
          std::exit(EXIT_FAILURE);
        }

        // Changed: from
        // Collect all points in the same box as p_i.
        // to:
        // Collect all points within distance d of p_i.
        std::vector<Index> local_ids;
        const Cell* c = cell(cx, cy, cz);
        if (c == nullptr || c->point_ids_within_distance_d.empty()) {
          continue;
        }

        for (Index pid : c->point_ids_within_distance_d) {
          // Distance check: only include points within distance d of p_i.
          Vector_3 diff = splat_position(pid) - splat_position(i);
          double dist2 = CGAL::to_double(diff * diff);
          double d_sq = 1.0 * grid_.get_box_size() * grid_.get_box_size();
          if (dist2 <= d_sq + 1e-6 || true) {
            local_ids.push_back(pid);
          }
        }

        // Need at least a few neighbors to define a local triangulation.
        if (local_ids.size() < 3) {
          continue;
        }

        ++computation_count;

        // Build a local tangent frame at p_i from its normal.
        std::vector<Vector_3> local_frame = compute_local_tangent_frame<Point_3, Vector_3, FT>(normals_[i]);
        Vector_3 u = local_frame[0];
        Vector_3 v = local_frame[1];

        // Project local neighbors onto the tangent plane and compute their 2D coordinates
        // Also input in the Delaunay Triangulation Kernel
        CGAL::Delaunay_triangulation_2<Kernel> DT;
        auto center_vh = DT.insert(Point_2(0, 0)); // keep the handle

        for (Index neighbor_id : local_ids) {
          if (neighbor_id == i) continue;

          Vector_3 vec = points_[neighbor_id] - points_[i];
          double x = CGAL::to_double(vec * u);
          double y = CGAL::to_double(vec * v);
          DT.insert(Point_2(x, y));
        }

        // Find the circumcenter of faces incident to the center vertex only.
        std::vector<FT> circumcenter_distances;
        FT max_r2 = FT(0);

        auto fc = DT.incident_faces(center_vh);
        if (fc != 0) {
          auto done = fc;

          do {
            if (!DT.is_infinite(fc)) {
              Point_2 a = fc->vertex(0)->point();
              Point_2 b = fc->vertex(1)->point();
              Point_2 cpt = fc->vertex(2)->point();

              Point_2 cc = CGAL::circumcenter(a, b, cpt);
              const FT r2 = cc.x() * cc.x() + cc.y() * cc.y();

              if (r2 > max_r2) {
                max_r2 = r2;
              }
            }
            ++fc;
          } while (fc != done);
        }

        if (max_r2 > FT(0)) {
          splat_sizes_[i] = CGAL::sqrt(max_r2);
          splat_sizes_[i] = (std::min)(splat_sizes_[i], global_spacing); // Clamp to global spacing to avoid outliers.
        }
      }

      std::cout << "Total local splat size computations: " << computation_count << " / " << splat_sizes_.size() << std::endl;

      max_splat_radius_ = *std::max_element(splat_sizes_.begin(), splat_sizes_.end());
      std::cout<<"Box size: " << grid_.get_box_size() << " Max splat radius: " << max_splat_radius_ << std::endl;
      return splat_sizes_;
    }


    /**
     * @brief Clears all per-cell accumulators.
     */
    void clear() {
      for (Cell& c : cells_) {
        c.point_ids_contained.clear();
        c.point_ids_within_distance_d.clear();
      }

      // Make sure all cells are cleared.
      for (const Cell& c : cells_) {
        CGAL_assertion(c.point_ids_contained.empty());
        CGAL_assertion(c.point_ids_within_distance_d.empty());
      }
    }

    /**
     * @brief Inserts a single point into the grid.
     *
     * @param point_id Index of the point in the input cloud.
     * @param p        Point position.
     */
    void insert(Index point_id, const Point_3& p) {
      int ix, iy, iz;
      if (!grid_.to_grid_coordinates(p, ix, iy, iz)) {
        return;
      }

      // Always add the point to the cell it falls into.
      const std::size_t idx = grid_.flat_index(ix, iy, iz);
      cells_[idx].add_point_contained(point_id);

      // And add it to all neighboring cells within distance d.
      for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
          for (int dz = -1; dz <= 1; ++dz) {
            int nx = ix + dx;
            int ny = iy + dy;
            int nz = iz + dz;
            if (grid_.valid_coordinates(nx, ny, nz) && grid_.distance_to_box(p, nx, ny, nz) <= grid_.get_box_size()) {
              const std::size_t nidx = grid_.flat_index(nx, ny, nz);
              cells_[nidx].add_point_within_distance_d(point_id);
            }
          }
        }
      }
    }



    /**
     * @brief Finds input points in grid cells intersecting a query neighborhood.
     *
     * @param p Query position.
     * @param radius Search radius.
     *
     * @return Unique input point identifiers in the neighborhood.
    */
    std::vector<Index> nearby_point_ids(const Point_3& p, FT radius) const {
      std::vector<Index> ids;

      int cx, cy, cz;
      if (!grid_.to_grid_coordinates(p, cx, cy, cz)) {
        return ids;
      }

      const int r = (std::max)(
        1,
        static_cast<int>(
          std::ceil(CGAL::to_double(radius) / CGAL::to_double(grid_.get_box_size()))
        )
      );

      for (int dx = -r; dx <= r; ++dx) {
        for (int dy = -r; dy <= r; ++dy) {
          for (int dz = -r; dz <= r; ++dz) {
            const int ix = cx + dx;
            const int iy = cy + dy;
            const int iz = cz + dz;

            const Cell* c = cell(ix, iy, iz);
            if (c == nullptr) {
              continue;
            }

            ids.insert(ids.end(), c->point_ids_contained.begin(), c->point_ids_contained.end());
          }
        }
      }

      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
      return ids;
    }

    /**
     * @brief Returns the point IDs of points within distance d of the specified cell's bounding box.
     * @param ix The x-coordinate of the cell in the grid.
     * @param iy The y-coordinate of the cell in the grid.
     * @param iz The z-coordinate of the cell in the grid.
     * @return A reference to the vector of point IDs within distance d of the cell's bounding box.
     */
    const std::vector<Index>& point_ids_within_distance_d_of_box(std::size_t ix, std::size_t iy, std::size_t iz) const {
      const Cell* c = cell(ix, iy, iz);
      CGAL_assertion(c != nullptr);

      return c->point_ids_within_distance_d;
    }

    /**
     * @brief Returns the point IDs of points contained within the specified cell's bounding box.
     * @param ix The x-coordinate of the cell in the grid.
     * @param iy The y-coordinate of the cell in the grid.
     * @param iz The z-coordinate of the cell in the grid.
     * @return A reference to the vector of point IDs contained within the cell's bounding box.
     */
    const std::vector<Index>& point_ids_in_box(std::size_t ix, std::size_t iy, std::size_t iz) const {
      const Cell* c = cell(ix, iy, iz);
      CGAL_assertion(c != nullptr);

      return c->point_ids_contained;
    }


    /**
     * @brief Determines if two splats intersect and returns the intersection points if they do.
     * @param splat_1 The index of the first splat.
     * @param splat_2 The index of the second splat.
     * @param intersection_point_1 The first intersection point (output).
     * @param intersection_point_2 The second intersection point (output).
     * @param minimum_angle The minimum allowed angle between splat normals.
     * @return True if the splats intersect, false otherwise.
     */
    bool splat_intersects_splat(Index splat_1, Index splat_2, Point_3& intersection_point_1, Point_3& intersection_point_2, FT minimum_angle) const {
      const Point_3& splat_1_center = splat_position(splat_1);
      const Vector_3 splat_1_normal = splat_normal(splat_1);
      const FT splat_1_radius = grid_.get_box_size();

      const Point_3& splat_2_center = splat_position(splat_2);
      const Vector_3 splat_2_normal = splat_normal(splat_2);
      const FT splat_2_radius = grid_.get_box_size();

      FT distance_between_centers = CGAL::sqrt(CGAL::squared_distance(splat_1_center, splat_2_center));
      if (distance_between_centers > splat_1_radius + splat_2_radius) {
        return false;
      }

      FT angle_between_normals = CGAL::approximate_angle(splat_1_normal, splat_2_normal);
      if (angle_between_normals < minimum_angle) {
        return false;
      }


      std::vector<Point_3> result_1 = circle_intersects_splat<Point_3, Vector_3, FT>(splat_1_center, splat_1_normal, splat_1_radius, splat_2_center, splat_2_normal, splat_2_radius);
      std::vector<Point_3> result_2 = circle_intersects_splat<Point_3, Vector_3, FT>(splat_2_center, splat_2_normal, splat_2_radius, splat_1_center, splat_1_normal, splat_1_radius);

      std::vector<Point_3> intersection_points = result_1;
      intersection_points.insert(intersection_points.end(), result_2.begin(), result_2.end());

      std::vector<Point_3> unique_intersection_points;

      for (const Point_3& p : intersection_points) {
        bool is_duplicate = false;
        for (const Point_3& q : unique_intersection_points) {
          if (CGAL::squared_distance(p, q) <= FT(1e-12)) {
            is_duplicate = true;
            break;
          }
        }
        if (!is_duplicate) {
          unique_intersection_points.push_back(p);
        }
      }

      // If there are exactly two unique intersection points, we have a valid intersection
      if (unique_intersection_points.size() == 2) {
        intersection_point_1 = unique_intersection_points[0];
        intersection_point_2 = unique_intersection_points[1];
        return true;
      }
      return false;
    }



    void set_all_splat_sizes(FT size) {
      splat_sizes_.assign(points_.size(), size);
      max_splat_radius_ = size;
    }

    void add_splat(CGAL::Surface_mesh<Point_3>& m,
                const Point_3& center,
                const Vector_3& normal,
                FT radius) const
    {
        const int num_vertices = 10;

        std::vector<Vector_3> tangent_frame_vectors = compute_local_tangent_frame<Point_3, Vector_3, FT>(normal);

        using Vertex_index = typename CGAL::Surface_mesh<Point_3>::Vertex_index;

        // Create the center vertex once.
        Vertex_index center_vertex = m.add_vertex(center);

        // Create the ring vertices once.
        std::vector<Vertex_index> ring;
        ring.reserve(num_vertices);

        for (int i = 0; i < num_vertices; ++i) {
            const double angle = 2.0 * M_PI * double(i) / double(num_vertices);

            Point_3 p = center + radius *
                (std::cos(angle) * tangent_frame_vectors[0] + std::sin(angle) * tangent_frame_vectors[1]);

            ring.push_back(m.add_vertex(p));
        }

        // Connect adjacent ring vertices to the shared center.
        for (int i = 0; i < num_vertices; ++i) {
            const int next = (i + 1) % num_vertices;
            m.add_face(center_vertex, ring[i], ring[next]);
        }
    }

    void add_triangle(CGAL::Surface_mesh<Point_3>& m, const Point_3& p1, const Point_3& p2, const Point_3& p3) const {
      typename CGAL::Surface_mesh<Point_3>::Vertex_index u = m.add_vertex(p1);
      typename CGAL::Surface_mesh<Point_3>::Vertex_index v = m.add_vertex(p2);
      typename CGAL::Surface_mesh<Point_3>::Vertex_index w = m.add_vertex(p3);
      m.add_face(u, v, w);
    }
    bool write_splat_sizes_to_file(const std::string& filename) const
    {
      std::ofstream out(filename);
      if (!out) {
        std::cerr << "Error: cannot open " << filename << " for writing.\n";
        return false;
      }

      for (const FT& size : splat_sizes_) {
        out << std::setprecision(17) << CGAL::to_double(size) << '\n';
      }

      return true;
    }

    /**
     * @brief A single grid cell.
     *
     * Each cell stores the indices of the points that fall inside it together
     * with the sum of their normals and the number of contributing normals.
     */
    struct Cell
    {
      // The indices of all input points within distance d of the cell center.
      std::vector<Index> point_ids_within_distance_d;
      std::vector<Index> point_ids_contained;

      void add_point_contained(Index pid) {
        point_ids_contained.push_back(pid);
      }
      void add_point_within_distance_d(Index pid) {
        point_ids_within_distance_d.push_back(pid);
      }
    };

    /**
     * @brief Returns the cell at integer grid coordinates.
     *
     * @param ix Cell index along x.
     * @param iy Cell index along y.
     * @param iz Cell index along z.
     *
     * @return Pointer to the requested cell, or `nullptr` if the coordinates
     *         are outside the grid.
     */
    const Cell* cell(int ix, int iy, int iz) const {
      if (!grid_.valid_coordinates(ix, iy, iz)) {
        return nullptr;
      }
      return &cells_[grid_.flat_index(ix, iy, iz)];
    }

    const Point_3& splat_position(Index idx) const {
      CGAL_assertion(idx < points_.size());
      return points_[idx];
    }

    const Vector_3& splat_normal(Index idx) const {
      CGAL_assertion(idx < normals_.size());
      return normals_[idx];
    }

    const FT& splat_size(Index idx) const {
      CGAL_assertion(idx < splat_sizes_.size());
      return splat_sizes_[idx];
    }

    std::vector<Point_3> points_;
    std::vector<Vector_3> normals_;
    std::vector<FT> splat_sizes_;
    const Box_grid<Kernel_>& grid_;
    std::vector<Cell> cells_;
    FT max_splat_radius_;
    std::vector<Vector_3> block_normals_;
};

} // namespace Splat_surface_reconstruction_3
} // namespace internal
} // namespace CGAL

#endif