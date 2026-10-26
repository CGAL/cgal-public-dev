// Copyright (c) 2026  GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pranav Jain

/**
 * @file splat_surface_reconstruction.h
 * @brief Splat surface reconstruction for point clouds with normals.
 *
 * This header provides the public reconstruction entry point together with
 * an internal dense box-grid acceleration structure used for preprocessing,
 * local normal estimation, and splat-size estimation.
 */
#ifndef CGAL_SPLAT_SURFACE_RECONSTRUCTION_H
#define CGAL_SPLAT_SURFACE_RECONSTRUCTION_H

// #include <CGAL/license/Splat_surface_reconstruction_3.h>
#include <CGAL/property_map.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/boost/graph/dijkstra_shortest_path.h>
#include <CGAL/number_utils.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Polygon_2.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <optional>
#include <queue>

namespace CGAL {

  template <typename PointRange, typename NormalRange, typename TriangleMesh>
  struct Splat_surface_reconstruction_3_test_access;

  template <typename Kernel_>
  struct Grid_test_access;

  /**
   * @brief Dense axis-aligned box grid over a fixed 3D domain.
   *
   * The grid stores points and accumulated normals per cell. It is used as a
   * local acceleration structure for neighborhood queries, cell-normal
   * estimation, and adaptive splat-size estimation.
   *
   * @tparam Kernel_ CGAL kernel type used for points, vectors, and scalars.
   */
  template <typename Kernel_>
  class Box_grid {
    private:
      friend struct Grid_test_access<Kernel_>;

    public:
    using Kernel   = Kernel_;
    using FT       = typename Kernel::FT;
    using Point_3  = typename Kernel::Point_3;
    using Point_2  = typename Kernel::Point_2;
    using Vector_3 = typename Kernel::Vector_3;
    using Index    = std::size_t;

    struct Initial_seed {
      Point_3 first;
      Point_3 second;
      Index splat_id;
    };

    public:
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
      if (!valid_coordinates(ix, iy, iz)) {
        return nullptr;
      }
      return &cells_[flat_index(ix, iy, iz)];
    }

    /**
     * @brief Constructs a dense 3D grid.
     *
     * The grid uses the same cell size along all three axes.
     *
     * @param box_size Cell side length.
     * @param min_x Minimum x-coordinate of the grid.
     * @param min_y Minimum y-coordinate of the grid.
     * @param min_z Minimum z-coordinate of the grid.
     * @param max_x Maximum x-coordinate of the grid.
     * @param max_y Maximum y-coordinate of the grid.
     * @param max_z Maximum z-coordinate of the grid.
    */
    Box_grid(FT box_size,
             FT min_x, FT min_y, FT min_z,
             FT max_x, FT max_y, FT max_z)
     : box_size_(box_size),
       min_x_(min_x), max_x_(max_x), min_y_(min_y),
       max_y_(max_y), min_z_(min_z), max_z_(max_z) {
      initialize_grid();
    }

    /// @brief Returns the grid cell side length.
    FT get_box_size() const {
      return box_size_;
    }

    /// @brief Returns the largest estimated splat radius.
    FT get_max_splat_radius() const {
      return max_splat_radius_;
    }

    /// @brief Returns the input point positions stored by the grid.
    const std::vector<Point_3>& points() const {
      return points_;
    }

    /// @brief Returns the input point normals stored by the grid.
    const std::vector<Vector_3>& normals() const {
      return normals_;
    }

    /// @brief Returns the estimated splat radius of every input point.
    const std::vector<FT>& splat_sizes() const {
      return splat_sizes_;
    }

    /// @brief Return the computed block normals.
    const std::vector<Vector_3>& block_normals() const {
      return block_normals_;
    }

    /**
     * @brief Builds the grid from input points and normals.
     *
     * Each point is inserted into its corresponding grid cell and its normal is
     * accumulated into that cell.
     *
     * @param points Input point positions.
     * @param normals Input point normals.
    */
    void build(const std::vector<Point_3>& points,
              const std::vector<Vector_3>& normals) {

      clear();

      points_ = points;
      normals_ = normals;

      CGAL_assertion(points_.size() == normals_.size());

      // Add each point to the grid boxes
      for (Index i = 0; i < points_.size(); ++i) {
        double len2 = CGAL::to_double(normals_[i].squared_length());
        insert(i, points_[i]);
      }

      filter_points_by_normals();

      compute_block_normals();
    }

    // TODO write description
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
          CGAL_assertion(pid < normals_.size());
          normal_sum = normal_sum + normals_[pid];
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
              return CGAL::to_double(normals_[pid] * normal_sum) < 0.0;
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
            CGAL_assertion(pid < normals_.size());

            const double dot = CGAL::to_double(
              normals_[pid] * candidate);

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
      FT global_spacing = 1.05*box_size_;
      splat_sizes_.resize(points_.size(), global_spacing);

      for (Index i = 0; i < points_.size(); ++i) {
        int cx, cy, cz;
        if (!to_grid_coordinates(points_[i], cx, cy, cz)) {
          std::cerr << "Warning: point " << i << " is out of grid bounds, skipping splat size estimation.\n";
          std::exit(EXIT_FAILURE);
        }

        // Collect all points in the same box as p_i.
        std::vector<Index> local_ids;
        const Cell* c = cell(cx, cy, cz);
        if (c == nullptr || c->point_ids_contained.empty()) {
          continue;
        }
        local_ids = c->point_ids_contained;

        // Need at least a few neighbors to define a local triangulation.
        if (local_ids.size() < 3) {
          continue;
        }

        // Build a local tangent frame at p_i from its normal.
        std::vector<Vector_3> local_frame = compute_local_tangent_frame(normals_[i]);
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

      max_splat_radius_ = *std::max_element(splat_sizes_.begin(), splat_sizes_.end());
      std::cout<<"Box size: " << box_size_ << " Max splat radius: " << max_splat_radius_ << std::endl;
      return splat_sizes_;
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
    Initial_seed get_initial_seed() const {
      std::cout<<"Searching for initial seed with splat size larger than box size..." << std::endl;
      for (Index i = 0; i < splat_sizes_.size(); ++i) {
        if (splat_sizes_[i] > box_size_ && i < normals_.size() && normals_[i] != CGAL::NULL_VECTOR) {
          // Define a central region: keep points that are at least 10% away from each face.
          const FT margin_x = FT(0.1) * (max_x_ - min_x_);
          const FT margin_y = FT(0.1) * (max_y_ - min_y_);
          const FT margin_z = FT(0.1) * (max_z_ - min_z_);
          const FT mid_min_x = min_x_ + margin_x;
          const FT mid_max_x = max_x_ - margin_x;
          const FT mid_min_y = min_y_ + margin_y;
          const FT mid_max_y = max_y_ - margin_y;
          const FT mid_min_z = min_z_ + margin_z;
          const FT mid_max_z = max_z_ - margin_z;

          auto is_in_middle = [&](const Point_3& p) -> bool {
            return (p.x() >= mid_min_x && p.x() <= mid_max_x &&
                    p.y() >= mid_min_y && p.y() <= mid_max_y &&
                    p.z() >= mid_min_z && p.z() <= mid_max_z);
          };

          // Require the seed point to be somewhere in the middle of the cloud.
          if (!is_in_middle(points_[i])) {
            continue;
          }

          auto frame = compute_local_tangent_frame(normals_[i]);
          const Vector_3& u = frame[0];

          const Point_3 seed1 = points_[i] + FT(0.5) * box_size_ * u;
          const Point_3 seed2 = points_[i] - FT(0.5) * box_size_ * u;

          // Make sure the two initial vertices are also not near the boundary.
          if (!is_in_middle(seed1) || !is_in_middle(seed2)) {
            continue;
          }

          std::cout<<"Found initial seed at point index " << i << " with splat size " << splat_sizes_[i] << "." << std::endl;
          return Initial_seed{seed1, seed2, i};
        }
      }

      std::cerr << "Warning: no suitable initial seed found with splat size larger than box size.\n";

      // Fallback to a default initial seed.
      auto frame = compute_local_tangent_frame(normals_[0]);
      const Vector_3& u = frame[0];
      const Point_3 seed1 = points_[1000] + FT(0.5) * box_size_ * u;
      const Point_3 seed2 = points_[1000] - FT(0.5) * box_size_ * u;
      return Initial_seed{seed1, seed2, 0};
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
      if (!to_grid_coordinates(p, cx, cy, cz)) {
        return ids;
      }

      const int r = (std::max)(
        1,
        static_cast<int>(
          std::ceil(CGAL::to_double(radius) / CGAL::to_double(box_size_))
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
     * @brief Intersects the construction circle of two parent vertices with a splat.
     *
     * The circle is defined by the parent pair and the global box size. Its
     * intersection with the splat plane is computed analytically and only
     * intersection points inside the splat disk are returned.
     *
     * @param parent_a First parent point.
     * @param parent_b Second parent point.
     * @param splat_center Center of the splat.
     * @param splat_normal Splat normal.
     * @param splat_radius Splat radius.
     *
     * @return Zero, one, or two valid intersection points.
    */
    std::vector<Point_3> intersect_circle_with_splat(const Point_3& parent_a,
                                                    const Point_3& parent_b,
                                                    const Point_3& splat_center,
                                                    const Vector_3& splat_normal,
                                                    FT splat_radius) const {
      std::vector<Point_3> result;

      const Vector_3 ab = parent_b - parent_a;
      const FT d2 = ab.squared_length();
      if (d2 <= FT(0)) {
        return result;
      }

      // Circle radius from the two-parent construction.
      const FT circle_r2 = box_size_ * box_size_ - d2 / FT(4);
      if (circle_r2 <= FT(0)) {
        return result;
      }

      const FT circle_r = CGAL::sqrt(circle_r2);

      // Circle center = midpoint of the two parents.
      const Point_3 mid(
        (parent_a.x() + parent_b.x()) / FT(2),
        (parent_a.y() + parent_b.y()) / FT(2),
        (parent_a.z() + parent_b.z()) / FT(2)
      );

      // Normal of the circle plane: along the parent-parent direction.
      const double ab_len2 = CGAL::to_double(d2);
      if (ab_len2 <= 0.0) {
        return result;
      }
      Vector_3 circle_normal = ab * (1.0 / std::sqrt(ab_len2));

      // Build an orthonormal basis (u, v) for the circle plane.
      std::vector<Vector_3> frame = compute_local_tangent_frame(circle_normal);
      const Vector_3& u = frame[0];
      const Vector_3& v = frame[1];

      // Normalize the splat normal.
      Vector_3 ns = splat_normal;
      const double ns_len2 = CGAL::to_double(ns.squared_length());
      if (ns_len2 <= 0.0) {
        return result;
      }
      ns = ns * (1.0 / std::sqrt(ns_len2));

      // Circle points are:
      //   x(t) = mid + circle_r * (cos(t) * u + sin(t) * v)
      //
      // Splat plane equation:
      //   dot(ns, x - splat_center) = 0
      //
      // This becomes:
      //   A cos(t) + B sin(t) = C
      const FT A = circle_r * (u * ns);
      const FT B = circle_r * (v * ns);
      const FT C = (splat_center - mid) * ns;

      const FT R2 = A * A + B * B;
      if (R2 <= FT(0)) {
        // Circle plane is parallel to splat plane.
        return result;
      }

      const FT R = CGAL::sqrt(R2);

      if (std::abs(CGAL::to_double(C)) > CGAL::to_double(R)) {
        // No intersection between circle and splat plane.
        return result;
      }

      const FT phi = std::atan2(CGAL::to_double(B), CGAL::to_double(A));
      const FT delta = std::acos(CGAL::to_double(C / R));

      const FT t1 = phi + delta;
      const FT t2 = phi - delta;

      const auto make_point = [&](FT t) -> Point_3 {
        const FT ct = std::cos(CGAL::to_double(t));
        const FT st = std::sin(CGAL::to_double(t));
        return Point_3(
          mid.x() + circle_r * (ct * u.x() + st * v.x()),
          mid.y() + circle_r * (ct * u.y() + st * v.y()),
          mid.z() + circle_r * (ct * u.z() + st * v.z())
        );
      };

      const FT splat_r2 = splat_radius * splat_radius;

      Point_3 p1 = make_point(t1);
      Point_3 p2 = make_point(t2);

      // Keep only points that lie inside the splat disk.
      if (CGAL::squared_distance(p1, splat_center) <= splat_r2) {
        result.push_back(p1);
      }

      // Avoid duplicate if the two solutions coincide numerically.
      if (CGAL::squared_distance(p2, splat_center) <= splat_r2 &&
          CGAL::squared_distance(p2, p1) > FT(1e-12))
      {
        result.push_back(p2);
      }

      return result;
    }

    /**
     * @brief Computes an orthonormal tangent frame for a normal.
     *
     * The returned vectors span the plane orthogonal to the input normal.
     *
     * @param n Surface normal.
     *
     * @return Two approximately unit-length tangent vectors.
    */
    std::vector<Vector_3> compute_local_tangent_frame(const Vector_3& n) const {
      // Find a vector that is not parallel to n to construct the tangent frame.
      Vector_3 temp;
      if (std::abs(CGAL::to_double(n.x())) < 0.9)
        temp = Vector_3(1,0,0);
      else
        temp = Vector_3(0,1,0);

      Vector_3 u = CGAL::cross_product(n, temp);
      Vector_3 v = CGAL::cross_product(n, u);

      // Each vector should be orthogonal to the others.
      CGAL_assertion(u*v < 1e-8);
      CGAL_assertion(n*u < 1e-8);
      CGAL_assertion(n*v < 1e-8);

      // None should be the zero vector.
      CGAL_assertion(u.squared_length() > 0);
      CGAL_assertion(v.squared_length() > 0);
      CGAL_assertion(n.squared_length() > 0);

      // Normalize u and v to have unit length.
      const double len_u = std::sqrt(CGAL::to_double(u.squared_length()));
      const double len_v = std::sqrt(CGAL::to_double(v.squared_length()));
      if (len_u > 0) u = u / len_u;
      if (len_v > 0) v = v / len_v;

      return {u, v};
    }

    /**
     * @brief Returns the averaged normal associated with the cell containing p.
     *
     * @param p Query point.
     *
     * @return Cell normal, or `CGAL::NULL_VECTOR` if the point is outside the grid.
    */
    Vector_3 box_normal(const Point_3& p) const {
      int ix, iy, iz;

      CGAL_assertion(to_grid_coordinates(p, ix, iy, iz));

      if (!to_grid_coordinates(p, ix, iy, iz)) {
        return CGAL::NULL_VECTOR;
      }

      return block_normals_[flat_index(ix,iy,iz)];
    }

    /**
     * @brief Converts a 3D point into integer grid coordinates.
     *
     * @param p Input point.
     * @param ix Output x index.
     * @param iy Output y index.
     * @param iz Output z index.
     *
     * @return `true` if the point is inside the grid.
    */
    bool to_grid_coordinates(const Point_3& p, int& ix, int& iy, int& iz) const {

      const double h = CGAL::to_double(box_size_);

      const double x = CGAL::to_double(p.x());
      const double y = CGAL::to_double(p.y());
      const double z = CGAL::to_double(p.z());

      const double minx = CGAL::to_double(min_x_);
      const double miny = CGAL::to_double(min_y_);
      const double minz = CGAL::to_double(min_z_);
      const double maxx = CGAL::to_double(max_x_);
      const double maxy = CGAL::to_double(max_y_);
      const double maxz = CGAL::to_double(max_z_);

      if (x < minx || x > maxx ||
          y < miny || y > maxy ||
          z < minz || z > maxz) {
        return false;
      }

      ix = static_cast<int>(std::floor((x - minx) / h));
      iy = static_cast<int>(std::floor((y - miny) / h));
      iz = static_cast<int>(std::floor((z - minz) / h));

      // Check that the indices are within the valid range.
      CGAL_assertion(ix >= 0 && iy >= 0 && iz >= 0);
      CGAL_assertion(static_cast<std::size_t>(ix) <= nx_);
      CGAL_assertion(static_cast<std::size_t>(iy) <= ny_);
      CGAL_assertion(static_cast<std::size_t>(iz) <= nz_);

      ix = (std::min)(ix, static_cast<int>(nx_) - 1);
      iy = (std::min)(iy, static_cast<int>(ny_) - 1);
      iz = (std::min)(iz, static_cast<int>(nz_) - 1);

      return valid_coordinates(ix, iy, iz);
    }

    /**
     * @brief Converts 3D grid coordinates into a flat array index.
     *
     * @param ix Cell index along x.
     * @param iy Cell index along y.
     * @param iz Cell index along z.
     *
     * @return Linearized index into the cell array.
    */
    std::size_t flat_index(int ix, int iy, int iz) const {
      return (static_cast<std::size_t>(ix) * ny_
            + static_cast<std::size_t>(iy)) * nz_
            + static_cast<std::size_t>(iz);
    }

    /**
     * @brief Checks whether integer grid coordinates are valid.
     *
     * @param ix Cell index along x.
     * @param iy Cell index along y.
     * @param iz Cell index along z.
     *
     * @return `true` if the coordinates lie inside the grid.
     */
    bool valid_coordinates(int ix, int iy, int iz) const {
      return ix >= 0 && iy >= 0 && iz >= 0 &&
            static_cast<std::size_t>(ix) < nx_ &&
            static_cast<std::size_t>(iy) < ny_ &&
            static_cast<std::size_t>(iz) < nz_;
    }

    /// @brief Returns the total number of grid cells.
    std::size_t number_of_cells() const {
      return nx_ * ny_ * nz_;
    }

    private:
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
     * @brief Allocates the 3D cell array from the current bounds and box size.
     */
    void initialize_grid() {
      const double h = CGAL::to_double(box_size_);
      const double extent_x = CGAL::to_double(max_x_ - min_x_);
      const double extent_y = CGAL::to_double(max_y_ - min_y_);
      const double extent_z = CGAL::to_double(max_z_ - min_z_);

      nx_ = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(extent_x / h)));
      ny_ = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(extent_y / h)));
      nz_ = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(extent_z / h)));

      cells_.clear();
      cells_.resize(number_of_cells());
    }

    /**
     * @brief Inserts a single point into the grid.
     *
     * @param point_id Index of the point in the input cloud.
     * @param p        Point position.
     */
    void insert(Index point_id, const Point_3& p) {
      int ix, iy, iz;
      if (!to_grid_coordinates(p, ix, iy, iz)) {
        return;
      }

      // Always add the point to the cell it falls into.
      const std::size_t idx = flat_index(ix, iy, iz);
      cells_[idx].add_point_contained(point_id);

      // And add it to all neighboring cells within distance d.
      for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
          for (int dz = -1; dz <= 1; ++dz) {
            int nx = ix + dx;
            int ny = iy + dy;
            int nz = iz + dz;
            if (valid_coordinates(nx, ny, nz) && distance_to_box(p, nx, ny, nz) <= box_size_) {
              const std::size_t nidx = flat_index(nx, ny, nz);
              cells_[nidx].add_point_within_distance_d(point_id);
            }
          }
        }
      }
    }

    FT distance_to_box(const Point_3& p, int ix, int iy, int iz) const {
      const Point_3 center = cell_center(ix, iy, iz);
      const Vector_3 diff = p - center;
      const FT half_box = box_size_ * FT(0.5);

      FT dx = std::max(FT(0), CGAL::abs(diff.x()) - half_box);
      FT dy = std::max(FT(0), CGAL::abs(diff.y()) - half_box);
      FT dz = std::max(FT(0), CGAL::abs(diff.z()) - half_box);

      return CGAL::sqrt(dx*dx + dy*dy + dz*dz);
    }

    /**
    * @brief Returns the center of a cell.
    */
    Point_3 cell_center(int ix, int iy, int iz) const {
      const FT x = min_x_ + (FT(ix) + FT(0.5)) * box_size_;
      const FT y = min_y_ + (FT(iy) + FT(0.5)) * box_size_;
      const FT z = min_z_ + (FT(iz) + FT(0.5)) * box_size_;
      return Point_3(x, y, z);
    }

    /**
     * @brief Returns the normalized average normal of a cell.
     */
    Vector_3 compute_cell_normal(int ix, int iy, int iz) const {
      const Cell& c = cells_[flat_index(ix, iy, iz)];
      if (c.normal_count == 0) {
        std::cout << "Warning: cell (" << ix << "," << iy << "," << iz << ") has no normals." << std::endl;
        return CGAL::NULL_VECTOR;
      }

      const FT inv = FT(1) / FT(c.normal_count);
      Vector_3 n = c.normal_sum * inv;

      const double len2 = CGAL::to_double(n.squared_length());
      if (len2 <= 0.0) {
        return CGAL::NULL_VECTOR;
      }

      return n * (1.0 / std::sqrt(len2));
    }

    const std::vector<Index>& point_ids_within_distance_d_of_box(std::size_t ix, std::size_t iy, std::size_t iz) const {
      CGAL_assertion(ix < nx_ && iy < ny_ && iz < nz_);
      const Cell* c = cell(ix, iy, iz);
      CGAL_assertion(c != nullptr);

      return c->point_ids_within_distance_d;
    }

    const std::vector<Index>& point_ids_in_box(std::size_t ix, std::size_t iy, std::size_t iz) const {
      CGAL_assertion(ix < nx_ && iy < ny_ && iz < nz_);
      const Cell* c = cell(ix, iy, iz);
      CGAL_assertion(c != nullptr);

      return c->point_ids_contained;
    }

    FT box_size_;
    FT min_x_ = 0, max_x_ = 0;
    FT min_y_ = 0, max_y_ = 0;
    FT min_z_ = 0, max_z_ = 0;

    std::size_t nx_ = 0;
    std::size_t ny_ = 0;
    std::size_t nz_ = 0;

    std::vector<Cell> cells_;
    std::vector<Point_3> points_;
    std::vector<Vector_3> normals_;
    std::vector<FT> splat_sizes_;
    std::vector<Vector_3> block_normals_;
    FT max_splat_radius_ = 0;

  //   ---------------------------------- DEBUG ---------------------------------------
  public:

  void set_all_splat_sizes(FT size) {
    splat_sizes_.assign(points_.size(), size);
    max_splat_radius_ = size;
  }

  void write_splat_mesh_obj(const char * filename) const {

    std::cout << "bbox: " << min_x_ << "," << min_y_ << "," << min_z_ << " - " << max_x_ << "," << max_y_ << "," << max_z_ << std::endl;

    typename CGAL::Surface_mesh<Point_3> mesh;

    for (std::size_t i = 0; i < points_.size(); ++i) {
      FT global_spacing = 1.05*box_size_;

      FT radius = splat_sizes_.size() == 0 ? global_spacing : splat_sizes_[i];

      add_splat(mesh, points_[i], normals_[i], radius);
    }

    CGAL::IO::write_polygon_mesh(filename, mesh);

    std::cout << "Splat mesh written to " << filename << std::endl;
  }

  Vector_3 normalize_vec(const Vector_3& n) const {
    const double len2 = CGAL::to_double(n.squared_length());
    if (len2 <= 0.0) {
      return CGAL::NULL_VECTOR;
    }
    const Vector_3 normalized = n / std::sqrt(len2);
    CGAL_assertion(std::abs(normalized.squared_length() - FT(1.0)) < 1e-8);
    return normalized;
  }

  void add_splat(CGAL::Surface_mesh<Point_3>& m,
               const Point_3& center,
               const Vector_3& normal,
               FT radius) const
  {
      const int num_vertices = 10;

      std::vector<Vector_3> tangent_frame_vectors = compute_local_tangent_frame(normal);

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

  
  //   public:
  //   /**
  //    * @brief Writes the input point cloud and normals to a PLY file.
  //    *
  //    * @param filename Output PLY filename.
  //    *
  //    * @return `true` if the file was written successfully.
  //    */
  //   bool write_point_cloud_ply(const std::string& filename) const
  //   {
  //     std::ofstream out(filename);
  //     if (!out) {
  //       std::cerr << "Error: cannot open " << filename << " for writing.\n";
  //       return false;
  //     }

  //     out << "ply\n";
  //     out << "format ascii 1.0\n";
  //     out << "element vertex " << points_.size() << "\n";
  //     out << "property float x\n";
  //     out << "property float y\n";
  //     out << "property float z\n";
  //     out << "property float nx\n";
  //     out << "property float ny\n";
  //     out << "property float nz\n";
  //     out << "end_header\n";

  //     for (std::size_t i = 0; i < points_.size(); ++i) {
  //       const Point_3& p = points_[i];
  //       const Vector_3 n = (i < normals_.size()) ? normals_[i] : CGAL::NULL_VECTOR;

  //       out << std::setprecision(17)
  //           << CGAL::to_double(p.x()) << ' '
  //           << CGAL::to_double(p.y()) << ' '
  //           << CGAL::to_double(p.z()) << ' '
  //           << CGAL::to_double(n.x()) << ' '
  //           << CGAL::to_double(n.y()) << ' '
  //           << CGAL::to_double(n.z()) << '\n';
  //     }

  //     return true;
  //   }

  //   /**
  //    * @brief Writes all grid vertices to a PLY file for visualization.
  //    *
  //    * @param filename Output PLY filename.
  //    *
  //    * @return `true` if the file was written successfully.
  //    */
  //   bool write_grid_vertices_ply(const std::string& filename) const
  //   {
  //     std::ofstream out(filename);
  //     if (!out) {
  //       std::cerr << "Error: cannot open " << filename << " for writing.\n";
  //       return false;
  //     }

  //     std::vector<Point_3> vertices;
  //     vertices.reserve((nx_ + 1) * (ny_ + 1) * (nz_ + 1));

  //     for (std::size_t ix = 0; ix <= nx_; ++ix) {
  //       for (std::size_t iy = 0; iy <= ny_; ++iy) {
  //         for (std::size_t iz = 0; iz <= nz_; ++iz) {
  //           const FT x = min_x_ + FT(ix) * box_size_;
  //           const FT y = min_y_ + FT(iy) * box_size_;
  //           const FT z = min_z_ + FT(iz) * box_size_;
  //           vertices.emplace_back(x, y, z);
  //         }
  //       }
  //     }

  //     out << "ply\n";
  //     out << "format ascii 1.0\n";
  //     out << "element vertex " << vertices.size() << "\n";
  //     out << "property float x\n";
  //     out << "property float y\n";
  //     out << "property float z\n";
  //     out << "end_header\n";

  //     for (const Point_3& p : vertices) {
  //       out << std::setprecision(17)
  //           << CGAL::to_double(p.x()) << ' '
  //           << CGAL::to_double(p.y()) << ' '
  //           << CGAL::to_double(p.z()) << '\n';
  //     }

  //     return true;
  //   }

  //   /**
  //    * @brief Writes occupied cell centers and their averaged normals to a PLY file.
  //    *
  //    * @param filename Output PLY filename.
  //    * @param normal_scale Scale factor for visualizing the normal direction.
  //    *
  //    * @return `true` if the file was written successfully.
  //    */
  //   bool write_cell_centers_and_normals_ply(const std::string& filename,
  //                                           double normal_scale = 0.25) const
  //   {
  //     std::ofstream out(filename);
  //     if (!out) {
  //       std::cerr << "Error: cannot open " << filename << " for writing.\n";
  //       return false;
  //     }

  //     std::vector<Point_3> centers;
  //     std::vector<Vector_3> normals;
  //     centers.reserve(cells_.size());
  //     normals.reserve(cells_.size());

  //     for (std::size_t ix = 0; ix < nx_; ++ix) {
  //       for (std::size_t iy = 0; iy < ny_; ++iy) {
  //         for (std::size_t iz = 0; iz < nz_; ++iz) {
  //           const Cell& c = cells_[flat_index(static_cast<int>(ix),
  //                                             static_cast<int>(iy),
  //                                             static_cast<int>(iz))];

  //           if (c.point_ids.empty()) {
  //             continue;
  //           }

  //           centers.push_back(cell_center(static_cast<int>(ix),
  //                                         static_cast<int>(iy),
  //                                         static_cast<int>(iz)));
  //           normals.push_back(compute_cell_normal(static_cast<int>(ix),
  //                                                 static_cast<int>(iy),
  //                                                 static_cast<int>(iz)));
  //         }
  //       }
  //     }

  //     out << "ply\n";
  //     out << "format ascii 1.0\n";
  //     out << "element vertex " << centers.size() * 2 << "\n";
  //     out << "property float x\n";
  //     out << "property float y\n";
  //     out << "property float z\n";
  //     out << "element edge " << centers.size() << "\n";
  //     out << "property int vertex1\n";
  //     out << "property int vertex2\n";
  //     out << "end_header\n";

  //     for (std::size_t i = 0; i < centers.size(); ++i) {
  //       const Point_3& c = centers[i];
  //       const Vector_3 n = normals[i];

  //       const Point_3 tip(
  //         c.x() + FT(normal_scale) * n.x(),
  //         c.y() + FT(normal_scale) * n.y(),
  //         c.z() + FT(normal_scale) * n.z()
  //       );

  //       out << std::setprecision(17)
  //           << CGAL::to_double(c.x()) << ' '
  //           << CGAL::to_double(c.y()) << ' '
  //           << CGAL::to_double(c.z()) << '\n';

  //       out << std::setprecision(17)
  //           << CGAL::to_double(tip.x()) << ' '
  //           << CGAL::to_double(tip.y()) << ' '
  //           << CGAL::to_double(tip.z()) << '\n';
  //     }

  //     for (std::size_t i = 0; i < centers.size(); ++i) {
  //       const int v0 = static_cast<int>(2 * i);
  //       const int v1 = static_cast<int>(2 * i + 1);
  //       out << v0 << ' ' << v1 << '\n';
  //     }

  //     return true;
  //   }
  //   ---------------------------------------------------------------------------

  };

  /**
   * @ingroup PkgSplatSurfaceReconstruction3Classes
   *
   * @brief Incrementally reconstructs a polygon mesh from splat candidates.
   *
   * The reconstruction grows a halfedge front from an initial seed. Candidate
   * vertices are generated from pairs of existing mesh vertices and accepted
   * according to proximity, projection, priority, and local topology tests.
   *
   * During growth, the class maintains an explicit halfedge connectivity graph.
   * Faces are created only after the growth phase by traversing the resulting
   * `next()` cycles.
   *
   * @tparam PointRange  Input point range type.
   * @tparam NormalRange Input normal range type.
   * @tparam TriangleMesh Output mesh type model of `HalfedgeListGraph`, `FaceListGraph`, and `MutableFaceGraph`
  */
  template <typename PointRange, typename NormalRange, typename TriangleMesh>
  class Splat_surface_reconstruction_3 {
    private:
      friend struct Splat_surface_reconstruction_3_test_access<PointRange, NormalRange, TriangleMesh>;

    public:
      using Point_3 = typename PointRange::value_type;
      using Kernel = typename Kernel_traits<Point_3>::Kernel;
      using Vector_3 = typename Kernel::Vector_3;
      using FT = typename Kernel::FT;
      using Grid = Box_grid<Kernel>;
      using Index = typename Grid::Index;

    private:
      using Vector_2 = typename Kernel::Vector_2;
      using Point_2  = typename Kernel::Point_2;

      using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
      using edge_descriptor = typename boost::graph_traits<TriangleMesh>::edge_descriptor;
      using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
      using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;

      struct Candidate {
      Point_3 position;
      Vector_3 normal;
      Index splat_id;
      vertex_descriptor first;
      vertex_descriptor second;
      int priority = 0;

      Candidate(const Point_3& position,
                const Vector_3& normal,
                Index splat_id,
                vertex_descriptor first,
                vertex_descriptor second)
        : position(position),
          normal(normal),
          splat_id(splat_id),
          first(first),
          second(second) {}
      };

    public:
      /**
       * @brief Constructs a reconstruction object and initializes mesh properties.
       *
       * Adds the per-vertex normal property map and per-edge property map,
       * initializes the spatial lookup structure, and creates the initial seed.
       *
       * @param grid Spatial grid containing input points, normals, and splat sizes.
       * @param output_mesh Mesh receiving the reconstructed surface.
      */
      Splat_surface_reconstruction_3(const Grid& grid, TriangleMesh& output_mesh)
      : grid_(grid),
        mesh_(output_mesh),
        points_pm_(get(CGAL::vertex_point, mesh_)) {

        bool created = false;
        std::tie(normals_pm_, created) = mesh_.template add_property_map<vertex_descriptor, Vector_3>("v:normal", CGAL::NULL_VECTOR);

        bool created_w = false;
        std::tie(edge_weight_pm_, created_w) = mesh_.template add_property_map<edge_descriptor, FT>("e:weight", FT(1));

        if (created) {
          std::cout << "Created vertex normal property map." << std::endl;
        }

        mesh_vertices_per_cell_.resize(grid_.number_of_cells());

        seed_from_grid();
      }

      /**
       * @brief Executes the complete reconstruction pipeline.
       *
       * Candidates are processed by priority until no candidates remain. Accepted
       * vertices are inserted into the halfedge graph and new candidates are
       * generated around them. Once growth stops, open halfedge cycles are converted
       * into polygonal faces.
      */
      void run() {
        if (!seeded_) {
          std::cerr << "No initial seed was found. Mesh growth did not start.\n";
          return;
        }
        else {
          std::cout << "Initial seed vertices created. Starting mesh growth...\n";
        }

        std::size_t rejected_near = 0;
        std::size_t rejected_proj = 0;
        std::size_t accepted = 0;

        int iterationIndex = 0;

        while (!candidate_queues_empty()) {
          Candidate cand(Point_3(0,0,0), CGAL::NULL_VECTOR, 0, vertex_descriptor(), vertex_descriptor()); // Dummy initialization
          if (!pop_candidate_from_queue(cand)) {
            break;
          }

          std::size_t target_parent_1 = 321;
          std::size_t target_parent_2 = 349;

          bool is_target = (cand.first.idx() == target_parent_1 || cand.first.idx() == target_parent_2)
                           && (cand.second.idx() == target_parent_1 || cand.second.idx() == target_parent_2);

          if (is_target) {
            std::cout << "Involved halfedge: " << halfedge(cand.first, cand.second, mesh_).first << "\n";
            std::cout << "Candidate involves target parents: " << cand.first << ", " << cand.second << std::endl;
            export_halfedge_mesh_json(mesh_, "involved_halfedge.json");
          }

          if (has_vertex_near(cand.position)) {
            rejected_near++;
            //std::cout << "Candidate rejected due to nearby vertex: " << cand.first << ", " << cand.second << std::endl;
            continue;
          }

          if (!projection_check(cand, grid_.box_normal(get(points_pm_, cand.first)))
           || !projection_check(cand, grid_.box_normal(get(points_pm_, cand.second)))
           || !projection_check(cand, grid_.box_normal(cand.position)) ) {
            rejected_proj++;
            //std::cout << "Candidate rejected due to projection check: " << cand.first << ", " << cand.second << std::endl;
            if (is_target) {
              CGAL_assertion(false);
            }
            continue;
          }

          int new_priority = compute_priority(cand);
          if (new_priority != cand.priority) {
            cand.priority = new_priority;
            push_candidate_in_queue(cand);
            continue;
          }

          vertex_descriptor nv = mesh_.add_vertex(cand.position);
          put(points_pm_, nv, cand.position);
          put(normals_pm_, nv, cand.normal);

          accepted++;
          if (!build_graph(cand.first, cand.second, nv, cand.normal)) {
            remove_vertex(nv, mesh_);
            continue;
          }

          insert_mesh_vertex(nv);

          push_candidates_from_vertex(nv);

          if (num_vertices(mesh_) % 1000 == 0) { // Print progress every 1000 vertices
            std::cout << num_vertices(mesh_) << " vertices in mesh." << "\n";
          }

          iterationIndex++;
        }

        assert_valid_triangle_loops();

        // Every half edge must have an opposite half edge (i.e., the mesh should be manifold).
        const auto null_h = boost::graph_traits<TriangleMesh>::null_halfedge();
        const auto null_v = boost::graph_traits<TriangleMesh>::null_vertex();
        for (const auto& he : halfedges(mesh_)) {
          CGAL_assertion(opposite(he, mesh_) != null_h);
          const auto u = source(he, mesh_);
          const auto v = target(he, mesh_);
          CGAL_assertion(u != null_v);
          CGAL_assertion(v != null_v);
          
          const auto [h_uv, found_uv] = halfedge(u, v, mesh_);
          const auto [h_vu, found_vu] = halfedge(v, u, mesh_);

          CGAL_assertion(found_uv && found_vu);
        }

        std::cout<<"Setting Faces..." << std::endl;
        fill_faces_from_next_cycles();

        assert_valid_triangle_faces();

        assert_final_mesh_manifold();

        std::cout << "Final mesh size: " << num_vertices(mesh_) << " vertices, " << num_edges(mesh_) << " edges, " << num_faces(mesh_) << " faces.\n";
        std::cout << "  Accepted: " << accepted << ", Rejected (near): " << rejected_near << ", Rejected (proj): " << rejected_proj << "\n";
      }

    private:

      bool is_halfedge_inside(halfedge_descriptor h, vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2) const {
        Vector_3 normal = normals_pm_[v0] + normals_pm_[v1] + normals_pm_[v2];
        normal = normal / std::sqrt(normal.squared_length());

        // Find the orientation of the triangle formed by the vertices with respect to the normal
        const auto& p0 = get(points_pm_, v0);
        const auto& p1 = get(points_pm_, v1);
        const auto& p2 = get(points_pm_, v2);

        return CGAL::orientation(p0, p1, p2, CGAL::ORIGIN + normal) == CGAL::POSITIVE;
      }

      // Asserts that all triangle faces in the mesh have valid halfedge connections
      void assert_valid_triangle_faces() const {
        size_t num_halfedges = halfedges(mesh_).size();

        std::vector<bool> halfedge_visited(num_halfedges, false);
        
        for (halfedge_descriptor h0 : halfedges(mesh_)) {
          const auto h0_target = target(h0, mesh_);
          const auto h0_source = source(h0, mesh_);

          for (halfedge_descriptor h1 : halfedges_around_source(h0_target, mesh_)) {

            // Not a halfedge that can complete the triangle
            if (source(h1, mesh_) != h0_target) {
              continue;
            }

            const auto h1_target = target(h1, mesh_);

            for (halfedge_descriptor h2 : halfedges_around_source(h1_target, mesh_)) {
              
              // Not a halfedge that can complete the triangle
              if (source(h2, mesh_) != h1_target) {
                continue;
              }

              const auto h2_target = target(h2, mesh_);

              if (!is_halfedge_inside(h0, h0_target, h1_target, h2_target) ||
                  !is_halfedge_inside(h1, h1_target, h2_target, h0_target) ||
                  !is_halfedge_inside(h2, h2_target, h0_target, h1_target)) {
                continue;
              }
              
              if (h2_target != h0_source) {
                continue;
              }

              bool loop_correct = (next(h0, mesh_) == h1 && next(h1, mesh_) == h2 && next(h2, mesh_) == h0);

              // Check that the edges are oriented correctly
              Point_3 p_a = points_pm_[h0_source];
              Point_3 p_b = points_pm_[h0_target];
              Point_3 p_c = points_pm_[h1_target];

              Vector_3 normal_a = normals_pm_[h0_source];
              Vector_3 normal_b = normals_pm_[h0_target];
              Vector_3 normal_c = normals_pm_[h1_target];
              Vector_3 final_normal = normal_a + normal_b + normal_c;
              final_normal = final_normal / std::sqrt(final_normal.squared_length());

              CGAL::Orientation orientation = orient3(p_a, p_b, p_c, final_normal);

              bool orientation_correct = (orientation == CGAL::LEFT_TURN);

              if (orientation_correct != loop_correct) {
                export_halfedge_mesh_json(mesh_, "invalid_orientation.json");
                std::cout << "Invalid orientation detected for triangle: h0=" << h0 << ", h1=" << h1 << ", h2=" << h2 << std::endl;
                std::cout << h0 << ": " << h0_source << "->" << h0_target << " (next: " << next(h0, mesh_) << ")" << std::endl;
                std::cout << h1 << ": " << source(h1, mesh_) << "->" << target(h1, mesh_) << " (next: " << next(h1, mesh_) << ")" << std::endl;
                std::cout << h2 << ": " << source(h2, mesh_) << "->" << target(h2, mesh_) << " (next: " << next(h2, mesh_) << ")" << std::endl;
                std::cout << "Orientation: " << orientation << std::endl;
                std::cout << "Points: p_a=" << p_a << ", p_b=" << p_b << ", p_c=" << p_c << std::endl;
                std::cout << "Final normal: " << final_normal << std::endl;
              }

              // If the orientation does not match the loop correctness, report an error.
              // A loop can have incorrect orientation but only if the loop is also incorrect.
              // However, if the loop is incorrect, this does not imply that the orientation is incorrect.
              CGAL_assertion(implies(loop_correct, orientation_correct));
            }
          }
        }
        std::cout << "Finished checking all triangle faces." << std::endl;
      }

      bool implies(bool b1, bool b2) const {
        return !b1 || b2;
      }

      void assert_valid_triangle_loops() const {
        size_t num_halfedges = halfedges(mesh_).size();

        std::vector<bool> halfedge_visited(num_halfedges, false);
        
        for (halfedge_descriptor h0 : halfedges(mesh_)) {
          const auto h0_target = target(h0, mesh_);
          const auto h0_source = source(h0, mesh_);

          for (halfedge_descriptor h1 : halfedges_around_source(h0_target, mesh_)) {

            // Not a halfedge that can complete the triangle
            if (source(h1, mesh_) != h0_target) {
              continue;
            }

            const auto h1_target = target(h1, mesh_);

            for (halfedge_descriptor h2 : halfedges_around_source(h1_target, mesh_)) {
              
              // Not a halfedge that can complete the triangle
              if (source(h2, mesh_) != h1_target) {
                continue;
              }

              const auto h2_target = target(h2, mesh_);

              if (!is_halfedge_inside(h0, h0_target, h1_target, h2_target) ||
                  !is_halfedge_inside(h1, h1_target, h2_target, h0_target) ||
                  !is_halfedge_inside(h2, h2_target, h0_target, h1_target)) {
                continue;
              }
              
              if (h2_target != h0_source) {
                continue;
              }

              bool loop_correct = (next(h0, mesh_) == h1 && next(h1, mesh_) == h2 && next(h2, mesh_) == h0);

              // Check that the edges are oriented correctly
              Point_3 p_a = points_pm_[h0_source];
              Point_3 p_b = points_pm_[h0_target];
              Point_3 p_c = points_pm_[h1_target];

              Vector_3 normal_a = normals_pm_[h0_source];
              Vector_3 normal_b = normals_pm_[h0_target];
              Vector_3 normal_c = normals_pm_[h1_target];
              Vector_3 final_normal = normal_a + normal_b + normal_c;
              final_normal = final_normal / std::sqrt(final_normal.squared_length());

              CGAL::Orientation orientation = orient3(p_a, p_b, p_c, final_normal);

              bool orientation_correct = (orientation == CGAL::LEFT_TURN);

              // if (orientation_correct != loop_correct) {
              //   export_halfedge_mesh_json(mesh_, "invalid_orientation.json");
              //   std::cout << "Invalid orientation detected for triangle: h0=" << h0 << ", h1=" << h1 << ", h2=" << h2 << std::endl;
              //   std::cout << h0 << ": " << h0_source << "->" << h0_target << " (next: " << next(h0, mesh_) << ")" << std::endl;
              //   std::cout << h1 << ": " << source(h1, mesh_) << "->" << target(h1, mesh_) << " (next: " << next(h1, mesh_) << ")" << std::endl;
              //   std::cout << h2 << ": " << source(h2, mesh_) << "->" << target(h2, mesh_) << " (next: " << next(h2, mesh_) << ")" << std::endl;
              //   std::cout << "Orientation: " << orientation << std::endl;
              //   std::cout << "Points: p_a=" << p_a << ", p_b=" << p_b << ", p_c=" << p_c << std::endl;
              //   std::cout << "Final normal: " << final_normal << std::endl;
              // }

              // The loop being correct implies that the orientation should also be correct.
              // If the loop is not correct, the orientation may or may not be correct.
              CGAL_assertion(implies(loop_correct, orientation_correct));
            }
          }
        }
        std::cout << "Finished checking all triangle loops." << std::endl;
      }

      void assert_final_mesh_manifold() const {
        // 1. Basic structural validity.
        CGAL_assertion(CGAL::is_valid_polygon_mesh(mesh_));

        // 2. The output should consist entirely of triangles.
        CGAL_assertion(CGAL::is_triangle_mesh(mesh_));

        // 3. No isolated vertices.
        for (const auto v : vertices(mesh_))
        {
          CGAL_assertion(halfedge(v, mesh_) != boost::graph_traits<TriangleMesh>
                                        ::null_halfedge());
        }

        // 4. No wire edges (edges with zero incident faces).
        for (const auto e : edges(mesh_))
        {
          const auto h = halfedge(e, mesh_);
          const bool has_face_1 = face(h, mesh_) !=
                                  boost::graph_traits<TriangleMesh>::null_face();
          const bool has_face_2 = face(opposite(h, mesh_), mesh_) !=
                                  boost::graph_traits<TriangleMesh>::null_face();

          CGAL_assertion(has_face_1 || has_face_2);
        }

        // All halfedges must be set correctly
        for (const auto h : halfedges(mesh_))
        {
          const auto h_next = next(h, mesh_);
          const auto h_prev = next(h_next, mesh_);
          CGAL_assertion(opposite(opposite(h, mesh_), mesh_) == h);

          const auto h_source = source(h, mesh_);
          const auto h_target = target(h, mesh_);
          const auto h_next_source = source(h_next, mesh_);

          CGAL_assertion(h_source == target(opposite(h, mesh_), mesh_));
          CGAL_assertion(h_target == source(opposite(h, mesh_), mesh_));

          CGAL_assertion(h_target == h_next_source);
        }

        std::size_t non_manifold_edges = 0;

        for (const auto h : halfedges(mesh_))
        {
          const auto a = source(h, mesh_);
          const auto b = target(h, mesh_);

          std::vector<std::size_t> incident_faces;

          // Find all faces containing the directed edge a -> b.
          for (const auto f : faces(mesh_))
          {
            auto hf = halfedge(f, mesh_);

            if (hf == TriangleMesh::null_halfedge())
              continue;

            const auto h0 = hf;
            std::size_t steps = 0;
            const std::size_t max_steps = num_halfedges(mesh_);

            do
            {
              if (steps++ >= max_steps)
              {
                std::cout << "Exceeded max steps in face cycle for face " << f.idx() << std::endl;
                CGAL_assertion(false);
                break;
              }

              if (source(hf, mesh_) == a && target(hf, mesh_) == b)
              {
                incident_faces.push_back(f.idx());
                break;
              }

              hf = next(hf, mesh_);

            } while (hf != h0);
          }

          if (incident_faces.size() > 2)
          {
            ++non_manifold_edges;

            std::cerr << "Directed edge "
                      << a.idx() << " -> " << b.idx()
                      << " belongs to " << incident_faces.size()
                      << " faces: ";

            for (const auto fi : incident_faces)
              std::cerr << fi << ' ';

            std::cerr << '\n';
          }
        }

        std::cerr << "Directed edges with more than two incident faces: "
                  << non_manifold_edges << '\n';

        assert(non_manifold_edges == 0 &&
              "An edge has more than two incident faces");
      }

      /**
       * @brief Inserts a mesh vertex into the spatial lookup structure.
       *
       * The vertex is stored in the grid cell containing its position.
      */
      void insert_mesh_vertex(vertex_descriptor vd) {
        int ix, iy, iz;

        const Point_3& p = get(points_pm_, vd);

        if (!grid_.to_grid_coordinates(p, ix, iy, iz))
            return;

        mesh_vertices_per_cell_[grid_.flat_index(ix, iy, iz)].push_back(vd);
      }

      /**
       * @brief Returns mesh vertices within a spatial neighborhood.
       *
       * @param p Query position.
       * @param radius Search radius.
       *
       * @return Mesh vertices in neighboring grid cells.
      */
      std::vector<vertex_descriptor> nearby_mesh_vertices(const Point_3& p, FT radius) const {
        std::vector<vertex_descriptor> out;

        int cx, cy, cz;

        if (!grid_.to_grid_coordinates(p, cx, cy, cz))
          return out;

        const int r = (std::max)(
          1,
          static_cast<int>(std::ceil(CGAL::to_double(radius) / CGAL::to_double(grid_.get_box_size()))));

        for (int dx = -r; dx <= r; ++dx)
          for (int dy = -r; dy <= r; ++dy)
            for (int dz = -r; dz <= r; ++dz) {
              if (!grid_.valid_coordinates(cx + dx, cy + dy, cz + dz))
                continue;

              const auto& cell = mesh_vertices_per_cell_[grid_.flat_index(cx + dx, cy + dy, cz + dz)];

              out.insert(out.end(), cell.begin(), cell.end());
            }

        return out;
      }

      /**
       * @brief Inserts a candidate into its priority queue.
       *
       * The candidate priority is clamped to the valid queue range before insertion.
      */
      void push_candidate_in_queue(Candidate cand) {
        CGAL_assertion(cand.priority >= 0);
        CGAL_assertion(cand.priority < NUM_PRIORITIES);

        candidate_queues_[cand.priority].push_back(cand);
      }

      /**
       * @brief Removes the highest-priority available candidate from the queues.
       *
       * @param cand Output candidate.
       *
       * @return `true` if a candidate was available.
      */
      bool pop_candidate_from_queue(Candidate& cand) {
        for (int p = NUM_PRIORITIES - 1; p >= 0; --p) {
          if (!candidate_queues_[p].empty()) {
            cand = candidate_queues_[p].front();
            candidate_queues_[p].pop_front();
            return true;
          }
        }

        return false;
      }

      /// @brief Returns whether every candidate queue is empty.
      bool candidate_queues_empty() const {
        for (const auto& q : candidate_queues_) {
          if (!q.empty())
            return false;
        }

        return true;
      }

      /**
       * @brief Initializes the mesh with the initial two-vertex seed.
       *
       * The seed creates a single undirected edge whose two halfedges form the
       * initial front. Candidate generation is then started from both seed vertices.
      */
      void seed_from_grid() {
        std::optional<typename Grid::Initial_seed> seed = grid_.get_initial_seed();
        if (!seed) {
          std::cerr << "No initial splat seed found.\n";
          return;
        }

        std::cout << "Adding initial seed vertices to mesh...\n";
        // v0_ = mesh_.add_vertex(seed->first);
        vertex_descriptor v0_ = add_vertex(mesh_);
        put(points_pm_, v0_, seed->first);
        put(normals_pm_, v0_, grid_.normals()[seed->splat_id]);

        vertex_descriptor v1_ = add_vertex(mesh_);
        put(points_pm_, v1_, seed->second);
        put(normals_pm_, v1_, grid_.normals()[seed->splat_id]);

        // Create the base edge between the two seed vertices.
        edge_descriptor e = add_edge(mesh_);
        halfedge_descriptor h = halfedge(e, mesh_);
        halfedge_descriptor ho = opposite(h, mesh_);

        set_target(h,  v1_, mesh_); // v0 -> v1
        set_target(ho, v0_, mesh_); // v1 -> v0
        set_halfedge(v0_, ho, mesh_);
        set_halfedge(v1_, h, mesh_);

        set_next(h, ho, mesh_);
        set_next(ho, h, mesh_);

        // Fully assert all properties of the initial vertices and edge
        CGAL_assertion(is_valid_edge_descriptor(e, mesh_));
        CGAL_assertion(is_valid_halfedge_descriptor(h, mesh_));
        CGAL_assertion(is_valid_halfedge_descriptor(ho, mesh_));
        CGAL_assertion(source(h, mesh_) == v0_);
        CGAL_assertion(source(ho, mesh_) == v1_);
        CGAL_assertion(target(h, mesh_) == v1_);
        CGAL_assertion(target(ho, mesh_) == v0_);
        CGAL_assertion(next(h, mesh_) == ho);
        CGAL_assertion(next(ho, mesh_) == h);
        CGAL_assertion(halfedge(v0_, mesh_) == ho);
        CGAL_assertion(halfedge(v1_, mesh_) == h);
        CGAL_assertion(is_valid_vertex_descriptor(v0_, mesh_));
        CGAL_assertion(is_valid_vertex_descriptor(v1_, mesh_));

        
        seeded_ = true;
        std::cout << "Mesh seeded with initial vertices. Starting growth...\n";

        insert_mesh_vertex(v0_);
        insert_mesh_vertex(v1_);

        push_candidates_from_vertex(v0_);
        push_candidates_from_vertex(v1_);

        std::cout << "seed_from_grid " << e << " added between " << v0_ << "  and " << v1_ << "\n";

      }

      /**
       * @brief Generates new candidates around a newly inserted mesh vertex.
       *
       * Nearby mesh vertices are paired with the supplied vertex and candidate
       * splat intersections are generated for each pair.
       *
       * @param nv Newly inserted mesh vertex.
      */
      void push_candidates_from_vertex(vertex_descriptor nv) {
        const Point_3 p = get(points_pm_, nv);

        auto nearby = nearby_mesh_vertices(p, FT(8.0) * grid_.get_box_size());

        for (vertex_descriptor other_v : nearby) {

          std::vector<Candidate> candidates = generate_candidates_from_parents(nv, other_v);

          for (const Candidate& c : candidates) {
            push_candidate_in_queue(c);
          }
        }
      }

      /**
       * @brief Generates candidate vertices from a pair of parent vertices.
       *
       * A construction circle derived from the parent pair is intersected with
       * nearby splats. Each valid intersection becomes a reconstruction candidate.
       *
       * @param parent_a First parent vertex.
       * @param parent_b Second parent vertex.
       *
       * @return Candidate vertices generated from the parent pair.
      */
      std::vector<Candidate> generate_candidates_from_parents(vertex_descriptor parent_a,
                                                              vertex_descriptor parent_b) const {
        std::vector<Candidate> candidates;

        const Point_3 pa = get(points_pm_, parent_a);
        const Point_3 pb = get(points_pm_, parent_b);

        const FT d = grid_.get_box_size();
        const FT d2 = CGAL::squared_distance(pa, pb);

        if (d2 < FT(0.9) * d * d) {
          return candidates;
        }

        // Circle radius from the two-parent construction.
        const FT circle_r = CGAL::sqrt((d * d) - (d2 / FT(4)));

        FT search_radius = circle_r + grid_.get_max_splat_radius();

        const Point_3 mid = midpoint(pa, pb);

        std::vector<Index> nearby = grid_.nearby_point_ids(mid, search_radius);

        const std::vector<Point_3>& points = grid_.points();
        const std::vector<Vector_3>& normals = grid_.normals();
        const std::vector<FT>& splat_sizes = grid_.splat_sizes();

        for (Index sid : nearby) {
          if (sid >= points.size() || sid >= normals.size()) {
            continue;
          }

          const FT splat_radius = splat_sizes[sid];

          const std::vector<Point_3> hits =
            grid_.intersect_circle_with_splat(pa, pb, points[sid], normals[sid], splat_radius);

          for (const Point_3& p : hits) {
            Candidate c(p, normals[sid], sid, parent_a, parent_b);
            c.priority = compute_priority(c);
            candidates.push_back(c);
          }
        }

        return candidates;
      }

      /**
       * @brief Tests whether a mesh vertex already exists near a candidate position.
       *
       * @param p Candidate position.
       *
       * @return `true` when an existing mesh vertex is closer than the reconstruction
       *         tolerance.
      */
      bool has_vertex_near(const Point_3& p) const {
        const FT tol = grid_.get_box_size() * FT(0.9);
        const FT tol2 = tol * tol;

        auto nearby = nearby_mesh_vertices(p, FT(2.0) * grid_.get_box_size());
        for (vertex_descriptor vd : nearby) {
          if (CGAL::squared_distance(get(points_pm_, vd), p) < tol2) {
            return true;
          }
        }
        return false;
      }

      /**
       * @brief Creates an oriented halfedge between two vertices.
       *
       * The new edge is initialized with the requested orientation `va -> vb`.
       * The opposite halfedge is initialized with orientation `vb -> va`.
       *
       * Vertex halfedge pointers are initialized when the corresponding vertex does
       * not already have a halfedge.
       *
       * @param va Source vertex.
       * @param vb Target vertex.
       *
       * @return The halfedge oriented from `va` to `vb`.
      */
      halfedge_descriptor create_halfedge(vertex_descriptor va,
                                          vertex_descriptor vb)
      {
        edge_descriptor e = add_edge(mesh_);
        halfedge_descriptor h_vavb = halfedge(e, mesh_);
        halfedge_descriptor h_vbva = opposite(h_vavb, mesh_);

        // h_vavb : va -> vb
        set_target(h_vavb, vb, mesh_);
        set_target(h_vbva, va, mesh_);

        if (halfedge(va, mesh_) == boost::graph_traits<TriangleMesh>::null_halfedge()) {
          set_halfedge(va, h_vbva, mesh_);
        }
        if (halfedge(vb, mesh_) == boost::graph_traits<TriangleMesh>::null_halfedge()) {
          set_halfedge(vb, h_vavb, mesh_);
        }

        return h_vavb;
      }

      // DEBUG: (TODO: remove)

      template <typename Point>
      bool export_halfedge_mesh_json(const CGAL::Surface_mesh<Point>& mesh, const std::string& filename) const
      {
        std::cout << "Start exporting halfedge mesh to " << filename << std::endl;
        using Mesh = CGAL::Surface_mesh<Point>;

        std::ofstream out(filename);
        if (!out) {
          std::cerr << "Could not open " << filename << '\n';
          return false;
        }

        out << "{\n";

        // --------------------------------------------------
        // Vertices: descriptor ID and 3D position
        // --------------------------------------------------
        out << "  \"vertices\": [\n";

        bool first = true;
        for (const auto v : mesh.vertices()) {
          if (!first) out << ",\n";
          first = false;

          const Point& p = mesh.point(v);
          const Vector_3 normal = normals_pm_[v];
          const Vector_3 box_normal = grid_.box_normal(p);

          out << "    {\"id\": " << v.idx()
              << ", \"point\": ["
              << CGAL::to_double(p.x()) << ", "
              << CGAL::to_double(p.y()) << ", "
              << CGAL::to_double(p.z()) << "], "
              << "\"normal\": ["
              << CGAL::to_double(normal.x()) << ", "
              << CGAL::to_double(normal.y()) << ", "
              << CGAL::to_double(normal.z()) << "], "
              << "\"box_normal\": ["
              << CGAL::to_double(box_normal.x()) << ", "
              << CGAL::to_double(box_normal.y()) << ", "
              << CGAL::to_double(box_normal.z()) << "]}";
        }

        out << "\n  ],\n";

        // --------------------------------------------------
        // Faces: ordered vertex IDs around each face
        // --------------------------------------------------
        out << "  \"faces\": [\n";

        first = true;
        for (const auto f : mesh.faces()) {
          if (!first) out << ",\n";
          first = false;

          out << "    {\"id\": " << f.idx() << ", \"vertices\": [";

          const auto h0 = halfedge(f, mesh);
          auto h = h0;
          bool first_vertex = true;

          do {
            if (!first_vertex) out << ", ";
            first_vertex = false;

            out << source(h, mesh).idx();
            h = next(h, mesh);
          } while (h != h0);

          out << "]}";
        }

        out << "\n  ],\n";

        // --------------------------------------------------
        // Splats
        // --------------------------------------------------
        out << "  \"splats\": [\n";

        const std::vector<FT> splat_sizes = grid_.splat_sizes();
        const std::vector<Vector_3> splat_normals = grid_.normals();
        const std::vector<Point> splat_positions = grid_.points();

        for (std::size_t i = 0; i < splat_sizes.size(); ++i) {
          if (i > 0) out << ",\n";

          out << "    {\"size\": " << CGAL::to_double(splat_sizes[i])
              << ", \"normal\": ["
              << CGAL::to_double(splat_normals[i].x()) << ", "
              << CGAL::to_double(splat_normals[i].y()) << ", "
              << CGAL::to_double(splat_normals[i].z()) << "], "
              << "\"position\": ["
              << CGAL::to_double(splat_positions[i].x()) << ", "
              << CGAL::to_double(splat_positions[i].y()) << ", "
              << CGAL::to_double(splat_positions[i].z()) << "]}";
        }
        out << "\n  ],\n";



        // --------------------------------------------------
        // Half-edges: complete directed connectivity
        // --------------------------------------------------
        out << "  \"halfedges\": [\n";

        first = true;
        for (const auto h : mesh.halfedges()) {
          if (!first) out << ",\n";
          first = false;

          const auto opp = opposite(h, mesh);
          const auto nxt = next(h, mesh);
          const auto f = face(h, mesh);

          // -1 represents the null face (a boundary half-edge).
          const long long face_id =
              (f == Mesh::null_face())
                  ? -1LL
                  : static_cast<long long>(f.idx());

          out << "    {\"id\": " << h.idx()
              << ", \"source\": " << source(h, mesh).idx()
              << ", \"target\": " << target(h, mesh).idx()
              << ", \"opposite\": " << opp.idx()
              << ", \"next\": " << nxt.idx()
              << ", \"face\": " << face_id
              << "}";
        }

        out << "\n  ]\n";
        out << "}\n";

        if (!out) {
          std::cerr << "Failed while writing " << filename << '\n';
          return false;
        }

        std::cout << "Exported mesh to " << filename << '\n';
        return true;
      }

      /**
       * @brief Inserts a new vertex into the local halfedge front.
       *
       * The two parent vertices define an existing front interval. The new vertex
       * is connected to both parents and the surrounding `next()` relationships are
       * rewired so that the old boundary interval is replaced by two new boundary
       * intervals.
       *
       * The local halfedge ordering around each parent is computed using the parent
       * box normals. The selected halfedges are first checked for compatibility
       * before any topology is modified.
       *
       * @param v0 First parent vertex.
       * @param v1 Second parent vertex.
       * @param nv Newly inserted vertex.
       * @param normal Normal associated with the new candidate.
       *
       * @return `true` if the local graph can be updated consistently, `false`
       *         otherwise.
      */
      bool build_graph(vertex_descriptor v0,
                       vertex_descriptor v1,
                       vertex_descriptor nv,
                       const Vector_3& normal) {

        const auto null_h = boost::graph_traits<TriangleMesh>::null_halfedge();
        const double two_pi = 2.0 * std::acos(-1.0);

        const Point_3 p0 = get(points_pm_, v0);
        const Point_3 p1 = get(points_pm_, v1);
        const Point_3 pn = get(points_pm_, nv);

        // The three vertices must be distinct.
        CGAL_assertion(v0 != v1);
        CGAL_assertion(v0 != nv);
        CGAL_assertion(v1 != nv);

        // The candidate triangle must not be degenerate.
        const Vector_3 e01 = p1 - p0;
        const Vector_3 e0n = pn - p0;
        CGAL_assertion(CGAL::cross_product(e01, e0n) != CGAL::NULL_VECTOR);

        // All three points should be distinct.
        CGAL_assertion(p0 != p1);
        CGAL_assertion(p0 != pn);
        CGAL_assertion(p1 != pn);

        auto normalize_vec = [](const Vector_3& n) -> Vector_3 {
          const double len2 = CGAL::to_double(n.squared_length());
          if (len2 <= 0.0) {
            return CGAL::NULL_VECTOR;
          }
          return n / std::sqrt(len2);
        };

        // Normalize the box normals before using them to construct tangent frames.
        Vector_3 n0 = normalize_vec(grid_.box_normal(p0));
        Vector_3 n1 = normalize_vec(grid_.box_normal(p1));
        Vector_3 n2 = normalize_vec(grid_.box_normal(pn));

        // // This normal can be very wrong sometimes
        // n0 = normals_pm_[v0];
        // n1 = normals_pm_[v1];

        //Vector_3 box_normal = normalize_vec(n0+n1+n2);

        // // Thus compute normal from the triangle (x impossible, v0 and v1 can be in any order)
        // Vector_3 n = normalize_vec(CGAL::cross_product(p1 - p0, pn - p0));
        // // Check that it agrees with the box normals.
        // if (n * box_normal < FT(0)) {
        //   n = -n;
        // }

        // n0 = n;
        // n1 = n;

        // CGAL_assertion(n0 != CGAL::NULL_VECTOR);
        // CGAL_assertion(n1 != CGAL::NULL_VECTOR);

        // if (n0 == CGAL::NULL_VECTOR || n1 == CGAL::NULL_VECTOR) {
        //   return false;
        // }

        // Construct an oriented tangent frame. The sign of v is adjusted so that
        // cross(u,v) agrees with the supplied normal.
        auto make_frame = [&](const Vector_3& n) -> std::pair<Vector_3, Vector_3> {
          std::vector<Vector_3> frame = grid_.compute_local_tangent_frame(n);
          Vector_3 u = frame[0];
          Vector_3 v = frame[1];

          if (u == CGAL::NULL_VECTOR || v == CGAL::NULL_VECTOR) {
            return {CGAL::NULL_VECTOR, CGAL::NULL_VECTOR};
          }

          if (CGAL::cross_product(u, v) * n < FT(0)) {
            v = -v;
          }

          CGAL_assertion(u != CGAL::NULL_VECTOR && v != CGAL::NULL_VECTOR);

          CGAL_assertion(u.squared_length() > FT(0));
          CGAL_assertion(v.squared_length() > FT(0));

          // Their cross product should agree with the normal
          CGAL_assertion(CGAL::cross_product(u, v) * n >= FT(0));

          // Check that u and v are perpendicular to the normal n.
          CGAL_assertion(
            std::abs(CGAL::to_double(u * n)) < 1e-8
          );
          CGAL_assertion(
            std::abs(CGAL::to_double(v * n)) < 1e-8
          );
          // Check that u and v are perpendicular to each other.
          CGAL_assertion(
            std::abs(CGAL::to_double(u * v)) < 1e-8
          );

          return {u, v};
        };

        // Return the polar angle of q around origin in the tangent frame.
        auto angle_of = [&](const Point_3& origin,
                            const Point_3& q,
                            const Vector_3& u,
                            const Vector_3& v) -> double {
          const Vector_3 d = q - origin;
          double a = std::atan2(CGAL::to_double(d * v), CGAL::to_double(d * u));
          if (a < 0.0) {
            a += two_pi;
          }
          CGAL_assertion(a >= 0.0);
          CGAL_assertion(a < two_pi);
          return a;
        };

        // Collect all incident halfedges around a parent and sort them by their
        // projected angle. These halfedges define the local front ordering.
        auto sorted_halfedges_around_vertex =
          [&](vertex_descriptor v, const Vector_3& local_normal)
            -> std::vector<std::pair<double, halfedge_descriptor>> {
          std::vector<std::pair<double, halfedge_descriptor>> out;
          const Point_3 pv = get(points_pm_, v);

          const auto [u, vaxis] = make_frame(local_normal);
          if (u == CGAL::NULL_VECTOR || vaxis == CGAL::NULL_VECTOR) {
            return out;
          }

          for (halfedge_descriptor h : halfedges_around_source(v, mesh_)) {
            CGAL_assertion(source(h, mesh_) == v);
            const vertex_descriptor t = target(h, mesh_);
            if (t == v) {
              continue;
            }
            out.emplace_back(angle_of(pv, get(points_pm_, t), u, vaxis), h);
          }

          std::sort(out.begin(), out.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });

          for (std::size_t i = 1; i < out.size(); ++i) {
            CGAL_assertion(out[i - 1].first <= out[i].first);
          }

          for (const auto& [angle, he] : out) {
            CGAL_assertion(he != null_h);
            CGAL_assertion(source(he, mesh_) == v);
            CGAL_assertion(angle >= 0.0 && angle < two_pi);
          }

          return out;
        };

        // Select the two front halfedges immediately clockwise and counter-clockwise
        // from the direction toward the new vertex.
        auto choose_closest_ccw_halfedges =
          [&](vertex_descriptor v,
              const Point_3& new_point,
              const Vector_3& local_normal)
            -> std::pair<halfedge_descriptor, halfedge_descriptor> {
          const Point_3 pv = get(points_pm_, v);
          const auto [u, vaxis] = make_frame(local_normal);
          if (u == CGAL::NULL_VECTOR || vaxis == CGAL::NULL_VECTOR) {
            return {null_h, null_h};
          }

          const double theta_new = angle_of(pv, new_point, u, vaxis);

          std::vector<std::pair<double, halfedge_descriptor>> sorted =
            sorted_halfedges_around_vertex(v, local_normal);

          if (sorted.empty()) {
            return {null_h, null_h};
          }

          auto it = std::lower_bound(
            sorted.begin(), sorted.end(), theta_new,
            [](const std::pair<double, halfedge_descriptor>& a, double theta) {
              return a.first < theta;
            });

          const auto& ccw = (it == sorted.end())   ? sorted.front() : *it;
          const auto& cw  = (it == sorted.begin()) ? sorted.back()  : *(it - 1);

          return {cw.second, ccw.second};
        };

        // Compute the local ordering using parent normals.
        
        auto [h1, h2] = choose_closest_ccw_halfedges(v0, pn, n0);
        auto [g2, g1] = choose_closest_ccw_halfedges(v1, pn, n1);

        // std::cout << "\n[BUILD_GRAPH ITERATION]" << std::endl;
        // std::cout << "v0: " << v0 << " v1: " << v1 << std::endl;
        // std::cout << "h1: " << h1 << " h2: " << h2 << " g1: " << g1 << " g2: " << g2 << std::endl;

        // Instead of using h2 and g2 from the local ordering, we use the edge connecting v0 and v1
        // First ensure that it exists in the mesh
        // this doesnt work in general but can be used to assert correctness in certain cases
        // const auto [h_v0v1, found_v0v1] = halfedge(v0, v1, mesh_);
        // const auto [h_v1v0, found_v1v0] = halfedge(v1, v0, mesh_);
        // CGAL_assertion(found_v0v1);
        // CGAL_assertion(found_v1v0);

        // g2 = h_v0v1;
        // h2 = h_v1v0;

        // // Then we can take h1 = opposite(next(h2))
        // h1 = opposite(next(h2, mesh_), mesh_);
        // // And g1 = previous(h2)

        // bool found_prev_h1 = false;
        // for (const auto& hi : mesh_.halfedges()) {
        //   if (next(hi, mesh_) == h2) {
        //     g1 = hi;
        //     found_prev_h1 = true;
        //     break;
        //   }
        // }
        
        // CGAL_assertion(found_prev_h1);

        // Ensure that the selected halfedges are oriented from the parent vertex to
        // the other vertex. If not, flip them to their opposite halfedge.
        if (target(h1, mesh_) != v0) h1 = opposite(h1, mesh_);
        CGAL_assertion(target(h1, mesh_) == v0);

        if (target(h2, mesh_) != v0) h2 = opposite(h2, mesh_);
        CGAL_assertion(target(h2, mesh_) == v0);

        if (target(g1, mesh_) != v1) g1 = opposite(g1, mesh_);
        CGAL_assertion(target(g1, mesh_) == v1);

        if (target(g2, mesh_) != v1) g2 = opposite(g2, mesh_);
        CGAL_assertion(target(g2, mesh_) == v1);

        // Keep only h2 and g1, and derive h1 and g2 from them.
        // h1 = opposite(next(h2, mesh_), mesh_);
        // g2 = opposite(next(g1, mesh_), mesh_);

        CGAL_assertion(h1 != null_h);
        CGAL_assertion(h2 != null_h);
        CGAL_assertion(g1 != null_h);
        CGAL_assertion(g2 != null_h);
        
        CGAL_assertion(target(h1, mesh_) == v0);
        CGAL_assertion(target(h2, mesh_) == v0);
        CGAL_assertion(target(g1, mesh_) == v1);
        CGAL_assertion(target(g2, mesh_) == v1);

        CGAL_assertion(source(h1, mesh_) != v0);
        CGAL_assertion(source(h2, mesh_) != v0);
        CGAL_assertion(source(g1, mesh_) != v1);
        CGAL_assertion(source(g2, mesh_) != v1);

        // std::cout << "after ensuring halfedge orientations:" << std::endl;
        // std::cout << "h1: " << h1 << " h2: " << h2 << " g1: " << g1 << " g2: " << g2 << std::endl;

        // If there exists an edge between v0 and v1, it should be represented by g2 and h2.
        if (halfedge(v0, v1, mesh_).second) {

          // std::cout << "halfedge between v0 and v1 exists." << std::endl;

          // std::cout << "actual: " << halfedge(v0, v1, mesh_).first << " " << halfedge(v1, v0, mesh_).first << std::endl;
          // std::cout << "g2: " << g2 << " h2: " << h2 << std::endl;
          // std::cout << "g1: " << g1 << " h1: " << h1 << std::endl;

          // It's either the pair g1, h1 or the pair g2, h2 that represents the edge between v0 and v1.
          if (g2 == halfedge(v0, v1, mesh_).first && h2 == halfedge(v1, v0, mesh_).first) {
            CGAL_assertion(g2 == opposite(h2, mesh_));
            CGAL_assertion(h2 == opposite(g2, mesh_));

            if (next(g1, mesh_) != h2) {
              std::cout << "Debugging halfedge mesh: next(h1, mesh_) != g2" << std::endl;
              std::cout << "h1: " << h1 << " next(h1, mesh_): " << next(h1, mesh_) << std::endl;
              std::cout << "g2: " << g2 << " next(g2, mesh_): " << next(g2, mesh_) << std::endl;
              std::cout << "g1: " << g1 << " next(g1, mesh_): " << next(g1, mesh_) << std::endl;
              std::cout << "h2: " << h2 << " next(h2, mesh_): " << next(h2, mesh_) << std::endl;
              export_halfedge_mesh_json(
                mesh_,
                "mesh_debug.json"
              );
            }

            // There must be continuity between g1 and h2
            CGAL_assertion(next(g1, mesh_) == h2);
          }
          else if (g1 == halfedge(v0, v1, mesh_).first && h1 == halfedge(v1, v0, mesh_).first) {
            CGAL_assertion(g1 == opposite(h1, mesh_));
            CGAL_assertion(h1 == opposite(g1, mesh_));

            if (next(h2, mesh_) != g1) {

              std::cout << "Debugging halfedge mesh: next(g2, mesh_) != h1" << std::endl;
              std::cout << "g2: " << g2 << " next(g2, mesh_): " << next(g2, mesh_) << std::endl;
              std::cout << "h1: " << h1 << " next(h1, mesh_): " << next(h1, mesh_) << std::endl;
              std::cout << "g1: " << g1 << " next(g1, mesh_): " << next(g1, mesh_) << std::endl;
              std::cout << "h2: " << h2 << " next(h2, mesh_): " << next(h2, mesh_) << std::endl;
              export_halfedge_mesh_json(
                mesh_,
                "mesh_debug.json"
              );
            }

            // There must be continuity between h2 and g1
            CGAL_assertion(next(h2, mesh_) == g1);
          }
          else {

            std::cout << "Debugging halfedge mesh: unexpected configuration" << std::endl;
            std::cout << "h1: " << h1 << " h2: " << h2 << " g1: " << g1 << " g2: " << g2 << std::endl;
            std::cout << "existing edge between v0 and v1: " << halfedge(v0, v1, mesh_).first << std::endl;
            std::cout << "existing edge between v1 and v0: " << halfedge(v1, v0, mesh_).first << std::endl;
            std::cout << "v0: " << v0 << std::endl;
            std::cout << "v1: " << v1 << std::endl;
            export_halfedge_mesh_json(
              mesh_,
              "mesh_debug.json"
            );
            CGAL_assertion(false);
          }

          // std::cout << "halfedge between v0 and v1 exists." << std::endl;

          // h2 = halfedge(v1, v0, mesh_).first;
          // g2 = halfedge(v0, v1, mesh_).first;
          // h1 = opposite(next(h2, mesh_), mesh_);
          // g1 = prev(h2, mesh_);

          // // g2 should be the edge from v0 to v1
          // CGAL_assertion(source(g2, mesh_) == v0 && target(g2, mesh_) == v1);
          // // h2 should be the edge from v1 to v0
          // CGAL_assertion(source(h2, mesh_) == v1 && target(h2, mesh_) == v0);

          // CGAL_assertion(g2 == halfedge(v0, v1, mesh_).first);
          // CGAL_assertion(h2 == halfedge(v1, v0, mesh_).first);

          // // g2 should be opposite of h2
          // CGAL_assertion(g2 == opposite(h2, mesh_));
          // CGAL_assertion(h2 == opposite(g2, mesh_));

          // // Next of g1 should be h2
          // CGAL_assertion(next(g1, mesh_) == h2);
        }

        // std::cout << "h1: " << h1 << " h2: " << h2 << " g1: " << g1 << " g2: " << g2 << std::endl;
        // std::cout << "next(h1, mesh_): " << next(h1, mesh_) << " next(h2, mesh_): " << next(h2, mesh_) << " next(g1, mesh_): " << next(g1, mesh_) << " next(g2, mesh_): " << next(g2, mesh_) << std::endl;


        // Halfedge: h1 source: vx target: v0
        // Halfedge: h2 source: v1 target: v0
        // Halfedge: g1 source: vx target: v1
        // Halfedge: g2 source: v0 target: v1

        // Halfedge: h1 source: vx target: v0
        // Halfedge: h2 source: vy target: v0
        // Halfedge: g1 source: v0 target: v1
        // Halfedge: g2 source: v1 target: v1

        // std::cout << "v0: " << v0 << " v1: " << v1 << std::endl;
        // int i = 0;
        // for (const auto& he : {h1, h2, g1, g2}) {
        //   std::cout << "Halfedge: " << he << "(" << i << ") source: " << source(he, mesh_) << " target: " << target(he, mesh_) << " next: " << next(he, mesh_) << " opposite: " << opposite(he, mesh_) << std::endl;
        //   i++;
        // }

        // std::cout << "Starting halfedge cycle check..." << std::endl;



        // Reject the candidate when the chosen halfedges aren't already connected.
        // CGAL_assertion(next(g1, mesh_) == opposite(g2, mesh_));
        // CGAL_assertion(next(h2, mesh_) == opposite(h1, mesh_));

        const auto null_f = boost::graph_traits<TriangleMesh>::null_face();

        CGAL_assertion(face(h1, mesh_) == null_f);
        CGAL_assertion(face(h2, mesh_) == null_f);
        CGAL_assertion(face(g1, mesh_) == null_f);
        CGAL_assertion(face(g2, mesh_) == null_f);

        // Check that the halfedge cycles around each parent vertex don't have repeated vertices.
        // This ensures that the local front is a simple cycle.
        // We limit the number of steps to avoid infinite loops.
        std::unordered_set<vertex_descriptor> visited;
        halfedge_descriptor h = opposite(h1, mesh_);
        bool reached_target = false;
        for (std::size_t steps = 0; steps <= 2000; ++steps) {

          CGAL_assertion(h != null_h);
          CGAL_assertion(opposite(h, mesh_) != null_h);
          CGAL_assertion(next(h, mesh_) != null_h);

          if (!visited.insert(source(h, mesh_)).second)
            return false;

          if (h == g1) {
            reached_target = true;
            break;
          }

          h = next(h, mesh_);
        }
        CGAL_assertion(reached_target);

        visited.clear();
        h = opposite(g2, mesh_);
        reached_target = false;
        for (std::size_t steps = 0; steps <= 2000; ++steps) {

          CGAL_assertion(h != null_h);
          CGAL_assertion(opposite(h, mesh_) != null_h);
          CGAL_assertion(next(h, mesh_) != null_h);

          if (!visited.insert(source(h, mesh_)).second)
            return false;

          if (h == h2) {
            reached_target = true;
            break;
          }

          h = next(h, mesh_);
        }

        CGAL_assertion(reached_target);

        // Create the two new undirected edges connecting the new vertex to its parents.
        halfedge_descriptor h_v0nv = create_halfedge(v0, nv);
        halfedge_descriptor h_nvv0 = opposite(h_v0nv, mesh_);
        halfedge_descriptor h_v1nv = create_halfedge(v1, nv);
        halfedge_descriptor h_nvv1 = opposite(h_v1nv, mesh_);

        // First boundary update:
        // g1 -> v1 -> nv -> v0 -> opposite(h1)
        set_next(g1, h_v1nv, mesh_);
        set_next(h_v1nv, h_nvv0, mesh_);
        set_next(h_nvv0, opposite(h1, mesh_), mesh_);

        // Second boundary update:
        // h2 -> v0 -> nv -> v1 -> opposite(g2)
        set_next(h2, h_v0nv, mesh_);
        set_next(h_v0nv, h_nvv1, mesh_);
        set_next(h_nvv1, opposite(g2, mesh_), mesh_);

        //assert_valid_triangle_faces();

        // Verify all six new successor relationships.
        CGAL_assertion(next(g1, mesh_) == h_v1nv);
        CGAL_assertion(next(h_v1nv, mesh_) == h_nvv0);
        CGAL_assertion(next(h_nvv0, mesh_) == opposite(h1, mesh_));

        CGAL_assertion(next(h2, mesh_) == h_v0nv);
        CGAL_assertion(next(h_v0nv, mesh_) == h_nvv1);
        CGAL_assertion(next(h_nvv1, mesh_) == opposite(g2, mesh_));

        // And verify the source and target vertices of the new halfedges.
        CGAL_assertion(source(h_v0nv, mesh_) == v0);
        CGAL_assertion(target(h_v0nv, mesh_) == nv);

        CGAL_assertion(source(h_nvv0, mesh_) == nv);
        CGAL_assertion(target(h_nvv0, mesh_) == v0);

        CGAL_assertion(source(h_v1nv, mesh_) == v1);
        CGAL_assertion(target(h_v1nv, mesh_) == nv);

        CGAL_assertion(source(h_nvv1, mesh_) == nv);
        CGAL_assertion(target(h_nvv1, mesh_) == v1);

        return true;
      }

      /**
       * @brief Converts closed `next()` cycles into polygonal faces.
       *
       * Every unassigned halfedge is used as the starting point of a boundary walk.
       * Closed cycles containing at least three halfedges are passed to the ear
       * clipping stage.
       *
       * Cycles that do not close within the safety limit are ignored.
      */
      void fill_faces_from_next_cycles()
      {
        const auto null_h = boost::graph_traits<TriangleMesh>::null_halfedge();
        const auto null_f = boost::graph_traits<TriangleMesh>::null_face();

        // Walk the current `next()` graph starting from an unassigned halfedge.
        for (halfedge_descriptor start : halfedges(mesh_)) {
          CGAL_assertion(start != null_h);

          // Only keep halfedges that are not yet assigned to a face.
          if (face(start, mesh_) != null_f) {
            continue;
          }

          std::vector<halfedge_descriptor> cycle;
          cycle.reserve(16);

          halfedge_descriptor h = start;
          const std::size_t max_steps = (std::min)((std::size_t)num_halfedges(mesh_), std::size_t(1000)); // If cycle is too long, it is considered a hole.
          bool closed = false;

          for (std::size_t i = 0; i < max_steps; ++i)
          {
            CGAL_assertion(h != null_h);

            if (face(h, mesh_) != null_f) {
              std::cout << "Encountered a halfedge that is already assigned to a face." << std::endl;
            }
            CGAL_assertion(face(h, mesh_) == null_f);

            cycle.push_back(h);

            h = next(h, mesh_);
            CGAL_assertion(h != null_h);

            if (h == start) {
              closed = true;
              break;
            }
          }

          CGAL_assertion(cycle.size() >= 3);

          // Only closed cycles with at least three halfedges can represent polygonal faces
          if (!closed || cycle.size() < 3) {
            continue;
          }

          // Verify that there are no preexisting halfedges between any non-consecutive vertices in the cycle.
          for (std::size_t i = 0; i < cycle.size(); ++i) {
            const auto u = source(cycle[i], mesh_);

            for (std::size_t j = 0; j < cycle.size(); ++j) {
              // The first and last vertices are consecutive in a closed cycle.
              std::size_t i_prev = (i + cycle.size() - 1) % cycle.size();
              std::size_t i_next = (i + 1) % cycle.size();

              if (j == i || j == i_prev || j == i_next) {
                continue;
              }

              const auto v_prev = source(cycle[i_prev], mesh_);
              const auto v = source(cycle[j], mesh_);
              const auto v_next = source(cycle[i_next], mesh_);


              const auto [h_uv, found_uv] = halfedge(u, v, mesh_);
              const auto [h_vu, found_vu] = halfedge(v, u, mesh_);

              if (found_uv || found_vu) {
                Vector_3 from_prev = points_pm_[v_prev] - points_pm_[v];
                Vector_3 to_next = points_pm_[v] - points_pm_[v_next];
                Vector_3 to_u = points_pm_[u] - points_pm_[v];
                const Vector_3 normal = normals_pm_[v];

                to_next = to_next / std::sqrt(CGAL::scalar_product(to_next, to_next));
                to_u = to_u / std::sqrt(CGAL::scalar_product(to_u, to_u));
                from_prev = from_prev / std::sqrt(CGAL::scalar_product(from_prev, from_prev));

                Vector_3 right_next = CGAL::cross_product(to_next, normal);
                right_next = right_next / std::sqrt(CGAL::scalar_product(right_next, right_next));

                Vector_3 right_prev = CGAL::cross_product(from_prev, normal);
                right_prev = right_prev / std::sqrt(CGAL::scalar_product(right_prev, right_prev));

                const double d_next = CGAL::scalar_product(to_u, right_next);
                const double d_prev = CGAL::scalar_product(to_u, right_prev);


                if (d_next > 0 || d_prev > 0) {
                  continue;
                }
                
                std::cout << "Found a preexisting halfedge between non-consecutive vertices with indices " << i << " and " << j << std::endl;
                std::cout << "Vertex u: " << u << std::endl;
                std::cout << "Vertex v: " << v << std::endl;
                std::cout << "Halfedge u->v: " << (found_uv ? "exists" : "does not exist") << std::endl;
                std::cout << "Halfedge v->u: " << (found_vu ? "exists" : "does not exist") << std::endl;
                if (found_uv) {
                  std::cout << "Halfedge u->v descriptor: " << h_uv << std::endl;
                  std::cout << "Halfedge u->v source: " << source(h_uv, mesh_) << ", target: " << target(h_uv, mesh_) << std::endl;
                  std::cout << "Halfedge u->v face: " << face(h_uv, mesh_) << std::endl;

                  std::cout << "All halfedges that go into vertex u:\n";
                  for (auto h : halfedges_around_target(u, mesh_)) {
                    std::cout << h << " " << source(h, mesh_) << "->" << target(h, mesh_) << ", next: " << next(h, mesh_) << "\n";
                  }
                  std::cout << std::endl;
                  std::cout << "All halfedges that go into vertex v:\n";
                  for (auto h : halfedges_around_target(v, mesh_)) {
                    std::cout << h << " " << source(h, mesh_) << "->" << target(h, mesh_) << ", next: " << next(h, mesh_) << "\n";
                  }
                  std::cout << std::endl;
                }
                if (found_vu) {
                  std::cout << "Halfedge v->u descriptor: " << h_vu << std::endl;
                  std::cout << "Halfedge v->u source: " << source(h_vu, mesh_) << ", target: " << target(h_vu, mesh_) << std::endl;
                  std::cout << "Halfedge v->u face: " << face(h_vu, mesh_) << std::endl;

                  std::cout << "All halfedges that go into vertex v:\n";
                  for (auto h : halfedges_around_target(v, mesh_)) {
                    std::cout << h << " " << source(h, mesh_) << "->" << target(h, mesh_) << ", next: " << next(h, mesh_) << "\n";
                  }
                  std::cout << std::endl;
                  std::cout << "All halfedges that go into vertex u:\n";
                  for (auto h : halfedges_around_target(u, mesh_)) {
                    std::cout << h << " " << source(h, mesh_) << "->" << target(h, mesh_) << ", next: " << next(h, mesh_) << "\n";
                  }
                  std::cout << std::endl;
                }
                std::cout << "In the cycle: ";
                for (std::size_t k = 0; k < cycle.size(); ++k) {
                  std::cout << cycle[k] << " [" << source(cycle[k], mesh_) << "]"<< (k + 1 < cycle.size() ? " -> " : "");
                }
                std::cout << std::endl;

                export_halfedge_mesh_json(
                  mesh_,
                  "mesh_debug.json"
                );
              }

              CGAL_assertion(!found_uv && !found_vu);
            }
          }

          // Ear clipping converts the polygon cycle into triangle faces while preserving
          // the existing halfedge connectivity.
          // TODO: should this have any reason to fail?
          bool success = ear_clip_and_add_faces(cycle);

          // Check that each halfedge in the cycle is now assigned to a non-null face
          if (success) {
            for (halfedge_descriptor h : cycle) {
              if (face(h, mesh_) == null_f) {
                std::cerr << "Halfedge " << h << " has null face!" << std::endl;
                export_halfedge_mesh_json(
                  mesh_,
                  "mesh_debug_null_face.json"
                );
              }
              CGAL_assertion(face(h, mesh_) != null_f);
            }
          }
        }
      }

      double triangle_angle(const Point_2& a, const Point_2& b, const Point_2& c) const {
        const double ux = CGAL::to_double(a.x() - b.x());
        const double uy = CGAL::to_double(a.y() - b.y());
        const double vx = CGAL::to_double(c.x() - b.x());
        const double vy = CGAL::to_double(c.y() - b.y());

        const double dot = ux * vx + uy * vy;
        const double nu = std::sqrt(ux * ux + uy * uy);
        const double nv = std::sqrt(vx * vx + vy * vy);
        if (nu <= 0.0 || nv <= 0.0)
          return 0.0;

        double cs = dot / (nu * nv);
        // cs = (std::max)(-1.0, (std::min)(1.0, cs));
        return std::acos(cs);
      }

      CGAL::Orientation orient3(const Point_3& a, const Point_3& b, const Point_3& c, const Vector_3& n) const {
        const Vector_3 ab = b - a;
        const Vector_3 ac = c - a;
        const Vector_3 cross = CGAL::cross_product(ab, ac);
        return (cross * n > 0) ? CGAL::LEFT_TURN : CGAL::RIGHT_TURN;
      }

      std::vector<std::pair<double, std::size_t>> find_ear_candidates(
          const std::vector<halfedge_descriptor>& cycle,
          const std::vector<vertex_descriptor>& verts,
          const std::vector<Point_2>& P,
          const Vector_3& n,
          const Point_3& origin,
          const Vector_3& u,
          const Vector_3& v,
          const std::vector<int>& idx,
          bool print_info)
      {
        // Verify inputs
        CGAL_assertion(cycle.size() >= 3);
        CGAL_assertion(verts.size() == cycle.size());
        CGAL_assertion(idx.size() == cycle.size());
        CGAL_assertion(P.size() == cycle.size());

        // The cycle must be well-formed (properly connected).
        for (std::size_t i = 0; i < cycle.size(); ++i) {
          const halfedge_descriptor h = cycle[i];
          const halfedge_descriptor h_next = cycle[(i+1) % cycle.size()];

          CGAL_assertion(next(h, mesh_) == h_next);
          CGAL_assertion(target(h, mesh_) == source(h_next, mesh_));
          CGAL_assertion(source(h, mesh_) == verts[i]);
          CGAL_assertion(target(h, mesh_) == verts[(i+1) % cycle.size()]);
        }

        // The vectors must be non-zero and orthogonal to each other.
        CGAL_assertion(u != CGAL::NULL_VECTOR);
        CGAL_assertion(v != CGAL::NULL_VECTOR);
        CGAL_assertion(n != CGAL::NULL_VECTOR);

        CGAL_assertion(std::abs(u * v) < 1e-8);
        CGAL_assertion(std::abs(u * n) < 1e-8);
        CGAL_assertion(std::abs(v * n) < 1e-8);

        // And must be normalized
        CGAL_assertion(std::abs(u * u - 1.0) < 1e-8);
        CGAL_assertion(std::abs(v * v - 1.0) < 1e-8);
        CGAL_assertion(std::abs(n * n - 1.0) < 1e-8);
        

        std::vector<std::pair<double, std::size_t>> ear_candidates;

        // 2D orientation test.
        const auto orient2 = [](const Point_2& a, const Point_2& b, const Point_2& c) {
          return CGAL::orientation(a, b, c);
        };

        for (std::size_t i = 0; i < idx.size(); ++i)
        {
          // Candidate ear = (prev,current,next).
          const std::size_t i0 = idx[(i + idx.size() - 1) % idx.size()];
          const std::size_t i1 = idx[i];
          const std::size_t i2 = idx[(i + 1) % idx.size()];

          const Point_2& a = P[i0];
          const Point_2& b = P[i1];
          const Point_2& c = P[i2];

          // Check if the triangle is oriented correctly (counter-clockwise).
          // if (orient2(a, b, c) != CGAL::LEFT_TURN) {

          //   if (print_info) {
          //     halfedge_descriptor h = halfedge(verts[i0], verts[i1], mesh_).first;
          //     std::cerr << "Halfedge with null face: " << h << std::endl;
          //     std::cerr << "Triangle not oriented correctly: " << i0 << ", " << i1 << ", " << i2 << std::endl;
          //     std::cerr << "Triangle vertices (2D): a=" << a << ", b=" << b << ", c=" << c << std::endl;
          //     std::cerr << "Orientation: " << orient2(a, b, c) << std::endl;
          //     export_halfedge_mesh_json(
          //       mesh_,
          //       "mesh_debug_null_face.json"
          //     );
          //     CGAL_assertion(false);
          //   }
          //   continue;
          // }

          bool wrong_orientation_3d = CGAL::cross_product(get(points_pm_, verts[i2]) - get(points_pm_, verts[i1]),
                                  get(points_pm_, verts[i0]) - get(points_pm_, verts[i1])) * n < FT(0);

          const Point_3& a3 = get(points_pm_, verts[i0]);
          const Point_3& b3 = get(points_pm_, verts[i1]);
          const Point_3& c3 = get(points_pm_, verts[i2]);

          const Vector_3 b_normal = grid_.box_normal(b3);

          bool wrong_orientation_3d_new = orient3(a3, b3, c3, b_normal) != CGAL::LEFT_TURN;

          if (wrong_orientation_3d && !wrong_orientation_3d_new)
          {
            if (print_info) {
              std::cerr << "Old check fails but new check passes" << std::endl;
            }
          }


          // Positive: left turn, negative: right turn in 3D projection.
          if (wrong_orientation_3d_new)
          {
            // This is fine since it can happen to concave polygons
            // However, if this happens without a preexisting connection between a and c, it is incorrect
            // Check if there is a preexisting connection between a and c
            //auto [he_ac, exists_ac] = halfedge(verts[i0], verts[i2], mesh_);
            // if (!exists_ac) {
            //   std::cout << "Halfedge between a and c exists: " << exists_ac << std::endl;
            //   std::cout << "Wrong orientation in 3D for triangle: " << i0 << ", " << i1 << ", " << i2 << std::endl;
            //   std::cout << "Triangle vertices (3D): a3=" << a3 << ", b3=" << b3 << ", c3=" << c3 << std::endl;
            //   std::cout << "Reference normal: " << b_normal << std::endl;
            //   std::cout << "Halfedge: " << halfedge(verts[i0], verts[i1], mesh_).first << std::endl;
            //   std::cout << "Orientation check: wrong_orientation_3d=" << wrong_orientation_3d << ", wrong_orientation_3d_new=" << wrong_orientation_3d_new << std::endl;
            //   // Wrong orientation in 3D, can happen on concave polygons
            //   export_halfedge_mesh_json(mesh_, "wrong_orientation.json");
            //   CGAL_assertion(false);
            // }
            continue;
          }

          bool intersecting = false;
          for (std::size_t j = 0; j < verts.size(); ++j) {
            for (halfedge_descriptor h : halfedges_around_source(verts[j], mesh_)) {
              vertex_descriptor t = target(h, mesh_);
              if (verts[j] == verts[i0] || verts[j] == verts[i2] ||
                  t == verts[i0] || t == verts[i2])
                continue;
              auto x = project(get(points_pm_, verts[j]), origin, u, v);
              auto y = project(get(points_pm_, t), origin, u, v);

              if (intersect_2d(a, c, x, y, false)) {
                // std::cout << "\nTriangle intersects with edge: " << i0 << ", " << i1 << ", " << i2 << " with edge " << h << std::endl;
                // std::cout << "Edge coordinates: (" << x << ") to (" << y << ")" << std::endl;
                // std::cout << "3d line " << get(points_pm_, verts[j]) << " to " << get(points_pm_, t) << std::endl;
                // std::cout << "Projection info: origin: " << origin << ", u: " << u << ", v: " << v << std::endl;
                // std::cout << "Normal: " << n << std::endl;
                // std::cout << "Projecting p1: " << get(points_pm_, verts[j]) << " to 2D: " << x << std::endl;
                // std::cout << "Projecting p2: " << get(points_pm_, t) << " to 2D: " << y << std::endl;
                // std::cout << "Triangle coordinates: (" << a << ") to (" << b << ") to (" << c << ")" << std::endl;

                // std::cout << "python data:" << std::endl;
                // std::cout << "cycle_points = [";
                // for (std::size_t k = 0; k < verts.size(); ++k) {
                //   Point_2 p = P[k];
                //   std::cout << "[" << CGAL::to_double(p.x()) << ", " << CGAL::to_double(p.y()) << "]";
                //   if (k + 1 < verts.size()) std::cout << ", ";
                // }
                // std::cout << "]" << std::endl;
                // std::cout << "x, y = ([" << x.x() << ", " << x.y() << "], [" << y.x() << ", " << y.y() << "])" << std::endl;
                // std::cout << "a, c = ([" << a.x() << ", " << a.y() << "], [" << c.x() << ", " << c.y() << "])" << std::endl;
                intersecting = true;
                break;
              }
            }
            if (intersecting)
              break;
          }

          if (intersecting) {
            if (print_info) {
              std::cerr << "Triangle intersects with another edge in 2D: " << i0 << ", " << i1 << ", " << i2 << std::endl;
            }
            continue;
          }

          const double ang = compute_triangle_angle(a3, b3, c3, b_normal);
          ear_candidates.emplace_back(ang, i); // Store candidate together with its interior angle.
        }

        return ear_candidates;
      }

      Point_2 project(const Point_3& p, const Point_3& origin, const Vector_3& u, const Vector_3& v) const {
        return Point_2(CGAL::to_double((p - origin) * u),
                       CGAL::to_double((p - origin) * v));
      }
      
      Vector_3 compute_cycle_normal(const std::vector<vertex_descriptor>& verts) const {
        Vector_3 n = CGAL::NULL_VECTOR;
        for (auto v: verts) {
          n = n + grid_.box_normal(get(points_pm_, v));
          //n = n + normals_pm_[v];
        }
        const double nlen2 = CGAL::to_double(n.squared_length());
        CGAL_assertion(nlen2 > 0.0);
        n = n / std::sqrt(nlen2);

        return n;
      }

      double compute_triangle_angle(const Point_3& a, const Point_3& b, const Point_3& c, const Vector_3& reference_normal) const {
        // First project the points onto the plane defined by the reference normal.
        const Vector_3 n = reference_normal;
        const Vector_3 u = (b - a) - ((b - a) * n) * n;
        const Vector_3 v = CGAL::cross_product(n, u);
        const Point_2 pa = project(a, a, u, v);
        const Point_2 pb = project(b, a, u, v);
        const Point_2 pc = project(c, a, u, v);

        // Compute the interior angle at vertex b in the projected 2D triangle.
        const double angle = triangle_angle(pa, pb, pc);
        return angle;
      }

      /**
       * @brief Triangulates a polygonal halfedge cycle using ear clipping.
       *
       * The cycle is projected once onto a local tangent plane. Valid ears are
       * identified using 2D orientation tests, containment tests, and a 3D normal
       * consistency test. Candidate ears are ordered by their interior angle.
       *
       * After an ear is committed, the halfedge cycle is updated and the remaining
       * polygon is processed again. When only three halfedges remain, the final
       * triangle is created directly.
       *
       * @param cycle Halfedge cycle representing a polygon boundary.
       * @return True if the polygon was successfully triangulated, false otherwise.
      */
      bool ear_clip_and_add_faces(std::vector<halfedge_descriptor> cycle) {
        CGAL_assertion(cycle.size() >= 3);

        // 2D orientation test.
        const auto orient2 = [](const Point_2& a, const Point_2& b, const Point_2& c) {
          return CGAL::orientation(a, b, c);
        };

        // Returns true if p lies on or inside triangle (a,b,c).
        const auto point_in_triangle =
          [&](const Point_2& p,
              const Point_2& a,
              const Point_2& b,
              const Point_2& c) -> bool {
          const auto o1 = orient2(a, b, p);
          const auto o2 = orient2(b, c, p);
          const auto o3 = orient2(c, a, p);
          return (o1 == CGAL::LEFT_TURN &&
                  o2 == CGAL::LEFT_TURN &&
                  o3 == CGAL::LEFT_TURN);
        };

        std::vector<vertex_descriptor> verts;
        verts.reserve(cycle.size());
        for (halfedge_descriptor h : cycle)
          verts.push_back(source(h, mesh_));

        // --------------------------------------------------------------------------
        // Extract the current polygon vertices.
        //
        // A polygon cycle must contain each vertex only once. Repeated vertices mean
        // that the halfedge walk is actually a pinched/multi-cycle boundary and cannot
        // be safely triangulated as a simple polygon.
        // --------------------------------------------------------------------------
        auto tmp = verts;
        std::sort(tmp.begin(), tmp.end());
        if (std::unique(tmp.begin(), tmp.end()) != tmp.end()) {
          std::cerr << "Warning: Ear clipping failed: cycle contains repeated vertices. Can result in holes.\n";
          for (halfedge_descriptor h : cycle) {
            vertex_descriptor v = source(h, mesh_);
            std::cerr << "  Vertex " << v << " at " << get(points_pm_, v) << "\n";
          }
          CGAL_assertion(false);
        }

        // Compute one supporting normal for the polygon. This normal is fixed for the
        // entire ear-clipping operation so that the projection does not change after
        // each ear removal.
        Vector_3 n = compute_cycle_normal(verts);

        // Project the polygon vertices into the tangent plane once. The resulting 2D
        // points are reused by the ear-clipping loop.
        std::vector<Vector_3> frame = grid_.compute_local_tangent_frame(n);
        Vector_3 u = frame[0];
        Vector_3 v = frame[1];
        CGAL_assertion(u != CGAL::NULL_VECTOR && v != CGAL::NULL_VECTOR);
        if (CGAL::cross_product(u, v) * n < FT(0))
          v = -v;
        const Point_3 origin = get(points_pm_, verts[0]);

        std::vector<Point_2> P;
        P.reserve(verts.size());

        // Iteratively remove one ear until only one triangle remains.
        while (cycle.size() > 3) {
          verts.clear();
          P.clear();

          verts.reserve(cycle.size());
          for (halfedge_descriptor h : cycle)
            verts.push_back(source(h, mesh_));

          for (vertex_descriptor vd : verts)
          {
            const Point_3 p = get(points_pm_, vd);
            P.emplace_back(project(p, origin, u, v));
          }

          std::vector<int> idx(P.size());
          std::iota(idx.begin(), idx.end(), 0);

          // ----------------------------------------------------------------------
          // Find every valid ear.
          // ----------------------------------------------------------------------

          std::vector<std::pair<double, std::size_t>> ear_candidates
            = find_ear_candidates(
                cycle,
                verts,
                P,
                n,
                origin,
                u,
                v,
                idx,
                false
              );
          
          // There should always be at least one valid ear unless the polygon has 3 or fewer vertices.
          bool found_ear = !ear_candidates.empty();

          if (!found_ear && cycle.size() > 3) {
            std::cerr << "No valid ear found for a polygon with more than 3 vertices.\n";
            std::vector<Point_3> points_in_cycle;
            for (halfedge_descriptor h : cycle) {
              points_in_cycle.push_back(mesh_.point(source(h, mesh_)));
              std::cout << "source in cycle: " << source(h, mesh_) << std::endl;
            }
            std::cout << "3d points: [";
            for (std::size_t i = 0; i < points_in_cycle.size(); ++i) {
              std::cout << "[" << points_in_cycle[i].x() << ", " << points_in_cycle[i].y() << ", " << points_in_cycle[i].z() << "]";
              if (i + 1 < points_in_cycle.size()) {
                std::cout << ", ";
              }
            }
            std::cout << "]" << std::endl;
            
            std::cout << "2d points: (projected on normal << " << n << ") [";
            for (std::size_t i = 0; i < P.size(); ++i) {
              std::cout << "[" << P[i].x() << ", " << P[i].y() << "]";
              if (i + 1 < P.size()) {
                std::cout << ", ";
              }
            }
            std::cout << "]" << std::endl;

            std::cout << "Find ear candidates output:" << std::endl;
            find_ear_candidates(
                cycle,
                verts,
                P,
                n,
                origin,
                u,
                v,
                idx,
                true
              );


            // int target_vertex = 349;
            // for (halfedge_descriptor h : cycle) {
            //   points_in_cycle.push_back(mesh_.point(source(h, mesh_)));
            //   std::cout << "source in cycle: " << source(h, mesh_) << std::endl;
            //   int target_vertex = 349;
            //   if (source(h, mesh_).idx() == target_vertex) {
            //     std::cout << "h: " << h << std::endl;
            //     export_halfedge_mesh_json(mesh_, "no_ear_vertex_349.json");
            //     std::cerr << "Found source vertex " << target_vertex << " in cycle." << std::endl;
            //     CGAL_assertion(false);
            //   }
            // }


            export_halfedge_mesh_json(mesh_, "no_ear.json");

            //CGAL_assertion(false);
            return false;
          }

          // There should always be at least one valid ear unless the polygon has 3 or fewer vertices.
          CGAL_assertion(found_ear || cycle.size() <= 3);

          // ----------------------------------------------------
          // commit_ear_step must:
          //  - add the ear triangle face,
          //  - update the mesh halfedge wiring,
          //  - erase the ear vertex from 'cycle'.
          // ----------------------------------------------------

          // Try ears from smallest angle to largest until one can be inserted topologically.
          std::sort(ear_candidates.begin(), ear_candidates.end());
          bool ear_committed = false;
          for (const auto& [angle, best_pos] : ear_candidates) {
            // commit_ear_step will now always return true, so no need to check its return value.
            CGAL_assertion(commit_ear_step(cycle, best_pos));
            break;
          }
          ear_committed = true;

          // This code really shouldn't be reached
          if (!ear_committed) {
            if (cycle.size() == 4) {
              for (std::size_t i = 0; i < cycle.size(); ++i) {
                commit_ear_step(cycle, i);
              }
            }
            return true;
          }
        }

        CGAL_assertion(cycle.size() == 3);

        // Commit the last remaining triangle.
        if (cycle.size() == 3)
        {
          CGAL_assertion(commit_ear_step(cycle, 1));
        }

        CGAL_assertion(cycle.empty());

        return true;
      }

      /**
       * @brief Commits one ear removal to the halfedge mesh.
       *
       * For an ear
       *
       *   a -> b -> c
       *
       * a new diagonal `a -> c` is created, the triangle `a-b-c` is assigned a face,
       * and the remaining polygon boundary is rewired to replace
       *
       *   a -> b -> c
       *
       * by
       *
       *   a -> c.
       *
       * The supplied cycle is updated to remove the ear vertex.
       *
       * @param cycle Current polygon halfedge cycle.
       * @param ear_pos Position of the ear vertex in `cycle`.
       *
       * @return `true` if the ear was successfully committed.
      */
      bool commit_ear_step(std::vector<halfedge_descriptor>& cycle,
                           std::size_t ear_pos) {
        const auto null_h = boost::graph_traits<TriangleMesh>::null_halfedge();
        const auto null_f = boost::graph_traits<TriangleMesh>::null_face();
        const auto null_v = boost::graph_traits<TriangleMesh>::null_vertex();

        const std::size_t n = cycle.size();

        // Ensure that the cycle has at least three vertices and the ear position is valid.
        CGAL_assertion(n >= 3);
        CGAL_assertion(ear_pos < n);

        // The final cycle already consists of exactly three halfedges, so no new
        // diagonal is needed. Convert the cycle directly into a triangular face.
        if (n == 3) {

          // Here we must have a valid cycle of three triangles
          CGAL_assertion(cycle[0] != null_h && cycle[1] != null_h && cycle[2] != null_h);
          CGAL_assertion(target(cycle[0], mesh_) == source(cycle[1], mesh_));
          CGAL_assertion(target(cycle[1], mesh_) == source(cycle[2], mesh_));
          CGAL_assertion(target(cycle[2], mesh_) == source(cycle[0], mesh_));
          CGAL_assertion(face(cycle[0], mesh_) == null_f);
          CGAL_assertion(face(cycle[1], mesh_) == null_f);
          CGAL_assertion(face(cycle[2], mesh_) == null_f);
          CGAL_assertion(next(cycle[0], mesh_) == cycle[1]);
          CGAL_assertion(next(cycle[1], mesh_) == cycle[2]);
          CGAL_assertion(next(cycle[2], mesh_) == cycle[0]);

          face_descriptor f = add_face(mesh_);
          CGAL_assertion(f != null_f);

          set_halfedge(f, cycle[0], mesh_);
          for (halfedge_descriptor h : cycle) {
            set_face(h, f, mesh_);
          }

          // This shouldn't be necessary with the assertions above
          // set_next(cycle[0], cycle[1], mesh_);
          // set_next(cycle[1], cycle[2], mesh_);
          // set_next(cycle[2], cycle[0], mesh_);

          cycle.clear();
          return true;
        }

        const std::size_t i_ab     = (ear_pos + n - 1) % n;
        const std::size_t i_bc     = ear_pos;
        const std::size_t i_before = (ear_pos + n - 2) % n;
        const std::size_t i_after  = (ear_pos + 1) % n;

        // Identify the two halfedges entering and leaving the ear vertex.
        halfedge_descriptor h_ab = cycle[i_ab];
        halfedge_descriptor h_bc = cycle[i_bc];

        CGAL_assertion(h_ab != null_h && h_bc != null_h);
        if (h_ab == null_h || h_bc == null_h) {
          return false;
        }

        const vertex_descriptor a = source(h_ab, mesh_);
        const vertex_descriptor b = target(h_ab, mesh_);
        const vertex_descriptor c = target(h_bc, mesh_);

        CGAL_assertion(a != null_v && b != null_v && c != null_v);
        CGAL_assertion(source(h_bc, mesh_) == b);

        if (source(h_bc, mesh_) != b) {
          return false;
        }

        // Vector_3 normal = normals_pm_[b];

        // std::vector<Vector_3> frame = grid_.compute_local_tangent_frame(normal);
        // const Vector_3& u = frame[0];
        // const Vector_3& v = frame[1];

        // auto project_to_2d = [&](const Point_3& p) -> Point_2 {
        //   const Vector_3 d = p - points_pm_[b];
        //   return Point_2(d * u, d * v);
        // };

        // const bool found = halfedge(a, c, mesh_).second;

        // if (found) {

        //   // Check if the midpoint of the diagonal (a, c) lies inside the projected polygon.
        //   // If it does, this is a real issue because it indicates that the diagonal intersects the interior of the polygon.


        //   // Project all vertices of the loop onto the plane defined by the normal at vertex b.
        //   CGAL::Polygon_2<Kernel> projected_polygon;

        //   for (halfedge_descriptor h : cycle) {
        //     const Point_3 p = points_pm_[target(h, mesh_)];
        //     projected_polygon.push_back(project_to_2d(p));
        //   }

        //   if (projected_polygon.is_simple()) {
        //     Point_3 pa = points_pm_[a];
        //     Point_3 pc = points_pm_[c];

        //     Point_3 mid_ac = CGAL::midpoint(pa, pc);

        //     Point_2 sample_point = project_to_2d(points_pm_[b]);
        //     const auto side = projected_polygon.bounded_side(sample_point);

        //     // If it is inside, that indicates a wrong connection
        //     CGAL_assertion(side != CGAL::ON_BOUNDED_SIDE);

        //     // Definitely no need to add another edge
        //     return false;
        //   }

        //   // TODO: add non simple check?
        // }
        

        // Create the diagonal a -> c. Its opposite halfedge c -> a is used to close
        // the new triangular face.
        halfedge_descriptor h_ac = create_halfedge(a, c);
        halfedge_descriptor h_ca = opposite(h_ac, mesh_);

        // The two boundary halfedges of the ear must not already belong to faces.
        CGAL_assertion(face(h_ac, mesh_) == null_f && face(h_ca, mesh_) == null_f);
        CGAL_assertion(face(h_ab, mesh_) == null_f && face(h_bc, mesh_) == null_f);
        if (face(h_ac, mesh_) != null_f && face(h_ca, mesh_) != null_f) {
          std::cout << "both new halfedges already belong to a face" << std::endl;
          return false;
        }
        if (face(h_ab, mesh_) != null_f || face(h_bc, mesh_) != null_f) {
          std::cout << "one of the boundary halfedges already belongs to a face" << std::endl;
          return false;
        }
        
        // Construct the new triangular face:
        //     a -> b -> c -> a
        face_descriptor f = add_face(mesh_);
        CGAL_assertion(f != null_f);
        if (f == null_f) {
          return false;
        }

        if (face(h_ca, mesh_) == null_f) {
          set_halfedge(f, h_ab, mesh_);

          set_face(h_ab, f, mesh_);
          set_face(h_bc, f, mesh_);
          set_face(h_ca, f, mesh_);

          set_next(h_ab, h_bc, mesh_);
          set_next(h_bc, h_ca, mesh_);
          set_next(h_ca, h_ab, mesh_);
        }

        // Rewire the remaining boundary:
        // ... -> a -> b -> c -> ...
        // becomes
        // ... -> a -> c -> ...
        halfedge_descriptor h_before = cycle[i_before];
        halfedge_descriptor h_after  = cycle[i_after];
        CGAL_assertion(h_before != null_h && h_after != null_h);

        set_next(h_before, h_ac, mesh_);
        set_next(h_ac, h_after, mesh_);

        // Update the explicit cycle representation so that the removed ear vertex
        // is no longer part of the polygon being triangulated.
        cycle[i_ab] = h_ac;
        cycle.erase(cycle.begin() + static_cast<std::ptrdiff_t>(i_bc));

        // Verify that the cycle is still valid after removing the ear vertex.
        for (std::size_t i = 0; i < cycle.size(); ++i) {
          CGAL_assertion(next(cycle[i], mesh_) == cycle[(i + 1) % cycle.size()]);
          CGAL_assertion(target(cycle[i], mesh_) == source(cycle[(i + 1) % cycle.size()], mesh_));
          CGAL_assertion(face(cycle[i], mesh_) == null_f);
        }

        // Verify that the newly created face has the correct halfedges.
        CGAL_assertion(face(h_ab, mesh_) == f);
        CGAL_assertion(face(h_bc, mesh_) == f);
        CGAL_assertion(face(h_ca, mesh_) == f);
        CGAL_assertion(next(h_ab, mesh_) == h_bc);
        CGAL_assertion(next(h_bc, mesh_) == h_ca);
        CGAL_assertion(next(h_ca, mesh_) == h_ab);
        CGAL_assertion(source(h_ab, mesh_) == a);
        CGAL_assertion(source(h_bc, mesh_) == b);
        CGAL_assertion(source(h_ca, mesh_) == c);
        CGAL_assertion(target(h_ab, mesh_) == b);
        CGAL_assertion(target(h_bc, mesh_) == c);
        CGAL_assertion(target(h_ca, mesh_) == a);

        return true;
      }

      /**
       * @brief Tests whether two 2D segments intersect in their interiors.
       *
       * Endpoint touching is not considered a strict intersection. When
       * `allow_overlap` is false, overlapping segments are considered intersecting.
       *
       * @param a First endpoint of the first segment.
       * @param b Second endpoint of the first segment.
       * @param c First endpoint of the second segment.
       * @param d Second endpoint of the second segment.
       * @param allow_overlap Whether overlapping segments should be ignored.
       *
       * @return `true` if the segments violate the local intersection constraint.
      */
      bool intersect_2d(const Point_2& a,
                        const Point_2& b,
                        const Point_2& c,
                        const Point_2& d,
                        bool allow_overlap = false) const {
        const CGAL::Segment_2<Kernel> s1(a, b);
        const CGAL::Segment_2<Kernel> s2(c, d);

        auto inter = CGAL::intersection(s1, s2);
        if (!inter) {
          return false;
        }

        Point_2 ip;
        if (CGAL::assign(ip, *inter)) {
          // Allow touching at endpoints.
          return !(ip == a || ip == b || ip == c || ip == d);
        }

        // If the intersection is not a point, it is an overlapping segment.
        if (!allow_overlap) {
          return true;
        }
        return false;
      };

      /**
       * @brief Checks whether a candidate is geometrically compatible with the
       *        surrounding mesh under a local planar projection.
       *
       * The parent-to-candidate edges are projected into the tangent plane defined
       * by the supplied normal. They are rejected when they intersect existing mesh
       * edges.
       *
       * @param cand Candidate reconstruction vertex.
       * @param n Projection normal.
       *
       * @return `true` when the candidate passes the local intersection test.
      */
      bool projection_check(const Candidate& cand, const Vector_3& n) const {
        CGAL_assertion(n != CGAL::NULL_VECTOR);

        const Point_3 pa = get(points_pm_, cand.first);
        const Point_3 pb = get(points_pm_, cand.second);
        const Point_3 pc = cand.position;

        auto local_frame = grid_.compute_local_tangent_frame(n);
        Vector_3 u = local_frame[0];
        Vector_3 v = local_frame[1];

        const Point_2 A = project(pa, pc, u, v);
        const Point_2 B = project(pb, pc, u, v);
        const Point_2 C = project(pc, pc, u, v);

        // New projected edges created if this candidate is accepted.
        const Point_2 new_e1_s = A;
        const Point_2 new_e1_t = C;
        const Point_2 new_e2_s = B;
        const Point_2 new_e2_t = C;

        //CHANGED: made the vicinity d instead of 4d
        auto nearby = nearby_mesh_vertices(pc, FT(1.0) * grid_.get_box_size());
        for (auto v : nearby) {
          const FT distance = std::sqrt(CGAL::squared_distance(pc, get(points_pm_, v)));
          // if (distance > FT(1.0) * grid_.get_box_size()) {
          //   continue;
          //   std::cout << "Vertex " << v << " is too far from candidate: distance = " << distance << " while allowed distance = " << FT(1.0) * grid_.get_box_size() << std::endl;
          // }
          //CGAL_assertion(distance <= FT(1.0) * grid_.get_box_size());

          for (halfedge_descriptor h : halfedges_around_target(v, mesh_)) {
            vertex_descriptor s = source(h, mesh_);
            vertex_descriptor t = target(h, mesh_);

            // Vector_3 n_s = get(normals_pm_, s);
            // Vector_3 n_t = get(normals_pm_, t);
            // Vector_3 n_h = n_s + n_t;
            // n_h = n_h / std::sqrt(n_h.squared_length());
            // if (n_h * n < 0) {
            //   continue;
            // }

            const Point_3 q0 = get(points_pm_, s);
            const Point_3 q1 = get(points_pm_, t);

            const Point_2 E0 = project(q0, pc, local_frame[0], local_frame[1]);
            const Point_2 E1 = project(q1, pc, local_frame[0], local_frame[1]);

            if (intersect_2d(new_e1_s, new_e1_t, E0, E1)) {
              return false;
            }
            if (intersect_2d(new_e2_s, new_e2_t, E0, E1)) {
              return false;
            }
          }
        }

        return true;
      }

      /**
       * @brief Searches the local halfedge neighborhood for another parent vertex.
       *
       * The search follows `next()` pointers for a bounded number of steps from
       * every incident halfedge.
       *
       * @param start_v Starting vertex.
       * @param target_v Vertex being searched for.
       *
       * @return `true` if the target vertex is encountered.
      */
      bool walk_finds_other_parent_from_vertex(vertex_descriptor start_v,
                                               vertex_descriptor target_v) const {

        for (halfedge_descriptor h0 : halfedges_around_source(start_v, mesh_)) {
          if (h0 == boost::graph_traits<TriangleMesh>::null_halfedge()) {
            return false;
          }

          halfedge_descriptor h = h0;
          for (std::size_t i = 0; i < window_size_; ++i) {
            const vertex_descriptor s = source(h, mesh_);
            const vertex_descriptor t = target(h, mesh_);

            if (s == target_v || t == target_v) {
              return true;
            }

            const halfedge_descriptor n = next(h, mesh_);
            if (n == boost::graph_traits<TriangleMesh>::null_halfedge()) {
              return false;
            }

            h = n;
          }
        }

        return false;
      }

      /**
       * @brief Determines whether a candidate connects two distinct boundary fronts.
       *
       * A candidate is considered to join two borders when neither parent can reach
       * the other through the local halfedge walks.
      */
      bool joins_two_borders(const Candidate& cand) const {
        const bool a_finds_b = walk_finds_other_parent_from_vertex(cand.first, cand.second);
        const bool b_finds_a = walk_finds_other_parent_from_vertex(cand.second, cand.first);

        return !a_finds_b && !b_finds_a;
      }

      /**
       * @brief Computes the processing priority of a reconstruction candidate.
       *
       * Higher priority is assigned to candidates adjacent to low-degree vertices,
       * followed by candidates that connect two separate boundary fronts.
       *
       * @param cand Candidate being prioritized.
       *
       * @return Integer priority in the range `[0, NUM_PRIORITIES-1]`.
      */
      int compute_priority(const Candidate& cand) const {

        // CHANGED: added this new priority rule
        // "At least one of the parent vertices has no edges incident to it"
        if (out_degree(cand.first, mesh_) == 0 || out_degree(cand.second, mesh_) == 0) {
          return 3; // highest priority
        }

        // "At least one of the parent vertices is a vertex with only one edge incident to it"
        if (out_degree(cand.first, mesh_) == 1 || out_degree(cand.second, mesh_) == 1) {
          return 2; // second highest priority
        }

        // "Inserting vc and its two edges joins two borders"
        if (joins_two_borders(cand)) {
          return 1; // medium priority
        }

        // "Inserting vc splits a border"
        return 0;   // lowest priority
      }

    private:

      // Constructor to set the state for a test
      Splat_surface_reconstruction_3(const Grid& grid, TriangleMesh& mesh, std::vector<Vector_3> normals, int i)
        : grid_(grid),
          mesh_(mesh),
          points_pm_(get(CGAL::vertex_point, mesh))
      {
        std::tie(normals_pm_, std::ignore) =
          mesh_.template add_property_map<vertex_descriptor, Vector_3>(
            "v:normal", CGAL::NULL_VECTOR);

        for (vertex_descriptor v : vertices(mesh_))
        {
          std::cout << "Setting normal for vertex " << v.idx() << " to " << normals[v.idx()] << std::endl;
          normals_pm_[v] = normals[v.idx()];
        }

        std::tie(edge_weight_pm_, std::ignore) =
          mesh_.template add_property_map<edge_descriptor, FT>(
            "e:weight", FT(1));

        mesh_vertices_per_cell_.resize(grid_.number_of_cells());
      }

      const Grid& grid_;
      TriangleMesh& mesh_;

      typename TriangleMesh:: template Property_map<vertex_descriptor,Point_3> points_pm_;
      typename TriangleMesh:: template Property_map<vertex_descriptor,Vector_3> normals_pm_;

      typename TriangleMesh::template Property_map<edge_descriptor, FT> edge_weight_pm_;

      vertex_descriptor v0_{boost::graph_traits<TriangleMesh>::null_vertex()};
      vertex_descriptor v1_{boost::graph_traits<TriangleMesh>::null_vertex()};
      bool seeded_ = false;

      static constexpr int NUM_PRIORITIES = 4; // 0 is lowest, NUM_PRIORITIES-1 is highest
      std::array<std::deque<Candidate>, NUM_PRIORITIES> candidate_queues_; // candidate_queues_[p] contains all candidates with priority p

      const std::size_t window_size_ = 8;

      std::vector<std::vector<vertex_descriptor>> mesh_vertices_per_cell_;
  };

} // namespace CGAL

#endif // CGAL_SPLAT_SURFACE_RECONSTRUCTION_H
