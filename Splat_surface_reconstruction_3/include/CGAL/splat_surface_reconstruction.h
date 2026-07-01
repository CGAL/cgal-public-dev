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
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <optional>
#include <queue>

namespace CGAL {

  /*!
    \ingroup PkgSplatSurfaceReconstruction3Ref

    Performs splat surface reconstruction as follows:

    - compute the ...
    - outputs the result in a polygon mesh

    This function relies mainly on the size parameter `spacing`.

    \tparam PointRange is a model of `Range`.
    \tparam NormalRange is a model of `Range`.
    \tparam PolygonMesh a model of `MutableFaceGraph` with an internal
    point property map.


    \param points the range of points.
    \param normals the range of normals.
    \param output_mesh where the reconstruction is stored.
    \param spacing size parameter.
    \return `true` if reconstruction succeeded, `false` otherwise.
  */

  /**
   * @brief Runs splat surface reconstruction on a point cloud with normals.
   *
   * The function takes a range of input points, a matching range of normals,
   * and reconstructs a polygon mesh using the splat-surface reconstruction
   * pipeline.
   *
   * @tparam PointRange  Model of Range whose value type is a CGAL point type.
   * @tparam NormalRange Model of Range whose value type is a CGAL vector type.
   * @tparam PolygonMesh Model of MutableFaceGraph with an internal point map.
   *
   * @param points       Input point range.
   * @param normals      Input normal range.
   * @param output_mesh  Output polygon mesh storing the reconstruction.
   * @param spacing      Global spacing parameter used by the algorithm.
   *
   * @return `true` if reconstruction succeeds, otherwise `false`.
   */
  template <typename PointRange,typename NormalRange, typename PolygonMesh>
  bool Splat_surface_reconstruction(const PointRange& points,
                               const NormalRange& normals,
                               PolygonMesh& output_mesh,
                               double spacing) {
    typedef typename PointRange::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel Kernel;
    typedef typename Kernel::Sphere_3 Sphere;
    typedef typename Kernel::FT FT;

    return true;
  }

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
      std::vector<Index> point_ids;
      Vector_3 normal_sum = CGAL::NULL_VECTOR;
      std::size_t normal_count = 0;

      /**
       * @brief Adds one point and its normal contribution to the cell.
       *
       * @param pid Point index in the input cloud.
       * @param n   Normal associated with the point.
       */
      void add_point(Index pid, const Vector_3& n)
      {
        point_ids.push_back(pid);
        if (n != CGAL::NULL_VECTOR) {
          normal_sum = normal_sum + n;
          ++normal_count;
        }
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
      if (!valid_coords(ix, iy, iz)) {
        return nullptr;
      }
      return &cells_[flat_index(ix, iy, iz)];
    }

    /**
     * @brief Constructs a box grid with explicit bounds.
     *
     * The same cell size is used along all axes, so each cell is cubic.
     *
     * @param box_size Grid cell side length.
     * @param min_x    Minimum x bound.
     * @param max_x    Maximum x bound.
     * @param min_y    Minimum y bound.
     * @param max_y    Maximum y bound.
     * @param min_z    Minimum z bound.
     * @param max_z    Maximum z bound.
     */
    Box_grid(FT box_size,
             FT min_x, FT max_x,
             FT min_y, FT max_y,
             FT min_z, FT max_z)
     : box_size_(box_size),
       min_x_(min_x), max_x_(max_x), min_y_(min_y),
       max_y_(max_y), min_z_(min_z), max_z_(max_z)
    {
      initialize_grid();
    }

    FT get_box_size() const {
      return box_size_;
    }

    FT get_max_splat_radius() const {
      return max_splat_radius_;
    }

    const std::vector<Point_3>& points() const {
      return points_;
    }

    const std::vector<Vector_3>& normals() const {
      return normals_;
    }

    const std::vector<FT>& splat_sizes() const {
      return splat_sizes_;
    }

    /**
     * @brief Builds the grid from input points and normals.
     *
     * Points are inserted into the cell containing their coordinates and each
     * cell accumulates the normals of the points it contains.
     *
     * @param points  Input 3D point cloud.
     * @param normals Input normals, one per point.
     */
    void build(const std::vector<Point_3>& points,
              const std::vector<Vector_3>& normals) {

      clear();

      points_ = points;
      normals_ = normals;

      if (normals_.size() != points_.size()) {
        std::cerr << "Warning: normals size does not match points size.\n";
        std::exit(EXIT_FAILURE);
      }

      for (Index i = 0; i < points_.size(); ++i) {
        insert(i, points_[i], normals_[i]);
      }
    }

    /**
     * @brief Computes one averaged unit normal per grid cell.
     *
     * Empty cells receive `CGAL::NULL_VECTOR`.
     *
     * @return A vector of block normals in grid order.
     */
    std::vector<Vector_3> compute_block_normals() const {
      std::vector<Vector_3> block_normals;
      block_normals.reserve(cells_.size());

      for (std::size_t ix = 0; ix < nx_; ++ix) {
        for (std::size_t iy = 0; iy < ny_; ++iy) {
          for (std::size_t iz = 0; iz < nz_; ++iz) {
            const Cell& c = cells_[flat_index(ix, iy, iz)];

            if (c.normal_count == 0) {
              // zero vector indicates empty cell, no normal information
              block_normals.push_back(Vector_3(0,0,0));
              continue;
            }

            const FT inv = FT(1) / FT(c.normal_count);
            Vector_3 n = c.normal_sum * inv;

            const double len2 = CGAL::to_double(n.squared_length());
            if (len2 <= 0.0) {
              block_normals.push_back(CGAL::NULL_VECTOR);
              continue;
            }

            Vector_3 block_normal = n * (1.0 / std::sqrt(len2));

            // Check consistency with every point normal in this cell.
            for (Index pid : c.point_ids) {
              if (pid >= normals_.size()) {
                continue;
              }

              const Vector_3& point_normal = normals_[pid];
              if (point_normal == CGAL::NULL_VECTOR) {
                continue;
              }

              if (CGAL::to_double(block_normal * point_normal) < 1e-12) {
                std::cout<<CGAL::to_double(block_normal * point_normal) << " at cell (" << ix << "," << iy << "," << iz << ") with point index " << pid << std::endl;
                std::cerr
                  << "Warning: block normal has negative dot product with a point normal.\n"
                  << "box_size is likely too large. Aborting.\n";
                std::exit(EXIT_FAILURE);
              }
            }

            block_normals.push_back(block_normal);
          }
        }
      }

      return block_normals;
    }

    /**
     * @brief Estimates one splat size per point using local 2D triangulation.
     *
     * For each point, the method:
     * - finds the point's cell,
     * - projects neighbors to the tangent plane,
     * - builds a local Delaunay triangulation in 2D,
     * - measures the farthest circumcenter distance.
     *
     * @param global_spacing Default splat size used as fallback.
     *
     * @return A vector of per-point splat radii.
     */
    std::vector<FT> estimate_individual_splat_sizes() {
    FT global_spacing = 1.05*box_size_;
    splat_sizes_.resize(points_.size(), global_spacing);

    for (Index i = 0; i < points_.size(); ++i) {
      int cx, cy, cz;
      if (!to_grid_coords(points_[i], cx, cy, cz)) {
        std::cerr << "Warning: point " << i << " is out of grid bounds, skipping splat size estimation.\n";
        std::exit(EXIT_FAILURE);
      }

      // Collect all points in the same box as p_i.
      std::vector<Index> local_ids;
      const Cell* c = cell(cx, cy, cz);
      if (c == nullptr || c->point_ids.empty()) {
        continue;
      }
      local_ids = c->point_ids;

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
     * @brief Returns a seed pair on the first splat whose radius exceeds box size.
     *
     * @return An optional seed with two seed points and the splat index.
     */
    std::optional<Initial_seed> get_initial_seed() const {
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
      std::exit(EXIT_FAILURE);
      return std::nullopt;
    }

    /**
     * @brief Returns point ids in nearby cells around a query point.
     *
     * @param p Query point.
     * @param radius Search radius.
     *
     * @return Unique point indices gathered from neighboring cells.
     */
    std::vector<Index> nearby_point_ids(const Point_3& p, FT radius) const {
      std::vector<Index> ids;

      int cx, cy, cz;
      if (!to_grid_coords(p, cx, cy, cz)) {
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

            ids.insert(ids.end(), c->point_ids.begin(), c->point_ids.end());
          }
        }
      }

      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
      return ids;
    }

    /**
     * @brief Intersects a circle with a splat.
     *
     * @param parent_a First parent point.
     * @param parent_b Second parent point.
     * @param splat_center Splat center.
     * @param splat_normal Splat normal.
     * @param splat_radius Splat radius.
     *
     * @return Vector of intersection points.
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
       * @brief Computes an orthonormal tangent frame at a point.
       *
       * @param n Input normal.
       *
       * @return Two orthonormal tangent directions spanning the tangent plane.
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

        // Normalize u and v to have unit length.
        const double len_u = std::sqrt(CGAL::to_double(u.squared_length()));
        const double len_v = std::sqrt(CGAL::to_double(v.squared_length()));
        if (len_u > 0) u = u / len_u;
        if (len_v > 0) v = v / len_v;

        return {u, v};
      }

      Vector_3 box_normal(const Point_3& p) const {
        int ix, iy, iz;

        if (!to_grid_coords(p, ix, iy, iz)) {
          return CGAL::NULL_VECTOR;
        }

        return compute_cell_normal(ix, iy, iz);
      }

    private:
    /**
     * @brief Clears all per-cell accumulators.
     */
    void clear() {
      for (Cell& c : cells_) {
        c.point_ids.clear();
        c.normal_sum = CGAL::NULL_VECTOR;
        c.normal_count = 0;
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
      cells_.resize(nx_ * ny_ * nz_);
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
    bool valid_coords(int ix, int iy, int iz) const {
      return ix >= 0 && iy >= 0 && iz >= 0 &&
            static_cast<std::size_t>(ix) < nx_ &&
            static_cast<std::size_t>(iy) < ny_ &&
            static_cast<std::size_t>(iz) < nz_;
    }

    /**
     * @brief Maps a 3D point to integer grid coordinates.
     *
     * @param p   Input point.
     * @param ix  Output cell index along x.
     * @param iy  Output cell index along y.
     * @param iz  Output cell index along z.
     *
     * @return `true` if the point is inside the grid bounds.
     */
    bool to_grid_coords(const Point_3& p, int& ix, int& iy, int& iz) const {

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

      ix = (std::min)(ix, static_cast<int>(nx_) - 1);
      iy = (std::min)(iy, static_cast<int>(ny_) - 1);
      iz = (std::min)(iz, static_cast<int>(nz_) - 1);

      return valid_coords(ix, iy, iz);
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
     * @brief Inserts a single point and its normal into the grid.
     *
     * @param point_id Index of the point in the input cloud.
     * @param p        Point position.
     * @param n        Point normal.
     */
    void insert(Index point_id, const Point_3& p, const Vector_3& n) {
      int ix, iy, iz;
      if (!to_grid_coords(p, ix, iy, iz)) {
        return;
      }

      const std::size_t idx = flat_index(ix, iy, iz);
      cells_[idx].add_point(point_id, n);
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
    FT max_splat_radius_ = 0;


    /*
    //////////////////// DEBUG /////////////////////
    */
    public:
    /**
     * @brief Writes the input point cloud and normals to a PLY file.
     *
     * @param filename Output PLY filename.
     *
     * @return `true` if the file was written successfully.
     */
    bool write_point_cloud_ply(const std::string& filename) const
    {
      std::ofstream out(filename);
      if (!out) {
        std::cerr << "Error: cannot open " << filename << " for writing.\n";
        return false;
      }

      out << "ply\n";
      out << "format ascii 1.0\n";
      out << "element vertex " << points_.size() << "\n";
      out << "property float x\n";
      out << "property float y\n";
      out << "property float z\n";
      out << "property float nx\n";
      out << "property float ny\n";
      out << "property float nz\n";
      out << "end_header\n";

      for (std::size_t i = 0; i < points_.size(); ++i) {
        const Point_3& p = points_[i];
        const Vector_3 n = (i < normals_.size()) ? normals_[i] : CGAL::NULL_VECTOR;

        out << std::setprecision(17)
            << CGAL::to_double(p.x()) << ' '
            << CGAL::to_double(p.y()) << ' '
            << CGAL::to_double(p.z()) << ' '
            << CGAL::to_double(n.x()) << ' '
            << CGAL::to_double(n.y()) << ' '
            << CGAL::to_double(n.z()) << '\n';
      }

      return true;
    }

    /**
     * @brief Writes all grid vertices to a PLY file for visualization.
     *
     * @param filename Output PLY filename.
     *
     * @return `true` if the file was written successfully.
     */
    bool write_grid_vertices_ply(const std::string& filename) const
    {
      std::ofstream out(filename);
      if (!out) {
        std::cerr << "Error: cannot open " << filename << " for writing.\n";
        return false;
      }

      std::vector<Point_3> vertices;
      vertices.reserve((nx_ + 1) * (ny_ + 1) * (nz_ + 1));

      for (std::size_t ix = 0; ix <= nx_; ++ix) {
        for (std::size_t iy = 0; iy <= ny_; ++iy) {
          for (std::size_t iz = 0; iz <= nz_; ++iz) {
            const FT x = min_x_ + FT(ix) * box_size_;
            const FT y = min_y_ + FT(iy) * box_size_;
            const FT z = min_z_ + FT(iz) * box_size_;
            vertices.emplace_back(x, y, z);
          }
        }
      }

      out << "ply\n";
      out << "format ascii 1.0\n";
      out << "element vertex " << vertices.size() << "\n";
      out << "property float x\n";
      out << "property float y\n";
      out << "property float z\n";
      out << "end_header\n";

      for (const Point_3& p : vertices) {
        out << std::setprecision(17)
            << CGAL::to_double(p.x()) << ' '
            << CGAL::to_double(p.y()) << ' '
            << CGAL::to_double(p.z()) << '\n';
      }

      return true;
    }

    /**
     * @brief Writes occupied cell centers and their averaged normals to a PLY file.
     *
     * @param filename Output PLY filename.
     * @param normal_scale Scale factor for visualizing the normal direction.
     *
     * @return `true` if the file was written successfully.
     */
    bool write_cell_centers_and_normals_ply(const std::string& filename,
                                            double normal_scale = 0.25) const
    {
      std::ofstream out(filename);
      if (!out) {
        std::cerr << "Error: cannot open " << filename << " for writing.\n";
        return false;
      }

      std::vector<Point_3> centers;
      std::vector<Vector_3> normals;
      centers.reserve(cells_.size());
      normals.reserve(cells_.size());

      for (std::size_t ix = 0; ix < nx_; ++ix) {
        for (std::size_t iy = 0; iy < ny_; ++iy) {
          for (std::size_t iz = 0; iz < nz_; ++iz) {
            const Cell& c = cells_[flat_index(static_cast<int>(ix),
                                              static_cast<int>(iy),
                                              static_cast<int>(iz))];

            if (c.point_ids.empty()) {
              continue;
            }

            centers.push_back(cell_center(static_cast<int>(ix),
                                          static_cast<int>(iy),
                                          static_cast<int>(iz)));
            normals.push_back(compute_cell_normal(static_cast<int>(ix),
                                                  static_cast<int>(iy),
                                                  static_cast<int>(iz)));
          }
        }
      }

      out << "ply\n";
      out << "format ascii 1.0\n";
      out << "element vertex " << centers.size() * 2 << "\n";
      out << "property float x\n";
      out << "property float y\n";
      out << "property float z\n";
      out << "element edge " << centers.size() << "\n";
      out << "property int vertex1\n";
      out << "property int vertex2\n";
      out << "end_header\n";

      for (std::size_t i = 0; i < centers.size(); ++i) {
        const Point_3& c = centers[i];
        const Vector_3 n = normals[i];

        const Point_3 tip(
          c.x() + FT(normal_scale) * n.x(),
          c.y() + FT(normal_scale) * n.y(),
          c.z() + FT(normal_scale) * n.z()
        );

        out << std::setprecision(17)
            << CGAL::to_double(c.x()) << ' '
            << CGAL::to_double(c.y()) << ' '
            << CGAL::to_double(c.z()) << '\n';

        out << std::setprecision(17)
            << CGAL::to_double(tip.x()) << ' '
            << CGAL::to_double(tip.y()) << ' '
            << CGAL::to_double(tip.z()) << '\n';
      }

      for (std::size_t i = 0; i < centers.size(); ++i) {
        const int v0 = static_cast<int>(2 * i);
        const int v1 = static_cast<int>(2 * i + 1);
        out << v0 << ' ' << v1 << '\n';
      }

      return true;
    }
  };


  /**
   * @brief Builds a splat mesh by growing a front from two initial seeds.
   *
   * This class keeps the halfedge structure in the output mesh, but the region
   * bookkeeping and priorities are intentionally not implemented yet.
   */
  template <typename PointRange, typename NormalRange, typename PolygonMesh>
  class Splat_surface_reconstruction_3 {
    public:
      using Point_3 = typename PointRange::value_type;
      using Kernel = typename Kernel_traits<Point_3>::Kernel;
      using Vector_3 = typename Kernel::Vector_3;
      using Vector_2 = typename Kernel::Vector_2;
      using Point_2  = typename Kernel::Point_2;
      using FT = typename Kernel::FT;
      using Grid = Box_grid<Kernel>;
      using Index = typename Grid::Index;

      using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
      using edge_descriptor = typename boost::graph_traits<PolygonMesh>::edge_descriptor;
      using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
      using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;

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
          second(second)
        {}
      };

      Splat_surface_reconstruction_3(const Grid& grid, PolygonMesh& output_mesh)
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

        seed_from_grid();
      }

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

        while (!candidate_queues_empty()) {
          Candidate cand(Point_3(0,0,0), CGAL::NULL_VECTOR, 0, vertex_descriptor(), vertex_descriptor()); // Dummy initialization
          if (!pop_candidate_from_queue(cand)) {
            break;
          }

          if (has_vertex_near(cand.position)) {
            rejected_near++;
            continue;
          }

          if (!projection_check(cand)) {
            rejected_proj++;
            continue;
          }

          int new_priority = compute_priority(cand);
          if (new_priority != cand.priority) {
            cand.priority = new_priority;
            push_candidate_in_queue(cand);
            continue;
          }

          accepted++;
          vertex_descriptor nv = mesh_.add_vertex(cand.position);
          put(points_pm_, nv, cand.position);
          put(normals_pm_, nv, cand.normal);

          if (!build_graph(cand.first, cand.second, nv, cand.normal)) {
            continue;
          }

          push_candidates_from_vertex(nv);
        }

        // CGAL::Polygon_mesh_processing::triangulate_faces(mesh_);

        std::cout << "Final mesh size: " << num_vertices(mesh_) << " vertices, " << num_edges(mesh_) << " edges,\n";
        std::cout << "  Accepted: " << accepted << ", Rejected (near): " << rejected_near << ", Rejected (proj): " << rejected_proj << "\n";
      }

    private:
      void push_candidate_in_queue(Candidate cand) {
        cand.priority =
          (std::max)(0,
          (std::min)(cand.priority, NUM_PRIORITIES - 1));

        candidate_queues_[cand.priority].push_back(cand);
      }

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

      bool candidate_queues_empty() const {
        for (const auto& q : candidate_queues_) {
          if (!q.empty())
            return false;
        }

        return true;
      }

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

        CGAL_assertion(is_valid_vertex_descriptor(v0_, mesh_));
        CGAL_assertion(is_valid_vertex_descriptor(v1_, mesh_));
        seeded_ = true;
        std::cout << "Mesh seeded with initial vertices. Starting growth...\n";

        push_candidates_from_vertex(v0_);
        push_candidates_from_vertex(v1_);
      }

      void push_candidates_from_vertex(vertex_descriptor nv) {
        const Point_3 p = get(points_pm_, nv);

        // std::vector<Index> nearby_ids = grid_.nearby_point_ids(p, FT(2.0) * grid_.get_box_size());

        for (auto other_v : vertices(mesh_)) {
          if (other_v == nv) {
            continue;
          }

          const Point_3 q = get(points_pm_, other_v);

          const FT d2 = CGAL::squared_distance(p, q);

          if (d2 > FT(4) * grid_.get_box_size() * grid_.get_box_size()) {
            continue;
          }

          std::vector<Candidate> candidates = generate_candidates_from_parents(nv, other_v);

          for (const Candidate& c : candidates) {
            push_candidate_in_queue(c);
          }
        }
      }

      std::vector<Candidate> generate_candidates_from_parents(vertex_descriptor parent_a,
                                                              vertex_descriptor parent_b) const {
        std::vector<Candidate> candidates;

        const Point_3 pa = get(points_pm_, parent_a);
        const Point_3 pb = get(points_pm_, parent_b);

        const FT d = grid_.get_box_size();
        const FT d2 = CGAL::squared_distance(pa, pb);

        // Circle radius from the two-parent construction.
        const FT circle_r = CGAL::sqrt((d * d) - (d2 / FT(4)));

        FT search_radius = circle_r + grid_.get_max_splat_radius();

        const Point_3 mid = Point_3(
          (pa.x() + pb.x()) / FT(2),
          (pa.y() + pb.y()) / FT(2),
          (pa.z() + pb.z()) / FT(2));

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

      bool has_vertex_near(const Point_3& p) const
      {
        const FT tol = grid_.get_box_size() * FT(0.5);
        const FT tol2 = tol * tol;

        for (auto vd : vertices(mesh_)) {
          if (CGAL::squared_distance(get(points_pm_, vd), p) < tol2) {
            return true;
          }
        }
        return false;
      }

      halfedge_descriptor create_edge_between(vertex_descriptor a,
                                              vertex_descriptor b) {
        edge_descriptor e = add_edge(mesh_);
        halfedge_descriptor h  = halfedge(e, mesh_);
        halfedge_descriptor ho = opposite(h, mesh_);

        set_target(h,  b, mesh_);   // a -> b
        set_target(ho, a, mesh_);   // b -> a

        if (halfedge(a, mesh_) == boost::graph_traits<PolygonMesh>::null_halfedge()) {
          set_halfedge(a, h, mesh_);
        }
        if (halfedge(b, mesh_) == boost::graph_traits<PolygonMesh>::null_halfedge()) {
          set_halfedge(b, ho, mesh_);
        }

        return h;
      }

      halfedge_descriptor ensure_halfedge(vertex_descriptor from,
                                          vertex_descriptor to) {

        halfedge_descriptor h = find_halfedge(from, to);
        if (h != boost::graph_traits<PolygonMesh>::null_halfedge()) {
          return h;
        }

        edge_descriptor e = add_edge(mesh_);
        halfedge_descriptor h0 = halfedge(e, mesh_);
        halfedge_descriptor h1 = opposite(h0, mesh_);

        // h0 : from -> to
        set_target(h0, to, mesh_);
        set_target(h1, from, mesh_);

        if (halfedge(from, mesh_) == boost::graph_traits<PolygonMesh>::null_halfedge()) {
          set_halfedge(from, h0, mesh_);
        }
        if (halfedge(to, mesh_) == boost::graph_traits<PolygonMesh>::null_halfedge()) {
          set_halfedge(to, h1, mesh_);
        }

        return h0;
      }

      halfedge_descriptor find_halfedge(vertex_descriptor from,
                                        vertex_descriptor to) const {
        for (halfedge_descriptor h : halfedges(mesh_)) {
          if (source(h, mesh_) == from &&
              target(h, mesh_) == to)
          {
            return h;
          }
        }

        return boost::graph_traits<PolygonMesh>::null_halfedge();
      }

      bool build_graph(const vertex_descriptor v0,
                        const vertex_descriptor v1,
                        const vertex_descriptor nv,
                        const Vector_3& normal) {

        const Point_3 p0 = get(points_pm_, v0);
        const Point_3 p1 = get(points_pm_, v1);
        const Point_3 p2 = get(points_pm_, nv);

        // Create or reuse v0-nv
        halfedge_descriptor h0 = ensure_halfedge(v0, nv);   // v0 -> nv
        halfedge_descriptor h0o = opposite(h0, mesh_);      // nv -> v0

        // Create or reuse nv-v1
        halfedge_descriptor h1o = ensure_halfedge(v1, nv);   // v1 -> nv
        halfedge_descriptor h1  = opposite(h1o, mesh_);      // nv -> v1

        // Make sure nv has a representative halfedge.
        if (halfedge(nv, mesh_) == boost::graph_traits<PolygonMesh>::null_halfedge()) {
          set_halfedge(nv, h0o, mesh_);
        }

        // Find directed edge v0 -> v1
        halfedge_descriptor h01 = find_halfedge(v0, v1);
        if (h01 == boost::graph_traits<PolygonMesh>::null_halfedge()) {
          // No base edge, so no face cycle yet.
          return true;
        }

        const halfedge_descriptor h10 = opposite(h01, mesh_);
        const Vector_3 tri_normal = CGAL::cross_product(p0 - p2, p1 - p2);

        auto f = add_face(mesh_);
        if (f == boost::graph_traits<PolygonMesh>::null_face()) {
          return false;
        }

        if (normal == CGAL::NULL_VECTOR || tri_normal * normal >= FT(0)) {
          // v0 -> nv -> v1 -> v0
          set_next(h0,  h1,  mesh_);
          set_next(h1,  h10, mesh_);
          set_next(h10, h0,  mesh_);

          set_face(h0,  f, mesh_);
          set_face(h1,  f, mesh_);
          set_face(h10, f, mesh_);
          set_halfedge(f, h0, mesh_);
        } else {
          // v0 -> v1 -> nv -> v0
          set_next(h01, h1o, mesh_);
          set_next(h1o, h0o, mesh_);
          set_next(h0o, h01, mesh_);

          set_face(h01, f, mesh_);
          set_face(h1o, f, mesh_);
          set_face(h0o, f, mesh_);
          set_halfedge(f, h01, mesh_);
        }

        return true;
      }

      bool projection_check(const Candidate& cand) const {
        const Point_3 pa = get(points_pm_, cand.first);
        const Point_3 pb = get(points_pm_, cand.second);
        const Point_3 pc = cand.position;

        // Use the box normal at the candidate position as the projection direction.
        Vector_3 n = grid_.box_normal(cand.position);

        auto local_frame = grid_.compute_local_tangent_frame(n);
        Vector_3 u = local_frame[0];
        Vector_3 v = local_frame[1];

        auto project = [&](const Point_3& p, const Point_3& origin) -> Point_2 {
          const Vector_3 d = p - origin;
          return Point_2(d * u, d * v);
        };

        const Point_2 A = project(pa, pc);
        const Point_2 B = project(pb, pc);
        const Point_2 C = project(pc, pc);

        // New projected edges created if this candidate is accepted.
        const Point_2 new_e1_s = A;
        const Point_2 new_e1_t = C;
        const Point_2 new_e2_s = B;
        const Point_2 new_e2_t = C;

        for (auto e : edges(mesh_)) {
          halfedge_descriptor h = halfedge(e, mesh_);
          vertex_descriptor s = source(h, mesh_);
          vertex_descriptor t = target(h, mesh_);

          // Ignore edges incident to the candidate's parents.
          if (s == cand.first || s == cand.second ||
              t == cand.first || t == cand.second) {
            continue;
          }

          const Point_3 q0 = get(points_pm_, s);
          const Point_3 q1 = get(points_pm_, t);

          if (CGAL::squared_distance(q0, pc) >= 4*grid_.get_box_size()*grid_.get_box_size() ||
              CGAL::squared_distance(q1, pc) >= 4*grid_.get_box_size()*grid_.get_box_size()) {
            continue;
          }

          const Point_2 E0 = project(q0, pc);
          const Point_2 E1 = project(q1, pc);

          const CGAL::Segment_2<Kernel> s1(new_e1_s, new_e1_t);
          const CGAL::Segment_2<Kernel> s2(E0, E1);

          auto strictly_intersect_2d = [&](const Point_2& a,
                                          const Point_2& b,
                                          const Point_2& c,
                                          const Point_2& d) -> bool {
            const CGAL::Segment_2<Kernel> s1(a, b);
            const CGAL::Segment_2<Kernel> s2(c, d);

            if (!CGAL::do_intersect(s1, s2)) {
              return false;
            }

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
            return true;
          };

          if (strictly_intersect_2d(new_e1_s, new_e1_t, E0, E1)) {
            return false;
          }
          if (strictly_intersect_2d(new_e2_s, new_e2_t, E0, E1)) {
            return false;
          }
        }

        return true;
      }

      std::size_t incident_vertex_degree(vertex_descriptor v) const {
        std::size_t deg = 0;

        for (auto e : edges(mesh_)) {
          halfedge_descriptor h = halfedge(e, mesh_);
          const vertex_descriptor t = target(h, mesh_);

          if (t == v) {
            ++deg;
          }
        }

        return deg;
      }

      bool walk_finds_other_parent_from_vertex(vertex_descriptor start_v,
                                               vertex_descriptor target_v) const {

        halfedge_descriptor h0 = halfedge(start_v, mesh_);
        if (h0 == boost::graph_traits<PolygonMesh>::null_halfedge()) {
          return false;
        }

        halfedge_descriptor h = h0;
        for (std::size_t i = 0; i < window_size_; ++i) {
          const vertex_descriptor s = source(h, mesh_);
          const vertex_descriptor t = target(h, mesh_);

          if (s == target_v || t == target_v) {
            return true;
          }

          const halfedge_descriptor o = opposite(h, mesh_);
          if (o == boost::graph_traits<PolygonMesh>::null_halfedge()) {
            return false;
          }

          const halfedge_descriptor n = next(o, mesh_);
          if (n == boost::graph_traits<PolygonMesh>::null_halfedge()) {
            return false;
          }

          h = n;
        }

        return false;
      }

      bool joins_two_borders(const Candidate& cand) const {
        const bool a_finds_b = walk_finds_other_parent_from_vertex(cand.first, cand.second);
        const bool b_finds_a = walk_finds_other_parent_from_vertex(cand.second, cand.first);

        return !a_finds_b && !b_finds_a;
      }

      int compute_priority(const Candidate& cand) const {
        if (incident_vertex_degree(cand.first) == 1 || incident_vertex_degree(cand.second) == 1) {
          return 2; // highest priority
        }

        else if (joins_two_borders(cand)) {
          return 1; // medium priority
        }

        return 0;   // default for now
      }

    private:
      const Grid& grid_;
      PolygonMesh& mesh_;

      typename PolygonMesh:: template Property_map<vertex_descriptor,Point_3> points_pm_;
      typename PolygonMesh:: template Property_map<vertex_descriptor,Vector_3> normals_pm_;

      typename PolygonMesh::template Property_map<edge_descriptor, FT> edge_weight_pm_;

      vertex_descriptor v0_{boost::graph_traits<PolygonMesh>::null_vertex()};
      vertex_descriptor v1_{boost::graph_traits<PolygonMesh>::null_vertex()};
      bool seeded_ = false;

      static constexpr int NUM_PRIORITIES = 3; // 0 is lowest, NUM_PRIORITIES-1 is highest
      std::array<std::deque<Candidate>, NUM_PRIORITIES> candidate_queues_; // candidate_queues_[p] contains all candidates with priority p

      const std::size_t window_size_ = 8;
      // in the class members
  };

} // namespace CGAL

#endif // CGAL_SPLAT_SURFACE_RECONSTRUCTION_H
