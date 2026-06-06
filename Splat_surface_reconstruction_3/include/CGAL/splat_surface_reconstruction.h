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
        std::exit(0);
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

            // If no points contributed to this cell, we leave the normal as NULL_VECTOR.
            if (c.normal_count == 0) {
              block_normals.push_back(CGAL::NULL_VECTOR);
              continue;
            }

            // Average the normals of points in this cell to get the block normal.
            const FT inv = FT(1) / FT(c.normal_count);
            Vector_3 n = c.normal_sum * inv;

            // Normalize the block normal to have unit length.
            const double len2 = CGAL::to_double(n.squared_length());
            block_normals.push_back(n * (1.0 / std::sqrt(len2)));
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
      splat_sizes_.resize(points_.size(), 2*box_size_);

      for (Index i = 0; i < points_.size(); ++i) {
        int cx, cy, cz;
        if (!to_grid_coords(points_[i], cx, cy, cz)) {
          std::cerr << "Warning: point " << i << " is out of grid bounds, skipping splat size estimation.\n";
          continue;
        }

        // Collect all points in the same box as p_i.
        std::vector<Index> local_ids;
        const Cell* c = cell(cx, cy, cz);
        if (c == nullptr || c->point_ids.empty()) {
          std::cerr << "Warning: no points found in the same cell for point " << i << ", skipping splat size estimation.\n";
          continue;
        }
        local_ids = c->point_ids;

        // Need at least a few neighbors to define a local triangulation.
        if (local_ids.size() < 3) {
          std::cerr << "Warning: not enough neighbors for point " << i << ", skipping splat size estimation.\n";
          continue;
        }

        // Build a local tangent frame at p_i from its normal.
        std::vector<Vector_3> local_frame = compute_local_tangent_frame(normals_[i]);
        Vector_3 u = local_frame[0];
        Vector_3 v = local_frame[1];

        // Project local neighbors onto the tangent plane and compute their 2D coordinates
        // Also input in the Delaunay Triangulation Kernel
        CGAL::Delaunay_triangulation_2<Kernel> DT;
        DT.insert(Point_2(0, 0)); // insert the center point itself at the origin of the local frame
        for (Index neighbor_id : local_ids) {
          if (neighbor_id == i) continue; // skip the center point
          Vector_3 vec = points_[neighbor_id] - points_[i];
          double x = CGAL::to_double(vec * u);
          double y = CGAL::to_double(vec * v);
          Point_2 point_2d(x, y);

          DT.insert(point_2d);
        }

        //Find the circumcenter of each triangle in Delaunay triangulation
        std::vector<FT> circumcenter_distances;
        for (auto face = DT.finite_faces_begin(); face != DT.finite_faces_end(); ++face) {
          Point_2 c = CGAL::circumcenter(face->vertex(0)->point(), face->vertex(1)->point(), face->vertex(2)->point());
          double dist = std::sqrt(CGAL::to_double(c.x()*c.x() + c.y()*c.y())); // distance from circumcenter to the origin (which is the projection of p_i)
          circumcenter_distances.push_back(FT(dist));
        }

        splat_sizes_[i] = *std::max_element(circumcenter_distances.begin(), circumcenter_distances.end());
      }

      // //// WARNING: HACK to prevent large splats 
      // for (Index i = 0; i < splat_sizes_.size(); ++i) {
      //   if (splat_sizes_[i] > FT(0.9)) {
      //     splat_sizes_[i] = 2*global_spacing;
      //   }
      // }

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
        if (splat_sizes_[i] > box_size_ && i < normals_.size() && normals_[i] != CGAL::NULL_VECTOR)
        {
          auto frame = compute_local_tangent_frame(normals_[i]);
          const Vector_3& u = frame[0];

          const Point_3 seed1 = points_[i] + FT(0.5) * box_size_ * u;
          const Point_3 seed2 = points_[i] - FT(0.5) * box_size_ * u;

          std::cout<<"Found initial seed at point index " << i << " with splat size " << splat_sizes_[i] << "." << std::endl;
          return Initial_seed{seed1, seed2, i};
        }
      }

      std::cerr << "Warning: no suitable initial seed found with splat size larger than box size.\n";
      exit(1);
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

      const int r = std::max(
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
     * @brief Intersects a sphere with one splat disk.
     *
     * @param sphere_center Sphere center.
     * @param sphere_radius Sphere radius.
     * @param splat_center Splat center.
     * @param splat_normal Splat normal.
     * @param splat_radius Splat disk radius.
     *
     * @return Up to two intersection points on the splat disk.
     */
    std::vector<Point_3> intersect_sphere_with_splat(const Point_3& sphere_center,
                                                    FT sphere_radius,
                                                    const Point_3& splat_center,
                                                    const Vector_3& splat_normal,
                                                    FT splat_radius) const {
      std::vector<Point_3> result;
      Vector_3 n = splat_normal;

      const double nlen2 = CGAL::to_double(n.squared_length());
      if (nlen2 <= 0.0) {
        return result;
      }
      n = n * (1.0 / std::sqrt(nlen2));

      const Vector_3 diff = sphere_center - splat_center;
      const FT signed_dist = diff * n;

      if (std::abs(CGAL::to_double(signed_dist)) > CGAL::to_double(sphere_radius)) {
        return result;
      }

      const FT circle_r2 = sphere_radius * sphere_radius - signed_dist * signed_dist;
      if (circle_r2 <= FT(0)) {
        return result;
      }

      const FT circle_r = CGAL::sqrt(circle_r2);

      const Point_3 circle_center = sphere_center - signed_dist * n;

      Vector_3 dir = splat_center - circle_center;
      dir = dir - (dir * n) * n;

      double dir_len2 = CGAL::to_double(dir.squared_length());
      if (dir_len2 <= 1e-18) {
        auto frame = compute_local_tangent_frame(splat_normal);
        dir = frame[0];
        dir_len2 = CGAL::to_double(dir.squared_length());
        if (dir_len2 <= 0.0) {
          return result;
        }
      }

      const Vector_3 u = dir * (1.0 / std::sqrt(dir_len2));

      const Point_3 p1 = circle_center + u * circle_r;
      const Point_3 p2 = circle_center - u * circle_r;

      const FT r2 = splat_radius * splat_radius;
      if (CGAL::squared_distance(p1, splat_center) <= r2) {
        result.push_back(p1);
      }
      if (CGAL::squared_distance(p2, splat_center) <= r2) {
        result.push_back(p2);
      }

      return result;
    }


    private:
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
    using FT = typename Kernel::FT;
    using Grid = Box_grid<Kernel>;
    using Index = typename Grid::Index;

    using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
    using edge_descriptor = typename boost::graph_traits<PolygonMesh>::edge_descriptor;
    using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;

    struct Candidate {
    Point_3 position;
    Vector_3 normal;
    Index splat_id;
    vertex_descriptor first;
    vertex_descriptor second;

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

      while (!candidate_queue_.empty()) {
        Candidate cand = candidate_queue_.front();
        candidate_queue_.pop();

        if (has_vertex_near(cand.position)) {
          continue;
        }

        vertex_descriptor nv = mesh_.add_vertex(cand.position);
        put(points_pm_, nv, cand.position);
        put(normals_pm_, nv, cand.normal);

        (void)connect_vertices(cand.first, nv);
        (void)connect_vertices(cand.second, nv);
        // Create the triangle (parent_a, parent_b, new vertex).
        auto f = mesh_.add_face(cand.first, cand.second, nv);
        if (f == boost::graph_traits<PolygonMesh>::null_face()) {
          // Try reversed orientation if needed.
          f = mesh_.add_face(cand.second, cand.first, nv);
        }


        push_candidates_from_parents(nv, cand.first);
        push_candidates_from_parents(nv, cand.second);

        std::cout << "Current mesh size: " << num_vertices(mesh_) << " vertices, " << num_edges(mesh_) << " edges,\n";
      }
    }

    private:
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

      (void)connect_vertices(v0_, v1_);

      seeded_ = true;
      std::cout << "Mesh seeded with initial vertices. Starting growth...\n";
      push_candidates_from_parents(v0_, v1_);
    }

    void push_candidates_from_parents(vertex_descriptor parent_a,
                                      vertex_descriptor parent_b) {
      std::vector<Candidate> candidates = generate_candidates_from_parents(parent_a, parent_b);
      std::cout << "Generated " << candidates.size() << " candidate points." << std::endl;
      for (const Candidate& c : candidates) {
        candidate_queue_.push(c);
      }
    }

    std::vector<Candidate> generate_candidates_from_parents(vertex_descriptor parent_a,
                                                            vertex_descriptor parent_b) const {
      std::vector<Candidate> candidates;

      const Point_3 pa = get(points_pm_, parent_a);
      const Point_3 pb = get(points_pm_, parent_b);

      const FT d = grid_.get_box_size();
      const FT d2 = CGAL::squared_distance(pa, pb);

      if (d2 >= FT(4) * d * d) {
        return candidates;
      }

      const FT sphere_r2 = d * d - d2 / FT(4);
      if (sphere_r2 <= FT(0)) {
        return candidates;
      }

      const FT sphere_r = CGAL::sqrt(sphere_r2);

      const Point_3 mid(
        (pa.x() + pb.x()) / FT(2),
        (pa.y() + pb.y()) / FT(2),
        (pa.z() + pb.z()) / FT(2));

      FT search_radius = sphere_r + grid_.get_max_splat_radius();
      search_radius = std::min(search_radius, grid_.get_box_size()); // WARNING: HACK

      std::vector<Index> nearby = grid_.nearby_point_ids(mid, search_radius);

      const std::vector<Point_3>& points = grid_.points();
      const std::vector<Vector_3>& normals = grid_.normals();
      const std::vector<FT>& splat_sizes = grid_.splat_sizes();

      for (Index sid : nearby) {
        if (sid >= points.size() || sid >= normals.size()) {
          continue;
        }

        const FT splat_radius = (sid < splat_sizes.size() && splat_sizes[sid] > FT(0))
                              ? splat_sizes[sid]
                              : 2*grid_.get_box_size();

        const FT broad_r = sphere_r + splat_radius;
        if (CGAL::squared_distance(mid, points[sid]) > broad_r * broad_r) {
          continue;
        }

        const std::vector<Point_3> hits =
          grid_.intersect_sphere_with_splat(mid,
                                            sphere_r,
                                            points[sid],
                                            normals[sid],
                                            splat_radius);

        for (const Point_3& p : hits) {
          candidates.push_back(Candidate(p, normals[sid], sid, parent_a, parent_b));
        }
      }

      return candidates;
    }

    bool has_vertex_near(const Point_3& p) const
    {
      const FT tol = grid_.get_box_size() * FT(0.8); // WARNING: HACK
      const FT tol2 = tol * tol;

      for (auto vd : vertices(mesh_)) {
        if (CGAL::squared_distance(get(points_pm_, vd), p) <= tol2) {
          return true;
        }
      }
      return false;
    }

    edge_descriptor connect_vertices(vertex_descriptor a, vertex_descriptor b) {
      edge_descriptor e = add_edge(mesh_);
      halfedge_descriptor h = halfedge(e, mesh_);

      set_target(h, b, mesh_);
      set_target(opposite(h, mesh_), a, mesh_);

      return e;
    }

  private:
    const Grid& grid_;
    PolygonMesh& mesh_;

    typename PolygonMesh:: template Property_map<vertex_descriptor,Point_3> points_pm_;
    typename PolygonMesh:: template Property_map<vertex_descriptor,Vector_3> normals_pm_;

    std::queue<Candidate> candidate_queue_;

    vertex_descriptor v0_{boost::graph_traits<PolygonMesh>::null_vertex()};
    vertex_descriptor v1_{boost::graph_traits<PolygonMesh>::null_vertex()};
    bool seeded_ = false;
  };

} // namespace CGAL

#endif // CGAL_SPLAT_SURFACE_RECONSTRUCTION_H
