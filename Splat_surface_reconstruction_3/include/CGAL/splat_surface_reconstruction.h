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

#ifndef CGAL_SPLAT_SURFACE_RECONSTRUCTION_H
#define CGAL_SPLAT_SURFACE_RECONSTRUCTION_H

#include <CGAL/license/Splat_surface_reconstruction_3.h>
#include <CGAL/property_map.h>
// #include <CGAL/Splat_surface_reconstruction_3/internal/Box_grid.h>
#include <CGAL/Kernel_traits.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

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
  template <typename PointRange,typename NormalRange, typename PolygonMesh>
  bool splat_surface_reconstruction(const PointRange& points,
                               const NormalRange& normals,
                               PolygonMesh& output_mesh,
                               double spacing) {
    typedef typename PointRange::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel Kernel;
    typedef typename Kernel::Sphere_3 Sphere;
    typedef typename Kernel::FT FT;

    return true;
  }

  template <typename Kernel_>
  class Box_grid {
  
    public:
    using Kernel   = Kernel_;
    using FT       = typename Kernel::FT;
    using Point_3  = typename Kernel::Point_3;
    using Vector_3 = typename Kernel::Vector_3;
    using Index    = std::size_t;

    public:
    struct Cell
    {
      std::vector<Index> point_ids;
      Vector_3 normal_sum = CGAL::NULL_VECTOR;
      std::size_t normal_count = 0;

      void add_point(Index pid, const Vector_3& n)
      {
        point_ids.push_back(pid);
        if (n != CGAL::NULL_VECTOR) {
          normal_sum = normal_sum + n;
          ++normal_count;
        }
      }
    };

    Box_grid(FT box_size,
             FT min_x, FT max_x,
             FT min_y, FT max_y,
             FT min_z, FT max_z) {
      
      box_size_ = box_size;
      min_x_ = min_x;
      max_x_ = max_x;
      min_y_ = min_y;
      max_y_ = max_y;
      min_z_ = min_z;
      max_z_ = max_z;
      initialize_grid();
    }

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

  private:
    void clear() {
      for (Cell& c : cells_) {
        c.point_ids.clear();
        c.normal_sum = CGAL::NULL_VECTOR;
        c.normal_count = 0;
      }
    }

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

    bool valid_coords(int ix, int iy, int iz) const {
      return ix >= 0 && iy >= 0 && iz >= 0 &&
            static_cast<std::size_t>(ix) < nx_ &&
            static_cast<std::size_t>(iy) < ny_ &&
            static_cast<std::size_t>(iz) < nz_;
    }

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

      ix = std::min(ix, static_cast<int>(nx_) - 1);
      iy = std::min(iy, static_cast<int>(ny_) - 1);
      iz = std::min(iz, static_cast<int>(nz_) - 1);

      return valid_coords(ix, iy, iz);
    }

    std::size_t flat_index(int ix, int iy, int iz) const {
      return (static_cast<std::size_t>(ix) * ny_
            + static_cast<std::size_t>(iy)) * nz_
            + static_cast<std::size_t>(iz);
    }

    void insert(Index point_id, const Point_3& p, const Vector_3& n) {
      int ix, iy, iz;
      if (!to_grid_coords(p, ix, iy, iz)) {
        return;
      }

      const std::size_t idx = flat_index(ix, iy, iz);
      cells_[idx].add_point(point_id, n);
    }

    private:
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
  };

}

#endif // CGAL_SPLAT_SURFACE_RECONSTRUCTION_H
