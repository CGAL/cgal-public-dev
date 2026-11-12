#ifndef CGAL_INTERNAL_SPLAT_SURFACE_RECONSTRUCTION_3_BOX_GRID_H
#define CGAL_INTERNAL_SPLAT_SURFACE_RECONSTRUCTION_3_BOX_GRID_H

#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Kernel/global_functions_3.h>

// Include the standard headers used by Box_grid.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

namespace CGAL {

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

  public:

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

    FT get_min_x() const {
      return min_x_;
    }

    FT get_min_y() const {
      return min_y_;
    }

    FT get_min_z() const {
      return min_z_;
    }

    FT get_max_x() const {
      return max_x_;
    }

    FT get_max_y() const {
      return max_y_;
    }

    FT get_max_z() const {
      return max_z_;
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

    bool is_point_inside_grid(const Point_3& p) const {
      int ix, iy, iz;
      return to_grid_coordinates(p, ix, iy, iz);
    }

    /// @brief Returns the total number of grid cells.
    std::size_t number_of_cells() const {
      return nx_ * ny_ * nz_;
    }

    /**
     * @brief Compute the distance from a point to the surface of a cell's bounding box.
     * @param p The point to measure the distance from.
     * @param ix The x-coordinate of the cell in the grid.
     * @param iy The y-coordinate of the cell in the grid.
     * @param iz The z-coordinate of the cell in the grid.
     * @return The distance from the point to the surface of the cell's bounding box.
     */
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

    Point_3 cell_min_corner(int ix, int iy, int iz) const {
      const FT x = min_x_ + FT(ix) * box_size_;
      const FT y = min_y_ + FT(iy) * box_size_;
      const FT z = min_z_ + FT(iz) * box_size_;
      return Point_3(x, y, z);
    }

    Point_3 cell_max_corner(int ix, int iy, int iz) const {
      const FT x = min_x_ + FT(ix + 1) * box_size_;
      const FT y = min_y_ + FT(iy + 1) * box_size_;
      const FT z = min_z_ + FT(iz + 1) * box_size_;
      return Point_3(x, y, z);
    }

  private:

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
    }

    FT box_size_;
    FT min_x_ = 0, max_x_ = 0;
    FT min_y_ = 0, max_y_ = 0;
    FT min_z_ = 0, max_z_ = 0;

    std::size_t nx_ = 0;
    std::size_t ny_ = 0;
    std::size_t nz_ = 0;
};


} // namespace CGAL

#endif