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

#include <CGAL/license/Splat_surface_reconstruction_3.h>
#include <CGAL/property_map.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

namespace CGAL {

  template <typename PointRange, typename NormalRange, typename PolygonMesh>
  class Splat_surface_reconstruction_3
  {
    public:
    using Point_3 = typename PointRange::value_type;
    using Kernel = typename Kernel_traits<Point_3>::Kernel;
    using Vector_3 = typename Kernel::Vector_3;
    using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
    using edge_descriptor = typename boost::graph_traits<PolygonMesh>::edge_descriptor;
    using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;

    struct Candidate {
      Point_3 position;
      Vector_3 normal;
      vertex_descriptor first, second;

      Candidate(const Point_3& position, const Vector_3& normal, vertex_descriptor first, vertex_descriptor second)
        : position(position), normal(normal), first(first), second(second)
        {}
    };

    Splat_surface_reconstruction_3(const PointRange& points,
                                   const NormalRange& normals,
                                   PolygonMesh& output_mesh,
                                   double spacing)
    : mesh(output_mesh), points_pm(get_property_map(boost::vertex_point, mesh))
    {
      bool created;
      std::tie(normals_pm, created) = mesh.add_property_map<vertex_descriptor, Vector_3>("v:normal", CGAL::NULL_VECTOR);
      if(created) {
        std::cout << "Created vertex normal property map." << std::endl;
      }
      // Assume that the first two points are the ones the algorithm starts with
      vertex_descriptor v0 = add_vertex(mesh);
      put(points_pm, v0, points[0]);
      put(normals_pm, v0, normals[0]);

      vertex_descriptor v1 = add_vertex(mesh);
      put(points_pm, v1, points[1]);
      put(normals_pm, v1, normals[1]);



      Candidate candidate(points[2], normals[2], v0, v1);

      vertex_descriptor nv = add_vertex(mesh);
      put(points_pm, nv, candidate.position);
      put(normals_pm, nv, candidate.normal);

      edge_descriptor fne = add_edge(mesh);
      halfedge_descriptor fnh = halfedge(fne, mesh);
      set_target(fnh, nv, mesh);
      set_target(opposite(fnh, mesh), v0, mesh);

      edge_descriptor nse = add_edge(mesh);
      halfedge_descriptor nsh = halfedge(nse, mesh);
      set_target(nsh, v1, mesh);
      set_target(opposite(nsh, mesh), nv, mesh);

      set_halfedge(nv, fnh, mesh);

      // if candidate.first ,  nv ,  and candidate.second  perform a left turn in the plane defined by the normals of the three points,
      // then we orient the face counterclockwise, otherwise we orient it clockwise

      if(CGAL::orientation(get(points_pm, candidate.first),
                           get(points_pm, nv),
                           get(points_pm, candidate.second),
                           get(points_pm, nv) + candidate.normal) == CGAL::LEFT_TURN) {
        set_next(fnh, nsh, mesh);
        set_next(opposite(nsh,mesh), opposite(fnh, mesh), mesh);
      }else{
        // todo
      }

      if(halfedge(candidate.first, mesh) == boost::graph_traits<PolygonMesh>::null_halfedge())
      {
        // the first edge becomes a circular list of two  halfedges
        set_halfedge(candidate.first, opposite(fnh,mesh), mesh);
        set_next(opposite(fnh, mesh), fnh, mesh);
      } else {
        // we add the new halfedge to the circular list of halfedges around candidate.first
        // Note that the edges have no particular order. We have to order them geometrically
        // around the vertex candidate.first to ensure the correct orientation of the faces.
        set_next(opposite(fnh, mesh), halfedge(candidate.first, mesh), mesh);
      }
    }


    PolygonMesh& mesh;
    typename PolygonMesh:: template Property_map<vertex_descriptor,Point_3> points_pm;
    typename PolygonMesh:: template Property_map<vertex_descriptor,Vector_3> normals_pm;
  };


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
    std::vector<FT> estimate_individual_splat_sizes(FT global_spacing) const {
      std::vector<FT> splat_sizes(points_.size(), global_spacing);

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
          splat_sizes[i] = global_spacing;
          continue;
        }

        // Build a local tangent frame at p_i from its normal.
        std::vector<Vector_3> local_frame = compute_local_tangent_frame(points_[i], normals_[i]);
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

        splat_sizes[i] = *std::max_element(circumcenter_distances.begin(), circumcenter_distances.end());
      }

      return splat_sizes;
    }

  private:
    /**
     * @brief Computes an orthonormal tangent frame at a point.
     *
     * @param n Input normal.
     *
     * @return Two orthonormal tangent directions spanning the tangent plane.
     */
    std::vector<Vector_3> compute_local_tangent_frame(const Point_3& p, const Vector_3& n) const {
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

    private:
    Point_3 cell_center(int ix, int iy, int iz) const
    {
      const FT x = min_x_ + (FT(ix) + FT(0.5)) * box_size_;
      const FT y = min_y_ + (FT(iy) + FT(0.5)) * box_size_;
      const FT z = min_z_ + (FT(iz) + FT(0.5)) * box_size_;
      return Point_3(x, y, z);
    }

    Vector_3 compute_cell_normal(int ix, int iy, int iz) const
    {
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
  };

} // namespace CGAL

#endif // CGAL_SPLAT_SURFACE_RECONSTRUCTION_H
