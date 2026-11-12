#ifndef CGAL_INTERNAL_SPLAT_SURFACE_RECONSTRUCTION_3_MESH_SURFACE_APPROXIMATION_H
#define CGAL_INTERNAL_SPLAT_SURFACE_RECONSTRUCTION_3_MESH_SURFACE_APPROXIMATION_H

#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

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

template <typename Kernel_, typename Mesh_>
class Mesh_surface_approximation : public Surface_approximation<Kernel_> {
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
    using Mesh = Mesh_;
    using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;
    using edge_descriptor = typename boost::graph_traits<Mesh>::edge_descriptor;
    using halfedge_descriptor = typename boost::graph_traits<Mesh>::halfedge_descriptor;
    using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;

    struct Triangle {
      Point_3 a, b, c;
      Vector_3 normal;

      Triangle(const Point_3& a_, const Point_3& b_, const Point_3& c_, const Vector_3& normal_)
        : a(a_), b(b_), c(c_), normal(normal_) {}
    };

    Mesh_surface_approximation(
      const Mesh& mesh,
      const Box_grid<Kernel>& grid
    ) :
      mesh_(mesh),
      grid_(grid)
    {
      cells_.resize(grid_.number_of_cells());

      std::vector<Point_3> vertices;

      // Add each face to the grid boxes
      for (const face_descriptor f : mesh_.faces())
      {
        // Get face normal
        const Vector_3 normal = CGAL::Polygon_mesh_processing::compute_face_normal(f, mesh_);

        vertices.clear();

        // Get face vertices and points
        for (const vertex_descriptor v :
            CGAL::vertices_around_face(mesh_.halfedge(f), mesh_))
        {
          const Point_3& p = mesh_.point(v);
          vertices.push_back(p);
        }

        if (vertices.size() == 3) {
          triangles_.emplace_back(vertices[0], vertices[1], vertices[2], normal);
        }
        else {
          std::cout << "Face with " << vertices.size() << " vertices encountered." << std::endl;
        }
      }

      for (std::size_t i = 0; i < triangles_.size(); ++i) {
        insert(i);
      }

      filter_triangles_by_normals();

      compute_block_normals();
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
    Initial_seed get_initial_seed() const override
    {
      const FT box_size = grid_.get_box_size();

      for (const Triangle& triangle : triangles_)
      {
        const Point_3 points[] = {
          triangle.a, triangle.b, triangle.c
        };

        // Find the longest edge.
        FT max_length_sq = FT(0);
        Point_3 p1, p2;

        for (int i = 0; i < 3; ++i)
        {
          for (int j = i + 1; j < 3; ++j)
          {
            const FT length_sq = CGAL::squared_distance(points[i], points[j]);

            if (length_sq > max_length_sq)
            {
              max_length_sq = length_sq;
              p1 = points[i];
              p2 = points[j];
            }
          }
        }

        if (max_length_sq < box_size * box_size)
          continue;

        // Place two points exactly box_size apart along the edge.
        const Vector_3 direction = (p2 - p1) / CGAL::sqrt(max_length_sq);
        const Point_3 first = p1;
        const Point_3 second = p1 + direction * box_size;

        return { first, second, triangle.normal };
      }

      std::cout << "No triangle contains an edge long enough for the initial seed." << std::endl;
      CGAL_assertion(false);
      return { Point_3(0,0,0), Point_3(0,0,0), CGAL::NULL_VECTOR };
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

      // Iterate over all nearby triangles
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

            const std::vector<Index>& box_triangle_ids = triangle_ids_in_box(
              static_cast<std::size_t>(nx),
              static_cast<std::size_t>(ny),
              static_cast<std::size_t>(nz)
            );

            for (const Index& triangle_idx : box_triangle_ids) {
              const Triangle& triangle = triangles_[triangle_idx];

              std::vector<Point_3> hits = circle_intersects_triangle<Point_3, Vector_3, FT>(
                circle_center,
                circle_normal,
                circle_radius, 
                triangle.a,
                triangle.b,
                triangle.c
              );

              for (const Point_3& hit : hits) {
                result.emplace_back(hit, triangle.normal);
              }
            }
          }
        }
      }

      return result;
    }


    std::vector<Feature_segment> find_feature_segments(FT minimum_angle) override {
      std::vector<Feature_segment> feature_segments;

      // Count feature edges incident to each vertex.
      std::vector<std::size_t> feature_edge_count(mesh_.number_of_vertices(), 0);

      // Accumulate the normals of feature edges at each vertex.
      std::vector<Vector_3> normal_sums(mesh_.number_of_vertices(), CGAL::NULL_VECTOR);

      // Iterate over all edges, find the angle for each and add it to the feature_segments if it exceeds the minimum_angle
      for (const edge_descriptor e : mesh_.edges())
      {
        const halfedge_descriptor h = mesh_.halfedge(e);

        const face_descriptor f1 = mesh_.face(h);
        const face_descriptor f2 = mesh_.face(mesh_.opposite(h));

        // Skip boundary edges (only one incident face).
        if (f1 == Mesh::null_face() || f2 == Mesh::null_face())
          continue;

        const Vector_3 n1 = CGAL::Polygon_mesh_processing::compute_face_normal(f1, mesh_);
        const Vector_3 n2 = CGAL::Polygon_mesh_processing::compute_face_normal(f2, mesh_);

        const FT angle = CGAL::approximate_angle(n1, n2);

        if (angle < minimum_angle)
          continue;

        const vertex_descriptor v1 = mesh_.source(h);
        const vertex_descriptor v2 = mesh_.target(h);
        const Point_3 source = mesh_.point(v1);
        const Point_3 target = mesh_.point(v2);
        const Vector_3 normal = normalize_vec<Vector_3, FT>(n1 + n2);

        feature_segments.emplace_back(source, target, normal);

        ++feature_edge_count[v1];
        ++feature_edge_count[v2];

        normal_sums[v1] = normal_sums[v1] + normal;
        normal_sums[v2] = normal_sums[v2] + normal;
      }

      // Collect vertices incident to at least three feature edges.
      for (const vertex_descriptor v : mesh_.vertices())
      {
        if (feature_edge_count[v] < 3)
          continue;

        const Point_3& point = mesh_.point(v);
        const Vector_3 normal = normalize_vec<Vector_3, FT>(normal_sums[v]);

        saved_initial_candidates_.emplace_back(point, normal);
      }

      return feature_segments;
    }


    std::vector<std::pair<Point_3, Vector_3>> find_initial_feature_vertex_candidates(FT minimum_angle) override {
      std::vector<std::pair<Point_3, Vector_3>> initial_candidates = saved_initial_candidates_;

      saved_initial_candidates_.clear();

      return initial_candidates;
    }


    void write_to_obj(const char * filename) const {

      typename CGAL::Surface_mesh<Point_3> mesh;

      for (const Triangle& tri : triangles_) {
        const Point_3& p0 = tri.a;
        const Point_3& p1 = tri.b;
        const Point_3& p2 = tri.c;

        typename CGAL::Surface_mesh<Point_3>::Vertex_index v0 = mesh.add_vertex(p0);
        typename CGAL::Surface_mesh<Point_3>::Vertex_index v1 = mesh.add_vertex(p1);
        typename CGAL::Surface_mesh<Point_3>::Vertex_index v2 = mesh.add_vertex(p2);

        mesh.add_face(v0, v1, v2);
      }

      CGAL::IO::write_polygon_mesh(filename, mesh);

      std::cout << "Splat mesh written to " << filename << std::endl;
    }

  private:


    /**
     * @brief Filters triangles within each grid cell based on their normals.
     */
    void filter_triangles_by_normals() {
      constexpr double min_normal_sum_length = 0.1;
      constexpr double min_normal_sum_length_squared =
        min_normal_sum_length * min_normal_sum_length;

      for (Cell& cell : cells_) {
        auto& ids = cell.triangle_ids_within_distance_d;

        if (ids.empty()) {
          continue;
        }

        // Step 1: Sum the associated triangle normals without normalization.
        Vector_3 normal_sum = CGAL::NULL_VECTOR;

        for (Index tid : ids) {
          normal_sum = normal_sum + triangle_normal(tid);
        }

        // Step 2: If the sum is too small, keep all associated triangles.
        const double length_squared =
          CGAL::to_double(normal_sum.squared_length());

        if (length_squared < min_normal_sum_length_squared) {
          continue;
        }

        // Step 3: Remove triangles whose normals oppose the summed normal.
        ids.erase(
          std::remove_if(
            ids.begin(),
            ids.end(),
            [&](Index tid) {
              return CGAL::to_double(triangle_normal(tid) * normal_sum) < 0.0;
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

        // Use the distance-associated triangle list, not cell membership.
        const auto& triangle_ids = cell.triangle_ids_within_distance_d;

        if (triangle_ids.empty()) {
          continue;
        }

        // Find the sampled direction with the largest minimum dot product.
        double best_score = -std::numeric_limits<double>::infinity();
        Vector_3 best_normal = CGAL::NULL_VECTOR;

        for (const Vector_3& candidate : sphere_samples) {
          double min_dot = std::numeric_limits<double>::infinity();

          for (Index tid : triangle_ids) {
            const double dot = CGAL::to_double(
              triangle_normal(tid) * candidate);

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
     * @brief Clears all per-cell accumulators.
     */
    void clear() {
      for (Cell& c : cells_) {
        c.triangle_ids_contained.clear();
        c.triangle_ids_within_distance_d.clear();
      }

      // Make sure all cells are cleared.
      for (const Cell& c : cells_) {
        CGAL_assertion(c.triangle_ids_contained.empty());
        CGAL_assertion(c.triangle_ids_within_distance_d.empty());
      }
    }

    /**
     * @brief Inserts a single triangle into the grid.
     *
     * @param triangle_id Index of the triangle in the input mesh.
     */
    void insert(Index triangle_id) {

      const Triangle& tri = triangle(triangle_id);

      const Point_3& a = tri.a;
      const Point_3& b = tri.b;
      const Point_3& c = tri.c;

      Point_3 triangle_min = Point_3(
        (std::min)((std::min)(a.x(), b.x()), c.x()),
        (std::min)((std::min)(a.y(), b.y()), c.y()),
        (std::min)((std::min)(a.z(), b.z()), c.z())
      );
      Point_3 triangle_max = Point_3(
        (std::max)((std::max)(a.x(), b.x()), c.x()),
        (std::max)((std::max)(a.y(), b.y()), c.y()),
        (std::max)((std::max)(a.z(), b.z()), c.z())
      );

      int min_ix, min_iy, min_iz;
      CGAL_assertion(grid_.to_grid_coordinates(triangle_min, min_ix, min_iy, min_iz));

      int max_ix, max_iy, max_iz;
      CGAL_assertion(grid_.to_grid_coordinates(triangle_max, max_ix, max_iy, max_iz));

      // And add it to all neighboring cells within distance d.
      for (int nx = min_ix; nx <= max_ix; ++nx) {
        for (int ny = min_iy; ny <= max_iy; ++ny) {
          for (int nz = min_iz; nz <= max_iz; ++nz) {
            if (!grid_.valid_coordinates(nx, ny, nz)) {
              continue;
            }

            // Check if the triangle actually intersects this cell's bounding box.  
            if (triangle_intersects_aabb<Point_3, Vector_3, FT>(
              a, b, c, grid_.cell_min_corner(nx, ny, nz), grid_.cell_max_corner(nx, ny, nz)))
            {
              const std::size_t nidx = grid_.flat_index(nx, ny, nz);
              cells_[nidx].add_triangle_contained(triangle_id);

              // It is then also within distance d of all cells around this cell that intersect the triangle's bounding box.
              // Add it to all neighboring cells around this cell.
              for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                  for (int dz = -1; dz <= 1; ++dz) {
                    const int nnx = nx + dx;
                    const int nny = ny + dy;
                    const int nnz = nz + dz;

                    if (!grid_.valid_coordinates(nnx, nny, nnz)) {
                      continue;
                    }

                    const std::size_t nnidx = grid_.flat_index(nnx, nny, nnz);
                    cells_[nnidx].add_triangle_within_distance_d(triangle_id);
                  }
                }
              }
            }
          }
        }
      }
    }



    /**
     * @brief Finds input triangles in grid cells intersecting a query neighborhood.
     *
     * @param p Query position.
     * @param radius Search radius.
     *
     * @return Unique input triangle identifiers in the neighborhood.
    */
    std::vector<Index> nearby_triangle_ids(const Point_3& p, FT radius) const {
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

            ids.insert(ids.end(), c->triangle_ids_contained.begin(), c->triangle_ids_contained.end());
          }
        }
      }

      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
      return ids;
    }

    /**
     * @brief Returns the triangle IDs of triangles within distance d of the specified cell's bounding box.
     * @param ix The x-coordinate of the cell in the grid.
     * @param iy The y-coordinate of the cell in the grid.
     * @param iz The z-coordinate of the cell in the grid.
     * @return A reference to the vector of triangle IDs within distance d of the cell's bounding box.
     */
    const std::vector<Index>& triangle_ids_within_distance_d_of_box(std::size_t ix, std::size_t iy, std::size_t iz) const {
      const Cell* c = cell(ix, iy, iz);
      CGAL_assertion(c != nullptr);

      return c->triangle_ids_within_distance_d;
    }

    /**
     * @brief Returns the point IDs of triangles contained within the specified cell's bounding box.
     * @param ix The x-coordinate of the cell in the grid.
     * @param iy The y-coordinate of the cell in the grid.
     * @param iz The z-coordinate of the cell in the grid.
     * @return A reference to the vector of triangle IDs contained within the cell's bounding box.
     */
    const std::vector<Index>& triangle_ids_in_box(std::size_t ix, std::size_t iy, std::size_t iz) const {
      const Cell* c = cell(ix, iy, iz);
      CGAL_assertion(c != nullptr);

      return c->triangle_ids_contained;
    }

    /**
     * @brief A single grid cell.
     *
     * Each cell stores the indices of the triangles that fall inside it.
     */
    struct Cell
    {
      // The indices of all input triangles within distance d of the cell's bounding box.
      std::vector<Index> triangle_ids_within_distance_d;
      std::vector<Index> triangle_ids_contained;

      void add_triangle_contained(Index tid) {
        bool exists_already = std::find(triangle_ids_contained.begin(), triangle_ids_contained.end(), tid) != triangle_ids_contained.end();

        if (exists_already) {
          return;
        }

        triangle_ids_contained.push_back(tid);
      }
      void add_triangle_within_distance_d(Index tid) {
        bool exists_already = std::find(triangle_ids_within_distance_d.begin(), triangle_ids_within_distance_d.end(), tid) != triangle_ids_within_distance_d.end();

        if (exists_already) {
          return;
        }

        triangle_ids_within_distance_d.push_back(tid);
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

    const Triangle& triangle(Index triangle_id) const {
      return triangles_[triangle_id];
    }

    const Vector_3& triangle_normal(Index triangle_id) const {
      return triangles_[triangle_id].normal;
    }

    const Mesh& mesh_;
    std::vector<Triangle> triangles_;
    const Box_grid<Kernel_>& grid_;
    std::vector<Cell> cells_;
    std::vector<Vector_3> block_normals_;
    std::vector<std::pair<Point_3, Vector_3>> saved_initial_candidates_;
};

} // namespace Splat_surface_reconstruction_3
} // namespace internal
} // namespace CGAL

#endif