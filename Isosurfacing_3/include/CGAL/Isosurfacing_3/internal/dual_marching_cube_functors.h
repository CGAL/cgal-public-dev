#ifndef CGAL_ISOSURFACING_3_INTERNAL_MANIFOLD_DUAL_CONTOURING_H
#define CGAL_ISOSURFACING_3_INTERNAL_MANIFOLD_DUAL_CONTOURING_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/tables.h>
#include <CGAL/Isosurfacing_3/internal/marching_cubes_functors.h>

#include <CGAL/assertions.h>
#include <CGAL/Container_helper.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/enumerable_thread_specific.h>
#else
#include <vector>
#endif

#include <array>
#include <bitset>
#include <unordered_map>

#include <CGAL/Eigen_matrix.h>
#include <CGAL/Eigen_vector.h>
#include <Eigen/SVD>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

// compute the dual vertices of each patch of a cell
template <typename Domain>
bool patch_position_QEM(
    const Domain& domain,
    const typename Domain::cell_descriptor& c,
    const std::vector<typename Domain::Geom_traits::Point_3>& patch_points,
    const std::vector<typename Domain::Geom_traits::Vector_3>& patch_normals,
    const bool constrain_to_cell,
    typename Domain::Geom_traits::Point_3& p) 
{
    using Geom_traits = typename Domain::Geom_traits;
    using FT = typename Geom_traits::FT;
    using Point_3 = typename Geom_traits::Point_3;
    using Vector_3 = typename Geom_traits::Vector_3;

    using Eigen_vector_3 = Eigen_vector<FT, 3>;
    using Eigen_matrix_3 = Eigen_matrix<FT, 3, 3>;
    using Eigen_vector_x = Eigen_vector<FT>;
    using Eigen_matrix_x = Eigen_matrix<FT>;

    typename Geom_traits::Compute_x_3 x_coord = domain.geom_traits().compute_x_3_object();
    typename Geom_traits::Compute_y_3 y_coord = domain.geom_traits().compute_y_3_object();
    typename Geom_traits::Compute_z_3 z_coord = domain.geom_traits().compute_z_3_object();
    typename Geom_traits::Construct_point_3 point = domain.geom_traits().construct_point_3_object();

    //
    FT x_min, y_min, z_min, x_max, y_max, z_max;
    x_min = y_min = z_min =   (std::numeric_limits<double>::max)();
    x_max = y_max = z_max = - (std::numeric_limits<double>::max)();
    FT x(0), y(0), z(0);

    if(constrain_to_cell)
    {
        typename Domain::Cell_vertices vertices = domain.cell_vertices(c);
        for(const auto& v : vertices)
        {
            const Point_3& cp = domain.point(v);
            x_min = std::min(x_min, x_coord(cp));
            y_min = std::min(y_min, y_coord(cp));
            z_min = std::min(z_min, z_coord(cp));

            x_max = std::max(x_max, x_coord(cp));
            y_max = std::max(y_max, y_coord(cp));
            z_max = std::max(z_max, z_coord(cp));
        }
    }

    std::size_t en = patch_points.size();
    if (en < 3) return false; // at least 3 points needed for a plane

    for(const auto& ep : patch_points)
    {
        x += x_coord(ep);
        y += y_coord(ep);
        z += z_coord(ep);
    }

    Point_3 com = point(x / FT(en), y / FT(en), z / FT(en));

    // SVD QEM 
    Eigen_matrix_3 A;
    A.setZero();
    Eigen_vector_3 rhs;
    rhs.setZero();
    for(std::size_t i=0; i<patch_points.size(); ++i)
    {
        Eigen::Matrix<FT, 3, 1> n_k;
        n_k(0) = x_coord(patch_normals[i]);
        n_k(1) = y_coord(patch_normals[i]);
        n_k(2) = z_coord(patch_normals[i]);

        Eigen::Matrix<FT, 3, 1> p_k;
        p_k(0) = x_coord(patch_points[i]);
        p_k(1) = y_coord(patch_points[i]);
        p_k(2) = z_coord(patch_points[i]);

        FT d_k = n_k.transpose() * p_k;
        A += n_k * n_k.transpose();
        rhs += d_k * n_k;
    }

    Eigen::JacobiSVD<Eigen::Matrix<FT, 3, 3>> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    svd.setThreshold(1e-3);

    Eigen::Matrix<FT, 3, 1> x_hat;
    x_hat << x_coord(com), y_coord(com), z_coord(com);

    Eigen::Matrix<FT, 3, 1> v_svd = x_hat + svd.solve(rhs - A * x_hat);

    if(constrain_to_cell)
    {
        v_svd[0] = std::clamp(v_svd[0], x_min, x_max);
        v_svd[1] = std::clamp(v_svd[1], y_min, y_max);
        v_svd[2] = std::clamp(v_svd[2], z_min, z_max);
    }

    p = point(v_svd[0], v_svd[1], v_svd[2]);
    return true;
}




// generate quads for dual vertices in DMC
template <typename Domain,
          typename EdgeDescriptor, 
          typename PolygonRange>
void generate_quad_dmc(const EdgeDescriptor& e,
                       const Domain& domain,
                       const std::unordered_map<EdgeDescriptor, std::vector<std::size_t>>& edge_to_dual_vertices,
                       PolygonRange& polygons)
{
    using FT = typename Domain::Geom_traits::FT;

    const auto& vertices = domain.incident_vertices(e);

    const FT val_0 = domain.value(vertices[0]);
    const FT val_1 = domain.value(vertices[1]);
    auto it = edge_to_dual_vertices.find(e);

    if (it == edge_to_dual_vertices.end())
        return; 

    const std::vector<std::size_t>& dual_vertices_idx = it->second;
    if (dual_vertices_idx.size()< 3)
        return; 

    // consistent global winding direction
    std::vector<std::size_t> temp = dual_vertices_idx; // make a local copy
    if(val_0>val_1) std::reverse(temp.begin(), temp.end());

    polygons.emplace_back(temp.begin(), temp.end());
}


// reorder the dual vertices of each edge to be in local cyclic order
template <typename Domain, 
          typename EdgeDescriptor>
void reorder_edge_to_dual_vertices(const Domain& domain,
                                   const std::unordered_map<typename Domain::cell_descriptor, std::vector<size_t>>& cell_to_dual_vertices,
                                   std::unordered_map<EdgeDescriptor, std::vector<std::size_t>>& edge_to_dual_vertices)
{
    for (auto& pair : edge_to_dual_vertices)
    {
        const EdgeDescriptor& e = pair.first;
        std::vector<std::size_t>& unordered_dual_vertices = pair.second;
        std::vector<std::size_t> ordered_dual_vertices;
        std::unordered_set<std::size_t> seen; // track unique duals for this edge

        const auto& cells = domain.incident_cells(e);
        for (const auto& cell : cells)
        {
            auto it = cell_to_dual_vertices.find(cell);
            if (it != cell_to_dual_vertices.end())
            {
                for (const auto& dual_idx : it->second)
                {
                    // add this dual if it actually appears in unordered_dual_vertices for this edge
                    if (std::find(unordered_dual_vertices.begin(), unordered_dual_vertices.end(), dual_idx) != unordered_dual_vertices.end()
                        && seen.insert(dual_idx).second) 
                    {
                        ordered_dual_vertices.push_back(dual_idx);
                    }
                }
            }
        }
        // replace the old by the ordered
        unordered_dual_vertices = std::move(ordered_dual_vertices);
    }
}


// dual marching cube
template <typename Domain, 
          typename PolygonRange>
void dual_marching_cubes(const Domain& domain,
                         const typename Domain::Geom_traits::FT isovalue,
                         std::vector<typename Domain::Geom_traits::Point_3>& points,
                         PolygonRange& polygons,
                         bool constrain_to_cell)
{
    using FT = typename Domain::Geom_traits::FT;
    using Point_3 = typename Domain::Geom_traits::Point_3;
    using Vector_3 = typename Domain::Geom_traits::Vector_3;
    using vertex_descriptor = typename Domain::vertex_descriptor;
    using edge_descriptor = typename Domain::edge_descriptor;
    using cell_descriptor = typename Domain::cell_descriptor;

    // edghe intersections and gradients 
    std::unordered_map<edge_descriptor, std::size_t> edge_to_point_id;
    std::vector<Point_3> edge_points;
    std::vector<Vector_3> edge_gradients;

    static_assert(Domain::VERTICES_PER_CELL == 8);

    auto edge_positioner = [&](const edge_descriptor& e) 
    {
        const auto& evs = domain.incident_vertices(e);
        const vertex_descriptor& v0 = evs[0];
        const vertex_descriptor& v1 = evs[1];
        const Point_3& p0 = domain.point(v0);
        const Point_3& p1 = domain.point(v1);
        const FT val0 = domain.value(v0);
        const FT val1 = domain.value(v1);

        Point_3 p;
        if (!domain.construct_intersection(p0, p1, val0, val1, isovalue, p))
            return;

        Vector_3 g = domain.gradient(p);

        edge_to_point_id[e] = edge_points.size();
        edge_points.push_back(p);
        edge_gradients.push_back(g);
    };

    domain.for_each_edge(edge_positioner);

    // compute the dual vertices for each polygon patch and output a edge to dual vertex map for the quad generation later
    std::unordered_map<edge_descriptor, std::vector<size_t>> edge_to_dual_vertices;
    std::unordered_map<Point_3, std::size_t> point_to_index; // to avoid duplicates
    std::unordered_map<cell_descriptor, std::vector<size_t>> cell_to_dual_vertices; // to keep track of local cyclic order

    auto cell_patch_positioner = [&](const cell_descriptor& c)
    {
        constexpr std::size_t vpc = Domain::VERTICES_PER_CELL;

        std::array<FT, vpc> values;
        std::array<Point_3, vpc> corners;
        const std::size_t i_case = get_cell_corners(domain, c, isovalue, corners, values, false);
        const auto& cell_edges = domain.cell_edges(c);

        const int* row = &Cube_table::polygon_cases[i_case * 16];
        int k = 0;
        std::vector<std::size_t> patch_dual_vertices_idx;
        std::vector<std::vector<int>> all_patch;

        // loop over edges in the lookup table
        while (row[k] != -1)
        {
            std::vector<int> current_patch_edge_local_idx;
            while (row[k] != -2 && row[k] != -1)
                current_patch_edge_local_idx.push_back(row[k++]); // collect local edge indices in one patch
            

            // collect intersections and gradients for this patch
            std::vector<Point_3> patch_points;
            std::vector<Vector_3> patch_grads;
            for (int local_edge_idx : current_patch_edge_local_idx) 
            {
                const edge_descriptor& global_edge = cell_edges[local_edge_idx]; 
                const auto it = edge_to_point_id.find(global_edge);
                if (it == edge_to_point_id.end())
                    continue; // skip if  not intersected
                auto global_idx = it->second;
                patch_points.push_back(edge_points[global_idx]);
                patch_grads.push_back(edge_gradients[global_idx]);
            }

            // compute the dual vertex for this patch
            if (patch_points.size() >= 3) {
                Point_3 dual_vertex;
                bool success = patch_position_QEM(domain, c, patch_points, patch_grads, constrain_to_cell, dual_vertex);
                if(success)
                {
                    std::size_t dual_vertex_idx;
                    auto u = point_to_index.find(dual_vertex);
                    if (u == point_to_index.end()) 
                    {
                        dual_vertex_idx = points.size();
                        points.push_back(dual_vertex);
                        point_to_index[dual_vertex] = dual_vertex_idx;
                    } 
                    else  dual_vertex_idx = u->second;

                    patch_dual_vertices_idx.push_back(dual_vertex_idx); // collect a dual vertex 
                    all_patch.push_back(current_patch_edge_local_idx); // collect a patch
                }
            }

            if (row[k] == -2) ++k;
        }

        CGAL_assertion(patch_dual_vertices_idx.size() == all_patch.size());

        // map edge to dual vertices
        // the ordering of patch_dual_vertices is consistent with all_patch
        int patch_size = patch_dual_vertices_idx.size();
        for (int patch_idx = 0; patch_idx < patch_size; patch_idx++)
        {
            // loop over all edges in this patch
            for (int local_edge_idx : all_patch[patch_idx]) 
            {
                auto global_edge = cell_edges[local_edge_idx];
                edge_to_dual_vertices[global_edge].push_back(patch_dual_vertices_idx[patch_idx]);
                cell_to_dual_vertices[c].push_back(patch_dual_vertices_idx[patch_idx]); 
            }
        }
    };

    domain.for_each_cell(cell_patch_positioner);

    // reorder the dual vertices of each edge to be in local cyclic order
    reorder_edge_to_dual_vertices(domain, cell_to_dual_vertices, edge_to_dual_vertices);

    // generate quads
    auto generate_quad = [&](const edge_descriptor& e) 
    {
        generate_quad_dmc(e, domain, edge_to_dual_vertices, polygons);
    };
    domain.for_each_edge(generate_quad);

};


} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif