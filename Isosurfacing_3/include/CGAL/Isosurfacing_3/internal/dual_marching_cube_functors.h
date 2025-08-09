#ifndef CGAL_ISOSURFACING_3_INTERNAL_DUAL_MARCHING_CUBE_H
#define CGAL_ISOSURFACING_3_INTERNAL_DUAL_MARCHING_CUBE_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/tables.h>
#include <CGAL/Isosurfacing_3/internal/tmc_dmc_functors.h>
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

    // vertex placement strategy
    enum class Vertex_strategy
    {
        QEM,
        Centroid
    };

namespace internal {

// generate dual vertex of each patch as the centroid of the patch
template <typename Domain>
bool patch_position_centroid(const Domain& domain,
                             const std::vector<typename Domain::Geom_traits::Point_3>& patch_points,
                             typename Domain::Geom_traits::Point_3& p)
{
    using Geom_traits = typename Domain::Geom_traits;
    using FT          = typename Geom_traits::FT;
    using Point_3     = typename Geom_traits::Point_3;

    typename Geom_traits::Compute_x_3 x_coord = domain.geom_traits().compute_x_3_object();
    typename Geom_traits::Compute_y_3 y_coord = domain.geom_traits().compute_y_3_object();
    typename Geom_traits::Compute_z_3 z_coord = domain.geom_traits().compute_z_3_object();
    typename Geom_traits::Construct_point_3 point = domain.geom_traits().construct_point_3_object();

    const std::size_t n = patch_points.size();
    if (n == 0)
        return false;

    FT x(0), y(0), z(0);
    for (const auto& q : patch_points)
    {
        x += x_coord(q);
        y += y_coord(q);
        z += z_coord(q);
    }

    p = point(x / n, y / n, z / n);
    return true;
}


// compute the dual vertex of each patch of a cell using qem
template <typename Domain>
bool patch_position_QEM(const Domain& domain,
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


// check if a point is nan
template <typename Domain>
bool is_valid_point(const typename Domain::Geom_traits::Point_3& p)
{
    return CGAL::is_finite(p.x()) &&
           CGAL::is_finite(p.y()) &&
           CGAL::is_finite(p.z());
}


// check duplicate patches
template <typename Point_3>
bool is_duplicate_patch(const std::vector<Point_3>& new_patch,
                        const std::vector<std::vector<Point_3>>& existing_patches,
                        double eps = 1e-8)
{
    auto are_points_equal = [eps](const Point_3& a, const Point_3& b) {
        return CGAL::squared_distance(a, b) < eps * eps;
    };

    auto point_cmp = [](const Point_3& a, const Point_3& b) {
        if (a.x() != b.x()) return a.x() < b.x();
        if (a.y() != b.y()) return a.y() < b.y();
        return a.z() < b.z();
    };

    std::vector<Point_3> sorted_new = new_patch;
    std::sort(sorted_new.begin(), sorted_new.end(), point_cmp);

    for (const auto& patch : existing_patches)
    {
        if (patch.size() != new_patch.size()) continue;

        std::vector<Point_3> sorted_existing = patch;
        std::sort(sorted_existing.begin(), sorted_existing.end(), point_cmp);

        bool equal = true;
        for (std::size_t i = 0; i < sorted_new.size(); ++i)
        {
            if (!are_points_equal(sorted_new[i], sorted_existing[i]))
            {
                equal = false;
                break;
            }
        }

        if (equal) return true;
    }

    return false;
}


// dual marching cube
template <typename Domain, 
          typename PolygonRange>
void dual_marching_cubes(const Domain& domain,
                         const typename Domain::Geom_traits::FT isovalue,
                         std::vector<typename Domain::Geom_traits::Point_3>& points,
                         PolygonRange& polygons,
                         bool constrain_to_cell,
                         Isosurfacing::Vertex_strategy strategy)
{
    using FT = typename Domain::Geom_traits::FT;
    using Point_3 = typename Domain::Geom_traits::Point_3;
    using Vector_3 = typename Domain::Geom_traits::Vector_3;
    using vertex_descriptor = typename Domain::vertex_descriptor;
    using edge_descriptor = typename Domain::edge_descriptor;
    using cell_descriptor = typename Domain::cell_descriptor;

    static_assert(Domain::VERTICES_PER_CELL == 8);

    // edghe intersections and gradients 
    // std::unordered_map<edge_descriptor, std::size_t> edge_to_point_id;
    // std::vector<Point_3> edge_points;
    // std::vector<Vector_3> edge_gradients;
    // auto edge_positioner = [&](const edge_descriptor& e) 
    // {
    //     const auto& evs = domain.incident_vertices(e);
    //     const vertex_descriptor& v0 = evs[0];
    //     const vertex_descriptor& v1 = evs[1];
    //     const Point_3& p0 = domain.point(v0);
    //     const Point_3& p1 = domain.point(v1);
    //     const FT val0 = domain.value(v0);
    //     const FT val1 = domain.value(v1);
    //     Point_3 p;
    //     if (!domain.construct_intersection(p0, p1, val0, val1, isovalue, p))
    //         return;
    //     Vector_3 g = domain.gradient(p);
    //     edge_to_point_id[e] = edge_points.size();
    //     edge_points.push_back(p);
    //     edge_gradients.push_back(g);
    // };
    // domain.for_each_edge(edge_positioner);

    // compute the dual vertices for each polygon patch and output a edge to dual vertex map for the quad generation later
    std::unordered_map<edge_descriptor, std::vector<size_t>> edge_to_dual_vertices;
    std::unordered_map<Point_3, std::size_t> point_to_index; // to avoid duplicates
    std::unordered_map<cell_descriptor, std::vector<size_t>> cell_to_dual_vertices; // to keep track of local cyclic order

    // debug
    // std::ofstream debug_out("patch_dual_vertices_log.txt");  
    // constexpr FT tol = 1e-5;
    // const Point_3 target_vertex(-0.295698, 2.88988, 3.21928);
    // const Point_3 target_vertex(-0.285791, 2.90516, 3.30311);


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

                // const auto it = edge_to_point_id.find(global_edge);
                // if (it == edge_to_point_id.end())
                //     continue; // skip if  not intersected
                // auto global_idx = it->second;
                // patch_points.push_back(edge_points[global_idx]);
                // patch_grads.push_back(edge_gradients[global_idx]);

                const auto& evs = domain.incident_vertices(global_edge);
                const vertex_descriptor& v0 = evs[0];
                const vertex_descriptor& v1 = evs[1];
                const Point_3& p0 = domain.point(v0);
                const Point_3& p1 = domain.point(v1);
                const FT val0 = domain.value(v0);
                const FT val1 = domain.value(v1);

                Point_3 p;
                if (!domain.construct_intersection(p0, p1, val0, val1, isovalue, p)) continue;

                // if (!is_valid_point<Domain>(p)) continue; // skip invalid points
                if (!is_valid_point<Domain>(p))
                {
                    std::cerr << "Invalid point detected: p0 (" 
                                << CGAL::to_double(p0.x()) << ", "
                                << CGAL::to_double(p0.y()) << ", "
                                << CGAL::to_double(p0.z()) << ")\n"
                                << "p1 ("
                                << CGAL::to_double(p1.x()) << ", "
                                << CGAL::to_double(p1.y()) << ", "
                                << CGAL::to_double(p1.z()) << ")\n";
                    continue;
                }


                Vector_3 g = domain.gradient(p);
                patch_points.push_back(p);
                patch_grads.push_back(g);
            }


            // compute the dual vertex for this patch
            if (patch_points.size() >= 3) 
            {
                Point_3 dual_vertex;

                bool success = false;
                switch(strategy)
                {
                    case Vertex_strategy::QEM:
                        success = patch_position_QEM(domain, c, patch_points, patch_grads, constrain_to_cell, dual_vertex);
                        break;
                    
                    case Vertex_strategy::Centroid:
                        success = patch_position_centroid(domain, patch_points, dual_vertex);
                        break;
                }
                
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

                    // debug
                    // auto close = [](const Point_3& a, const Point_3& b) 
                    // {
                    //     return (CGAL::abs(a.x() - b.x()) < tol &&
                    //             CGAL::abs(a.y() - b.y()) < tol &&
                    //             CGAL::abs(a.z() - b.z()) < tol);
                    // };
                    // if (close(dual_vertex, target_vertex)) 
                    // {
                    //     std::cout << "Matched target dual vertex!\n";
                    //     std::cout << "i_case = " << i_case << "\n";
                    //     std::ofstream ofs("matched_cell_vertices.off");
                    //     ofs << "OFF\n";
                    //     ofs << "8 12 0\n"; 
                    //     for (std::size_t i = 0; i < vpc; ++i)
                    //     {
                    //         const Point_3& p = corners[i];
                    //         ofs << p.x() << " " << p.y() << " " << p.z() << "\n";
                    //     }
                    //     ofs.close();
                    //     std::exit(0);
                    // }
                    // for (const auto& pt : patch_points)
                    //     debug_out << "(" << pt.x() << " " << pt.y() << " " << pt.z() << ") ";
                    // debug_out << "=> (" << dual_vertex.x() << " " << dual_vertex.y() << " " << dual_vertex.z() << ")\n";
                    // if (close(dual_vertex, target_vertex)) {
                    //     debug_out << "Target dual vertex generated, stopping.\n";
                    //     debug_out.close();
                    //     std::exit(0);
                    // }

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

}


// check if a point p is on the edge e
template <typename Domain>
bool point_lies_on_segment(const Domain& domain,
                           const typename Domain::Geom_traits::Point_3& p, 
                           const typename Domain::edge_descriptor& e,
                           double tol = 1e-6)
{
    using Point_3 = typename Domain::Geom_traits::Point_3;
    using vertex_descriptor = typename Domain::vertex_descriptor;
    using Vector_3 = typename Domain::Geom_traits::Vector_3;

    const auto& evs = domain.incident_vertices(e);
    const vertex_descriptor& v0 = evs[0];
    const vertex_descriptor& v1 = evs[1];
    const Point_3& p0 = domain.point(v0);
    const Point_3& p1 = domain.point(v1);

    Vector_3 v = p1 - p0;
    Vector_3 w = p - p0;

    Vector_3 cross = CGAL::cross_product(v, w);
    if (CGAL::squared_length(cross) > tol * tol)
        return false;

    // check 0 <= dot(w, v) <= ||v||^2
    auto dot = w * v;
    auto v_len2 = CGAL::squared_length(v);
    if (dot < -tol || dot > v_len2 + tol)
        return false;

    return true;
}


// topologically correct dual marching cube
// dual marching cube
template <typename Domain, 
          typename PolygonRange>
void dual_marching_cubes_tmc(const Domain& domain,
                         const typename Domain::Geom_traits::FT isovalue,
                         std::vector<typename Domain::Geom_traits::Point_3>& points,
                         PolygonRange& polygons,
                         bool constrain_to_cell,
                         Isosurfacing::Vertex_strategy strategy)
{
    using FT = typename Domain::Geom_traits::FT;
    using Point_3 = typename Domain::Geom_traits::Point_3;
    using Vector_3 = typename Domain::Geom_traits::Vector_3;
    using vertex_descriptor = typename Domain::vertex_descriptor;
    using edge_descriptor = typename Domain::edge_descriptor;
    using cell_descriptor = typename Domain::cell_descriptor;

    static_assert(Domain::VERTICES_PER_CELL == 8);

    // compute the dual vertices for each polygon patch and output a edge to dual vertex map for the quad generation later
    std::unordered_map<edge_descriptor, std::vector<size_t>> edge_to_dual_vertices;
    std::unordered_map<Point_3, std::size_t> point_to_index; // to avoid duplicates
    std::unordered_map<cell_descriptor, std::vector<size_t>> cell_to_dual_vertices; // to keep track of local cyclic order

    // document the source of each dual vertex for nonmanifold edge tracking
    enum class DualSrc : uint8_t { Regular, TMC_NoTunnel, TMC_Tunnel };

    // prefer Tunnel > NoTunnel > Regular when the same point is reused
    auto upgrade_src = [](DualSrc oldv, DualSrc newv) {
        if (oldv == DualSrc::TMC_Tunnel || newv == DualSrc::TMC_Tunnel) return DualSrc::TMC_Tunnel;
        if (oldv == DualSrc::TMC_NoTunnel || newv == DualSrc::TMC_NoTunnel) return DualSrc::TMC_NoTunnel;
        return DualSrc::Regular;
    };

    std::unordered_map<std::size_t, DualSrc> dual_src;

    // debug
    // std::array<Point_3, 8> target = {{
    //                                 Point_3(-2.10145, -2.10145, 2.10145),
    //                                 Point_3(-1.95652, -2.10145, 2.10145),
    //                                 Point_3(-2.10145, -1.95652, 2.10145),
    //                                 Point_3(-1.95652, -1.95652, 2.10145),
    //                                 Point_3(-2.10145, -2.10145, 2.24638),
    //                                 Point_3(-1.95652, -2.10145, 2.24638),
    //                                 Point_3(-2.10145, -1.95652, 2.24638),
    //                                 Point_3(-1.95652, -1.95652, 2.24638),
    //                             }};

    // // debug 
    // const Point_3 target(185.299, 179.63, 72.2626);
    // const Point_3 target(183.596, 178.42400000000001, 73.252499999999998);
    // const Point_3 target(64.6465, 134.465, 93.0505);
    // const Point_3 target(64.6465, 139.636, 93.0505);
    // const Point_3 target(175.83799999999999, 87.919200000000004, 92.060599999999994);
    // const Point_3 target(176.34899999999999, 88.960099999999997, 92.653599999999997);
    // const Point_3 target(82.830799999999996, 90.669200000000004, 61.5822);
    // const Point_3 target(82.747500000000002, 90.505099999999999, 61.373699999999999);
    // const Point_3 target(189.779, 154.18299999999999, 53.939999999999998);
    // const Point_3 target(190.131, 155.852, 54.178600000000003);
    // const Point_3 target(109.822, 58.737900000000003, 13.5657);
    // const Point_3 target(109.726, 58.836799999999997, 14.1892);
    // const Point_3 target(63.229100000000003, 127.131, 88.099199999999996);
    // const Point_3 target(61.230699999999999, 140.911, 84.399600000000007);
    // const Point_3 target(169.55600000000001, 95.151700000000005, 49.012700000000002);
    // const Point_3 target(166.16800000000001, 93.595299999999995, 46.805700000000002);
    // const Point_3 target(109.907, 108.431, 35.1872);
    // const Point_3 target(109.661, 110.134, 34.3311);
    // const Point_3 target(110.302, 109.476, 35.0793);
    // const Point_3 target(110.331, 112.404, 34.3367);

    auto is_cell_tmc_tunnel = [&](const auto& triangles, const auto& cell_edges) -> bool
    {
        // Any triangle vertex that is NOT on any of the 12 cell edges ⇒ interior (hexagon) vertex.
        int interior_count = 0;

        for (const auto& tri : triangles) 
        {
            for (const auto& p : tri) 
            {
                bool on_any_edge = false;
                for (int i = 0; i < 12; ++i) 
                {
                    if (point_lies_on_segment(domain, p, cell_edges[i])) { on_any_edge = true; break; }
                }
                if (!on_any_edge) { interior_count++; }
            }
        }
        // Tunnel cases (popcount(q_sol)==6) create multiple interior hexagon uses.
        // One interior vertex might appear in many tris, so we just need a threshold > 0.
        return interior_count > 0;
    };

    auto cell_patch_positioner = [&](const cell_descriptor& c)
    {

        constexpr std::size_t vpc = Domain::VERTICES_PER_CELL;

        const auto& cell_edges = domain.cell_edges(c);

        std::array<FT, vpc> values;
        std::array<Point_3, vpc> corners;
        const std::size_t i_case = get_cell_corners(domain, c, isovalue, corners, values, true);

        // skip empty / full cells
        constexpr std::size_t ones = (1 << 8) - 1;
        if((i_case & ones) == ones || // all bits set
            (i_case & ones) == 0) // no bits set
            return;
        
        // // debug
        // for (int i = 0; i < 8; ++i) 
        // {
        //     if (CGAL::squared_distance(corners[i], target[i]) < 1e-8)
        //         return;
        // }

        // ambiguous case 
        const int tcm = Cube_table::t_ambig[i_case];
        if(tcm == 105)
        {
            std::vector<std::array<Point_3, 3>> triangles;
            std::vector<std::vector<Point_3>> all_patch_ambi;

            bool success = CGAL::Isosurfacing::internal::p_slice(domain, c, isovalue, corners, values, i_case, true, triangles);
            if(!success) return;
            
            const bool tmc_tunnel = is_cell_tmc_tunnel(triangles, cell_edges); // check if a tunnel is generated in the tmc case

            // then merge triangles for patches
            for (const auto& tri : triangles) 
            {
                bool merged = false;

                for (auto& patch : all_patch_ambi) 
                {
                    for (const auto& p1 : tri) 
                    {
                        if (std::find(patch.begin(), patch.end(), p1) != patch.end()) 
                        {
                            for (const auto& pt : tri) 
                            {
                                if (std::find(patch.begin(), patch.end(), pt) == patch.end()) patch.push_back(pt);
                            }
                            merged = true;
                            break;
                        }
                    }
                    if (merged) break;
                }

                if (!merged) all_patch_ambi.emplace_back(tri.begin(), tri.end()); // new patch with this triangle’s vertices
            }
            
            // // Deduplicate patches
            // std::vector<std::vector<Point_3>> unique_patches;
            // for (const auto& patch : all_patch_ambi)
            // {
            //     if (!is_duplicate_patch(patch, unique_patches))
            //         unique_patches.push_back(patch);
            // }
            // all_patch_ambi = std::move(unique_patches);

            // //debug
            // if ([&]() 
            //     {
            //         std::array<Point_3, 8> target = {{
            //             Point_3(-2.10145, -2.10145, 2.10145),
            //             Point_3(-1.95652, -2.10145, 2.10145),
            //             Point_3(-2.10145, -1.95652, 2.10145),
            //             Point_3(-1.95652, -1.95652, 2.10145),
            //             Point_3(-2.10145, -2.10145, 2.24638),
            //             Point_3(-1.95652, -2.10145, 2.24638),
            //             Point_3(-2.10145, -1.95652, 2.24638),
            //             Point_3(-1.95652, -1.95652, 2.24638),
            //         }};
            //         for (int i = 0; i < 8; ++i) {
            //             if (CGAL::squared_distance(corners[i], target[i]) > 1e-8)
            //                 return false;
            //         }
            //         return true;
            //     }())
            // {
            //     std::cout << "all_patch_ambi has patches of size ";
            //     for (const auto& patch: all_patch_ambi)
            //     {
            //         std::cout << patch.size() << " | "; 
            //         for (const auto& p: patch)
            //         {
            //             std::cout << "vertex: (" 
            //                     << CGAL::to_double(p.x()) << ", "
            //                     << CGAL::to_double(p.y()) << ", "
            //                     << CGAL::to_double(p.z()) << ") \n";
            //         }
            //     }
            //     std::cout<< "\n";
            // }

            // compute the dual vertex for each patch
            for (const auto& patch_points: all_patch_ambi )
            {
                if (patch_points.size() >= 3) 
                {
                    std::vector<Vector_3> patch_grads; 
                    std::vector<edge_descriptor> current_patch_edges;
                    std::vector<Point_3> filtered_patch_points; // filter out invalid points

                    for (const auto& p: patch_points)
                    {
                        // if (!is_valid_point<Domain>(p)) continue; // skip invalid points
                        if (!is_valid_point<Domain>(p))
                        {
                            std::cerr << "Invalid point detected: p (" 
                                      << CGAL::to_double(p.x()) << ", "
                                      << CGAL::to_double(p.y()) << ", "
                                      << CGAL::to_double(p.z()) << ")\n";
                            continue;
                        }

                        filtered_patch_points.push_back(p);
                        patch_grads.push_back(domain.gradient(p)); // collect gradients

                        // loop over each cell edge and check if it is one that the current patch point intersect with
                        for (int i = 0; i < 12; ++i) 
                        {
                            const auto& e = cell_edges[i];

                            // // debug
                            // if (d < 1e-4)
                            // {
                            //     std::cout << "Matching candidate:\n";
                            //     std::cout << "patch point: (" << CGAL::to_double(p.x()) << ", "
                            //                                 << CGAL::to_double(p.y()) << ", "
                            //                                 << CGAL::to_double(p.z()) << ")\n";
                            //     std::cout << "edge interp : (" << CGAL::to_double(pe.x()) << ", "
                            //                                 << CGAL::to_double(pe.y()) << ", "
                            //                                 << CGAL::to_double(pe.z()) << ")\n";
                            //     std::cout << "dist = " << d << "\n";
                            // }
                            if (point_lies_on_segment(domain, p,e))
                            {
                                current_patch_edges.push_back(e);
                                break;
                            }
                            
                        }
                    }

                    if (filtered_patch_points.size() < 3) continue; 
                    // //debug
                    // if ([&]() 
                    //     {
                    //         std::array<Point_3, 8> target = {{
                    //             Point_3(-2.10145, -2.10145, 2.10145),
                    //             Point_3(-1.95652, -2.10145, 2.10145),
                    //             Point_3(-2.10145, -1.95652, 2.10145),
                    //             Point_3(-1.95652, -1.95652, 2.10145),
                    //             Point_3(-2.10145, -2.10145, 2.24638),
                    //             Point_3(-1.95652, -2.10145, 2.24638),
                    //             Point_3(-2.10145, -1.95652, 2.24638),
                    //             Point_3(-1.95652, -1.95652, 2.24638),
                    //         }};
                    //         for (int i = 0; i < 8; ++i) {
                    //             if (CGAL::squared_distance(corners[i], target[i]) > 1e-8)
                    //                 return false;
                    //         }
                    //         return true;
                    //     }())
                    // {
                    //     std::cout << "current_patch_edges has patches of size "
                    //      << current_patch_edges.size(); 
                    //     std::cout<< "\n";
                    //     for (const auto& e: current_patch_edges)
                    //     {
                    //         const auto& evs = domain.incident_vertices(e);
                    //         const vertex_descriptor& p0 = evs[0];
                    //         const vertex_descriptor& p1 = evs[1];
                    //         const Point_3& v0 = domain.point(p0);
                    //         const Point_3& v1 = domain.point(p1);
                    //         std::cout << "Edge: (" 
                    //             << CGAL::to_double(v0.x()) << ", "
                    //             << CGAL::to_double(v0.y()) << ", "
                    //             << CGAL::to_double(v0.z()) << ") -- ("
                    //             << CGAL::to_double(v1.x()) << ", "
                    //             << CGAL::to_double(v1.y()) << ", "
                    //             << CGAL::to_double(v1.z()) << ") \n";
                    //     }
                    // }
            
                    Point_3 dual_vertex;

                    bool success = false;
                    switch(strategy)
                    {
                        case Vertex_strategy::QEM:
                            success = patch_position_QEM(domain, c, filtered_patch_points, patch_grads, constrain_to_cell, dual_vertex);
                            // if (success) std::cout << "using QEM \n";
                            // else std::cout << "QEM failed\n";

                            break;
                        
                        case Vertex_strategy::Centroid:
                            success = patch_position_centroid(domain, filtered_patch_points, dual_vertex);
                            // if (success) std::cout << "using centroid \n";
                            // else std::cout << "centroid failed\n";

                            break;
                    }
                    
                    // //debug
                    // if ([&]() 
                    //     {
                    //         std::array<Point_3, 8> target = {{
                    //             Point_3(-2.10145, -2.10145, 2.10145),
                    //             Point_3(-1.95652, -2.10145, 2.10145),
                    //             Point_3(-2.10145, -1.95652, 2.10145),
                    //             Point_3(-1.95652, -1.95652, 2.10145),
                    //             Point_3(-2.10145, -2.10145, 2.24638),
                    //             Point_3(-1.95652, -2.10145, 2.24638),
                    //             Point_3(-2.10145, -1.95652, 2.24638),
                    //             Point_3(-1.95652, -1.95652, 2.24638),
                    //         }};
                    //         for (int i = 0; i < 8; ++i) {
                    //             if (CGAL::squared_distance(corners[i], target[i]) > 1e-8)
                    //                 return false;
                    //         }
                    //         return true;
                    //     }())
                    // {
                    //     std::cout << "print patch edges (v0,v1) \n";
                    //     if (current_patch_edges.empty()) {
                    //         std::cout << "No edges found for current patch.\n";
                    //     }
                    //     for (const auto& e: current_patch_edges)
                    //     {
                    //         const auto& evs = domain.incident_vertices(e);
                    //         const vertex_descriptor& p0 = evs[0];
                    //         const vertex_descriptor& p1 = evs[1];
                    //         const Point_3& v0 = domain.point(p0);
                    //         const Point_3& v1 = domain.point(p1);
                    //         std::cout << "Edge: (" 
                    //             << CGAL::to_double(v0.x()) << ", "
                    //             << CGAL::to_double(v0.y()) << ", "
                    //             << CGAL::to_double(v0.z()) << ") -- ("
                    //             << CGAL::to_double(v1.x()) << ", "
                    //             << CGAL::to_double(v1.y()) << ", "
                    //             << CGAL::to_double(v1.z()) << ") \n";
                    //     }
                    // }

                    // debug
                    // if (CGAL::squared_distance(dual_vertex, target) < 1e-6) 
                    // {
                    //     std::cout << "generated from tmc cases \n";
                    //     std::cout << "Found matching dual vertex!\n";
                    //     std::cout << "dual vertex: (" 
                    //             << CGAL::to_double(dual_vertex.x()) << ", "
                    //             << CGAL::to_double(dual_vertex.y()) << ", "
                    //             << CGAL::to_double(dual_vertex.z()) << ")\n";
                    //     std::cout << "  Cell vertices:\n";
                    //     for (std::size_t i = 0; i < vpc; ++i) 
                    //     {
                    //         const Point_3& corner = domain.point(domain.cell_vertices(c)[i]);
                    //         std::cout << "    (" 
                    //                 << CGAL::to_double(corner.x()) << ", "
                    //                 << CGAL::to_double(corner.y()) << ", "
                    //                 << CGAL::to_double(corner.z()) << ")\n";
                    //     }
                    //     std::cout << "  Patch points:\n";
                    //     for (const auto& pt : patch_points) 
                    //     {
                    //         std::cout << "    (" 
                    //                 << CGAL::to_double(pt.x()) << ", "
                    //                 << CGAL::to_double(pt.y()) << ", "
                    //                 << CGAL::to_double(pt.z()) << ")\n";
                    //     }
                    //     std::cout << "  Triangles:\n";
                    //     for (std::size_t i = 0; i < triangles.size(); ++i)
                    //     {
                    //         const auto& tri = triangles[i];
                    //         std::cout << "  Triangle " << i << ":\n";
                    //         for (const auto& pt : tri)
                    //         {
                    //             std::cout << "    (" 
                    //                     << CGAL::to_double(pt.x()) << ", "
                    //                     << CGAL::to_double(pt.y()) << ", "
                    //                     << CGAL::to_double(pt.z()) << ")\n";
                    //         }
                    //     }
                    //     // std::exit(0);
                    // }

                    if(success)
                    {
                        std::size_t dual_vertex_idx;
                        auto u = point_to_index.find(dual_vertex);
                        if (u == point_to_index.end()) // if no duplicates
                        {
                            dual_vertex_idx = points.size();
                            points.push_back(dual_vertex);
                            point_to_index[dual_vertex] = dual_vertex_idx;

                            // label the source
                            dual_src[dual_vertex_idx] = tmc_tunnel ? DualSrc::TMC_Tunnel
                                               : DualSrc::TMC_NoTunnel;
                        } 
                        else  
                        {
                            dual_vertex_idx = u->second;

                            // if the point is duplicate, find the old tag and update the new one
                            auto it = dual_src.find(dual_vertex_idx);
                            DualSrc prev = (it == dual_src.end()) // if we cant find the previous label, default to the new tag
                                ? (tmc_tunnel ? DualSrc::TMC_Tunnel: DualSrc::TMC_NoTunnel) 
                                : it->second;
                            dual_src[dual_vertex_idx] = upgrade_src(prev, tmc_tunnel ? DualSrc::TMC_Tunnel
                                                                                    : DualSrc::TMC_NoTunnel);
                        }

                        cell_to_dual_vertices[c].push_back(dual_vertex_idx); 

                        for (const auto& e: current_patch_edges)
                        {
                            edge_to_dual_vertices[e].push_back(dual_vertex_idx);
                        }
                        
                    }
                }
            }
            // debug
            // if ([&]() 
            //     {
            //         std::array<Point_3, 8> target = {{
            //             Point_3(-0.362319, 2.82609, 3.11594),
            //             Point_3(-0.217391, 2.82609, 3.11594),
            //             Point_3(-0.362319, 2.97101, 3.11594),
            //             Point_3(-0.217391, 2.97101, 3.11594),
            //             Point_3(-0.362319, 2.82609, 3.26087),
            //             Point_3(-0.217391, 2.82609, 3.26087),
            //             Point_3(-0.362319, 2.97101, 3.26087),
            //             Point_3(-0.217391, 2.97101, 3.26087),
            //         }};
            //         for (int i = 0; i < 8; ++i) {
            //             if (CGAL::squared_distance(corners[i], target[i]) > 1e-8)
            //                 return false;
            //         }
            //         return true;
            //     }())
            // {
            //     // std::cout << "points of triangulation \n";
            //     // for(const auto& t: triangles)
            //     // {
            //     //     for(const auto& p : t)
            //     //     {
            //     //         std::cout << "("
            //     //                 << CGAL::to_double(p.x()) << ", "
            //     //                 << CGAL::to_double(p.y()) << ", "
            //     //                 << CGAL::to_double(p.z()) << ") ";
            //     //     }
            //     //     std::cout << "\n";
            //     // }
            //     // std::cout << "points of patches \n";
            //     // for (const auto& patch: all_patch_ambi)
            //     // {
            //     //     for(const auto& p : patch)
            //     //     {
            //     //         std::cout << "("
            //     //                 << CGAL::to_double(p.x()) << ", "
            //     //                 << CGAL::to_double(p.y()) << ", "
            //     //                 << CGAL::to_double(p.z()) << ") ";
            //     //     }
            //     //     std::cout << "\n";
            //     // }
            //     // // print cell corners
            //     // for (std::size_t i = 0; i < vpc; ++i)
            //     // {
            //     //     const Point_3& p = corners[i];
            //     //     std::cout << p.x() << " " << p.y() << " " << p.z() << "\n";
            //     // }
            //     std::exit(0);
            // }

        }

        else // regular cases
        {
            const int* row = &Cube_table::polygon_cases[i_case * 16];
            int k = 0;

            std::vector<std::vector<int>> all_patch;
            std::vector<std::size_t> patch_dual_vertices_idx;

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

                    const auto& evs = domain.incident_vertices(global_edge);
                    const vertex_descriptor& v0 = evs[0];
                    const vertex_descriptor& v1 = evs[1];
                    const Point_3& p0 = domain.point(v0);
                    const Point_3& p1 = domain.point(v1);
                    const FT val0 = domain.value(v0);
                    const FT val1 = domain.value(v1);

                    Point_3 p;
                    if (!domain.construct_intersection(p0, p1, val0, val1, isovalue, p)) continue;

                    // if (!is_valid_point<Domain>(p)) continue; // skip invalid points
                    if (!is_valid_point<Domain>(p))
                        {
                            std::cerr << "Invalid point detected: p0 (" 
                                      << CGAL::to_double(p0.x()) << ", "
                                      << CGAL::to_double(p0.y()) << ", "
                                      << CGAL::to_double(p0.z()) << ")\n"
                                      << "p1 ("
                                      << CGAL::to_double(p1.x()) << ", "
                                      << CGAL::to_double(p1.y()) << ", "
                                      << CGAL::to_double(p1.z()) << ")\n";
                            continue;
                        }


                    Vector_3 g = domain.gradient(p);
                    patch_points.push_back(p);
                    patch_grads.push_back(g);
                }


                // compute the dual vertex for this patch
                if (patch_points.size() >= 3) 
                {
                    Point_3 dual_vertex;

                    bool success = false;
                    switch(strategy)
                    {
                        case Vertex_strategy::QEM:
                            success = patch_position_QEM(domain, c, patch_points, patch_grads, constrain_to_cell, dual_vertex);
                            break;
                        
                        case Vertex_strategy::Centroid:
                            success = patch_position_centroid(domain, patch_points, dual_vertex);
                            break;
                    }
                    
                    // debug
                    // if (CGAL::squared_distance(dual_vertex, target) < 1e-6)
                    // {
                    //     std::cout << "generated from regular cases\n";
                    //     std::cout << "Found matching dual vertex!\n";
                    //     std::cout << "dual vertex: (" 
                    //             << CGAL::to_double(dual_vertex.x()) << ", "
                    //             << CGAL::to_double(dual_vertex.y()) << ", "
                    //             << CGAL::to_double(dual_vertex.z()) << ")\n";
                    //     std::cout << "  Cell vertices:\n";
                    //     for (std::size_t i = 0; i < vpc; ++i) 
                    //     {
                    //         const Point_3& corner = domain.point(domain.cell_vertices(c)[i]);
                    //         std::cout << "    (" 
                    //                 << CGAL::to_double(corner.x()) << ", "
                    //                 << CGAL::to_double(corner.y()) << ", "
                    //                 << CGAL::to_double(corner.z()) << ")\n";
                    //     }
                    //     std::cout << "  Patch points:\n";
                    //     for (const auto& pt : patch_points) 
                    //     {
                    //         std::cout << "    (" 
                    //                 << CGAL::to_double(pt.x()) << ", "
                    //                 << CGAL::to_double(pt.y()) << ", "
                    //                 << CGAL::to_double(pt.z()) << ")\n";
                    //     }
                    //     std::cout << "case number = " << i_case << "\n";
                    //     // std::exit(0);
                    // }

                    if(success)
                    {
                        std::size_t dual_vertex_idx;
                        auto u = point_to_index.find(dual_vertex);
                        if (u == point_to_index.end()) 
                        {
                            dual_vertex_idx = points.size();
                            points.push_back(dual_vertex);
                            point_to_index[dual_vertex] = dual_vertex_idx;
                            
                            // label the source
                            dual_src[dual_vertex_idx] = DualSrc::Regular;

                        } 
                        else  
                        {
                            dual_vertex_idx = u->second;

                            // if the point is duplicate, find the old tag and update the new one
                            auto it = dual_src.find(dual_vertex_idx);
                            DualSrc prev = (it == dual_src.end()) ? DualSrc::Regular : it->second;
                            dual_src[dual_vertex_idx] = upgrade_src(prev, DualSrc::Regular);
                        }


                        patch_dual_vertices_idx.push_back(dual_vertex_idx); // collect a dual vertex 
                        all_patch.push_back(current_patch_edge_local_idx); // collect a patch


                    }
                }

                if (row[k] == -2) ++k;
            }

            CGAL_assertion(patch_dual_vertices_idx.size() == all_patch.size());

            // map edge to dual vertices
            // the ordering of patch_dual_vertices_idx is consistent with all_patch
            int patch_size = patch_dual_vertices_idx.size();
            for (int patch_idx = 0; patch_idx < patch_size; patch_idx++)
            {
                // loop over all edges in this patch
                for (int local_edge_idx : all_patch[patch_idx]) 
                {
                    auto global_edge = cell_edges[local_edge_idx];
                    edge_to_dual_vertices[global_edge].push_back(patch_dual_vertices_idx[patch_idx]);
                }
                cell_to_dual_vertices[c].push_back(patch_dual_vertices_idx[patch_idx]); 
            }
        }
        
    };

    domain.for_each_cell(cell_patch_positioner);

    // reorder the dual vertices of each edge to be in local cyclic order
    reorder_edge_to_dual_vertices(domain, cell_to_dual_vertices, edge_to_dual_vertices);

    // // debug
    // const Point_3 v1(-2.10145, -2.10145, 2.10145);
    // const Point_3 v2(-2.10145, -2.10145, 2.24638);

    // generate quads
    auto generate_quad = [&](const edge_descriptor& e) 
    {
        // // debug
        // auto it = edge_to_dual_vertices.find(e);
        // if (it != edge_to_dual_vertices.end())
        // {
        //     const std::vector<std::size_t>& dual_vertices_idx = it->second;
        //     std::cout << "Edge with " << dual_vertices_idx.size() << " dual vertices:\n";
        //     const auto& evs = domain.incident_vertices(e);
        //     const Point_3& p0 = domain.point(evs[0]);
        //     const Point_3& p1 = domain.point(evs[1]);
        //      if ((CGAL::squared_distance(v1, p0) < 1e-8 && CGAL::squared_distance(v2, p1) < 1e-8) ||
        //     (CGAL::squared_distance(v1, p1) < 1e-8 && CGAL::squared_distance(v2, p0) < 1e-8))
        //     {
        //         const std::vector<std::size_t>& dual_vertices_idx = it->second;
        //         std::cout << "Matched edge with " << dual_vertices_idx.size() << " dual vertices:\n";
        //         for (std::size_t idx : dual_vertices_idx) 
        //         {
        //             const auto& pt = points[idx];
        //             std::cout << "    Dual vertex: (" 
        //                     << CGAL::to_double(pt.x()) << ", " 
        //                     << CGAL::to_double(pt.y()) << ", " 
        //                     << CGAL::to_double(pt.z()) << ")\n";
        //         }
        //         std::exit(0); 
        //     }   
        // }

        generate_quad_dmc(e, domain, edge_to_dual_vertices, polygons);
    };
    domain.for_each_edge(generate_quad);

    // count how many polygons share each dual edge 
    std::unordered_map<
        std::pair<std::size_t,std::size_t>, 
        int, 
        boost::hash<std::pair<std::size_t,std::size_t>>
    > dual_edge_use;

    // normalize orientation
    auto norm_edge = [](std::size_t a, std::size_t b) 
    {
        return (a < b) ? std::pair<std::size_t,std::size_t>(a,b)
                    : std::pair<std::size_t,std::size_t>(b,a);
    };

    for (const auto& poly : polygons) 
    {
        const int m = static_cast<int>(poly.size());   // triangles or quads (or n-gons)
        for (int i = 0; i < m; ++i) {
            auto e = norm_edge(poly[i], poly[(i+1)%m]); //%m to wrap around
            dual_edge_use[e] += 1;
        }
    }

    // report nonmanifold edgesz
    for (const auto& [e, cnt] : dual_edge_use) 
    {
        auto [a,b] = e;

        if (cnt == 1) std::cout << "[DUAL-EDGE boundary] (" << a << "," << b << ") faces=" << cnt << "\n";
        else if (cnt >= 3) 
        {
            std::cout << "[DUAL-EDGE NONMANIFOLD] ...\n";

            // the the source of the dual vertex
            auto get_src = [&](std::size_t dual_vertex_index) 
            {
                if (auto it = dual_src.find(dual_vertex_index); it != dual_src.end()) return it->second;

                // if somehow cannot find the tag, warn and default to regular
                std::cout << "dual vertex index" << dual_vertex_index << "does not have a tag \n";
                return DualSrc::Regular;
            };

            auto src_a = get_src(a);
            auto src_b = get_src(b);
            std::cout << "[DUAL-EDGE NONMANIFOLD] (" << a << "," << b << ") incident_faces=" << cnt << "\n";
            // 0 = Regular; 1 = TMC_NoTunnel; 2 = TMC_Tunnel
            std::cout << "  sources: a=" << int(src_a) << " b=" << int(src_b) << "\n";
        }
    }

    
}


} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif