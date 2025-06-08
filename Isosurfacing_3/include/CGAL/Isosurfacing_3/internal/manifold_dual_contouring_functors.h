#ifndef CGAL_ISOSURFACING_3_INTERNAL_MANIFOLD_DUAL_CONTOURING_H
#define CGAL_ISOSURFACING_3_INTERNAL_MANIFOLD_DUAL_CONTOURING_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/tables.h>

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

namespace CGAL {
namespace Isosurfacing {
namespace internal {

// compute the dual vertices of each patch of a cell
template <typename Domain, 
          typename Vertices,
          typename EdgeToPointIDMap>
std::vector<typename Domain::Geom_traits::Point_3>
patch_position_per_cell_QEM(const Domain& domain,
                            const typename Domain::cell_descriptor& c,
                            std::size_t i_case,
                            const Vertices& vertices,
                            const EdgeToPointIDMap& edge_to_point_id, 
                            bool constrain_to_cell)
{
    using Point_3 = typename Domain::Geom_traits::Point_3;
    using Vector_3 = typename Domain::Geom_traits::Vector_3;
    std::vector<Point_3> patch_vertices;

    std::vector<std::size_t> current_patch_edges;

    for (int k = 0; k < 16; ++k) 
    {
        int e_id =Cube_table::polygon_cases[i_case*16 + k];
        if (e_id == -1) break;
        if (e_id == -2) 
        {
            if (!current_patch_edges.empty()) 
            {
                // collect the points and gradients for the current patch
                std::vector<Point_3> patch_points;
                std::vector<Vector_3> patch_gradients;
                for (std::size_t i = 0; i < current_patch_edges.size(); i++) {
                    patch_points.push_back(vertices[e_id]);
                    patch_gradients.push_back(domain.gradient(vertices[e_id]));
                }
                Point_3 dual_vertex;
                // compute the dual vertex
                bool success = cell_position_QEM(c, domain, constrain_to_cell,
                    edge_to_point_id, patch_points, patch_gradients, dual_vertex);
                if (success) patch_vertices.push_back(dual_vertex);
                current_patch_edges.clear();
            }
            continue;
        }
        current_patch_edges.push_back(e_id);
    }

    // compute the vertex for the last patch
    if (!current_patch_edges.empty()) 
    {
        std::vector<Point_3> patch_points;
        std::vector<Vector_3> patch_gradients;
        for (std::size_t i = 0; i < current_patch_edges.size(); i++) {
            auto e_id = current_patch_edges[i];
            patch_points.push_back(vertices[e_id]);
            patch_gradients.push_back(domain.gradient(vertices[e_id]));
        }
        Point_3 dual_vertex;
        bool success = cell_position_QEM(c, domain, constrain_to_cell,
            edge_to_point_id, patch_points, patch_gradients, dual_vertex);
        if (success) patch_vertices.push_back(dual_vertex);
    }
    return patch_vertices;
}


} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL