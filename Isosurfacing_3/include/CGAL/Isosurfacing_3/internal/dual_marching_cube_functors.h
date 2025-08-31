#ifndef CGAL_ISOSURFACING_3_INTERNAL_DUAL_MARCHING_CUBE_H
#define CGAL_ISOSURFACING_3_INTERNAL_DUAL_MARCHING_CUBE_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/tables.h>
#include <CGAL/Isosurfacing_3/internal/tmc_dmc_functors.h>
#include <CGAL/Isosurfacing_3/internal/marching_cubes_functors.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>

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
#include <unordered_set>
#include <algorithm>
#include <numeric>

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

    #include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/number_utils.h> // CGAL::to_double
#include <fstream>
#include <string>

// Writes ALL voxel cells as independent hexahedra (6 quads each).
// Grid must provide xdim(), ydim(), zdim() (vertex counts) and point(i,j,k).
template <typename Grid>
bool write_all_cells_off(const Grid& grid, const std::string& path)
{
    const std::size_t X = grid.xdim(); // vertex dims
    const std::size_t Y = grid.ydim();
    const std::size_t Z = grid.zdim();
    if (X < 2 || Y < 2 || Z < 2) return false;

    const std::size_t nx = X - 1;      // cell dims
    const std::size_t ny = Y - 1;
    const std::size_t nz = Z - 1;

    const std::size_t n_cells  = nx * ny * nz;
    const std::size_t n_verts  = n_cells * 8;  // per cell (dup vertices)
    const std::size_t n_faces  = n_cells * 6;  // 6 quads per cell

    std::ofstream out(path);
    if (!out) return false;

    out << "OFF\n";
    out << n_verts << " " << n_faces << " 0\n";

    // Emit vertices: for each cell, write its 8 corners in fixed order.
    // Corner order (local indices 0..7):
    // 0:(i,j,k)   1:(i+1,j,k)   2:(i+1,j+1,k)   3:(i,j+1,k)
    // 4:(i,j,k+1) 5:(i+1,j,k+1) 6:(i+1,j+1,k+1) 7:(i,j+1,k+1)
    for (std::size_t k = 0; k < nz; ++k)
      for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t i = 0; i < nx; ++i)
        {
            const auto p000 = grid.point(i    , j    , k    );
            const auto p100 = grid.point(i + 1, j    , k    );
            const auto p110 = grid.point(i + 1, j + 1, k    );
            const auto p010 = grid.point(i    , j + 1, k    );
            const auto p001 = grid.point(i    , j    , k + 1);
            const auto p101 = grid.point(i + 1, j    , k + 1);
            const auto p111 = grid.point(i + 1, j + 1, k + 1);
            const auto p011 = grid.point(i    , j + 1, k + 1);

            const auto emit = [&](const auto& P){
                out << CGAL::to_double(P.x()) << " "
                    << CGAL::to_double(P.y()) << " "
                    << CGAL::to_double(P.z()) << "\n";
            };

            emit(p000); emit(p100); emit(p110); emit(p010);
            emit(p001); emit(p101); emit(p111); emit(p011);
        }

    // Emit faces: 6 quads per cell, referencing the 8 local verts.
    // Outward orientation for an axis-aligned right-handed grid.
    // Face order: z- (bottom), z+ (top), x- , x+ , y- , y+
    std::size_t base = 0;
    for (std::size_t c = 0; c < n_cells; ++c, base += 8)
    {
        // bottom z- : 0 1 2 3
        out << "4 " << base+0 << " " << base+1 << " " << base+2 << " " << base+3 << "\n";
        // top z+    : 4 7 6 5
        out << "4 " << base+4 << " " << base+7 << " " << base+6 << " " << base+5 << "\n";
        // x-        : 0 3 7 4
        out << "4 " << base+0 << " " << base+3 << " " << base+7 << " " << base+4 << "\n";
        // x+        : 1 5 6 2
        out << "4 " << base+1 << " " << base+5 << " " << base+6 << " " << base+2 << "\n";
        // y-        : 0 4 5 1
        out << "4 " << base+0 << " " << base+4 << " " << base+5 << " " << base+1 << "\n";
        // y+        : 3 2 6 7
        out << "4 " << base+3 << " " << base+2 << " " << base+6 << " " << base+7 << "\n";
    }

    return true;
}

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


// hash for std::pair<size_t,size_t>
struct PairHash 
{
    std::size_t operator()(const std::pair<std::size_t,std::size_t>& p) const noexcept
    {
        return std::hash<std::size_t>{}(p.first)^(std::hash<std::size_t>{}(p.second) << 1);
    } 
};

struct FaceRec 
{
    std::array<std::size_t,4> quad_verts;   // {v0,v1,v2,v3} from the original face
    std::array<std::size_t,4> ring_centers; // {c01,c12,c23,c30} in the same order
};

// subdivide each quad into 5 quads 
// remove all incident quads
template <typename Domain,
          typename PolygonRange>
void subdivide_quads_skip_incident(const Domain& domain,
                                    const std::vector<std::size_t>& face_ids,
                                    const std::vector<std::pair<std::size_t,std::size_t>>& nme, 
                                    std::vector<typename Domain::Geom_traits::Point_3>& points,
                                    PolygonRange& polygons,
                                    std::unordered_map<std::size_t, FaceRec>& old_fid_to_FaceRec, // for vertex merging later
                                    std::vector<std::vector<std::size_t>>& to_add)
{
    using Point_3 = typename Domain::Geom_traits::Point_3;

    auto norm_edge = [](std::size_t a, std::size_t b) 
    {
        return (a < b) ? std::pair<std::size_t,std::size_t>(a,b)
                    : std::pair<std::size_t,std::size_t>(b,a);
    };

    // hashset for fast lookup
    std::unordered_set<std::pair<std::size_t,std::size_t>, PairHash> nme_set;
    nme_set.reserve(nme.size());
    for (const auto& e: nme) nme_set.insert(e);

    auto is_nme = [&](std::size_t a, std::size_t b) -> bool {return nme_set.count(norm_edge(a,b));}; //if found, return 1 else 0

    // compute the centroid and return its index
    auto centroid3 = [&](std::size_t a, std::size_t b, std::size_t c) -> std::size_t 
    {
        std::vector<Point_3> patch_points = {points[a], points[b], points[c]};
        Point_3 p;
        bool success = patch_position_centroid(domain, patch_points, p);
        if (!success) 
        {
            std::cout << "WARN: patch_position_centroid not success\n";
            return 0;
        }
        points.push_back(p);
        return points.size() - 1;
    };

    std::vector<std::size_t> to_erase;
    to_erase.reserve(face_ids.size());

    // std::vector<std::vector<std::size_t>> to_add;
    // to_add.reserve(face_ids.size() * 4);

    for(std::size_t fid: face_ids)
    {
        if (fid >= polygons.size()) {
            std::cout << "WARN: face id " << fid << " out of range; skipping\n";
            continue;
        }

        const auto& face = polygons[fid];
        if (face.size() != 4)
        {
            std::cout << "WARN: face " << fid << " is not a quad; skipping subdivision\n";
            continue;
        }

        const std::size_t v0 = face[0];
        const std::size_t v1 = face[1];
        const std::size_t v2 = face[2];
        const std::size_t v3 = face[3];

        const std::size_t c01 = centroid3(v0, v1, v2);
        const std::size_t c12 = centroid3(v1, v2, v3);
        const std::size_t c23 = centroid3(v2, v3, v0);
        const std::size_t c30 = centroid3(v3, v0, v1);
        
        old_fid_to_FaceRec.emplace(fid, FaceRec{{v0,v1,v2,v3},{c01,c12,c23,c30}});

        // center quad
        to_add.push_back({ c01, c12, c23, c30 });

        // add each side quad only if its outer edge is NOT an NME
        if (!is_nme(v0, v1)) to_add.push_back({ v0, v1, c01, c30 });
        if (!is_nme(v1, v2)) to_add.push_back({ v1, v2, c12, c01 });
        if (!is_nme(v2, v3)) to_add.push_back({ v2, v3, c23, c12 });
        if (!is_nme(v3, v0)) to_add.push_back({ v3, v0, c30, c23 });

        to_erase.push_back(fid);
    }

    // erase originals in descending order 
    std::sort(to_erase.begin(), to_erase.end());
    to_erase.erase(std::unique(to_erase.begin(), to_erase.end()), to_erase.end());
    for (auto it = to_erase.rbegin(); it != to_erase.rend(); ++it) {
         // convert index to iterator and erase that polygon
        auto poly_it = polygons.begin() + static_cast<std::ptrdiff_t>(*it);
        polygons.erase(poly_it);
    }

    // // append new quads
    // for (auto& f : to_add)
    //     polygons.push_back(std::move(f));
}

// get the vertices of the associated primal edge given the 4 cells of the 4 dual vertices of a quad
template <typename Domain>
std::pair<typename Domain::vertex_descriptor,
          typename Domain::vertex_descriptor>
primal_edge_for_quad(const Domain& domain,
                     const typename Domain::cell_descriptor& c0,
                     const typename Domain::cell_descriptor& c1,
                     const typename Domain::cell_descriptor& c2,
                     const typename Domain::cell_descriptor& c3)
{
    using vertex_descriptor = typename Domain::vertex_descriptor;

    const auto& v0 = domain.cell_vertices(c0);
    const auto& v1 = domain.cell_vertices(c1);
    const auto& v2 = domain.cell_vertices(c2);
    const auto& v3 = domain.cell_vertices(c3);

    auto contains = [&](const auto& arr, const auto& v) -> bool 
    {
        for (const auto& a : arr) if (a == v) return true;
        return false; 
    };

    std::vector<vertex_descriptor> common;
    common.reserve(2);

    // pick vertices that are present in all 4 cells
    for (const vertex_descriptor& a : v0)
    {
        if (contains(v1, a) && contains(v2, a) && contains(v3, a))
        {
            // check duplicate
            bool duplicate = false;
            for (const auto& c : common) if (c == a) { duplicate = true; break; }
            if (!duplicate) common.push_back(a);
        }
    }

    // Expect exactly 2 primal vertices (a primal edge)
    CGAL_assertion(common.size() == 2);

    return {common[0], common[1]};
}

// groups  : vector of disjoint sets  (e.g. {{A,B}, {D,E,F}, …})
// u, v    : two indices that must end up in the same set
static void merge_into(std::vector<std::vector<std::size_t>>& groups,std::size_t u, std::size_t v)
{
    if (u == v) return;                     

    int igu = -1, igv = -1;                     // group indices containing u and v
    for (int i = 0; i < (int)groups.size(); ++i)
    {
        const auto& G = groups[i];
        if (igu == -1 && std::find(G.begin(), G.end(), u) != G.end()) igu = i;
        if (igv == -1 && std::find(G.begin(), G.end(), v) != G.end()) igv = i;
        if (igu != -1 && igv != -1) break;
    }

    if (igu == -1 && igv == -1)            // neither in a group → create new
        groups.push_back({u, v});

    else if (igu != -1 && igv == -1)       // u in group, v not
        groups[igu].push_back(v);

    else if (igu == -1 && igv != -1)       // v in group, u not
        groups[igv].push_back(u);

    else if (igu != igv)                   // both in different groups then merge
    {
        if (groups[igu].size() < groups[igv].size()) std::swap(igu, igv);
        auto& big  = groups[igu];
        auto& small = groups[igv];
        big.insert(big.end(), small.begin(), small.end()); // alwatys merge the small size group into the big one to save total moves
        groups.erase(groups.begin() + igv);
    }
}


// compute primal faces for each nme
template <typename Domain>
std::vector<typename Domain::vertex_descriptor>
primal_face(const Domain& domain,
                const typename Domain::cell_descriptor& cell_a,
                const typename Domain::cell_descriptor& cell_b)
{
    using vertex_descriptor = typename Domain::vertex_descriptor;

    const auto& verts_a = domain.cell_vertices(cell_a);
    const auto& verts_b = domain.cell_vertices(cell_b);

    std::array<vertex_descriptor,4> common{};

    int k = 0;

    // find the 4 common primal vertices (the separating face)
    for (const vertex_descriptor& va : verts_a) 
    {
        for (const vertex_descriptor& vb : verts_b) 
        {
            if (va == vb) 
            {
                if (k < 4) common[k++] = va;
                break;
            }
        }
    }

    CGAL_assertion(k==4);

    auto in_common = [&](const vertex_descriptor& v)
    {
        for (int i=0;i<4;++i) if (common[i] == v) return true;
        return false;
    };

    std::vector<vertex_descriptor> ordered;
    ordered.reserve(4);
    ordered.push_back(common[0]);

    const auto& cell_edges = domain.cell_edges(cell_a);

    for (int added = 1; added < 4;) 
    {
        const auto curr = ordered.back();
        bool advanced = false;

        for (const auto& ce : cell_edges) 
        {
            const auto& ev = domain.incident_vertices(ce);
            // check both directions; stay on the separating face using in_common
            if (ev[0] == curr && in_common(ev[1]) &&
                std::find(ordered.begin(), ordered.end(), ev[1]) == ordered.end())
            {
                ordered.push_back(ev[1]);
                ++added;
                advanced = true;
                break;
            }
            if (ev[1] == curr && in_common(ev[0]) &&
                std::find(ordered.begin(), ordered.end(), ev[0]) == ordered.end())
            {
                ordered.push_back(ev[0]);
                ++added;
                advanced = true;
                break;
            }
        }

        CGAL_assertion(advanced); // prevents infinite loop 
        if (!advanced) break;     
    }

    CGAL_assertion(ordered.size() == 4);

    return ordered;
}



// pairing direction for asymptotic decider
enum class AD_pair{D02, D13};

// decide the correct pairing based on asymptotic decider
// S = (f0 - iso)*(f2 - iso) - (f1 - iso)*(f3 - iso)
// if S > 0 → AD_pair::D02  
//      – choose diagonal pf[0]–pf[2]  
//      – edge groupings: {E01, E12} and {E23, E30}  
// if S < 0 → AD_pair::D13  
//      – choose diagonal pf[1]–pf[3]  
//      – edge groupings: {E01, E30} and {E12, E23}  
//  if S ≈ 0 → tie case, resolved deterministically (here defaults to D02)
template <typename Domain>
AD_pair asymptotic_decider(const Domain& domain,
                           const std::vector<typename Domain::vertex_descriptor>& pf,
                           const typename Domain::Geom_traits::FT iso)
{
    using FT = typename Domain::Geom_traits::FT;

    const FT f0 = domain.value(pf[0]); // pf order must be a 0-1-2-3 loop
    const FT f1 = domain.value(pf[1]);
    const FT f2 = domain.value(pf[2]);
    const FT f3 = domain.value(pf[3]);

    const FT S = (f0 - iso)*(f2 - iso) - (f1 - iso)*(f3 - iso);

    const FT tol = FT(0);
    if (S >  tol) return AD_pair::D02; // diagonal (0,2) → groups {E01,E12} & {E23,E30}
    if (S < -tol) return AD_pair::D13; // diagonal (1,3) → groups {E01,E30} & {E12,E23}
    return AD_pair::D02;               // deterministic tiebreak
}



// for each nme, merge vertices of incident quads after subdvision according to asymptotic decider
template <typename Domain, typename PairHash>
void centers_to_merge_for_nme(const Domain& domain,
                        const std::pair<std::size_t, std::size_t> e, //nme; [a,b] must be normalized already: a<b
                        const std::unordered_map<std::pair<std::size_t,std::size_t>, std::vector<std::size_t>, PairHash>& edge_to_faces,
                        const std::unordered_map<std::size_t, FaceRec>& face_rec,  
                        const std::vector<typename Domain::cell_descriptor>& dual_vertex_cell,
                        const typename Domain::Geom_traits::FT isovalue,
                        std::vector<std::vector<std::size_t>>& to_be_merged_points,
                    std::vector<typename Domain::Geom_traits::Point_3>& points)
{
    using vertex_descriptor = typename Domain::vertex_descriptor;

    auto [a,b] = e;

    // incident faces to this nme e
    auto it = edge_to_faces.find(e);
    CGAL_assertion(it != edge_to_faces.end()); 
    const auto& faces = it->second;

    if (faces.size()!=4) return; // if less then 4 incident faces, then the mesh is inherently nonmanifold, nothing to do

    // find the index of a value 
    auto find_idx = [](const std::array<std::size_t, 4>& v, std::size_t key)->int
    {
        for (int i=0;i<(int)v.size();++i) if (v[i]==key) return i;
        return -1;
    };

    // For each incident face, locate which edge is (a,b) in that face’s v0..v3 cycle,
    // then pick the two adjacent side centers (indices k and (k+3)%4).
    struct FaceSide 
    {
        std::array<std::size_t,2> side_centers;    // the two centers along the (a,b) edge of this face
        std::array<vertex_descriptor, 2>  primal_edge;     // the primal edge this face sits on (2 primal verts)
    };
    std::vector<FaceSide> per_face; per_face.reserve(4);

    for (std::size_t fid : faces)
    {
        const auto iterat = face_rec.find(fid);
        CGAL_assertion(iterat != face_rec.end());
        const auto& rec = iterat->second;

        // original dual-vertex cycle of that quad
        const std::array<std::size_t,4>& quad_verts = rec.quad_verts;            // size 4: {v0,v1,v2,v3}
        const std::array<std::size_t, 4>& centroid_verts = rec.ring_centers;         // size 4: {c01,c12,c23,c30}

        // locate the position where edge (V[k], V[k+1]) == (a,b) (up to orientation)
        int ia = find_idx(quad_verts, a);
        int ib = find_idx(quad_verts, b);
        CGAL_assertion(ia>=0 && ib>=0);
        // ensure they are adjacent on the 4-cycle
        bool adjab = ((ia+1)%4 == ib); 
        bool adjba = ((ib+1)%4 == ia);
        CGAL_assertion(adjab || adjba);

        int k = adjab ? ia : ib; // edge is (V[k], V[k+1]) but orientation can go either way

        // the two side centers adjacent to that edge on this face:
        // consistent with the face vertices ordering
        std::array<std::size_t,2> sides = { centroid_verts[k], centroid_verts[ (k+3) % 4 ] };

        //debug print
        const auto q0 = points[sides[0]];
        const auto q1 = points[sides[1]];
        std::cout << "[face " << fid << "] k=" << k
            << " sides={" << sides[0] << "," << sides[1] << "}\n"
            << "   s[0]: (" << CGAL::to_double(q0.x()) << ", "
                            << CGAL::to_double(q0.y()) << ", "
                            << CGAL::to_double(q0.z()) << ")\n"
            << "   s[1]: (" << CGAL::to_double(q1.x()) << ", "
                            << CGAL::to_double(q1.y()) << ", "
                            << CGAL::to_double(q1.z()) << ")\n";

        const auto c0 = dual_vertex_cell[quad_verts[0]];
        const auto c1 = dual_vertex_cell[quad_verts[1]];
        const auto c2 = dual_vertex_cell[quad_verts[2]];
        const auto c3 = dual_vertex_cell[quad_verts[3]];
        auto pe = primal_edge_for_quad(domain, c0,c1,c2,c3); // get the primal edge associated with this quad

        //debug
        // const auto& v0 = pe.first;
        // const auto& v1 = pe.second;
        // const auto& P0 = domain.point(v0);
        // const auto& P1 = domain.point(v1);
        // std::cout << "Primal edge : \n"
        //         << "  v0 -> (" << CGAL::to_double(P0.x()) << ", "
        //                         << CGAL::to_double(P0.y()) << ", "
        //                         << CGAL::to_double(P0.z()) << ")\n"
        //         << "  v1 -> (" << CGAL::to_double(P1.x()) << ", "
        //                         << CGAL::to_double(P1.y()) << ", "
        //                         << CGAL::to_double(P1.z()) << ")\n";

        per_face.push_back({sides, {pe.first, pe.second}});
    }

    const auto& cell_a = dual_vertex_cell[a];
    const auto& cell_b = dual_vertex_cell[b];
    auto pf = primal_face(domain, cell_a, cell_b); 

    AD_pair pairing = asymptotic_decider(domain, pf, isovalue); 

    // debug
    using FT = typename Domain::Geom_traits::FT;
    const auto P0 = domain.point(pf[0]); const FT V0 = domain.value(pf[0]);
    const auto P1 = domain.point(pf[1]); const FT V1 = domain.value(pf[1]);
    const auto P2 = domain.point(pf[2]); const FT V2 = domain.value(pf[2]);
    const auto P3 = domain.point(pf[3]); const FT V3 = domain.value(pf[3]);
    std::cout << "[AD audit] pf corners (id -> coord, value):\n";
    auto pv = [&](int i, const auto& P, FT V){
        std::cout << "  pf["<<i<<"] -> (" << CGAL::to_double(P.x()) << ", "
                                         << CGAL::to_double(P.y()) << ", "
                                         << CGAL::to_double(P.z()) << "), "
                  << "val=" << CGAL::to_double(V) << "\n";
    };
    pv(0,P0,V0); pv(1,P1,V1); pv(2,P2,V2); pv(3,P3,V3);
    std::cout << "Asymptotic decider pairing: " << ((pairing == AD_pair::D02) ? "D02" : "D13") << "\n";
    std::cout << "isovalue = " << CGAL::to_double(isovalue) << "\n";

    // helper to compare unordered primal edges
    auto same_edge = [](const std::array<vertex_descriptor,2>& e, vertex_descriptor x, vertex_descriptor y)
    {
        return ( (e[0]==x && e[1]==y) || (e[0]==y && e[1]==x) );
    };

    //    E0=(pf[0],pf[1]), E1=(pf[1],pf[2]), E2=(pf[2],pf[3]), E3=(pf[3],pf[0])
    std::array<std::size_t,2> side_centers_E0{SIZE_MAX, SIZE_MAX}; // store face-centers associated with E0
    std::array<std::size_t,2> side_centers_E1{SIZE_MAX, SIZE_MAX}; // store face-centers associated with E1
    std::array<std::size_t,2> side_centers_E2{SIZE_MAX, SIZE_MAX}; // store face-centers associated with E2
    std::array<std::size_t,2> side_centers_E3{SIZE_MAX, SIZE_MAX}; // store face-centers associated with E3


    for (const auto& fs : per_face) // fs = FaceSide struct
    {
        const auto& e = fs.primal_edge;
        if (same_edge(e, pf[0],pf[1])) side_centers_E0 = fs.side_centers;// E0 
        else if (same_edge(e, pf[1],pf[2])) side_centers_E1 = fs.side_centers;// E1 
        else if (same_edge(e, pf[2],pf[3])) side_centers_E2 = fs.side_centers;// E2
        else if (same_edge(e, pf[3],pf[0])) side_centers_E3 = fs.side_centers;// E3
        else {CGAL_assertion(false);}
    }
    CGAL_assertion(side_centers_E0[0]!=SIZE_MAX && side_centers_E0[1]!=SIZE_MAX && side_centers_E1[0]!=SIZE_MAX
    && side_centers_E1[1]!=SIZE_MAX && side_centers_E2[0]!=SIZE_MAX && side_centers_E2[1]!=SIZE_MAX && side_centers_E3[0]!=SIZE_MAX && side_centers_E3[1]!=SIZE_MAX);

    // vertex pairing based on the result from asymptotic decider
    if (pairing == AD_pair::D02) // {E0, E1}, {E2, E3}
    { 
        // since the ordering (either CCW/CW) is consistent, the side_centers goes around and thus has flipped pairs
        // we need to merge e.g. E0[0] with E1[1], E0[1] with E1[0], etc.
        merge_into(to_be_merged_points, side_centers_E0[0], side_centers_E1[1]);
        merge_into(to_be_merged_points, side_centers_E0[1], side_centers_E1[0]);
        merge_into(to_be_merged_points, side_centers_E2[0], side_centers_E3[1]);
        merge_into(to_be_merged_points, side_centers_E2[1], side_centers_E3[0]);
    } 
    else 
    { // AD_pair::D13 // {E0, E3}, {E1, E2}
        merge_into(to_be_merged_points, side_centers_E0[0], side_centers_E3[1]);
        merge_into(to_be_merged_points, side_centers_E0[1], side_centers_E3[0]);
        merge_into(to_be_merged_points, side_centers_E2[0], side_centers_E1[1]);
        merge_into(to_be_merged_points, side_centers_E2[1], side_centers_E1[0]);
    }
}




// update polygons and points to resolve manifold edges
// merge points at their centroid
template <typename Domain, typename PolygonRange>
void update_quads(const Domain& domain,
                  const std::vector<std::vector<std::size_t>> to_be_merged_points,
                  std::vector<typename Domain::Geom_traits::Point_3>& points,
                  std::vector<std::vector<std::size_t>>& to_add,
                  PolygonRange& polygons)
{
    using Point_3  = typename Domain::Geom_traits::Point_3;

    const std::size_t n = points.size();

    //build a representative map with identity e.g. rep_of[4] = 4
    std::vector<std::size_t> rep_of(n);
    std::iota(rep_of.begin(), rep_of.end(), 0);

    // place each merge group at its centroid, update points[rep]
    for (const auto& group: to_be_merged_points)
    {
        if (group.empty()) continue;

        // the group is reprsented by the smallest index rep
        const std::size_t rep = *std::min_element(group.begin(), group.end());

        // gather patch points
        std::vector<Point_3> patch_pts;
        patch_pts.reserve(group.size());
        for (std::size_t v : group) {
            patch_pts.push_back(points[v]);
            rep_of[v] = rep;               // redirect every group member to rep
        }

        // compute centroid; if it ever returns false, keep current rep position
        Point_3 p = points[rep];
        if (patch_position_centroid(domain, patch_pts, p)) points[rep] = p; // write centroid into representative
    }

    std::vector<std::size_t> new_index(n, static_cast<std::size_t>(-1));
    std::vector<Point_3> new_points;
    new_points.reserve(n);

    for (int i = 0; i < n; ++i) {
        if (rep_of[i] == i) // identity representative, i.e. the unaffected
        {                      
            new_index[i] = static_cast<std::size_t>(new_points.size());
            new_points.push_back(points[i]);       // copy unaffected coord
        }
    }
    // give affected the index of their representative
    for (int i = 0; i < n; ++i) 
    {
        if (new_index[i] == static_cast<std::size_t>(-1)) new_index[i] = new_index[ rep_of[i] ];
    }

    points = std::move(new_points); // replace new_points by points and discard new_ponts

    // helper to remap old index to new index
    auto remap_poly = [&](std::vector<std::size_t>& q) 
    {
        for (auto& v : q) v = new_index[v];
        CGAL_assertion(q[0]!=q[1] && q[1]!=q[2] && q[2]!=q[3] && q[3]!=q[0] && q[0]!=q[2] && q[1]!=q[3]);
    };

    for (auto& q : polygons) remap_poly(q);
    for (auto& q : to_add)   remap_poly(q);

    // append to_add to polygons
    polygons.insert(polygons.end(),
                    std::make_move_iterator(to_add.begin()),
                    std::make_move_iterator(to_add.end()));
    to_add.clear(); // dont need to_add anymore
}



// use orient_polygon_soup to return all nonmanifold vertices
// not optimal here since we create copies of points and polygons
struct NMVVisitor : CGAL::Polygon_mesh_processing::Default_orientation_visitor 
{
  std::unordered_set<std::size_t>& out;  
  explicit NMVVisitor(std::unordered_set<std::size_t>& out_) : out(out_) {}
  void non_manifold_vertex(std::size_t v, std::size_t) { out.insert(v); }
};

template <typename PointRange, typename PolygonRange>
std::unordered_set<std::size_t>
find_nonmanifold_vertices(const PointRange& points, const PolygonRange& polygons)
{
  auto pts   = points;   // run on copies; the routine mutates inputs
  auto polys = polygons;

  std::unordered_set<std::size_t> nmv;
  NMVVisitor vis(nmv);
  CGAL::Polygon_mesh_processing::orient_polygon_soup(pts, polys, CGAL::parameters::visitor(vis));

  return nmv;            
}


// dimension of grid
struct gridDim
{
    std::size_t nx, ny, nz;
};

// check if a cell is on the boundary of the grid
template<typename Domain>
inline bool is_boundary_cell(const typename Domain::cell_descriptor& c, const gridDim& d) 
{
    return (c[0] == 0 || c[0] == d.nx-1 ||
            c[1] == 0 || c[1] == d.ny-1 ||
            c[2] == 0 || c[2] == d.nz-1);
}

// to determine which sides of the cell should be padded
enum cellSide : uint8_t 
{
  XMIN=1<<0, XMAX=1<<1, YMIN=1<<2,
  YMAX=1<<3, ZMIN=1<<4, ZMAX=1<<5
};

template<typename Domain>
uint8_t which_cellSide(const typename Domain::cell_descriptor& c, const gridDim& d) 
{
    uint8_t m = 0;
    if (c[0] == 0)      m |= XMIN;
    if (c[0] == d.nx-1) m |= XMAX;
    if (c[1] == 0)      m |= YMIN;
    if (c[1] == d.ny-1) m |= YMAX;
    if (c[2] == 0)      m |= ZMIN;
    if (c[2] == d.nz-1) m |= ZMAX;
    return m;
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

struct PostProcessOff{};
struct PostProcessOn{};

// topologically correct dual marching cube
// dual marching cube
template <typename Domain, 
          typename PolygonRange>
void dual_marching_cubes_tmc(const Domain& domain,
                             const typename Domain::Geom_traits::FT isovalue,
                             std::vector<typename Domain::Geom_traits::Point_3>& points,
                             PolygonRange& polygons,
                             bool constrain_to_cell,
                             Isosurfacing::Vertex_strategy strategy,
                             PostProcessOff)
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
    std::vector<cell_descriptor> dual_vertex_cell; // use later for track primal face for asymptotic decider


    // document the source of each dual vertex for nonmanifold edge tracking
    enum class DualSrc : uint8_t { Regular, TMC_NoTunnel, TMC_Tunnel };

    // prefer Tunnel > NoTunnel > Regular when the same point is reused
    auto upgrade_src = [](DualSrc oldv, DualSrc newv) {
        if (oldv == DualSrc::TMC_Tunnel || newv == DualSrc::TMC_Tunnel) return DualSrc::TMC_Tunnel;
        if (oldv == DualSrc::TMC_NoTunnel || newv == DualSrc::TMC_NoTunnel) return DualSrc::TMC_NoTunnel;
        return DualSrc::Regular;
    };

    std::unordered_map<std::size_t, DualSrc> dual_src;

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

                            if (point_lies_on_segment(domain, p,e))
                            {
                                current_patch_edges.push_back(e);
                                break;
                            }
                            
                        }
                    }

                    if (filtered_patch_points.size() < 3) continue; 

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

                    if(success)
                    {
                        std::size_t dual_vertex_idx;
                        auto u = point_to_index.find(dual_vertex);
                        if (u == point_to_index.end()) // if no duplicates
                        {
                            dual_vertex_idx = points.size();
                            points.push_back(dual_vertex);
                            point_to_index[dual_vertex] = dual_vertex_idx;
                            dual_vertex_cell.push_back(c);

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

                    if(success)
                    {
                        std::size_t dual_vertex_idx;
                        auto u = point_to_index.find(dual_vertex);
                        if (u == point_to_index.end()) 
                        {
                            dual_vertex_idx = points.size();
                            points.push_back(dual_vertex);
                            point_to_index[dual_vertex] = dual_vertex_idx;
                            dual_vertex_cell.push_back(c);
                            
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

    // generate quads
    auto generate_quad = [&](const edge_descriptor& e) 
    {

        generate_quad_dmc(e, domain, edge_to_dual_vertices, polygons);
    };
    domain.for_each_edge(generate_quad);

}


// function overload
// postprocessing on to resolve nonmanifoldness
template <typename Domain, 
          typename PolygonRange>
void dual_marching_cubes_tmc(const Domain& domain,
                             const typename Domain::Geom_traits::FT isovalue,
                             std::vector<typename Domain::Geom_traits::Point_3>& points,
                             PolygonRange& polygons,
                             bool constrain_to_cell,
                             Isosurfacing::Vertex_strategy strategy,
                             PostProcessOn,
                             gridDim gridDims)
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
    std::vector<cell_descriptor> dual_vertex_cell; // use later for track primal face for asymptotic decider


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
    // const Point_3 target(109.907, 108.431, 35.1872);
    // const Point_3 target(109.661, 110.134, 34.3311);
    // const Point_3 target(110.302, 109.476, 35.0793);
    // const Point_3 target(110.331, 112.404, 34.3367);
    const Point_3 target1(65.8811, 135.8070, 93.4010);
    const Point_3 target2(65.8587, 138.4150, 93.4168);
    // const Point_3 target1(110.331, 112.404, 34.3367);
    // const Point_3 target2(109.661, 110.134, 34.3311);

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
                            dual_vertex_cell.push_back(c);

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
                            dual_vertex_cell.push_back(c);
                            
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

    // -------------------------------------------------------------------------------------
    // below is for resolving nonmanifold edges and nonmanifold vertices
    // -------------------------------------------------------------------------------------

    // 1. resolve nonmanifold edges

    // count how many polygons share each dual edge 
    std::unordered_map<
        std::pair<std::size_t,std::size_t>, 
        int, 
        boost::hash<std::pair<std::size_t,std::size_t>>
    > dual_edge_use;

    std::unordered_map<
        std::pair<std::size_t,std::size_t>, 
        std::vector<std::size_t>,
        boost::hash<std::pair<std::size_t,std::size_t>>
    > edge_to_faces;

    // normalize orientation
    auto norm_edge = [](std::size_t a, std::size_t b) 
    {
        return (a < b) ? std::pair<std::size_t,std::size_t>(a,b)
                    : std::pair<std::size_t,std::size_t>(b,a);
    };

    std::size_t fid = 0;
    for (const auto& poly : polygons) 
    {
        const int m = static_cast<int>(poly.size());   
        for (int i = 0; i < m; ++i) 
        {
            auto e = norm_edge(poly[i], poly[(i+1)%m]); //%m to wrap around
            dual_edge_use[e] += 1;

            auto& faces = edge_to_faces[e];
            if (std::find(faces.begin(), faces.end(), fid) == faces.end()) faces.push_back(fid);
        }
        ++fid;
    }

    // get all nonmanifold edges; nme if > 2 incident faces 
    std::vector<std::pair<std::size_t, std::size_t>> nme;
    for (const auto& [e, cnt] : dual_edge_use) 
    {
        auto [a,b] = e;

        // if (cnt == 1) std::cout << "[DUAL-EDGE boundary] (" << a << "," << b << ") faces=" << cnt << "\n";
        if (cnt == 1) continue; // boundary edge
        else if (cnt >= 3) 
        {
            nme.emplace_back(a,b);
            // std::cout << "[DUAL-EDGE NONMANIFOLD] ...\n";

            // the the source of the dual vertex
            auto get_src = [&](std::size_t dual_vertex_index) 
            {
                if (auto it = dual_src.find(dual_vertex_index); it != dual_src.end()) return it->second;

                // if somehow cannot find the tag, warn and default to regular
                // std::cout << "dual vertex index" << dual_vertex_index << "does not have a tag \n";
                return DualSrc::Regular;
            };

            auto src_a = get_src(a);
            auto src_b = get_src(b);
            // std::cout << "[DUAL-EDGE NONMANIFOLD] (" << a << "," << b << ") incident_faces=" << cnt << "\n";
            // 0 = Regular; 1 = TMC_NoTunnel; 2 = TMC_Tunnel
            // std::cout << "  sources: a=" << int(src_a) << " b=" << int(src_b) << "\n";
        }
    }

    if (nme.empty()) return;

    std::unordered_map<std::size_t, FaceRec> old_fid_to_FaceRec;
    std::vector<std::vector<std::size_t>> to_add;

    // collect all unique faces incident to nme
    std::unordered_set<std::size_t> uniq;
    uniq.reserve(nme.size() * 4);
    for (auto& e: nme)
    {
        auto it = edge_to_faces.find(e);
        if (it == edge_to_faces.end()) continue;
        uniq.insert(it->second.begin(), it->second.end());
    }
    std::vector<std::size_t> face_ids(uniq.begin(), uniq.end());
    std::sort(face_ids.begin(), face_ids.end());

    subdivide_quads_skip_incident(domain, face_ids, nme, points, polygons, old_fid_to_FaceRec, to_add);

    // return all sets of points id need to be merged 
    std::vector<std::vector<std::size_t>> to_be_merged_points; 
    for (auto& e: nme)
    {
        const Point_3 p1 = points[e.first];
        const Point_3 p2 = points[e.second];
        if ((CGAL::squared_distance(p1, target1) < 1e-6 && CGAL::squared_distance(p2, target2) < 1e-6) || (CGAL::squared_distance(p1, target2) < 1e-6 && CGAL::squared_distance(p2, target1) < 1e-6))
        {
            std::cout << "Matched target nonmanifold edge \n";

        }
        centers_to_merge_for_nme(domain, e, edge_to_faces, old_fid_to_FaceRec, dual_vertex_cell, isovalue, to_be_merged_points, points);

    }

    update_quads(domain, to_be_merged_points, points, to_add, polygons);

    // 2. resolve nonmanifold vertices
    std::unordered_set<std::size_t> nmv = find_nonmanifold_vertices(points, polygons);
    if (nmv.empty()) 
    {
        std::cout << "no nonmanifold vertices found\n; exit...";
        return;
    }
    // for (auto v_id: nmv)
    // {
    //     const auto c = dual_vertex_cell[v_id];
    //     if (!is_boundary_cell(con))
    // }

}


} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif