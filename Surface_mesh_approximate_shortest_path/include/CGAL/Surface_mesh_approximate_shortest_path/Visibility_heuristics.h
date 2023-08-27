#ifndef CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_VISIBILITY_HEURISTICS_H
#define CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_VISIBILITY_HEURISTICS_H

#include <algorithm>

class Visibility_heuristics
{
public:
    Visibility_heuristics() {};

    double get_heuristic_parameter(double e0, double e1, double e2, double h)
    {
        // get longest and shortest edges
        auto max_e = std::max(e0, std::max(e1, e2));
        auto min_e = std::min(e0, std::min(e1, e2));

        // get heuristic blending parameter
        float threshold_1 = 5.1424f;
        float threshold_2 = 4.20638f;
        float threshold_3 = 0.504201f;
        float threshold_4 = 2.84918f;
        std::array<float,16> lambda = {0.320991f, 0.446887f, 0.595879f,  0.270094f,  0.236679f, 0.159685f,  0.0872932f, 0.434132f,
                                        1.0f,      0.726262f, 0.0635997f, 0.0515979f, 0.56903f,  0.0447586f, 0.0612103f, 0.718198f};

        auto b0 = max_e > threshold_1 * e0;
        auto b1 = max_e > threshold_2 * min_e;
        auto b2 = h < threshold_3 * e0;     // the publication and the supplemental implementation
        auto b3 = h < threshold_4 * max_e;  // enumerate the parameters differently
        int idx = b0 + b1 * 2 + b2 * 4 + b3 * 8;
        double l = lambda[idx];

        return l;
    }

    double operator() (Surface_mesh_approximate_shortest_path& shopa, halfedge_descriptor h, Point_2 P, Point_2 C)
    {
        // get edge lengths
        Surface_mesh& mesh = shopa.m_mesh;
        auto edge_lengths = shopa.Get_edge_length_map();
        FT e0 = edge_lengths[mesh.edge(h)];
        FT e1 = edge_lengths[mesh.edge(mesh.next(h))];
        FT e2 = edge_lengths[mesh.edge(mesh.prev(h))];
        FT height = P.y();

        // look up the blending weight lambda
        FT lambda = shopa.m_Visibility_heuristic(e0, e1, e2, height);

        // compute heuristic point coordinates
        FT Qx = lambda * C.x() + (1-lambda) * P.x();
        FT Qy = lambda * C.y() + (1-lambda) * P.y();

        return Point_2(Qx, Qy);
    }
};

#endif CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_VISIBILITY_HEURISTICS_H
