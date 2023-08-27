// Copyright (c) 2023 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Robert Piel

//#ifndef CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_TRAITS_H
//#define CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_TRAITS_H

#include <CGAL/squared_distance_3.h>

#include <CGAL/boost/graph/helpers.h>
#include <boost/array.hpp>

namespace CGAL {

/*
// forward declarations
class Kernel;
class Surface_mesh;

template<class Kernel, class Surface_mesh>
class Surface_mesh_approximate_shortest_path_traits
{
    class Visibility_heuristic;
};
typedef Surface_mesh_approximate_shortest_path_traits<Kernel, Surface_mesh> Traits;

class Skip_condition;

class Enqueue_policy;

template<class Traits,
         class Visibility_heuristic,
         class Skip_condition,
         class Enqueue_policy>
class Surface_mesh_approximate_shortest_path;
typedef Surface_mesh_approximate_shortest_path<Traits, Visibility_heuristic, Skip_condition, Enqueue_policy> Approximate_shortest_path;
*/

template <class K, class SurfaceMesh>
class Surface_mesh_approximate_shortest_path_traits : public K
{
public:

    typedef SurfaceMesh Surface_mesh;

    typedef K Kernel;
    typedef typename Kernel::FT FT;

    typedef typename Kernel::Point_2 Point_2;

    typedef boost::graph_traits<Surface_mesh> Graph_traits;
    typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
    typedef typename Graph_traits::edge_descriptor edge_descriptor;

    typedef typename Surface_mesh::template Property_map<edge_descriptor, FT> Edge_property_map;

private:
    Kernel m_kernel;

public:
    Surface_mesh_approximate_shortest_path_traits() {}

    Surface_mesh_approximate_shortest_path_traits(const Kernel& kernel)
        : m_kernel(kernel) {}

    // visility heuristics
    class Visibility_heuristic
    {
    public:
        Visibility_heuristic() {};

        Point_2 operator() (Surface_mesh& mesh, Edge_property_map& edge_lengths, halfedge_descriptor h, Point_2 P, Point_2 C)
        {
            // get edge lengths
            FT e0 = edge_lengths[mesh.edge(h)];
            FT e1 = edge_lengths[mesh.edge(mesh.next(h))];
            FT e2 = edge_lengths[mesh.edge(mesh.prev(h))];
            FT height = P.y();

            // look up the blending weight lambda
            FT lambda = get_heuristic_parameter(e0, e1, e2, height);

            // compute heuristic point coordinates
            FT Qx = lambda * C.x() + (1-lambda) * P.x();
            FT Qy = lambda * C.y() + (1-lambda) * P.y();

            return Point_2(Qx, Qy);
        }

        FT get_heuristic_parameter(FT e0, FT e1, FT e2, FT h)
        {
            // get longest and shortest edges
            auto max_e = std::max(e0, std::max(e1, e2));
            auto min_e = std::min(e0, std::min(e1, e2));

            // get heuristic blending parameter
            FT threshold_1 = 5.1424f;
            FT threshold_2 = 4.20638f;
            FT threshold_3 = 0.504201f;
            FT threshold_4 = 2.84918f;
            std::array<FT,16> lambda = {0.320991f, 0.446887f, 0.595879f,  0.270094f,  0.236679f, 0.159685f,  0.0872932f, 0.434132f,
                                         1.0f,      0.726262f, 0.0635997f, 0.0515979f, 0.56903f,  0.0447586f, 0.0612103f, 0.718198f};

            auto b0 = max_e > threshold_1 * e0;
            auto b1 = max_e > threshold_2 * min_e;
            auto b2 = h < threshold_3 * e0;     // the publication and the supplemental implementation
            auto b3 = h < threshold_4 * max_e;  // enumerate the parameters differently
            int idx = b0 + b1 * 2 + b2 * 4 + b3 * 8;
            FT l = lambda[idx];

            return l;
        }
    };

};

}

//#endif // CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_TRAITS_H
