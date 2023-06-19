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

#ifndef CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_TRAITS_H
#define CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_TRAITS_H

#include <CGAL/squared_distance_3.h>

#include <CGAL/Surface_mesh_approximate_shortest_path/function_objects.h>

#include <CGAL/boost/graph/helpers.h>
#include <boost/array.hpp>

namespace CGAL {

template <class K, class TriangleMesh>
class Surface_mesh_approximate_shortest_path_traits : public K
{
// typedefs
    public:

        typedef TriangleMesh Triangle_mesh;

        typedef K Kernel;
        typedef typename Kernel::FT FT;

        typedef typename Surface_mesh_approximate_shortest_path_3::Compute_squared_edge_length<Kernel, Triangle_mesh>                         Compute_squared_edge_length;
        typedef typename Surface_mesh_approximate_shortest_path_3::Unfold_triangle_3_along_halfedge<Kernel, Triangle_mesh>                    Unfold_triangle_3_along_halfedge;
        typedef typename Surface_mesh_approximate_shortest_path_3::Reconstruct_source_point_in_triangle_tangent_space<Kernel, Triangle_mesh>  Reconstruct_source_point_in_triangle_tangent_space;
        typedef typename Surface_mesh_approximate_shortest_path_3::Construct_triangle_centroid_2<Kernel>                                      Construct_triangle_centroid_2;
        typedef typename Surface_mesh_approximate_shortest_path_3::Construct_heuristic_point_2<Kernel, Triangle_mesh>                         Construct_heuristic_point_2;
        typedef typename Surface_mesh_approximate_shortest_path_3::Edge_intersection_test<Kernel>                                             Edge_intersection_test;

// Predicates (queries that can typically be answered by yes or no)
    // is the heuristic point visible from the virtual source?
    // is vertex a saddle vertex?

// Constructions
public:

    //class Construct_heuristic_point
    //{
        // this needs Construct_Barycenter and the LUT from the paper

        // a linear interpolation between the barycenter and the vertex
        // away from the halfedge needs to be done
    //};

    // in Surface_mesh_shortest_path, there are some typedefs here
    // they involve Surface_mesh_shortest_paths_3 which I have not seen yet

private:
    Kernel m_kernel;

public:
    Surface_mesh_approximate_shortest_path_traits() {}

    Surface_mesh_approximate_shortest_path_traits(const Kernel& kernel)
        : m_kernel(kernel) {}

    // function to give outside access by just returning an instance of the class
    Compute_squared_edge_length compute_squared_edge_length_object() const { return Compute_squared_edge_length(); };
    Unfold_triangle_3_along_halfedge unfold_triangle_3_along_halfedge_object() const { return Unfold_triangle_3_along_halfedge(); };
    Reconstruct_source_point_in_triangle_tangent_space reconstruct_source_point_in_triangle_tangent_space_object() const { return Reconstruct_source_point_in_triangle_tangent_space(); };
    Construct_triangle_centroid_2 construct_centroid_2_object() const { return Construct_triangle_centroid_2(); };
    Construct_heuristic_point_2 construct_heuristic_point_2_object() const { return Construct_heuristic_point_2(); };
    Edge_intersection_test edge_intersection_test_object() const { return Edge_intersection_test(); };
};

}

#endif // CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_TRAITS_H
