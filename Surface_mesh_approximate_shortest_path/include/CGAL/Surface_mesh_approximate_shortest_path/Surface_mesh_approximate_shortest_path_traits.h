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

#include <CGAL/boost/graph/helpers.h>
#include <boost/array.hpp>

template <
    class K,
    class TriangleMesh >
class Surface_mesh_approximate_shortest_path_traits : public K
{
// typedefs
    public:

        typedef K Kernel;

        typedef TriangleMesh Triangle_mesh;

        typedef boost::graph_traits<Triangle_mesh> Graph_traits;
        typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
        typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
        typedef typename Graph_traits::face_descriptor face_descriptor;

        typedef typename Kernel::FT FT; // FT is the FieldNumberType (Field as in algebraic field)
        // barycentric coordinates type
        typedef typename std::array<FT, 3> barycentric_coordinates;

        //=== result types of triangle unfolding ===
        // 1. triangle_3_edge_2_parameterization stores the three (helf)edge lengths.
        // Let (e_0, e_1, e_2) be the three edge lengths
        // corresponding to the halfedge h, next(h) and next(next(h)).
        typedef typename std::array<FT, 3> triangle_3_unfolded_edge_2_parameterization;
        // 2. point_3_unfolded_tangent_space_parameterization stores the (x,y) coordinates
        // of points in the tangent plane of the unfolded triangle and the source point coordinates.
        typedef typename std::array<FT, 2> point_3_unfolded_tangent_space_parameterization;

        struct unfolded_triangle_in_tangent_plane {
            point_3_unfolded_tangent_space_parameterization A;
            point_3_unfolded_tangent_space_parameterization B;
            point_3_unfolded_tangent_space_parameterization P;
        };


// Predicates (queries that can typically be answered by yes or no)
    // is the heuristic point visible from the virtual source?
    // is vertex a saddle vertex?
    //typedef typename Surface_mesh_approximate_shortest_paths_3

// Constructions
public:
    class Unfold_along_halfedge
    {
        // this class unfolds two adjacent 3D-triangles into a 2D domain while preserving their edge lengths
    public:
        typedef unfolded_triangle_in_tangent_plane result_type;

        // compute the triangle points, edge lengths and coordinates of virtual source (naming according to the paper)
        result_type operator() (halfedge_descriptor h) const
        {
            unfolded_triangle_in_tangent_plane unfolded_triangle;

            // first point
            unfolded_triangle.A = {0., 0.};

            // length of first edge
            triangle_3_unfolded_edge_2_parameterization squared_edge_lenghts;
            squared_edge_lenghts[0] = CGAL::squared_distance(source(h), target(h));

            // second point
            unfolded_triangle.B = {CGAL::sqrt(squared_edge_lenghts[0]), 0.};

            // second edge length
            halfedge_descriptor nexth = next(h);
            squared_edge_lenghts[1] = CGAL::squared_distance(source(nexth), target(nexth));

            // third edge length
            halfedge_descriptor nextnexth = next(nexth);
            squared_edge_lenghts[2] = CGAL::squared_distance(source(nextnexth), target(nextnexth));

            // third point
            point_3_unfolded_tangent_space_parameterization P;
            P[0] = squared_edge_lenghts[0] + (squared_edge_lenghts[2] - squared_edge_lenghts[1]);
            P[0] /= 2*squared_edge_lenghts[0];
            P[1] = CGAL::sqrt(squared_edge_lenghts[2] - CGAL::square(P[0]));
            unfolded_triangle.P = P;

            return unfolded_triangle;
        }

    };

    class Construct_barycenter
    {
        // this should be implemented in the kernel already
    };

    class Construct_heuristic_point
    {
        // this needs Construct_Barycenter and the LUT from the paper

        // a linear interpolation between the barycenter and the vertex
        // away from the halfedge needs to be done
    };

    // in Surface_mesh_shortest_path, there are some typedefs here
    // they involve Surface_mesh_shortest_paths_3 which I have not seen yet

private:
    Kernel m_kernel;
    Unfold_along_halfedge m_unfold_along_halfedge_object;
    Construct_barycenter m_construct_barycenter_object;
    Construct_heuristic_point m_construct_heuristic_point_object;

public:
    Surface_mesh_approximate_shortest_path_traits() {}

    Surface_mesh_approximate_shortest_path_traits(const Kernel& kernel)
        : m_kernel(kernel)
        , m_unfold_along_halfedge_object()
    {}

};

#endif // CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_TRAITS_H
