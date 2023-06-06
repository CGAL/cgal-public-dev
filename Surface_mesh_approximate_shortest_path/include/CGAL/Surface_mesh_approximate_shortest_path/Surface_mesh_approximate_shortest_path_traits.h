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
        //typedef typename Kernel::ConstructBarycenter_3 CounstructBarycenter_3;

        typedef typename Surface_mesh_approximate_shortest_path_3::Unfold_triangle_3_along_halfedge<Kernel, Triangle_mesh> Unfold_triangle_3_along_halfedge;


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
    Unfold_triangle_3_along_halfedge m_unfold_triangle_3_along_halfedge_object() const { return m_unfold_triangle_3_along_halfedge_object; }

};

}

#endif // CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_TRAITS_H
