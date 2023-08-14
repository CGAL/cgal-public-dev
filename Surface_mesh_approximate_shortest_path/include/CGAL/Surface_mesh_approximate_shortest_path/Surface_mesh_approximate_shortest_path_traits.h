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

namespace CGAL {

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

private:
    Kernel m_kernel;

public:
    Surface_mesh_approximate_shortest_path_traits() {}

    Surface_mesh_approximate_shortest_path_traits(const Kernel& kernel)
        : m_kernel(kernel) {}

    // visility heuristics
    class VisibilityHeuristic
    {
    public:
        VisibilityHeuristic() {};

        virtual Point_2 operator() ();
    };

};

}

#endif // CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_TRAITS_H
