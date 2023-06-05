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

#ifndef CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_H
#define CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_H

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>

template<class Traits>
class Surface_mesh_approximate_shortest_path
{
public:
    typedef typename Traits::Triangle_mesh Triangle_mesh;
    typedef boost::graph_traits<Triangle_mesh> Graph_traits;
    typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
    typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
    typedef typename Graph_traits::face_descriptor face_descriptor;

    // unfolding is a construction that should be done in the traits class
    typedef typename Traits::Unfold_along_halfedge Unfold_along_halfedge;

private:
    const Triangle_mesh& mesh;

public:
    Surface_mesh_approximate_shortest_path(const Triangle_mesh& mesh,
                               const Traits& traits = Traits())
        : mesh(mesh) {};

};

#endif // CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_H
