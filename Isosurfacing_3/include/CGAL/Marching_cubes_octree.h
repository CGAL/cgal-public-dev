// Copyright (c) 2022  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ágoston Sipos

#ifndef CGAL_MARCHINGCUBESOCTREE_H
#define CGAL_MARCHINGCUBESOCTREE_H

#include <CGAL/Octree_domain.h>
#include <CGAL/Octree_mesh_extractor.h>

namespace CGAL {

template<class Domain_, class PointRange, class PolygonRange>
void make_polygon_mesh_using_marching_cubes_on_octree(const Domain_& domain, const typename Domain_::FT iso_value,
                                            PointRange& points, PolygonRange& polygons) {

    if constexpr(std::is_same_v<Domain_, CGAL::Octree_domain<typename Domain_::Geom_traits>>) {
        const typename Domain_::Octree& octree = domain.getOctree();

        CGAL::Octree_mesh_extractor<typename Domain_::Geom_traits> extractor (octree, iso_value);

        domain.iterate_voxels(extractor);

        points = extractor.get_vertices();
        polygons = extractor.get_faces();
    }
    else {
        throw CGAL::Precondition_exception("Octree_marching_cubes"
                    , "std::is_same_v<Domain_, CGAL::Octree_domain<typename Domain_::Geom_traits>>"
                    , "Marching_cubes_octree.h", 35, "Octree isosurface extraction is only available on an Octree_domain");
    }
}

}

#endif