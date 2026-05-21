// Copyright (c) 2026  GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pranav Jain

#ifndef CGAL_SPLAT_SURFACE_RECONSTRUCTION_H
#define CGAL_SPLAT_SURFACE_RECONSTRUCTION_H

#include <CGAL/license/Splat_surface_reconstruction_3.h>


#include <CGAL/property_map.h>
#include <CGAL/Splat_surface_reconstruction_3/internal/Box_grid.h>


namespace CGAL {

  /*!
    \ingroup PkgSplatSurfaceReconstruction3Ref

    Performs splat surface reconstruction as follows:

    - compute the ...
    - outputs the result in a polygon mesh

    This function relies mainly on the size parameter `spacing`.

    \tparam PointRange is a model of `Range`.
    \tparam NormalRange is a model of `Range`.
    \tparam PolygonMesh a model of `MutableFaceGraph` with an internal
    point property map.


    \param points the range of points.
    \param normals the range of normals.
    \param output_mesh where the reconstruction is stored.
    \param spacing size parameter.
    \return `true` if reconstruction succeeded, `false` otherwise.
  */
  template <typename PointRange,typename NormalRange,
            typename PolygonMesh>
  bool
  splat_surface_reconstruction(const PointRange& points,
                               const NormalRange& normals,
                               PolygonMesh& output_mesh,
                               double spacing)
  {
    typedef typename PointRange::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel Kernel;
    typedef typename Kernel::Sphere_3 Sphere;
    typedef typename Kernel::FT FT;


    return true;
  }


}


#endif // CGAL_SPLAT_SURFACE_RECONSTRUCTION_H
