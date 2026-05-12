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

    \tparam PointInputIterator is a model of `InputIterator`.

    \tparam PointMap is a model of `ReadablePropertyMap` with value
    type `Point_3<Kernel>`.

    \tparam NormalMap is a model of `ReadablePropertyMap` with value
    type `Vector_3<Kernel>`.

    \tparam PolygonMesh a model of `MutableFaceGraph` with an internal
    point property map.


    \param begin iterator on the first point of the sequence.
    \param end past the end iterator of the point sequence.
    \param point_map property map: value_type of `InputIterator` -> Point_3.
    \param normal_map property map: value_type of `InputIterator` -> Vector_3.
    \param output_mesh where the reconstruction is stored.
    \param spacing size parameter.
    \return `true` if reconstruction succeeded, `false` otherwise.
  */
  template <typename PointInputIterator,
            typename PointMap,
            typename NormalMap,
            typename PolygonMesh>
  bool
  splat_surface_reconstruction(PointInputIterator begin,
                               PointInputIterator end,
                               PointMap point_map,
                               NormalMap normal_map,
                               PolygonMesh& output_mesh,
                               double spacing)
  {
    typedef typename boost::property_traits<PointMap>::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel Kernel;
    typedef typename Kernel::Sphere_3 Sphere;
    typedef typename Kernel::FT FT;


    return true;
  }


}


#endif // CGAL_SPLAT_SURFACE_RECONSTRUCTION_H
