// Copyright (c) 2017  GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot, Pierre Alliez, Tong Zhao, Hongyi Liu

#ifndef CGAL_POISSON_SURFACE_RECONSTRUCTION_H
#define CGAL_POISSON_SURFACE_RECONSTRUCTION_H

#include <CGAL/license/Implicit_surface_reconstruction_3.h>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Implicit_reconstruction_function.h>
#include <CGAL/property_map.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Implicit_contouring.h>

/*!
\file Poisson_surface_reconstruction.h
*/

namespace CGAL {


    /*!
    \ingroup PkgImplicitSurfaceReconstruction3Ref

    Performs surface reconstruction as follows:

    - compute the Poisson implicit function, through a conjugate
    gradient solver, represented as a piecewise linear function
    stored on a 3D Delaunay mesh generated via octree based discretization 
    or Delaunay refinement
    - meshes the function with a user-defined precision using another
    round of Delaunay refinement: it contours the isosurface
    corresponding to the isovalue of the median of the function
    values at the input points
    - outputs the result in a polygon mesh

    
    If `use_octree` is true, the function creates an octree based discretization with predefined 
    depth, where the total size of the discretization vertices is less sensitive to input size.
    If octree based discretization is not used, Delaunay refinement will be used, where the total size of the 
    discretization vertices in general grows with input point size.

    If `use_marching_tets` is true, the function will create surface based on the same discretization
    used for solver. The isocontouring process with marching tetrehedra is faster but result in lower 
    quality mesh. If marching tetrehedron is not used, this function use implicit surface generation, 
    which relies mainly on the size parameter `average_spacing_ratio`. The function computes the 
    spacing using the multiple of `compute_average_spacing()`. Smaller values increase the precision 
    of the output mesh at the cost of higher computation time.

    Parameters `sm_angle`, `sm_radius` and `sm_distance` work
    similarly to the parameters of `SurfaceMeshFacetsCriteria_3`. The
    latest two are defined with respect to spacing computed by `compute_average_spacing()` 
    and `average_spacing_ratio`.

    \tparam PointInputIterator is a model of `InputIterator`.

    \tparam PointMap is a model of `ReadablePropertyMap` with value
    type `Point_3<Kernel>`.

    \tparam NormalMap is a model of `ReadablePropertyMap` with value
    type `Vector_3<Kernel>`.

    \tparam PolygonMesh a model of `MutableFaceGraph` with an internal
    point property map.

    \tparam Tag is a tag whose type affects the behavior of the
    meshing algorithm (see `make_surface_mesh()`).

    \param PointList the container of input points
    \param point_map property map: value_type of `InputIterator` -> Point_3.
    \param normal_map property map: value_type of `InputIterator` -> Vector_3.
    \param output_mesh where the reconstruction is stored.
    \param fitting data fitting term weight. If the result is over-smoothed, increasing this weighting generally produce a more detailed result.
    \param use_octree use octree based discretization if true, else use Delaunay refinement 
    \param use_marching_tets use marching tetrehedra on the solver's discretization for fast mesh generation, else use implicit surface generation
    \param average_spacing_ratio size ratio parameter. Ignored if use marching tetrehedra
    \param sm_angle bound for the minimum facet angle in degrees. Ignored if use marching tetrehedra
    \param sm_radius bound for the radius of the surface Delaunay balls (relatively to the `average_spacing`). Ignored if use marching tetrehedra
    \param sm_distance bound for the center-center distances (relatively to the `average_spacing`). Ignored if use marching tetrehedra
    \param tag surface mesher tag.
    \return `true` if reconstruction succeeded, `false` otherwise.
    */
    template <typename PointList,
        typename PointMap,
        typename NormalMap,
        typename PolygonMesh,
        typename Tag = CGAL::Manifold_with_boundary_tag>
        bool
        poisson_surface_reconstruction_delaunay (PointList points,
            PointMap point_map,
            NormalMap normal_map,
            PolygonMesh& output_mesh,
            double fitting = 1,
            bool use_octree = true,
            bool use_marching_tets = false,
            double average_spacing_ratio = 6,
            double sm_angle = 20.0,
            double sm_radius = 100.0,
            double sm_distance = 0.025,
            Tag tag = Tag())
    {
        typedef typename boost::property_traits<PointMap>::value_type Point;
        typedef typename Kernel_traits<Point>::Kernel Kernel;
        typedef typename Kernel::Sphere_3 Sphere;
        typedef typename Kernel::FT FT;

        typedef CGAL::Implicit_reconstruction_function<Kernel, PointList, NormalMap> Implicit_reconstruction_function;
        typedef typename CGAL::Surface_mesher::Surface_mesh_default_triangulation_3_generator<Kernel>::Type STr;
        typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
        typedef CGAL::Implicit_surface_3<Kernel, Implicit_reconstruction_function> Surface_3;

        Implicit_reconstruction_function function;
        function.initialize_point_map(points, point_map, normal_map, use_octree);

        if ( ! function.compute_poisson_implicit_function_new(fitting) )
            return false;

        if (! implicit_contouring(points, 
            point_map, 
            function, 
            output_mesh, 
            use_marching_tets, 
            average_spacing_ratio,
            sm_angle,
            sm_radius,
            sm_distance,
            tag))
            return false;

        return true;
    }

}


#endif // CGAL_POISSON_SURFACE_RECONSTRUCTION_H
