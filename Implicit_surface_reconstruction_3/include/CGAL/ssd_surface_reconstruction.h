
#ifndef CGAL_SSD_SURFACE_RECONSTRUCTION_H
#define CGAL_SSD_SURFACE_RECONSTRUCTION_H

#include <CGAL/license/Implicit_surface_reconstruction_3.h>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Implicit_reconstruction_function.h>
#include <CGAL/property_map.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Implicit_contouring.h>

namespace CGAL {


    /*!
    \ingroup PkgPoissonSurfaceReconstruction3Ref

    Performs surface reconstruction as follows:

    - compute the Poisson implicit function, through a conjugate
    gradient solver, represented as a piecewise linear function
    stored on a 3D Delaunay mesh generated via Delaunay refinement
    - meshes the function with a user-defined precision using another
    round of Delaunay refinement: it contours the isosurface
    corresponding to the isovalue of the median of the function
    values at the input points
    - outputs the result in a polygon mesh

    This function relies mainly on the size parameter `spacing`. A
    reasonable solution is to use the average spacing of the input
    point set (using `compute_average_spacing()` for example). Smaller
    values increase the precision of the output mesh at the cost of
    higher computation time.

    Parameters `sm_angle`, `sm_radius` and `sm_distance` work
    similarly to the parameters of `SurfaceMeshFacetsCriteria_3`. The
    latest two are defined with respect to `spacing`.

    \tparam PointInputIterator is a model of `InputIterator`.

    \tparam PointMap is a model of `ReadablePropertyMap` with value
    type `Point_3<Kernel>`.

    \tparam NormalMap is a model of `ReadablePropertyMap` with value
    type `Vector_3<Kernel>`.

    \tparam PolygonMesh a model of `MutableFaceGraph` with an internal
    point property map.

    \tparam Tag is a tag whose type affects the behavior of the
    meshing algorithm (see `make_surface_mesh()`).

    \param begin iterator on the first point of the sequence.
    \param end past the end iterator of the point sequence.
    \param point_map property map: value_type of `InputIterator` -> Point_3.
    \param normal_map property map: value_type of `InputIterator` -> Vector_3.
    \param output_mesh where the reconstruction is stored.
    \param spacing size parameter.
    \param sm_angle bound for the minimum facet angle in degrees.
    \param sm_radius bound for the radius of the surface Delaunay balls (relatively to the `average_spacing`).
    \param sm_distance bound for the center-center distances (relatively to the `average_spacing`).
    \param tag surface mesher tag.
    \return `true` if reconstruction succeeded, `false` otherwise.
    */
    template <typename PointList,
        typename PointMap,
        typename NormalMap,
        typename PolygonMesh,
        typename Tag = CGAL::Manifold_with_boundary_tag>
        bool
        ssd_surface_reconstruction_delaunay (PointList points,
            PointMap point_map,
            NormalMap normal_map,
            PolygonMesh& output_mesh,
            double fitting = 10,
            double laplacian = 1,
            double hessian = 1e-3,
            bool use_octree = true,
            bool use_marching_tets = false,
            double spacing_ratio = 6,
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

        if ( ! function.compute_ssd_implicit_function(fitting, laplacian, hessian) )
            return false;

        if (! implicit_contouring(points, 
            point_map, 
            function, 
            output_mesh, 
            use_marching_tets, 
            spacing_ratio,
            sm_angle,
            sm_radius,
            sm_distance,
            tag))
            return false;

        return true;
    }

}


#endif // CGAL_POISSON_SURFACE_RECONSTRUCTION_H
