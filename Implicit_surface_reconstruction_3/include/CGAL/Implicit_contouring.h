#ifndef CGAL_IMPLICIT_CONTOURING_H
#define CGAL_IMPLICIT_CONTOURING_H

namespace CGAL {

    template <
        typename PointList,
        typename PointMap,
        typename Implicit_reconstruction_function,
        typename Polyhedron,
        typename Tag = CGAL::Manifold_with_boundary_tag
    >
    bool
	implicit_contouring(
        PointList& points,
        PointMap point_map,
        Implicit_reconstruction_function & function,
        Polyhedron & output_mesh,
        bool use_marching_tets,
        double spacing_ratio,
        double sm_angle,
        double sm_radius,
        double sm_distance,
        Tag tag = Tag())
	{
        typedef typename Kernel::FT FT;
        typedef typename Kernel::Sphere_3 Sphere;
        typedef typename boost::property_traits<PointMap>::value_type Point;

        typedef typename CGAL::Surface_mesher::Surface_mesh_default_triangulation_3_generator<Kernel>::Type STr;
        typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
        typedef CGAL::Implicit_surface_3<Kernel, Implicit_reconstruction_function> Surface_3;

        if (use_marching_tets)
            return function.marching_tetrahedra(0, output_mesh); // default isovalue is 0


        FT spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(points, spacing_ratio, CGAL::parameters::point_map(point_map));
        Point inner_point = function.get_inner_point();
        Sphere bsphere = function.bounding_sphere();
        FT radius = CGAL::approximate_sqrt(bsphere.squared_radius());

        FT sm_sphere_radius = 5.0 * radius;
        FT sm_dichotomy_error = sm_distance * spacing / 1000.0;

        Surface_3 surface(function,
            Sphere (inner_point, sm_sphere_radius * sm_sphere_radius),
            sm_dichotomy_error / sm_sphere_radius);

        CGAL::Surface_mesh_default_criteria_3<STr> criteria (sm_angle,
            sm_radius * spacing,
            sm_distance * spacing);

        STr tr;
        C2t3 c2t3(tr);

        CGAL::make_surface_mesh(c2t3,
            surface,
            criteria,
            tag);

        if(tr.number_of_vertices() == 0)
            return false;

        CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, output_mesh);

        return true;
	}

}

#endif // CGAL_IMPLICIT_CONTOURING_H