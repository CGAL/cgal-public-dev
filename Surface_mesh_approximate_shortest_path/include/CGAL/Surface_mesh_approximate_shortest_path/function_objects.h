#ifndef CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_FUNCTION_OBJECTS_H
#define CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_FUNCTION_OBJECTS_H

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/helpers.h>

namespace CGAL {

namespace Surface_mesh_approximate_shortest_path_3 {

template <class Kernel, class Surface_mesh>
class Unfold_triangle_3_along_halfedge{
public:
    typedef typename Kernel::FT         FT;
    typedef typename Kernel::Point_2    Point_2; // the kernel is 3D - does that lead to problems?
    typedef typename Kernel::Vector_2   Vector_2;

    struct unfolded_triangle_2 {
        //Point_2 A; // note that the first point is always (0,0)
        Point_2 B;
        Point_2 P;
    };

    typedef Surface_mesh                                SM;
    typedef typename boost::graph_traits<SM>            graph_traits;
    typedef typename graph_traits::halfedge_descriptor  halfedge_descriptor;

    typedef unfolded_triangle_2 result_type;

public:
    Unfold_triangle_3_along_halfedge()
    {
    }

    //Unfold_triangle_3_along_halfedge(const Kernel& kernel)
    //    : m_

    result_type operator() (halfedge_descriptor h)
    {
        unfolded_triangle_2 unfolded_triangle;

        // first point
        // unfolded_triangle.A = {0., 0.};

        // edge length
        FT e0 = squared_distance(SM::source(h), SM::target(h));

        // second point
        unfolded_triangle.B = {sqrt(e0), 0.};

        // second edge length
        halfedge_descriptor nexth = SM::next(h);
        FT e1 = squared_distance(SM::source(nexth), SM::target(nexth));

        // third edge length
        FT e2 = squared_distance(SM::target(nexth), SM::source(h));

        // third point
        Point_2 P;
        P[0] = e0 + (e2 - e1);
        P[0] /= 2*e0;
        P[1] = sqrt(e2 - square(P[0]));
        unfolded_triangle.P = P;

        return unfolded_triangle;
    }

};

}

}

#endif // CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_FUNCTION_OBJECTS_H
