#ifndef CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_FUNCTION_OBJECTS_H
#define CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_FUNCTION_OBJECTS_H

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/helpers.h>

#include <array>

namespace CGAL {

namespace Surface_mesh_approximate_shortest_path_3 {

template <class Kernel, class Surface_mesh>
class Compute_squared_edge_length {
public:
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_3 Point_3;

    typedef Surface_mesh                                SM;
    typedef typename boost::graph_traits<SM>            graph_traits;
    typedef typename graph_traits::vertex_descriptor    vertex_descriptor;
    typedef typename graph_traits::halfedge_descriptor  halfedge_descriptor;

public:
    Compute_squared_edge_length() {};

    FT operator() (SM& tmesh, halfedge_descriptor h)
    {
        vertex_descriptor v1 = tmesh.source(h);
        vertex_descriptor v2 = tmesh.target(h);

        return operator() (tmesh, v1, v2);
    }

    FT operator() (SM& tmesh, vertex_descriptor v1, vertex_descriptor v2)
    {
        Point_3 v1_point = tmesh.point(v1);
        Point_3 v2_point = tmesh.point(v2);

        double length = squared_distance(v1_point, v2_point);
        return length;
    }
};

template <class Kernel, class Surface_mesh>
class Find_edge_length_and_update_property_map
{
public:
    typedef typename Kernel::FT FT;

    typedef Surface_mesh                    SM;
    typedef boost::graph_traits<SM>         Graph_traits;
    typedef typename Graph_traits::edge_descriptor edge_descriptor;
    typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;

    typedef typename SM::template Property_map<edge_descriptor, FT> Edge_property_map;

    typedef Compute_squared_edge_length<Kernel, SM> Compute_squared_edge_length;

    typedef FT result_type;

public:
    Find_edge_length_and_update_property_map() {};

    result_type operator() (SM& mesh, halfedge_descriptor h)
    {
        return Compute_squared_edge_length()(mesh, h);
    }

    result_type operator() (SM& mesh, halfedge_descriptor h, Edge_property_map& edge_lengths)
    {
        edge_descriptor e = mesh.edge(h);
        FT edge_length = edge_lengths[e];

        if (edge_length == 0)
        {
            edge_length = Compute_squared_edge_length()(mesh, h);
            edge_lengths[e] = edge_length;
        }

        return edge_length;
    }
};

template <class Kernel, class Surface_mesh>
class Unfold_triangle_3_along_halfedge{
public:
    typedef typename Kernel::FT           FT;
    typedef typename Kernel::Point_2      Point_2;
    typedef typename Kernel::Point_3      Point_3;

    struct unfolded_triangle_2 {
        //Point_2 A; // note that the first point is always (0,0)
        Point_2 B;
        Point_2 P;
    };

    typedef Surface_mesh                                SM;
    typedef typename boost::graph_traits<SM>            graph_traits;
    typedef typename graph_traits::edge_descriptor      edge_descriptor;
    typedef typename graph_traits::halfedge_descriptor  halfedge_descriptor;

    typedef typename SM::template Property_map<edge_descriptor, FT> Edge_property_map;

    typedef Compute_squared_edge_length<Kernel, SM> Compute_squared_edge_length;
    typedef Find_edge_length_and_update_property_map<Kernel, SM> Find_edge_length_and_update_property_map;

    typedef unfolded_triangle_2 result_type;

public:
    Unfold_triangle_3_along_halfedge() {}

    result_type operator() (SM& tmesh, halfedge_descriptor h)
    {
        unfolded_triangle_2 unfolded_triangle;

        // edge length
        FT e0 = Compute_squared_edge_length()(tmesh, h);

        // second point
        unfolded_triangle.B = Point_2(sqrt(e0), 0.);

        // second edge length
        halfedge_descriptor nexth = tmesh.next(h);
        FT e1 = Compute_squared_edge_length()(tmesh, nexth);

        // third edge length
        FT e2 = Compute_squared_edge_length()(
            tmesh, tmesh.target(nexth), tmesh.source(h));

        // third point
        FT Px = (e0 + (e2 - e1)) / (2 * e0);
        FT Py = sqrt(e2 - square(Px));
        unfolded_triangle.P = Point_2(Px, Py);

        std::cout << "unfolded triangle" << std::endl;
        std::cout << "B: " << unfolded_triangle.B.x() << "," << unfolded_triangle.B.y() << std::endl;
        std::cout << "P: " << unfolded_triangle.P.x() << "," << unfolded_triangle.P.y() << std::endl;

        return unfolded_triangle;
    }

    result_type operator() (SM& tmesh, halfedge_descriptor h, Edge_property_map& edge_lengths)
    {
        unfolded_triangle_2 unfolded_triangle;

        // edge length
        FT e0 = Find_edge_length_and_update_property_map()(tmesh, h, edge_lengths);

        // second point
        unfolded_triangle.B = Point_2(sqrt(e0), 0.);

        // second edge length
        halfedge_descriptor nexth = tmesh.next(h);
        FT e1 = Find_edge_length_and_update_property_map()(tmesh, nexth, edge_lengths);

        // third edge length
        halfedge_descriptor prevh = tmesh.next(nexth);
        FT e2 = Find_edge_length_and_update_property_map()(tmesh, prevh, edge_lengths);

        // third point
        FT Px = (e0 + (e2 - e1)) / (2 * e0);
        FT Py = sqrt(e2 - square(Px));
        unfolded_triangle.P = Point_2(Px, Py);

        std::cout << "unfolded triangle" << std::endl;
        std::cout << "B: " << unfolded_triangle.B.x() << "," << unfolded_triangle.B.y() << std::endl;
        std::cout << "P: " << unfolded_triangle.P.x() << "," << unfolded_triangle.P.y() << std::endl;
        std::cout << "with (squared) edge lengths" << std::endl;
        std::cout << "A -> B: " << edge_lengths[tmesh.edge(h)] << std::endl;
        std::cout << "B -> P: " << edge_lengths[tmesh.edge(nexth)] << std::endl;
        std::cout << "P -> A: " << edge_lengths[tmesh.edge(prevh)] << std::endl;

        return unfolded_triangle;
    }
};

template <class Kernel, class Surface_mesh>
class Reconstruct_source_point_in_triangle_tangent_space{
public:
    typedef typename Kernel::FT           FT;
    typedef typename Kernel::Point_2      Point_2;
    typedef typename Kernel::Point_3      Point_3;

    typedef Surface_mesh                                SM;
    typedef typename boost::graph_traits<SM>            Graph_traits;
    typedef typename Graph_traits::vertex_descriptor    vertex_descriptor;
    typedef typename Graph_traits::edge_descriptor      edge_descriptor;
    typedef typename Graph_traits::halfedge_descriptor  halfedge_descriptor;
    typedef typename Graph_traits::face_descriptor      face_descriptor;

    typedef typename SM::template Property_map<edge_descriptor, FT> Edge_property_map;

    typedef Compute_squared_edge_length<Kernel, SM> Compute_squared_edge_length;
    typedef Find_edge_length_and_update_property_map<Kernel, Surface_mesh> Find_edge_length_and_update_property_map;

    typedef Point_2 result_type;

public:
    Reconstruct_source_point_in_triangle_tangent_space() {}

    result_type operator() (SM& tmesh, halfedge_descriptor h)
    {
        // edge length
        FT e0 = Compute_squared_edge_length()(tmesh, h);

        // just for now until the proper data structures are implemented:
        FT d2_AU = 1.; // squared distances
        FT d2_BU = 1.;

        // first coordinate of the virtual geodesic source S
        FT Sx = (e0 + (d2_AU - d2_BU)) / (2.*e0);
        FT Sy = -sqrt(d2_AU - square(Sx));

        // Source point in triangle tangent plane
        Point_2 S = {Sx, Sy};

        std::cout << "reconstructed source" << std::endl;
        std::cout << "S: " << S.x() << "," << S.y() << std::endl;

        return S;
    }

    /* result_type operator() (SM& tmesh,
                halfedge_descriptor h,
                Edge_property_map& edge_lengths,
                Face_values_map& face_values)
    {
        // edge length
        FT e0 = Find_edge_length_and_update_property_map()(tmesh, h, edge_lengths);

        // find the correct entries in the face_values_map
        face_descriptor f = tmesh.face(h);
        FT d2_AU = face_values[f].d2v0; // squared distances
        FT d2_BU = face_values[f].d2v1; // this is NOT yet correct. One needs to find a global to
        // local relationship for the vertices in the respective triangle here.
        // it always needs to be the distances corresponding to the vertices source(h) and target(h)
        // but how do I respect that in the face_values and then here?
        // Is there something like a local vertex_property_map for the local vertices?

        // first coordinate of the virtual geodesic source S
        FT Sx = (e0 + (d2_AU - d2_BU)) / (2.*e0);
        FT Sy = -sqrt(d2_AU - square(Sx));

        // Source point in triangle tangent plane
        Point_2 S = {Sx, Sy};

        std::cout << "reconstructed source:" << std::endl;
        std::cout << "B: " << S.x() << "," << S.y() << std::endl;

        return S;
    } */
};

template <class Kernel>
class Construct_triangle_centroid_2
{
public:
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_2 Point_2;

    typedef Point_2 result_type;

public:
    Construct_triangle_centroid_2() {}

    result_type operator() (Point_2 B, Point_2 P)
    {
        FT cx = (B.x() + P.x()) / 3.;
        FT cy = P.y() / 3;

        return Point_2(cx, cy);
    }
};

template <class Kernel>
class Get_heuristic_parameter
{
public:
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_2 Point_2;
    typedef FT result_type;

public:
    Get_heuristic_parameter() {};

    result_type operator() (
        FT e0, // e0 is always the length of the halfedge along which we unfolded
        FT e1,
        FT e2,
        Point_2 P)
    {
        FT h = P.y();

        // get longest and shortest edges
        auto max_e = std::max(e0, std::max(e1, e2));
        auto min_e = std::min(e0, std::min(e1, e2));

        // compute the 4 parameters listed in 3.2 in the paper
        //FT tau1 = max_e / min_e;
        //FT tau2 = max_e / e0;
        //FT tau3 = h / max_e;
        //FT tau4 = h / e0;

        // get heuristic blending parameter
        FT threshold_1 = 5.1424f;
        FT threshold_2 = 4.20638f;
        FT threshold_3 = 0.504201f;
        FT threshold_4 = 2.84918f;
        std::array<FT,16> lambda = {0.320991f, 0.446887f, 0.595879f,  0.270094f,  0.236679f, 0.159685f,  0.0872932f, 0.434132f,
                                    1.0f,      0.726262f, 0.0635997f, 0.0515979f, 0.56903f,  0.0447586f, 0.0612103f, 0.718198f};

        auto b0 = max_e > threshold_1 * e0;
        auto b1 = max_e > threshold_2 * min_e;
        auto b2 = h < threshold_3 * e0;     // the publication and the supplemental implementation
        auto b3 = h < threshold_4 * max_e;  // enumerate the parameters differently
        int idx = b0 + b1 * 2 + b2 * 4 + b3 * 8;
        FT l = lambda[idx];

        return l;
    }
};

template <class Kernel, class Surface_mesh>
class Construct_heuristic_point_2
{
public:
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_2 Point_2;

    typedef Surface_mesh                                SM;
    typedef typename boost::graph_traits<SM>            Graph_traits;
    typedef typename Graph_traits::edge_descriptor      edge_descriptor;
    typedef typename Graph_traits::halfedge_descriptor  halfedge_descriptor;

    typedef typename SM::template Property_map<edge_descriptor, FT> Edge_property_map;

    typedef Find_edge_length_and_update_property_map<Kernel, Surface_mesh> Find_edge_length_and_update_property_map;
    typedef Unfold_triangle_3_along_halfedge<Kernel, Surface_mesh> Unfold_triangle_3_along_halfedge;

    typedef Point_2 result_type;

public:
    Construct_heuristic_point_2() {}

    result_type operator() (SM& mesh, halfedge_descriptor h, Edge_property_map& edge_lengths)
    {
        // get barycenter
        auto unfolded_triangle = Unfold_triangle_3_along_halfedge()(mesh, h, edge_lengths);
        Point_2 C = Construct_triangle_centroid_2<Kernel>()(unfolded_triangle.B, unfolded_triangle.P);
        std::cout << "centroid has coordinates" << std::endl;
        std::cout << C.x() << "," << C.y() << std::endl;

        // get edge lengths
        FT e0 = edge_lengths[mesh.edge(h)];
        FT e1 = edge_lengths[mesh.edge(mesh.next(h))];
        FT e2 = edge_lengths[mesh.edge(mesh.prev(h))];

        // look up to blending weight lambda
        FT lambda = Get_heuristic_parameter<Kernel>()(e0, e1, e2, unfolded_triangle.P);

        // compute heuristic point coordinates
        FT Qx = lambda * C.x() + (1-lambda) * unfolded_triangle.P.x();
        FT Qy = lambda * C.y() + (1-lambda) * unfolded_triangle.P.y();

        std::cout << "heuristic point has coordinates" << std::endl;
        std::cout << Qx << "," << Qy << std::endl;

        return Point_2(Qx, Qy);
    }
};

}

}

#endif // CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_FUNCTION_OBJECTS_H
