#ifndef CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_FUNCTION_OBJECTS_H
#define CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_FUNCTION_OBJECTS_H

#include "CGAL/Surface_mesh_approximate_shortest_path/Surface_mesh_approximate_shortest_path.h" // this is not supposed to be here! structure problem
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
        FT Px = (e0 + (e2 - e1)) / (2*sqrt(e0));
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
        halfedge_descriptor oppo_halfedge = h;//tmesh.opposite(h);

        // edge length
        FT e0 = Find_edge_length_and_update_property_map()(tmesh, oppo_halfedge, edge_lengths);

        // second point
        unfolded_triangle.B = Point_2(sqrt(e0), 0.);

        // second edge length
        halfedge_descriptor next_halfedge = tmesh.next(oppo_halfedge);
        FT e1 = Find_edge_length_and_update_property_map()(tmesh, next_halfedge, edge_lengths);

        // third edge length
        halfedge_descriptor prev_halfedge = tmesh.next(next_halfedge);
        FT e2 = Find_edge_length_and_update_property_map()(tmesh, prev_halfedge, edge_lengths);

        // third point
        FT Px = (e0 + (e2 - e1)) / (2*sqrt(e0));
        FT Py = sqrt(e2 - square(Px));
        unfolded_triangle.P = Point_2(Px, Py);

        //std::cout << "unfolded triangle" << std::endl;
        //std::cout << "B: " << unfolded_triangle.B.x() << "," << unfolded_triangle.B.y() << std::endl;
        //std::cout << "P: " << unfolded_triangle.P.x() << "," << unfolded_triangle.P.y() << std::endl;
        //std::cout << "with (squared) edge lengths" << std::endl;
        //std::cout << "A -> B: " << edge_lengths[tmesh.edge(oppo_halfedge)] << std::endl;
        //std::cout << "B -> P: " << edge_lengths[tmesh.edge(next_halfedge)] << std::endl;
        //std::cout << "P -> A: " << edge_lengths[tmesh.edge(prev_halfedge)] << std::endl;

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
    typedef typename CGAL::Face_values<Kernel> Face_values; // that only works because of the bad inclusion, see top of this file.
    typedef typename Surface_mesh::template Property_map<face_descriptor, Face_values> Face_values_map;

    typedef Compute_squared_edge_length<Kernel, SM> Compute_squared_edge_length;
    typedef Find_edge_length_and_update_property_map<Kernel, Surface_mesh> Find_edge_length_and_update_property_map;

    typedef Point_2 result_type;

public:
    Reconstruct_source_point_in_triangle_tangent_space() {}

    result_type operator() (SM& tmesh,
                halfedge_descriptor h,
                Edge_property_map& edge_lengths,
                Face_values_map& face_values)
    {
        // edge length
        FT e0 = Find_edge_length_and_update_property_map()(tmesh, h, edge_lengths);

        halfedge_descriptor oppo_halfedge = h;//tmesh.opposite(h);
        // find the correct entries in the face_values_map
        if (is_border(oppo_halfedge, tmesh)) {
            std::cerr << "halfedge opposite to " << h << " is on border and hence there is no way to reconstruct the source" << std::endl;
        }
        face_descriptor opposite_face = tmesh.face(oppo_halfedge);

        vertex_descriptor A = tmesh.target(oppo_halfedge); // this is swapped because target(h) == source(opposite(h))
        vertex_descriptor B = tmesh.source(oppo_halfedge);
        //std::cout << opposite_face << std::endl;
        //std::cout << A << std::endl;
        //std::cout << B << std::endl;

        int A_loc = vertex_index_in_face(A, opposite_face, tmesh);
        int B_loc = vertex_index_in_face(B, opposite_face, tmesh);

        FT d2A = face_values[opposite_face].d2verts[A_loc];
        FT d2B = face_values[opposite_face].d2verts[B_loc];

        // first coordinate of the virtual geodesic source S
        FT Sx = (e0 + (d2A - d2B)) / (2.*sqrt(e0));
        FT d2ASx = d2A - square(Sx);
        FT Sy;

        if (-1e-10 < d2ASx && d2ASx < 0.)  // there should never be negative numbers unless d2A == Sx and numerical errors hit
        {
            //d2ASx = 0.;
            Sy = 0;
        }
        else
        {
            assert(d2ASx >= 0.);
            Sy = -sqrt(d2ASx);
        }

        // Source point in triangle tangent plane
        Point_2 S = {Sx, Sy};

        //std::cout << "reconstructed source:" << std::endl;
        //std::cout << "B: " << S.x() << "," << S.y() << std::endl;

        return S;
    }
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
        //std::cout << std::fixed;
        //std::cout << std::setprecision(14);
        //std::cout << "centroid has coordinates" << std::endl;
        //std::cout << C.x() << "," << C.y() << std::endl;

        return operator() (mesh, h, unfolded_triangle.P, C, edge_lengths);
    }

    result_type operator() (SM& mesh, halfedge_descriptor h, Point_2 P, Point_2 C, Edge_property_map& edge_lengths)
    {
        // get edge lengths
        FT e0 = edge_lengths[mesh.edge(h)];
        FT e1 = edge_lengths[mesh.edge(mesh.next(h))];
        FT e2 = edge_lengths[mesh.edge(mesh.prev(h))];

        // look up the blending weight lambda
        FT lambda = Get_heuristic_parameter<Kernel>()(e0, e1, e2, P);

        // compute heuristic point coordinates
        FT Qx = lambda * C.x() + (1-lambda) * P.x();
        FT Qy = lambda * C.y() + (1-lambda) * P.y();

        //std::cout << "heuristic point has coordinates" << std::endl;
        //std::cout << Qx << "," << Qy << std::endl;

        return Point_2(Qx, Qy);
    }
};

template <class Kernel>
class Edge_intersection_test
{
public:
    typedef typename Kernel::Point_2        Point_2;

    typedef typename Kernel::Left_turn_2    Left_turn_2;

    typedef std::pair<bool, bool> result_type; // the first bool indicates whether we have an intersection,
    // the secong bool tells us whether the source is to the left (false) or to the right (true) (only relevant if the first bool is false)

public:
    Edge_intersection_test() {};

    result_type operator() (Point_2 S, Point_2 Q, Point_2 A, Point_2 B)
    {
        bool left_of_edge = Left_turn_2()(S, A, Q);
        bool not_right_of_edge = Left_turn_2()(S, B, Q);

        result_type intersection_result;

        // check all the possible cases
        if (left_of_edge && !not_right_of_edge) {
            std::cerr << "Intersection test with edge failed. The source was identified to be both left and right of the edge (which is not possible)";
        }
        else if (left_of_edge && not_right_of_edge) {
            // source is to the left of the edge
            intersection_result.first = false;
            intersection_result.second = false;
        }
        else if (!left_of_edge && !not_right_of_edge) {
            // source is to the right of the edge
            intersection_result.first = false;
            intersection_result.second = true;
        }
        else if (!left_of_edge && not_right_of_edge) {
            // source is below the edge such that we get an intersection!
            intersection_result.first = true;
            intersection_result.second = false;
        }

        return intersection_result;
    }

};

}

}

#endif // CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_FUNCTION_OBJECTS_H
