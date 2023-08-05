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

#include <CGAL/Polygon_mesh_processing/locate.h>

namespace CGAL {

template<class Kernel>
struct Face_values {
    typedef typename Kernel::FT FT;

    FT sigma;
    FT d;
    std::array<FT,3> d2verts;

    Face_values(FT _sigma=0., FT _d=std::numeric_limits<FT>::max(), std::array<FT,3> _d2verts = {-1., -1., -1.})
        : sigma(_sigma), d(_d), d2verts(_d2verts) {}

    friend std::ostream & operator <<(std::ostream& stream, const Face_values vals)
    {
        return ( stream << vals.sigma << "\t" << vals.d << "\t"
                       << vals.d2verts[0]  << "\t" << vals.d2verts[1] << "\t" << vals.d2verts[2] << std::endl);
    };
};

template<class Kernel>
struct Unfolded_triangle_2 {
    typedef typename Kernel::Point_2 Point_2;

    //Point_2 A; // note that the first point is always (0,0)
    Point_2 B;
    Point_2 P;
};

template<class Traits>
class Surface_mesh_approximate_shortest_path
{
public:
    typedef typename Traits::Kernel Kernel;
    typedef typename Kernel::FT FT;

    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Point_3 Point_3;

    typedef typename Traits::Surface_mesh Surface_mesh;
    typedef boost::graph_traits<Surface_mesh> Graph_traits;

    typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
    typedef typename Graph_traits::edge_descriptor edge_descriptor;
    typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
    typedef typename Graph_traits::face_descriptor face_descriptor;

    typedef typename Surface_mesh::template Property_map<edge_descriptor, FT> Edge_property_map;

    typedef Unfolded_triangle_2<Kernel>     Unfolded_triangle;
    typedef Face_values<Kernel>             Face_values;

public:
    typedef typename Surface_mesh::template Property_map<face_descriptor, Face_values> Face_values_map;

    typedef typename Traits::SingleQueue    SingleQueue;

    typedef Polygon_mesh_processing::Barycentric_coordinates<FT> Barycentric_coordinates;
    typedef Polygon_mesh_processing::Face_location<Surface_mesh, FT> Face_location;

private:
    const Traits m_traits;
    Surface_mesh& m_mesh;

    Edge_property_map m_edge_lengths;
    Face_values_map m_face_values;
    std::vector<Face_location> m_target_face_locations;

    typedef std::queue<halfedge_descriptor> Halfedge_queue;
    Halfedge_queue m_A, m_B;

public:
    Surface_mesh_approximate_shortest_path(Surface_mesh& mesh,
                                const Traits& traits = Traits())
        : m_mesh(mesh), m_A(), m_B()
        {
            bool created_edge_property_map, created_face_property_map;
            boost::tie(m_edge_lengths, created_edge_property_map) = m_mesh.template add_property_map<edge_descriptor, FT>("edge_lengths");
            assert(created_edge_property_map);

            boost::tie(m_face_values, created_face_property_map) = m_mesh.template add_property_map<face_descriptor, Face_values>("face_values");
            assert(created_face_property_map);
        };

    Edge_property_map& Get_edge_length_map() { return m_edge_lengths; };

    class Compute_squared_edge_length {
    public:
        Compute_squared_edge_length() {};

        FT operator() (Surface_mesh& tmesh, halfedge_descriptor h)
        {
            vertex_descriptor v1 = tmesh.source(h);
            vertex_descriptor v2 = tmesh.target(h);

            return operator() (tmesh, v1, v2);
        }

        FT operator() (Surface_mesh& tmesh, vertex_descriptor v1, vertex_descriptor v2)
        {
            Point_3 v1_point = tmesh.point(v1);
            Point_3 v2_point = tmesh.point(v2);

            double length = squared_distance(v1_point, v2_point);
            return length;
        }
    };

    class Find_edge_length_and_update_property_map
    {
    public:
        typedef FT result_type;

    public:
        Find_edge_length_and_update_property_map() {};

        result_type operator() (Surface_mesh& mesh, halfedge_descriptor h)
        {
            return Compute_squared_edge_length()(mesh, h);
        }

        result_type operator() (Surface_mesh& mesh, halfedge_descriptor h, Edge_property_map& edge_lengths)
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
    typedef Find_edge_length_and_update_property_map Find_edge_length;

    class Unfold_triangle_3_along_halfedge
    {
    public:
        typedef Unfolded_triangle result_type;

        public:
            Unfold_triangle_3_along_halfedge() {}

            result_type operator() (Surface_mesh& tmesh, halfedge_descriptor h)
            {
                Unfolded_triangle unfolded_triangle;

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

                return unfolded_triangle;
            }

            result_type operator() (Surface_mesh& tmesh, halfedge_descriptor h, Edge_property_map& edge_lengths)
            {
                Unfolded_triangle unfolded_triangle;

                // edge length
                FT e0 = Find_edge_length_and_update_property_map()(tmesh, h, edge_lengths);

                // second point
                unfolded_triangle.B = Point_2(sqrt(e0), 0.);

                // second edge length
                halfedge_descriptor next_halfedge = tmesh.next(h);
                FT e1 = Find_edge_length_and_update_property_map()(tmesh, next_halfedge, edge_lengths);

                // third edge length
                halfedge_descriptor prev_halfedge = tmesh.next(next_halfedge);
                FT e2 = Find_edge_length_and_update_property_map()(tmesh, prev_halfedge, edge_lengths);

                // third point
                FT Px = (e0 + (e2 - e1)) / (2*sqrt(e0));
                FT Py = sqrt(e2 - square(Px));
                unfolded_triangle.P = Point_2(Px, Py);

                return unfolded_triangle;
            }
        };
    typedef Unfold_triangle_3_along_halfedge unfold_triangle_3;

    class Reconstruct_source_point_in_triangle_tangent_space{
    public:
        typedef Point_2 result_type;

    public:
        Reconstruct_source_point_in_triangle_tangent_space() {}

        result_type operator() (Surface_mesh& tmesh,
                                halfedge_descriptor h,
                                Edge_property_map& edge_lengths,
                                Face_values_map& face_values)
        {
            FT e0 = Find_edge_length()(tmesh, h, edge_lengths);
            if (is_border(h, tmesh)) {
                std::cerr << "halfedge opposite to " << h << " is on border and hence there is no way to reconstruct the source" << std::endl;
            }

            // find the correct entries in the face_values_map
            face_descriptor face = tmesh.face(h);
            vertex_descriptor A = tmesh.target(h); // this is swapped because target(h) == source(opposite(h))
            vertex_descriptor B = tmesh.source(h);

            int A_loc = vertex_index_in_face(A, face, tmesh);
            int B_loc = vertex_index_in_face(B, face, tmesh);

            FT d2A = face_values[face].d2verts[A_loc];
            FT d2B = face_values[face].d2verts[B_loc];

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
            return S;
        }
    };
    typedef Reconstruct_source_point_in_triangle_tangent_space Reconstruct_source_point;

    class Construct_triangle_centroid_2
    {
    public:
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

    class Get_heuristic_parameter
    {
    public:
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

    class Construct_heuristic_point_2
    {
    public:
        typedef Point_2 result_type;

    public:
        Construct_heuristic_point_2() {}

        result_type operator() (Surface_mesh& mesh, halfedge_descriptor h, Edge_property_map& edge_lengths)
        {
            // get barycenter
            Unfolded_triangle unfolded_triangle = unfold_triangle_3()(mesh, h, edge_lengths);
            Point_2 C = Construct_triangle_centroid_2()(unfolded_triangle.B, unfolded_triangle.P);

            return operator() (mesh, h, unfolded_triangle.P, C, edge_lengths);
        }

        result_type operator() (Surface_mesh& mesh, halfedge_descriptor h, Point_2 P, Point_2 C, Edge_property_map& edge_lengths)
        {
            // get edge lengths
            FT e0 = edge_lengths[mesh.edge(h)];
            FT e1 = edge_lengths[mesh.edge(mesh.next(h))];
            FT e2 = edge_lengths[mesh.edge(mesh.prev(h))];

            // look up the blending weight lambda
            FT lambda = Get_heuristic_parameter()(e0, e1, e2, P);

            // compute heuristic point coordinates
            FT Qx = lambda * C.x() + (1-lambda) * P.x();
            FT Qy = lambda * C.y() + (1-lambda) * P.y();

            return Point_2(Qx, Qy);
        }
    };

    class Edge_intersection_test
    {
    public:
        typedef typename Kernel::Left_turn_2    Left_turn_2;

        // TO BE CHANGED!!!
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


    void add_source_point(Point_3 source_point)
    {
        // locate source_point on mesh
        // NOTE: if the point is not on the mesh, it will be projected onto the mesh
        Face_location source_point_location = Polygon_mesh_processing::locate(source_point, m_mesh);

        add_source_point(source_point_location);
    };

    void add_source_point(Face_location source_point)
    {
        // get coordinates of source_point
        Point_3 source_point_coordinates = Polygon_mesh_processing::construct_point(source_point, m_mesh);
        face_descriptor face = source_point.first;
        halfedge_descriptor h0 = m_mesh.halfedge(face);

        // set sigma and d for the source triangle
        m_face_values[face].sigma = 0.;
        Face_location barycenter;
        barycenter.first = face;
        barycenter.second = {1./3, 1./3, 1./3};
        Point_3 face_barycenter = Polygon_mesh_processing::construct_point(barycenter, m_mesh);
        m_face_values[face].d = sqrt(squared_distance(source_point_coordinates, face_barycenter));

        // set distances to face vertices
        for (halfedge_descriptor h : m_mesh.halfedges_around_face(h0))
        {
            Find_edge_length()(m_mesh, h, m_edge_lengths);

            vertex_descriptor v_idx = m_mesh.source(h);
            //std::cout << "vertex index " << v_idx << std::endl;
            int local_v_idx = vertex_index_in_face(v_idx, face, m_mesh);
            //std::cout << "local vertex index " << local_v_idx << std::endl;
            FT dist = squared_distance(source_point_coordinates, m_mesh.point(v_idx));
            m_face_values[face].d2verts[local_v_idx] = dist;
            //std::cout << std::endl;

            m_A.push(h);
        }

        //std::cout << face << " with values " << m_face_values[face] << std::endl;
    };

    void add_target(Point_3 target_point)
    {
        // locate source_point on mesh
        // NOTE: if the point is not on the mesh, it will be projected onto the mesh
        Face_location target_face_location = Polygon_mesh_processing::locate(target_point, m_mesh);

        add_target(target_face_location);
    };

    void add_target(Face_location target_point)
    {
        m_target_face_locations.push_back(target_point);
    }

    std::pair<bool, int> belongs_to_target_face(halfedge_descriptor h)
    {
        std::pair<bool, int> is_in_target_face = { false, -1 };
        face_descriptor face = m_mesh.face(h);
        for (int i = 0; i < m_target_face_locations.size(); i++)
        {
            if (m_target_face_locations[i].first == face)
            {
                is_in_target_face.first = true;
                is_in_target_face.second = i;
            }
        }

        return is_in_target_face;
    }

    //=== PROPAGATION ===//

    typedef std::pair<bool, bool> intersection_result;

    std::pair<FT, FT> get_new_dist(intersection_result intersection,
                                   halfedge_descriptor h,
                                   Point_2 C,
                                   Point_2 S)
    {
        face_descriptor face = m_mesh.face(h);
        FT e0 = sqrt(m_edge_lengths[m_mesh.edge(h)]);

        std::pair<FT, FT> prev_geodesic_dist;
        prev_geodesic_dist.first = m_face_values[face].sigma;
        prev_geodesic_dist.second = m_face_values[face].d;

        std::pair<FT, FT> new_geodesic_dist;
        if (intersection.first)
        {
            // we have an intersection
            new_geodesic_dist.first = prev_geodesic_dist.first;
            new_geodesic_dist.second = new_geodesic_dist.first + sqrt( squared_distance(C, S) );
        }
        else
        {
            if (intersection.second)
            {
                // right turn (source is to the right of B)
                new_geodesic_dist.first = sqrt( square(S.x()-e0) + square(S.y()) ) + prev_geodesic_dist.first;
                new_geodesic_dist.second = new_geodesic_dist.first + sqrt( square(C.x()-e0) + square(C.y()) );
            }
            else
            {
                // left turn (source is to the left of A, t<0 case)
                new_geodesic_dist.first = sqrt(square(S.x()) + square(S.y())) + prev_geodesic_dist.first;
                new_geodesic_dist.second = new_geodesic_dist.first + sqrt(square(C.x()) + square(C.y()));
            }
        }

        return new_geodesic_dist;
    }

    std::array<int, 3> get_local_vertex_indices(halfedge_descriptor h)
    {
        face_descriptor face = m_mesh.face(h);
        vertex_descriptor A = m_mesh.source(h);
        vertex_descriptor B = m_mesh.target(h);
        vertex_descriptor P = m_mesh.target(m_mesh.next(h));

        int local_idx_A_in_face = vertex_index_in_face(A, face, m_mesh);
        int local_idx_B_in_face = vertex_index_in_face(B, face, m_mesh);
        int local_idx_P_in_face = vertex_index_in_face(P, face, m_mesh);

        std::array<int, 3> local_vertex_indices = {local_idx_A_in_face, local_idx_B_in_face, local_idx_P_in_face};

        return local_vertex_indices;
    }

    void set_squared_vertex_distances_intersection(halfedge_descriptor oppo_h,
                                                   Point_2 P_coords,
                                                   Point_2 S_coords)
    {
        face_descriptor face = m_mesh.face(oppo_h);
        halfedge_descriptor h = m_mesh.opposite(oppo_h);
        face_descriptor prev_face = m_mesh.face(h);

        Face_values prev_face_values = m_face_values[prev_face];
        std::array<FT,3> prev_face_vertex_distances = prev_face_values.d2verts;

        vertex_descriptor A = m_mesh.target(oppo_h);
        int local_idx_A_in_face = vertex_index_in_face(A, face, m_mesh);
        int local_idx_A_in_prev_face = vertex_index_in_face(A, prev_face, m_mesh);
        m_face_values[face].d2verts[local_idx_A_in_face] = prev_face_vertex_distances[local_idx_A_in_prev_face];

        vertex_descriptor B = m_mesh.source(oppo_h);
        int local_idx_B_in_face = vertex_index_in_face(B, face, m_mesh);
        int local_idx_B_in_prev_face = vertex_index_in_face(B, prev_face, m_mesh);
        m_face_values[face].d2verts[local_idx_B_in_face] = prev_face_vertex_distances[local_idx_B_in_prev_face];

        vertex_descriptor P = m_mesh.target(m_mesh.next(oppo_h));
        int local_idx_P_in_face = vertex_index_in_face(P, face, m_mesh);
        m_face_values[face].d2verts[local_idx_P_in_face] = squared_distance(P_coords, S_coords);
    }

    void set_squared_vertex_distances_left_turn(halfedge_descriptor h)
    {
        face_descriptor face = m_mesh.face(h);
        std::array<int, 3> local_vertex_indices = get_local_vertex_indices(h);

        m_face_values[face].d2verts[local_vertex_indices[0]] = 0.;
        m_face_values[face].d2verts[local_vertex_indices[1]] = m_edge_lengths[m_mesh.edge(h)];
        m_face_values[face].d2verts[local_vertex_indices[2]] = m_edge_lengths[m_mesh.edge(m_mesh.prev(h))];
    }

    void set_squared_vertex_distances_right_turn(halfedge_descriptor h)
    {
        face_descriptor face = m_mesh.face(h);
        std::array<int, 3> local_vertex_indices = get_local_vertex_indices(h);

        m_face_values[face].d2verts[local_vertex_indices[0]] = m_edge_lengths[m_mesh.edge(h)];
        m_face_values[face].d2verts[local_vertex_indices[1]] = 0.;
        m_face_values[face].d2verts[local_vertex_indices[2]] = m_edge_lengths[m_mesh.edge(m_mesh.next(h))];
    }

    void set_squared_vertex_distances(intersection_result intersection,
                                      halfedge_descriptor h,
                                      std::pair<FT, FT> new_geodesic_dist,
                                      Point_2 P,
                                      Point_2 S)
    {
        halfedge_descriptor oppo_h = m_mesh.opposite(h);
        face_descriptor next_face = m_mesh.face(oppo_h);
        m_face_values[next_face].sigma = new_geodesic_dist.first;
        m_face_values[next_face].d = new_geodesic_dist.second;

        if (intersection.first)
        {
            // we have an intersection
            set_squared_vertex_distances_intersection(oppo_h, P, S);
        }
        else
        {
            if (intersection.second)
            {
                // right turn (source is to the right of B)
                set_squared_vertex_distances_right_turn(oppo_h);
            }
            else
            {
                // left turn (source is to the left of A, t<0 case)
                set_squared_vertex_distances_left_turn(oppo_h);
            }
        }
    }

    std::pair<bool, FT> update_face_values(intersection_result intersection,
                            halfedge_descriptor h,
                            Point_2 C,
                            Point_2 P,
                            Point_2 S)
    {
        std::pair<FT, FT> new_geodesic_dist = get_new_dist(intersection, h, C, S);
        face_descriptor face = m_mesh.face(m_mesh.opposite(h));

        if (face.idx() == 19)
        {
            std::cout << "BEFORE UPDATE: " << m_face_values[face] << std::endl;
            std::cout << "new geodesic distances: " << new_geodesic_dist.first << ", " << new_geodesic_dist.second << std::endl;
        }

        if (new_geodesic_dist.second < m_face_values[face].d)
        {
            set_squared_vertex_distances(intersection, h, new_geodesic_dist, P, S);

            return { true, new_geodesic_dist.second };
        }

        if (face.idx() == 19)
        {
            std::cout << "AFTER UPDATE: " << m_face_values[face] << std::endl;
        }

        return { false, FT(-1.) };
    }

    template<class SkipCondition>
    void enqueue_new_halfedges(halfedge_descriptor h, FT new_geodesic_dist, int iter)
    {
        FT dist_i(iter);
        bool enqueue_in_A = SkipCondition()();

        if (enqueue_in_A)
        {
            m_A.push(m_mesh.next(h));
            m_A.push(m_mesh.prev(h));
        }
        else
        {
            m_B.push(m_mesh.next(h));
            m_B.push(m_mesh.prev(h));
        }
    }

    void propagate_over_halfedge(halfedge_descriptor halfedge, int iter)
    {
        halfedge_descriptor opposite_halfedge = m_mesh.opposite(halfedge);
        if (m_mesh.is_border(opposite_halfedge)) { return; }
        FT e0 = sqrt(Find_edge_length()(m_mesh, halfedge, m_edge_lengths));

        // constructions
        auto unfolded_triangle = unfold_triangle_3()(m_mesh, opposite_halfedge, m_edge_lengths);
        Point_2 A(FT(0.), FT(0.));
        Point_2 B = unfolded_triangle.B;
        Point_2 P = unfolded_triangle.P;

        Point_2 C, Q;
        std::pair<bool, int> is_target_face = belongs_to_target_face(opposite_halfedge);
        if (is_target_face.first)
        {
            Barycentric_coordinates coords = m_target_face_locations[is_target_face.second].second;
            std::array<int, 3> local_vert_idx = get_local_vertex_indices(opposite_halfedge);
            // C becomes the target point and we also test visibility with respect to the target point (C = Q)
            FT C1 = coords[local_vert_idx[1]] * B.x() + coords[local_vert_idx[2]] * P.x();
            FT C2 = coords[local_vert_idx[2]] * P.y();
            C = { C1, C2 };
            Q = C;
        }
        else
        {
            C = Construct_triangle_centroid_2()(B, P);
            Q = Construct_heuristic_point_2()(m_mesh, opposite_halfedge, P, C, m_edge_lengths);
        }
        Point_2 S = Reconstruct_source_point()(m_mesh, halfedge, m_edge_lengths, m_face_values);

        // intersection test
        auto intersection = Edge_intersection_test()(S, Q, A, B);
        std::pair<bool, FT> update_result = update_face_values(intersection, halfedge, C, P, S);
        if (update_result.first) { enqueue_new_halfedges<SingleQueue>(opposite_halfedge, update_result.second, iter); }
    }

    void propagate_geodesic_source(Face_location source)
    {
        // initialise the algorithm
        add_source_point(source);

        int iter = 1;
        while (!m_A.empty())
        {
            std::cout << "iteration " << iter << std::endl;
            std::cout << "updating faces: ";
            // iterate over halfedges in queue and update face values
            while (!m_A.empty()) {
                halfedge_descriptor h = m_A.front();
                m_A.pop();

                propagate_over_halfedge(h, iter);
                std::cout << m_mesh.face(h) << "\t" << std::endl;
            }

            m_A = m_B;
            std::queue<halfedge_descriptor> empty;
            std::swap( m_B, empty );
            iter++;
        }

        std::cout << "algorithm has concluded and the following face values were obtained:" << std::endl;
        for (face_descriptor f : faces(m_mesh))
        {
            std::cout << f.idx() << ":\t" <<  m_face_values[f] << std::endl;
        }
    }

    void propagate_geodesic_source(Point_3 source)
    {
        Face_location source_location = Polygon_mesh_processing::locate(source, m_mesh);

        propagate_geodesic_source(source_location);
    }

    // output queries
    Face_values& get_face_values(face_descriptor face)
    {
        return m_face_values[face];
    }
};

}

#endif // CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_H
