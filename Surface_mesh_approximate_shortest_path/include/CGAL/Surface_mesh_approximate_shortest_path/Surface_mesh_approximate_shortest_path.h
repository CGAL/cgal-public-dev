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

template <class Kernel>
struct Face_values {
    typedef typename Kernel::FT FT;

    FT sigma;
    FT d;
    std::array<FT,3> d2verts;

    Face_values(FT _sigma=0., FT _d=std::numeric_limits<FT>::max(), std::array<FT,3> _d2verts = {-1., -1., -1.}) // we need some default values along the lines of CGAL::infty
        : sigma(_sigma), d(_d), d2verts(_d2verts) {}       // so that the comparator with any real number says that it is larger

    friend std::ostream & operator <<(std::ostream& stream, const Face_values vals)
    {
        return ( stream << vals.sigma << "\t" << vals.d << "\t"
                       << vals.d2verts[0]  << "\t" << vals.d2verts[1] << "\t" << vals.d2verts[2] << std::endl);
    };
};

template<class Traits>
class Surface_mesh_approximate_shortest_path
{
public:
    typedef typename Traits::Kernel Kernel;
    typedef typename Kernel::FT FT;

    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Point_3 Point_3;

    typedef typename Traits::Triangle_mesh Triangle_mesh;
    typedef boost::graph_traits<Triangle_mesh> Graph_traits;

    typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
    typedef typename Graph_traits::edge_descriptor edge_descriptor;
    typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
    typedef typename Graph_traits::face_descriptor face_descriptor;

    typedef typename Triangle_mesh::template Property_map<edge_descriptor, FT> Edge_property_map;

    typedef Face_values<Kernel> Face_values;

public:
    typedef typename Triangle_mesh::template Property_map<face_descriptor, Face_values> Face_values_map;

    typedef typename Traits::Compute_squared_edge_length                        Compute_squared_edge_length;
    typedef typename Traits::Unfold_triangle_3_along_halfedge                   Unfold_triangle_3_along_halfedge;
    typedef typename Traits::Reconstruct_source_point_in_triangle_tangent_space Reconstruct_source_point_in_triangle_tangent_space;
    typedef typename Traits::Construct_triangle_centroid_2                      Construct_triangle_centroid_2;
    typedef typename Traits::Construct_heuristic_point_2                        Construct_heuristic_point_2;
    typedef typename Traits::Edge_intersection_test                             Edge_intersection_test;

    typedef Polygon_mesh_processing::Barycentric_coordinates<FT> Barycentric_coordinates;
    typedef Polygon_mesh_processing::Face_location<Triangle_mesh, FT> Face_location;

private:
    const Traits m_traits;
    Triangle_mesh& m_mesh;

    Edge_property_map m_edge_lengths;
    Face_values_map m_face_values;
    std::vector<Face_location> m_target_face_locations;

    typedef std::queue<halfedge_descriptor> Halfedge_queue;
    Halfedge_queue m_A, m_B;

public:
    Surface_mesh_approximate_shortest_path(Triangle_mesh& mesh,
                                const Traits& traits = Traits())
        : m_mesh(mesh), m_A(), m_B()
        {
            //std::cout << mesh.number_of_faces() << std::endl;

            bool created_edge_property_map, created_face_property_map;
            boost::tie(m_edge_lengths, created_edge_property_map) = m_mesh.template add_property_map<edge_descriptor, FT>("edge_lengths");
            assert(created_edge_property_map);

            boost::tie(m_face_values, created_face_property_map) = m_mesh.template add_property_map<face_descriptor, Face_values>("face_values");
            assert(created_face_property_map);

            // test initialization of face_value_map
            //std::cout << "face values for face 0:" << std::endl << m_face_values[face_descriptor(0)] << std::endl;
        };

    Edge_property_map& Get_edge_length_map() { return m_edge_lengths; };

    Unfold_triangle_3_along_halfedge unfold_triangle_3_along_halfedge_object()
        { return m_traits.unfold_triangle_3_along_halfedge_object(); }
    Reconstruct_source_point_in_triangle_tangent_space reconstruct_source_point_in_triangle_tangent_space_object()
        { return m_traits.reconstruct_source_point_in_triangle_tangent_space_object(); }
    Construct_triangle_centroid_2 construct_centroid_object()
        { return m_traits.construct_centroid_2_object(); };
    Construct_heuristic_point_2 construct_heuristic_point_object()
        { return m_traits.construct_heuristic_point_2_object(); };
    Edge_intersection_test edge_intersection_test_object()
        { return m_traits.edge_intersection_test_object(); };

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
            m_traits.find_edge_length_and_update_property_map_object()(m_mesh, h, m_edge_lengths);

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

    bool update_face_values(intersection_result intersection,
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

            return true;
        }

        if (face.idx() == 19)
        {
            std::cout << "AFTER UPDATE: " << m_face_values[face] << std::endl;
        }

        return false;
    }

    void enqueue_new_halfedges(halfedge_descriptor h)
    {
        m_A.push(m_mesh.next(h));
        m_A.push(m_mesh.prev(h));
    }

    void propagate_over_halfedge(halfedge_descriptor halfedge)
    {
        halfedge_descriptor opposite_halfedge = m_mesh.opposite(halfedge);
        if (m_mesh.is_border(opposite_halfedge)) { return; }
        FT e0 = sqrt(m_traits.find_edge_length_and_update_property_map_object()(m_mesh, halfedge, m_edge_lengths));

        // constructions
        auto unfolded_triangle = unfold_triangle_3_along_halfedge_object()(m_mesh, opposite_halfedge, m_edge_lengths);
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
            C = construct_centroid_object()(B, P);
            Q = construct_heuristic_point_object()(m_mesh, opposite_halfedge, P, C, m_edge_lengths);
        }
        Point_2 S = reconstruct_source_point_in_triangle_tangent_space_object()(m_mesh, halfedge, m_edge_lengths, m_face_values);

        // intersection test
        auto intersection = edge_intersection_test_object()(S, Q, A, B);
        auto update_result = update_face_values(intersection, halfedge, C, P, S);
        if (update_result) { enqueue_new_halfedges(opposite_halfedge); }
    }

    void propagate_geodesic_source(Face_location source)
    {
        // initialise the algorithm
        add_source_point(source);

        // iterate over halfedges in queue and update face values
        while (!m_A.empty()) {
            halfedge_descriptor h = m_A.front();
            m_A.pop();

            propagate_over_halfedge(h);
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
