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

#include <CGAL/Polygon_mesh_processing/locate.h>>

namespace CGAL {

template <class Kernel>
struct Face_values {
    typedef typename Kernel::FT FT;

    FT sigma;
    FT d;
    std::array<FT,3> d2verts;

    Face_values(FT _sigma=1., FT _d=2., std::array<FT,3> _d2verts = {3., 4., 5.}) // we need some default values along the lines of CGAL::infty
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

    typedef typename Traits::Triangle_mesh Triangle_mesh;
    typedef boost::graph_traits<Triangle_mesh> Graph_traits;

    //typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
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

private:
    const Traits m_traits;
    Triangle_mesh& mesh;

    Edge_property_map m_edge_lengths;
    Face_values_map m_face_values;

public:
    Surface_mesh_approximate_shortest_path(Triangle_mesh& mesh,
                                const Traits& traits = Traits())
        : mesh(mesh)
        {
            //std::cout << mesh.number_of_faces() << std::endl;

            bool created_edge_property_map, created_face_property_map;
            boost::tie(m_edge_lengths, created_edge_property_map) = mesh.template add_property_map<edge_descriptor, FT>("edge_lengths");
            assert(created_edge_property_map);

            boost::tie(m_face_values, created_face_property_map) = mesh.template add_property_map<face_descriptor, Face_values>("face_values");
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

        typedef CGAL::Polygon_mesh_processing::Barycentric_coordinates<FT> Barycentric_coordinates;

    class Add_source_point
    {
    public:
        typedef std::pair<face_descriptor, Barycentric_coordinates> face_location;

    public:
        Add_source_point() {};

        void operator() (face_location source_point)
        {
            m_face_values[source_point.first]
        }
    };

};

}

#endif // CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_H
