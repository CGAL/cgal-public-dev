// Copyright (c) 2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is a part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY, AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Keyu Chen， Dmitry Anisimov.

/*!
  \file Harmonic_mesh.h
*/

#ifndef CGAL_HARMONIC_MESH_H
#define CGAL_HARMONIC_MESH_H

#include <CGAL/license/Barycentric_coordinates_2.h>

#include <CGAL/disable_warnings.h>

// CGAL headers.
#include <CGAL/assertions.h>
#include <CGAL/Polygon_2_algorithms.h>

// Barycentric coordinates headers.
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>

// Boost headers.
#include <boost/optional/optional.hpp>

// Number utils headers
#include <CGAL/number_utils.h>

// Delaunay_triangulation headers
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>



// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Introduction of Harmonic_mesh_2

template<class Traits>
    class Harmonic_mesh_2
{

public:

    /// \name Types
    /// @{

    /// Number type.
    typedef typename Traits::FT         FT;
    typedef typename std::vector<FT>    FT_vector;

    /// Point type.
    typedef typename Traits::Point_2 Point_2;
    typedef typename std::vector<Point_2> Point_vector;



    /// @}

    // \name Creation
    Harmonic_mesh_2(const std::vector<typename Traits::Point_2> &vertices) :
        vertex(vertices),
        number_of_vertices(vertex.size())
    {
        insert_constraint(cdt, vertex);
        detect_shape_scale(shape_scale, vertex);
        // Initialize some private parameters here.
    }

    // Create partition using assigned number of mesh vertices.
    void create_mesh(FT &partition_constraint)
    {
        FT max_edge_length = shape_scale * partition_constraint;
        create_denaulay_mesh(max_edge_length);
    }

    Point_vector locate_point(Point_2 &query)
    {
        Point_vector triangle_vertex;
        locate_triangle_face(cdt, query, triangle_vertex);
        return triangle_vertex;
    }

    Point_vector all_vertices()
    {
        Point_vector all_mesh_vertices;
        list_all_vertices(cdt, all_mesh_vertices);
        return all_mesh_vertices;
    }

    Point_vector get_neighbor(Point_2 &vertex)
    {
        Point_vector neighbors;
        list_all_neighbors(cdt, neighbors, vertex);
        return neighbors;
    }

    void print_information()
    {
        //std::cout<<"Prior class function available."<<std::endl;
    }

private:

    /// Constrained_Delaunay_triangulation type
    typedef CGAL::Triangulation_vertex_base_2<Traits> Vb;
    typedef CGAL::Delaunay_mesh_face_base_2<Traits> Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
    typedef CGAL::Constrained_Delaunay_triangulation_2<Traits, Tds> CDT;
    typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
    //typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Delaunay_mesher;

    typedef typename CDT::Vertex_iterator Vertex_iterator;
    typedef typename CDT::Vertex_circulator Vertex_circulator;
    typedef typename CDT::Vertex_handle Vertex_handle;
    typedef typename CDT::Vertex Vertex;
    typedef typename CDT::Face_handle Face_handle;
    typedef typename CDT::Face Face;
    typedef typename CDT::Point CDT_Point;

    typedef typename Traits::Vector_2 Vector_2;
    //typedef typename std::vector<Point_2> Point_vector;

    const Point_vector &vertex;

    const size_t number_of_vertices;

    FT shape_scale;

    CDT cdt;
    //Delaunay_mesher delaunay_mesher;

    void insert_constraint(CDT &cdt, const Point_vector &vertex)
    {

        for (size_t i = 0; i < number_of_vertices; ++i) {
            size_t ip = (i + 1) % number_of_vertices;

            Vertex_handle va = cdt.insert(CDT_Point(vertex[i]));
            Vertex_handle vb = cdt.insert(CDT_Point(vertex[ip]));

            cdt.insert_constraint(va, vb);
        }

        //delaunay_mesher(cdt);
    }

    void detect_shape_scale(FT &shape_scale, const Point_vector &vertex)
    {
        FT max_x, min_x, max_y, min_y;
        if(vertex.size()) {
            max_x = vertex[0].x();
            min_x = vertex[0].x();
            max_y = vertex[0].y();
            min_y = vertex[0].y();
        }
        for (size_t i = 1; i < number_of_vertices; ++i) {
            max_x = std::max<FT>(max_x, vertex[i].x());
            max_y = std::max<FT>(max_y, vertex[i].y());
            min_x = std::min<FT>(min_x, vertex[i].x());
            min_y = std::min<FT>(min_y, vertex[i].y());
        }
        shape_scale = std::min<FT>(max_x - min_x, max_y - min_y);
    }

    void create_denaulay_mesh(FT &max_edge_length)
    {
        //delaunay_mesher.set_criteria(Criteria(0.125, max_edge_length));
        std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
        //CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, 0.5));
        CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, max_edge_length));
        std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
    }

    void locate_triangle_face(CDT &cdt, Point_2 query, Point_vector &triangle_vertex)
    {
        Face_handle triangle_face_handle = cdt.locate(query);
        Face triangle_face = *triangle_face_handle;

        Vertex_handle first_vertex_handle = triangle_face.vertex(0);
        Vertex_handle second_vertex_handle = triangle_face.vertex(1);
        Vertex_handle third_vertex_handle = triangle_face.vertex(2);

        Vertex first_vertex = *first_vertex_handle;
        Vertex second_vertex = *second_vertex_handle;
        Vertex third_vertex = *third_vertex_handle;

        Point_2 first_vertex_location = first_vertex.point();
        Point_2 second_vertex_location = second_vertex.point();
        Point_2 third_vertex_location = third_vertex.point();

        triangle_vertex.push_back(first_vertex_location);
        triangle_vertex.push_back(second_vertex_location);
        triangle_vertex.push_back(third_vertex_location);
    }

    void list_all_vertices(CDT &cdt, Point_vector &all_mesh_vertices)
    {
        Vertex_iterator all_vertices_begin = cdt.finite_vertices_begin();

        while(all_vertices_begin != cdt.finite_vertices_end())
        {
            Vertex_handle current_vertex_handle = all_vertices_begin;
            Vertex current_vertex = *current_vertex_handle;
            all_mesh_vertices.push_back(current_vertex.point());

            all_vertices_begin++;
        }
    }

    void list_all_neighbors(CDT &cdt, Point_vector &neighbors, Point_2 &vertex)
    {
        Vertex query = Vertex(vertex);
        Vertex_circulator all_neighbors_begin = cdt.incident_vertices(&query);
        Vertex_circulator all_neighbors_end = all_neighbors_begin;
        Vertex_handle first_neighbor_handle = all_neighbors_begin;
        Vertex first_neighbor = *first_neighbor_handle;
        neighbors.push_back(first_neighbor.point());
        all_neighbors_begin++;
        while(all_neighbors_begin != all_neighbors_end) {
            Vertex_handle current_neighbor_handle = all_neighbors_begin;
            Vertex current_neighbor = *current_neighbor_handle;
            neighbors.push_back(current_neighbor.point());
            all_neighbors_begin++;
        }
    }
};

}  // namespace Barycentric_coordinates

}  // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_HARMONIC_MESH_H
