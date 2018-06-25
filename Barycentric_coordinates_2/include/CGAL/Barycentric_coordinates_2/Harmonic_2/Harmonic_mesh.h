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
    Harmonic_mesh_2(const std::vector<typename Traits::Point_2> &vertices, const Traits &b_traits) :
        vertex(vertices),
        barycentric_traits(b_traits),
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

    typedef typename CDT::Vertex_handle Vertex_handle;
    typedef typename CDT::Point CDT_Point;

    typedef typename Traits::Vector_2 Vector_2;
    //typedef typename std::vector<Point_2> Point_vector;

    const Point_vector &vertex;

    const Traits &barycentric_traits;

    const size_t number_of_vertices;

    FT shape_scale;

    CDT cdt;
    //Delaunay_mesher delaunay_mesher;

    void insert_constraint(CDT &cdt, const Point_vector &vertex) {

        for (size_t i = 0; i < number_of_vertices; ++i) {
            size_t ip = (i + 1) % number_of_vertices;

            Vertex_handle va = cdt.insert(CDT_Point(vertex[i]));
            Vertex_handle vb = cdt.insert(CDT_Point(vertex[ip]));

            cdt.insert_constraint(va, vb);
        }

        //delaunay_mesher(cdt);
    }

    void detect_shape_scale(FT &shape_scale, const Point_vector &vertex) {
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

    void create_denaulay_mesh(FT &max_edge_length) {
        //delaunay_mesher.set_criteria(Criteria(0.125, max_edge_length));
        std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
        //CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, 0.5));
        CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, max_edge_length));
        std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
    }
};

}  // namespace Barycentric_coordinates

}  // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_HARMONIC_MESH_H
