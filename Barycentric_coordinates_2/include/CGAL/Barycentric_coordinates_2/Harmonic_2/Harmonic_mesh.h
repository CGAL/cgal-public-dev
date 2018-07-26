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
#include <CGAL/property_map.h>

// Barycentric coordinates headers.
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>

// Boost headers.
#include <boost/optional/optional.hpp>

// Number utils headers
#include <CGAL/number_utils.h>

// Delaunay_triangulation headers
#include <CGAL/Triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <CGAL/IO/File_poly.h>



// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Introduction of Harmonic_mesh_2

template<class Traits>
    class Harmonic_mesh
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
    Harmonic_mesh(const std::vector<typename Traits::Point_2> &vertices, const Traits &barycentric_traits) :
        m_vertex(vertices),
        m_barycentric_traits(barycentric_traits)
    {
        insert_constraint(cdt);
        detect_shape_scale(shape_scale);
        // Initialize some private parameters here.
    }

    // Create partition using assigned number of mesh vertices.
    void create_mesh(FT &partition_constraint)
    {
        FT max_edge_length = shape_scale * partition_constraint;
        create_denaulay_mesh(max_edge_length);
    }

    Point_vector get_all_vertices()
    {
        Point_vector all_mesh_vertices;
        list_all_vertices(cdt, all_mesh_vertices);
        return all_mesh_vertices;
    }

    std::vector<int> get_neighbor(int i)
    {
        std::vector<int> neighbors;
        list_all_neighbors(cdt, neighbors, i);
        return neighbors;
    }

    std::vector<int> get_boundary_vertices()
    {
        std::vector<int> boundary;
        list_all_boundary_vertices(cdt, boundary);
        return boundary;
    }

    std::vector<int> locate_point(Point_2 query_point)
    {
        Point_2 query = query_point;
        std::vector<int> triangle_indices;
        triangle_indices.resize(3);
        locate_triangle_face(cdt, query, triangle_indices);

        return triangle_indices;
    }

    void save_mesh(Point_vector &mesh_vertices, std::ostream &f)
    {
        write_triangle_OBJ_file(cdt, mesh_vertices, f);
    }

private:

    /// Constrained_Delaunay_triangulation type
    struct VertexInfo2
    {
        VertexInfo2(){}

        int index;
        std::vector<int> neighbor;

        bool is_on_boundary;
    };

    typedef CGAL::Triangulation_vertex_base_with_info_2<VertexInfo2, Traits> Vbb;
    typedef CGAL::Triangulation_vertex_base_2<Traits, Vbb> Vb;
    typedef CGAL::Delaunay_mesh_face_base_2<Traits> Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
    typedef CGAL::Constrained_Delaunay_triangulation_2<Traits, Tds> CDT;
    typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;

    typedef typename CDT::Vertex_iterator Vertex_iterator;
    typedef typename CDT::Vertex_circulator Vertex_circulator;
    typedef typename CDT::Vertex_handle Vertex_handle;
    typedef typename CDT::Vertex Vertex;
    typedef typename CDT::Edge Edge;
    typedef typename CDT::Face_handle Face_handle;
    typedef typename CDT::Face Face;
    typedef typename CDT::Point CDT_Point;
    typedef typename CDT::Constrained_edges_iterator Constrained_edges_iterator;

    typedef typename Traits::Vector_2 Vector_2;

    const Point_vector &m_vertex;

    const Traits &m_barycentric_traits;

    FT shape_scale;

    CDT cdt;

    std::vector<Vertex_handle> Mesh_handles;

    void insert_constraint(CDT &cdt)
    {
        const size_t number_of_vertices = m_vertex.size();

        for (size_t i = 0; i < number_of_vertices; ++i) {
            size_t ip = (i + 1) % number_of_vertices;

            Vertex_handle va = cdt.insert(CDT_Point(m_vertex[i]));
            Vertex_handle vb = cdt.insert(CDT_Point(m_vertex[ip]));

            cdt.insert_constraint(va, vb);
        }
    }

    void detect_shape_scale(FT &shape_scale)
    {
        FT max_x, min_x, max_y, min_y;
        if(m_vertex.size()) {
            max_x = m_vertex[0].x();
            min_x = m_vertex[0].x();
            max_y = m_vertex[0].y();
            min_y = m_vertex[0].y();
        }

        const size_t number_of_vertices = m_vertex.size();

        for (size_t i = 1; i < number_of_vertices; ++i) {
            max_x = std::max<FT>(max_x, m_vertex[i].x());
            max_y = std::max<FT>(max_y, m_vertex[i].y());
            min_x = std::min<FT>(min_x, m_vertex[i].x());
            min_y = std::min<FT>(min_y, m_vertex[i].y());
        }
        shape_scale = std::min<FT>(max_x - min_x, max_y - min_y);
    }

    void create_denaulay_mesh(FT &max_edge_length)
    {
        const size_t number_of_vertices = m_vertex.size();
        //delaunay_mesher.set_criteria(Criteria(0.125, max_edge_length));
        std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
        //CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, 0.5));
        CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, max_edge_length));
        std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
    }

    void list_all_vertices(CDT &cdt, Point_vector &mesh_vertices)
    {
        int i = 0;
        for (Vertex_iterator vertex_handle = cdt.finite_vertices_begin(); vertex_handle != cdt.finite_vertices_end(); ++vertex_handle)
        {
            if(!cdt.is_infinite(vertex_handle)){
                Mesh_handles.push_back(vertex_handle);

                Vertex v = *vertex_handle;
                mesh_vertices.push_back(v.point());
                vertex_handle->info().index = i;
                i++;
            }
        }
    }

    void locate_triangle_face(const CDT &cdt, Point_2 query, std::vector<int> &triangle_indices)
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

        int indice1 = first_vertex_handle->info().index;
        int indice2 = second_vertex_handle->info().index;
        int indice3 = third_vertex_handle->info().index;

        triangle_indices[0] = indice1;
        triangle_indices[1] = indice2;
        triangle_indices[2] = indice3;
    }

    void list_all_neighbors(CDT &cdt, std::vector<int> &neighbors, int i)
    {
        Vertex_handle query_handle = Mesh_handles[i];

        Vertex_circulator all_neighbors_begin = cdt.incident_vertices(query_handle);
        Vertex_circulator all_neighbors_end = all_neighbors_begin;
        Vertex_handle first_neighbor_handle = all_neighbors_begin;
        if(!cdt.is_infinite(first_neighbor_handle)){
            neighbors.push_back(first_neighbor_handle->info().index);
        }

        all_neighbors_begin++;
        while(all_neighbors_begin != all_neighbors_end) {
            Vertex_handle current_neighbor_handle = all_neighbors_begin;
            if(!cdt.is_infinite(current_neighbor_handle)){
                neighbors.push_back(current_neighbor_handle->info().index);
            }
            all_neighbors_begin++;
        }
    }

    void list_all_boundary_vertices(CDT &cdt, std::vector<int> &boundary)
    {
        for(Constrained_edges_iterator start = cdt.constrained_edges_begin(); start != cdt.constrained_edges_end(); ++start){
            Edge constrained_edge = *start;
            Face_handle face_handle = constrained_edge.first;
            int i = constrained_edge.second;
            Face face = *face_handle;
            Vertex_handle boundary_vertex_first = face.vertex(CDT::cw(i));
            Vertex_handle boundary_vertex_second = face.vertex(CDT::ccw(i));
            boundary.push_back(boundary_vertex_first->info().index);
            boundary.push_back(boundary_vertex_second->info().index);
        }
    }

    void write_triangle_OBJ_file(CDT &cdt, Point_vector &mesh_vertices, std::ostream &f)
    {
        typedef typename CDT::Vertex_handle Vertex_handle;
        typedef typename CDT::Finite_vertices_iterator
          Finite_vertices_iterator;
        typedef typename CDT::Finite_faces_iterator
          Finite_faces_iterator;


        //for(Finite_vertices_iterator vit = cdt.finite_vertices_begin();
        //    vit != cdt.finite_vertices_end();
        //    ++vit)
        //{
        //    f << "v " << vit->point() << " 0.000000" << std::endl;
        //}
        for(size_t i = 0; i < mesh_vertices.size(); ++i)
        {
            f << "v " << mesh_vertices[i] << " 0.000000" << std::endl;
        }

        f << std::endl;

        for(Finite_faces_iterator fit = cdt.finite_faces_begin();
            fit != cdt.finite_faces_end();
            ++fit)
        {
            Face face = *fit;
            Vertex_handle first_vertex = face.vertex(0);
            Vertex_handle second_vertex = face.vertex(1);
            Vertex_handle third_vertex = face.vertex(2);
            f << "f "
              << first_vertex->info().index + 1 << " "
              << second_vertex->info().index + 1 << " "
              << third_vertex->info().index + 1
              << std::endl;
        }

        //f << std::endl;
    }
};

}  // namespace Barycentric_coordinates

}  // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_HARMONIC_MESH_H
