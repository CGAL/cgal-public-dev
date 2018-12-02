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
#include <CGAL/lloyd_optimize_mesh_2.h>
#include <CGAL/Triangulation_conformer_2.h>

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
        detect_minimal_edge_length(minimal_edge_length);
        // Initialize some private parameters here.
    }

    // Create partition using assigned number of mesh vertices.
    void create_mesh(FT &partition_constraint)
    {
        shape_bound = minimal_edge_length * partition_constraint;
        insert_constraint(cdt);
        create_denaulay_mesh();
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
    typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;

    typedef typename CDT::Vertex_iterator Vertex_iterator;
    typedef typename CDT::Vertex_circulator Vertex_circulator;
    typedef typename CDT::Vertex_handle Vertex_handle;
    typedef typename CDT::Vertex Vertex;
    typedef typename CDT::Edge Edge;
    typedef typename CDT::Edge_iterator Edge_iterator;
    typedef typename CDT::Face_handle Face_handle;
    typedef typename CDT::Face Face;
    typedef typename CDT::Point CDT_Point;
    typedef typename CDT::Constrained_edges_iterator Constrained_edges_iterator;

    typedef typename Traits::Vector_2 Vector_2;

    const Point_vector &m_vertex;

    const Traits &m_barycentric_traits;

    FT minimal_edge_length;
    FT shape_bound;

    CDT cdt;

    std::vector<Vertex_handle> Mesh_handles;



    void detect_minimal_edge_length(FT &minimal_edge_length)
    {
        const size_t number_of_vertices = m_vertex.size();

        Vector_2 edge(m_vertex[0], m_vertex[1]);
        FT edge_length = edge.squared_length();
        minimal_edge_length = std::sqrt(edge_length);
        for (size_t i = 1; i < number_of_vertices; i++)
        {
            Vector_2 edge(m_vertex[i % number_of_vertices], m_vertex[(i + 1) % number_of_vertices]);
            FT edge_length = edge.squared_length();
            edge_length = std::sqrt(edge_length);

            if (edge_length < minimal_edge_length || i == 0)
                minimal_edge_length = edge_length;
        }
    }

    void insert_constraint(CDT &cdt)
    {
        const size_t number_of_vertices = m_vertex.size();

        for (size_t i = 0; i < number_of_vertices; ++i)
        {
            size_t ip = (i + 1) % number_of_vertices;

            Point_2 sampled_boundary_point_a;
            Point_2 sampled_boundary_point_b;

            Vector_2 polygon_edge = Vector_2(m_vertex[i], m_vertex[ip]);
            FT edge_length = polygon_edge.squared_length();
            edge_length = std::sqrt(edge_length);

            size_t step_num = edge_length / shape_bound;
            for (size_t num = 0; num < step_num; num++)
            {
                FT a = (FT)num / (FT)step_num;
                FT b = (FT)(num + 1) / (FT)step_num;
                FT a_x = m_vertex[i].x() * a + m_vertex[ip].x() * (1 - a);
                FT a_y = m_vertex[i].y() * a + m_vertex[ip].y() * (1 - a);
                FT b_x = m_vertex[i].x() * b + m_vertex[ip].x() * (1 - b);
                FT b_y = m_vertex[i].y() * b + m_vertex[ip].y() * (1 - b);
                sampled_boundary_point_a = Point_2(a_x, a_y);
                sampled_boundary_point_b = Point_2(b_x, b_y);

                cdt.insert(sampled_boundary_point_a);

                Vertex_handle va = cdt.insert(CDT_Point(sampled_boundary_point_a));
                Vertex_handle vb = cdt.insert(CDT_Point(sampled_boundary_point_b));

                cdt.insert_constraint(va, vb);
            }
        }
    }

    void create_denaulay_mesh()
    {
        const size_t number_of_vertices = m_vertex.size();
        std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
        Mesher mesher(cdt);
        mesher.set_criteria(Criteria(0.125, shape_bound));
        mesher.refine_mesh();
        CGAL::make_conforming_Delaunay_2(cdt);
        CGAL::make_conforming_Gabriel_2(cdt);
        //CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, shape_bound));
        //CGAL::lloyd_optimize_mesh_2(cdt, CGAL::parameters::max_iteration_number = 10);
        std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
    }

    void list_all_vertices(CDT &cdt, Point_vector &mesh_vertices)
    {
        int i = 0;
        for (Vertex_iterator vertex_handle = cdt.finite_vertices_begin(); vertex_handle != cdt.finite_vertices_end(); ++vertex_handle)
        {
            if(!cdt.is_infinite(vertex_handle))
            {
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
        Vertex query_vertex = *query_handle;
        Point_2 query_point = query_vertex.point();

        Vertex_circulator all_neighbors_begin = cdt.incident_vertices(query_handle);
        Vertex_circulator all_neighbors_end = all_neighbors_begin;
        Vertex_handle first_neighbor_handle = all_neighbors_begin;
        Vertex first_neighbor_vertex = *first_neighbor_handle;
        Point_2 first_neighbor_point = first_neighbor_vertex.point();
        Point_2 mid_point((query_point.x() + first_neighbor_point.x()) / FT(2), (query_point.y() + first_neighbor_point.y()) / FT(2));
        if (!cdt.is_infinite(first_neighbor_handle) &&
            CGAL::bounded_side_2(m_vertex.begin(), m_vertex.end(), query_point, m_barycentric_traits) == CGAL::ON_BOUNDED_SIDE)
        {
            neighbors.push_back(first_neighbor_handle->info().index);
        }

        all_neighbors_begin++;
        while (all_neighbors_begin != all_neighbors_end)
        {
            Vertex_handle current_neighbor_handle = all_neighbors_begin;
            Vertex current_neighbor_vertex = *current_neighbor_handle;
            Point_2 current_neighbor_point = current_neighbor_vertex.point();
            Point_2 mid_point((query_point.x() + current_neighbor_point.x()) / FT(2), (query_point.y() + current_neighbor_point.y()) / FT(2));
            if (!cdt.is_infinite(current_neighbor_handle) &&
                CGAL::bounded_side_2(m_vertex.begin(), m_vertex.end(), query_point, m_barycentric_traits) == CGAL::ON_BOUNDED_SIDE)
            {
                neighbors.push_back(current_neighbor_handle->info().index);
            }
            all_neighbors_begin++;
        }
    }

    void list_all_boundary_vertices(CDT &cdt, std::vector<int> &boundary)
    {
        for (Constrained_edges_iterator start = cdt.constrained_edges_begin(); start != cdt.constrained_edges_end(); ++start)
        {
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

        for (size_t i = 0; i < mesh_vertices.size(); ++i)
        {
            f << "v " << mesh_vertices[i] << " 0.000000" << std::endl;
        }

        f << std::endl;

        for (Finite_faces_iterator fit = cdt.finite_faces_begin();
             fit != cdt.finite_faces_end();
             ++fit)
        {
            Face face = *fit;
            Vertex_handle first_vertex = face.vertex(0);
            Vertex_handle second_vertex = face.vertex(1);
            Vertex_handle third_vertex = face.vertex(2);

            Vertex v1 = *first_vertex;
            Vertex v2 = *second_vertex;
            Vertex v3 = *third_vertex;
            Point_2 p1 = v1.point();
            Point_2 p2 = v2.point();
            Point_2 p3 = v3.point();

            FT q_x = (p1.x() + p2.x() + p3.x()) / FT(3);
            FT q_y = (p1.y() + p2.y() + p3.y()) / FT(3);
            Point_2 query_point(q_x, q_y);
            if (CGAL::bounded_side_2(m_vertex.begin(), m_vertex.end(), query_point, m_barycentric_traits) == CGAL::ON_BOUNDED_SIDE)
            {
                f << "f "
                  << first_vertex->info().index + 1 << " "
                  << second_vertex->info().index + 1 << " "
                  << third_vertex->info().index + 1
                  << std::endl;
            }
        }
    }
};

}  // namespace Barycentric_coordinates

}  // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_HARMONIC_MESH_H
