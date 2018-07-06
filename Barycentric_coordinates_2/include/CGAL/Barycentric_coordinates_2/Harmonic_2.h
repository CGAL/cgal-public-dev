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
  \file Harmonic_2.h
*/

#ifndef CGAL_HARMONIC_2_H
#define CGAL_HARMONIC_2_H

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
#include <boost/tuple/tuple.hpp>

// Mesh Weights and Solver headers
#include <CGAL/Barycentric_coordinates_2/Harmonic_2/Harmonic_mesh.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_2/Harmonic_interpolator.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_2/Harmonic_solver.h>




// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Brief introduction about Harmonic coordinates




template<class Traits, class Mesh, class Interpolator, class Solver >
    class Harmonic_2
{

public:

    /// \name Types
    /// @{

    /// Number type.
    typedef typename Traits::FT      FT;

    /// Point type.
    typedef typename Traits::Point_2 Point_2;

    /// @}



    // \name Creation

    // Brief introduction of Maximum_entropy_2 class, its constructor, input and output etc.
    Harmonic_2(const std::vector<typename Traits::Point_2> &vertices, const Traits &b_traits) :
        vertex(vertices),
        barycentric_traits(b_traits),
        number_of_vertices(vertex.size()),
        dense_mesher(Mesh(vertices, barycentric_traits)),
        sparse_mesher(Mesh(vertices, barycentric_traits)),
        fast_solver(Solver(vertices, barycentric_traits)),
        precise_solver(Solver(vertices, barycentric_traits)),
        interpolator(Interpolator(barycentric_traits)),
        is_sparse_mesh_created(false),
        is_dense_mesh_created(false)
    {
        // Initialize some private parameters here.

    }

    // Computation of Maximum Entropy Weight Functions, to keep this interface the same with other coordinates.

    // This function computes weights for single query point, but in this case, Maximum Entropy coordinate can
    // not provide such weights.
    // We keep this interface, leave the content empty except an assertion.
    template<class OutputIterator>
        inline boost::optional<OutputIterator> weights(const Point_2 &query_point, OutputIterator &output)
    {
        //assertion
        std::cout<<"Currently this function is not available."<<std::endl;
    }

    // Computation of Maximum Entropy Basis Functions

    // This function computes Maximum Entropy barycentric coordinates for a chosen query point on the bounded side of an arbitrary polygon.
    // \pre The provided polygon is arbitrary one.
    template<class OutputIterator>
        inline boost::optional<OutputIterator> coordinates_on_bounded_side(const Point_2 &query_point, OutputIterator &output, const Type_of_algorithm type_of_algorithm)
    {
        switch(type_of_algorithm)
        {
            case PRECISE:
            return coordinates_on_bounded_side_precise_2(query_point, output);
            break;

            case FAST:
            return coordinates_on_bounded_side_fast_2(query_point, output);
            break;
        }

        // Pointer cannot be here. Something went wrong.
        const bool type_of_algorithm_failure = true;
        CGAL_postcondition( !type_of_algorithm_failure );
        if(!type_of_algorithm_failure) return boost::optional<OutputIterator>(output);
        else return boost::optional<OutputIterator>();
    }

    template<class OutputIterator>
        inline boost::optional<OutputIterator> coordinates_on_unbounded_side(const Point_2 &query_point, OutputIterator &output, const Type_of_algorithm type_of_algorithm, const bool warning_tag = true)
    {
        //assertion
        std::cout<<"Currently compute coordinates on unbounded side is not available."<<std::endl;
    }

    // Information Functions

    // This function prints some information about harmonic coordinates.
    void print_coordinates_information(std::ostream &output_stream) const
    {
        return print_coordinates_information_2(output_stream);
    }

private:

    // Some convenient typedefs.
    typedef typename Traits::Vector_2              Vector_2;
    typedef typename std::vector<FT>               FT_vector;
    typedef typename std::vector<int>              Index_vector;
    typedef typename std::vector<Point_2>          Point_vector;

    // Nametype for mesh vertex storation, type is <index, location, neighbor, coordinates, boundary>
    typedef boost::tuple<int, Point_2, Index_vector, FT_vector, bool> Indexed_mesh_vertex;

    std::vector<Indexed_mesh_vertex> dense_mesh_vertices;
    std::vector<Indexed_mesh_vertex> sparse_mesh_vertices;


    // Internal global variables.
    const Point_vector &vertex;

    const Traits &barycentric_traits;

    const size_t number_of_vertices;

    bool is_sparse_mesh_created;
    bool is_dense_mesh_created;

    // Mesh class
    Mesh dense_mesher;
    Mesh sparse_mesher;

    // Solver class
    Solver fast_solver;
    Solver precise_solver;

    // Interpolator class
    Interpolator interpolator;

    template<class OutputIterator>
        boost::optional<OutputIterator> coordinates_on_bounded_side_precise_2(const Point_2 &query_point, OutputIterator &output)
    {
        // First check whether the partition has been created.
        // This part could be refined as singleton or static behavior. We will improve that later.
        if(!is_dense_mesh_created){
            // Here we set up a dense constraint for dense partion, the max edge length should be less than polygon_scale*dense_partition_constraint.
            FT dense_partition_constraint = FT(1)/FT(20);
            dense_mesher.create_mesh(dense_partition_constraint);

            // Compute harmonic coordinates at each mesh vertices, store them in the property map.
            compute_harmonic_coordinates(dense_mesher, precise_solver, dense_mesh_vertices);

            is_dense_mesh_created = true;
        }

        // Locate query_point in the created partition, return three triangle vertices (interpolate coordinates).
        Point_2 query = query_point;
        std::vector<int> triangle_indices = dense_mesher.locate_point(query);

        std::vector<Indexed_mesh_vertex> indexed_triangle_vertices;
        for(size_t i = 0; i < 3; ++i)
            indexed_triangle_vertices.push_back(dense_mesh_vertices[triangle_indices[i]]);


        Point_vector triangle_vertices;
        for(size_t i = 0; i < 3; ++i)
            triangle_vertices.push_back(indexed_triangle_vertices[i].get<1>());

        FT_vector triangle_coordinates = interpolator.interpolate(triangle_vertices, query_point);

        FT_vector coordinates;
        FT C(0);
        for(size_t i = 0; i < number_of_vertices; ++i) {
            FT c(0);
            for(size_t j = 0; j < 3; ++j) {
                c += triangle_coordinates[j] * indexed_triangle_vertices[j].get<3>()[i];
            }
            coordinates.push_back(c);
            C += c;
        }

        for(size_t i = 0; i < number_of_vertices; ++i) {
            *output = coordinates[i]/C;
            ++output;
        }

        return boost::optional<OutputIterator>(output);
    }

    template<class OutputIterator>
        boost::optional<OutputIterator> coordinates_on_bounded_side_fast_2(const Point_2 &query_point, OutputIterator &output)
    {
        // First check whether the partition has been created.
        // This part could be refined as singleton or static behavior. We will improve that later.
        if(!is_sparse_mesh_created){
            // Here we set up a dense constraint for dense partion, the max edge length should be less than polygon_scale*dense_partition_constraint.
            FT sparse_partition_constraint = FT(1)/FT(2);
            sparse_mesher.create_mesh(sparse_partition_constraint);

            // Compute harmonic coordinates at each mesh vertices, store them in the property map.
            compute_harmonic_coordinates(sparse_mesher, fast_solver, sparse_mesh_vertices);

            is_sparse_mesh_created = true;
        }

        // Locate query_point in the created partition, return three triangle vertices (interpolate coordinates).
        Point_2 query = query_point;
        std::vector<int> triangle_indices = sparse_mesher.locate_point(query);

        std::vector<Indexed_mesh_vertex> indexed_triangle_vertices;
        for(size_t i = 0; i < 3; ++i)
            indexed_triangle_vertices.push_back(sparse_mesh_vertices[triangle_indices[i]]);


        Point_vector triangle_vertices;
        for(size_t i = 0; i < 3; ++i)
            triangle_vertices.push_back(indexed_triangle_vertices[i].get<1>());

        FT_vector triangle_coordinates = interpolator.interpolate(triangle_vertices, query_point);

        FT_vector coordinates;
        FT C(0);
        for(size_t i = 0; i < number_of_vertices; ++i) {
            FT c(0);
            for(size_t j = 0; j < 3; ++j) {
                c += triangle_coordinates[j] * indexed_triangle_vertices[j].get<3>()[i];
            }
            coordinates.push_back(c);
            C += c;
        }

        for(size_t i = 0; i < number_of_vertices; ++i) {
            *output = coordinates[i]/C;
            ++output;
        }

        return boost::optional<OutputIterator>(output);
    }

    // OTHER FUNCTIONS.

    // Print some information about harmonic coordinates.
    void print_coordinates_information_2(std::ostream &output_stream) const
    {

    }

    void compute_harmonic_coordinates(Mesh &mesher, Solver &solver, std::vector<Indexed_mesh_vertex> &all_mesh_vertices)
    {
        /// 1. Store all mesh vertices and index at tuple_property_map all_mesh_vertices;
        Point_vector mesh_vertices = mesher.get_all_vertices();
        all_mesh_vertices.resize(mesh_vertices.size());
        for(size_t i = 0; i < mesh_vertices.size(); ++i)
        {
            all_mesh_vertices[i].get<0>() = i;
            all_mesh_vertices[i].get<1>() = mesh_vertices[i];
            all_mesh_vertices[i].get<4>() = false;
        }

        /// 2. Get 1-Ring neighbor's location, store neighbors' index at tuple_property_map;
        for(size_t i = 0; i < mesh_vertices.size(); ++i)
        {
            std::vector<int> neighbors = mesher.get_neighbor(i);

            all_mesh_vertices[i].get<2>().resize(neighbors.size());
            all_mesh_vertices[i].get<2>() = neighbors;
        }

        /// 3. Pass the neighbor connectivity and all mesh vertex location to Solver class, solve and store the coordinates at each mesh vertices.
        std::vector<int> boundary_id = mesher.get_boundary_vertices();
        for(size_t i = 0; i < boundary_id.size(); ++i)
        {
            all_mesh_vertices[boundary_id[i]].get<4>() = true;
        }
        std::vector<bool> is_on_boundary_info;
        for(size_t i = 0; i < all_mesh_vertices.size(); ++i)
        {
            is_on_boundary_info.push_back(all_mesh_vertices[i].get<4>());
        }

        solver.set_mesh(mesh_vertices);
        solver.compute_boundary_coordinates(is_on_boundary_info);

        for(size_t i = 0; i < mesh_vertices.size(); ++i)
        {
            solver.set_connection(all_mesh_vertices[i].get<0>(), all_mesh_vertices[i].get<2>());
        }

        solver.solve();

        for(size_t i = 0; i <mesh_vertices.size(); ++i)
        {
            all_mesh_vertices[i].get<3>() = solver.get_coordinates(i);
        }


    }

};

} // namespace Barycentric_coordinates

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_HARMONIC_2_H
