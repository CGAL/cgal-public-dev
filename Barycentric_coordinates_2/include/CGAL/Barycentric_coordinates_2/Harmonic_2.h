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

// Barycentric coordinates headers.
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>

// Boost headers.
#include <boost/optional/optional.hpp>

// Mesh Weights and Solver headers
#include <CGAL/Barycentric_coordinates_2/Harmonic_2/Harmonic_mesh.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_2/Harmonic_weights.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_2/Harmonic_solver.h>




// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Brief introduction about Harmonic coordinates




template<class Traits, class Mesh/*, class Weights, class Solver*/ >
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
        mesher(Mesh(vertices, barycentric_traits)),
        is_sparse_mesh_created(false),
        is_dense_mesh_created(false)
        //interpolator(Weights(barycentric_traits)),
        //solver(Solver(vertices, barycentric_traits))
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

    // This function computes Maximum Entropy barycentric coordinates for a chosen query point on the unbounded side of an arbitrary polygon.
    // Due to the constraint of Maximum Entropy coordinate, we can not compute coordinates for unbounded side points.
    // We keep the interface and leave this function empty except an assertion.
    template<class OutputIterator>
        inline boost::optional<OutputIterator> coordinates_on_unbounded_side(const Point_2 &query_point, OutputIterator &output, const Type_of_algorithm type_of_algorithm, const bool warning_tag = true)
    {
        //assertion
        std::cout<<"Currently compute coordinates on unbounded side is not available."<<std::endl;
    }

    // Information Functions

    // This function prints some information about maximum entropy coordinates.
    void print_coordinates_information(std::ostream &output_stream) const
    {
        return print_coordinates_information_2(output_stream);
    }

private:

    // Some convenient typedefs.
    typedef typename Traits::Vector_2              Vector_2;
    typedef typename std::vector<FT>               FT_vector;
    typedef typename std::vector<Point_2>          Point_vector;

    // Internal global variables.
    const Point_vector &vertex;

    const Traits &barycentric_traits;

    const size_t number_of_vertices;

    bool is_sparse_mesh_created;
    bool is_dense_mesh_created;

    // Mesh class
    Mesh mesher;

    // Weights class
    //Weights interpolator;

    // Solver class
    //Solver solver;

    template<class OutputIterator>
        boost::optional<OutputIterator> coordinates_on_bounded_side_precise_2(const Point_2 &query_point, OutputIterator &output)
    {
        // First check whether the partition has been created.
        // This part could be refined as singleton or static behavior. We will improve that later.
        if(!is_dense_mesh_created){
            // Here we set up a dense constraint for dense partion, the max edge length would be less than polygon_scale*dense_partition_constraint.
            FT dense_partition_constraint = FT(1)/FT(20);
            mesher.create_mesh(dense_partition_constraint);

            is_dense_mesh_created = true;
        }

        //// Locate query_point in the created partition, return one single point (perfect condition) or three triangle vertices (interpolate coordinates).
        //Point_vector location = mesher.locate_point(query_point);

        //switch (location.size()) {
        //    // query_point perfectly locates on location[0].
        //    case 1:
        //    FT_vector coordinates = mesher.get_coordinates(location[0]);
        //    break;

        //    // query_point locates inside a triangle.
        //    case 3:
        //    FT_vector neighbor1 = mesher.get_coordinates(location[0]);
        //    FT_vector neighbor2 = mesher.get_coordinates(location[1]);
        //    FT_vector neighbor3 = mesher.get_coordinates(location[2]);
        //    FT_vector coordinates = interpolator.interpolate(neighbor1, neighbor2, neighbor3, location);
        //    break;
        //}

        //for(size_t i = 0; i < number_of_vertices; ++i) {
        //    *output = coordinates[0];
        //    ++output;
        //}

        return boost::optional<OutputIterator>(output);
    }

    template<class OutputIterator>
        boost::optional<OutputIterator> coordinates_on_bounded_side_fast_2(const Point_2 &query_point, OutputIterator &output)
    {
        // First check whether the partition has been created.
        // This part could be refined as singleton or static behavior. We will improve that later.
        if(!is_sparse_mesh_created){
            // Here we set up a sparse constraint for sparse partion, the max edge length would be less than polygon_scale*sparse_partition_constraint.
            FT sparse_partition_constraint = FT(1)/FT(5);
            mesher.create_mesh(sparse_partition_constraint);

            is_sparse_mesh_created = true;
        }

        //// Locate query_point in the created partition, return one single point (perfect condition) or three triangle vertices (interpolate coordinates).
        //Point_vector location = mesher.locate_point(query_point);

        //switch (location.size()) {
        //    // query_point perfectly locates on location[0].
        //    case 1:
        //    FT_vector coordinates = mesher.get_coordinates(location[0]);
        //    break;

        //    // query_point locates inside a triangle.
        //    case 3:
        //    FT_vector neighbor1 = mesher.get_coordinates(location[0]);
        //    FT_vector neighbor2 = mesher.get_coordinates(location[1]);
        //    FT_vector neighbor3 = mesher.get_coordinates(location[2]);
        //    FT_vector coordinates = interpolator.interpolate(neighbor1, neighbor2, neighbor3, location);
        //    break;
        //}

        //for(size_t i = 0; i < number_of_vertices; ++i) {
        //    *output = coordinates[0];
        //    ++output;
        //}

        return boost::optional<OutputIterator>(output);
    }

    // OTHER FUNCTIONS.

    // Print some information about discrete harmonic coordinates.
    void print_coordinates_information_2(std::ostream &output_stream) const
    {
        // More detail information about Maximum Entropy coordinates.
        //output_stream << std::endl << "THIS FUNCTION IS UNDER CONSTRUCTION." << std::endl << std::endl;
    }

};

} // namespace Barycentric_coordinates

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_HARMONIC_2_H
