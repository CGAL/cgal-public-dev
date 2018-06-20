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
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Delaunay_mesher_2.h>



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

    /// @}

    // \name Creation
    Harmonic_mesh_2(const std::vector<typename Traits::Point_2> &vertices, const Traits &b_traits) :
        vertex(vertices),
        barycentric_traits(b_traits),
        number_of_vertices(vertex.size()),
        squared_distance_2(barycentric_traits.compute_squared_distance_2_object())
    {
        // Initialize some private parameters here.
    }

    // This function computes prior functions for each query point.
    void compute_prior_functions(typename Traits::Point_2 &query_point, FT_vector &m)
    {
        // Call MEC1 prior function.
        compute_prior_functions_type_one(query_point, m);

    }

    void print_information()
    {
        //std::cout<<"Prior class function available."<<std::endl;
    }

private:

    typedef typename Traits::Vector_2 Vector_2;
    //typedef typename std::vector<Vector_2> Vector_vector;
    typedef typename std::vector<Point_2> Point_vector;

    const Point_vector &vertex;

    const Traits &barycentric_traits;

    const size_t number_of_vertices;

    typename Traits::Compute_squared_distance_2 squared_distance_2;

    // Some details and private compute functions
    void compute_prior_functions_type_one(typename Traits::Point_2 &query_point, FT_vector &m)
    {
        FT_vector r(number_of_vertices), e(number_of_vertices), ro(number_of_vertices);
        FT PItilde = 0.0;

        for (int i = 0; i < number_of_vertices; ++i ) {
            Vector_2 r_vector, e_vector;
            size_t ip = (i + 1) % number_of_vertices;

            r[i] = static_cast<FT >(sqrt(CGAL::to_double(squared_distance_2(vertex[i], query_point))) );
            e[i] = static_cast<FT >(sqrt(CGAL::to_double(squared_distance_2(vertex[ip], vertex[i]))) );
        }


        for (int i = 0; i < number_of_vertices; ++i ) {
            size_t ip = (i + 1) % number_of_vertices;
            ro[i] = r[i] + r[ip] - e[i];
        }

        for (int i = 0; i < number_of_vertices; ++i ) {
            size_t im = (i + number_of_vertices - 1) % number_of_vertices;
            FT denom = ro[im] * ro[i];

            m[i] = 1.0 / denom;
            PItilde += m[i];
        }

        for (int i = 0; i < number_of_vertices; ++i)
            m[i] /= PItilde;

    }
};

}  // namespace Barycentric_coordinates

}  // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_HARMONIC_MESH_H
