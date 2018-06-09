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
  \file Maximum_entropy_prior_function.h
*/

#ifndef CGAL_MAXIMUM_ENTROPY_PRIOR_FUNCTION_H
#define CGAL_MAXIMUM_ENTROPY_PRIOR_FUNCTION_H

#include <CGAL/license/Barycentric_coordinates_2.h>

#include <CGAL/disable_warnings.h>

// CGAL headers.
#include <CGAL/assertions.h>
#include <CGAL/Polygon_2_algorithms.h>

// Barycentric coordinates headers.
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>

// Boost headers.
#include <boost/optional/optional.hpp>

// Add Eigen headers later
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Eigen_vector.h>
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Vector_2.h>



// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Introduction of Maximum_entropy_prior_function_type_one_2

template<class Traits>
    class Maximum_entropy_prior_function_type_one_2
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
    Maximum_entropy_prior_function_type_one_2(const std::vector<typename Traits::Point_2> &vertices, const Traits &b_traits) :
        vertex(vertices),
        barycentric_traits(b_traits),
        number_of_vertices(vertex.size()),
        squared_distance_2(barycentric_traits.compute_squared_distance_2_object())
    {
        // For many query points, we will create a Prior class just once and reuse(call) function compute_prior_functions() for each query point.
        //vertex = vertices;
        //number_of_vertices = vertex.size();
        //std::cout<<"Prior class created."<<std::endl;
    }

    // This function computes prior functions for each query point.
    void compute_prior_functions(typename Traits::Point_2 &query_point, FT_vector &m)
    {
        // Implement MEC1 prior function.
        //FT_vector m(number_of_vertices);
        //std::cout<<"compute_prior_functions called"<<std::endl;
        compute_prior_functions_type_one(query_point, m);
        //std::cout<<"compute_prior_functions finished"<<std::endl;

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
        // Vector_2 r_vector, e_vector.
        //std::cout<<"compute_prior_functions_type_one called"<<std::endl;

        FT_vector r(number_of_vertices), e(number_of_vertices), ro(number_of_vertices);
        FT PItilde = 0.0;

        //for(int i=0;i<number_of_vertices;i++){
        //    std::cout<<"vertex "<<i<<" : "<<vertex[i].x()<<" "<<vertex[i].y()<<std::endl;
        //}
        //std::cout<<"query_point "<<" : "<<query_point.x()<<" "<<query_point.y()<<std::endl;

        for (int i = 0; i < number_of_vertices; ++i ) {
            //Vector_2 r_vector, e_vector;
            size_t ip = (i + 1) % number_of_vertices;
            //std::cout<<"ip "<<ip<<std::endl;
            //std::cout<<"assertion"<<std::endl;
            //r_vector = Vector_2(vertex[i], query_point);
            r[i] = FT(sqrt(squared_distance_2(vertex[i], query_point)));
            //std::cout<<"r[i] "<<r[i]<<std::endl;
            //std::cout<<"r[i]"<<r[i]<<std::endl;
            //e_vector = Vector_2(vertex[ip], vertex[i]);
            e[i] =FT(sqrt(squared_distance_2(vertex[ip], vertex[i])));
            //std::cout<<"e[i] "<<e[i]<<std::endl;
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

        //std::cout<<"compute_prior_functions_type_one finished"<<std::endl;
    }
};


// Introduction of Maximum_entropy_prior_function_type_two_2

template<class Traits>
    class Maximum_entropy_prior_function_type_two_2
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
    Maximum_entropy_prior_function_type_two_2(const Traits &b_traits)
    {
        // For many query points, we will create a Prior class just once and reuse(call) function compute_prior_functions() for each query point.
    }

    // This function computes prior functions for each query point.
    void compute_prior_functions(typename Traits::Point_2 &query_point, const FT_vector &m)
    {
        // Implement MEC2 prior function.
    }

private:

    // Some details and private compute functions
};

}  // namespace Barycentric_coordinates

}  // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MAXIMUM_ENTROPY_PRIOR_FUNCTION_H
