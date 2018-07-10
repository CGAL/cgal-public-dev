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

// Number utils headers
#include <CGAL/number_utils.h>

// Add Eigen headers later
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Eigen_vector.h>
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Vector_2.h>

// Property map headers.
#include <CGAL/property_map.h>



// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Introduction of Maximum_entropy_prior_function_type_one_2

template<class Traits, class Element, class Point_map >
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

    /// Element type.
    typedef std::vector<Element> Element_range;


    /// @}

    // \name Creation
    Maximum_entropy_prior_function_type_one_2(const Element_range &elements, const Point_map &point_map, const Traits &barycentric_traits) :
        m_elements(elements),
        m_point_map(point_map),
        m_barycentric_traits(barycentric_traits),
        number_of_vertices(m_elements.size()),
        squared_distance_2(m_barycentric_traits.compute_squared_distance_2_object())
    {
        // Initialize some private parameters here.
    }

    // This function computes prior functions for each query point.
    template<class Range>
        void compute_prior_functions(typename Traits::Point_2 &query_point, Range &m)
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
    typedef typename std::vector<Point_2> Point_vector;

    //const Point_vector &vertex;
    const Element_range m_elements;

    const Point_map m_point_map;

    const Traits &m_barycentric_traits;

    const size_t number_of_vertices;

    typename Traits::Compute_squared_distance_2 squared_distance_2;

    // Some details and private compute functions
    template<class Range>
        void compute_prior_functions_type_one(typename Traits::Point_2 &query_point, Range &m)
    {
        FT_vector r(number_of_vertices), e(number_of_vertices), ro(number_of_vertices), output(number_of_vertices);
        FT PItilde(0);

        for (size_t i = 0; i < number_of_vertices; ++i ) {
            Vector_2 r_vector, e_vector;
            size_t ip = (i + 1) % number_of_vertices;

            r[i] = static_cast<FT >(sqrt(CGAL::to_double(squared_distance_2(get(m_point_map, m_elements[i]), query_point))) );
            e[i] = static_cast<FT >(sqrt(CGAL::to_double(squared_distance_2(get(m_point_map, m_elements[ip]), get(m_point_map, m_elements[i])))) );
        }


        for (size_t i = 0; i < number_of_vertices; ++i ) {
            size_t ip = (i + 1) % number_of_vertices;
            ro[i] = r[i] + r[ip] - e[i];
        }

        for (size_t i = 0; i < number_of_vertices; ++i ) {
            size_t im = (i + number_of_vertices - 1) % number_of_vertices;
            FT denom = ro[im] * ro[i];

            output[i] = 1.0 / denom;
            PItilde += output[i];
        }

        for (size_t i = 0; i < number_of_vertices; ++i)
            output[i] /= PItilde;

        m.clear();
        m.resize(number_of_vertices);
        for(size_t i = 0; i < number_of_vertices; ++i)
            m[i] = output[i];
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
