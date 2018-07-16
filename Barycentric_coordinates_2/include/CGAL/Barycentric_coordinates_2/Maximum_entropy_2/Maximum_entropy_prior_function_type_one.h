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

// Finds a square root of the provided value of the type `Kernel::FT` by first converting it to the double type and then taking the square root using the `CGAL::sqrt()` function.
template<class Traits>
    class Default_sqrt
{
    typedef typename Traits::FT FT;

public:
    FT operator()(const FT &value) const
    {
        return FT(CGAL::sqrt(CGAL::to_double(value)));
    }
};

BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Sqrt, Sqrt, false)

// Case: do_not_use_default = false.
template<class Traits, bool do_not_use_default = Has_nested_type_Sqrt<Traits>::value>
    class Get_sqrt
{
public:
    typedef Default_sqrt<Traits> Sqrt;

    static Sqrt sqrt_object(const Traits&)
    {
        return Sqrt();
    }
};

// Case: do_not_use_default = true.
template<class Traits>
    class Get_sqrt<Traits, true>
{
public:
    typedef typename Traits::Sqrt Sqrt;

    static Sqrt sqrt_object(const Traits &traits)
    {
        return traits.sqrt_object();
    }
};

// Introduction of Maximum_entropy_prior_function_type_one_2

template<class Traits >
    class Maximum_entropy_prior_function_type_one
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
    Maximum_entropy_prior_function_type_one(const std::vector<typename Traits::Point_2> &vertices, const Traits &barycentric_traits) :
        m_vertex(vertices),
        m_barycentric_traits(barycentric_traits),
        squared_distance_2(m_barycentric_traits.compute_squared_distance_2_object()),
        sqrt(Get_sqrt<Traits>::sqrt_object(m_barycentric_traits))
    {
        // Initialize some private parameters here.
        const size_t number_of_vertices = m_vertex.size();
        r.resize(number_of_vertices);
        e.resize(number_of_vertices);
        ro.resize(number_of_vertices);
        output.resize(number_of_vertices);
    }

    // This function computes prior functions for each query point.
    template<class Range>
        void compute(typename Traits::Point_2 query_point, Range &m)
    {
        FT PItilde = FT(0);

        const size_t number_of_vertices = m_vertex.size();

        for (size_t i = 0; i < number_of_vertices; ++i ) {
            Vector_2 r_vector, e_vector;
            size_t ip = (i + 1) % number_of_vertices;

            r[i] = sqrt(CGAL::to_double(squared_distance_2(m_vertex[i], query_point )) );
            e[i] = sqrt(CGAL::to_double(squared_distance_2(m_vertex[ip], m_vertex[i] )) );
        }


        for (size_t i = 0; i < number_of_vertices; ++i ) {
            size_t ip = (i + 1) % number_of_vertices;
            ro[i] = r[i] + r[ip] - e[i];
        }

        for (size_t i = 0; i < number_of_vertices; ++i ) {
            size_t im = (i + number_of_vertices - 1) % number_of_vertices;
            FT denom = ro[im] * ro[i];

            output[i] = FT(1) / denom;
            PItilde += output[i];
        }

        for (size_t i = 0; i < number_of_vertices; ++i)
            output[i] = output[i] / PItilde;

        m.clear();
        m.resize(number_of_vertices);
        for(size_t i = 0; i < number_of_vertices; ++i)
            m[i] = output[i];
    }

    void print_information()
    {
        //std::cout<<"Prior class function available."<<std::endl;
    }

private:

    typedef typename Traits::Vector_2 Vector_2;
    typedef typename std::vector<Point_2> Point_vector;

    const Point_vector &m_vertex;

    const Traits &m_barycentric_traits;

    typename Traits::Compute_squared_distance_2 squared_distance_2;

    typename Get_sqrt<Traits>::Sqrt sqrt;

    FT_vector r, e, ro, output;

};

}  // namespace Barycentric_coordinates

}  // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MAXIMUM_ENTROPY_PRIOR_FUNCTION_H
