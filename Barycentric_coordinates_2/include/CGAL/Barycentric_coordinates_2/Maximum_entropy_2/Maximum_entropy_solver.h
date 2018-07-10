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
  \file Maximum_entropy_solver.h
*/

/*
    TODO_list:
    1. Add Eigen solver to function solve_linear_system().
    2. Test different max_num_iter and tol, find some proper value for PRECISE and FAST mode.
    3. Check if the exp() using (in function partition()) is correct, currently we receive CGAL ERROR with exact kernel.
*/

#ifndef CGAL_MAXIMUM_ENTROPY_SOLVER_H
#define CGAL_MAXIMUM_ENTROPY_SOLVER_H

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

// Property map headers.
#include <CGAL/property_map.h>

namespace CGAL {

namespace Barycentric_coordinates{

template<class Traits, class Element, class Point_map >
    class Maximum_entropy_newton_solver
{

public:

    /// \name Types
    /// @{

    /// Number type.
    typedef typename Traits::FT         FT;
    typedef typename std::vector<FT>    FT_vector;
    //typedef typename CGAL::Eigen_solver_traits<> Eigen_solver;
    typedef typename CGAL::Eigen_matrix<FT>    Matrix;
    //typedef CGAL::cpp11::array<FT, 6>   Eigen_matrix;

    /// Point type.
    typedef typename Traits::Point_2 Point_2;

    /// Element type.
    typedef std::vector<Element> Element_range;


    /// @}

    // \name Creation
    Maximum_entropy_newton_solver(const Element_range &elements, const Point_map &point_map, const Traits &barycentric_traits) :
        m_elements(elements),
        m_point_map(point_map),
        m_barycentric_traits(barycentric_traits),
        number_of_vertices(m_elements.size())
    {
        // Initialize some private parameters here.
    }

    // Main function, solve the Newton iteration problem with a user determined type_of_algorithm(max_num_iter and tol).
    template<class Range>
        void solve(Range &lambda, Matrix vtilde, Range &m, const Type_of_algorithm type_of_algorithm)
    {
        switch (type_of_algorithm)
        {
            case PRECISE :
            max_number_iter = 1000;
            tol = 1.0e-12;
            optimize_parameters(lambda, vtilde, m, max_number_iter, tol);

            case FAST :
            max_number_iter = 500;
            tol = 1.0e-6;
            optimize_parameters(lambda, vtilde, m, max_number_iter, tol);
        }
    }

private:

    typedef typename std::vector<Point_2> Point_vector;

    const Element_range m_elements;

    const Point_map m_point_map;

    const Traits &m_barycentric_traits;

    const size_t number_of_vertices;

    int max_number_iter;
    FT tol;

    template<class Range>
        void optimize_parameters(Range &lambda, Matrix vtilde, Range &m, int max_number_iter, FT tol)
    {
        const FT alpha(1);
        FT_vector output_lambda(lambda), output_m(m);
        for (size_t iter = 0; iter < max_number_iter; ++iter) {
            FT_vector g(2);
            compute_gradient(output_lambda, vtilde, output_m, g);


            const FT g_norm = g[0] * g[0] + g[1] * g[1];
            if (g_norm < tol) break;

            Matrix H(2, 2);
            compute_hessian(output_lambda, vtilde, output_m, H);

            FT_vector delta_lambda(2);
            solve_linear_system(g, H, vtilde, delta_lambda);

            output_lambda[0] = output_lambda[0] + alpha * delta_lambda[0];
            output_lambda[1] = output_lambda[1] + alpha * delta_lambda[1];
        }

        lambda.clear();
        lambda.resize(2);
        for(size_t i = 0; i < lambda.size(); ++i)
            lambda[i] = output_lambda[i];

        m.clear();
        m.resize(number_of_vertices);
        for(size_t i = 0; i < m.size(); ++i)
            m[i] = output_m[i];
    }

    // Implement details.
    // Compute first derivative.
    inline void compute_gradient(const FT_vector lambda, Matrix vtilde, FT_vector m, FT_vector &g)
    {
        FT dZ1(0), dZ2(0);
        for (size_t i = 0; i < number_of_vertices; ++i) {

            const FT Zival = partition(vtilde, m, lambda, i);
            dZ1 += Zival * vtilde(i, 0);
            dZ2 += Zival * vtilde(i, 1);
        }

        g[0] = -dZ1;
        g[1] = -dZ2;
    }

    // Compute second derivative.
    inline void compute_hessian(const FT_vector lambda, Matrix vtilde, FT_vector m, Matrix &H)
    {
        FT dZ11(0), dZ12(0), dZ22(0);
        for (size_t i = 0; i < number_of_vertices; ++i) {

            const FT Zival = partition(vtilde, m, lambda, i);

            dZ11 += Zival * vtilde(i, 0) * vtilde(i, 0);
            dZ12 += Zival * vtilde(i, 0) * vtilde(i, 1);
            dZ22 += Zival * vtilde(i, 1) * vtilde(i, 1);
        }

        H.set(0, 0, dZ11);
        H.set(0, 1, dZ12);
        H.set(1, 0, dZ12);
        H.set(1, 1, dZ22);
    }

    // Later we will add another eigen solver to this function.
    inline void solve_linear_system(const FT_vector g, const Matrix H, Matrix vtilde, FT_vector &delta_lambda)
    {
        FT denom0 = H(0, 0);
        FT denom1 = H(1, 1) * H(0, 0) - H(1, 0) * H(0, 1);


        delta_lambda[1] = (H(1, 0) * g[0] - g[1] * H(0, 0)) / denom1;
        delta_lambda[0] = (-g[0] - H(0, 1) * delta_lambda[1]) / denom0;
    }


    inline FT partition(const Matrix vtilde, const FT_vector m, const FT_vector lambda, const int index)
    {
        assert(index >= 0);
        FT dot_product = lambda[0] * vtilde(index, 0) + lambda[1] * vtilde(index, 1);

        FT exponent = static_cast<FT >(exp(CGAL::to_double(-dot_product)) );

        return m[index] * exponent;
    }
};

} // namespace Barycentric_coordinates

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MAXIMUM_ENTROPY_SOLVER_H
