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

// Add partition headers.
#include <CGAL/Barycentric_coordinates_2/Maximum_entropy_2/Partition.h>

namespace CGAL {

namespace Barycentric_coordinates{

/*!
 * \ingroup PkgBarycentric_coordinates_2
 * The class `Maximum_entropy_newton_solver` implements numerical newton solver for maximum entropy coordinates ( \cite cgal:bc:hs-mecap-08 ).
 * This class is parameterized by a traits class `Traits`.

\tparam Traits must be a model of the concept `BarycentricTraits_2`.

\cgalModels `MaximumEntropySolver`

*/

template<class Traits >
    class Maximum_entropy_newton_solver
{

public:

    /// \name Types
    /// @{

    /// Number type.
    typedef typename Traits::FT         FT;
    typedef typename std::vector<FT>    FT_vector;

    /// Point type.
    typedef typename Traits::Point_2 Point_2;

    /// Matrix type.
    typedef typename CGAL::Eigen_matrix<FT> Matrix;


    /// @}

    // \name Creation
    Maximum_entropy_newton_solver(const std::vector<typename Traits::Point_2> &vertices, const Traits &barycentric_traits) :
        m_vertex(vertices),
        m_barycentric_traits(barycentric_traits),
        partition(Partition())
    {

    }

    // Main function, solve the Newton iteration problem with a user determined type_of_algorithm(max_num_iter and tol).
    void solve(FT_vector &lambda, const Matrix &vtilde, const FT_vector &m, const size_t max_number_iter, const FT tol)
    {
        const FT alpha = FT(1);

        for (size_t iter = 0; iter < max_number_iter; ++iter) {
            FT_vector g(2);
            compute_gradient(lambda, vtilde, m, g);

            const FT g_norm = g[0] * g[0] + g[1] * g[1];
            if (g_norm < tol) break;

            Matrix H(2, 2);
            compute_hessian(lambda, vtilde, m, H);

            FT_vector delta_lambda(2);

            solve_linear_system(g, H, delta_lambda);

            lambda[0] = lambda[0] + alpha * delta_lambda[0];
            lambda[1] = lambda[1] + alpha * delta_lambda[1];
        }
    }

private:

    typedef typename std::vector<Point_2> Point_vector;

    typedef typename CGAL::Eigen_solver_traits<> Eigen_solver;
    typedef typename CGAL::Eigen_vector<FT>    Vector;

    // Internal global variables.
    const Point_vector &m_vertex;

    const Traits &m_barycentric_traits;

    Partition partition;

    // Implement details.
    // Compute first derivative.
    inline void compute_gradient(const FT_vector &lambda, const Matrix &vtilde, const FT_vector &m, FT_vector &g)
    {
        FT dZ1 = FT(0);
        FT dZ2 = FT(0);

        const size_t number_of_vertices = m_vertex.size();

        for (size_t i = 0; i < number_of_vertices; ++i) {

            const FT Zival = partition(vtilde, m, lambda, i);
            dZ1 += Zival * vtilde(i, 0);
            dZ2 += Zival * vtilde(i, 1);
        }

        g[0] = -dZ1;
        g[1] = -dZ2;
    }

    // Compute second derivative.
    inline void compute_hessian(const FT_vector &lambda, const Matrix &vtilde, const FT_vector &m, Matrix &H)
    {
        FT dZ11 = FT(0);
        FT dZ12 = FT(0);
        FT dZ22 = FT(0);

        const size_t number_of_vertices = m_vertex.size();

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

    // Solve a 2x2 linear system.
    inline void solve_linear_system(const FT_vector &g, const Matrix &H, FT_vector &delta_lambda)
    {
        // Closed-form solver, may be not stable. (Deprecated)
        /*
        FT denom0 = H(0, 0);
        FT denom1 = H(1, 1) * H(0, 0) - H(1, 0) * H(0, 1);

        delta_lambda[1] = (H(1, 0) * g[0] - g[1] * H(0, 0)) / denom1;
        delta_lambda[0] = (-g[0] - H(0, 1) * delta_lambda[1]) / denom0;
        */


        // Eigen solver.
        Eigen::SparseMatrix<FT> m_H;
        m_H.resize(2,2);
        m_H.insert(0,0) = H(0,0);
        m_H.insert(1,0) = H(1,0);
        m_H.insert(0,1) = H(0,1);
        m_H.insert(1,1) = H(1,1);

        Eigen::Vector2d v_g(g[0], g[1]);

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<FT> > ldlt;
        ldlt.compute(m_H);
        Eigen::Vector2d v_delta_lambda = ldlt.solve(-v_g);
        delta_lambda.clear();
        delta_lambda.resize(2);
        delta_lambda[0] = v_delta_lambda[0];
        delta_lambda[1] = v_delta_lambda[1];
    }

};

} // namespace Barycentric_coordinates

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MAXIMUM_ENTROPY_SOLVER_H
