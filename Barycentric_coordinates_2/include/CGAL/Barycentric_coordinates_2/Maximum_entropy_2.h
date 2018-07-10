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
  \file Maximum_entropy_2.h
*/

#ifndef CGAL_MAXIMUM_ENTROPY_2_H
#define CGAL_MAXIMUM_ENTROPY_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

#include <CGAL/disable_warnings.h>

// CGAL headers.
#include <CGAL/assertions.h>
#include <CGAL/Polygon_2_algorithms.h>

// Barycentric coordinates headers.
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>

// Boost headers.
#include <boost/optional/optional.hpp>

// Eigen headers.
#include <CGAL/Eigen_vector.h>
#include <CGAL/Eigen_matrix.h>

// Property map headers.
#include <CGAL/property_map.h>

// Solver and Prior headers
#include <CGAL/Barycentric_coordinates_2/Maximum_entropy_2/Maximum_entropy_solver.h>
#include <CGAL/Barycentric_coordinates_2/Maximum_entropy_2/Maximum_entropy_prior_function.h>



// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Brief introduction about Maximum Entropy coordinates




template<class Traits, class Prior, class Solver, class Element, class Point_map >
    class Maximum_entropy_2
{

public:

    /// \name Types
    /// @{

    /// Number type.
    typedef typename Traits::FT      FT;

    /// Point type.
    typedef typename Traits::Point_2 Point_2;

    /// Element type.
    typedef std::vector<Element> Element_range;

    /// @}



    // \name Creation

    // Brief introduction of Maximum_entropy_2 class, its constructor, input and output etc.
    Maximum_entropy_2(const Element_range &elements, const Point_map &point_map, const Traits &b_traits) :
        //vertex(vertices),
        m_elements(elements),
        m_point_map(point_map),
        barycentric_traits(b_traits),
        number_of_vertices(m_elements.size()),
        prior(Prior(m_elements, m_point_map, barycentric_traits)),
        solver(Solver(m_elements, m_point_map, barycentric_traits))
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
    typedef typename Traits::Vector_2 Vector_2;
    typedef typename std::vector<FT>      FT_vector;
    typedef typename std::vector<Point_2> Point_vector;
    //typedef typename CGAL::Eigen_solver_traits<> Eigen_solver;
    typedef typename CGAL::Eigen_matrix<FT>        Matrix;

    // Internal global variables.
    const Element_range m_elements;

    const Point_map m_point_map;

    const Traits &barycentric_traits;

    const size_t number_of_vertices;

    //Point_vector &vertex;

    // Prior class
    Prior prior;

    // Solver class
    Solver solver;

    template<class OutputIterator>
        boost::optional<OutputIterator> coordinates_on_bounded_side_precise_2(const Point_2 &query_point, OutputIterator &output)
    {
        // Implementation of precise mec computing.
        // For precise edition, we set up smaller tolerance or more max iteration steps in solver.solve() function.
        Vector_2 s;
        Point_2 query_point_2 = query_point;

        Matrix vtilde(number_of_vertices, 2);

        for(size_t i = 0; i < number_of_vertices; ++i) {
            s = Vector_2(get(m_point_map, m_elements[i]),query_point);

            vtilde.set(i, 0, s.x());
            vtilde.set(i, 1, s.y());
        }

        FT_vector m(number_of_vertices);
        prior.compute_prior_functions(query_point_2, m);


        FT_vector lambda(2);
        solver.solve(lambda, vtilde, m, PRECISE);

        FT_vector z(number_of_vertices);
        FT Z(0);
        for(size_t i = 0; i < number_of_vertices; ++i) {
            z[i] = partition(vtilde, m, lambda, (int) i);
            Z += z[i];
        }

        for(size_t i = 0; i < number_of_vertices; ++i) {
            *output = z[i] / Z;
            output++;
        }

    }

    template<class OutputIterator>
        boost::optional<OutputIterator> coordinates_on_bounded_side_fast_2(const Point_2 &query_point, OutputIterator &output)
    {
        // Implementation of fast mec computing.
        // For fast edition, we set up larger tolerance and less iteration steps in solver.solve() function.
        Vector_2 s;
        Point_2 query_point_2 = query_point;

        Matrix vtilde(number_of_vertices, 2);

        for(size_t i = 0; i < number_of_vertices; ++i) {
            s = Vector_2(get(m_point_map, m_elements[i]),query_point);

            vtilde.set(i, 0, s.x());
            vtilde.set(i, 1, s.y());
        }

        FT_vector m(number_of_vertices);
        prior.compute_prior_functions(query_point_2, m);


        FT_vector lambda(2);

        solver.solve(lambda, vtilde, m, FAST);


        FT_vector z(number_of_vertices);
        FT Z(0);
        for(size_t i = 0; i < number_of_vertices; ++i) {
            z[i] = partition(vtilde, m, lambda, (int) i);
            Z += z[i];
            //std::cout<<"z["<<i<<"] : "<<z[i]<<std::endl;
        }

        for(size_t i = 0; i < number_of_vertices; ++i) {
            *output = z[i] / Z;
            output++;
        }

        return boost::optional<OutputIterator>(output);

    }

    // OTHER FUNCTIONS.

    // Print some information about discrete harmonic coordinates.
    void print_coordinates_information_2(std::ostream &output_stream) const
    {
        // More detail information about Maximum Entropy coordinates.
        //output_stream << std::endl << "THIS FUNCTION IS UNDER CONSTRUCTION." << std::endl << std::endl;
    }

    // Compute partition values using dot product of matrix vtilde and vector lambda.
    inline FT partition(const Matrix &vtilde, const FT_vector &m, const FT_vector &lambda, const int index) const {
        assert(index >= 0);
        FT dot_product = lambda.at(0) * vtilde(index, 0) + lambda.at(1) * vtilde(index, 1);

        FT exponent = static_cast<FT >(exp(CGAL::to_double(-dot_product)) );

        return m[index] * exponent;
    }


};

} // namespace Barycentric_coordinates

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MAXIMUM_ENTROPY_2_H
