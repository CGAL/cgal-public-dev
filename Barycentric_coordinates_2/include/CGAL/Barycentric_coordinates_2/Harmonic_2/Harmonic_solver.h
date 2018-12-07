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
  \file Harmonic_solver.h
*/
#ifndef CGAL_HARMONIC_SOLVER_H
#define CGAL_HARMONIC_SOLVER_H

#include <CGAL/license/Barycentric_coordinates_2.h>

#include <CGAL/disable_warnings.h>

// CGAL headers.
#include <CGAL/assertions.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/property_map.h>
#include <CGAL/Kernel/global_functions.h>

// Barycentric coordinates headers.
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>

// Segment coordinates headers.
#include <CGAL/Barycentric_coordinates_2/Segment_coordinates_2.h>

// Boost headers.
#include <boost/optional/optional.hpp>

// Eigen headers.
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Eigen_vector.h>




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

// Introduction of Harmonic_solver_2
template<class Traits>
    class Harmonic_solver
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

    /// Vector type.
    typedef typename Traits::Vector_2 Vector_2;

    /// Segment type.
    typedef typename Traits::Segment_2 Segment_2;

    typedef CGAL::cpp11::array<FT,2> Pair;





    /// @}

    // \name Creation
    Harmonic_solver(const std::vector<typename Traits::Point_2> &vertices, const Traits &barycentric_traits) :
        m_vertex(vertices),
        m_barycentric_traits(barycentric_traits),
        squared_distance_2(m_barycentric_traits.compute_squared_distance_2_object()),
        sqrt(Get_sqrt<Traits>::sqrt_object(m_barycentric_traits))
    {
        // Initialize some private parameters here.
    }

    void set_mesh(const Point_vector all_mesh_vertices)
    {
        mesh_vertices = all_mesh_vertices;
    }

    void compute_boundary_coordinates(const std::vector<bool> is_on_boundary_info)
    {
        boundary_info = is_on_boundary_info;
        indices.resize(mesh_vertices.size());
        int numB = 0, numI = 0;
        for(size_t i = 0; i < boundary_info.size(); ++i)
        {
            if(boundary_info[i]) indices[i] = numB++;
            else indices[i] = numI++;
        }
        const size_t number_of_vertices = m_vertex.size();
        boundary = Eigen::MatrixXd::Zero(numB, number_of_vertices);
        b = Eigen::MatrixXd::Zero(numI, number_of_vertices);
        x = Eigen::MatrixXd::Zero(numI, number_of_vertices);
        A.resize(numI, numI);
        tripletList.reserve(numI * 7);

        for(size_t i = 0; i < is_on_boundary_info.size(); ++i)
        {
            if(boundary_info[i]) {
                compute_segment_coordinates(indices[i], i);
            }
        }
    }

    void set_connection(int index, std::vector<int> neighbor_index)
    {
        if(!boundary_info[index] && neighbor_index.size() != 0)
        {
            FT_vector alphaCot(neighbor_index.size());
            FT_vector betaCot(neighbor_index.size());

            for(size_t j = 0; j < neighbor_index.size(); ++j)
            {
                const size_t jp = (j + 1) % neighbor_index.size();

                Vector_2 s1(mesh_vertices[index], mesh_vertices[neighbor_index[j]]);
                Vector_2 s2(mesh_vertices[neighbor_index[jp]], mesh_vertices[neighbor_index[j]]);
                alphaCot[j] = cotangent(s2, s1);

                Vector_2 s3(mesh_vertices[neighbor_index[j]], mesh_vertices[neighbor_index[jp]]);
                Vector_2 s4(mesh_vertices[index], mesh_vertices[neighbor_index[jp]]);
                betaCot[j] = cotangent(s4, s3);
            }

            FT W(0);
            for(size_t j = 0; j < neighbor_index.size(); ++j)
            {
                const size_t jp = (j + 1) % neighbor_index.size();
                const size_t idx = neighbor_index[jp];

                const FT w = -(alphaCot[j] + betaCot[jp]);
                W -= w;

                const size_t number_of_vertices = m_vertex.size();

                if(boundary_info[idx])
                {
                    for(size_t k = 0; k < number_of_vertices; ++k)
                        b(indices[index], k) -= boundary(indices[idx], k) * w;
                }
                else
                {
                    tripletList.push_back(T(indices[index], indices[idx], w));
                }
            }
            tripletList.push_back(T(indices[index], indices[index], W));
        }
    }

    void solve()
    {
        A.setFromTriplets(tripletList.begin(), tripletList.end());
        A.makeCompressed();

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<FT> > ldlt;
        ldlt.compute(A);
        x = ldlt.solve(b);
    }

    FT_vector get_coordinates(size_t i)
    {
        FT_vector computed_coordinates;
        const size_t number_of_vertices = m_vertex.size();
        computed_coordinates.resize(number_of_vertices);

        if(boundary_info[i])
        {
            for(size_t k = 0; k < number_of_vertices; ++k)
            {
                computed_coordinates[k] = boundary(indices[i], k);
            }
        }
        else
        {
            for(size_t k = 0; k < number_of_vertices; ++k)
            {
                computed_coordinates[k] = x(indices[i], k);
            }
        }
        return computed_coordinates;
    }

private:

    const Point_vector &m_vertex;

    const Traits &m_barycentric_traits;

    typename Traits::Compute_squared_distance_2 squared_distance_2;

    typename Get_sqrt<Traits>::Sqrt sqrt;

    Point_vector mesh_vertices;

    Eigen::MatrixXd boundary;

    Eigen::SparseMatrix<FT> A;
    Eigen::MatrixXd b;
    Eigen::MatrixXd x;

    typedef Eigen::Triplet<FT> T;
    std::vector<T> tripletList;


    std::vector<bool> boundary_info;
    std::vector<int> indices;

    void compute_segment_coordinates(int matrix_index, int boundary_index)
    {
        Point_2 boundary_point = mesh_vertices[boundary_index];
        const size_t number_of_vertices = m_vertex.size();
        for(size_t i = 0; i < number_of_vertices; i++) {
            size_t ip = (i + 1) % number_of_vertices;
            if(CGAL::are_strictly_ordered_along_line(m_vertex[i], boundary_point, m_vertex[ip]))
            {
                const Pair segment_coordinates = CGAL::Barycentric_coordinates::compute_segment_coordinates_2(m_vertex[i], m_vertex[ip], boundary_point, Traits());
                boundary(matrix_index, i) = segment_coordinates[0];
                boundary(matrix_index, ip) = segment_coordinates[1];
                break;
            }
            else if(i==number_of_vertices)
            {
                std::cout << "segment coordinates computing error for boundary points!!!" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    FT cotangent(const Vector_2 s2, const Vector_2 s1)
    {
        FT scalar_product = s2.x() * s1.x() + s2.y() * s1.y();
        FT cross_product = s2.x() * s1.y() - s2.y() * s1.x();

        FT cotangent_value = scalar_product / CGAL::abs<FT>(cross_product);
    }

};

}  // namespace Barycentric_coordinates

}  // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_HARMONIC_SOLVER_H
