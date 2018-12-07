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
  \file Harmonic_interpolator.h
*/

#ifndef CGAL_HARMONIC_INTERPOLATOR_H
#define CGAL_HARMONIC_INTERPOLATOR_H

#include <CGAL/license/Barycentric_coordinates_2.h>

#include <CGAL/disable_warnings.h>

// CGAL headers.
#include <CGAL/assertions.h>
#include <CGAL/Polygon_2_algorithms.h>

// Barycentric coordinates headers.
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>

// Boost headers.
#include <boost/optional/optional.hpp>

#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>


// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Introduction of Harmonic_interpolator

template<class Traits>
    class Harmonic_interpolator
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

    typedef typename CGAL::cpp11::array<FT,3> Coords;



    /// @}

    // \name Creation
    Harmonic_interpolator(const Traits &barycentric_traits) :
        m_barycentric_traits(barycentric_traits)
    {

    }

    FT_vector interpolate(Point_vector &triangle_vertices, Point_2 query_point)
    {
        FT_vector triangle_coordinates;
        triangle_coordinates.reserve(3);

        const Coords coords = CGAL::Barycentric_coordinates::compute_triangle_coordinates_2(triangle_vertices[0], triangle_vertices[1], triangle_vertices[2], query_point, Traits());

        for(size_t i = 0; i < 3; ++i)
        {
            triangle_coordinates[i] = coords[i];
        }
        return triangle_coordinates;

    }

private:

    const Traits &m_barycentric_traits;


};


}  // namespace Barycentric_coordinates

}  // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_HARMONIC_INTERPOLATOR_H
