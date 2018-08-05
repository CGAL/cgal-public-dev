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
  \file Maximum_entropy_parameters.h
*/


#ifndef CGAL_MAXIMUM_ENTROPY_PARAMETERS_H
#define CGAL_MAXIMUM_ENTROPY_PARAMETERS_H

#include <CGAL/license/Barycentric_coordinates_2.h>

#include <CGAL/disable_warnings.h>

// CGAL headers.
#include <CGAL/assertions.h>

// Barycentric coordinates headers.
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>


namespace CGAL {

namespace Barycentric_coordinates{

/*!
 * \ingroup PkgBarycentric_coordinates_2
 * The class `Maximum_entropy_parameters` implements different parameters including `max_number_of_iterations` and `tolerance` for `Maximum_entropy_solver` class.
 * This class is parameterized by a traits class `Traits` and completed `Parameters` for `Maximum_entropy_2` class.

\tparam Traits must be a model of the concept `BarycentricTraits_2`.

\cgalModels `MaximumEntropyParameters`

*/

template<class Traits >
    class Maximum_entropy_parameters
{

public:

    /// \name Types
    /// @{

    /// Number type.
    typedef typename Traits::FT         FT;

    /// @}

    // \name Creation
    Maximum_entropy_parameters(const Traits &barycentric_traits) :
        m_barycentric_traits(barycentric_traits)
    {
        precise_max_iter_num = 1000;
        fast_max_iter_num = 500;

        precise_tolerance = FT(1)/FT(1000000000000); // 1e-12
        fast_tolerance = FT(1)/FT(1000000); //1e-6
    }

    Maximum_entropy_parameters(const size_t max_iter_num, const FT tolerance, const Traits &barycentric_traits) :
        m_barycentric_traits(barycentric_traits)
    {
        precise_max_iter_num = max_iter_num;
        fast_max_iter_num = max_iter_num;

        precise_tolerance = tolerance;
        fast_tolerance = tolerance;
    }

    size_t max_number_of_iterations(const Type_of_algorithm type_of_algorithm)
    {
        switch (type_of_algorithm)
        {
            case PRECISE :
            return precise_max_iter_num;

            case FAST :
            return fast_max_iter_num;

            default:
            break;
        }
    }

    FT tolerance(const Type_of_algorithm type_of_algorithm)
    {
        switch (type_of_algorithm)
        {
            case PRECISE :
            return precise_tolerance;

            case FAST :
            return fast_tolerance;

            default:
            break;
        }
    }

    void set_max_number_of_iterations(const size_t max_num_iter_)
    {
        precise_max_iter_num = max_num_iter_;
        fast_max_iter_num = max_num_iter_;
    }

    void set_tolerance(const FT tolerance_)
    {
        precise_tolerance = tolerance_;
        fast_tolerance = tolerance_;
    }


private:

    const Traits &m_barycentric_traits;

    size_t precise_max_iter_num;
    size_t fast_max_iter_num;

    FT precise_tolerance;
    FT fast_tolerance;





};

} // namespace Barycentric_coordinates

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MAXIMUM_ENTROPY_PARAMETERS_H
