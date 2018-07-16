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
  \file Partition.h
*/

#ifndef CGAL_PARTITION_H
#define CGAL_PARTITION_H

#include <CGAL/license/Barycentric_coordinates_2.h>

#include <CGAL/disable_warnings.h>

// CGAL headers.
#include <CGAL/assertions.h>

// Boost headers.
#include <boost/optional/optional.hpp>

// Add CGAL/Eigen headers.
#include <CGAL/Eigen_matrix.h>


    class Partition
{

public:

    // \name Creation
    Partition()
    {

    }

    template<typename FT>
        inline FT operator()(const typename CGAL::Eigen_matrix<FT> &vtilde, const typename std::vector<FT> &m, const typename std::vector<FT> &lambda, const size_t index)
    {
        FT dot_product = lambda[0] * vtilde(index, 0) + lambda[1] * vtilde(index, 1);

        FT exponent = static_cast<FT >(exp(CGAL::to_double(-dot_product)) );

        return m[index] * exponent;
    }

private:


};

#include <CGAL/enable_warnings.h>

#endif // CGAL_MAXIMUM_ENTROPY_SOLVER_H
