//============================================================================
//
// Copyright (c) 2001-2008 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
//----------------------------------------------------------------------------
//
// Library       : QdX
// File          : demos/xtriangulate/include/includes_common.h
// QdX_release   : $Name:  $
// Revision      : $Revision: 1.4 $
// Revision_date : $Date: 2009-07-24 13:21:30 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
//============================================================================

#ifndef XTRIANGULATE_INCLUDES_COMMON_H
#define XTRIANGULATE_INCLUDES_COMMON_H

#include <CGAL/basic.h>
#include <CGAL/algorithm.h>

#include<boost/array.hpp>
#include <iostream>
 
template <class NT, std::size_t N>
std::ostream& operator <<(std::ostream& os, 
        const boost::array<NT, N>& v) {

    CGAL::output_range(os, v.begin(), v.end(), ", ", "[", "]");
    return os;
}

typedef boost::array< float, 3 >  Point_3f;
typedef boost::array< double, 3 > Point_3d; 
typedef boost::array< double, 2 > Point_2d;

typedef std::vector< Point_3f > Point_vec_3f;
typedef std::vector< Point_3d > Point_vec_3d;

typedef boost::array< unsigned, 3 > Tri_index;
typedef std::vector< Tri_index > Triangles_vector;

#endif // XTRIANGULATE_INCLUDES_COMMON_H

