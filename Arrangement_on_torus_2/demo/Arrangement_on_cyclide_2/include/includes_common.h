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
// File          : demos/xsurface/include/includes_common.h
// QdX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
//============================================================================

#ifndef XSURFACE_INCLUDES_COMMON_H
#define XSURFACE_INCLUDES_COMMON_H

// various debug outputs
//#define CGAL_SL_VERBOSE 1
//#define CGAL_ARR_CONSTRUCTION_SL_VISITOR_VERBOSE 1
//#define CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE 1

//#define CGAL_ARR_TOPOLOGY_TRAITS_VERBOSE 1
//#define CGAL_ARR_SIGN_OF_SUBPATH_VERBOSE 1

#ifndef CYCLIDE_DEMO_USE_CORE
#define CYCLIDE_DEMO_USE_CORE 1 
#endif

#ifndef CGAL_FILTERED_CKvA_2
#define CGAL_FILTERED_CKvA_2 0
#endif

#include <CGAL/basic.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/algorithm.h>

#include<boost/array.hpp>
#include <iostream>

#if CYCLIDE_DEMO_USE_CORE
typedef CGAL::CORE_arithmetic_kernel Arithmetic_kernel;
#else
typedef CGAL::LEDA_arithmetic_kernel Arithmetic_kernel;
#endif

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

#endif // XSURFACE_INCLUDES_COMMON_H 
