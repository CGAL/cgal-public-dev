//============================================================================
//
// Copyright (c) 2001-2010 Max-Planck-Institut Saarbruecken (Germany).
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
// File          : demos/xsurface/include/includes_ckva.h
// QdX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
//============================================================================

#ifndef XSURFACE_INCLUDES_CKVA_H
#define XSURFACE_INCLUDES_CKVA_H

//!@file includes_ckva.h
//! collects CGAL::CKvA_2 related includes for arrangement computation

// TODO filtered versions?
#include <CGAL/Arr_surfaces_intersecting_dupin_cyclide_traits_2.h>


typedef Arithmetic_kernel::Integer Integer;

typedef CGAL::Arr_surfaces_intersecting_dupin_cyclide_traits_2< Integer > 
Geo_traits;

typedef Geo_traits::Dupin_cyclide_3 Base_surface_3;

typedef Geo_traits::Surface_3 Surface_3;
typedef Geo_traits::Point_2 Point_3;
typedef Geo_traits::X_monotone_curve_2 X_monotone_arc_3;

typedef Surface_3::Polynomial_3 Polynomial_3;

#define SURFACE_CONSTRUCTION Surface_3::surface_cache()

typedef std::vector< Point_3 > Points_3;
typedef std::vector< X_monotone_arc_3 > XArcs_3;

typedef std::vector< Base_surface_3 > Base_surfaces;
typedef std::vector< Surface_3 > Surface_set;

#endif // XSURFACE_INCLUDES_CKVA_H

