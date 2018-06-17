// Copyright (c) 2018  Liangliang Nan
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Liangliang Nan

#ifndef CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_MIP_SOLVER_CONFIG_H
#define CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_MIP_SOLVER_CONFIG_H


#define HAS_SCIP	0
#define HAS_GLPK	1


#if (HAS_SCIP == 0 && HAS_GLPK == 0)
#error No MIP solver available. 
#endif


#endif	// CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_MIP_SOLVER_CONFIG_H
