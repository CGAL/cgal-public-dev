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
// Author(s)     : Liangliang Nan

#ifndef CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_MIP_SOLVER_H
#define CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_MIP_SOLVER_H

#include <CGAL/mip/mip_model.h>

/*!
\file mip_solver.h
*/

namespace CGAL {

	/** \ingroup PkgPolygonalSurfaceReconstruction
	*
	*	A Mixed Integer Program solver (It can also be used to solve general
	*   linear programs.
	*   Currently it is a wrapper encapsulating SCIP.
	*/
	class MIP_solver
	{
	public:
		MIP_solver() {}
		~MIP_solver() {}

		/// Solves the problem and returns false if fails.
		virtual bool solve(const MIP_model* model);

		/// Returns the result. 
		/// The result can also be retrieved using Variable::solution_value().
		/// NOTE: (1) result is valid only if the solver succeeded.
		///       (2) each entry in the result corresponds to the variable with the
		///			 same index in the model.
		const std::vector<double>& solution() const { return result_; }

	private:
		std::vector<double>	result_;
	};


} //namespace CGAL


#include <CGAL/mip/mip_solver_interface_SCIP.h>


#endif	// CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_MIP_SOLVER_H
