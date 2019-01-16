// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri

#ifndef CGAL_LEVELS_OF_DETAIL_ENUMERATIONS_H
#define CGAL_LEVELS_OF_DETAIL_ENUMERATIONS_H

namespace CGAL {

  namespace Levels_of_detail {

    /*!
      \ingroup PkgLevelsOfDetailRef
      \brief Different enumerations used by the Levels Of Detail algorithm.
    */

    /// \name Semantic Label
    /// @{

    /// This label shows a semantic class a given point belongs to.
    enum class Semantic_label { 
			
      /// Any class that is not handled by the algorithm.
			UNASSIGNED = 0,

      /// Ground points.
			GROUND = 1,

      /// Points treated as a building boundary, e.g. walls.
			BUILDING_BOUNDARY = 2,

      /// Points treated as a building interior, e.g. roofs.
      BUILDING_INTERIOR = 3, 

      /// Vegetation points.
      VEGETATION = 4
		};

    /// @}

  } // Levels_of_detail

} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_ENUMERATIONS_H
