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

  /// This label represents a semantic class a given point belongs to.
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

  /// \name Visibility Label
  /// @{

	enum class Visibility_label {

    // All elements that do not belong to any object.
    OUTSIDE = 0,

    // All elements that belong to a given object.
    INSIDE = 1
	};

  /// @}

  /// \name Extrusion Type
  /// @{

  /// This enumeration allows to choose a type of extrusion.
  enum class Extrusion_type { 
  
    /// Extrudes the top of the object to its minimum height.
    MIN = 0,

    /// Extrudes the top of the object to its average height.
    AVERAGE = 1,

    /// Extrudes the top of the object to its maximum height.
    MAX = 2
  };  

  /// @}

  /// \name Reconstruction Type
  /// @{

  /// This enumeration allows to choose a type of reconstruction.
  enum class Reconstruction_type { 
			
    /// Only ground represented as a plane.
    PLANAR_GROUND = 0,

    /// Only ground represented as a smooth surface.
    SMOOTH_GROUND = 1,

    /// Only buildings as footprints.
    BUILDINGS0 = 2,

    /// Only buildings as boxes.
    BUILDINGS1 = 3,

    /// Only buildings.
    BUILDINGS2 = 4,

    /// Only trees as footprints.
    VEGETATION0 = 5,

    /// Only trees as cylinders.
    VEGETATION1 = 6,

    /// Only trees.
    VEGETATION2 = 7,

    /// All objects with level of detail 0.
    LOD0 = 8,

    /// All objects with level of detail 1.
    LOD1 = 9,

    /// All objects with level of detail 2.
    LOD2 = 10
	};

  /// This enumeration represents different type of urban objects.
  enum class Urban_object_type {

    /// Ground.
    GROUND = 0,

    /// Building.
    BUILDING = 1,

    /// Tree.
    TREE = 2
  };

  /// @}

} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_ENUMERATIONS_H
