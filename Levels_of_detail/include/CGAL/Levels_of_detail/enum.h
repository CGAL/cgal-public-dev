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
//
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri
//

#ifndef CGAL_LEVELS_OF_DETAIL_ENUM_H
#define CGAL_LEVELS_OF_DETAIL_ENUM_H

#include <CGAL/license/Levels_of_detail.h>

namespace CGAL {
namespace Levels_of_detail {

  /*!
    \ingroup PkgLevelsOfDetailRef
      
    \brief Various enums used by the `CGAL::Levels_of_detail::Levels_of_detail`.
  */

  /// \name Semantic Label
  /// @{

  /// This label represents a semantic class of a point.
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
    VEGETATION = 4,

    // Unspecified.
    UNSPECIFIED = 5

	}; // Semantic_label

  /// @}

  /// \name Visibility Label
  /// @{

  /// This label represents a position of an item with respect to an object.
	enum class Visibility_label {

    // Outside the object.
    OUTSIDE = 0,

    // Inside the object.
    INSIDE = 1,

    // Unspecified.
    UNSPECIFIED = 2

	}; // Visibility_label

  /// @}

  /// \name Extrusion Type
  /// @{

  /// This enum enables to choose a type of extrusion for an object.
  enum class Extrusion_type { 

    /// Extrudes the footprint of the object to its average height.
    AVG = 0,

    /// Extrudes the footprint of the object to its maximum height.
    MAX = 1,

    // Unspecified.
    UNSPECIFIED = 2

  }; // Extrusion_type

  /// @}

  /// \name Reconstruction Type
  /// @{

  /// This enum enables to choose a type of reconstruction.
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
    TREES0 = 5,

    /// Only trees as cylinders.
    TREES1 = 6,

    /// Only trees.
    TREES2 = 7,

    /// All objects with the level of detail 0.
    LOD0 = 8,

    /// All objects with the level of detail 1.
    LOD1 = 9,

    /// All objects with the level of detail 2.
    LOD2 = 10,

    // Unspecified.
    UNSPECIFIED = 11

	}; // Reconstruction_type

  /// \name Urban Object Type
  /// @{

  /// This enum represents different types of urban objects.
  enum class Urban_object_type {

    /// Ground.
    GROUND = 0,

    /// Building.
    BUILDING = 1,

    /// Tree.
    TREE = 2,

    /// Unspecified.
    UNSPECIFIED = 3

  }; // Urban_object_type

  /// @}

  /// \name Intermediate Step
  /// @{

  /// This enum enables to choose different intermediate steps of the algorithm.
  enum class Intermediate_step {

    /// Tree clusters.
    TREE_CLUSTERS = 0,

    /// Tree points.
    TREE_POINTS = 1,

    /// Tree boundaries.
    TREE_BOUNDARIES = 2,

    /// Tree footprints.
    TREE_FOOTPRINTS = 3,

    /// Extruded tree boundaries.
    EXTRUDED_TREE_BOUNDARIES = 4,

    /// Extruded tree footprints.
    EXTRUDED_TREE_FOOTPRINTS = 5,

    /// Tree trunks.
    TREE_TRUNKS = 6,

    /// Tree crowns.
    TREE_CROWNS = 7,

    /// Building clusters.
    BUILDING_CLUSTERS = 8,

    /// Building boundary points.
    BUILDING_BOUNDARY_POINTS = 9,

    /// Building wall points.
    BUILDING_WALL_POINTS = 10,

    /// Building approximate boundaries.
    BUILDING_APPROXIMATE_BOUNDARIES = 11,

    /// Building partitioning 2.
    BUILDING_PARTITIONING_2 = 12,

    /// Building points.
    BUILDING_POINTS = 13,

    /// Building boundaries.
    BUILDING_BOUNDARIES = 14,

    /// Building footprints.
    BUILDING_FOOTPRINTS = 15,

    /// Extruded building boundaries.
    EXTRUDED_BUILDING_BOUNDARIES = 16,

    /// Extruded building footprints.
    EXTRUDED_BUILDING_FOOTPRINTS = 17,

    /// Building roof points.
    BUILDING_ROOF_POINTS = 18,

    /// Approximate building bounds.
    APPROXIMATE_BUILDING_BOUNDS = 19,

    /// Building partitioning 3.
    BUILDING_PARTITIONING_3 = 20,

    /// Building walls.
    BUILDING_WALLS = 21,

    /// Building roofs.
    BUILDING_ROOFS = 22

  }; // Intermediate_step

  /// @}

} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_ENUM_H