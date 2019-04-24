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

    /// Input ground points.
    INPUT_GROUND_POINTS = 0,

    /// Input vegetation points.
    INPUT_VEGETATION_POINTS = 1,

    /// Tree clusters.
    TREE_CLUSTERS = 2,

    /// Tree points.
    TREE_POINTS = 3,

    /// Tree boundaries.
    TREE_BOUNDARIES = 4,

    /// Tree footprints.
    TREE_FOOTPRINTS = 5,

    /// Extruded tree boundaries.
    EXTRUDED_TREE_BOUNDARIES = 6,

    /// Extruded tree footprints.
    EXTRUDED_TREE_FOOTPRINTS = 7,

    /// Tree trunks.
    TREE_TRUNKS = 8,

    /// Tree crowns.
    TREE_CROWNS = 9,

    /// Input building boundary points.
    INPUT_BUILDING_BOUNDARY_POINTS = 10,

    /// Input building interior points.
    INPUT_BUILDING_INTERIOR_POINTS = 11,

    /// Building clusters.
    BUILDING_CLUSTERS = 12,

    /// Building boundary points.
    BUILDING_BOUNDARY_POINTS = 13,

    /// Building wall points.
    BUILDING_WALL_POINTS = 14,

    /// Building approximate boundaries.
    BUILDING_APPROXIMATE_BOUNDARIES = 15,

    /// Building partitioning 2.
    BUILDING_PARTITIONING_2 = 16,

    /// Building points.
    BUILDING_POINTS = 17,

    /// Building boundaries.
    BUILDING_BOUNDARIES = 18,

    /// Building footprints.
    BUILDING_FOOTPRINTS = 19,

    /// Extruded building boundaries.
    EXTRUDED_BUILDING_BOUNDARIES = 20,

    /// Extruded building footprints.
    EXTRUDED_BUILDING_FOOTPRINTS = 21,

    /// Building roof points.
    BUILDING_ROOF_POINTS = 22,

    /// Approximate building bounds.
    APPROXIMATE_BUILDING_BOUNDS = 23,

    /// Building partitioning 3.
    BUILDING_PARTITIONING_3 = 24,

    /// Building walls.
    BUILDING_WALLS = 25,

    /// Building roofs.
    BUILDING_ROOFS = 26

  }; // Intermediate_step

  /// @}

  /// \cond SKIP_IN_MANUAL
  enum class Wire_type {

    /// Planar ground wire.
    PLANAR_GROUND_WIRE = 0,

    /// Smooth ground wire.
    SMOOTH_GROUND_WIRE = 1,

    /// Trees wire 0.
    TREES_WIRE0 = 2,

    /// Trees wire 1.
    TREES_WIRE1 = 3,

    /// Trees wire 2.
    TREES_WIRE2 = 4,

    /// Buildings wire 0.
    BUILDINGS_WIRE0 = 5,

    /// Buildings wire 1.
    BUILDINGS_WIRE1 = 6,

    /// Buildings wire 2.
    BUILDINGS_WIRE2 = 7,

    /// LOD wire 0.
    LOD_WIRE0 = 8,

    /// LOD wire 1.
    LOD_WIRE1 = 9,

    /// LOD wire 2.
    LOD_WIRE2 = 10
  };
  /// \endcond

} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_ENUM_H
