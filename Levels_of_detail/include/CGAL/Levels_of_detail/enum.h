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
    \ingroup PkgLevelsOfDetailEnumerations

    \brief This enum enables to choose a type of input.
  */
  enum class Input_type {

    /// Input is a set of points with labels.
		POINT_SET = 0,

    /// Input is a soup of polygons with labels.
    POLYGON_SOUP = 1

	}; // Input_type

  /*!
    \ingroup PkgLevelsOfDetailEnumerations

    \brief This label represents a semantic class of an item.
  */
  enum class Semantic_label {

    /// Ground points.
		GROUND = 0,

    /// Items treated as vegetation.
    VEGETATION = 1,

    /// Items treated as a building boundary, e.g. walls.
		BUILDING_BOUNDARY = 2,

    /// Items treated as a building interior, e.g. roofs.
    BUILDING_INTERIOR = 3,

    /// Any item that is not handled by the algorithm.
    UNCLASSIFIED = 4

	}; // Semantic_label

  /*!
    \ingroup PkgLevelsOfDetailEnumerations

    \brief This label represents a position of an item with respect to an object.
  */
	enum class Visibility_label {

    /// Outside the object.
    OUTSIDE = 0,

    /// Inside the object.
    INSIDE = 1,

	}; // Visibility_label

  /*!
    \ingroup PkgLevelsOfDetailEnumerations

    \brief This enum enables to choose a type of extrusion for an object.
  */
  enum class Extrusion_type {

    /// Extrudes the footprint of an object to its average height.
    AVG = 0,

    /// Extrudes the footprint of an object to its maximum height.
    MAX = 1

  }; // Extrusion_type

  /*!
    \ingroup PkgLevelsOfDetailEnumerations

    \brief This enum enables to choose a type of object reconstruction.
  */
  enum class Reconstruction_type {

    /// Ground is represented as a plane.
    PLANAR_GROUND = 0,

    /// Ground is represented as a smooth surface.
    SMOOTH_GROUND = 1,

    /// Trees are represented as discs.
    TREES0 = 2,

    /// Trees are represented as cylinders.
    TREES1 = 3,

    /// Complete trees.
    TREES2 = 4,

    /// Buildings are represented as polygons.
    BUILDINGS0 = 5,

    /// Buildings are represented as extruded polygons.
    BUILDINGS1 = 6,

    /// Complete buildings.
    BUILDINGS2 = 7,

    /// All available objects with the level of detail 0.
    LOD0 = 8,

    /// All available objects with the level of detail 1.
    LOD1 = 9,

    /// All available objects with the level of detail 2.
    LOD2 = 10

	}; // Reconstruction_type

  /*!
    \ingroup PkgLevelsOfDetailEnumerations

    \brief This enum represents different types of objects.
  */
  enum class Urban_object_type {

    /// Ground.
    GROUND = 0,

    /// Tree trunk.
    TREE_TRUNK = 1,

    /// Tree crown.
    TREE_CROWN = 2,

    /// Building wall.
    BUILDING_WALL = 3,

    /// Building roof.
    BUILDING_ROOF = 4

  }; // Urban_object_type

  /*!
    \ingroup PkgLevelsOfDetailEnumerations

    \brief This enum enables to choose an intermediate step of the algorithm to
    retreive the corresponding data.
  */
  enum class Intermediate_step {

    /// Input ground points.
    INPUT_GROUND_POINTS = 0,

    /// Input vegetation points.
    INPUT_VEGETATION_POINTS = 1,

    /// All vegetation points are returned with the index of the corresponding cluster.
    TREE_CLUSTERS = 2,

    /// All vegetation points are returned with the index of the corresponding tree.
    TREE_POINTS = 3,

    /// Tree boundaries are returned as a set of polylines with the index of the corresponding tree.
    TREE_BOUNDARIES = 4,

    /// Tree footprints are returned as a polygon soup with the index of the corresponding tree.
    TREE_FOOTPRINTS = 5,

    /// Extruded tree boundaries are returned as a polygon soup with the index of the corresponding tree.
    EXTRUDED_TREE_BOUNDARIES = 6,

    /// Extruded tree footprints are returned as a polygon soup with the index of the corresponding tree.
    EXTRUDED_TREE_FOOTPRINTS = 7,

    /// Tree trunks are returned as a polygon soup with the index of the corresponding tree.
    TREE_TRUNKS = 8,

    /// Tree crowns are returned as a polygon soup with the index of the corresponding tree.
    TREE_CROWNS = 9,

    /// Input building boundary points.
    INPUT_BUILDING_BOUNDARY_POINTS = 10,

    /// Input building interior points.
    INPUT_BUILDING_INTERIOR_POINTS = 11,

    /// All building points are returned with the index of the corresponding cluster.
    BUILDING_CLUSTERS = 12,

    /// Building boundary points.
    BUILDING_BOUNDARY_POINTS = 13,

    /// Building wall points are returned with the index of the corresponding wall.
    BUILDING_WALL_POINTS = 14,

    /// Building approximate boundaries are returned as a set of polylines.
    BUILDING_APPROXIMATE_BOUNDARIES = 15,

    /// Building partitioning is returned as a planar polygon soup with faces labeled inside/outside.
    BUILDING_PARTITIONING_2 = 16,

    /// Building points are returned with the index of the corresponding building.
    BUILDING_POINTS = 17,

    /// Building boundaries are returned as a set of polylines with the index of the corresponding building.
    BUILDING_BOUNDARIES = 18,

    /// Building footprints are returned as a polygon soup with the index of the corresponding building.
    BUILDING_FOOTPRINTS = 19,

    /// Extruded building boundaries are returned as a polygon soup with the index of the corresponding building.
    EXTRUDED_BUILDING_BOUNDARIES = 20,

    /// Extruded building footprints are returned as a polygon soup with the index of the corresponding building.
    EXTRUDED_BUILDING_FOOTPRINTS = 21,

    /// Building roof points are returned with the index of the corresponding roof.
    BUILDING_ROOF_POINTS = 22,

    /// Approximate building bounds are returned as a quadrilateral soup.
    APPROXIMATE_BUILDING_BOUNDS = 23,

    /// Building partitioning is returned as a polygon soup that bounds the corresponding building's interior.
    BUILDING_PARTITIONING_3 = 24,

    /// Building walls are returned as a polygon soup with the index of the corresponding building.
    BUILDING_WALLS = 25,

    /// Building roofs are returned as a polygon soup with the index of the corresponding building.
    BUILDING_ROOFS = 26

  }; // Intermediate_step

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
