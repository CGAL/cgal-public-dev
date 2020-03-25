// Copyright (c) 2020 SARL GeometryFactory (France).
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

#ifndef CGAL_URBAN_AREA_PROCESSING_ENUM_H
#define CGAL_URBAN_AREA_PROCESSING_ENUM_H

// #include <CGAL/license/Urban_area_processing.h>

namespace CGAL {
namespace Urban_area_processing {

  /*!
    \ingroup PkgUrbanAreaProcessingRefEnum
      
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
    \ingroup PkgUrbanAreaProcessingRefEnum
      
    \brief This label represents a position of an item with respect to an object.
  */
  enum class Visibility_label {

    /// Outside the object.
    OUTSIDE = 0,

    /// Inside the object.
    INSIDE = 1,

  }; // Visibility_label

  /*!
    \ingroup PkgUrbanAreaProcessingRefEnum
      
    \brief This enum enables to choose a type of extrusion for an object.
  */
  enum class Extrusion_type { 

    /// Extrudes the footprint of an object to its average height.
    AVG = 0,

    /// Extrudes the footprint of an object to its maximum height.
    MAX = 1

  }; // Extrusion_type

  /*!
    \ingroup PkgUrbanAreaProcessingRefEnum
      
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
    \ingroup PkgUrbanAreaProcessingRefEnum
      
    \brief This enum enables to choose a type of object reconstruction.
  */ 
  enum class Reconstruction_type { 
      
    /// Level of detail 0.
    LOD0 = 0,

    /// Level of detail 1.
    LOD1 = 1,

    /// Level of detail 2.
    LOD2 = 2

  }; // Reconstruction_type

} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_ENUM_H
