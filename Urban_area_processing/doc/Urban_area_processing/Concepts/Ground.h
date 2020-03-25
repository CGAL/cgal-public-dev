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

#ifndef CGAL_URBAN_AREA_PROCESSING_GROUND_CONCEPT_H
#define CGAL_URBAN_AREA_PROCESSING_GROUND_CONCEPT_H

// #include <CGAL/license/Urban_area_processing.h>

/*!
\ingroup PkgUrbanAreaProcessingRefConcepts
\cgalConcept

A concept that describes the set of methods required by all classes that represent 
a ground object. The ground object is an important building block of the urban reconstruction
pipeline since it represents a base where all available urban objects are placed.

\cgalHasModel 
- `CGAL::Urban_area_processing::Planar_ground`,
- `CGAL::Urban_area_processing::Smooth_ground`
*/
class Ground {

public:
  
  /*!  
    reconstructs the ground.
  */
  void reconstruct() {

  }

  /*!  
    places an urban object onto the ground.

    \tparam UrbanObject 
    must be a model of `Urban_object` that is a tree, a building, etc.
  */
  template<typename UrbanObject>
  void place(const UrbanObject& urban_object) {

  }

  /*!  
    updates all internal data structures, if necessary, after placing different 
    urban objects.
  */
  void update() {

  }

  /*!  
    outputs the ground object to `os` in City GML format.
  */
  void output_to_city_gml(std::ofstream& os) {

  }
};

#endif // CGAL_URBAN_AREA_PROCESSING_GROUND_CONCEPT_H
