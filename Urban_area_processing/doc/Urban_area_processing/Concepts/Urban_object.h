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

#ifndef CGAL_URBAN_AREA_PROCESSING_URBAN_OBJECT_CONCEPT_H
#define CGAL_URBAN_AREA_PROCESSING_URBAN_OBJECT_CONCEPT_H

// #include <CGAL/license/Urban_area_processing.h>

/*!
\ingroup PkgUrbanAreaProcessingRefConcepts
\cgalConcept

A concept that describes the set of methods required by all classes that represent 
different urban objects. Such objects may be trees, buildings, etc.

\cgalHasModel 
- `CGAL::Urban_area_processing::Tree`,
- `CGAL::Urban_area_processing::Building`
*/
class Urban_object {

public:

  /*!  
    returns a polygon with holes, which represents the footprint of the urban object.

    \tparam OutputIterator 
    must be an output iterator whose value type is `std::vector< std::pair<Point_3, std::size_t> >`,
    where the first item in the pair is a point and second item is the contour index. 
    If the latter is `std::size_t(-1)` then this contour is outer, otherwise it is a hole
    and the stored index is the index of the corresponding outer contour.
  */
  template<typename OutputIterator>
  void footprint(OutputIterator contours) {

  }

  /*!  
    returns a set of triangles, which represent the building blocks of the object, 
    which constitute to its middle part such as walls in a building and trunk in a tree.

    \tparam OutputIterator 
    must be an output iterator whose value type is `Triangle_3`.
  */
  template<typename OutputIterator>
  void boundary(OutputIterator triangles) {

  }

  /*!  
    returns a set of triangles, which represent the building blocks of the object, 
    which constitute to its upper part such as roofs in a building and crown in a tree.

    \tparam OutputIterator 
    must be an output iterator whose value type is `Triangle_3`.
  */
  template<typename OutputIterator>
  void top(OutputIterator triangles) {

  }

  /*!  
    outputs the urban object to os in City GML format.
  */
  void output_to_city_gml(std::ofstream& os) {

  }
};

#endif // CGAL_URBAN_AREA_PROCESSING_URBAN_OBJECT_CONCEPT_H
