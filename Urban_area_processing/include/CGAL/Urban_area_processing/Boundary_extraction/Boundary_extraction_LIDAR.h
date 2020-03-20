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

#ifndef CGAL_URBAN_AREA_PROCESSING_BOUNDARY_EXTRACTION_LIDAR_H
#define CGAL_URBAN_AREA_PROCESSING_BOUNDARY_EXTRACTION_LIDAR_H

// #include <CGAL/license/Urban_area_processing.h>

// CGAL includes.
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Urban_area_processing/internal/Boundary_from_triangulation_2.h>

namespace CGAL {
namespace Urban_area_processing {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Boundary_extraction_LIDAR {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;

    using Fi  = Face_info<Traits>;
    using Fbi = CGAL::Triangulation_face_base_with_info_2<Fi, Traits>;
    using Fb  = CGAL::Alpha_shape_face_base_2<Traits, Fbi>;

    using Vi  = Vertex_info<Traits>;
    using Vbi = CGAL::Triangulation_vertex_base_with_info_2<Vi, Traits>;
    using Vb  = CGAL::Alpha_shape_vertex_base_2<Traits, Vbi>;
    
    using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
    using Delaunay = CGAL::Delaunay_triangulation_2<Traits, Tds>;
    using Alpha_shape_2 = CGAL::Alpha_shape_2<Delaunay>;

    using Boundary_extractor = 
      internal::Boundary_from_triangulation_2<Traits, Alpha_shape_2>;

    struct Triangulation_wrapper {
      Delaunay delaunay;

      bool empty() const {
        return delaunay.number_of_faces() == 0;
      }
    };

    Boundary_extraction_LIDAR(
      const InputRange& input_range,
      const PointMap point_map,
      const FT scale) : 
    m_input_range(input_range),
    m_point_map(point_map),
    m_scale(scale) { 

      CGAL_precondition(scale > FT(0));
      CGAL_precondition(input_range.size() > 0);
      insert_in_triangulation();
    }

    template<typename OutputIterator>
    void extract(OutputIterator boundaries) {

      std::cout << "* extracting LIDAR boundary... " << std::endl;
      CGAL_precondition(!m_triangulation.empty());
      if (m_triangulation.empty()) return;

      Alpha_shape_2 alpha_shape(
        m_triangulation.delaunay, m_scale, Alpha_shape_2::GENERAL);
      count_and_mark_faces(alpha_shape);

      Boundary_extractor extractor(alpha_shape);
      extractor.extract(boundaries);
      std::cout << std::endl;
    }

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const FT m_scale;

    Triangulation_wrapper m_triangulation;

    void insert_in_triangulation() {

      auto& delaunay = m_triangulation.delaunay;
      delaunay.clear();
      for (auto item = m_input_range.begin(); 
      item != m_input_range.end(); ++item)
        delaunay.insert(internal::point_2_from_point_3(
          get(m_point_map, *item)));
    }

    void count_and_mark_faces(
      Alpha_shape_2& alpha_shape) const {

      std::size_t count = 0;
      for (auto face = alpha_shape.finite_faces_begin();
      face != alpha_shape.finite_faces_end(); ++face, ++count) {
        face->info().index = count;
        if (alpha_shape.classify(face) == Alpha_shape_2::INTERIOR)
          face->info().label = 0;
        else 
          face->info().label = std::size_t(-1);
      }
    }
  };

} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_BOUNDARY_EXTRACTION_LIDAR_H
