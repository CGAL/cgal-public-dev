// All rights reserved.
// Copyright (c) 2020 SARL GeometryFactory (France).
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

#ifndef CGAL_URBAN_AREA_PROCESSING_INTERNAL_FIND_HOLES_2_H
#define CGAL_URBAN_AREA_PROCESSING_INTERNAL_FIND_HOLES_2_H

// #include <CGAL/license/Urban_area_processing.h>

// STL includes.
#include <vector>
#include <utility>
#include <algorithm>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// Internal includes.
#include <CGAL/Urban_area_processing/internal/Contouring/Boundary_from_triangulation_2.h>

// TODO:
// 1. I should merge this class with some functions from the class Boundary_from_triangulation_2.h.

namespace CGAL {
namespace Urban_area_processing {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap,
  typename InputTriangulation>
  class Find_holes_2 {
  
  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;
    using Triangulation = InputTriangulation;

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
    
    using Face_handle = typename Delaunay::Face_handle;
    using Face_handles = std::vector<Face_handle>;

    Find_holes_2(
      const Input_range& input_range,
      const Point_map point_map,
      Triangulation& triangulation,
      const FT scale,
      const bool verbose) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_triangulation(triangulation),
    m_scale(scale),
    m_verbose(verbose)
    { }

    std::size_t mark_holes() {

      Delaunay delaunay;
      insert_in_triangulation(delaunay);
      Alpha_shape_2 alpha_shape(
        delaunay, m_scale, Alpha_shape_2::GENERAL);
      count_and_mark_faces(alpha_shape);
      const std::size_t num_holes = add_holes(alpha_shape);
      return num_holes;
    }

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;
    Triangulation& m_triangulation;
    const FT m_scale;
    const bool m_verbose;

    void insert_in_triangulation(
      Delaunay& delaunay) const {
      
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

    std::size_t add_holes(
      const Alpha_shape_2& alpha_shape) {

      std::vector<Face_handles> holes;
      Boundary_extractor extractor(alpha_shape, true);
      extractor.get_holes(holes);
      update_triangulation(holes);
      return holes.size();
    }

    void update_triangulation(
      const std::vector<Face_handles>& holes) {

      Point_2 center;
      for (const auto& hole : holes) {
        for (const auto face : hole) {
          internal::compute_barycenter(face, center);
          const auto fh = m_triangulation.locate(center);
          fh->info().label = 0; // outside
          fh->info().tagged = false; // outside
        }
      }
    }
  };

} // internal
} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_INTERNAL_FIND_HOLES_2_H
