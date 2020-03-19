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
#include <CGAL/Shape_detection/Region_growing.h>

// Internal includes.
#include <CGAL/Urban_area_processing/struct.h>

#include "../../../../test/Urban_area_processing/include/Saver.h"

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

    using Face_handle = typename Delaunay::Face_handle;
    using Vertex_handle = typename Delaunay::Vertex_handle;
    
    using Neighbor_query = 
      internal::Triangulation_neighbor_query_2<Traits, Face_handle>;
    using Connected_region = 
      internal::Triangulation_connected_region_2<Traits, Face_handle>;
    using Region_growing = 
      CGAL::Shape_detection::Region_growing<
      std::vector<Face_handle>, Neighbor_query, Connected_region>;

    using Indices = std::vector<std::size_t>;

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

      std::cout << "LIDAR boundary extraction: " << std::endl;
      CGAL_precondition(!m_triangulation.empty());
      if (m_triangulation.empty()) return;

      Alpha_shape_2 alpha_shape(
        m_triangulation.delaunay, m_scale, Alpha_shape_2::GENERAL);
      traverse_boundary(alpha_shape, boundaries);
    }

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const FT m_scale;

    Triangulation_wrapper m_triangulation;

    void insert_in_triangulation() {

      auto& delaunay = m_triangulation.delaunay;
      delaunay.clear();
      for (auto it = m_input_range.begin(); it != m_input_range.end(); ++it)
        delaunay.insert(internal::point_2_from_point_3(
          get(m_point_map, *it)));
    }

    template<typename OutputIterator>
    void traverse_boundary(
      Alpha_shape_2& alpha_shape,
      OutputIterator boundaries) const {

      count_and_mark_faces(alpha_shape);
      const std::vector<Face_handle> faces(
        alpha_shape.finite_face_handles().begin(), 
        alpha_shape.finite_face_handles().end());

      std::vector<Indices> regions;
      find_connected_regions(1, faces, regions);
      std::cout << "* regions detected: " << regions.size() << std::endl;

      const std::size_t num_outer_faces = mark_outer_faces(alpha_shape);
      std::cout << "* num outer faces: " << num_outer_faces << std::endl;

      const std::size_t num_holes = mark_holes(faces);
      std::cout << "* num holes: " << num_holes << std::endl;
      
      Saver<Traits> saver;
      for (std::size_t i = 0; i < regions.size(); ++i) {
        const auto& region = regions[i];

        saver.export_polygon_soup(faces, region, 
        "/Users/monet/Documents/gf/urban-area-processing/logs/region_" + std::to_string(i));
        extract_boundaries(
          alpha_shape, faces, region, boundaries);
      }
    }

    void count_and_mark_faces(
      Alpha_shape_2& alpha_shape) const {

      std::size_t count = 0;
      for (auto fit = alpha_shape.finite_faces_begin();
      fit != alpha_shape.finite_faces_end(); ++fit, ++count) {
        fit->info().object_index = count;
        if (alpha_shape.classify(fit) == Alpha_shape_2::INTERIOR)
          fit->info().label = 1;
        else 
          fit->info().label = std::size_t(-1);
      }
    }

    void find_connected_regions(
      const std::size_t ref_label,
      const std::vector<Face_handle>& faces,
      std::vector<Indices>& regions) const {

      Neighbor_query neighbor_query(ref_label, faces);
      Connected_region connected_region(ref_label, faces);
      Region_growing region_growing(
        faces, neighbor_query, connected_region);

      regions.clear();
      region_growing.detect(std::back_inserter(regions));

      if (ref_label == 0) {
        for (std::size_t i = 0; i < regions.size(); ++i) {
          for (const std::size_t idx : regions[i]) {
            const auto fh = *(faces.begin() + idx);
            fh->info().tagged = true;
          }
        }
      }
    }

    std::size_t mark_outer_faces(
      const Alpha_shape_2& alpha_shape) const {

      std::size_t num_outer_faces = 0;
      set_default_state(alpha_shape);
      for (auto fit = alpha_shape.all_faces_begin(); 
      fit != alpha_shape.all_faces_end(); ++fit) {
        if (alpha_shape.is_infinite(fit)) {
          fit->info().label = 0;
          propagate(alpha_shape, fit, num_outer_faces);
        }
      }
      set_default_state(alpha_shape);
      return num_outer_faces;
    }

    void set_default_state(
      const Alpha_shape_2& alpha_shape) const {
      for (auto fit = alpha_shape.all_faces_begin(); 
      fit != alpha_shape.all_faces_end(); ++fit)
        fit->info().used = false;
    }

    void propagate(
      const Alpha_shape_2& alpha_shape,
      const Face_handle fh,
      std::size_t& num_outer_faces) const {

      bool stop = false;
      do {
        stop = true;
        for (std::size_t k = 0; k < 3; ++k) {
          const auto fhn = fh->neighbor(k);

          if (alpha_shape.is_infinite(fhn)) continue;
          if (fhn->info().used) continue;
          if (fhn->info().label == 1) continue;

          num_outer_faces += 1;
          fhn->info().used = true;
          fhn->info().label = 0;
          propagate(alpha_shape, fhn, num_outer_faces);
          stop = false;
        }
      } while (!stop);
    }

    std::size_t mark_holes(
      const std::vector<Face_handle>& faces) const {

      std::vector<Indices> regions;
      find_connected_regions(std::size_t(-1), faces, regions);

      for (std::size_t i = 0; i < regions.size(); ++i) {
        for (const std::size_t idx : regions[i]) {
          const auto fh = *(faces.begin() + idx);
          fh->info().label = i + 2;
        }
      }
      return regions.size();
    }

    template<typename OutputIterator>
    void extract_boundaries(
      const Alpha_shape_2& alpha_shape,
      const std::vector<Face_handle>& faces,
      const Indices& region,
      OutputIterator boundaries) const {

      std::set<std::size_t> labels;
      for (auto fit = alpha_shape.all_faces_begin(); 
      fit != alpha_shape.all_faces_end(); ++fit) {
        const std::size_t label = fit->info().label;
        if (label != 1) labels.insert(label);
      }

      std::vector<Point_2> contour;
      for (const std::size_t label : labels) {
        const auto edge = get_boundary_edge(label, faces, region);
        if (edge.second != std::size_t(-1)) {

          traverse(label, alpha_shape, faces, region, edge, contour);
          std::cout << "- found contour (size): " << contour.size() << std::endl;
          *(++boundaries) = std::make_pair(contour, std::size_t(-1));
        }
      }
    }

    std::pair<Face_handle, std::size_t> get_boundary_edge(
      const std::size_t ref_label,
      const std::vector<Face_handle>& faces,
      const Indices& region) const {

      for (const std::size_t idx : region) {
        const auto fh = *(faces.begin() + idx);
        if (fh->info().label != 1) continue;

        for (std::size_t k = 0; k < 3; ++k) {
          const auto fhn = fh->neighbor(k);
          if (fhn->info().label == ref_label)
            return std::make_pair(fh, k);
        }
      }
      return std::make_pair(Face_handle(), std::size_t(-1));
    }

    void traverse(
      const std::size_t ref_label,
      const Alpha_shape_2& alpha_shape,
      const std::vector<Face_handle>& faces,
      const Indices& region,
      const std::pair<Face_handle, std::size_t>& edge,
      std::vector<Point_2>& contour) const {

      contour.clear();
      if (region.size() == 1) {
        contour.push_back(edge.first->vertex(0)->point());
        contour.push_back(edge.first->vertex(1)->point());
        contour.push_back(edge.first->vertex(2)->point());
      }

      const auto start = edge;
      auto next = edge;
      CGAL_precondition(start.first != Face_handle());
      
      do {
        auto fh = next.first;
        std::size_t k = next.second;
        CGAL_precondition(k != std::size_t(-1));

        const std::size_t kadd = (k + 1) % 3;
        const auto& point = fh->vertex(kadd)->point();
        contour.push_back(point);

        std::size_t knext = (k + 2) % 3;
        auto vh = fh->vertex(knext);
        next = get_next_edge(alpha_shape, vh, fh, ref_label);
        
        if (next.second == std::size_t(-1)) {
          contour.push_back(vh->point());
          knext = (knext + 1) % 3;
          vh = fh->vertex(knext);
          next = get_next_edge(alpha_shape, vh, fh, ref_label);

          if (next.second == std::size_t(-1)) {
            contour.push_back(vh->point());
            knext = (knext + 1) % 3;
            vh = fh->vertex(knext);
            next = get_next_edge(alpha_shape, vh, fh, ref_label);
          }
        }
      } while (next != start);
    }

    std::pair<Face_handle, std::size_t> get_next_edge(
      const Alpha_shape_2& alpha_shape,
      const Vertex_handle vh,
      const Face_handle ref_fh,
      const std::size_t ref_label) const {

      auto fh = alpha_shape.incident_faces(vh);
      const auto start = fh;
      do {
        if (fh->info().label == 1 && fh != ref_fh) {
          for (std::size_t k = 0; k < 3; ++k) {
            const auto fhn = fh->neighbor(k);
            if (fhn->info().label == ref_label)
              return std::make_pair(fh, k);
          }
        }
        ++fh;
      } while (fh != start);
      return std::make_pair(Face_handle(), std::size_t(-1));
    }
  };

} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_BOUNDARY_EXTRACTION_LIDAR_H
