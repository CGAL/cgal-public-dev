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

#ifndef CGAL_URBAN_AREA_PROCESSING_INTERNAL_BOUNDARY_FROM_TRIANGULATION_H
#define CGAL_URBAN_AREA_PROCESSING_INTERNAL_BOUNDARY_FROM_TRIANGULATION_H

// #include <CGAL/license/Urban_area_processing.h>

// STL includes.
#include <vector>
#include <utility>
#include <algorithm>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Shape_detection/Region_growing.h>

// Internal includes.
#include <CGAL/Urban_area_processing/struct.h>
#include <CGAL/Urban_area_processing/internal/utils.h>
#include <CGAL/Urban_area_processing/internal/Shape_detection/Triangulation_labeled_region_2.h>
#include <CGAL/Urban_area_processing/internal/Shape_detection/Triangulation_neighbor_query_2.h>

// Utils.
#include "../../../../../test/Urban_area_processing/include/Saver.h"

// TODO:
// 1. Make this algorithm more reliable. Now, it fails on data sets with multiple
// little holes and subregions, which are connected by only a vertex rather than edge.
// The problem is how to identify in this case what is the main region and what is 
// its hole. Another issue is when a part of one region is only one triangle away
// from the other region. Fix could be either create an alternative region growing
// that includes faces connected by a vertex at the same region or create a separate thread 
// for each region and work locally around it that is creating its holes and traversing them.
// 2. Should I change indices from 0 - inside, 1 - outside to 1 - inside, 0 - outside?

namespace CGAL {
namespace Urban_area_processing {
namespace internal {

  template<
  typename GeomTraits,
  typename InputTriangulation>
  class Boundary_from_triangulation_2 {

  public:
    using Traits = GeomTraits;
    using Triangulation = InputTriangulation;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;

    using Face_handle = typename Triangulation::Face_handle;
    using Vertex_handle = typename Triangulation::Vertex_handle;
    
    using Neighbor_query = 
      internal::Triangulation_neighbor_query_2<Traits, Face_handle>;
    using Labeled_region = 
      internal::Triangulation_labeled_region_2<Traits, Face_handle>;
    using Region_growing = 
      CGAL::Shape_detection::Region_growing<
      std::vector<Face_handle>, Neighbor_query, Labeled_region>;

    using Indices = std::vector<std::size_t>;
    using Points_2 = std::vector<Point_2>;
    using Face_handles = std::vector<Face_handle>;

    Boundary_from_triangulation_2(
      const Triangulation& triangulation,
      const bool verbose = true) : 
    m_triangulation(triangulation),
    m_faces(
      triangulation.finite_face_handles().begin(), 
      triangulation.finite_face_handles().end()),
      m_verbose(verbose) { 
      
      CGAL_precondition(triangulation.number_of_faces() > 0);
      CGAL_precondition(m_faces.size() == triangulation.number_of_faces());
    }

    void get_holes(
      std::vector<Face_handles>& holes) {

      const std::size_t ref_label = 0; // all interior regions are labeled 0
      const std::size_t min_faces = 2; // min number of faces in the region
      
      std::vector<Indices> regions;
      find_connected_regions(ref_label, min_faces, regions);
      if (m_verbose)
        std::cout << "- interior regions are detected: " << regions.size() << std::endl;

      tag_interior_faces(regions);
      if (m_verbose)
        std::cout << "- interior faces are tagged" << std::endl;

      const std::size_t num_outer_faces = mark_outer_faces();
      if (m_verbose)
        std::cout << "- outer faces are marked: " << num_outer_faces << std::endl;

      regions.clear();
      mark_holes(regions);
      if (m_verbose)
        std::cout << "- holes are found: " << regions.size() << std::endl;

      extract_holes(regions, holes);
      if (m_verbose)
        std::cout << "- holes are extracted: " << holes.size() << std::endl;
    }

    template<typename OutputIterator>
    void extract(OutputIterator boundaries) {

      const std::size_t ref_label = 0; // all interior regions are labeled 0
      const std::size_t min_faces = 2; // min number of faces in the region
      
      std::vector<Indices> regions;
      find_connected_regions(ref_label, min_faces, regions);
      if (m_verbose)
        std::cout << "- interior regions are detected: " << regions.size() << std::endl;

      tag_interior_faces(regions);
      if (m_verbose)
        std::cout << "- interior faces are tagged" << std::endl;

      save_regions(regions);
      if (m_verbose)
        std::cout << "- regions are saved" << std::endl;

      const std::size_t num_outer_faces = mark_outer_faces();
      if (m_verbose)
        std::cout << "- outer faces are marked: " << num_outer_faces << std::endl;
      
      /*
      save_triangulation(1, "outer");
      if (m_verbose)
        std::cout << "- outer triangulation is saved" << std::endl; */

      std::vector<Indices> holes;
      mark_holes(holes);
      const std::size_t num_holes = holes.size();
      if (m_verbose)
        std::cout << "- holes are found: " << num_holes << std::endl;

      /*
      for (std::size_t i = 2; i < num_holes + 2; ++i)
        save_triangulation(i, "holes/hole_" + std::to_string(i));
      if (m_verbose)
        std::cout << "- holes are saved" << std::endl; */

      Indices labels;
      const std::size_t num_labels = num_holes + 2;
      create_reference_labels(num_labels, labels);
      if (m_verbose)
        std::cout << "- reference labels are created: " << labels.size() << std::endl;

      std::size_t num_boundaries = 0;
      for (std::size_t i = 0; i < regions.size(); ++i) {
        const auto& region = regions[i];
        num_boundaries += extract_boundaries(region, i, labels, boundaries);
      }
      if (m_verbose)
        std::cout << "- boundaries are extracted: " << num_boundaries << std::endl;
    }

  private:
    const Triangulation& m_triangulation;
    const std::vector<Face_handle> m_faces;
    const bool m_verbose;

    void find_connected_regions(
      const std::size_t ref_label,
      const std::size_t min_faces,
      std::vector<Indices>& regions) const {

      Neighbor_query neighbor_query(m_faces);
      Labeled_region labeled_region(m_faces, 
        ref_label, min_faces);
      Region_growing region_growing(m_faces, 
        neighbor_query, labeled_region);

      regions.clear();
      region_growing.detect(std::back_inserter(regions));
    }

    void tag_interior_faces(
      const std::vector<Indices>& regions) const {

      for (const auto face : m_faces) {
        face->info().tagged = false;
        face->info().label = std::size_t(-1);
      }

      for (std::size_t i = 0; i < regions.size(); ++i) {
        const auto& region = regions[i];
        for (const std::size_t face_index : region) {
          const auto face = *(m_faces.begin() + face_index);
          face->info().tagged = true;
          face->info().object_index = i;
          face->info().label = 0;
        }
      }
    }

    void save_regions(
      const std::vector<Indices>& regions) const {

      Saver<Traits> saver;
      for (std::size_t i = 0; i < regions.size(); ++i) {
        const auto& region = regions[i];
        saver.export_polygon_soup(m_faces, region, 
        "/Users/monet/Documents/gf/urban-area-processing/logs/region_" + 
        std::to_string(i));
      }
    }

    std::size_t mark_outer_faces() {

      std::size_t num_outer_faces = 0;
      set_default_state();
      for (auto face = m_triangulation.all_faces_begin(); 
      face != m_triangulation.all_faces_end(); ++face) {
        if (m_triangulation.is_infinite(face)) {
          face->info().label = 1; // we use 1 to mark all exterior faces
          propagate(face, num_outer_faces);
        }
      }
      set_default_state();
      return num_outer_faces;
    }

    void set_default_state() {
      for (auto face = m_triangulation.all_faces_begin(); 
      face != m_triangulation.all_faces_end(); ++face)
        face->info().used = false;
    }

    void propagate(
      const Face_handle face,
      std::size_t& num_outer_faces) const {

      bool stop = false;
      do {
        stop = true;
        for (std::size_t k = 0; k < 3; ++k) {
          const auto nface = face->neighbor(k);

          if (m_triangulation.is_infinite(nface)) continue; // infinite
          if (nface->info().used) continue; // already traversed
          if (nface->info().label == 0) continue; // interior faces

          num_outer_faces += 1;
          nface->info().used = true;
          nface->info().label = 1;
          propagate(nface, num_outer_faces);
          stop = false;
        }
      } while (!stop);
    }

    void save_triangulation(
      const std::size_t ref_label,
      const std::string name) const {

      Saver<Traits> saver;
      saver.export_polygon_soup(m_triangulation, ref_label, 
      "/Users/monet/Documents/gf/urban-area-processing/logs/" + name);
    }

    void mark_holes(
      std::vector<Indices>& regions) const {

      const std::size_t ref_label = std::size_t(-1); // hole faces
      const std::size_t min_faces = 1; // min number of faces in the region

      regions.clear();
      find_connected_regions(ref_label, min_faces, regions);

      for (std::size_t i = 0; i < regions.size(); ++i) {
        const auto& region = regions[i];
        for (const std::size_t face_index : region) {
          const auto face = *(m_faces.begin() + face_index);
          face->info().label = i + 2;
        }
      }
    }

    void create_reference_labels(
      const std::size_t num_labels,
      Indices& labels) const {

      labels.clear();
      labels.reserve(num_labels - 1); // zero label is interior, skip it
      for (std::size_t i = 1; i < num_labels; ++i)
        labels.push_back(i);
    }

    template<typename OutputIterator>
    std::size_t extract_boundaries(
      const Indices& region,
      const std::size_t region_index,
      const Indices& labels,
      OutputIterator boundaries) const {

      Points_2 contour;
      std::vector<Points_2> contours;

      // Compute contours.
      for (const std::size_t label : labels) {
        const auto edge = get_boundary_edge(label, region);
        if (edge.second != std::size_t(-1)) {
          traverse(label, region, region_index, edge, contour);
          if (contour.size() >= 3)
            contours.push_back(contour);
        }
      }

      // Set main contour and its holes.
      std::sort(contours.begin(), contours.end(), 
      [](const Points_2& a, const Points_2& b) -> bool { 

        FT a_perimeter = FT(0);
        for (std::size_t i = 0; i < a.size(); ++i) {
          const std::size_t ip = (i + 1) % a.size();
          a_perimeter += CGAL::squared_distance(a[i], a[ip]);
        }

        FT b_perimeter = FT(0);
        for (std::size_t i = 0; i < b.size(); ++i) {
          const std::size_t ip = (i + 1) % b.size();
          b_perimeter += CGAL::squared_distance(b[i], b[ip]);
        }

        return a_perimeter > b_perimeter;
      });

      if (contours.size() == 0) return 0;
      *(++boundaries) = std::make_pair(contours[0], std::size_t(-1));
      for (std::size_t i = 1; i < contours.size(); ++i)
        *(++boundaries) = std::make_pair(contours[i], region_index);
      return contours.size();
    }

    std::pair<Face_handle, std::size_t> get_boundary_edge(
      const std::size_t ref_label,
      const Indices& region) const {

      for (const std::size_t face_index : region) {
        const auto face = *(m_faces.begin() + face_index);
        CGAL_assertion(face->info().label == 0); // interior face

        for (std::size_t k = 0; k < 3; ++k) {
          const auto nface = face->neighbor(k);
          if (nface->info().label == ref_label)
            return std::make_pair(face, k);
        }
      }
      return std::make_pair(Face_handle(), std::size_t(-1));
    }

    void traverse(
      const std::size_t ref_label,
      const Indices& region,
      const std::size_t region_index,
      const std::pair<Face_handle, std::size_t>& edge,
      Points_2& contour) const {

      CGAL_assertion(
        edge.second != std::size_t(-1));

      contour.clear();
      if (region.size() == 1) {
        contour.push_back(edge.first->vertex(0)->point());
        contour.push_back(edge.first->vertex(1)->point());
        contour.push_back(edge.first->vertex(2)->point());
      }
      
      std::size_t count = 0;
      const auto start = edge;
      auto next = edge;
      do {
        auto face = next.first;
        std::size_t k = next.second;

        CGAL_assertion(k != std::size_t(-1));
        const auto& point = face->vertex((k + 1) % 3)->point();
        contour.push_back(point);

        std::size_t knext = (k + 2) % 3;
        auto vertex = face->vertex(knext);
        next = get_next_edge(vertex, face, ref_label, region_index);
        
        k = next.second;
        if (k == std::size_t(-1)) {
          contour.push_back(vertex->point());

          knext = (knext + 1) % 3;
          vertex = face->vertex(knext);
          next = get_next_edge(vertex, face, ref_label, region_index);

          k = next.second;
          if (k == std::size_t(-1)) {
            contour.clear(); return;
          }
        }
        ++count;
      } while (next != start && count <= 100000);

      if (count >= 100000) {
        contour.clear(); return;
      }
    }

    std::pair<Face_handle, std::size_t> get_next_edge(
      const Vertex_handle vertex,
      const Face_handle ref_face,
      const std::size_t ref_label,
      const std::size_t region_index) const {

      auto face = m_triangulation.incident_faces(vertex);
      if (face.is_empty()) 
        return std::make_pair(Face_handle(), std::size_t(-1));

      std::size_t count = 0;
      const auto start = face;
      do {
        if (
          face != ref_face && 
          face->info().object_index == region_index) {
          
          for (std::size_t k = 0; k < 3; ++k) {
            const auto nface = face->neighbor(k);
            if (
              face->vertex(k) != vertex &&
              nface->info().label == ref_label) {
              
              return std::make_pair(face, k);
            }
          }
        }
        ++face; ++count;
      } while (face != start && count <= 100);
      return std::make_pair(Face_handle(), std::size_t(-1));
    }

    void extract_holes(
      const std::vector<Indices>& regions,
      std::vector<Face_handles>& holes) const {

      holes.clear();
      holes.resize(regions.size());

      for (std::size_t i = 0; i < regions.size(); ++i) {
        const auto& region = regions[i];
        for (const std::size_t face_index : region)
          holes[i].push_back(m_faces[face_index]);
      }
    }
  };

} // internal
} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_INTERNAL_BOUNDARY_FROM_TRIANGULATION_H
