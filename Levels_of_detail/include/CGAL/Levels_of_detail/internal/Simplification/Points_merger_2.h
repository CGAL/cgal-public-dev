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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_POINTS_MERGER_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_POINTS_MERGER_2_H

// STL includes.
#include <map>
#include <utility>
#include <iostream>
#include <vector>
#include <string>
#include <list>
#include <limits>
#include <set>
#include <algorithm>
#include <iterator>

// CGAL includes.
#include <CGAL/barycenter.h>
#include <CGAL/property_map.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/point_generators_2.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Spatial search.
#include <CGAL/Levels_of_detail/internal/Spatial_search/K_neighbor_query.h>
#include <CGAL/Levels_of_detail/internal/Spatial_search/Sphere_neighbor_query.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"
#include "../../../../../test/Levels_of_detail/include/Utilities.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Points_merger_2 {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Segment_2 = typename Traits::Segment_2;
    using Line_2 = typename Traits::Line_2;

    using Pair = std::pair<Point_2, std::size_t>;
    using Point_map = CGAL::First_of_pair_property_map<Pair>;

    using K_neighbor_query =
      internal::K_neighbor_query<Traits, std::vector<Pair>, Point_map>;
    using Sphere_neighbor_query =
      internal::Sphere_neighbor_query<Traits, std::vector<Pair>, Point_map>;

    using Indices = std::vector<std::size_t>;

    using Triangulation = internal::Triangulation<Traits>;
    
    using FBI = typename Triangulation::CFB;
    using FB  = CGAL::Alpha_shape_face_base_2<Traits, FBI>;
    using VBI = typename Triangulation::VBI;
    using VB  = CGAL::Alpha_shape_vertex_base_2<Traits, VBI>;
    using TDS = CGAL::Triangulation_data_structure_2<VB, FB>;
    using TAG = CGAL::Exact_predicates_tag;
    using Delaunay = CGAL::Constrained_Delaunay_triangulation_2<Traits, TDS, TAG>;
    using Alpha_shape_2 = CGAL::Alpha_shape_2<Delaunay>;

    using Location_type = typename Alpha_shape_2::Locate_type;
    using Face_handle = typename Alpha_shape_2::Face_handle;
    using Vertex_handle = typename Alpha_shape_2::Vertex_handle;

    Points_merger_2(
      const FT noise_level,
      const FT alpha) :
    m_noise_level(noise_level),
    m_alpha(alpha) { 

      CGAL_precondition(m_alpha > FT(0));
    }

    template<typename Point_map>
    void merge(
      const Indices& indices,
      const Point_map& point_map,
      const std::vector< std::vector<Segment_2> >& contours,
      Triangulation& result) {
      
      Delaunay triangulation;
      insert_points(indices, point_map, triangulation);
      Alpha_shape_2 alpha_shape(
        triangulation, m_alpha, Alpha_shape_2::GENERAL);
      tag_faces(alpha_shape);
      save_triangulation(alpha_shape, 
        "/Users/monet/Documents/lod/logs/buildings/tmp/alpha_shape-original", false);
      
      auto& delaunay = result.delaunay;
      convert(indices, point_map, contours, alpha_shape, delaunay);
      save_triangulation(delaunay, 
        "/Users/monet/Documents/lod/logs/buildings/tmp/delaunay-original", false);

      use_graphcut(delaunay);
      save_triangulation(delaunay, 
        "/Users/monet/Documents/lod/logs/buildings/tmp/delaunay-clean", false);
    }

  private:
    const FT m_noise_level;
    const FT m_alpha;

    template<
    typename Range, 
    typename Point_map,
    typename Base>
    void insert_points(
      const Range& range, 
      const Point_map& point_map,
      Base& base) {

      for (auto it = range.begin(); it != range.end(); ++it)
        base.insert(
          internal::point_2_from_point_3(get(point_map, *it)));
    }

    template<typename Base>
    void insert_constraints(
      const std::vector< std::vector<Segment_2> >& contours,
      Base& base) {

      for (const auto& contour : contours) {
        for (const auto& segment : contour) {
          const auto& s = segment.source();
          const auto& t = segment.target();

          const auto vh1 = base.insert(s);
          const auto vh2 = base.insert(t);

          if (vh1 != vh2)
            base.insert_constraint(vh1, vh2);
        }
      }
    }

    void tag_faces(
      Alpha_shape_2& alpha_shape) {

      for (auto fit = alpha_shape.finite_faces_begin();
      fit != alpha_shape.finite_faces_end(); ++fit)
        fit->info().tagged = true;

      for (auto fh = alpha_shape.finite_faces_begin();
      fh != alpha_shape.finite_faces_end(); ++fh) {
        
        bool found = false;
        for (std::size_t k = 0; k < 3; ++k) {
          const auto fhn = fh->neighbor(k);
          if (alpha_shape.is_infinite(fhn)) {
            found = true; break;
          }
        }
        if (!found) continue;
        if (!fh->info().tagged) continue;

        Face_handle seed = static_cast<Face_handle>(fh);
        propagate(alpha_shape, seed);
      }
    }

    void propagate(
      const Alpha_shape_2& alpha_shape, 
      Face_handle fh) {

      if (alpha_shape.classify(fh) == Alpha_shape_2::INTERIOR)
        return;
      fh->info().tagged = false;

      for (std::size_t k = 0; k < 3; ++k) {
        auto fhn = fh->neighbor(k);
        if (!alpha_shape.is_infinite(fhn) && fhn->info().tagged)
          propagate(alpha_shape, fhn);
      }
    }

    template<typename Base>
    void save_triangulation(
      const Base& base,
      const std::string path,
      const bool out_labels) {

      const FT z = FT(0);
      std::size_t num_vertices = 0;
      internal::Indexer<Point_3> indexer;

      std::vector<Point_3> vertices; 
      std::vector<Indices> faces; 
      std::vector<Color> fcolors;

      Polygon_inserter<Traits> inserter(faces, fcolors);
      auto output_vertices = std::back_inserter(vertices);
      auto output_faces = boost::make_function_output_iterator(inserter);

      output_triangulation(
        base, indexer, num_vertices, 
        output_vertices, output_faces, z, out_labels);
      
      Saver<Traits> saver;
      saver.export_polygon_soup(vertices, faces, fcolors, path);
    }

    template<
    typename Base,
    typename Indexer,
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_triangulation(
      const Base& base,
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const FT z,
      const bool out_labels) const {

      std::vector<std::size_t> face(3);
      for (auto fh = base.finite_faces_begin(); 
      fh != base.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;
        
        for (std::size_t k = 0; k < 3; ++k) {
          const Point_2& q = fh->vertex(k)->point();
          const Point_3 p = Point_3(q.x(), q.y(), z);
          const std::size_t idx = indexer(p);
          if (idx == num_vertices) {
            *(vertices++) = p; 
            ++num_vertices;
          }
          face[k] = idx;
        }
        if (out_labels)
          *(faces++) = std::make_pair(face, fh->info().label);
        else 
          *(faces++) = std::make_pair(face, 1);
      }
    }

    template<
    typename Point_map,
    typename Base>
    void convert(
      const Indices& indices,
      const Point_map& point_map,
      const std::vector< std::vector<Segment_2> >& contours,
      const Alpha_shape_2& alpha_shape,
      Base& base) {

      base.clear();
      insert_points(indices, point_map, base);
      insert_constraints(contours, base);
      update_tagged_faces(alpha_shape, base);
    }

    template<typename Base>
    void update_tagged_faces(
      const Alpha_shape_2& alpha_shape,
      Base& base) {

      for (auto fh = base.finite_faces_begin();
      fh != base.finite_faces_end(); ++fh) {
        
        const Point_2 b = CGAL::barycenter(
          fh->vertex(0)->point(), FT(1),
          fh->vertex(1)->point(), FT(1),
          fh->vertex(2)->point(), FT(1));
        
        Location_type type; int stub;
        const auto bh = alpha_shape.locate(b, type, stub);

        if (bh->info().tagged)
          fh->info().tagged = true;
        else 
          fh->info().tagged = false;
      }
    }

    template<typename Base>
    void use_graphcut(
      Base& base) {

    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_POINTS_MERGER_2_H
