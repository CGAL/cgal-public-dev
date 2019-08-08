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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_ALPHA_SHAPES_FILTERING_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_ALPHA_SHAPES_FILTERING_2_H

// CGAL includes.
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/point_generators_2.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Alpha_shapes_filtering_2 {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Triangle_2 = typename Traits::Triangle_2;

    using Vb = CGAL::Alpha_shape_vertex_base_2<Traits>;
    using Fb = CGAL::Alpha_shape_face_base_2<Traits>;
    using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
    using Triangulation_2 = CGAL::Delaunay_triangulation_2<Traits, Tds>;
    using Alpha_shape_2 = CGAL::Alpha_shape_2<Triangulation_2>;
    using Points_in_triangle = CGAL::Random_points_in_triangle_2<Point_2>;

    Alpha_shapes_filtering_2(const FT alpha) : 
    m_alpha(alpha) 
    { }

    template<
    typename Range, 
    typename Point_map>
    void add_points(const Range& range, Point_map point_map) {
      insert_in_triangulation(range, point_map);
    }

    void get_filtered_points(
      const FT edge_sampling_2, 
      std::vector<Point_2>& result) {
      
      if (m_triangulation.number_of_faces() == 0) return;
      CGAL_precondition(m_alpha > FT(0));
      CGAL_precondition(m_triangulation.number_of_faces() != 0);

      Alpha_shape_2 alpha_shape(m_triangulation, m_alpha, Alpha_shape_2::GENERAL);
      sample_edges(alpha_shape, edge_sampling_2, result);
    }

    void get_samples(
      const FT edge_sampling_2,
      const std::size_t face_sampling_2, 
      std::vector<Point_2>& result) {
      
      if (m_triangulation.number_of_faces() == 0) return;
      CGAL_precondition(m_alpha > FT(0));
      CGAL_precondition(m_triangulation.number_of_faces() != 0);

      Alpha_shape_2 alpha_shape(m_triangulation, m_alpha, Alpha_shape_2::GENERAL);
      sample_faces(alpha_shape, edge_sampling_2, face_sampling_2, result);
    }

  private:
    const FT m_alpha;
    Triangulation_2 m_triangulation;

    template<
    typename Range, 
    typename Point_map>
    void insert_in_triangulation(
      const Range& range, 
      Point_map point_map) {
                  
      for (auto it = range.begin(); it != range.end(); ++it)
        m_triangulation.insert(
          internal::point_2_from_point_3(get(point_map, *it)));
    }

    void sample_edges(
      const Alpha_shape_2& alpha_shape,
      const FT edge_sampling_2,
      std::vector<Point_2>& result) const {

      for (auto eit = alpha_shape.alpha_shape_edges_begin(); 
        eit != alpha_shape.alpha_shape_edges_end(); ++eit) {

        const Point_2& source = eit->first->vertex((eit->second + 1) % 3)->point();
        const Point_2& target = eit->first->vertex((eit->second + 2) % 3)->point();
        sample_edge(source, target, edge_sampling_2, result);
      }
    }

    void sample_edge(
      const Point_2& source, 
      const Point_2& target,
      const FT edge_sampling_2,
      std::vector<Point_2>& result) const {

      const FT distance = internal::distance(source, target);
        
      CGAL_precondition(edge_sampling_2 > FT(0));
      const std::size_t nb_pts = 
      static_cast<std::size_t>(CGAL::to_double(distance / edge_sampling_2)) + 1;
        
      CGAL_precondition(nb_pts > 0);
      for (std::size_t i = 0; i <= nb_pts; ++i) {

        const FT ratio = static_cast<FT>(i) / static_cast<FT>(nb_pts);
        result.push_back(
          Point_2(
            source.x() * (FT(1) - ratio) + target.x() * ratio,
            source.y() * (FT(1) - ratio) + target.y() * ratio));
      }
    }

    void sample_faces(
      const Alpha_shape_2& alpha_shape,
      const FT edge_sampling_2,
      const std::size_t face_sampling_2,
      std::vector<Point_2>& result) const {

      for (auto fit = alpha_shape.finite_faces_begin(); 
        fit != alpha_shape.finite_faces_end(); ++fit) {

        const auto type = alpha_shape.classify(fit);
        if (type == Alpha_shape_2::INTERIOR) {
          
          const auto& p0 = fit->vertex(0)->point();
          const auto& p1 = fit->vertex(1)->point();
          const auto& p2 = fit->vertex(2)->point();

          // Add edges.
          sample_edge(p0, p1, edge_sampling_2, result);
          sample_edge(p1, p2, edge_sampling_2, result);
          sample_edge(p2, p0, edge_sampling_2, result);

          // Add face.
          const Triangle_2 triangle(p0, p1, p2);
          Points_in_triangle generator(triangle);
          std::copy_n(generator, face_sampling_2, 
          std::back_inserter(result));
        }
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_ALPHA_SHAPES_FILTERING_2_H
