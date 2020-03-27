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

#ifndef CGAL_URBAN_AREA_PROCESSING_INTERNAL_MERGE_CONTOURS_2_H
#define CGAL_URBAN_AREA_PROCESSING_INTERNAL_MERGE_CONTOURS_2_H

// #include <CGAL/license/Urban_area_processing.h>

// STL includes.
#include <vector>
#include <utility>
#include <algorithm>

// Boost includes.
#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#include <CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// Internal includes.
#include <CGAL/Urban_area_processing/internal/utils.h>

// Utils.
#include "../../../../../test/Urban_area_processing/include/Saver.h"

// TODO:
// 1. Should I remove alpha shape computation at the beginning of the algorithm?
// 2. I should probably extract a function that tags all exterior faces of the alpha shape.
// It is also required in the class Boundary_from_triangulation_2.h.
// 3. Clean up stupid face handle types.
// 4. Try to merge update_tags() with mark_outer_faces().

namespace CGAL {
namespace Urban_area_processing {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename SegmentMap>
  class Merge_contours_2 {
  
  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;
    using Triangle_2 = typename Traits::Triangle_2;

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

    using Indices = std::vector<std::size_t>;
    using Doubles = std::vector<double>;
    using Size_pair = std::pair<std::size_t, std::size_t>;

    using Alpha_expansion = CGAL::internal::Alpha_expansion_graph_cut_boost;

    Merge_contours_2(
      const Input_range& input_range,
      const Segment_map segment_map,
      const FT scale,
      const FT noise,
      const bool verbose = true) :
    m_input_range(input_range),
    m_segment_map(segment_map),
    m_scale(scale),
    m_noise(noise),
    m_verbose(verbose) { 

      CGAL_precondition(input_range.size() > 0);
    }

    template<
    typename PointRange,
    typename PointMap,
    typename Triangulation>
    void merge(
      const PointRange& point_range,
      const PointMap point_map,
      Triangulation& triangulation) const {

      CGAL_precondition(point_range.size() >= 3);
      create_initial_triangulation(
        point_range, point_map, triangulation);
      if (m_verbose)
        std::cout << "- initial triangulation is created" << std::endl;

      save_triangulation(triangulation, "triangulation_initial");
      if (m_verbose)
        std::cout << "- initial triangulation is saved" << std::endl;

      use_graphcut(triangulation);
      if (m_verbose)
        std::cout << "- final triangulation is created" << std::endl;

      save_triangulation(triangulation, "triangulation_final");
      if (m_verbose)
        std::cout << "- final triangulation is saved" << std::endl;
    }

  private:
    const Input_range& m_input_range;
    const Segment_map m_segment_map;
    const FT m_scale;
    const FT m_noise;
    const bool m_verbose;

    template<
    typename PointRange,
    typename PointMap,
    typename Triangulation>
    void create_initial_triangulation(
      const PointRange& point_range,
      const PointMap point_map,
      Triangulation& triangulation) const {
      
      Delaunay delaunay;
      insert_in_triangulation(point_range, point_map, delaunay);
      Alpha_shape_2 alpha_shape(
        delaunay, m_scale, Alpha_shape_2::GENERAL);
      mark_outer_faces(alpha_shape);
      convert(point_range, point_map, alpha_shape, 
      triangulation);
    }

    template<
    typename PointRange,
    typename PointMap,
    typename Base>
    void insert_in_triangulation(
      const PointRange& point_range,
      const PointMap point_map,
      Base& base) const {
      
      base.clear();
      for (auto item = point_range.begin(); 
      item != point_range.end(); ++item)
        base.insert(internal::point_2_from_point_3(
          get(point_map, *item)));
    }

    void mark_outer_faces(
      Alpha_shape_2& alpha_shape) const {

      for (auto face = alpha_shape.finite_faces_begin();
      face != alpha_shape.finite_faces_end(); ++face)
        face->info().tagged = true; // means inside

      for (auto face = alpha_shape.finite_faces_begin();
      face != alpha_shape.finite_faces_end(); ++face) {
        
        bool found = false;
        for (std::size_t k = 0; k < 3; ++k) {
          const auto nface = face->neighbor(k);
          if (alpha_shape.is_infinite(nface)) {
            found = true; break;
          }
        }
        if (!found) continue;
        if (!face->info().tagged) continue;
        propagate(alpha_shape, face);
      }
    }

    void propagate(
      const Alpha_shape_2& alpha_shape, 
      const Face_handle face) const {

      if (alpha_shape.classify(face) == Alpha_shape_2::INTERIOR) 
        return;
      face->info().tagged = false;

      for (std::size_t k = 0; k < 3; ++k) {
        const auto nface = face->neighbor(k);
        if (!alpha_shape.is_infinite(nface) && nface->info().tagged)
          propagate(alpha_shape, nface);
      }
    }

    template<
    typename PointRange,
    typename PointMap,
    typename Triangulation>
    void convert(
      const PointRange& point_range,
      const PointMap point_map,
      const Alpha_shape_2& alpha_shape,
      Triangulation& triangulation) const {

      triangulation.clear();
      insert_in_triangulation(
        point_range, point_map, triangulation);
      insert_constraints(triangulation);
      update_tagged_faces(alpha_shape, triangulation);
    }

    template<typename Base>
    void insert_constraints(
      Base& base) const {

      for (const auto& item : m_input_range) {
        const auto& segment = get(m_segment_map, item);
        const auto& source = segment.source();
        const auto& target = segment.target();

        const auto vh1 = base.insert(source);
        const auto vh2 = base.insert(target);
        if (vh1 != vh2)
          base.insert_constraint(vh1, vh2);
      }
    }

    template<typename Triangulation>
    void update_tagged_faces(
      const Alpha_shape_2& alpha_shape,
      Triangulation& triangulation) const {

      Point_2 center;
      for (auto face = triangulation.finite_faces_begin();
      face != triangulation.finite_faces_end(); ++face) {
        internal::compute_barycenter(face, center);
        const auto fh = alpha_shape.locate(center);
        if (fh->info().tagged)
          face->info().tagged = true;
        else 
          face->info().tagged = false;
      }
    }

    template<typename Triangulation>
    void save_triangulation(
      const Triangulation& triangulation,
      const std::string name) const {

      Saver<Traits> saver;
      saver.export_polygon_soup(triangulation,
      "/Users/monet/Documents/gf/urban-area-processing/logs/" + name);
    }

    template<typename Triangulation>
    void use_graphcut(
      Triangulation& triangulation) const {

      preprocess_triangulation(triangulation);
      save_triangulation(triangulation, "delaunay_original");

      estimate_labels(triangulation);
      save_triangulation(triangulation, "delaunay_approx");

      apply_graphcut(triangulation);
      save_triangulation(triangulation, "delaunay_graphcut");
    }

    template<typename Triangulation>
    void preprocess_triangulation(
      Triangulation& triangulation) const {

      set_object_indices(triangulation);
      clear_labels(triangulation);
      resize_probabilities(triangulation);
    }

    template<typename Triangulation>
    void set_object_indices(
      Triangulation& triangulation) const {

      std::size_t count = 0;
      for (auto face = triangulation.finite_faces_begin();
      face != triangulation.finite_faces_end(); ++face) {
        if (face->info().tagged) {
          face->info().index = count; ++count;
        }
      }
    }

    template<typename Triangulation>
    void clear_labels(
      Triangulation& triangulation) const {

      for (auto face = triangulation.finite_faces_begin();
      face != triangulation.finite_faces_end(); ++face) {
        if (face->info().tagged) 
          face->info().label = 1; // inside
        else 
          face->info().label = 0; // outside
      }
    }

    template<typename Triangulation>
    void resize_probabilities(
      Triangulation& triangulation) const {

      for (auto face = triangulation.finite_faces_begin();
      face != triangulation.finite_faces_end(); ++face) {
        if (!face->info().tagged) continue;
        face->info().probabilities.clear();
        face->info().probabilities.resize(2, FT(0)); // two labels
      }
    }

    template<typename Triangulation>
    void estimate_labels(
      Triangulation& triangulation) const {

      estimate_in_out(triangulation);
      update_labels(triangulation);
    }

    template<typename Triangulation>
    void estimate_in_out(
      Triangulation& triangulation) const {

      estimate_in_out_from_boundary_faces(triangulation);
      normalize_probabilities(triangulation);
    }

    template<typename Triangulation>
    void estimate_in_out_from_boundary_faces(
      Triangulation& triangulation) const {

      using FaceHandle = typename Triangulation::Face_handle;
      for (auto face = triangulation.finite_faces_begin();
      face != triangulation.finite_faces_end(); ++face) {
        if (!face->info().tagged) continue;
        
        for (std::size_t k = 0; k < 3; ++k) {
          const auto nface = face->neighbor(k);
          if (nface->info().tagged) continue;
          add_statistics_from_boundary_face(
            static_cast<FaceHandle>(face), triangulation);
          break;
        }
      }
    }

    template<
    typename FaceHandle,
    typename Triangulation>
    void add_statistics_from_boundary_face(
      const FaceHandle face, Triangulation& triangulation) const {

      Point_2 center;
      const FT radius = FT(1);
      const std::size_t num_samples = 24;
      internal::compute_barycenter(face, center);

      std::vector<Point_2> samples1;
      internal::create_points_on_circle(
        center, radius, FT(0), num_samples, samples1);

      std::vector<Point_2> samples2;
      internal::create_points_on_circle(
        center, radius, FT(180), num_samples, samples2);

      for (std::size_t k = 0; k < num_samples / 2; ++k) {
        const auto& q1 = samples1[k];
        const auto& q2 = samples2[k];
        apply_line_walk_from_boundary_face(
          center, q1, face, triangulation);
        apply_line_walk_from_boundary_face(
          center, q2, face, triangulation);
      }
    }

    template<
    typename FaceHandle,
    typename Triangulation>
    void apply_line_walk_from_boundary_face(
      const Point_2& p,
      const Point_2& q,
      const FaceHandle ref,
      Triangulation& triangulation) const {

      std::vector< std::vector<FaceHandle> > regions;
      const bool success = apply_line_walk(
        p, q, ref, triangulation, regions);
      if (success)
        add_region_statistics(p, q, regions);
    }

    template<
    typename FaceHandle,
    typename Triangulation>
    bool apply_line_walk(
      const Point_2& p,
      const Point_2& q,
      const FaceHandle ref,
      const Triangulation& triangulation,
      std::vector< std::vector<FaceHandle> >& regions) const {

      if (triangulation.oriented_side(ref, p) == CGAL::ON_NEGATIVE_SIDE)
        return false;

      auto circ = triangulation.line_walk(p, q, ref);
      const auto end = circ;

      regions.clear();
      std::vector<FaceHandle> region;
      region.push_back(circ);
      do {

        auto f1 = circ; ++circ;
        auto f2 = circ;

        if (triangulation.is_infinite(f2) || !f2->info().tagged) {
          if (!region.empty()) 
            regions.push_back(region); 
          break;
        }

        const bool success = internal::are_neighbors(f1, f2);
        if (success) {
          const std::size_t idx = f1->index(f2);
          const auto edge = std::make_pair(f1, idx);

          if (triangulation.is_constrained(edge)) {
            if (!region.empty())
              regions.push_back(region); 
            region.clear();
          }
        }
        region.push_back(f2);
      } while (circ != end);
      return true;
    }

    template<typename FaceHandle>
    void add_region_statistics(
      const Point_2& p1, const Point_2& p2,
      const std::vector< std::vector<FaceHandle> >& regions) const {

      // >= 0 regions.
      const std::size_t num_regions = regions.size();
      if (num_regions == 0) return;

      // >= 1 regions.
      std::vector<FT> weights;
      compute_angle_weights(p1, p2, regions, weights);

      if (num_regions == 1) {
        std::size_t count = 0;
        const auto& region = regions[0];
        const std::size_t num_faces = region.size();
        
        Point_2 q1, q2;
        internal::compute_barycenter(region[0], q1);
        internal::compute_barycenter(region[num_faces - 1], q2);

        const FT distance = internal::distance(q1, q2);
        if (distance < m_noise) {         
          for (std::size_t i = 0; i < region.size(); ++i) {
            region[i]->info().probabilities[0] += FT(1) * weights[count]; // outside
            ++count;
          }
        } else {
          for (std::size_t i = 0; i < region.size(); ++i) {
            region[i]->info().probabilities[1] += FT(1) * weights[count]; // inside
            ++count;
          }
        } return;
      } 

      // >= 2 regions.
      std::size_t count = 0;
      for (std::size_t k = 0; k < num_regions - 1; ++k) {
        const auto& region = regions[k];
        for (std::size_t i = 0; i < region.size(); ++i) {
          region[i]->info().probabilities[1] += FT(1) * weights[count]; // inside
          ++count;
        }
      }

      const auto& region = regions[num_regions - 1];
      for (std::size_t i = 0; i < region.size(); ++i) {
        region[i]->info().probabilities[0] += FT(1) * weights[count]; // outside
        ++count;
      }
    }

    template<typename FaceHandle>
    void compute_angle_weights(
      const Point_2& p1, const Point_2& p2,
      const std::vector< std::vector<FaceHandle> >& regions,
      std::vector<FT>& weights) const {

      // >= 0 regions.
      std::size_t num_weights = 0;
      for (const auto& region : regions)
        num_weights += region.size();
      
      const std::size_t num_regions = regions.size();

      if (num_regions == 0)
        return;

      // >= 1 regions.
      if (num_regions == 1) {
        weights.clear();

        const auto& region = regions[0];
        const std::size_t num_faces = region.size();

        Point_2 q1, q2;
        internal::compute_barycenter(region[0], q1);
        internal::compute_barycenter(region[num_faces - 1], q2);

        const FT distance = internal::distance(q1, q2);
        if (distance < m_noise)
          weights.resize(num_weights, FT(1));
        else
          weights.resize(num_weights, FT(0));
        return;
      }

      // >= 2 regions.
      FT product = FT(1);
      const Segment_2 segment1 = Segment_2(p1, p2);
      
      std::vector<FT> angles;
      for (std::size_t i = 0; i < num_regions - 1; ++i) {
        const std::size_t ip = i + 1;

        const auto& region0 = regions[i];
        const auto& region1 = regions[ip];

        const auto face  = region0[region0.size() - 1];
        const auto nface = region1[0];

        const std::size_t idx = face->index(nface);
        const auto& vh1 = face->vertex( (idx + 1) % 3 );
        const auto& vh2 = face->vertex( (idx + 2) % 3 );

        const auto& q1 = vh1->point();
        const auto& q2 = vh2->point();

        const Segment_2 segment2 = Segment_2(q1, q2);

        const FT angle_d = internal::angle_degree_2(segment1, segment2);
        const FT angle_2 = CGAL::abs(internal::get_angle_2(angle_d));
        
        const FT angle = angle_2 * static_cast<FT>(CGAL_PI) / FT(180);
        angles.push_back(angle);
      }

      FT denom = -FT(1);
      for (const auto& angle : angles)
        denom = CGAL::max(denom, angle);

      for (const auto& angle : angles)
        product *= internal::smooth_step(angle, denom);

      weights.clear();
      weights.resize(num_weights, product);
    }

    template<typename Triangulation>
    void normalize_probabilities(
      Triangulation& triangulation) const {

      for (auto face = triangulation.finite_faces_begin();
      face != triangulation.finite_faces_end(); ++face) {
        if (!face->info().tagged) continue;

        FT& inside  = face->info().probabilities[1];
        FT& outside = face->info().probabilities[0];

        if (inside == FT(0) && outside == FT(0)) {
          inside = FT(1) / FT(2); outside = FT(1) / FT(2);
          continue;
        }

        const FT sum = inside + outside;
        inside /= sum; outside /= sum;
      }
    }

    template<typename Triangulation>
    void update_labels(
      Triangulation& triangulation) const {

      for (auto face = triangulation.finite_faces_begin();
      face != triangulation.finite_faces_end(); ++face) {
        if (!face->info().tagged) continue;

        const auto& probabilities = face->info().probabilities;
        if (probabilities[1] >= FT(1) / FT(2))
          face->info().label = 1; // inside
        else 
          face->info().label = 0; // outside
      }
    }

    template<typename Triangulation>
    void apply_graphcut(
      Triangulation& triangulation) const {

      Indices labels;
      set_initial_labels(triangulation, labels);

      std::vector<Size_pair> edges;
      Doubles edge_weights;
      set_graphcut_edges(triangulation, edges, edge_weights);

      std::vector<Doubles> cost_matrix;
      set_cost_matrix(triangulation, cost_matrix);

      Alpha_expansion graphcut;
      graphcut(edges, edge_weights, cost_matrix, labels);

      set_new_labels(labels, triangulation);
      update_tags(triangulation);
    }

    template<typename Triangulation>
    void set_initial_labels(
      const Triangulation& triangulation,
      Indices& labels) const {

      labels.clear();
      for (auto face = triangulation.finite_faces_begin();
      face != triangulation.finite_faces_end(); ++face) {
        if (!face->info().tagged) continue;
        labels.push_back(face->info().label);
      }
    }

    template<typename Triangulation>
    void set_graphcut_edges(
      const Triangulation& triangulation,
      std::vector<Size_pair>& edges,
      Doubles& edge_weights) const {

      edges.clear();
      edge_weights.clear();

      FT max_value = -FT(1);
      for (auto eh = triangulation.finite_edges_begin();
      eh != triangulation.finite_edges_end(); ++eh) {
        const auto& edge = *eh;

        const auto face = edge.first;
        const auto idx = edge.second;
        const auto nface = face->neighbor(idx);

        if (face->info().tagged && nface->info().tagged) {
          
          const std::size_t idxi = face->info().index;
          const std::size_t idxj = nface->info().index;
          edges.push_back(std::make_pair(idxi, idxj));

          const auto& p1 = face->vertex((idx + 1) % 3)->point();
          const auto& p2 = face->vertex((idx + 2) % 3)->point();
          const FT distance = internal::distance(p1, p2); // should I change it to the squared_distance?

          FT edge_weight = distance;
          max_value = CGAL::max(max_value, edge_weight);
          if (triangulation.is_constrained(edge))
            edge_weight = -FT(1);
          edge_weights.push_back(CGAL::to_double(edge_weight));
        }
      }

      const FT beta = FT(1);
      for (auto& edge_weight : edge_weights) {
        if (edge_weight == -1.0)
          edge_weight = CGAL::to_double(max_value);

        CGAL_assertion(max_value != FT(0));
        edge_weight /= CGAL::to_double(max_value);
        edge_weight *= CGAL::to_double(beta);
      }
    }

    template<typename Triangulation>
    void set_cost_matrix(
      const Triangulation& triangulation,
      std::vector<Doubles>& cost_matrix) const {

      cost_matrix.clear();
      cost_matrix.resize(2);

      FT max_value = -FT(1);
      std::vector<FT> weights;

      for (auto face = triangulation.finite_faces_begin();
      face != triangulation.finite_faces_end(); ++face) {
        if (!face->info().tagged) continue;
        
        const auto& p0 = face->vertex(0)->point();
        const auto& p1 = face->vertex(1)->point();
        const auto& p2 = face->vertex(2)->point();
        
        const Triangle_2 triangle = Triangle_2(p0, p1, p2);
        const FT area = CGAL::abs(triangle.area());

        const FT weight = area;
        max_value = CGAL::max(max_value, weight);
        weights.push_back(weight);
      }

      CGAL_assertion(max_value != FT(0));
      for (auto& weight : weights)
        weight /= max_value;

      std::size_t count = 0;
      for (auto face = triangulation.finite_faces_begin();
      face != triangulation.finite_faces_end(); ++face) {
        if (!face->info().tagged) continue;

        for (std::size_t k = 0; k < 2; ++k)
          cost_matrix[k].push_back(get_cost(
            weights[count], face->info().probabilities[k]));
        ++count;
      }
    }

    double get_cost(
      const FT weight, const FT probability) const {
      return CGAL::to_double((FT(1) - probability) * weight);
    }

    template<typename Triangulation>
    void set_new_labels(
      const Indices& labels,
      Triangulation& triangulation) const {
      
      std::size_t count = 0;
      for (auto face = triangulation.finite_faces_begin();
      face != triangulation.finite_faces_end(); ++face) {
        if (!face->info().tagged) continue;
        face->info().label = labels[count];
        ++count;
      }
    }

    template<typename Triangulation>
    void update_tags(
      Triangulation& triangulation) const {

      for (auto face = triangulation.finite_faces_begin();
      face != triangulation.finite_faces_end(); ++face) {
        if (!face->info().tagged) continue;
        
        bool found = false;
        for (std::size_t k = 0; k < 3; ++k) {
          const auto nface = face->neighbor(k);
          if (!nface->info().tagged || triangulation.is_infinite(nface)) {
            found = true; break;
          }
        }
        if (!found) continue;

        if (face->info().label == 0) // outside
          propagate_tag(triangulation, face);
      }
    }

    template<
    typename Triangulation,
    typename FaceHandle>
    void propagate_tag(
      const Triangulation& triangulation, 
      FaceHandle face) const {

      if (face->info().label == 0)
        face->info().tagged = false;
      else return;

      for (std::size_t k = 0; k < 3; ++k) {
        const auto nface = face->neighbor(k);
        if (!triangulation.is_infinite(nface) && nface->info().tagged)
          propagate_tag(triangulation, nface);
      }
    }
  };

} // internal
} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_INTERNAL_MERGE_CONTOURS_2_H
