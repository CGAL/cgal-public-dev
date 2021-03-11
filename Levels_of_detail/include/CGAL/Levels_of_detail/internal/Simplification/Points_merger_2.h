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

// Boost includes.
#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#include <CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h>

// CGAL includes.
#include <CGAL/enum.h>
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
    using Vector_2 = typename Traits::Vector_2;
    using Line_2 = typename Traits::Line_2;
    using Triangle_2 = typename Traits::Triangle_2;

    using Pair = std::pair<Point_2, std::size_t>;
    using Point_map = CGAL::First_of_pair_property_map<Pair>;

    using K_neighbor_query =
      internal::K_neighbor_query<Traits, std::vector<Pair>, Point_map>;
    using Sphere_neighbor_query =
      internal::Sphere_neighbor_query<Traits, std::vector<Pair>, Point_map>;

    using Indices = std::vector<std::size_t>;

    using Triangulation = internal::Triangulation<Traits>;
    using BaseTri = typename Triangulation::Delaunay;
    using LF_circulator = typename BaseTri::Line_face_circulator;
    using F_handle = typename BaseTri::Face_handle;

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

    using Size_pair = std::pair<std::size_t, std::size_t>;
    using Alpha_expansion = CGAL::internal::Alpha_expansion_graph_cut_boost;
    using Edge = std::pair<F_handle, std::size_t>;

    Points_merger_2(
      const FT noise_level,
      const FT alpha) :
    m_noise_level(noise_level),
    m_alpha(alpha),
    m_pi(static_cast<FT>(CGAL_PI)) {

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
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/alpha_shape-original", false);

      BaseTri& delaunay = result.delaunay;
      convert(indices, point_map, contours, alpha_shape, delaunay);
      use_graphcut(delaunay);
      save_triangulation(delaunay,
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/delaunay-clean", false);
    }

  private:
    const FT m_noise_level;
    const FT m_alpha;
    const FT m_pi;

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
        propagate_alpha_shape(alpha_shape, fh);
      }
    }

    void propagate_alpha_shape(
      const Alpha_shape_2& alpha_shape,
      Face_handle fh) {

      if (alpha_shape.classify(fh) == Alpha_shape_2::INTERIOR)
        return;
      fh->info().tagged = false;

      for (std::size_t k = 0; k < 3; ++k) {
        auto fhn = fh->neighbor(k);
        if (!alpha_shape.is_infinite(fhn) && fhn->info().tagged)
          propagate_alpha_shape(alpha_shape, fhn);
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
    typename Point_map>
    void convert(
      const Indices& indices,
      const Point_map& point_map,
      const std::vector< std::vector<Segment_2> >& contours,
      const Alpha_shape_2& alpha_shape,
      BaseTri& base) {

      base.clear();
      insert_points(indices, point_map, base);
      insert_constraints(contours, base);
      update_tagged_faces(alpha_shape, base);
    }

    void update_tagged_faces(
      const Alpha_shape_2& alpha_shape,
      BaseTri& base) {

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

    void use_graphcut(
      BaseTri& base) {

      set_object_indices(base);

      clear_labels(base);

      resize_probabilities(base);

      save_triangulation(base,
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/delaunay-original", true);

      compute_in_out(base);

      update_labels(base);

      save_triangulation(base,
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/delaunay-approx", true);

      const FT beta = FT(1);

      std::vector<std::size_t> labels;
      set_initial_labels(base, labels);

      std::vector<Size_pair> edges;
      std::vector<double> edge_weights;
      set_graphcut_edges(beta, base, edges, edge_weights);

      std::vector< std::vector<double> > cost_matrix;
      set_cost_matrix(base, cost_matrix);

      Alpha_expansion graphcut;
      graphcut(edges, edge_weights, cost_matrix, labels);

      set_new_labels(labels, base);

      save_triangulation(base,
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/delaunay-graphcut", true);

      update_tags(base);
    }

    void set_object_indices(
      BaseTri& base) {

      std::size_t count = 0;
      for (auto fh = base.finite_faces_begin();
      fh != base.finite_faces_end(); ++fh) {
        if (fh->info().tagged) {
          fh->info().object_index = count;
          ++count;
        }
      }
    }

    void clear_labels(
      BaseTri& base) {

      for (auto fh = base.finite_faces_begin();
      fh != base.finite_faces_end(); ++fh) {
        if (fh->info().tagged)
          fh->info().label = 1;
        else
          fh->info().label = 0;
      }
    }

    void resize_probabilities(
      BaseTri& base) {

      for (auto fh = base.finite_faces_begin();
      fh != base.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;

        fh->info().probabilities.clear();
        fh->info().probabilities.resize(2, FT(0));
      }
    }

    void compute_in_out(
      BaseTri& base) {

      /*
      compute_in_out_stable(base);
      return; */

      add_in_out_from_alpha_shape_boundary(base);
      /* add_in_out_from_detected_boundary(base); */
      normalize_probabilities(base);
    }

    void normalize_probabilities(
      BaseTri& base) {

      for (auto fh = base.finite_faces_begin();
      fh != base.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;

        FT& inside  = fh->info().probabilities[1];
        FT& outside = fh->info().probabilities[0];

        if (inside == FT(0) && outside == FT(0)) {
          inside = FT(1) / FT(2); outside = FT(1) / FT(2);
          continue;
        }

        const FT sum = inside + outside;
        inside /= sum; outside /= sum;
      }
    }

    void add_in_out_from_alpha_shape_boundary(
      BaseTri& base) {

      for (auto fh = base.finite_faces_begin();
      fh != base.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;

        for (std::size_t k = 0; k < 3; ++k) {
          const auto fhn = fh->neighbor(k);
          if (fhn->info().tagged) continue;

          const std::size_t idx = fh->index(fhn);
          const Edge edge = std::make_pair(fh, idx);
          add_statistics_from_alpha_shape_boundary(edge, base);
        }
      }
    }

    void add_in_out_from_detected_boundary(
      BaseTri& base) {

      for (auto fh = base.finite_faces_begin();
      fh != base.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;

        for (std::size_t k = 0; k < 3; ++k) {
          const auto fhn = fh->neighbor(k);
          if (!fhn->info().tagged) continue;

          const std::size_t idx = fh->index(fhn);
          const Edge edge = std::make_pair(fh, idx);
          add_statistics_from_detected_boundary(edge, base);
        }
      }
    }

    void add_statistics_from_alpha_shape_boundary(
      const Edge& edge,
      BaseTri& base) {

      const auto fh = edge.first;
      const FT radius = FT(1);
      const std::size_t num_samples = 24;
      const Point_2 center = CGAL::barycenter(
        fh->vertex(0)->point(), FT(1),
        fh->vertex(1)->point(), FT(1),
        fh->vertex(2)->point(), FT(1));

      std::vector<Point_2> samples1;
      create_points_on_circle(
        center, radius, FT(0), num_samples, samples1);

      std::vector<Point_2> samples2;
      create_points_on_circle(
        center, radius, FT(180), num_samples, samples2);

      for (std::size_t k = 0; k < num_samples / 2; ++k) {
        const auto& p1 = samples1[k];
        const auto& p2 = samples2[k];

        apply_line_walk_from_alpha_shape_boundary(
          center, p1, fh, base);
        apply_line_walk_from_alpha_shape_boundary(
          center, p2, fh, base);
      }
    }

    void add_statistics_from_detected_boundary(
      const Edge& edge,
      BaseTri& base) {

      if (!base.is_constrained(edge))
        return;

      const auto fh = edge.first;
      const FT radius = FT(1);
      const std::size_t num_samples = 24;
      const Point_2 center = CGAL::barycenter(
        fh->vertex(0)->point(), FT(1),
        fh->vertex(1)->point(), FT(1),
        fh->vertex(2)->point(), FT(1));

      std::vector<Point_2> samples1;
      create_points_on_circle(
        center, radius, FT(0), num_samples, samples1);

      std::vector<Point_2> samples2;
      create_points_on_circle(
        center, radius, FT(180), num_samples, samples2);

      for (std::size_t k = 0; k < num_samples / 2; ++k) {
        const auto& p1 = samples1[k];
        const auto& p2 = samples2[k];

        apply_line_walk_from_detected_boundary(
          center, p1, fh, base);
        apply_line_walk_from_detected_boundary(
          center, p2, fh, base);
      }
    }

    void apply_line_walk_from_alpha_shape_boundary(
      const Point_2& p,
      const Point_2& q,
      const F_handle ref,
      BaseTri& base) {

      std::vector< std::vector<F_handle> > regions;
      const bool success = add_line_walk(p, q, ref, base, regions);
      if (success)
        add_region_alpha_shape_v1(p, q, regions);
    }

    void add_region_alpha_shape_v1(
      const Point_2& p1, const Point_2& p2,
      const std::vector< std::vector<F_handle> >& regions) {

      const std::size_t num_regions = regions.size();
      if (num_regions == 0) return;

      std::vector<FT> weights;
      compute_angle_weights(p1, p2, regions, weights);

      if (num_regions == 1) {

        std::size_t count = 0;
        const auto& region = regions[0];
        const std::size_t num_faces = region.size();

        const auto q1 = get_point(region[0]);
        const auto q2 = get_point(region[num_faces - 1]);

        const FT distance = internal::distance(q1, q2);
        if (distance < m_noise_level) {
          for (std::size_t i = 0; i < region.size(); ++i) {
            region[i]->info().probabilities[0] += FT(1) * weights[count];
            ++count;
          }
        } else {
          for (std::size_t i = 0; i < region.size(); ++i) {
            region[i]->info().probabilities[1] += FT(1) * weights[count];
            ++count;
          }
        }

        return;
      }

      if (num_regions == 2) {

        std::size_t rg_idx0 = 0, rg_idx1 = 1, count = 0;
        if (regions[0].size() > regions[1].size()) {
          rg_idx0 = 1; rg_idx1 = 0;
        }

        const auto& region0 = regions[0];
        for (std::size_t i = 0; i < region0.size(); ++i) {
          region0[i]->info().probabilities[rg_idx0] += FT(1) * weights[count];
          ++count;
        }

        const auto& region1 = regions[1];
        for (std::size_t i = 0; i < region1.size(); ++i) {
          region1[i]->info().probabilities[rg_idx1] += FT(1) * weights[count];
          ++count;
        }

        return;
      }

      if (num_regions >= 3) {

        std::size_t count = 0;
        const auto& region0 = regions[0];
        for (std::size_t i = 0; i < region0.size(); ++i) {
          region0[i]->info().probabilities[0] += FT(1) * weights[count];
          ++count;
        }

        for (std::size_t k = 1; k < num_regions - 1; ++k) {
          const auto& regionk = regions[k];
          for (std::size_t i = 0; i < regionk.size(); ++i) {
            regionk[i]->info().probabilities[1] += FT(1) * weights[count];
            ++count;
          }
        }

        const auto& region1 = regions[num_regions - 1];
        for (std::size_t i = 0; i < region1.size(); ++i) {
          region1[i]->info().probabilities[0] += FT(1) * weights[count];
          ++count;
        }

        return;
      }
    }

    void compute_angle_weights(
      const Point_2& p1, const Point_2& p2,
      const std::vector< std::vector<F_handle> >& regions,
      std::vector<FT>& weights) {

      std::size_t num_weights = 0;
      for (const auto& region : regions)
        num_weights += region.size();

      const std::size_t num_regions = regions.size();

      if (num_regions == 0)
        return;

      if (num_regions == 1) {
        weights.clear();

        const auto& region = regions[0];
        const std::size_t num_faces = region.size();

        const auto q1 = get_point(region[0]);
        const auto q2 = get_point(region[num_faces - 1]);

        const FT distance = internal::distance(q1, q2);
        if (distance < m_noise_level)
          weights.resize(num_weights, FT(1));
        else
          weights.resize(num_weights, FT(0));

        return;
      }

      FT product = FT(1);
      const Segment_2 segment1 = Segment_2(p1, p2);

      std::vector<FT> angles;
      for (std::size_t i = 0; i < num_regions - 1; ++i) {
        const std::size_t ip = i + 1;

        const auto& region0 = regions[i];
        const auto& region1 = regions[ip];

        const auto& fh  = region0[region0.size() - 1];
        const auto& fhn = region1[0];

        const std::size_t idx = fh->index(fhn);
        const auto& vh1 = fh->vertex( (idx + 1) % 3 );
        const auto& vh2 = fh->vertex( (idx + 2) % 3 );

        const auto& q1 = vh1->point();
        const auto& q2 = vh2->point();

        const Segment_2 segment2 = Segment_2(q1, q2);

        const FT angle_d = angle_degree_2(segment1, segment2);
        const FT angle_2 = CGAL::abs(get_angle_2(angle_d));

        const FT angle = angle_2 * static_cast<FT>(CGAL_PI) / FT(180);
        angles.push_back(angle);
      }

      FT denom = -FT(1);
      for (const auto& angle : angles)
        denom = CGAL::max(denom, angle);

      for (const auto& angle : angles)
        product *= smooth_step(angle, denom);

      weights.clear();
      weights.resize(num_weights, product);
    }

    FT smooth_step(
      const FT num, const FT denom) {

      /*
      const FT x = num / denom;
      CGAL_assertion(x >= FT(0) && x <= FT(1)); */

      /* const FT y = x; */
      /* const FT y = 3 * x * x - 2 * x * x * x; */

      const FT y = static_cast<FT>(std::sin(CGAL::to_double(num)));
      return y;
    }

    FT angle_degree_2(
      const Segment_2& s1, const Segment_2& s2) {

      const Vector_2 v1 =  s1.to_vector();
      const Vector_2 v2 = -s2.to_vector();

		  const FT det = CGAL::determinant(v1, v2);
		  const FT dot = CGAL::scalar_product(v1, v2);
      const FT angle_rad = static_cast<FT>(
        std::atan2(CGAL::to_double(det), CGAL::to_double(dot)));
      const FT angle_deg = angle_rad * FT(180) / m_pi;
      return angle_deg;
    }

    FT get_angle_2(const FT angle) {

      FT angle_2 = angle;
      if (angle_2 > FT(90)) angle_2 = FT(180) - angle_2;
      else if (angle_2 < -FT(90)) angle_2 = FT(180) + angle_2;
      return angle_2;
    }

    void compute_distance_weights(
      const std::vector< std::vector<F_handle> >& regions,
      std::vector<FT>& weights) {

      std::size_t num_weights = 0;
      for (const auto& region : regions)
        num_weights += region.size();

      weights.clear();
      weights.reserve(num_weights);

      std::vector<FT> distances;
      distances.reserve(num_weights);

      const auto q1 = get_point(regions[0][0]);
      for (const auto& region : regions) {
        for (const auto fh : region) {
          const auto q2 = get_point(fh);
          const FT distance = internal::distance(q1, q2);

          if (distance < m_noise_level)
            distances.push_back(-FT(1));
          else
            distances.push_back(distance);
        }
      }

      const std::size_t num_regions = regions.size();
      const std::size_t num_faces   = regions[num_regions - 1].size();

      const auto q3 = get_point(regions[num_regions - 1][num_faces - 1]);
      const FT sum_distance  = internal::distance(q1, q3);

      if (sum_distance != FT(0)) {
        for (const auto& distance : distances) {

          FT weight = FT(1);
          if (distance != -FT(1))
            weight = FT(0);
          weights.push_back(weight);
        }
      } else {
        weights.clear();
        weights.resize(num_weights, FT(1));
      }
    }

    void add_region_alpha_shape_v2(
      const std::vector< std::vector<F_handle> >& regions) {

      const std::size_t num_regions = regions.size();
      if (num_regions == 0) return;

      if (num_regions > 0) {

        const auto& region = regions[0];
        const auto q1 = get_point(region[0]);
        for (auto fh : region) {
          const auto q2 = get_point(fh);

          if (internal::distance(q1, q2) < m_noise_level)
            fh->info().probabilities[0] += FT(1);
          else
            fh->info().probabilities[1] += FT(1);
        }
      }

      if (num_regions > 1) {

        const auto& region = regions[num_regions - 1];
        const std::size_t num_faces = region.size();
        const auto q1 = get_point(region[num_faces - 1]);
        for (auto fh : region) {
          const auto q2 = get_point(fh);

          if (internal::distance(q1, q2) < m_noise_level)
            fh->info().probabilities[0] += FT(1);
          else
            fh->info().probabilities[1] += FT(1);
        }
      }

      if (num_regions > 2) {
        for (std::size_t i = 1; i < num_regions - 1; ++i)
          for (auto fh : regions[i])
            fh->info().probabilities[1] += FT(1);
      }
    }

    void apply_line_walk_from_detected_boundary(
      const Point_2& p,
      const Point_2& q,
      const F_handle ref,
      BaseTri& base) {

      std::vector< std::vector<F_handle> > regions;
      const bool success = add_line_walk(p, q, ref, base, regions);
      if (success)
        add_region_detected_v1(p, q, regions);
    }

    void add_region_detected_v1(
      const Point_2& p1, const Point_2& p2,
      const std::vector< std::vector<F_handle> >& regions) {

      const std::size_t num_regions = regions.size();
      if (num_regions == 0) return;

      std::vector<FT> weights;
      compute_angle_weights(p1, p2, regions, weights);

      if (num_regions == 1) {

        std::size_t count = 0;
        const auto& region = regions[0];
        const std::size_t num_faces = region.size();

        const auto q1 = get_point(region[0]);
        const auto q2 = get_point(region[num_faces - 1]);

        const FT distance = internal::distance(q1, q2);
        if (distance < m_noise_level) {
          for (std::size_t i = 0; i < region.size(); ++i) {
            region[i]->info().probabilities[0] += FT(1) * weights[count];
            ++count;
          }
        } else {
          for (std::size_t i = 0; i < region.size(); ++i) {
            region[i]->info().probabilities[1] += FT(1) * weights[count];
            ++count;
          }
        }

        return;
      }

      if (num_regions >= 2) {

        std::size_t count = 0;
        for (std::size_t k = 0; k < num_regions - 1; ++k) {
          const auto& regionk = regions[k];
          for (std::size_t i = 0; i < regionk.size(); ++i) {
            regionk[i]->info().probabilities[1] += FT(1) * weights[count];
            ++count;
          }
        }

        const auto& region1 = regions[num_regions - 1];
        for (std::size_t i = 0; i < region1.size(); ++i) {
          region1[i]->info().probabilities[0] += FT(1) * weights[count];
          ++count;
        }

        return;
      }
    }

    void add_region_detected_v2(
      const std::vector< std::vector<F_handle> >& regions) {

      const std::size_t num_regions = regions.size();
      if (num_regions == 0) return;

      const auto& region = regions[num_regions - 1];
      const std::size_t num_faces = region.size();
      const auto q1 = get_point(region[num_faces - 1]);
      for (auto fh : region) {
        const auto q2 = get_point(fh);

        if (internal::distance(q1, q2) < m_noise_level)
          fh->info().probabilities[0] += FT(1);
        else
          fh->info().probabilities[1] += FT(1);
      }
    }

    bool add_line_walk(
      const Point_2& p,
      const Point_2& q,
      const F_handle ref,
      const BaseTri& base,
      std::vector< std::vector<F_handle> >& regions) {

      if (base.oriented_side(ref, p) == CGAL::ON_NEGATIVE_SIDE)
        return false;

      LF_circulator circ = base.line_walk(p, q, ref);
      const LF_circulator end = circ;

      regions.clear();
      std::vector<F_handle> region;
      region.push_back(circ);
      do {

        LF_circulator f1 = circ; ++circ;
        LF_circulator f2 = circ;

        if (base.is_infinite(f2) || !f2->info().tagged) {
          if (!region.empty())
            regions.push_back(region);
          break;
        }

        const bool success = are_neighbors(f1, f2);
        if (success) {
          const std::size_t idx = f1->index(f2);
          const auto edge = std::make_pair(f1, idx);

          if (base.is_constrained(edge)) {
            if (!region.empty())
              regions.push_back(region);
            region.clear();
          }
        }
        region.push_back(f2);
      } while (circ != end);

      return true;
    }

    Point_2 get_point(const F_handle fh) {
      return CGAL::barycenter(
        fh->vertex(0)->point(), FT(1),
        fh->vertex(1)->point(), FT(1),
        fh->vertex(2)->point(), FT(1));
    }

    void compute_in_out_stable(
      BaseTri& base) {

      for (auto fh = base.finite_faces_begin();
      fh != base.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;
        compute_statistics(base, fh);
      }
    }

    void compute_statistics(
      const BaseTri& base,
      F_handle fh) {

      const FT radius = FT(1);
      const std::size_t num_samples = 24;
      const Point_2 center = CGAL::barycenter(
        fh->vertex(0)->point(), FT(1),
        fh->vertex(1)->point(), FT(1),
        fh->vertex(2)->point(), FT(1));

      std::vector<Point_2> samples1;
      create_points_on_circle(
        center, radius, FT(0), num_samples, samples1);

      std::vector<Point_2> samples2;
      create_points_on_circle(
        center, radius, FT(180), num_samples, samples2);

      FT inside = FT(0), outside = FT(0);
      for (std::size_t k = 0; k < num_samples / 2; ++k) {
        const auto& p1 = samples1[k];
        const auto& p2 = samples2[k];

        const auto pair = get_in_out_value(
          base, center, p1, p2, fh);

        inside  += pair.first;
        outside += pair.second;
      }

      if (inside == FT(0) && outside == FT(0)) {
        inside = FT(1); outside = FT(1);
      }

      const FT sum = inside + outside;
      inside /= sum; outside /= sum;

      auto& probabilities = fh->info().probabilities;

      probabilities[0] = outside;
      probabilities[1] = inside;
    }

    void create_points_on_circle(
      const Point_2& center,
      const FT radius,
      const FT start,
      const std::size_t num_samples,
      std::vector<Point_2>& samples) {

      samples.clear();
      samples.reserve(num_samples);

      FT factor = FT(360) / static_cast<FT>(num_samples);
      factor *= static_cast<FT>(CGAL_PI); factor /= FT(180);

      FT init = start;
      init *= static_cast<FT>(CGAL_PI); init /= FT(180);

      for (std::size_t i = 0; i < num_samples / 2; ++i) {
        const double angle =
          CGAL::to_double(init) + double(i) * CGAL::to_double(factor);

        const FT cosa = static_cast<FT>(std::cos(angle));
        const FT sina = static_cast<FT>(std::sin(angle));

        const FT x = center.x() + radius * cosa;
        const FT y = center.y() + radius * sina;

        samples.push_back(Point_2(x, y));
      }
    }

    std::pair<FT, FT> get_in_out_value(
      const BaseTri& base,
      const Point_2& p,
      const Point_2& q1,
      const Point_2& q2,
      const F_handle ref) {

      const auto pair1 = apply_line_walk(base, p, q1, ref);
      const auto pair2 = apply_line_walk(base, p, q2, ref);

      const FT inside  = pair1.first  + pair2.first;
      const FT outside = pair1.second + pair2.second;

      return std::make_pair(inside, outside);
    }

    std::pair<FT, FT> apply_line_walk(
      const BaseTri& base,
      const Point_2& p,
      const Point_2& q,
      const F_handle ref) {

      LF_circulator circ = base.line_walk(p, q, ref);
      const LF_circulator end = circ;

      std::size_t inter = 0;
      do {

        LF_circulator f1 = circ; ++circ;
        LF_circulator f2 = circ;

        const bool success = are_neighbors(f1, f2);
        if (!success) break;

        const std::size_t idx = f1->index(f2);
        const auto edge = std::make_pair(f1, idx);
        if (base.is_constrained(edge)) ++inter;
        if (base.is_infinite(f2) || !f2->info().tagged) break;

      } while (circ != end);

      if (inter % 2 == 0) return std::make_pair(FT(0), FT(1));
      else return std::make_pair(FT(1), FT(0));
    }

    bool are_neighbors(
      LF_circulator f1, LF_circulator f2) {

      for (std::size_t i = 0; i < 3; ++i) {
        const std::size_t ip = (i + 1) % 3;

        const auto p1 = f1->vertex(i);
        const auto p2 = f1->vertex(ip);

        for (std::size_t j = 0; j < 3; ++j) {
          const std::size_t jp = (j + 1) % 3;

          const auto q1 = f2->vertex(j);
          const auto q2 = f2->vertex(jp);

          if (
            ( p1 == q1 && p2 == q2) ||
            ( p1 == q2 && p2 == q1) ) {

            return true;
          }
        }
      }
      return false;
    }

    void update_labels(
      BaseTri& base) {

      for (auto fh = base.finite_faces_begin();
      fh != base.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;

        const auto& probabilities = fh->info().probabilities;
        if (probabilities[1] >= FT(1) / FT(2)) // inside
          fh->info().label = 1;
        else
          fh->info().label = 0;
      }
    }

    void set_initial_labels(
      const BaseTri& base,
      std::vector<std::size_t>& labels) {

      labels.clear();
      for (auto fh = base.finite_faces_begin();
      fh != base.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;
        labels.push_back(fh->info().label);
      }
    }

    void set_graphcut_edges(
      const FT beta,
      const BaseTri& base,
      std::vector<Size_pair>& edges,
      std::vector<double>& edge_weights) {

      edges.clear();
      edge_weights.clear();

      FT max_value = -FT(1);
      for (auto eh = base.finite_edges_begin();
      eh != base.finite_edges_end(); ++eh) {
        const auto& edge = *eh;

        const auto  fh = edge.first;
        const auto idx = edge.second;
        const auto fhn = fh->neighbor(idx);

        if (fh->info().tagged && fhn->info().tagged) {

          const std::size_t idxi =  fh->info().object_index;
          const std::size_t idxj = fhn->info().object_index;
          edges.push_back(std::make_pair(idxi, idxj));

          const auto& p1 = fh->vertex((idx + 1) % 3)->point();
          const auto& p2 = fh->vertex((idx + 2) % 3)->point();
          const FT distance = internal::distance(p1, p2);

          FT edge_weight = distance;
          max_value = CGAL::max(max_value, edge_weight);
          if (base.is_constrained(edge))
            edge_weight = -FT(1);
          edge_weights.push_back(CGAL::to_double(edge_weight));
        }
      }

      for (auto& edge_weight : edge_weights) {
        if (edge_weight == -1.0)
          edge_weight = CGAL::to_double(max_value);
        edge_weight /= CGAL::to_double(max_value);
        edge_weight *= CGAL::to_double(beta);
      }
    }

    void set_cost_matrix(
      const BaseTri& base,
      std::vector< std::vector<double> >& cost_matrix) {

      cost_matrix.clear();
      cost_matrix.resize(2);

      FT max_value = -FT(1);
      std::vector<FT> weights;

      for (auto fh = base.finite_faces_begin();
      fh != base.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;

        const auto& p0 = fh->vertex(0)->point();
        const auto& p1 = fh->vertex(1)->point();
        const auto& p2 = fh->vertex(2)->point();

        const Triangle_2 triangle = Triangle_2(p0, p1, p2);
        const FT area = CGAL::abs(triangle.area());

        const FT weight = area;
        max_value = CGAL::max(max_value, weight);
        weights.push_back(weight);
      }

      for (auto& weight : weights)
        weight /= max_value;

      std::size_t count = 0;
      for (auto fh = base.finite_faces_begin();
      fh != base.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;

        for (std::size_t k = 0; k < 2; ++k)
          cost_matrix[k].push_back(get_cost(
            weights[count], fh->info().probabilities[k]));
        ++count;
      }
    }

    double get_cost(
      const FT weight,
      const FT probability) {
      return CGAL::to_double((FT(1) - probability) * weight);
    }

    void set_new_labels(
      const std::vector<std::size_t>& labels,
      BaseTri& base) {

      std::size_t count = 0;
      for (auto fh = base.finite_faces_begin();
      fh != base.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;

        fh->info().label = labels[count];
        ++count;
      }
    }

    void update_tags(
      BaseTri& base) {

      for (auto fh = base.finite_faces_begin();
      fh != base.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;

        bool found = false;
        for (std::size_t k = 0; k < 3; ++k) {
          const auto fhn = fh->neighbor(k);
          if (!fhn->info().tagged || base.is_infinite(fhn)) {
            found = true; break;
          }
        }
        if (!found) continue;

        if (fh->info().label == 0)
          propagate_base(base, fh);
      }
    }

    void propagate_base(
      const BaseTri& base,
      F_handle fh) {

      if (fh->info().label == 0)
        fh->info().tagged = false;
      else return;

      for (std::size_t k = 0; k < 3; ++k) {
        auto fhn = fh->neighbor(k);
        if (!base.is_infinite(fhn) && fhn->info().tagged)
          propagate_base(base, fhn);
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_POINTS_MERGER_2_H
