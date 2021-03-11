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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_KINETIC_PARTITIONING_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_KINETIC_PARTITIONING_2_H

// STL includes.
#include <list>
#include <vector>
#include <utility>
#include <unordered_map>

// CGAL includes.
#include <CGAL/Barycentric_coordinates_2/Segment_coordinates_2.h>

// Kinetic includes.
#include "kinetic2/kinetic_model.h"
#include "kinetic2/propagation.h"

// Internal includes.
#include <CGAL/Levels_of_detail/internal/struct.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<class GeomTraits>
  class Kinetic_partitioning_2 {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;
    using Plane_3 = typename Traits::Plane_3;

    using Partition_2 = internal::Partition_2<Traits>;
    using Partition_edge_2 = typename Partition_2::Edge;
    using Partition_face_2 = typename Partition_2::Face;

    using Kinetic_face = Face;
    using Kinetic_edge = Edge;
    using Kinetic_hedge = HalfEdge;
    using Kinetic_segments = std::vector<Segment*>;
    using Kinetic_faces = std::list<Kinetic_face*>;

    struct Kinetic_model_2 {

      Kinetic_Model* model;
      Kinetic_model_2() {
        model = new Kinetic_Model();
        model->reinit();
      }
      ~Kinetic_model_2() {
        segments().clear();
        delete model;
      }
      Kinetic_segments& segments() {
        return model->segments;
      }
      const Kinetic_faces& faces() const {
        return model->graph->faces;
      }
      void set_min_face_width(const FT min_face_width) {
        model->set_prop_merge_min_thinness(CGAL::to_double(min_face_width));
      }
      void set_max_intersections(const std::size_t max_intersections) {
        model->set_prop_ttl(static_cast<int>(max_intersections));
      }
      Kinetic_Model* get() {
        return model;
      }
      Kinetic_edge* edges_begin() {
        return model->graph->edges_head;
      }
    };

    Kinetic_partitioning_2(
      const FT kinetic_min_face_width_2,
      const std::size_t kinetic_max_intersections_2) :
    m_min_face_width(kinetic_min_face_width_2),
    m_max_intersections(kinetic_max_intersections_2),
    m_bbox_scale(FT(2)),
    m_aspacing_scale(FT(4)),
    m_num_neighbors(6)
    { }

    void compute(
      std::vector<Segment_2>& segments,
      Partition_2& partition) const {

      // Initialize the model.
      Kinetic_model_2 model;

      // We first compute translation factor and size of the bounding box and
      // then translate all segments such that the bottom left corner of the
      // bounding box becomes (0, 0);
      Point_2 translation;
      std::pair<FT, FT> bbox_size;
      compute_translation_and_bbox_size(segments, translation, bbox_size);
      translate_segments(translation, segments);

      // Scale segments.
      const FT average_spacing =
      internal::average_spacing_2(segments, m_num_neighbors) / m_aspacing_scale;
      std::pair<FT, FT> scale;
      get_scale(bbox_size, average_spacing, scale);
      scale_segments(scale, segments);

      // Compute partition.
      set_kinetic_segments(segments, model);
      compute_partition(bbox_size, scale, model);

      // Get back all polygon faces.
      translation = Point_2(-translation.x(), -translation.y());
      create_partition(translation, scale, model, partition);
    }

  private:

    // External parameters.
    const FT m_min_face_width;
    const std::size_t m_max_intersections;

    // Internal parameters.
    const FT m_bbox_scale;
    const FT m_aspacing_scale;
    const std::size_t m_num_neighbors;

    void compute_translation_and_bbox_size(
      const std::vector<Segment_2>& segments,
      Point_2& translation,
      std::pair<FT, FT>& bbox_size) const {

      std::vector<Point_2> bbox;
      internal::bounding_box_2(segments, bbox);

      const FT minx = bbox[0].x(); const FT miny = bbox[0].y();
      const FT maxx = bbox[2].x(); const FT maxy = bbox[2].y();

      const FT lengthx = CGAL::abs(maxx - minx);
      const FT lengthy = CGAL::abs(maxy - miny);

      // Position segments at the center of the bounding box.
      translation =
      Point_2(minx - lengthx / FT(2), miny - lengthy / FT(2));

      const FT bbox_width = lengthx * m_bbox_scale;
      const FT bbox_height = lengthy * m_bbox_scale;
      bbox_size = std::make_pair(bbox_width, bbox_height);
    }

    void translate_segments(
      const Point_2& translation,
      std::vector<Segment_2>& segments) const {

      for (auto& segment : segments) {
        const Point_2& source = segment.source();
        const Point_2& target = segment.target();

        const FT x1 = source.x() - translation.x();
        const FT y1 = source.y() - translation.y();

        const FT x2 = target.x() - translation.x();
        const FT y2 = target.y() - translation.y();

        Point_2 new_source = Point_2(x1, y1);
        Point_2 new_target = Point_2(x2, y2);
        segment = Segment_2(new_source, new_target);
      }
    }

    void get_scale(
      const std::pair<FT, FT>& bbox_size,
      const FT average_spacing,
      std::pair<FT, FT>& scale) const {

      const FT x = bbox_size.first;
      const FT y = bbox_size.second;

      CGAL_assertion(
        x > FT(0) && y > FT(0) &&
        average_spacing > FT(0));

      const FT x_num = x / average_spacing;
      const FT y_num = y / average_spacing;
      scale = std::make_pair(x_num, y_num);
    }

    void scale_segments(
      const std::pair<FT, FT>& scale,
      std::vector<Segment_2>& segments) const {

      const FT scale_x = scale.first;
      const FT scale_y = scale.second;

      for (auto& segment : segments) {
        const Point_2& source = segment.source();
        const Point_2& target = segment.target();

        const FT x1 = scale_x * source.x();
        const FT y1 = scale_y * source.y();

        const FT x2 = scale_x * target.x();
        const FT y2 = scale_y * target.y();

        Point_2 new_source = Point_2(x1, y1);
        Point_2 new_target = Point_2(x2, y2);
        segment = Segment_2(new_source, new_target);
      }
    }

    void set_kinetic_segments(
      const std::vector<Segment_2>& segments,
      Kinetic_model_2& model) const {

      CGAL_assertion(segments.size() > 0);
      auto& model_segments = model.segments();

      model_segments.clear();
      model_segments.reserve(segments.size());
      for (std::size_t i = 0; i < segments.size(); ++i) {

        const Segment_2& segment = segments[i];
        const double width =
        CGAL::sqrt(CGAL::to_double(segment.squared_length()));

        const Point_2& source = segment.source();
        const Point_2& target = segment.target();

        const double x1 = CGAL::to_double(source.x());
        const double y1 = CGAL::to_double(source.y());

        const double x2 = CGAL::to_double(target.x());
        const double y2 = CGAL::to_double(target.y());

        model_segments.push_back(
          new Segment(i, x1, y1, x2, y2, width, 0.0, 0.0, 0.0, false));
      }
    }

    void compute_partition(
      const std::pair<FT, FT>& bbox_size,
      const std::pair<FT, FT>& scale,
      Kinetic_model_2& model) const {

      Propagation propagation;

      const std::size_t rows =
      static_cast<std::size_t>(
        std::ceil(CGAL::to_double(scale.second * bbox_size.second)));
      const std::size_t cols =
      static_cast<std::size_t>(
        std::ceil(CGAL::to_double(scale.first * bbox_size.first)));

      propagation.dmitry_size_rows = rows;
      propagation.dmitry_size_cols = cols;

      model.set_min_face_width(m_min_face_width);
      model.set_max_intersections(m_max_intersections);
      propagation.propagate(model.get());
    }

    void create_partition(
      const Point_2& translation,
      const std::pair<FT, FT>& scale,
      Kinetic_model_2& model,
      Partition_2& partition) const {

      partition.clear();
      std::unordered_map<int, int> fmap;
      create_partition_faces(
        translation, scale, model.faces(), fmap, partition);
      create_partition_neighbors(
        translation, scale, model.faces(), fmap, partition);
      create_partition_edges(
        translation, scale, fmap, model.edges_begin(), partition);
    }

    void create_partition_faces(
      const Point_2& translation,
      const std::pair<FT, FT>& scale,
      const Kinetic_faces& faces,
      std::unordered_map<int, int>& fmap,
      Partition_2& partition) const {

      CGAL_assertion(faces.size() > 0);
      partition.faces.reserve(faces.size());

      int id = 0; std::vector<Point_2> polygon;
      for (const auto& fh : faces) {
        const auto& face = *fh;

        const int face_id = face.id_face;
        CGAL_assertion(face_id >= 0 && face_id < faces.size());
        fmap[face_id] = id; ++id;

        create_polygon(translation, scale, face, polygon);
        partition.faces.push_back(Partition_face_2(polygon));
      }
      CGAL_assertion(partition.faces.size() == faces.size());
    }

    void create_polygon(
      const Point_2& translation,
      const std::pair<FT, FT>& scale,
      const Kinetic_face& face,
      std::vector<Point_2>& polygon) const {

      polygon.clear();
      const FT scale_x = scale.first;
      const FT scale_y = scale.second;

      const auto& vertices = face.vertices;
      CGAL_assertion(vertices.size() >= 3);

      for (const auto& vertex : vertices) {
        const auto& p = vertex.first->pt;
        const FT x = FT(p.x) / scale_x - translation.x();
        const FT y = FT(p.y) / scale_y - translation.y();
        polygon.push_back(Point_2(x, y));
      }
      CGAL_assertion(polygon.size() == vertices.size());
    }

    void create_partition_neighbors(
      const Point_2& translation,
      const std::pair<FT, FT>& scale,
      const Kinetic_faces& kinetic_faces,
      const std::unordered_map<int, int>& fmap,
      Partition_2& partition) const {

      CGAL_assertion(partition.faces.size() == kinetic_faces.size());
      CGAL_assertion(fmap.size() > 0);

      std::size_t i = 0;
      for (auto fit = kinetic_faces.begin();
      fit != kinetic_faces.end(); ++fit, ++i) {
        const auto& kinetic_face = **fit;
        auto& partition_face = partition.faces[i];
        create_partition_face_neighbors(
          translation, scale, fmap, kinetic_face, partition_face);
      }
      CGAL_assertion(i == kinetic_faces.size());
    }

    void create_partition_face_neighbors(
      const Point_2& translation,
      const std::pair<FT, FT>& scale,
      const std::unordered_map<int, int>& fmap,
      const Kinetic_face& kinetic_face,
      Partition_face_2& partition_face) const {

      const auto& edges = kinetic_face.edges;

      auto& neighbors = partition_face.neighbors;
      auto& nedges = partition_face.edges;
      auto& nconstr = partition_face.constraints;
      neighbors.clear();
      neighbors.reserve(edges.size());
      nedges.clear();
      nedges.reserve(edges.size());
      nconstr.clear();

      Segment_2 segment;
      for (auto& eh : edges) {
        auto& hedge = *eh;

        const auto& opposite = *(hedge.opposite());
        if (opposite.f == nullptr) {
          neighbors.push_back(-1);
          get_segment(translation, scale, hedge, segment);
          nedges.push_back(segment);
          nconstr[neighbors.back()] = true;
          continue;
        }

        const auto& face = *(opposite.f);
        const int face_id = face.id_face;
        CGAL_assertion(face_id >= 0);
        CGAL_assertion(fmap.find(face_id) != fmap.end());
        const int fidx = fmap.at(face_id);
        neighbors.push_back(fidx);
        get_segment(translation, scale, hedge, segment);
        nedges.push_back(segment);
        nconstr[fidx] = get_constraint(hedge);
      }
      CGAL_assertion(neighbors.size() == edges.size());
      CGAL_assertion(nedges.size() == edges.size());
      CGAL_assertion(nconstr.size() > 0);
    }

    void get_segment(
      const Point_2& translation,
      const std::pair<FT, FT>& scale,
      const Kinetic_hedge& hedge,
      Segment_2& segment) const {

      const FT scale_x = scale.first;
      const FT scale_y = scale.second;

      const auto& p1 = (hedge.e)->v1->pt;
      const auto& p2 = (hedge.e)->v2->pt;

      const FT x1 = FT(p1.x) / scale_x - translation.x();
      const FT y1 = FT(p1.y) / scale_y - translation.y();
      const FT x2 = FT(p2.x) / scale_x - translation.x();
      const FT y2 = FT(p2.y) / scale_y - translation.y();

      segment = Segment_2(Point_2(x1, y1), Point_2(x2, y2));
    }

    bool get_constraint(
      const Kinetic_hedge& hedge) const {

      const auto& edge = hedge.e;
      const auto& p1 = edge->v1->pt;
      const auto& p2 = edge->v2->pt;

      if (edge->type == INNER_EDGE) {
        auto inner = static_cast<Inner_Edge *const>(edge);

        const auto& rays = inner->rays;
        for (const auto& ray : rays) {
          const auto& seg = ray->parent;

          const auto& q1 = seg->end1;
          const auto& q2 = seg->end2;
          return overlap(p1, p2, q1, q2);
        }
      }
      return false;
    }

    bool overlap(
      const Point2d& a, const Point2d& b,
      const Point2d& c, const Point2d& d) const {

      const Point_2 p1 = Point_2(a.x, a.y);
      const Point_2 p2 = Point_2(b.x, b.y);
      const Point_2 q1 = Point_2(c.x, c.y);
      const Point_2 q2 = Point_2(d.x, d.y);

      if (is_inside(p1, p2, q1, q2))
        return true;
      if (is_inside(q1, q2, p1, p2))
        return true;
      return false;
    }

    bool is_inside(
      const Point_2& p1, const Point_2& p2,
      const Point_2& q1, const Point_2& q2) const {

      const Traits traits;
      const auto res1 =
      Barycentric_coordinates::compute_segment_coordinates_2(p1, p2, q1, traits);
      const auto res2 =
      Barycentric_coordinates::compute_segment_coordinates_2(p1, p2, q2, traits);

      const FT bval = -FT(1) / FT(10);
      const FT tval = FT(11) / FT(10);

      if (res1[0] > bval && res1[1] > bval &&
          res1[0] < tval && res1[1] < tval) return true;

      if (res2[0] > bval && res2[1] > bval &&
          res2[0] < tval && res2[1] < tval) return true;

      // if (res1[0] > bval && res1[1] < tval) return true;
      // if (res2[0] > bval && res2[1] < tval) return true;

      return false;
    }

    void create_partition_edges(
      const Point_2& translation,
      const std::pair<FT, FT>& scale,
      const std::unordered_map<int, int>& fmap,
      Kinetic_edge* edge,
      Partition_2& partition) const {

      const FT scale_x = scale.first;
      const FT scale_y = scale.second;

      CGAL_assertion(fmap.size() > 0);
      while (edge != NULL) {
        const auto& p = edge->v1->pt;
        const auto& q = edge->v2->pt;

        const FT x1 = FT(p.x) / scale_x - translation.x();
        const FT y1 = FT(p.y) / scale_y - translation.y();
        const FT x2 = FT(q.x) / scale_x - translation.x();
        const FT y2 = FT(q.y) / scale_y - translation.y();

        const auto& hei = *(edge->v1_v2);
        const auto& hej = *(edge->v2_v1);

        int face_id_i, face_id_j, fidxi, fidxj;

        if (hei.f == nullptr) {
          fidxi = -1;
        } else {
          face_id_i = (hei.f)->id_face;
          CGAL_assertion(face_id_i >= 0);
          CGAL_assertion(fmap.find(face_id_i) != fmap.end());
          fidxi = fmap.at(face_id_i);
        }

        if (hej.f == nullptr) {
          fidxj = -1;
        } else {
          face_id_j = (hej.f)->id_face;
          CGAL_assertion(face_id_j >= 0);
          CGAL_assertion(fmap.find(face_id_j) != fmap.end());
          fidxj = fmap.at(face_id_j);
        }

        partition.edges.push_back(Partition_edge_2(
          Point_2(x1, y1), Point_2(x2, y2), fidxi, fidxj));
        edge = edge->e_next;
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_KINETIC_PARTITIONING_2_H
