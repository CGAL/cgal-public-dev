#ifndef CGAL_LEVELS_OF_DETAIL_KINETIC_PARTITIONING_2_H
#define CGAL_LEVELS_OF_DETAIL_KINETIC_PARTITIONING_2_H

// STL includes.
#include <list>
#include <vector>
#include <utility>
#include <unordered_map>

// Kinetic includes.
#include "kinetic_model.h"
#include "propagation.h"

// Internal includes.
#include <CGAL/Levels_of_detail/internal/structures.h>

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

    using Polygon_face_2 = Polygon_face_2<Traits>;

    using Kinetic_face = Face;
    using Kinetic_faces = std::list<Kinetic_face *>;

    Kinetic_partitioning_2(
      const std::size_t max_intersections, 
      const FT min_face_width) :
    m_max_intersections(max_intersections), 
    m_min_face_width(min_face_width), 
    m_bbox_scale(FT(2)),
    m_aspacing_scale(FT(4)),
    m_num_neighbors(6)
    { }

    void compute(
      std::vector<Segment_2>& segments,
      std::vector<Polygon_face_2>& polygon_faces) const {
      
      // Initialize the model.
      Kinetic_Model* model = new Kinetic_Model();
      model->reinit();

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

      // Compute partitioning.
      set_kinetic_segments(segments, model);
      compute_partitioning(bbox_size, scale, model);

      // Get back all polygon faces.
      translation = Point_2(-translation.x(), -translation.y());
      create_output(translation, scale, model, polygon_faces);

      // Clear.
      model->segments.clear();
      delete model;
    }

  private:

    // External parameters.
    const std::size_t m_max_intersections;
    const FT m_min_face_width;

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

      // Positions segments at the center of the bounding box.
      translation = 
      Point_2(minx - lengthx / FT(2), miny - lengthy / FT(2));

      const FT bbox_width = lengthx * m_bbox_scale;
      const FT bbox_height = lengthy * m_bbox_scale;

      bbox_size = 
      std::make_pair(bbox_width, bbox_height);
    }

    void translate_segments(
      const Point_2& translation,
      std::vector<Segment_2>& segments) const {

      for (std::size_t i = 0; i < segments.size(); ++i) {
        Segment_2& segment = segments[i];
                      
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

      CGAL_precondition(
        x > FT(0) && 
        y > FT(0) && 
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

      for (std::size_t i = 0; i < segments.size(); ++i) {
        Segment_2& segment = segments[i];
                      
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
      Kinetic_Model* model) const {
                  
      CGAL_precondition(segments.size() > 0);
      auto& model_segments = model->segments;

      model_segments.reserve(segments.size());
      for (std::size_t i = 0; i < segments.size(); ++i) {
        
        const Segment_2& segment = segments[i];
        const double width = 
        CGAL::sqrt(
          CGAL::to_double(
            segment.squared_length()));

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

    void compute_partitioning(
      const std::pair<FT, FT>& bbox_size, 
      const std::pair<FT, FT>& scale, 
      Kinetic_Model* model) const {

      Propagation propagation;

      const std::size_t rows = 
      static_cast<std::size_t>(
        std::ceil(
          CGAL::to_double(
            scale.second * bbox_size.second)));

      const std::size_t cols = 
      static_cast<std::size_t>(
        std::ceil(
          CGAL::to_double(
            scale.first * bbox_size.first)));

      propagation.dmitry_size_rows = rows;
      propagation.dmitry_size_cols = cols;

      model->set_prop_ttl(static_cast<int>(m_max_intersections));
      model->set_prop_merge_min_thinness(CGAL::to_double(m_min_face_width));

      propagation.propagate(model);
    }

    void create_output(
      const Point_2& translation, 
      const std::pair<FT, FT>& scale, 
      const Kinetic_Model* model, 
      std::vector<Polygon_face_2>& polygon_faces) const {

      const auto* graph = model->graph;
      const auto& kinetic_faces = graph->faces;
      create_polygon_faces(translation, scale, kinetic_faces, polygon_faces);
    }

    void create_polygon_faces(
      const Point_2& translation, 
      const std::pair<FT, FT>& scale, 
      const Kinetic_faces& kinetic_faces, 
      std::vector<Polygon_face_2>& polygon_faces) const {

      CGAL_precondition(kinetic_faces.size() > 0);
      
      polygon_faces.clear();
      polygon_faces.resize(kinetic_faces.size());

      std::size_t i = 0;
      std::unordered_map<std::size_t, std::size_t> fmap;

      for (auto fit = kinetic_faces.begin(); 
      fit != kinetic_faces.end(); ++fit, ++i) {
        const auto& kinetic_face = **fit;

        const int id_face = kinetic_face.id_face;
        if (id_face < 0) 
          std::cerr << 
            "Error (create_polygon_faces): Negative face id!" 
          << std::endl;

        const std::size_t face_id = static_cast<std::size_t>(id_face);
        fmap[face_id] = i;

        Polygon_face_2 polygon_face;
        create_face_vertices(translation, scale, kinetic_face, polygon_face);
        polygon_faces[i] = polygon_face;
      }

      i = 0;
      for (auto fit = kinetic_faces.begin(); 
      fit != kinetic_faces.end(); ++fit, ++i) {
        const auto& kinetic_face = **fit;
        
        auto& polygon_face = polygon_faces[i];
        create_face_neighbors(fmap, kinetic_face, polygon_face);
      }
    }

    void create_face_vertices(
      const Point_2& translation,
      const std::pair<FT, FT>& scale,
      const Kinetic_face& kinetic_face,
      Polygon_face_2& polygon_face) const {

      const FT scale_x = scale.first;
      const FT scale_y = scale.second;

      const auto& kinetic_vertices = kinetic_face.vertices;
      auto& face_vertices = polygon_face.vertices;

      CGAL_precondition(kinetic_vertices.size() >= 3);
      face_vertices.reserve(kinetic_vertices.size());

      for (auto vit = kinetic_vertices.begin(); 
      vit != kinetic_vertices.end(); ++vit) {  
        
        const auto& point = vit->first->pt;

        const FT x = FT(point.x) / scale_x - translation.x();
        const FT y = FT(point.y) / scale_y - translation.y();

        face_vertices.push_back(Point_2(x, y));
      }
    }

    void create_face_neighbors(
      const std::unordered_map<std::size_t, std::size_t>& fmap,
      const Kinetic_face& kinetic_face,
      Polygon_face_2& polygon_face) const {

      const auto& kinetic_edges = kinetic_face.edges;
      auto& face_neighbors = polygon_face.neighbors;

      for (auto eit = kinetic_edges.begin();
      eit != kinetic_edges.end(); ++eit) {

        const auto& opposite = *((*eit)->opposite());
        if (opposite.f == nullptr) 
          continue;

        const auto& kinetic_face = *(opposite.f);
        const int id_face = kinetic_face.id_face;

        if (id_face < 0) 
          std::cerr << 
            "Error (create_face_neighbors): Negative face id!" 
          << std::endl;

        const std::size_t face_id = 
        static_cast<std::size_t>(id_face);

        CGAL_precondition(fmap.find(face_id) != fmap.end());
        face_neighbors.push_back(fmap.at(face_id));
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_KINETIC_PARTITIONING_2_H
