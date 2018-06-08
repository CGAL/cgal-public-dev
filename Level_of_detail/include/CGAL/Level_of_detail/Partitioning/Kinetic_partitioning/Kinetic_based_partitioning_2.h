#ifndef CGAL_LEVEL_OF_DETAIL_KINETIC_BASED_PARTITIONING_H
#define CGAL_LEVEL_OF_DETAIL_KINETIC_BASED_PARTITIONING_H

// STL includes.
#include <map>
#include <list>
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Level_of_detail/Tools/Transformations/Translator.h>
#include <CGAL/Level_of_detail/Tools/Estimations/Bounding_box_estimator.h>
#include <CGAL/Level_of_detail/Tools/Estimations/Average_spacing_estimator.h>

// Kinetic includes.
#include "kinetic_model.h"
#include "propagation.h"

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

        template<class InputKernel, class PartitionFace>
        class Kinetic_based_partitioning_2 {

        public:
            using Kernel         = InputKernel;
            using Partition_face = PartitionFace;

            using FT        = typename Kernel::FT;
            using Point_2   = typename Kernel::Point_2;
            using Segment_2 = typename Kernel::Segment_2;

            using Bounding_box           = std::vector<Point_2>;
            using Bounding_box_estimator = LOD::Bounding_box_estimator<Kernel>;
            
            using Bounding_box_size = std::pair<FT, FT>;
            using Scale             = std::pair<FT, FT>;

            using Internal_point_map   = CGAL::Identity_property_map<Point_2>;
            using Internal_segment_map = CGAL::Identity_property_map<Segment_2>;

            using Segments             = std::vector<Segment_2>;
            using Model_segments       = std::vector<Segment *>; // Segment is the class from the kinetic framework

            using Translator                = LOD::Translator<Kernel>;
            using Average_spacing_estimator = LOD::Average_spacing_estimator<Kernel>;
            
            using Face_vertices = std::vector<Point_2>;

            // Kinetic related types.
            using Kinetic_face     = Face;
            using Kinetic_faces    = std::list<Kinetic_face *>;
            using Kinetic_vertex   = std::pair<Vertex *, Vertex_Type>;
            using Kinetic_vertices = std::list<Kinetic_vertex>;
            using Kinetic_point    = Point2d;

            using Const_kinetic_faces_iterator    = typename Kinetic_faces::const_iterator;
            using Const_kinetic_vertices_iterator = typename Kinetic_vertices::const_iterator;

            Kinetic_based_partitioning_2(const size_t num_intersections, const FT min_face_width) :
            m_num_intersections(num_intersections), 
            m_min_face_width(min_face_width), 
            m_bbox_scale(FT(2)),
            m_aspacing_scale(FT(4)),
            m_num_neighbours(6)
            { }

            template<class Elements, class Segment_map, class Output>
            void compute(const Elements &elements, const Segment_map &segment_map, Output &output) const {
            
                // Initialize the model.
                Kinetic_Model* model = new Kinetic_Model(); // Kinetic_Model is the class from the kinetic framework
                model->reinit();

                // Create segments.
                Segments segments;
                create_segments(elements, segment_map, segments);

                // Compute translation factor and size of the data bounding box.
                Point_2           translation; 
                Bounding_box_size bbox_size;
                compute_translation_and_bbox_size(segments, translation, bbox_size);

                // Translate all segments such that the bottom left corner becomes (0, 0);
                const Translator translator;
                translator.translate_segments_2(translation, m_internal_segment_map, segments);

                // Compute average spacing.
                const Average_spacing_estimator average_spacing_estimator(m_num_neighbours);
                const FT average_spacing = average_spacing_estimator.compute_on_segments_2(segments, m_internal_segment_map) / m_aspacing_scale;
                
                // Compute the proper scale.
                Scale scale;
                get_scale(bbox_size, average_spacing, scale);

                // Set input segments to the kinetic model.
                set_kinetic_model_segments(segments, scale, model);

                // Compute partition.
                compute_partition(bbox_size, scale, model);

                // Get reverse translation factor.
                translation = Point_2(-translation.x(), -translation.y());

                // Get back all partition faces.
                create_output(translation, scale, model, output);

                // Clear.
                model->segments.clear();
                delete model;
            }

        private:
            // Important.
            const size_t m_num_intersections;
            const FT     m_min_face_width;

            // Internal parameters.
            const FT     m_bbox_scale;
            const FT     m_aspacing_scale;
            const size_t m_num_neighbours;

            const Internal_point_map   m_internal_point_map;
            const Internal_segment_map m_internal_segment_map;

            template<class Elements, class Segment_map>
            void create_segments(const Elements &elements, const Segment_map &segment_map, Segments &segments) const {
                using Const_elements_iterator = typename Elements::const_iterator;

                segments.clear();
                segments.resize(elements.size());

                size_t i = 0;
                for (Const_elements_iterator ce_it = elements.begin(); ce_it != elements.end(); ++ce_it, ++i)
                    segments[i] = get(segment_map, *ce_it);
            }

            void compute_translation_and_bbox_size(const Segments &segments, Point_2 &translation, Bounding_box_size &bbox_size) const {
                
                Bounding_box bbox;
                const Bounding_box_estimator bounding_box_estimator;
                bounding_box_estimator.compute_bounding_box_2(segments, m_internal_segment_map, bbox);

                const FT minx = bbox[0].x(); const FT miny = bbox[0].y();
                const FT maxx = bbox[2].x(); const FT maxy = bbox[2].y();

                const FT lengthx = CGAL::abs(maxx - minx);
                const FT lengthy = CGAL::abs(maxy - miny);

                const FT bbox_width  = lengthx * m_bbox_scale;
                const FT bbox_height = lengthy * m_bbox_scale;
                    
                const FT scale = FT(2) * m_bbox_scale;

                const FT half_box_x = bbox_width  / scale;
                const FT half_box_y = bbox_height / scale;

                translation = Point_2(minx - half_box_x, miny - half_box_y);
                bbox_size = std::make_pair(bbox_width, bbox_height);
            }

            void get_scale(const Bounding_box_size &bbox_size, const FT average_spacing, Scale &scale) const {

                const FT x = bbox_size.first;
                const FT y = bbox_size.second;

                CGAL_precondition(x > FT(0) && y > FT(0) && average_spacing > FT(0));

                const FT x_num = x / average_spacing;
                const FT y_num = y / average_spacing;

                scale = std::make_pair(x_num, y_num);
            }

            void set_kinetic_model_segments(const Segments &segments, const Scale &scale, Kinetic_Model *model) const {
                
                CGAL_precondition(segments.size() > 0);
                Model_segments &model_segments = model->segments;
                
                const double scale_x = CGAL::to_double(scale.first);
                const double scale_y = CGAL::to_double(scale.second);

                model_segments.reserve(segments.size());
                for (size_t i = 0; i < segments.size(); ++i) {

                    const Point_2 &source = segments[i].source();
                    const Point_2 &target = segments[i].target();

                    const double x1 = scale_x * CGAL::to_double(source.x());
                    const double y1 = scale_y * CGAL::to_double(source.y());

                    const double x2 = scale_x * CGAL::to_double(target.x());
                    const double y2 = scale_y * CGAL::to_double(target.y());

                    const Segment_2 scaled_segment = Segment_2(Point_2(FT(x1), FT(y1)), Point_2(FT(x2), FT(y2)));

                    const double width = CGAL::sqrt(CGAL::to_double(scaled_segment.squared_length()));
                    model_segments.push_back(new Segment(i, x1, y1, x2, y2, width, 0.0, 0.0, 0.0, false));
                }
            }

            void compute_partition(const Bounding_box_size &bbox_size, const Scale &scale, Kinetic_Model *model) const {

                Propagation propagation; // Propagation is the class from the kinetic framework

                const size_t rows = static_cast<size_t>(std::ceil(CGAL::to_double(scale.second * bbox_size.second)));
                const size_t cols = static_cast<size_t>(std::ceil(CGAL::to_double(scale.first  * bbox_size.first)));

                propagation.dmitry_size_rows = rows;
                propagation.dmitry_size_cols = cols;

                model->set_prop_ttl(static_cast<int>(m_num_intersections));
                model->set_prop_merge_min_thinness(CGAL::to_double(m_min_face_width));

                propagation.propagate(model);
            }

            template<class Output>
            void create_output(const Point_2 &translation, const Scale &scale, const Kinetic_Model *model, Output &output) const {

                const Partition *graph = model->graph; // Partition is the class from the kinetic framework
                const Kinetic_faces &kinetic_faces = graph->faces;
                create_partition_faces(translation, scale, kinetic_faces, output);
            }

            template<class Output>
            void create_partition_faces(const Point_2 &translation, const Scale &scale, const Kinetic_faces &kinetic_faces, Output &output) const {

                output.clear();
                CGAL_precondition(kinetic_faces.size() > 0);
                
                for (Const_kinetic_faces_iterator kf_it = kinetic_faces.begin(); kf_it != kinetic_faces.end(); ++kf_it)
                    create_partition_face(translation, scale, *kf_it, output);
            }

            template<class Output>
            void create_partition_face(const Point_2 &translation, const Scale &scale, const Kinetic_face *kinetic_face, Output &output) const {
                
                const Kinetic_vertices &kinetic_face_vertices = kinetic_face->vertices;
                CGAL_precondition(kinetic_face_vertices.size() > 2);
                Face_vertices face_vertices(kinetic_face_vertices.size());

                const FT scale_x = scale.first;
                const FT scale_y = scale.second;

                size_t i = 0;
                for (Const_kinetic_vertices_iterator kv_it = kinetic_face_vertices.begin(); kv_it != kinetic_face_vertices.end(); ++kv_it, ++i) {
                    const Kinetic_point &point = kv_it->first->pt;

                    const FT x = FT(point.x) / scale_x - translation.x();
                    const FT y = FT(point.y) / scale_y - translation.y();

                    face_vertices[i] = Point_2(x, y);
                }
                output.push_back(Partition_face(face_vertices, m_internal_point_map));
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_KINETIC_BASED_PARTITIONING_H