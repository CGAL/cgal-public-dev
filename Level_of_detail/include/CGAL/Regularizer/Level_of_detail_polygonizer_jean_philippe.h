#ifndef CGAL_LEVEL_OF_DETAIL_POLYGONIZER_JEAN_PHILIPPE_H
#define CGAL_LEVEL_OF_DETAIL_POLYGONIZER_JEAN_PHILIPPE_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\"
#else 
#define PSR "/" 
#endif 

// STL includes.
#include <map>
#include <list>
#include <cmath>
#include <vector>
#include <utility>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/compute_average_spacing.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>

// Jean Philippe includes.
#include <CGAL/Regularizer/jean_philippe/kinetic_model.h>
#include "CGAL/Regularizer/jean_philippe/propagation.h"

namespace CGAL {

	namespace LOD {

        template<class KernelTraits, class DataStructure>
        class Level_of_detail_polygonizer_jean_philippe { 

        public:
            typedef KernelTraits  Kernel;
            typedef DataStructure Data_structure;
            
            typedef typename Kernel::FT        FT;
            typedef typename Kernel::Point_2   Point_2;
            typedef typename Kernel::Point_3   Point_3;
            typedef typename Kernel::Segment_2 Segment_2;

            using Segments       = std::vector<Segment_2>;
            using Model_segments = std::vector<Segment *>;

            using Bbox_size = std::pair<FT, FT>;
            using Scale     = std::pair<FT, FT>;

            using Log = CGAL::LOD::Mylog;

            using JP_Face     = Face;
            using JP_Faces    = std::list<JP_Face *>;
            using JP_Vertex   = std::pair<Vertex *, Vertex_Type>;
            using JP_Vertices = std::list<JP_Vertex>;
            using JP_Point    = Point2d;

            using Container  = typename Data_structure::Container;
            using Containers = typename Data_structure::Containers;

            using Polygon = typename Container::Polygon;
            using Points  = std::vector<Point_2>;

            Level_of_detail_polygonizer_jean_philippe() :
            m_silent(false), m_debug(false), m_num_intersections(2), m_min_face_width(FT(1)), 
            m_num_neighbours(6), m_local_scaling(FT(4)) { }

            void polygonize(Segments &segments, Data_structure &data_structure) const {

                if (m_debug) {
                    const std::string stub = "";
                    Log segments_exporter; segments_exporter.export_segments_as_obj("tmp" + std::string(PSR) + "polygonizer_input_segments_jean_philippe", segments, stub);
                }

                Kinetic_Model* model = new Kinetic_Model();
                
                // Initialize the model.
                initialize_kinetic_model(model);

                // Find bottom left corner of the segments' bounding box.
                Point_2 bl; Bbox_size bbox_size;
                find_bottom_left_corner(segments, bl, bbox_size);

                // Translate all segments so that bottom left corner becomes (0, 0);
                translate_segments(bl, segments, true);

                // Compute average spacing and scale.
                Scale scale;

                const FT average_spacing = get_average_spacing(segments);
                get_scale(bbox_size, average_spacing, scale);

                // Set input segments to the model.
                set_input_segments(segments, scale, model);
                
                // Compute partition.
                compute_partition(model, bbox_size, scale);

                // Get back segments.
                get_back_segments(segments, model, scale);

                // Get reverse translation factor.
                bl = Point_2(-bl.x(), -bl.y());

                // Translate back all segments.
                translate_segments(bl, segments, false);

                // Save results.
                if (!m_silent) save_partition(segments);

                // Built data structure.
                built_data_structure(bl, scale, model, data_structure);

                delete model;
            }

            void make_silent(const bool new_state) {
                m_silent = new_state;
            }

            void set_number_of_intersections(const size_t new_value) {
                assert(new_value > 0);
                m_num_intersections = new_value;
            }

            void set_min_face_width(const FT new_value) {
                assert(new_value > FT(0));
                m_min_face_width = new_value;
            }

        private:
            bool m_silent;
            bool m_debug;

            size_t m_num_intersections;
            FT     m_min_face_width;

            const size_t m_num_neighbours;
            const FT     m_local_scaling;

            void initialize_kinetic_model(Kinetic_Model *model) const {
                model->reinit();
            }

            void find_bottom_left_corner(const Segments &segments, Point_2 &bl, Bbox_size &bbox_size) const {
                assert(segments.size() > 0);

                FT minx = FT(1000000000000);
                FT miny = FT(1000000000000);

                FT maxx = -FT(1000000000000);
                FT maxy = -FT(1000000000000);

                for (size_t i = 0; i < segments.size(); ++i) {

                    const Point_2 &source = segments[i].source();
                    const Point_2 &target = segments[i].target();

                    // Min.
                    minx = CGAL::min(minx, source.x());
                    minx = CGAL::min(minx, target.x());

                    miny = CGAL::min(miny, source.y());
                    miny = CGAL::min(miny, target.y());

                    // Max.
                    maxx = CGAL::max(maxx, source.x());
                    maxx = CGAL::max(maxx, target.x());

                    maxy = CGAL::max(maxy, source.y());
                    maxy = CGAL::max(maxy, target.y());
                }

                bl = Point_2(minx, miny);

                const FT lengthx = CGAL::abs(maxx - minx);
                const FT lengthy = CGAL::abs(maxy - miny);

                const FT bbox_width  = lengthx;
                const FT bbox_height = lengthy;

                bbox_size = std::make_pair(bbox_width, bbox_height);
            }

            void translate_segments(const Point_2 &bl, Segments &segments, const bool save) const {
                assert(segments.size() > 0);

                Point_2 new_source, new_target;
                for (size_t i = 0; i < segments.size(); ++i) {
                    
                    const Point_2 &source = segments[i].source();
                    const Point_2 &target = segments[i].target();

                    const FT x1 = source.x() - bl.x();
                    const FT y1 = source.y() - bl.y();

                    const FT x2 = target.x() - bl.x();
                    const FT y2 = target.y() - bl.y();

                    new_source = Point_2(x1, y1);
                    new_target = Point_2(x2, y2);

                    segments[i] = Segment_2(new_source, new_target);
                }

                if (m_debug && save) {
                    const std::string stub = "";
                    Log segments_exporter; segments_exporter.export_segments_as_obj("tmp" + std::string(PSR) + "polygonizer_translated_input_segments_jean_philippe", segments, stub);
                }
            }

            FT get_average_spacing(const Segments &segments) const {

                using Local_Kernel = CGAL::Simple_cartesian<double>;
				using Point_3ft    = Local_Kernel::Point_3;

                assert(segments.size() > 0);
                std::vector<Point_3ft> points(segments.size() * 2);

                size_t count = 0;
                for (size_t i = 0; i < segments.size(); ++i) {

                    const Point_2 &source = segments[i].source();
                    const Point_2 &target = segments[i].target();

                    points[count++] = Point_3ft(CGAL::to_double(source.x()), CGAL::to_double(source.y()), 0.0);
                    points[count++] = Point_3ft(CGAL::to_double(target.x()), CGAL::to_double(target.y()), 0.0);
                }

                const double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(points.begin(), points.end(), CGAL::Identity_property_map<Point_3ft>(), m_num_neighbours, Local_Kernel());
                return static_cast<FT>(average_spacing) / m_local_scaling;
            }

            void get_scale(const Bbox_size &bbox_size, const FT average_spacing, Scale &scale) const {

                const FT x = bbox_size.first;
                const FT y = bbox_size.second;

                assert(x > FT(0) && y > FT(0) && average_spacing > FT(0));

                const FT x_num = x / average_spacing;
                const FT y_num = y / average_spacing;

                scale = std::make_pair(x_num, y_num);
            }

            void set_input_segments(const Segments &segments, const Scale &scale, Kinetic_Model *model) const {
                
                assert(segments.size() > 0);
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

            void compute_partition(Kinetic_Model *model, const Bbox_size &bbox_size, const Scale &scale) const {

                Propagation propagation;

                const size_t rows = static_cast<size_t>(std::ceil(CGAL::to_double(scale.second * bbox_size.second)));
                const size_t cols = static_cast<size_t>(std::ceil(CGAL::to_double(scale.first  * bbox_size.first)));

                propagation.dmitry_size_rows = rows;
                propagation.dmitry_size_cols = cols;

                model->set_prop_ttl(static_cast<int>(m_num_intersections));
                model->set_prop_merge_min_thinness(CGAL::to_double(m_min_face_width));

                propagation.propagate(model);
            }

            void get_back_segments(Segments &segments, Kinetic_Model *model, const Scale &scale) const {
                segments.clear();

                const FT scale_x = scale.first;
                const FT scale_y = scale.second;

	            Edge* e = model->graph->edges_head;
                while (e != NULL) {

                    const Point2d &pt_1 = e->v1->pt;
                    const Point2d &pt_2 = e->v2->pt;

                    segments.push_back(Segment_2(Point_2(FT(pt_1.x) / scale_x, FT(pt_1.y) / scale_y),  Point_2(FT(pt_2.x) / scale_x, FT(pt_2.y) / scale_y)));
                    e = e->e_next;
                }
            }

            void save_partition(const Segments &segments) const {

                const std::string stub = "";
                Log segments_exporter; segments_exporter.export_segments_as_obj("tmp" + std::string(PSR) + "polygonizer_partition_jean_philippe", segments, stub);
            }

            void built_data_structure(const Point_2 &bl, const Scale &scale, Kinetic_Model *model, Data_structure &data_structure) const {

                Partition *graph = model->graph;
                const JP_Faces &jp_faces = graph->faces;

                built_polygons(bl, jp_faces, scale, data_structure);
                if (false) save_polygons(data_structure);
            }

            void built_polygons(const Point_2 &bl, const JP_Faces &jp_faces, const Scale &scale, Data_structure &data_structure) const {
                assert(jp_faces.size() > 0);

                Containers &containers = data_structure.containers();
                
                containers.clear();
                containers.resize(jp_faces.size());

                size_t i = 0;
                for (typename JP_Faces::const_iterator fit = jp_faces.begin(); fit != jp_faces.end(); ++fit, ++i)
                    create_polygon(bl, *fit, scale, i, containers);
            }

            void create_polygon(const Point_2 &bl, const JP_Face *jp_face, const Scale &scale, const size_t container_index, Containers &containers) const {

                const JP_Vertices &jp_face_vertices = jp_face->vertices;
                
                assert(jp_face_vertices.size() > 2);
                Points points(jp_face_vertices.size());

                const FT scale_x = scale.first;
                const FT scale_y = scale.second;

                size_t i = 0;
                for (typename JP_Vertices::const_iterator vit = jp_face_vertices.begin(); vit != jp_face_vertices.end(); ++vit, ++i) {
                 
                    const JP_Point &point = vit->first->pt;

                    const FT x = FT(point.x) / scale_x - bl.x();
                    const FT y = FT(point.y) / scale_y - bl.y();

                    points[i] = Point_2(x, y);
                }
                containers[container_index].polygon = Polygon(points.begin(), points.end());
            }

            void save_polygons(const Data_structure &data_structure) const {

                const Containers &containers = data_structure.containers();
                assert(containers.size() > 0);

                Log exporter;
                exporter.save_polygons<Containers, Polygon, Kernel>(containers, "tmp" + std::string(PSR) + "polygonizer_polygons_jean_philippe");
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_POLYGONIZER_JEAN_PHILIPPE_H