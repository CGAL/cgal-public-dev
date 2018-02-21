#ifndef CGAL_LEVEL_OF_DETAIL_POLYGONIZER_JEAN_PHILIPPE_H
#define CGAL_LEVEL_OF_DETAIL_POLYGONIZER_JEAN_PHILIPPE_H

#if defined(WIN32) || defined(_WIN32) 
#define PS "\\"
#else 
#define PS "/" 
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
            typedef typename Kernel::Segment_2 Segment_2;

            using Segments       = std::vector<Segment_2>;
            using Model_segments = std::vector<Segment *>;

            using Bbox_size = std::pair<size_t, size_t>;
            using Log = CGAL::LOD::Mylog;

            Level_of_detail_polygonizer_jean_philippe() :
            m_silent(false), m_debug(false), m_num_intersections(2), m_min_face_width(FT(3)) { }

            void polygonize(Segments &segments) const {

                if (m_debug) {
                    const std::string stub = "";
                    Log segments_exporter; segments_exporter.export_segments_as_obj("tmp" + std::string(PS) + "polygonizer_input_segments_jean_philippe", segments, stub);
                }

                Kinetic_Model* model = new Kinetic_Model();
                
                // Initialize the model.
                initialize_kinetic_model(model);

                // Find bottom left corner of the segments' bounding box.
                Point_2 bl; Bbox_size bbox_size;
                find_bottom_left_corner(segments, bl, bbox_size);

                // Translate all segments so that bottom left corner becomes (0, 0);
                translate_segments(bl, segments, true);

                // Set input segments to the model.
                set_input_segments(segments, model);
                
                // Compute partition.
                compute_partition(model, bbox_size);

                // Get back segments.
                get_back_segments(segments, model);

                // Get reverse translation factor.
                bl = Point_2(-bl.x(), -bl.y());

                // Translate back all segments.
                translate_segments(bl, segments, false);

                // Save results.
                if (!m_silent) save_partition(segments);

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

            void built_data(Data_structure & /* data_structure */) {

            }

        private:
            bool m_silent;
            bool m_debug;

            size_t m_num_intersections;
            FT     m_min_face_width;

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

                const FT lengthx = maxx - minx;
                const FT lengthy = maxy - miny;

                const size_t bbox_width  = static_cast<size_t>(std::ceil(CGAL::to_double(lengthx)));
                const size_t bbox_height = static_cast<size_t>(std::ceil(CGAL::to_double(lengthy)));

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
                    Log segments_exporter; segments_exporter.export_segments_as_obj("tmp" + std::string(PS) + "polygonizer_translated_input_segments_jean_philippe", segments, stub);
                }
            }

            void set_input_segments(const Segments &segments, Kinetic_Model *model) const {
                
                assert(segments.size() > 0);
                Model_segments &model_segments = model->segments;
                
                model_segments.reserve(segments.size());
                for (size_t i = 0; i < segments.size(); ++i) {

                    const Point_2 &source = segments[i].source();
                    const Point_2 &target = segments[i].target();

                    const double x1 = CGAL::to_double(source.x());
                    const double y1 = CGAL::to_double(source.y());

                    const double x2 = CGAL::to_double(target.x());
                    const double y2 = CGAL::to_double(target.y());

                    const double width = CGAL::sqrt(CGAL::to_double(segments[i].squared_length()));
                    model_segments.push_back(new Segment(i, x1, y1, x2, y2, width, 0.0, 0.0, 0.0, false));
                }
            }

            void compute_partition(Kinetic_Model *model, const Bbox_size &bbox_size) const {

                Propagation propagation;

                const size_t rows = bbox_size.second;
                const size_t cols = bbox_size.first;

                propagation.dmitry_size_rows = rows;
                propagation.dmitry_size_cols = cols;

                model->set_prop_ttl(static_cast<int>(m_num_intersections));
                model->set_prop_merge_min_thinness(CGAL::to_double(m_min_face_width));

                propagation.propagate(model);
            }

            void get_back_segments(Segments &segments, Kinetic_Model *model) const {
                segments.clear();

	            Edge* e = model->graph->edges_head;
                while (e != NULL) {
                    const Point2d &pt_1 = e->v1->pt;
                    const Point2d &pt_2 = e->v2->pt;

                    segments.push_back(Segment_2(Point_2(pt_1.x, pt_1.y),  Point_2(pt_2.x, pt_2.y)));
                    e = e->e_next;
                }
            }

            void save_partition(const Segments &segments) const {
                const std::string stub = "";
                Log segments_exporter; segments_exporter.export_segments_as_obj("tmp" + std::string(PS) + "polygonizer_partition_jean_philippe", segments, stub);
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_POLYGONIZER_JEAN_PHILIPPE_H