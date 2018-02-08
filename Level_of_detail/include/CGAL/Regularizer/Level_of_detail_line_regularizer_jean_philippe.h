#ifndef CGAL_LEVEL_OF_DETAIL_LINE_REGULARIZER_JEAN_PHILIPPE_H
#define CGAL_LEVEL_OF_DETAIL_LINE_REGULARIZER_JEAN_PHILIPPE_H

#if defined(WIN32) || defined(_WIN32) 
#define PS "\\"
#else 
#define PS "/" 
#endif 

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <utility>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/enum.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>

// Jean Philippe includes.
#include <CGAL/Regularizer/jean_philippe/kinetic_model.h>
#include <CGAL/Regularizer/jean_philippe/regularization_angles.h>
#include <CGAL/Regularizer/jean_philippe/regularization_angles_quadratic.h>
#include <CGAL/Regularizer/jean_philippe/regularization_ordinates.h>
#include <CGAL/Regularizer/jean_philippe/regularization_ordinates_quadratic.h>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits, class BoundaryData, class ProjectedPoints>
        class Level_of_detail_line_regularizer_jean_philippe { 

        public:
            typedef KernelTraits    Kernel;
            typedef BoundaryData 	Boundary_data;
			typedef ProjectedPoints Projected_points;

			typedef typename Kernel::Point_2   Point_2;
			typedef typename Kernel::Line_2    Line_2;
            typedef typename Kernel::Segment_2 Segment_2;
			typedef typename Kernel::FT 	   FT;

            typedef std::vector<Line_2>           Lines;
            typedef std::vector<Segment_2>        Segments;

            using Log = CGAL::LOD::Mylog;

            Level_of_detail_line_regularizer_jean_philippe() : 
            m_silent(false), 
            m_debug(false),
            m_add_ordinates(false)
            { }

			void make_silent(const bool state) {
				m_silent = state;
			}

            void add_ordinates(const bool state) {
                m_add_ordinates = state;
            }

            void process(const Boundary_data &, const Projected_points &, const Segments &segments, Lines &lines) {

                if (m_debug) {
                    const std::string stub = "";
                    Log segments_exporter; segments_exporter.export_segments_as_obj("tmp" + std::string(PS) + "input_segments_jean_philippe", segments, stub);
                }

                Kinetic_Model* model = new Kinetic_Model();
                
                initialize_kinetic_model(model);
                set_segments(segments, model);

                regularize_segments(model);
                get_back_lines(model, lines);

                delete model;
            }

            inline Segments &get_regularized_segments() {
                return m_regularized_segments;
            }

        private:
            bool m_silent;
            bool m_debug;
            bool m_add_ordinates;

            Segments m_regularized_segments;

            void initialize_kinetic_model(Kinetic_Model *model) {
                model->reinit();
            }

            void set_segments(const Segments &segments, Kinetic_Model *model) {
                std::vector<Segment *> &model_segments = model->segments;
                
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

            void regularize_segments(Kinetic_Model *model) {

                regularize_angles(model);
                if (m_add_ordinates) regularize_ordinates(model);
            }

            void regularize_angles(Kinetic_Model *model) {
                
                Regularization_Angles* m_rega = nullptr;
                m_rega = new Regularization_Angles_Quadratic();
                m_rega->regularize(model);
                delete m_rega;
            }

            void regularize_ordinates(Kinetic_Model *model) {
                
                Regularization_Ordinates* m_regp = nullptr;
                m_regp = new Regularization_Ordinates_Quadratic();
                m_regp->regularize(model);
		        delete m_regp;
            }

            void get_back_lines(Kinetic_Model *model, Lines &lines) {
                get_segments(model, m_regularized_segments);

                if (!m_silent) {
                    const std::string stub = "";
                    Log segments_exporter; segments_exporter.export_segments_as_obj("tmp" + std::string(PS) + "regularized_segments_jean_philippe", m_regularized_segments, stub);
                }

                get_lines(m_regularized_segments, lines);
            }

            void get_segments(Kinetic_Model *model, Segments &segments) {
                
                std::vector<Segment *> &model_segments = model->segments;
                segments.clear();

                segments.reserve(model_segments.size());
                for (size_t i = 0; i < model_segments.size(); ++i) {

                    const Point2d source = model_segments[i]->finalEnd1;
                    const Point2d target = model_segments[i]->finalEnd2;

                    const FT x1 = static_cast<FT>(source.x);
                    const FT y1 = static_cast<FT>(source.y);

                    const FT x2 = static_cast<FT>(target.x);
                    const FT y2 = static_cast<FT>(target.y);

                    segments.push_back(Segment_2(Point_2(x1, y1), Point_2(x2, y2)));
                }
            }

            void get_lines(const Segments &segments, Lines &lines) {

                lines.clear();
                lines.resize(segments.size());

                for (size_t i = 0; i < segments.size(); ++i) {

                    const Point_2 &source = segments[i].source();
                    const Point_2 &target = segments[i].target();

                    lines[i] = Line_2(source, target);
                }

                if (m_debug) {
                    const std::string stub = "";
                    Log lines_exporter; lines_exporter.export_lines_as_obj("tmp" + std::string(PS) + "regularized_lines_jean_philippe", lines, stub);
                }
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_LINE_REGULARIZER_JEAN_PHILIPPE_H