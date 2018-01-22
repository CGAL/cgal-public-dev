#ifndef CGAL_LEVEL_OF_DETAIL_LINE_REGULARIZER_CGAL_H
#define CGAL_LEVEL_OF_DETAIL_LINE_REGULARIZER_CGAL_H

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
#include <CGAL/Iterator_range.h>
#include <CGAL/property_map.h>
#include <CGAL/enum.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Regularizer/Level_of_detail_regularize_planes.h>
#include <CGAL/Regularizer/Level_of_detail_regularize_property_maps.h>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits, class BoundaryData, class ProjectedPoints>
        class Level_of_detail_line_regularizer_cgal { 

        public:
            typedef KernelTraits    Kernel;
            typedef BoundaryData 	Boundary_data;
			typedef ProjectedPoints Projected_points;

			typedef typename Kernel::Point_2   Point_2;
            typedef typename Kernel::Point_3   Point_3;
			typedef typename Kernel::Line_2    Line_2;
            typedef typename Kernel::Line_3    Line_3;
            typedef typename Kernel::Segment_2 Segment_2;
            typedef typename Kernel::Plane_3   Plane_3;
            typedef typename Kernel::Vector_2  Vector_2;
            typedef typename Kernel::Vector_3  Vector_3;
			typedef typename Kernel::FT 	   FT;

            typedef std::pair<Point_3, int>                Point_with_plane;
            typedef typename Boundary_data::const_iterator Boundary_iterator;

            typedef typename CGAL::First_of_pair_property_map<Point_with_plane>  Point_map;
            typedef typename CGAL::Second_of_pair_property_map<Point_with_plane> Index_map;
            typedef typename CGAL::Identity_property_map<Plane_3>                Plane_map;

            typedef std::vector<Point_with_plane> Points;
            typedef std::vector<Plane_3>          Planes;
            typedef std::vector<Line_2>           Lines;
            typedef std::vector<Segment_2>        Segments;

            typedef std::vector<Point_3>    Debug_quad;
            typedef std::vector<Debug_quad> Debug_quads;

            using Log = CGAL::LOD::Mylog;

            Level_of_detail_line_regularizer_cgal() : 
            m_silent(false), 
            m_debug(false),
            m_parallelism(false),
            m_orthogonality(false),
            m_coplanarity(false),
            m_z_symmetry(false),
            m_angle(FT(25))
            { }

			void make_silent(const bool silent) {
				m_silent = silent;
			}

            void process(const Boundary_data &building_boundaries, const Projected_points &building_boundaries_projected, const Segments &segments, Lines &lines) const {
                
                assert(!"This method is not working! Do not use it!");
                assert(building_boundaries.size() == segments.size());

                Points points;
                create_points(building_boundaries, building_boundaries_projected, points);

                Planes planes;
                create_planes(segments, planes);

                regularize(points, planes);
                get_lines(planes, lines);
            }

        private:
            bool m_silent;
            bool m_debug;

            const bool m_parallelism;
            const bool m_orthogonality;
            const bool m_coplanarity;
            const bool m_z_symmetry;

            const FT m_angle;

            void create_points(const Boundary_data &building_boundaries, const Projected_points &building_boundaries_projected, Points &points) const {

                points.clear();
                int plane_index = 0;

                for (Boundary_iterator bit = building_boundaries.begin(); bit != building_boundaries.end(); ++bit, ++plane_index) {
                    
                    const std::vector<int> &indices = (*bit).second;
                    for (size_t i = 0; i < indices.size(); ++i) {
                        
                        const Point_2 &point_2 = building_boundaries_projected.at(indices[i]);
                        const Point_3  point_3 = Point_3(point_2.x(), point_2.y(), FT(0));
                            
                        const Point_with_plane pwp = std::make_pair(point_3, plane_index);
                        points.push_back(pwp);
                    }
                }

                if (m_debug) {
                    Log saver; saver.export_points_with_planes_as_xyz("tmp" + std::string(PS) + "regularization_points", points);
                }
            }

            void create_planes(const Segments &segments, Planes &planes) const {
                Debug_quads debug_quads(segments.size());

                planes.clear();
                planes.resize(segments.size());

                int plane_index = 0;
                for (size_t i = 0; i < segments.size(); ++i, ++plane_index) {
                 
                    const Plane_3 plane = create_plane_from_segment(segments[i], debug_quads[i]);
                    planes[plane_index] = plane;
                }

                if (m_debug) {
                    Log saver; saver.save_quads_as_ply<Debug_quads, Point_3>(debug_quads, "tmp" + std::string(PS) + "regularization_planes");
                }
            }

            Plane_3 create_plane_from_segment(const Segment_2 &segment, Debug_quad &debug_quad) const {

                const Point_2 &p1 = segment.source();
                const Point_2 &p2 = segment.target();

                const Point_3 a = Point_3(p1.x(), p1.y(), FT(0));
                const Point_3 b = Point_3(p2.x(), p2.y(), FT(0));
                const Point_3 c = Point_3(p2.x(), p2.y(), FT(10));

                if (m_debug) {
                    const Point_3 d = Point_3(p1.x(), p1.y(), FT(10));    
                    debug_quad.resize(4);

                    debug_quad[0] = a;
                    debug_quad[1] = b;
                    debug_quad[2] = c;
                    debug_quad[3] = d;
                }

                return Plane_3(a, b, c);
            }

            void regularize(const Points &points, Planes &planes) const {
                
                CGAL::LOD::regularize_planes(points, Point_map(), planes, Plane_map(), CGAL::LOD::Point_to_shape_index_map<Kernel>(points),
                m_parallelism, m_orthogonality, m_coplanarity, m_z_symmetry);
            }

            void get_lines(const Planes &planes, Lines &lines) const {

                lines.clear();
                lines.resize(planes.size());

                for (size_t i = 0; i < planes.size(); ++i) {

                    const  Point_3 point_3      = planes[i].point();
                    const   Line_3 pnd_line_3   = planes[i].perpendicular_line(point_3);
                    const Vector_3 pnd_vector_3 = pnd_line_3.to_vector();
                    const Vector_2 pnd_vector_2 = Vector_2(pnd_vector_3.x(), pnd_vector_3.y());
                    const Vector_2 vector_2     = pnd_vector_2.perpendicular(CGAL::COUNTERCLOCKWISE);
                    const  Point_2 point_2      = Point_2(point_3.x(), point_3.y());
                    const   Line_2 line_2       = Line_2(point_2, vector_2);

                    lines[i] = line_2;
                }
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_LINE_REGULARIZER_CGAL_H