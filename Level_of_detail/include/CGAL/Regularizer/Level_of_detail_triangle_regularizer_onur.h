#ifndef CGAL_LEVEL_OF_DETAIL_TRIANGLE_REGULARIZER_ONUR_H
#define CGAL_LEVEL_OF_DETAIL_TRIANGLE_REGULARIZER_ONUR_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\"
#else 
#define PSR "/" 
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

namespace CGAL {

	namespace LOD {

        template<class KernelTraits, class BoundaryData, class ProjectedPoints>
        class Level_of_detail_triangle_regularizer_onur { 

        public:
            typedef KernelTraits    Kernel;
            typedef BoundaryData 	Boundary_data;
			typedef ProjectedPoints Projected_points;

			typedef typename Kernel::Point_2   Point_2;
            typedef typename Kernel::Point_3   Point_3;
            typedef typename Kernel::Vector_2  Vector_2;
            typedef typename Kernel::Vector_3  Vector_3;
			typedef typename Kernel::Line_2    Line_2;
            typedef typename Kernel::Line_3    Line_3;
            typedef typename Kernel::Segment_2 Segment_2;
            typedef typename Kernel::Plane_3   Plane_3;
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

            Level_of_detail_triangle_regularizer_onur() : 
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

            void process(const Boundary_data & /* building_boundaries */, const Projected_points & /* building_boundaries_projected */, const Segments & /* segments */, Lines & /* lines */) const {

                assert(!"This method is not yet implemented!");
            }

        private:
            bool m_silent;
            bool m_debug;

            const bool m_parallelism;
            const bool m_orthogonality;
            const bool m_coplanarity;
            const bool m_z_symmetry;

            const FT m_angle;
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_TRIANGLE_REGULARIZER_ONUR_H