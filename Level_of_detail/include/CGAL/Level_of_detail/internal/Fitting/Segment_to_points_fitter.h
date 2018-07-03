#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_TO_POINTS_FITTER_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_TO_POINTS_FITTER_H

// CGAL includes.
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Level_of_detail/internal/Fitting/Line_to_points_fitter.h>
#include <CGAL/Level_of_detail/internal/Estimations/End_points_estimator.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

		template<class InputKernel>
		class Segment_to_points_fitter {

		public:
			using Kernel = InputKernel;

            using FT        = typename Kernel::FT;
            using Line_2    = typename Kernel::Line_2;
            using Point_2   = typename Kernel::Point_2;
            using Vector_2  = typename Kernel::Vector_2;
            using Segment_2 = typename Kernel::Segment_2;

            using Line_to_points_fitter    = LOD::Line_to_points_fitter<Kernel>;
            using End_points_estimator     = LOD::End_points_estimator<Kernel>;

            using Projected_points   = std::list<Point_2>;
            using Identity_point_map = CGAL::Identity_property_map<Point_2>;

			template<class Elements, class Point_map>
			FT fit_segment_2(const Elements &elements, const Point_map &point_map, Segment_2 &segment) const {
                CGAL_precondition(elements.size() > 1);

                // Fit line to all points.
                Line_2 line;
                const Line_to_points_fitter line_to_points_fitter;
				const FT quality = line_to_points_fitter.fit_line_2(elements, point_map, line);

                // Project points onto the line.
                Projected_points projected_points;
                for (typename Elements::const_iterator ce_it = elements.begin(); ce_it != elements.end(); ++ce_it)
                  projected_points.push_back(line.projection(get(point_map, *ce_it)));
                
                // Find end points of the segment;
                Point_2 a, b;
                Identity_point_map identity_point_map;

                const End_points_estimator end_points_estimator;
                end_points_estimator.find_end_points_wrt_barycentre_2(projected_points, identity_point_map, a, b);

                segment = Segment_2(a, b);
                return quality;
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_TO_POINTS_FITTER_H
