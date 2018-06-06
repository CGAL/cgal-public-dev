#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_TO_POINTS_FITTER_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_TO_POINTS_FITTER_H

// LOD includes.
#include <CGAL/Level_of_detail/Tools/Fitters/Line_to_points_fitter.h>

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

            using Line_to_points_fitter = LOD::Line_to_points_fitter<Kernel>;

            Segment_to_points_fitter(const FT tolerance = FT(1) / FT(1000000)) :
            m_big_value(FT(100000000000000)),
            m_tolerance(tolerance)
            { }

			template<class Elements, class Point_map>
			FT fit_segment_2(const Elements &elements, const Point_map &point_map, Segment_2 &segment) const {
				
                CGAL_precondition(elements.size() > 0);
                using Const_elements_iterator = typename Elements::const_iterator;

                // Fit line to all points.
                Line_2 line;
                const Line_to_points_fitter line_to_points_fitter;
				const FT quality = line_to_points_fitter.fit_line_2(elements, point_map, line);

                // Create segment.
                FT minx =  m_big_value, miny =  m_big_value;
			    FT maxx = -m_big_value, maxy = -m_big_value;

                for (Const_elements_iterator element = elements.begin(); element != elements.end(); ++element) {
					
                    const Point_2 &point    = get(point_map, *element);
                    const Point_2 projected = line.projection(point);

                    minx = CGAL::min(minx, projected.x());
                    maxx = CGAL::max(maxx, projected.x());

                    miny = CGAL::min(miny, projected.y());
					maxy = CGAL::max(maxy, projected.y());
                }
                segment = Segment_2(Point_2(minx, miny), Point_2(maxx, maxy));

                // Rotate segment if needed.
                const Vector_2 v1 = line.to_vector();
			    const Vector_2 v2 = segment.to_vector();
			
				if ((v1.y() < FT(0) && v2.y() >= FT(0) && CGAL::abs(v1.y() - v2.y()) > m_tolerance) ||
				    (v2.y() < FT(0) && v1.y() >= FT(0) && CGAL::abs(v1.y() - v2.y()) > m_tolerance) ||
					(v1.x() < FT(0) && v2.x() >= FT(0) && CGAL::abs(v1.x() - v2.x()) > m_tolerance) ||
					(v2.x() < FT(0) && v1.x() >= FT(0) && CGAL::abs(v1.x() - v2.x()) > m_tolerance)) {

					segment = Segment_2(Point_2(minx, maxy), Point_2(maxx, miny));
				}
                return quality;
            }

        private:
            const FT m_big_value;
            const FT m_tolerance;
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_TO_POINTS_FITTER_H