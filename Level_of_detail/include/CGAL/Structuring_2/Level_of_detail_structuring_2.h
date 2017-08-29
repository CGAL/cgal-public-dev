#ifndef CGAL_LEVEL_OF_DETAIL_STRUCTURING_2_H
#define CGAL_LEVEL_OF_DETAIL_STRUCTURING_2_H

// STL includes.
#include <map>
#include <vector>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits>
		class Level_of_detail_structuring_2 {

		public:
			typedef KernelTraits Traits;

			typedef typename Traits::Point_2   Point;
			typedef typename Traits::Line_2    Line;
			typedef typename Traits::Segment_2 Segment;

			typedef std::map<int, Point>             Projected_points;
			typedef std::map<int, std::vector<int> > Connected_components;
			typedef std::vector<Line> 				 Lines;
			typedef std::vector<Segment> 			 Segments;

			// Remove later.
			using Log = CGAL::LOD::Mylog;

			Level_of_detail_structuring_2(const Projected_points &points, const Connected_components &components, const Lines &lines) :
			m_points(points), m_cc(components), m_lines(lines) { }

			int structure_point_set(Segments &) {

				// (START) Create log.
				Log log; log.out << "START EXECUTION\n\n\n";
				
				// -----------------------------------------

				const auto number_of_structured_segments = -1;

				// to be implemented!


				// -------------------------------

				// (END) Save log.
				log.out << "\n\nFINISH EXECUTION";
				log.save("structuring_2");

				return number_of_structured_segments;
			}

		private:
			const Projected_points &m_points;
			const Connected_components &m_cc;
			const Lines &m_lines;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_STRUCTURING_2_H