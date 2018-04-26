#ifndef CGAL_LEVEL_OF_DETAIL_LINEAR_REGION_GROWING_H
#define CGAL_LEVEL_OF_DETAIL_LINEAR_REGION_GROWING_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits>
		class Level_of_detail_linear_region_growing {

		public:
			typedef KernelTraits Kernel;

			using FT 	     = typename Kernel::FT;
			using Point_2    = typename Kernel::Point_2;
			using Triangle_2 = typename Kernel::Triangle_2;
			using Segment_2  = typename Kernel::Segment_2;
	
			using Index   = int;
			using Indices = std::vector<Index>;

			using States   = std::vector<bool>;
			using Segments = std::vector<Segment_2>;
			using Output   = std::vector<Indices>;

			Level_of_detail_linear_region_growing() : 
			m_tolerance(FT(1) / FT(100000)) 
			{ }

			void find_connected_segments(const Segments &segments, const States &to_be_used, Output &output) const {
				
				output.clear();
				assert(segments.size() == to_be_used.size());

				if (segments.size() == 0) return;
				States was_used(segments.size(), false);

				for (size_t i = 0; i < segments.size(); ++i) {
					if (to_be_used[i] && !was_used[i]) {
					
						Indices indices;
						grow_region(segments, to_be_used, i, was_used, indices);

						assert(indices.size() != 0);
						output.push_back(indices);
					}
				}
			}

		private:
			const FT m_tolerance;

			void grow_region(const Segments &segments, const States &to_be_used, const size_t segment_index, States &was_used, Indices &indices) const {
				
				assert(was_used.size() == segments.size());
				assert(segment_index >= 0 && segment_index < segments.size());
				
				was_used[segment_index] = true;
				indices.push_back(segment_index);

				for (size_t i = 0; i < segments.size(); ++i) {
					if (to_be_used[i] && !was_used[i]) {

						if (go_together(segments[segment_index], segments[i]))
							grow_region(segments, to_be_used, i, was_used, indices);
					}
				}
			}

			bool go_together(const Segment_2 &seg1, const Segment_2 &seg2) const {

				if (are_close_by(seg1, seg2) && are_collinear(seg1, seg2)) return true;
				return false;
			}

			bool are_close_by(const Segment_2 &seg1, const Segment_2 &seg2) const {

				const Point_2 &p1 = seg1.source();
				const Point_2 &q1 = seg1.target();

				const Point_2 &p2 = seg2.source();
				const Point_2 &q2 = seg2.target();

				if (is_close_by(p1, p2)) return true;
				if (is_close_by(p1, q2)) return true;
				if (is_close_by(q1, p2)) return true;
				if (is_close_by(q1, q2)) return true;

				return false;
			}

			bool is_close_by(const Point_2 &p, const Point_2 &q) const {

				const FT eps = m_tolerance;
				if (CGAL::abs(p.x() - q.x()) < eps && CGAL::abs(p.y() - q.y()) < eps) return true;
				return false;
			}

			bool are_collinear(const Segment_2 &seg1, const Segment_2 &seg2) const {
				
				const Point_2 &p1 = seg1.source();
				const Point_2 &q1 = seg1.target();

				const Point_2 &p2 = seg2.source();
				const Point_2 &q2 = seg2.target();

				if (are_quasi_collinear(p1, q1, p2) && are_quasi_collinear(p1, q1, q2)) return true;
				return false;
			}

			bool are_quasi_collinear(const Point_2 &p, const Point_2 &q, const Point_2 &r) const {
				if (is_close_by(p, r) || is_close_by(q, r)) return true;

				const FT eps = m_tolerance;
				const Triangle_2 triangle = Triangle_2(p, q, r);

				if (CGAL::abs(triangle.area()) < eps) return true;
				return false;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_LINEAR_REGION_GROWING_H