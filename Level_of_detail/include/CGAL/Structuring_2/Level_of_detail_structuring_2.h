#ifndef CGAL_LEVEL_OF_DETAIL_STRUCTURING_2_H
#define CGAL_LEVEL_OF_DETAIL_STRUCTURING_2_H

// STL includes.
#include <map>
#include <vector>
#include <cmath>

// CGAL includes.
#include <CGAL/number_utils.h>
#include <CGAL/utils.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>

// NOTES:
// 1. Clutter has to be handled separately. See Section 2 Clutter.
// 2. All points associated with a line should come with epsilon parameter.

namespace CGAL {

	namespace LOD {

		template<class KernelTraits>
		class Level_of_detail_structuring_2 {

		public:
			typedef KernelTraits Traits;

			typedef typename Traits::Point_2   Point;
			typedef typename Traits::Line_2    Line;
			typedef typename Traits::Segment_2 Segment;
			typedef typename Traits::FT 	   FT;

			typedef std::map<int, Point>             Points;
			typedef std::map<int, std::vector<int> > Connected_components;
			typedef std::vector<Line> 				 Lines;
			typedef std::vector<Segment> 			 Segments;

			typedef typename Connected_components::const_iterator CC_iterator;

			// Remove later or change.
			enum class Structured_label { LINEAR, CORNER };
			using Log = CGAL::LOD::Mylog;

			using Structured_points  = std::vector<Point>; 			   // new structured points
			using Structured_labels  = std::vector<Structured_label>;  // for each structured point store a label: linear or corner
			using Structured_anchors = std::vector<std::vector<int> >; // for each structured point store all associated primitives: indices of lines (may be 1 or 2)

			typename Traits::Compute_squared_distance_2 squared_distance;

			Level_of_detail_structuring_2(const Points &points, const Connected_components &components, const Lines &lines) :
			m_points(points), m_cc(components), m_lines(lines) { 

				assert(components.size() == lines.size());
			}

			// This is a 2D version of the algorithm in the paper:
			// Surface Reconstruction Through Point Set Structuring, F. Lafarge, P. Alliez.
			int structure_point_set(Segments &) {

				// (START) Create log.
				Log log; log.out << "START EXECUTION\n\n\n";
				
				// -----------------------------------------

				const auto number_of_structured_segments = -1;


				// (1) Project all points onto the given lines.
				project();

				log.out << "(1) All points are projected onto the given lines. The results are saved in tmp/projected" << std::endl;


				// (2) Compute min, average, and max distances from all points to the given lines.
				compute_distances();

				log.out << "(2) Min, avg, and max distances for each set of points are computed. The results are saved in tmp/distances" << std::endl;


				// (3) Find epsilon for each set of points.
				// Alternatively it can be set as one unique value using the corresponding function.


				// (4) Find one unique segment for each set of points.


				// -------------------------------

				// (END) Save log.
				log.out << "\n\nFINISH EXECUTION";
				log.save("structuring_2");

				return number_of_structured_segments;
			}

			const Structured_points& get_structured_points() const {
				return m_str_points;
			}

			const Structured_labels& get_structured_labels() const {
				return m_str_labels;
			}

			const Structured_anchors& get_structured_anchors() const {
				return m_str_anchors;
			}

			void set_epsilon(const FT value) {

				assert(value > FT(0));

				m_eps.clear();
				m_eps.resize(m_cc.size(), value);
			}

		private:
			const Points 			   &m_points;
			const Connected_components &m_cc;
			const Lines                &m_lines;

			std::vector<FT> m_min_dist, m_avg_dist, m_max_dist;
			std::vector<FT> m_eps;

			Structured_points  m_str_points;
			Structured_labels  m_str_labels;
			Structured_anchors m_str_anchors;

			Points m_projected;

			void project() {

				clear_projected();
				assert(m_projected.size() == m_points.size());
			}

			void clear_projected() {

				m_projected.clear();
				m_projected.resize(m_points.size());
			}

			void compute_distances() {
				
				clear_distances();

				Log log;

				assert(m_min_dist.size() == m_cc.size());
				assert(m_avg_dist.size() == m_cc.size());
				assert(m_max_dist.size() == m_cc.size());

				assert(m_projected.size() == m_points.size());

				const FT big_value = FT(1000000);
				size_t number_of_components = 0;

				for (CC_iterator it = m_cc.begin(); it != m_cc.end(); ++it, ++number_of_components) {
					const auto num_points = (*it).second.size();

					m_min_dist[number_of_components] =  big_value; 
					m_avg_dist[number_of_components] =      FT(0); 
					m_max_dist[number_of_components] = -big_value;

					for (size_t i = 0; i < num_points; ++i) {
						
						const auto index = (*it).second[i];
						assert(index >= 0);

						const Point &p = m_points.at(index);
						const Point &q = m_projected.at(index);

						const auto dist = CGAL::sqrt(squared_distance(p, q));

						m_min_dist[number_of_components] = CGAL::min(m_min_dist[number_of_components], dist);
						m_avg_dist[number_of_components] += dist;
						m_max_dist[number_of_components] = CGAL::max(m_max_dist[number_of_components], dist);
					}
					m_avg_dist[number_of_components] /= static_cast<FT>(num_points);

					log.out << "  min: " << m_min_dist[number_of_components] << 
					           "; avg: " << m_avg_dist[number_of_components] <<
					           "; max: " << m_max_dist[number_of_components] << std::endl;

					log.save("tmp/distances");
				}
			}

			void clear_distances() {

				m_min_dist.clear(); m_avg_dist.clear(); m_max_dist.clear();
				m_min_dist.resize(m_cc.size(), FT(0)); m_avg_dist.resize(m_cc.size(), FT(0)); m_max_dist.resize(m_cc.size(), FT(0));
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_STRUCTURING_2_H