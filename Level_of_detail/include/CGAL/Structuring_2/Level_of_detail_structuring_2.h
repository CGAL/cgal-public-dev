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
			typename Traits::Compute_scalar_product_2   dot_product;

			Level_of_detail_structuring_2(const Points &points, const Connected_components &components, const Lines &lines) :
			m_points(points), m_cc(components), m_lines(lines), m_tol(FT(1) / FT(1000000)), m_big_value(FT(1000000)) { 

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
				compute_epsilon_values();

				log.out << "(3) Epsilon values are computed for each component. The results are saved in tmp/epsilons" << std::endl;


				// (4) Find one unique segment for each set of points.
				find_segments();

				log.out << "(4) Segments are extracted for each set of points. The results are saved in tmp/segments" << std::endl;


				// (5) Find side length Lp used in the occupancy grid for each segment.
				find_lp();

				log.out << "(5) Side lengths are found for each occupancy grid. The results are saved in tmp/lp" << std::endl;


				// (6) Fill in the occupancy grid for each segment.
				const auto fill_all = true;
				fill_in_occupancy_grid(fill_all);

				log.out << "(6) The occupancy grid is projected and filled! All new points are created. The results are saved in tmp/occupancy" << std::endl;


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
				set_default_epsilon_values(value);
			}

		private:
			const Points 			   &m_points;
			const Connected_components &m_cc;
			const Lines                &m_lines;

			std::vector<FT> m_min_dist, m_avg_dist, m_max_dist;
			std::vector<FT> m_eps;
			
			std::vector<FT>  m_lp;
			std::vector<int> m_times;

			Structured_points  m_str_points;
			Structured_labels  m_str_labels;
			Structured_anchors m_str_anchors;

			Points   m_projected;
			Segments m_segments;

			const FT m_tol;
			const FT m_big_value;

			void project() {

				clear_projected();

				Log log;

				assert(!m_lines.empty());
				assert(m_lines.size() == m_cc.size());
				assert(m_projected.empty());

				size_t number_of_lines = 0;
				for (CC_iterator it = m_cc.begin(); it != m_cc.end(); ++it, ++number_of_lines) {
					
					const auto num_points = (*it).second.size();
					for (size_t i = 0; i < num_points; ++i) {

						const auto index = (*it).second[i];
						assert(index >= 0);

						const Point &p = m_points.at(index);
						m_projected[index] = project_onto_a_line(p, m_lines[number_of_lines]);

						assert(!std::isnan(m_projected.at(index).x()) && !std::isnan(m_projected.at(index).y()));
						assert(!std::isinf(m_projected.at(index).x()) && !std::isinf(m_projected.at(index).y()));

						log.out << "index: " << index << ", " << m_projected.at(index) << std::endl;
					}
					log.out << std::endl;
				}
				log.save("tmp/projected");
			}

			void clear_projected() {

				m_projected.clear();
			}

			// Improve this function!
			Point project_onto_a_line(const Point &p, const Line &line) const {

				const auto a = line.point(0);
				const auto b = line.point(1);

				const auto projected = a + dot_product(p - a, b - a) / dot_product(b - a, b - a) * (b - a);
				
				if (std::isnan(projected.x()) || std::isnan(projected.y())) return line.projection(p);
				else return projected;
			}

			void compute_distances() {
				
				clear_distances();

				Log log;

				assert(!m_cc.empty());
				assert(!m_points.empty());

				assert(m_min_dist.size() == m_cc.size());
				assert(m_avg_dist.size() == m_cc.size());
				assert(m_max_dist.size() == m_cc.size());

				assert(!m_projected.empty());
				assert(m_projected.size() == m_points.size());

				size_t number_of_components = 0;

				for (CC_iterator it = m_cc.begin(); it != m_cc.end(); ++it, ++number_of_components) {
					const auto num_points = (*it).second.size();

					m_min_dist[number_of_components] =  m_big_value; 
					m_avg_dist[number_of_components] =      FT(0); 
					m_max_dist[number_of_components] = -m_big_value;

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

					assert(m_max_dist[number_of_components] <  m_big_value);
					assert(m_min_dist[number_of_components] > -m_big_value);

					log.out << "  min: " << m_min_dist[number_of_components] << 
					           "; avg: " << m_avg_dist[number_of_components] <<
					           "; max: " << m_max_dist[number_of_components] << std::endl;
				}
				log.save("tmp/distances");
			}

			void clear_distances() {

				m_min_dist.clear(); m_avg_dist.clear(); m_max_dist.clear();
				m_min_dist.resize(m_cc.size(), FT(0)); m_avg_dist.resize(m_cc.size(), FT(0)); m_max_dist.resize(m_cc.size(), FT(0));
			}

			void compute_epsilon_values() {

				set_default_epsilon_values(FT(0));

				Log log;

				assert(!m_eps.empty());
				assert(m_eps.size() == m_cc.size());
				assert(!m_max_dist.empty());
				assert(m_max_dist.size() == m_eps.size());

				const auto number_of_components = m_max_dist.size();
				const FT default_eps = 0.025;

				for (size_t i = 0; i < number_of_components; ++i) {

					if (m_max_dist[i] < m_tol) m_eps[i] = default_eps;
					else m_eps[i] = m_max_dist[i];

					log.out << "eps: " << m_eps[i] << std::endl;
				}
				log.save("tmp/epsilons");
			}

			void set_default_epsilon_values(const FT value) {

				m_eps.clear();
				m_eps.resize(m_cc.size(), value);
			}

			// Improve this function.
			void find_segments() {

				clear_segments();

				Log log;

				assert(!m_segments.empty());
				assert(m_segments.size() == m_lines.size() && m_segments.size() == m_cc.size());
				assert(!m_projected.empty());

				size_t number_of_segments = 0;
				for (CC_iterator it = m_cc.begin(); it != m_cc.end(); ++it, ++number_of_segments) {

					const auto num_points = (*it).second.size();

					auto minx =  m_big_value, miny =  m_big_value;
					auto maxx = -m_big_value, maxy = -m_big_value;				

					for (size_t i = 0; i < num_points; ++i) {
						const auto index = (*it).second[i];

						const Point &p = m_projected.at(index);

						minx = CGAL::min(minx, p.x());
						maxx = CGAL::max(maxx, p.x());

						miny = CGAL::min(miny, p.y());
						maxy = CGAL::max(maxy, p.y());
					}
					m_segments[number_of_segments] = Segment(Point(minx, miny), Point(maxx, maxy));

					auto v1 = m_lines[number_of_segments].to_vector();
					auto v2 = m_segments[number_of_segments].to_vector();

					// Rotate segments if needed.
					if ((v1.y() < FT(0) && v2.y() >= FT(0) && std::fabs(v1.y() - v2.y()) > m_tol) ||
						(v2.y() < FT(0) && v1.y() >= FT(0) && std::fabs(v1.y() - v2.y()) > m_tol) ||
						(v1.x() < FT(0) && v2.x() >= FT(0) && std::fabs(v1.x() - v2.x()) > m_tol) ||
						(v2.x() < FT(0) && v1.x() >= FT(0) && std::fabs(v1.x() - v2.x()) > m_tol)) {

						m_segments[number_of_segments] = Segment(Point(minx, maxy), Point(maxx, miny));

						log.out << "seg: " << m_segments[number_of_segments].source() << " ----- " << m_segments[number_of_segments].target() << std::endl;
						continue;
					}
					log.out << "seg: " << m_segments[number_of_segments].source() << " ----- " << m_segments[number_of_segments].target() << std::endl;
				}

				assert(number_of_segments == m_lines.size());
				log.save("tmp/segments");
			}

			void clear_segments(){

				m_segments.clear();
				m_segments.resize(m_lines.size());
			}

			void find_lp() {

				clear_lp();

				Log log;

				assert(!m_lp.empty());
				assert(!m_times.empty());

				assert(m_lp.size() == m_cc.size() && m_lp.size() == m_eps.size());
				assert(m_lp.size() == m_segments.size());

				assert(m_times.size() == m_lp.size());

				const auto number_of_segments = m_segments.size();
				for (size_t i = 0; i < number_of_segments; ++i) {

					const auto seg_length = CGAL::sqrt(m_segments[i].squared_length());
					
					const auto upper_bound = CGAL::sqrt(FT(2)) * m_eps[i];

					const auto initial = upper_bound / FT(2);

					const auto times = std::floor(seg_length / initial);

					const auto rest = seg_length - times * initial;

					const auto extra = rest / times;

					m_lp[i]    = initial + extra;
					m_times[i] = static_cast<int>(times);

					assert(m_lp[i] < upper_bound);
					assert(std::fabs(seg_length - m_lp[i] * m_times[i]) < m_tol);

					log.out << "lp: " << m_lp[i] << "; times: " << m_times[i] << std::endl;
				}
				log.save("tmp/lp");
			}

			void clear_lp() {

				m_lp.clear();
				m_lp.resize(m_cc.size());

				m_times.clear();
				m_times.resize(m_cc.size());
			}

			void fill_in_occupancy_grid(const bool) {

			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_STRUCTURING_2_H