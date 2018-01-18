#ifndef CGAL_LEVEL_OF_DETAIL_CLUTTER_FILTERING_H
#define CGAL_LEVEL_OF_DETAIL_CLUTTER_FILTERING_H

#if defined(WIN32) || defined(_WIN32) 
#define PS "\\"
#define PN "\r\n"
#else 
#define PS "/" 
#define PN "\n"
#endif

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>
#include <unordered_set>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Simple_cartesian.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Utils/Level_of_detail_utils_simple.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class BoundaryData, class ProjectedPoints>
		class Level_of_detail_clutter_filtering {

		public:
			typedef KernelTraits 	Kernel;
			typedef BoundaryData 	Boundary_data;
			typedef ProjectedPoints Projected_points;

			typedef typename Kernel::FT 	 FT;
			typedef typename Kernel::Point_2 Point_2;
			typedef typename Kernel::Line_2  Line_2;

			typedef std::pair<int, Point_2> 					 				Projected_point;
			typedef typename CGAL::Second_of_pair_property_map<Projected_point> Point_map;

			typedef CGAL::Search_traits_2<Kernel>                       					 Search_traits_2;
			typedef CGAL::Search_traits_adapter<Projected_point, Point_map, Search_traits_2> Search_traits;
			typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> 					     Neighbor_search;
			typedef CGAL::Fuzzy_sphere<Search_traits>                    					 Fuzzy_circle;
			typedef CGAL::Kd_tree<Search_traits>					                         Fuzzy_tree;

			typedef typename Projected_points::const_iterator Point_iterator;
			typedef typename Neighbor_search::iterator  	  Neighbour_iterator;

			typename Kernel::Compute_squared_distance_2 squared_distance;
			typename Kernel::Compute_scalar_product_3   dot_product;

            typedef Level_of_detail_utils_simple<Kernel> Simple_utils;
			using Log = CGAL::LOD::Mylog;

            Level_of_detail_clutter_filtering() : m_scale(-FT(1)), m_mean(-FT(1)), m_silent(false) { }

            void set_scale(const FT scale) {
                
                assert(scale > FT(0));
                m_scale = scale;
            }

            void set_mean(const FT mean) {

                assert(mean > FT(0));
                m_mean = mean;
            }

            void make_silent(const bool silent) {
                m_silent = silent;
            }

            int filter(Boundary_data &, Projected_points &boundary_clutter_projected) const {
                assert(!boundary_clutter_projected.empty());

                Fuzzy_tree tree;
				create_tree(tree, boundary_clutter_projected);

                Projected_points filtered_points;
				filter_points(filtered_points, tree, boundary_clutter_projected);
				boundary_clutter_projected = filtered_points;

                const int number_of_new_points = static_cast<int>(filtered_points.size());
				assert(number_of_new_points >= 0);

                if (!m_silent) {
                    Log exporter; exporter.export_projected_points_as_xyz("tmp" + std::string(PS) + "filtered_clutter", boundary_clutter_projected, "unused path");
                }

                return number_of_new_points;
            }

        private:
            FT m_scale;
            FT m_mean;

            bool m_silent;
            Simple_utils m_simple_utils;

            void create_tree(Fuzzy_tree &tree, const Projected_points &boundary_clutter_projected) const {

				assert(!boundary_clutter_projected.empty());
				tree.insert(boundary_clutter_projected.begin(), boundary_clutter_projected.end());
			}

            void filter_points(Projected_points &filtered_points, const Fuzzy_tree &tree, const Projected_points &boundary_clutter_projected) const {

                assert(!boundary_clutter_projected.empty());
				filtered_points.clear();

				for (Point_iterator pit = boundary_clutter_projected.begin(); pit != boundary_clutter_projected.end(); ++pit) {
					
					const Projected_point &projected = *pit;
					handle_projected_point(filtered_points, tree, projected);
				}
            }

            void handle_projected_point(Projected_points &filtered_points, const Fuzzy_tree &tree, const Projected_point &query) const {
				
                Projected_points neighbours;
				find_nearest_neighbours(neighbours, tree, query);

                const bool preserve_query = apply_filter(neighbours, query);
                if (preserve_query) filtered_points[query.first] = query.second;
            }

            void find_nearest_neighbours(Projected_points &neighbours, const Fuzzy_tree &tree, const Projected_point &query) const {
				neighbours.clear();

				Fuzzy_circle circle;
				compute_circle(query.second, circle);

				tree.search(std::inserter(neighbours, neighbours.end()), circle);
				neighbours[query.first] = query.second;
            }

			void compute_circle(const Point_2 &centre, Fuzzy_circle &circle) const {

				assert(m_scale > FT(0));
				circle = Fuzzy_circle(centre, m_scale);
			}

            bool apply_filter(const Projected_points &neighbours, const Projected_point &query) const {
                
                assert(!neighbours.empty());
                bool preserve_query = false;

                Line_2 line;
                fit_line_to_points(line, neighbours);

                std::vector<FT> distances(neighbours.size());
                
                const int num_good_points = compute_distances(distances, line, neighbours);
                const FT query_dist       = compute_distance(line, query.second);
                const FT local_mean       = compute_mean(distances);

                if (should_query_be_preserved(num_good_points, query_dist, local_mean)) preserve_query = true;
                return preserve_query;
            }

			void fit_line_to_points(Line_2 &line, const Projected_points &points) const {

      			using Local_Kernel = CGAL::Simple_cartesian<double>;
				using Point_2ft    = Local_Kernel::Point_2;
				using Line_2ft     = Local_Kernel::Line_2;

				const size_t num_points = points.size();
				std::vector<Point_2ft> tmp_points(num_points);

				size_t count = 0;
				for (Point_iterator pit = points.begin(); pit != points.end(); ++pit, ++count) {
				
					const Point_2 &p = (*pit).second;

					const double x = CGAL::to_double(p.x());
					const double y = CGAL::to_double(p.y());

					tmp_points[count] = Point_2ft(x, y);
				}
				assert(num_points == count);

				Line_2ft tmp_line;
				CGAL::linear_least_squares_fitting_2(tmp_points.begin(), tmp_points.end(), tmp_line, CGAL::Dimension_tag<0>());
				line = Line_2(static_cast<FT>(tmp_line.a()), static_cast<FT>(tmp_line.b()), static_cast<FT>(tmp_line.c()));
			}

            int compute_distances(std::vector<FT> &distances, const Line_2 &line, const Projected_points &points) const {
                assert(distances.size() == points.size());

                int num_good_points = 0;
                size_t count = 0;

                for (Point_iterator pit = points.begin(); pit != points.end(); ++pit, ++count) {
                    const Point_2 &point = (*pit).second;

                    distances[count] = compute_distance(line, point);
                    if (is_within_tolerance(distances[count])) ++num_good_points;
                }

                assert(points.size() == count);
                return num_good_points;
            }

            FT compute_distance(const Line_2 &line, const Point_2 &point) const {
                
                const Point_2 projected = m_simple_utils.project_onto_line(line, point);
                return static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance(point, projected))));
            }

            FT compute_mean(const std::vector<FT> &values) const {
                assert(!values.empty());

                FT mean = FT(0);
                for (size_t i = 0; i < values.size(); ++i) mean += values[i];
                mean /= static_cast<FT>(values.size());

                return mean;
            }

            bool is_within_tolerance(const FT value) const {
                
                if (value < m_mean) return true;
                return false;
            }

            bool should_query_be_preserved(const int num_good_points, const FT query_dist, const FT local_mean) const {
                bool preserve_query = true;

                if (!is_within_tolerance(local_mean)) preserve_query = preserved_state_on();
                if (num_good_points < 3)              preserve_query = preserved_state_on();
                if (!is_within_tolerance(query_dist)) preserve_query = preserved_state_on();

                return preserve_query;
            }

            bool preserved_state_on() const {
                return false;
            }

            bool preserved_state_off() const {
                return true;
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_CLUTTER_FILTERING_H