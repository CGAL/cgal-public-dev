#ifndef CGAL_LEVEL_OF_DETAIL_POINTS_BASED_REGION_GROWING_2_H
#define CGAL_LEVEL_OF_DETAIL_POINTS_BASED_REGION_GROWING_2_H

// STL includes.
#include <map>
#include <list>
#include <cmath>
#include <algorithm>
#include <unordered_set>

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel, class PointIdentifier>
		class Points_based_region_growing_2 {

		public:
		    using Kernel           = InputKernel;
            using Point_identifier = PointIdentifier;

			using FT 	   = typename Kernel::FT;
			using Point_2  = typename Kernel::Point_2;
			using Vector_2 = typename Kernel::Vector_2;

			Points_based_region_growing_2(const FT epsilon, const FT cluster_epsilon, const FT normal_threshold, const FT min_points) :
			m_epsilon(epsilon),
			m_cluster_epsilon(cluster_epsilon),
			m_normal_threshold(normal_threshold),
			m_min_points(min_points)
			{ }

			template<class Elements, class Point_map, class Normal_map, class Range>
			void detect(const Elements &elements, const Point_map &point_map, const Normal_map &normal_map, std::list<Range> &output) const {
				CGAL_precondition(elements.size() > 1);

			}

		private:
			const FT m_epsilon;
			const FT m_cluster_epsilon;
			const FT m_normal_threshold;
			const FT m_min_points;
		};

	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_POINTS_BASED_REGION_GROWING_2_H