#ifndef CGAL_LEVEL_OF_DETAIL_END_POINTS_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_END_POINTS_ESTIMATOR_H

// LOD includes.
#include <CGAL/Level_of_detail/Tools/Estimations/Barycentre_estimator.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

        template<class InputKernel>
        class End_points_estimator {

        public:
            using Kernel = InputKernel;
            
            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;

            using Barycentre_estimator = LOD::Barycentre_estimator<Kernel>;
            typename Kernel::Compute_squared_distance_2 squared_distance_2;

            template<class Elements, class Point_map>
            void find_end_points_wrt_barycentre_2(const Elements &elements, const Point_map &point_map, Point_2 &a, Point_2 &b) const {
                CGAL_precondition(elements.size() > 1);

                // Find barycentre.
                Point_2 barycentre;
                const Barycentre_estimator barycentre_estimator;
                barycentre_estimator.compute_barycentre_2(elements, point_map, barycentre);

                // Find the first end point.
                find_furthest_point_2(elements, point_map, barycentre, a);

                // Find the second end point.
                find_furthest_point_2(elements, point_map, a, b);
            }

            template<class Elements, class Point_map>
            void find_furthest_point_2(const Elements &elements, const Point_map &point_map, const Point_2 &reference_point, Point_2 &resulting_point) const {
    
                CGAL_precondition(elements.size() > 0);
                using Const_elements_iterator = typename Elements::const_iterator;

                FT max_distance = FT(0);
                Const_elements_iterator max_it;

                for (Const_elements_iterator ce_it = elements.begin(); ce_it != elements.end(); ++ce_it) {
                    const Point_2 &point = get(point_map, *ce_it);

                    const FT distance = squared_distance_2(point, reference_point);
                    if (distance > max_distance) {
                        
                        max_distance = distance;
                        max_it       = ce_it;
                    }
                }
                resulting_point = get(point_map, *max_it);
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_END_POINTS_ESTIMATOR_H