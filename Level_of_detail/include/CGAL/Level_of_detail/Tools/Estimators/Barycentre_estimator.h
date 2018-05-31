#ifndef CGAL_LEVEL_OF_DETAIL_BARYCENTRE_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_BARYCENTRE_ESTIMATOR_H

namespace CGAL {

	namespace Level_of_detail {

        template<class InputKernel>
        class Barycentre_estimator {

        public:
            using Kernel = InputKernel;
            
            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;          
            
            template<class Elements, class Point_map>
			void compute_barycentre_2(const Elements &elements, const Point_map &point_map, Point_2 &barycentre) const {

                CGAL_precondition(elements.size() > 0);
                FT x = FT(0), y = FT(0), size = FT(0);

				for (typename Elements::const_iterator element = elements.begin(); element != elements.end(); ++element, size += FT(1)) {
                    const Point_2 &point = get(point_map, *element);

                    x += point.x();
                    y += point.y();
                }

				x /= size;
				y /= size;

				barycentre = Point_2(x, y);
			}
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BARYCENTRE_ESTIMATOR_H