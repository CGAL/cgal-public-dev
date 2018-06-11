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
                using Const_elements_iterator = typename Elements::const_iterator;
                
                FT x = FT(0), y = FT(0), size = FT(0);

				for (Const_elements_iterator ce_it = elements.begin(); ce_it != elements.end(); ++ce_it, size += FT(1)) {
                    const Point_2 &point = get(point_map, *ce_it);

                    x += point.x();
                    y += point.y();
                }

				x /= size;
				y /= size;

				barycentre = Point_2(x, y);
			}

            template<class Face_handle>
            void compute_triangulation_face_barycentre_2(const Face_handle &face_handle, Point_2 &barycentre) const {

				FT x = FT(0), y = FT(0);
				for (size_t i = 0; i < 3; ++i) {

					x += face_handle->vertex(i)->point().x();
					y += face_handle->vertex(i)->point().y();
				}

				x /= FT(3);
				y /= FT(3);

				barycentre = Point_2(x, y);
			}
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BARYCENTRE_ESTIMATOR_H