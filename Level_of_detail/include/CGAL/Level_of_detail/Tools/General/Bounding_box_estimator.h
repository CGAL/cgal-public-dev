#ifndef CGAL_LEVEL_OF_DETAIL_BOUNDING_BOX_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_BOUNDING_BOX_ESTIMATOR_H

namespace CGAL {

	namespace Level_of_detail {

        template<class InputKernel>
        class Bounding_box_estimator {

        public:
            using Kernel = InputKernel;
            
            using FT      = typename Kernel::FT;
            using Point_3 = typename Kernel::Point_3;
            using Plane_3 = typename Kernel::Plane_3;

            Bounding_box_estimator() : m_big_value(FT(100000000000000)) { }

            template<class Elements, class Point_map, class Bounding_box_range_3>
            void compute_horizontal_bounding_box_3(const Elements &elements, const Point_map &point_map, const Plane_3 &plane, Bounding_box_range_3 &bounding_box) const {

				FT minx =  m_big_value, miny =  m_big_value;
                FT maxx = -m_big_value, maxy = -m_big_value;
                
                FT z = FT(0); size_t i = 0;
				for (typename Elements::const_iterator element = elements.begin(); element != elements.end(); ++element, ++i) {
					
                    const Point_3 &point    = get(point_map, *element);
                    const Point_3 projected = plane.projection(point);

                    minx = CGAL::min(minx, projected.x());
					miny = CGAL::min(miny, projected.y());

					maxx = CGAL::max(maxx, projected.x());
					maxy = CGAL::max(maxy, projected.y());

                    z += projected.z();
                }
                z /= static_cast<FT>(i);

                bounding_box.clear();
                bounding_box.push_back(Point_3(minx, miny, z));
                bounding_box.push_back(Point_3(maxx, miny, z));
                bounding_box.push_back(Point_3(maxx, maxy, z));
                bounding_box.push_back(Point_3(minx, maxy, z));
            }

        private:
            const FT m_big_value;
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BOUNDING_BOX_ESTIMATOR_H