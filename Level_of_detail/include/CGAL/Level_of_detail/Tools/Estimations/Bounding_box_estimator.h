#ifndef CGAL_LEVEL_OF_DETAIL_BOUNDING_BOX_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_BOUNDING_BOX_ESTIMATOR_H

namespace CGAL {

	namespace Level_of_detail {

        template<class InputKernel>
        class Bounding_box_estimator {

        public:
            using Kernel = InputKernel;
            
            using FT        = typename Kernel::FT;
            using Point_2   = typename Kernel::Point_2;
            using Segment_2 = typename Kernel::Segment_2;
            using Point_3   = typename Kernel::Point_3;
            using Plane_3   = typename Kernel::Plane_3;

            Bounding_box_estimator() : 
            m_big_value(FT(100000000000000)) 
            { }

            template<class Elements, class Segment_map, class Bounding_box_range_2>
            void compute_bounding_box_2(const Elements &elements, const Segment_map &segment_map, Bounding_box_range_2 &bounding_box) const {
                
                CGAL_precondition(elements.size() > 0);
                using Const_elements_iterator = typename Elements::const_iterator;

                FT minx =  m_big_value, miny =  m_big_value;
                FT maxx = -m_big_value, maxy = -m_big_value;

				for (Const_elements_iterator ce_it = elements.begin(); ce_it != elements.end(); ++ce_it) {
                    const Segment_2 &segment = get(segment_map, *ce_it);
                    
                    const Point_2 &source = segment.source();
                    const Point_2 &target = segment.target();

                    minx = CGAL::min(minx, source.x()); minx = CGAL::min(minx, target.x());
                    miny = CGAL::min(miny, source.y()); miny = CGAL::min(miny, target.y());

                    maxx = CGAL::max(maxx, source.x()); maxx = CGAL::max(maxx, target.x());
                    maxy = CGAL::max(maxy, source.y()); maxy = CGAL::max(maxy, target.y());
                }

                bounding_box.clear();
                bounding_box.push_back(Point_2(minx, miny));
                bounding_box.push_back(Point_2(maxx, miny));
                bounding_box.push_back(Point_2(maxx, maxy));
                bounding_box.push_back(Point_2(minx, maxy));
            }

            template<class Elements, class Point_map, class Bounding_box_range_3>
            void compute_bounding_box_3(const Elements &elements, const Point_map &point_map, const Plane_3 &plane, Bounding_box_range_3 &bounding_box) const {
                
                CGAL_precondition(elements.size() > 0);
                using Const_elements_iterator = typename Elements::const_iterator;

				FT minx =  m_big_value, miny =  m_big_value;
                FT maxx = -m_big_value, maxy = -m_big_value;
                
                FT z = FT(0), size = FT(0);
				for (Const_elements_iterator ce_it = elements.begin(); ce_it != elements.end(); ++ce_it, size += FT(1)) {
					
                    const Point_3 &point    = get(point_map, *ce_it);
                    const Point_3 projected = plane.projection(point);

                    minx = CGAL::min(minx, projected.x());
					miny = CGAL::min(miny, projected.y());

					maxx = CGAL::max(maxx, projected.x());
					maxy = CGAL::max(maxy, projected.y());

                    z += projected.z();
                }
                z /= size;

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