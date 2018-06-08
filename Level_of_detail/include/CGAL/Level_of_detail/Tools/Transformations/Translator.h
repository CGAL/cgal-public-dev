#ifndef CGAL_LEVEL_OF_DETAIL_TRANSLATOR_H
#define CGAL_LEVEL_OF_DETAIL_TRANSLATOR_H

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel>
		class Translator {

		public:
			using Kernel = InputKernel;
            
            using FT        = typename Kernel::FT;
            using Point_2   = typename Kernel::Point_2;
            using Segment_2 = typename Kernel::Segment_2;

            template<class Segment_map, class Elements>
            void translate_segments_2(const Point_2 &translation, const Segment_map &segment_map, Elements &elements) const {
                
                CGAL_precondition(elements.size() > 0);
                using Elements_iterator = typename Elements::iterator;

                Point_2 new_source, new_target;
                Segment_2 new_segment;

                for (Elements_iterator ce_it = elements.begin(); ce_it != elements.end(); ++ce_it) {
                    const Segment_2 &segment = get(segment_map, *ce_it);
                    
                    const Point_2 &source = segment.source();
                    const Point_2 &target = segment.target();

                    const FT x1 = source.x() - translation.x();
                    const FT y1 = source.y() - translation.y();

                    const FT x2 = target.x() - translation.x();
                    const FT y2 = target.y() - translation.y();

                    new_source = Point_2(x1, y1);
                    new_target = Point_2(x2, y2);

                    new_segment = Segment_2(new_source, new_target);
                    put(segment_map, *ce_it, new_segment);
                }
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_TRANSLATOR_H