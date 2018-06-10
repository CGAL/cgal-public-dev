#ifndef CGAL_LEVEL_OF_DETAIL_CONSTRAINED_TRIANGULATION_CREATOR_H
#define CGAL_LEVEL_OF_DETAIL_CONSTRAINED_TRIANGULATION_CREATOR_H

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel>
		class Constrained_triangulation_creator {

		public:
			using Kernel = InputKernel;

            template<class Faces_range, class Point_map, class Triangulation>
            void make_triangulation_with_info(const Faces_range &faces_range, const Point_map &point_map, Triangulation &triangulation) const {


            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_CONSTRAINED_TRIANGULATION_CREATOR_H