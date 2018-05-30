#ifndef CGAL_LEVEL_OF_DETAIL_GRID_BASED_FILTERING_H
#define CGAL_LEVEL_OF_DETAIL_GRID_BASED_FILTERING_H

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel, class PointIdentifier>
		class Grid_based_filtering {

        public:
            using Kernel           = InputKernel;
            using Point_identifier = PointIdentifier;

            template<class Elements, class Point_map, class Output>
            void apply(const Elements &elements, const Point_map &point_map, Output &output) const {
            
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_GRID_BASED_FILTERING_H