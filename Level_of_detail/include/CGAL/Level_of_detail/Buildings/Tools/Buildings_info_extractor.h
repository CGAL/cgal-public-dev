#ifndef CGAL_LEVEL_OF_DETAIL_BUILDINGS_INFO_EXTRACTOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDINGS_INFO_EXTRACTOR_H

namespace CGAL {

	namespace Level_of_detail {

        template<class InputKernel>
        class Buildings_info_extractor {

        public:
            using Kernel    = InputKernel;
            using Segment_2 = typename Kernel::Segment_2;

            template<class Buildings, class Segments>
            void get_wall_segments(const Buildings &buildings, Segments &segments) const {
                segments.clear();

                for (auto cb_it = buildings.begin(); cb_it != buildings.end(); ++cb_it) {
                    const auto &building = *cb_it;

                    for (auto edge = building.floor_edges().begin(); edge != building.floor_edges().end(); ++edge)
                        segments.push_back(*edge);
                }
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDINGS_INFO_EXTRACTOR_H