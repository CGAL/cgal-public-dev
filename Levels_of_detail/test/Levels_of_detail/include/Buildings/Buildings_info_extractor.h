#ifndef CGAL_LEVEL_OF_DETAIL_BUILDINGS_INFO_EXTRACTOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDINGS_INFO_EXTRACTOR_H

namespace CGAL {

	namespace Level_of_detail {

        template<class InputKernel>
        class Buildings_info_extractor {

        public:
            using Kernel = InputKernel;
            
            using Point_2   = typename Kernel::Point_2;
            using Segment_2 = typename Kernel::Segment_2;
            using Point_3   = typename Kernel::Point_3;

            template<class Buildings, class Segments>
            void get_wall_segments(const Buildings &buildings, Segments &segments) const {
                segments.clear();

                for (auto cb_it = buildings.begin(); cb_it != buildings.end(); ++cb_it) {
                    const auto &building = *cb_it;

                    for (auto edge = building.floor_edges().begin(); edge != building.floor_edges().end(); ++edge)
                        segments.push_back(*edge);
                }
            }

            template<class Buildings, class Triangle>
            void get_flat_roof_triangles(const Buildings &buildings, std::list<Triangle> &triangles) const {
                triangles.clear();

                for (auto cb_it = buildings.begin(); cb_it != buildings.end(); ++cb_it) {
                    const auto &building = *cb_it;

                    for (auto face = building.floor_faces().begin(); face != building.floor_faces().end(); ++face) {

                        const Point_2 &p1 = face->vertex(0);
                        const Point_2 &p2 = face->vertex(1);
                        const Point_2 &p3 = face->vertex(2);

                        Triangle triangle;
                        triangle.push_back(Point_3(p1.x(), p1.y(), building.local_ground_height() + building.height()));
                        triangle.push_back(Point_3(p2.x(), p2.y(), building.local_ground_height() + building.height()));
                        triangle.push_back(Point_3(p3.x(), p3.y(), building.local_ground_height() + building.height()));
                        
                        triangles.push_back(triangle);
                    }
                }
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDINGS_INFO_EXTRACTOR_H