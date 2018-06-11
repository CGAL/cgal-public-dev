#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_H

// STL includes.
#include <list>

namespace CGAL {

	namespace Level_of_detail {

        template<class InputKernel>
        class Building {
        
        public:
            using Kernel = InputKernel;
            
            using FT         = typename Kernel::FT;
            using Point_2    = typename Kernel::Point_2;
            using Segment_2  = typename Kernel::Segment_2;
            using Triangle_2 = typename Kernel::Triangle_2;

            using Floor_edges = std::list<Segment_2>;
            using Floor_faces = std::list<Triangle_2>;

            inline void add_floor_face(const Point_2 &p1, const Point_2 &p2, const Point_2 &p3) {
                m_floor_faces.push_back(Triangle_2(p1, p2, p3));
            }

        private:
            Floor_edges m_floor_edges;
            Floor_faces m_floor_faces;
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_H