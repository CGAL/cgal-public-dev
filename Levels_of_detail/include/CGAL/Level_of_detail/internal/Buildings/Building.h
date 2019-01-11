#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_H

// STL includes.
#include <list>
#include <vector>

namespace CGAL {

	namespace Level_of_detail {

        template<class InputKernel, class FloorFaceHandle>
        class Building {
        
        public:
            using Kernel            = InputKernel;
            using Floor_face_handle = FloorFaceHandle;
            
            using FT         = typename Kernel::FT;
            using Point_2    = typename Kernel::Point_2;
            using Segment_2  = typename Kernel::Segment_2;
            using Triangle_2 = typename Kernel::Triangle_2;

            using Edge    = Segment_2;
            using Polygon = Triangle_2;

            using Floor_edges        = std::list<Edge>;
            using Floor_faces        = std::list<Polygon>;
            using Floor_face_handles = std::vector<Floor_face_handle>;
            
            using Const_floor_faces_iterator        = typename Floor_faces::const_iterator;
            using Const_floor_edges_iterator        = typename Floor_edges::const_iterator;
            using Const_floor_face_handles_iterator = typename Floor_face_handles::const_iterator;

            Building() :
            m_index(-1),
            m_is_valid(true),
            m_height(FT(0)),
            m_local_ground_height(FT(0))
            { }

            inline int &index() {
                return m_index;
            }

            inline const int &index() const {
                return m_index;
            }

            inline bool &is_valid() {
                return m_is_valid;
            }

            inline const bool &is_valid() const {
                return m_is_valid;
            }


            // Extra functions.
            inline void add_floor_edge(const Point_2 &p1, const Point_2 &p2) {
                m_floor_edges.push_back(Edge(p1, p2));
            }

            inline void add_floor_face(const Point_2 &p1, const Point_2 &p2, const Point_2 &p3) {
                m_floor_faces.push_back(Polygon(p1, p2, p3));
            }

            inline void add_floor_face_handle(const Floor_face_handle &floor_face_handle) {
                m_floor_face_handles.push_back(floor_face_handle);
            }


            // Heights.
            inline FT &height() {
                return m_height;
            }

            inline const FT &height() const {
                return m_height;
            }

            inline FT &local_ground_height() {
                return m_local_ground_height;
            }

            inline const FT &local_ground_height() const {
                return m_local_ground_height;
            }


            // Floor.
            inline Floor_edges &floor_edges() {
                return m_floor_edges;
            }

            inline const Floor_edges &floor_edges() const {
                return m_floor_edges;
            }

            inline Floor_faces &floor_faces() {
                return m_floor_faces;
            }

            inline const Floor_faces &floor_faces() const {
                return m_floor_faces;
            }

            inline Floor_face_handles &floor_face_handles() {
                return m_floor_face_handles;
            }

            inline const Floor_face_handles &floor_face_handles() const {
                return m_floor_face_handles;
            }

        private:
            int  m_index;
            bool m_is_valid;

            FT m_height;
            FT m_local_ground_height;

            Floor_edges m_floor_edges;
            Floor_faces m_floor_faces;

            Floor_face_handles m_floor_face_handles;       
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_H