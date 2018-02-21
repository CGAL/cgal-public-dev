#ifndef CGAL_LEVEL_OF_DETAIL_CONTAINER_VERTEX_H
#define CGAL_LEVEL_OF_DETAIL_CONTAINER_VERTEX_H

// STL includes.
#include <map>
#include <list>
#include <cmath>
#include <vector>
#include <cassert>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
        class Level_of_detail_container_face;

        template<class KernelTraits>
        class Level_of_detail_container_vertex {

        public:
            typedef KernelTraits Kernel;

            using FT    = typename Kernel::FT;
            using Face  = CGAL::LOD::Level_of_detail_container_face<Kernel>;
            using Faces = std::vector<Face>;

            Level_of_detail_container_vertex(const FT x, const FT y) : m_x(x), m_y(y) { }

            // Coordinates of the vertex.
            FT x() const {
                return m_x;
            }

            FT y() const {
                return m_y;
            }

            // Adjacent faces.
            inline const Faces &adjacent_faces() const {
                return m_faces;
            }

            inline Faces &adjacent_faces() {
                return m_faces;
            }

        private:
            FT m_x, m_y;
            Faces m_faces;
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_CONTAINER_VERTEX_H