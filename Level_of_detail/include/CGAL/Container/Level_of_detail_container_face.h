#ifndef CGAL_LEVEL_OF_DETAIL_CONTAINER_FACE_H
#define CGAL_LEVEL_OF_DETAIL_CONTAINER_FACE_H

// STL includes.
#include <map>
#include <list>
#include <cmath>
#include <vector>
#include <cassert>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
        class Level_of_detail_container_vertex;

        template<class KernelTraits>
        class Level_of_detail_container_edge;

        template<class KernelTraits>
        class Level_of_detail_container_face {

        public:
            typedef KernelTraits Kernel;
            using FT = typename Kernel::FT;

            using Vertex = CGAL::LOD::Level_of_detail_container_vertex<Kernel>;
            using Edge   = CGAL::LOD::Level_of_detail_container_edge<Kernel>;

            using Vertices = std::vector<Vertex>;
            using Edges    = std::vector<Edge>;

            Level_of_detail_container_face(const Vertices &vertices, const Edges &edges) : m_vertices(vertices), m_edges(edges) { }

            // Face vertices.
            inline const Vertices &vertices() const {
                return m_vertices;
            }

            inline Vertices &vertices() {
                return m_vertices;
            }

            // Face edges.
            inline const Edges &edges() const {
                return m_edges;
            }

            inline Edges &edges() {
                return m_edges;
            }

        private:
            Vertices m_vertices;
            Edges    m_edges;
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_CONTAINER_FACE_H