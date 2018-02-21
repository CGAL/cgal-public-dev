#ifndef CGAL_LEVEL_OF_DETAIL_CONTAINER_EDGE_H
#define CGAL_LEVEL_OF_DETAIL_CONTAINER_EDGE_H

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
        class Level_of_detail_container_halfedge;

        template<class KernelTraits>
        class Level_of_detail_container_edge {

        public:
            typedef KernelTraits Kernel;
            using FT = typename Kernel::FT;
            
            using Vertex   = CGAL::LOD::Level_of_detail_container_vertex<Kernel>;
            using Halfedge = CGAL::LOD::Level_of_detail_container_halfedge<Kernel>;

            Level_of_detail_container_edge(const Vertex &source, const Vertex &target) : m_source(source), m_target(target) { }

            // Vertices of the edge.
            inline const Vertex &source() const {
                return m_source;
            }

            inline const Vertex &target() const {
                return m_target;
            }

            // Halfedge.
            inline const Halfedge &halfedge() const {
                return m_halfedge;
            }

            inline Halfedge &halfedge() {
                return m_halfedge;
            }

            // Opposite halfedge.
            inline const Halfedge &opposite() const {
                return m_opposite;
            }

            inline Halfedge &opposite() {
                return m_opposite;
            }

        private:
            Vertex m_source;
            Vertex m_target;

            Halfedge m_halfedge;
            Halfedge m_opposite;
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_CONTAINER_EDGE_H