#ifndef CGAL_LEVEL_OF_DETAIL_CONTAINER_HALFEDGE_H
#define CGAL_LEVEL_OF_DETAIL_CONTAINER_HALFEDGE_H

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
        class Level_of_detail_container_halfedge {

        public:
            typedef KernelTraits Kernel;
            using FT = typename Kernel::FT;

            using Face = CGAL::LOD::Level_of_detail_container_face<Kernel>;

            Level_of_detail_container_halfedge(const Face &adjacent_face) : m_adjacent_face(adjacent_face) { }

            inline const Face &adjacent_face() const {
                return m_adjacent_face;
            }

        private:
            Face m_adjacent_face;
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_CONTAINER_HALFEDGE_H