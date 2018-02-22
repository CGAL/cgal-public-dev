#ifndef CGAL_LEVEL_OF_DETAIL_CONTAINER_H
#define CGAL_LEVEL_OF_DETAIL_CONTAINER_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <cassert>

// CGAL includes.
#include <CGAL/IO/Color.h>
#include <CGAL/Polygon_2.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Container/Level_of_detail_container_edge.h>
#include <CGAL/Container/Level_of_detail_container_face.h>
#include <CGAL/Container/Level_of_detail_container_vertex.h>
#include <CGAL/Container/Level_of_detail_container_halfedge.h>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
        struct MyContainer {
        
        public:
            typedef KernelTraits Kernel;
            typedef CGAL::Polygon_2<Kernel> Polygon;
            
            using FT          = typename Kernel::FT;
            using Colour      = CGAL::Color;
            using Constraints = std::vector<bool>;

            Polygon polygon;
            Constraints constraints;

            FT inside    = FT(1);
            Color colour = Color(55, 255, 55);
        };

        template<class KernelTraits>
        class Level_of_detail_container {

        public:
            typedef KernelTraits Kernel;

            using Point_2 = typename Kernel::Point_2;

            using Container  = CGAL::LOD::MyContainer<Kernel>;
            using Containers = std::vector<Container>;

            Level_of_detail_container() { }

            inline size_t number_of_containers() const {
                return m_containers.size();
            }

            inline const Containers &containers() const {
                return m_containers;
            }

            inline Containers &containers() {
                return m_containers;
            }

            int locate(const Point_2 &query) const {
                assert(m_containers.size() > 0);

                for (size_t i = 0; i < m_containers.size(); ++i)
                    if (m_containers[i].polygon.has_on_bounded_side(query))
                        return static_cast<int>(i);

                return -1;
            }

        private:
            Containers m_containers;
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_CONTAINER_H