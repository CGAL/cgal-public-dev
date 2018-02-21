#ifndef CGAL_LEVEL_OF_DETAIL_CONTAINER_H
#define CGAL_LEVEL_OF_DETAIL_CONTAINER_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <cassert>

// CGAL includes.
#include <CGAL/Polygon_2.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Container/Level_of_detail_container_edge.h>
#include <CGAL/Container/Level_of_detail_container_face.h>
#include <CGAL/Container/Level_of_detail_container_halfedge.h>
#include <CGAL/Container/Level_of_detail_container_vertex.h>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
        struct MyPolygon {
        
        public:
            typedef KernelTraits Kernel;
            typedef CGAL::Polygon_2<Kernel> Polygon;
            
            using FT = typename Kernel::FT;

            Polygon polygon;
            FT inside = -FT(1);
        };

        template<class KernelTraits>
        class Level_of_detail_container {

        public:
            typedef KernelTraits Kernel;

            using Edge     = CGAL::LOD::Level_of_detail_container_edge<Kernel>;
            using Face     = CGAL::LOD::Level_of_detail_container_face<Kernel>;
            using Halfedge = CGAL::LOD::Level_of_detail_container_halfedge<Kernel>;
            using Vertex   = CGAL::LOD::Level_of_detail_container_vertex<Kernel>;

            using Polygon = CGAL::LOD::MyPolygon<Kernel>;

            using Polygons = std::vector<Polygon>;
            using Faces    = std::vector<Face>;

            Level_of_detail_container() { }

            inline size_t number_of_polygons() const {
                return m_polygons.size();
            }

            inline const Polygons &polygons() const {
                return m_polygons;
            }

            inline Polygons &polygons() {
                return m_polygons;
            }

        private:
            Polygons m_polygons;
            Faces    m_faces;
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_CONTAINER_H