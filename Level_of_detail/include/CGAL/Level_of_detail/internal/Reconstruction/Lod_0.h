#ifndef CGAL_LEVEL_OF_DETAIL_LOD_0_H
#define CGAL_LEVEL_OF_DETAIL_LOD_0_H

// STL includes.
#include <vector>
#include <string>

// LOD includes.
#include <CGAL/Level_of_detail/internal/Reconstruction/Lod_0_builder.h>

namespace CGAL {

	namespace Level_of_detail {

        template<class InputKernel, class InputBuilding, class OutputMesh>
        class Lod_0 {
        
        public:
            using Kernel   = InputKernel;
            using Building = InputBuilding;
            using Mesh     = OutputMesh;

            using FT      = typename Kernel::FT;
            using Point_3 = typename Kernel::Point_3;

            using Points_3 = std::vector<Point_3>;

            using Ground_face = Points_3;
            using Wall_faces  = std::vector<Points_3>;
            using Roof_faces  = std::vector<Points_3>;

            using HDS = typename Mesh::HalfedgeDS;

            template<class Buildings, class Ground>
            void reconstruct(const Buildings &buildings, const Ground &ground) {
                
                clear();
                reconstruct_lod(buildings, ground);
            }

            inline std::string name() const {
                return "LOD_0";
            }

            inline size_t number_of_faces() const {
                return 1 + m_wall_faces.size() + m_roof_faces.size();
            }

            size_t number_of_vertices() const {
                size_t num_vertices = m_ground_face.size();
                
                for (auto wall_face = m_wall_faces.begin(); wall_face != m_wall_faces.end(); ++wall_face)
                    num_vertices += wall_face->size();

                for (auto roof_face = m_roof_faces.begin(); roof_face != m_roof_faces.end(); ++roof_face)
                    num_vertices += roof_face->size();

                return num_vertices;
            }

            inline const Ground_face &ground_face() const {
                return m_ground_face;
            }

            inline const Wall_faces &wall_faces() const {
                return m_wall_faces;
            }

            inline const Roof_faces &roof_faces() const {
                return m_roof_faces;
            }

            inline const Mesh &mesh() const {
                return m_mesh;
            }

            void clear() {
                m_ground_face.clear();
                m_wall_faces.clear();
                m_roof_faces.clear();
            }

        private:
            Ground_face m_ground_face;
            Wall_faces  m_wall_faces;
            Roof_faces  m_roof_faces;

            Mesh m_mesh;

            template<class Buildings, class Ground>
            void reconstruct_lod(const Buildings &buildings, const Ground &ground) {
                using Mesh_builder = Lod_0_builder<Kernel, HDS, Building, Buildings, Ground, Ground_face, Wall_faces, Roof_faces>;

                Mesh_builder mesh_builder(buildings, ground, m_ground_face, m_wall_faces, m_roof_faces);
                m_mesh.delegate(mesh_builder);
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_LOD_0_H
