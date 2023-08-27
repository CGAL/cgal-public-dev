// Copyright (c) 2023
// INRIA Sophia-Antipolis (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jeffrey Cochran

#ifndef COLLISION_SCENE_3_H
#define COLLISION_SCENE_3_H

#include <CGAL/Continuous_collision_detection_3/internal/AABB_triangle_trajectory_primitive.h>

#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Mesh_index.h>
#include <CGAL/Trajectories.h>
#include <CGAL/Collision_mesh_3.h>

#include <vector>
#include <utility>


namespace CGAL {



    template <class K>
    class Collision_scene {

        public:
        
            template<class Local_index> struct Scene_index;

            typedef          K                              Kernel;
            typedef typename ::CGAL::Mesh_index             Mesh_index;
            typedef typename K::Point_3                     Point;
            typedef typename K::Vector_3                    Vector;
            typedef          Collision_mesh<K>              Mesh;
            typedef typename Mesh::Vertex_index             Vertex_index;
            typedef typename Mesh::Edge_index               Edge_index;
            typedef typename Mesh::Halfedge_index           Halfedge_index;
            typedef typename Mesh::Face_index               Face_index;
            typedef typename Mesh::Face_range               Face_range;
            typedef typename Mesh::Vertex_range             Vertex_range;

            typedef          Scene_index<Vertex_index>                                  Scene_vertex_index;
            typedef          Scene_index<Face_index>                                    Scene_face_index;
            typedef          std::vector<Scene_vertex_index>                            Scene_vertex_range;
            typedef          std::vector<Scene_face_index>                              Scene_face_range;
            typedef          Triangle_trajectory_observer<K, Scene_face_index>          Trajectory;
            typedef          std::vector<Trajectory>                                    Trajectory_range;
            typedef          AABB_Triangle_trajectory_primitive<K, Scene_face_index>    Trajectory_primitive;
            typedef          ::CGAL::AABB_traits<K, Trajectory_primitive>               AABB_traits;
            typedef          AABB_tree<AABB_traits>                                     Tree;
            typedef typename Tree::Primitive_id                                         Primitive_id;


        private:
            Trajectory make_trajectory(const Mesh & mesh, const Scene_face_index& scene_face_index);

            Scene_vertex_range  vertices_;
            Scene_face_range    faces_;
            Trajectory_range    trajectories_;
            std::vector<Mesh>&  meshes_;
            Tree                trajectory_tree_;

        public:

            explicit Collision_scene(std::vector<Mesh> & meshes);

            const Scene_vertex_range& vertices() const;
            const Scene_face_range& faces() const;
            const Tree& tree() const;

            Collision_mesh<K> joined_meshes();

            void color(const Scene_face_index& ti, CGAL::IO::Color c);

            template <class UpdateFunctor>
            void update_state(UpdateFunctor& update_functor, bool update_tree=false);

            void update_tree();
    };


    // ============
    // Constructors
    // ============
    template <class K>
    Collision_scene<K>::Collision_scene(std::vector<Mesh> & meshes) : meshes_{meshes}
    {
        size_t i{0};
        for( const auto & mesh_ : meshes_ )
        {
            Mesh_index current_mesh_index{i};

            vertices_.reserve(vertices_.size() + mesh_.num_vertices());
            for( const auto & current_vertex_index : mesh_.vertices()) { 
                vertices_.push_back( 
                    Scene_vertex_index{
                        current_mesh_index,
                        current_vertex_index
                    } 
                );
            }

            size_t k{faces_.size()};
            faces_.reserve(faces_.size() + mesh_.num_faces());
            trajectories_.reserve(trajectories_.size() + mesh_.num_faces());
            for( const auto & current_face_index : mesh_.faces()) { 
                faces_.push_back( 
                    Scene_face_index{
                        current_mesh_index,
                        current_face_index
                    }
                );
                trajectories_.push_back(
                    make_trajectory(mesh_, faces_[k])
                );
                ++k;
            }

            ++i;
        }

        trajectory_tree_ = Tree(trajectories_.begin(),trajectories_.end());
    };
    

    // ========================
    // Member Class definitions
    // ========================
    template <class K> 
    template <class Local_index>
    struct Collision_scene<K>::Scene_index : std::pair<Mesh_index, Local_index> {

        typedef std::pair<Mesh_index, Local_index> Base;
        typedef std::size_t size_type;

        Scene_index() {}

        Scene_index(Mesh_index mi, Local_index li)
        : Base(mi,li)
        {}

        Scene_index(const Base& index_pair)
        : Base(index_pair)
        {}

        friend std::ostream& operator<<(std::ostream& os, Scene_index const& Scene_index)
        {
            return (os << "(" << Scene_index.first << ", " << Scene_index.second << ")" );
        } 

        bool operator<(Scene_index const& other_Scene_index) const {
            bool less_than = this->first == other_Scene_index.first 
                ? this->second < other_Scene_index.second
                : this->first < other_Scene_index.first;
            return less_than;
        }

        Mesh_index mesh_index() const {
          return this->first;
        }
        
        Local_index local_index() const {
          return this->second;
        }

    };

    // ===============
    // Member routines
    // ===============    
    template <class K>
    const typename Collision_scene<K>::Tree& Collision_scene<K>::tree() const
    {
        return this->trajectory_tree_;
    }

    template <class K>
    const typename Collision_scene<K>::Scene_vertex_range& Collision_scene<K>::vertices() const
    {
        return this->vertices_;
    }

    template <class K>
    const typename Collision_scene<K>::Scene_face_range& Collision_scene<K>::faces() const
    {
        return this->faces_;
    }

    template <class K>
    auto Collision_scene<K>::make_trajectory( 
        const Mesh& mesh, 
        const Scene_face_index& scene_face_index 
    ) -> Trajectory 
    {

        size_t j{0};
        std::vector<const Point*> trajectory_points(6);

        for( const auto & vert_index : mesh.vertices_around_face(mesh.halfedge(scene_face_index.local_index())) )
        {
            trajectory_points.at(j) = & mesh.point(vert_index);
            trajectory_points.at(j+1) = & mesh.next_point(vert_index);

            j += 2;
        }

        return Trajectory(
            trajectory_points[0], // v0_current
            trajectory_points[2], // v0_next
            trajectory_points[4], // v1_current
            trajectory_points[1], // v1_next
            trajectory_points[3], // v2_current
            trajectory_points[5], // v2_next
            scene_face_index
        );
    };

    template <class K>
    Collision_mesh<K> Collision_scene<K>::joined_meshes()
    {
        Mesh joined_mesh{meshes_[0]};

        auto begin = meshes_.begin();
        auto end = meshes_.end();
        begin++;

        std::for_each(
            begin,
            end,
            [&joined_mesh](const auto& m){
                joined_mesh += m;
            }
        );

        return joined_mesh;
    };

    template <class K>
    void Collision_scene<K>::color(const Scene_face_index& ti, CGAL::IO::Color c) {
        meshes_[ti.mesh_index().idx()].color(ti.local_index(), c);
        return;
    };
    
    template <class K>
    template <class UpdateFunctor>
    void Collision_scene<K>::update_state(UpdateFunctor& update_functor, bool update_tree_)
    {
        Mesh* mesh_ptr = nullptr;
        for( auto& scene_vertex_index : vertices_) 
        {
            mesh_ptr = &(meshes_[scene_vertex_index.mesh_index()]);
            update_functor(mesh_ptr, scene_vertex_index);
        }
        
        if( update_tree_ ) {
            update_tree();
        }
    }
    
    template <class K>
    void Collision_scene<K>::update_tree()
    {
        for( auto& t : trajectories_ ) 
        { 
            t.update(); 
        }
        // TODO: use the update functor that sloriot provided
        trajectory_tree_ = Tree(trajectories_.begin(),trajectories_.end());
    }


}

#endif