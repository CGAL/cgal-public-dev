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

#include <Collision_mesh_3.h>
#include <AABB_triangle_trajectory_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <vector>
#include <utility>

namespace CGAL {

    template <class K>
    class Collision_scene {

        public:
            typedef typename K::Point_3                     Point;
            typedef typename K::Vector_3                    Vector;
            typedef          Collision_mesh<K>              Mesh;
            typedef typename Mesh::Vertex_index             Vertex_index;
            typedef typename Mesh::Edge_index               Edge_index;
            typedef typename Mesh::Halfedge_index           Halfedge_index;
            typedef typename Mesh::Face_index               Face_index;
            typedef typename Mesh::Face_range               Face_range;
            typedef typename Mesh::Vertex_around_face_range Vertex_around_face_range;

            struct Mesh_index {

                std::size_t idx_;

                explicit Mesh_index(std::size_t idx) : idx_{ idx } {}
                
                friend std::ostream& operator<<(std::ostream& os, Mesh_index const& mesh_index)
                {
                    return (os << "m" << mesh_index.idx_ );
                } 

            };

            struct Trajectory_index {
                Face_index face_index;
                Mesh_index mesh_index;

                Trajectory_index(Mesh_index mi, Face_index face_index) : mesh_index{ mi }, face_index{ face_index }  {}
                
                friend std::ostream& operator<<(std::ostream& os, Trajectory_index const& trajectory_index)
                {
                    return (os << "(" << trajectory_index.mesh_index << ", " << trajectory_index.face_index << ")" );
                } 

            };

            typedef          Triangle_trajectory<K, Trajectory_index>                   Trajectory;
            typedef          AABB_Triangle_trajectory_primitive<K, Trajectory_index>    Trajectory_primitive;
            typedef          ::CGAL::AABB_traits<K, Trajectory_primitive>               AABB_traits;
            typedef          AABB_tree<AABB_traits>                                     Tree;
            typedef typename Tree::Primitive_id                                         Primitive_id;

            std::vector<Mesh>&     meshes_;
            std::vector<Trajectory> trajectories_;
            Tree                    trajectory_tree;

        private:

            std::vector<Trajectory> get_trajectories(const Mesh & mesh, const Mesh_index& mesh_index);

        public:

            explicit Collision_scene(std::vector<Mesh> & meshes) : meshes_{meshes}
            {
                int i{0};
                for( const auto & mesh_ : meshes_ )
                {
                    std::vector<Trajectory> tmp = get_trajectories(mesh_, Mesh_index(i));
                    trajectories_.reserve(trajectories_.size() + tmp.size());
                    trajectories_.insert(trajectories_.end(), tmp.begin(), tmp.end());
                    ++i;
                }

                trajectory_tree = Tree(trajectories_.begin(),trajectories_.end());
            }

    };
    
    template <class K>
    std::vector<typename Collision_scene<K>::Trajectory> Collision_scene<K>::get_trajectories( 
        const Mesh& mesh, 
        const Mesh_index& mesh_index 
    ){

        std::vector<Trajectory> trajectories;
        trajectories.reserve(mesh.num_faces());

        int j;
        std::vector<const Point*> trajectory_points(6);
        for( const auto& face_index : mesh.faces())
        {
            j = 0;

            for( const auto & vert_index : mesh.vertices_around_face(mesh.halfedge(face_index)) )
            {
                trajectory_points.at(j) = & mesh.point(vert_index);
                trajectory_points.at(j+1) = & mesh.next_point(vert_index);

                j += 2;
            }

            trajectories.push_back(
                Trajectory(
                    trajectory_points[0],
                    trajectory_points[2],
                    trajectory_points[4],
                    trajectory_points[1],
                    trajectory_points[3],
                    trajectory_points[5],
                    Trajectory_index(mesh_index, face_index)
                )
            );

        }

        return trajectories;
    };

}

#endif