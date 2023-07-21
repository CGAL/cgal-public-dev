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
            typedef typename K::Point_3                                                 Point;
            typedef typename K::Vector_3                                                Vector;
            typedef          Collision_mesh<K>                                          CMesh;
            typedef typename CMesh::Vertex_index                                        Vertex_index;
            typedef typename CMesh::Edge_index                                          Edge_index;
            typedef typename CMesh::Halfedge_index                                      Halfedge_index;
            typedef typename CMesh::Face_index                                          Face_index;
            typedef typename CMesh::Face_range                                          Face_range;
            typedef          std::pair<int, Face_index>                                 Trajectory_index;
            typedef typename CMesh::Vertex_around_face_range                            Vertex_around_face_range;
            typedef          Triangle_trajectory<K, Trajectory_index>                   Trajectory;
            typedef          AABB_Triangle_trajectory_primitive<K, Trajectory_index>    Trajectory_primitive;
            typedef          ::CGAL::AABB_traits<K, Trajectory_primitive>               AABB_traits;
            typedef          AABB_tree<AABB_traits>                                     Tree;
            typedef typename Tree::Primitive_id                                         Primitive_id;

            std::vector<CMesh>&     meshes_;
            std::vector<Trajectory> trajectories_;
            Tree                    trajectory_tree;

        private:

            void append_trajectories(const CMesh& mesh, const int i);

        public:

            Collision_scene(std::vector<CMesh> & meshes) : meshes_{meshes}
            {
                int i{0};
                for( const auto & m : meshes_ )
                {
                    append_trajectories(m, i);
                    ++i;
                }

                trajectory_tree = Tree(trajectories_.begin(),trajectories_.end());
            }

    };

    template <class K>
    void Collision_scene<K>::append_trajectories(const Collision_mesh<K>& mesh, const int mesh_index)
    {

        std::vector<Trajectory> tmp = get_trajectories( mesh, mesh_index );

        trajectories_.reserve(trajectories_.size() + tmp.size());
        trajectories_.insert(trajectories_.end(), tmp.begin(), tmp.end());

        return;
    }

}

#endif