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

#ifndef DRAW_COLLISION_SCENE_H
#define DRAW_COLLISION_SCENE_H

#include <iostream>
#include <set>
#include <vector>

#include <Update_functors.h>
#include <CGAL/draw_surface_mesh.h>

namespace CGAL {

template <class CollisionScene>
void draw_collision_scene( CollisionScene& scene, bool draw_next=false )
{
    using K             = typename CollisionScene::Kernel;
    using Swap_functor  = Swap_current_next_functor<CollisionScene>;
    
    for( const auto& fi : scene.faces())
    {
        scene.color(fi, CGAL::IO::blue());
    }

    if( draw_next ) {

        // Color collisions white
        auto collisions = CGAL::get_collisions(scene);
        auto collision_indices = convert_candidates_to_indices(collisions);
        for( const auto& fi : collision_indices )
        {
            scene.color(fi, CGAL::IO::white());
        }

        Swap_functor scf = Swap_functor();
        scene.update_state(scf);
        ::CGAL::draw_color(scene.joined_meshes());
        scene.update_state(scf);
    }
    else {
        ::CGAL::draw_color(scene.joined_meshes());
    }

}

template <class CollisionCandidate, class Index=typename CollisionCandidate::Index>
std::set<Index>
convert_candidates_to_indices(
    std::vector<CollisionCandidate> candidates
)
{
    std::set<typename CollisionCandidate::Index> collision_indices;
    std::for_each(
        candidates.begin(), 
        candidates.end(), 
        [&collision_indices](const auto& candidate){ 
            collision_indices.insert(candidate.first->index);
            collision_indices.insert(candidate.second->index);
        }
    );
    return collision_indices;
}

}

#endif